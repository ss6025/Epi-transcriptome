#!/usr/bin/env python3
"""
intersect_te.py — Intersect A-to-I editing sites with a loci BED file.

Reads all *.AGTC.bed files from --agtc-dir, intersects with --loci-bed
(TEs, known ADAR targets, repeat elements, or any genomic regions), and
writes per-sample and combined output tables.

Optional --metadata CSV (sample, genotype, condition, treatment) enables
subsetting and stratified summary plots.

Usage:
  python3 intersect_te.py \\
      --agtc-dir  filtered_agtc/ \\
      --loci-bed  mm10_te.bed \\
      --output    te_intersect/ \\
      [--metadata sample_metadata.csv] \\
      [--genotype WT,KO] \\
      [--condition tumor,immune] \\
      [--treatment IFNAR,IgG] \\
      [--loci-one-based] \\
      [--threads 4]
"""
from __future__ import annotations

import argparse
import bisect
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd


# ============================================================
# Interval helpers (pure Python, no bedtools dependency)
# ============================================================

def merge_intervals(ivs: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not ivs:
        return []
    ivs = sorted(ivs)
    out = [ivs[0]]
    for s, e in ivs[1:]:
        if s < out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def site_overlaps(start: int, end: int, merged: list[tuple[int, int]]) -> bool:
    """Binary-search test: does [start, end) overlap any interval in *sorted* merged list?"""
    idx = bisect.bisect_right(merged, (end, end)) - 1
    if idx < 0:
        return False
    while idx >= 0 and merged[idx][0] < end:
        if merged[idx][1] > start:
            return True
        idx -= 1
    return False


def overlap_width(s1: int, e1: int, s2: int, e2: int) -> int:
    return max(0, min(e1, e2) - max(s1, s2))


# ============================================================
# Chromosome name normalisation
# ============================================================

def norm_chr(c) -> str:
    c = str(c).strip()
    if c.startswith("chr"):
        return c
    if c.isdigit() or c in ("X", "Y", "M", "MT"):
        return "chr" + c
    return c


# ============================================================
# Load loci BED
# ============================================================

def load_loci_bed(
    path: Path, one_based: bool = False
) -> dict[str, list[tuple[int, int, str]]]:
    """
    Load a BED-like loci file. Columns expected: chrom, start, end[, name, ...].
    Returns { chrom: [(start0, end0, name), ...] } sorted by start.
    Supports plain BED (0-based half-open) or 1-based coords (--loci-one-based).
    Also handles SQuIRE output with TE_chr / TE_start / TE_stop / TE_name columns.
    """
    comp = "gzip" if str(path).endswith(".gz") else None
    try:
        df = pd.read_csv(path, sep="\t", comment="#", low_memory=False, compression=comp)
    except Exception as exc:
        sys.exit(f"[FATAL] Cannot read loci BED {path}: {exc}")

    df.columns = [str(c).strip() for c in df.columns]
    cols_lc = {c.lower(): c for c in df.columns}

    # SQuIRE layout
    squire = all(k in cols_lc for k in ("te_chr", "te_start", "te_stop"))
    if squire:
        c_chr  = cols_lc["te_chr"]
        c_s    = cols_lc["te_start"]
        c_e    = cols_lc["te_stop"]
        c_nm   = cols_lc.get("te_name", cols_lc.get("te_class", None))
    else:
        # Standard BED: first three columns
        c_chr, c_s, c_e = df.columns[0], df.columns[1], df.columns[2]
        c_nm = df.columns[3] if len(df.columns) > 3 else None

    by_chrom: dict[str, list[tuple[int, int, str]]] = {}
    for i in range(len(df)):
        ch = norm_chr(df[c_chr].iloc[i])
        if not ch:
            continue
        try:
            ts = int(pd.to_numeric(df[c_s].iloc[i], errors="coerce"))
            te = int(pd.to_numeric(df[c_e].iloc[i], errors="coerce"))
        except (TypeError, ValueError):
            continue
        bs, be = (ts - 1, te) if one_based else (ts, te)
        if bs < 0 or be <= bs:
            continue
        name = str(df[c_nm].iloc[i]).strip() if c_nm else "locus"
        by_chrom.setdefault(ch, []).append((bs, be, name))

    for ch in by_chrom:
        by_chrom[ch].sort(key=lambda x: (x[0], x[1]))
    return by_chrom


# ============================================================
# Load AGTC BED
# ============================================================

_AGTC_COLS = ["Chrom", "Start", "End", "Name", "Score", "Strand",
               "AllSubs", "Coverage", "Frequency", "CellType"]


def load_agtc_bed(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep="\t", low_memory=False)
    except Exception as exc:
        print(f"[WARN] Cannot read {path}: {exc}")
        return pd.DataFrame()
    df.columns = [str(c).strip() for c in df.columns]
    # Normalise column names
    rename = {}
    for col in df.columns:
        cl = col.lower().replace("-", "").replace("_", "")
        if cl in ("chrom", "chr", "region", "chromosome"):
            rename[col] = "Chrom"
        elif cl in ("start",):
            rename[col] = "Start"
        elif cl in ("end", "stop"):
            rename[col] = "End"
        elif cl in ("coverage", "coverageq30"):
            rename[col] = "Coverage"
        elif cl in ("frequency", "freq"):
            rename[col] = "Frequency"
        elif cl in ("allsubs", "subs", "sub"):
            rename[col] = "AllSubs"
        elif cl in ("celltype", "cell_type", "cell"):
            rename[col] = "CellType"
    df = df.rename(columns=rename)
    df["_c"] = df["Chrom"].apply(norm_chr)
    df["_s"] = pd.to_numeric(df.get("Start", 0), errors="coerce").fillna(0).astype(int)
    df["_e"] = pd.to_numeric(df.get("End",   0), errors="coerce").fillna(0).astype(int)
    return df


# ============================================================
# Intersection
# ============================================================

def classify_sites(
    agtc: pd.DataFrame,
    loci_by_chrom: dict[str, list[tuple[int, int, str]]],
) -> tuple[pd.Series, pd.Series]:
    """Return (in_loci: bool Series, locus_name: str Series)."""
    merged_by = {
        ch: merge_intervals([(t[0], t[1]) for t in rows])
        for ch, rows in loci_by_chrom.items()
    }
    in_loci: list[bool] = []
    names:   list[str]  = []

    for i in range(len(agtc)):
        ch = agtc["_c"].iloc[i]
        s  = int(agtc["_s"].iloc[i])
        e  = int(agtc["_e"].iloc[i])
        loci = loci_by_chrom.get(ch)
        if not loci:
            in_loci.append(False)
            names.append("nonLoci")
            continue
        merged = merged_by[ch]
        if not site_overlaps(s, e, merged):
            in_loci.append(False)
            names.append("nonLoci")
            continue
        # Find the best-overlapping named locus
        best_w, best_name = -1, "Loci"
        for lo, hi, nm in loci:
            if lo >= e or hi <= s:
                continue
            w = overlap_width(s, e, lo, hi)
            if w > best_w:
                best_w, best_name = w, nm
        in_loci.append(True)
        names.append(best_name if best_w > 0 else "Loci")

    return (
        pd.Series(in_loci, index=agtc.index, dtype=bool),
        pd.Series(names,   index=agtc.index, dtype=str),
    )


# ============================================================
# Process a single AGTC BED file
# ============================================================

def process_one(args_tuple) -> dict:
    agtc_path, loci_by_chrom, output_dir = args_tuple
    sample = agtc_path.stem.replace(".AGTC", "")
    out_dir = Path(output_dir)

    df = load_agtc_bed(agtc_path)
    if df.empty:
        return {"sample": sample, "total": 0, "in_loci": 0, "not_loci": 0}

    in_loci, locus_names = classify_sites(df, loci_by_chrom)
    df["InLoci"]    = in_loci
    df["LocusName"] = locus_names

    # Write all sites with loci annotation
    all_out = out_dir / f"{sample}.loci_annotated.bed"
    df.drop(columns=["_c", "_s", "_e"], errors="ignore").to_csv(
        all_out, sep="\t", index=False
    )

    # Write only the in-loci sites
    loci_out = out_dir / f"{sample}.in_loci.bed"
    df[df["InLoci"]].drop(columns=["_c", "_s", "_e"], errors="ignore").to_csv(
        loci_out, sep="\t", index=False
    )

    n_in  = int(in_loci.sum())
    n_tot = len(df)
    print(f"[OK] {sample}: {n_in}/{n_tot} sites in loci -> {loci_out.name}")
    return {"sample": sample, "total": n_tot, "in_loci": n_in, "not_loci": n_tot - n_in}


# ============================================================
# Metadata subsetting
# ============================================================

def load_sample_metadata(path: Path) -> pd.DataFrame:
    """CSV: sample, genotype, condition, treatment (header required)."""
    try:
        df = pd.read_csv(path)
        df.columns = [c.strip().lower() for c in df.columns]
        if "sample" not in df.columns:
            print(f"[WARN] metadata missing 'sample' column — subsetting disabled")
            return pd.DataFrame()
        return df
    except Exception as exc:
        print(f"[WARN] Cannot read metadata {path}: {exc}")
        return pd.DataFrame()


def filter_samples_by_metadata(
    agtc_files: list[Path],
    meta: pd.DataFrame,
    genotypes: list[str] | None,
    conditions: list[str] | None,
    treatments: list[str] | None,
) -> list[Path]:
    if meta.empty:
        return agtc_files
    keep_samples = set(meta["sample"].astype(str))
    if genotypes and "genotype" in meta.columns:
        g = meta[meta["genotype"].isin(genotypes)]["sample"].astype(str)
        keep_samples &= set(g)
    if conditions and "condition" in meta.columns:
        c = meta[meta["condition"].isin(conditions)]["sample"].astype(str)
        keep_samples &= set(c)
    if treatments and "treatment" in meta.columns:
        t = meta[meta["treatment"].isin(treatments)]["sample"].astype(str)
        keep_samples &= set(t)
    filtered = [
        f for f in agtc_files
        if any(s in f.stem for s in keep_samples)
    ]
    print(f"[INFO] Metadata filter: {len(filtered)}/{len(agtc_files)} samples kept")
    return filtered


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Intersect AGTC editing sites with a loci BED file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--agtc-dir",  required=True, type=Path,
                        help="Directory containing *.AGTC.bed files (from step 4)")
    parser.add_argument("--loci-bed",  required=True, type=Path,
                        help="Loci BED file (TEs, ADAR targets, repeats, etc.)")
    parser.add_argument("--output",    required=True, type=Path,
                        help="Output directory")
    parser.add_argument("--metadata",  type=Path, default=None,
                        help="Sample metadata CSV: sample,genotype,condition,treatment")
    parser.add_argument("--genotype",  type=str, default=None,
                        help="Comma-separated genotypes to include (e.g. WT,KO)")
    parser.add_argument("--condition", type=str, default=None,
                        help="Comma-separated conditions to include (e.g. tumor,immune)")
    parser.add_argument("--treatment", type=str, default=None,
                        help="Comma-separated treatments to include (e.g. IFNAR,IgG)")
    parser.add_argument("--loci-one-based", action="store_true",
                        help="Loci BED uses 1-based coordinates (default: 0-based)")
    parser.add_argument("--threads",   type=int, default=4)
    args = parser.parse_args()

    if not args.agtc_dir.is_dir():
        sys.exit(f"[FATAL] --agtc-dir not found: {args.agtc_dir}")
    if not args.loci_bed.exists():
        sys.exit(f"[FATAL] --loci-bed not found: {args.loci_bed}")

    args.output.mkdir(parents=True, exist_ok=True)

    # Parse comma-separated filter lists
    genotypes  = [g.strip() for g in args.genotype.split(",")]  if args.genotype  else None
    conditions = [c.strip() for c in args.condition.split(",")] if args.condition else None
    treatments = [t.strip() for t in args.treatment.split(",")] if args.treatment else None

    # Load loci
    print(f"[INFO] Loading loci BED: {args.loci_bed}")
    loci_by_chrom = load_loci_bed(args.loci_bed, one_based=args.loci_one_based)
    n_chroms = len(loci_by_chrom)
    n_loci   = sum(len(v) for v in loci_by_chrom.values())
    print(f"[INFO] Loaded {n_loci:,} loci across {n_chroms} chromosomes")

    # Discover AGTC BED files
    agtc_files = sorted(args.agtc_dir.glob("*.AGTC.bed"))
    if not agtc_files:
        sys.exit(f"[FATAL] No *.AGTC.bed files in {args.agtc_dir}")
    print(f"[INFO] Found {len(agtc_files)} AGTC BED file(s)")

    # Optional metadata subsetting
    meta = pd.DataFrame()
    if args.metadata and args.metadata.exists():
        meta = load_sample_metadata(args.metadata)
    if not meta.empty or any([genotypes, conditions, treatments]):
        agtc_files = filter_samples_by_metadata(
            agtc_files, meta, genotypes, conditions, treatments
        )

    if not agtc_files:
        sys.exit("[FATAL] No AGTC files remain after metadata filtering")

    # Process in parallel
    task_args = [(f, loci_by_chrom, args.output) for f in agtc_files]
    results = []
    with ProcessPoolExecutor(max_workers=args.threads) as pool:
        for r in pool.map(process_one, task_args):
            results.append(r)

    # Summary table
    summary_df = pd.DataFrame(results).sort_values("sample")
    summary_df["pct_in_loci"] = (
        summary_df["in_loci"] / summary_df["total"].replace(0, float("nan")) * 100
    ).round(2)

    # Attach metadata columns if available
    if not meta.empty:
        summary_df = summary_df.merge(meta, on="sample", how="left")

    summary_path = args.output / "intersection_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\n[OK] Summary -> {summary_path}")

    # Combined in-loci BED across all samples
    combined_parts = []
    for f in agtc_files:
        sample = f.stem.replace(".AGTC", "")
        loci_out = args.output / f"{sample}.in_loci.bed"
        if loci_out.exists() and loci_out.stat().st_size > 0:
            part = pd.read_csv(loci_out, sep="\t", low_memory=False)
            if not part.empty:
                part["Sample"] = sample
                if not meta.empty and "sample" in meta.columns:
                    row = meta[meta["sample"] == sample]
                    for col in ["genotype", "condition", "treatment"]:
                        if col in meta.columns:
                            part[col.capitalize()] = (
                                row[col].values[0] if len(row) > 0 else "NA"
                            )
                combined_parts.append(part)

    if combined_parts:
        combined = pd.concat(combined_parts, ignore_index=True)
        combined_path = args.output / "all_samples_in_loci.bed"
        combined.to_csv(combined_path, sep="\t", index=False)
        print(f"[OK] Combined in-loci BED -> {combined_path}")
        print(f"     Total sites: {len(combined):,}")

    # Print summary
    print("\n--- Intersection Summary ---")
    print(summary_df[["sample", "total", "in_loci", "pct_in_loci"]].to_string(index=False))


if __name__ == "__main__":
    main()
