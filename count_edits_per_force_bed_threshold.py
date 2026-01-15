#!/usr/bin/env python3
import argparse, sys, re, csv
from pathlib import Path
import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree

def eprint(*a): print(*a, file=sys.stderr)

def normalize_chrom(chrom: str, contigs_set: set[str]) -> str | None:
    chrom = str(chrom).strip()
    if chrom in contigs_set:
        return chrom
    if chrom.startswith("chr") and chrom[3:] in contigs_set:
        return chrom[3:]
    if ("chr" + chrom) in contigs_set:
        return "chr" + chrom
    return None

def load_fai_contigs(fai: Path):
    contigs = []
    with fai.open() as f:
        for line in f:
            if line.strip():
                contigs.append(line.split("\t", 1)[0].strip())
    return set(contigs), contigs

def detect_split(header_line: str):
    return ("\t" if "\t" in header_line else None)  # None => split on whitespace

def read_force_to_trees(force_tsv: Path, contigs_set: set[str], force_one_based: bool):
    """
    Builds IntervalTrees keyed by chrom.
    Expects columns (names can vary but must exist):
      chr/start/end/max_mean_dsRNA_force/class/family/rep.id (or rep_id)
    Coordinates:
      - If force_one_based=True: treats start/end as 1-based inclusive and converts to 0-based half-open.
      - Else: treats start/end as already 0-based half-open (BED-like).
    """
    trees: dict[str, IntervalTree] = {}

    with force_tsv.open() as f:
        header = f.readline().rstrip("\n")
        if not header:
            raise SystemExit(f"[FATAL] Force TSV empty: {force_tsv}")
        sep = detect_split(header)
        cols = header.split(sep) if sep else header.split()
        cols_lc = [c.strip().lower() for c in cols]

        def find(*names):
            for nm in names:
                nm = nm.lower()
                if nm in cols_lc:
                    return cols_lc.index(nm)
            return None

        ci = find("chr", "chrom", "chromosome", "seqname")
        si = find("start", "chromstart")
        ei = find("end", "chromend")
        fi = find("max_mean_dsrna_force", "max_mean_dsRNA_force".lower(), "dsrna_force", "force")
        cls_i = find("class")
        fam_i = find("family")
        rep_i = find("rep.id", "rep_id", "repid", "repeat", "name")

        missing = [k for k,v in {
            "chr":ci,"start":si,"end":ei,"max_mean_dsRNA_force":fi,"class":cls_i,"family":fam_i,"rep_id":rep_i
        }.items() if v is None]
        if missing:
            raise SystemExit(f"[FATAL] Force TSV missing columns {missing}. Header={cols}")

        n_kept = 0
        for line in f:
            if not line.strip(): continue
            p = (line.rstrip("\n").split(sep) if sep else line.split())
            if len(p) <= max(ci,si,ei,fi,cls_i,fam_i,rep_i):
                continue

            chrom = normalize_chrom(p[ci], contigs_set)
            if chrom is None:
                continue

            try:
                s = int(float(p[si])); e = int(float(p[ei]))
                force = float(p[fi])
            except:
                continue

            if force_one_based:
                # 1-based inclusive -> 0-based half-open
                s0 = s - 1
                e0 = e
            else:
                # already BED-like
                s0 = s
                e0 = e

            if e0 <= s0:
                continue

            cls = p[cls_i].strip() or "NA"
            fam = p[fam_i].strip() or "NA"
            rep = p[rep_i].strip() or "NA"

            tree = trees.get(chrom)
            if tree is None:
                tree = IntervalTree()
                trees[chrom] = tree
            tree.add(Interval(s0, e0, (force, cls, fam, rep)))
            n_kept += 1

    eprint(f"[INFO] Loaded force intervals: {n_kept}")
    return trees

def iter_reditools_chunks(path: Path, chunksize: int):
    return pd.read_csv(path, sep="\t", dtype=str, chunksize=chunksize)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force_tsv", required=True)
    ap.add_argument("--reditools_table", required=True)
    ap.add_argument("--genome_fai", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--min_cov", type=float, default=None)
    ap.add_argument("--min_frequency", type=float, default=None)
    ap.add_argument("--chunksize", type=int, default=500000)
    ap.add_argument("--force_one_based", action="store_true",
                    help="If force TSV start/end are 1-based inclusive, convert to BED. Try this if you get 0 overlaps.")
    args = ap.parse_args()

    contigs_set, _ = load_fai_contigs(Path(args.genome_fai))
    trees = read_force_to_trees(Path(args.force_tsv), contigs_set, args.force_one_based)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # We will write CSV with proper quoting (BaseCount has commas!)
    header = [
        "Region","Position","Reference","Strand","Coverage-q30","MeanQ",
        "BaseCount[A,C,G,T]","AllSubs","Frequency",
        "max_mean_dsRNA_force","Class","Family","rep_id"
    ]

    wrote_header = False
    n_out = 0

    with out_csv.open("w", newline="") as wf:
        writer = csv.writer(wf)

        for chunk in iter_reditools_chunks(Path(args.reditools_table), args.chunksize):
            cols = {c.lower(): c for c in chunk.columns}

            need = ["region","position","reference","strand","coverage-q30","meanq","allsubs","frequency"]
            missing = [x for x in need if x not in cols]
            if missing:
                raise SystemExit(f"[FATAL] REDItools missing columns {missing}. Found={list(chunk.columns)}")

            bc_col = next((c for c in chunk.columns if "basecount" in c.lower()), None)
            if bc_col is None:
                raise SystemExit("[FATAL] Could not find BaseCount column.")

            region_norm = chunk[cols["region"]].map(lambda x: normalize_chrom(x, contigs_set))
            pos1 = pd.to_numeric(chunk[cols["position"]], errors="coerce")
            cov  = pd.to_numeric(chunk[cols["coverage-q30"]], errors="coerce")
            freq = pd.to_numeric(chunk[cols["frequency"]], errors="coerce")

            ok = region_norm.notna() & pos1.notna()
            if args.min_cov is not None:
                ok = ok & cov.notna() & (cov >= args.min_cov)
            if args.min_frequency is not None:
                ok = ok & freq.notna() & (freq >= args.min_frequency)

            sub = chunk.loc[ok].copy()
            if sub.empty:
                continue

            sub["_chr"] = region_norm[ok]
            sub["_pos1"] = pos1[ok].astype(int)
            # 0-based coordinate for point query
            sub["_p0"] = sub["_pos1"] - 1

            if not wrote_header:
                writer.writerow(header)
                wrote_header = True

            # For each row: query intervals containing point, keep max force
            for _, r in sub.iterrows():
                chrom = r["_chr"]
                p0 = int(r["_p0"])
                tree = trees.get(chrom)
                if tree is None:
                    continue
                hits = tree.at(p0)
                if not hits:
                    continue

                best = None  # (force, cls, fam, rep)
                for itv in hits:
                    force, cls, fam, rep = itv.data
                    if best is None or force > best[0]:
                        best = (force, cls, fam, rep)

                force, cls, fam, rep = best
                writer.writerow([
                    chrom,
                    str(int(r["_pos1"])),
                    r[cols["reference"]],
                    r[cols["strand"]],
                    r[cols["coverage-q30"]],
                    r[cols["meanq"]],
                    r[bc_col],
                    r[cols["allsubs"]],
                    r[cols["frequency"]],
                    f"{force:.6g}",
                    cls, fam, rep
                ])
                n_out += 1

    if not wrote_header:
        raise SystemExit("[FATAL] No output rows written (0 overlaps). Try --force_one_based or check chr naming/build.")

    eprint(f"[OK] Wrote {n_out} overlapping sites to {out_csv}")

if __name__ == "__main__":
    main()

