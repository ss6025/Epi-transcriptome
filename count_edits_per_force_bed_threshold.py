#!/usr/bin/env python3
import argparse
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional, Tuple, List, Set

import math

# ----------------------------- utils -----------------------------

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def run(cmd: List[str], stdout_path: Optional[Path] = None) -> None:
    if stdout_path is None:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            eprint("CMD FAILED:\n  " + " ".join(map(str, cmd)))
            eprint(p.stderr)
            raise RuntimeError("Command failed")
        return
    else:
        with stdout_path.open("w") as w:
            p = subprocess.run(cmd, stdout=w, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            eprint("CMD FAILED:\n  " + " ".join(map(str, cmd)))
            eprint(p.stderr)
            raise RuntimeError("Command failed")
        return

def which_or_die(prog: str):
    from shutil import which
    if which(prog) is None:
        raise SystemExit(f"[FATAL] Required program not found in PATH: {prog}")

def sanitize_chrom(s: str) -> str:
    return s.replace("\r", "").strip()

def load_fai_contigs_and_rank(genome_fai: Path) -> Tuple[Set[str], List[str]]:
    contigs = []
    with genome_fai.open() as f:
        for line in f:
            if not line.strip():
                continue
            c = line.split("\t", 1)[0].strip()
            if c:
                contigs.append(c)
    if not contigs:
        raise SystemExit(f"[FATAL] No contigs read from fai: {genome_fai}")
    return set(contigs), contigs

# ---------------------- REDItools sites -> BEDs -------------------

def detect_header_and_delim(path: Path) -> Tuple[bool, str, Optional[List[str]]]:
    with path.open() as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            delim = "\t" if line.count("\t") >= line.count(",") else ","
            fields = line.rstrip("\n").split(delim)
            lower = [x.strip().lower() for x in fields]
            if ("region" in lower and "position" in lower) or ("chr" in lower and "start" in lower and "end" in lower):
                return True, delim, fields
            return False, delim, None
    return False, "\t", None

def looks_like_int(x: str) -> bool:
    try:
        int(x)
        return True
    except Exception:
        return False

_sub_re = re.compile(r"^(?:[ACGT][ACGT]|-|\.)$", re.IGNORECASE)

def infer_site_schema_headerless(path: Path, delim: str, max_lines: int = 2000):
    coord_conf_bed = 0
    coord_conf_pos = 0
    sub_counts = []

    with path.open() as f:
        n = 0
        for line in f:
            if n >= max_lines:
                break
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(delim)
            n += 1

            if len(parts) >= 3 and looks_like_int(parts[1]) and looks_like_int(parts[2]):
                s = int(parts[1]); e = int(parts[2])
                if e > s:
                    coord_conf_bed += 1
            if len(parts) >= 2 and looks_like_int(parts[1]):
                coord_conf_pos += 1

            while len(sub_counts) < len(parts):
                sub_counts.append(0)

            for i, v in enumerate(parts):
                v = v.strip()
                if _sub_re.match(v):
                    sub_counts[i] += 1

    coord_mode = "bed" if coord_conf_bed >= max(10, coord_conf_pos // 2) else "pos"

    sub_idx = None
    candidates = [(c, i) for i, c in enumerate(sub_counts) if i not in (0, 1, 2)]
    if candidates:
        candidates.sort(reverse=True)
        if candidates[0][0] >= 10:
            sub_idx = candidates[0][1]

    return coord_mode, sub_idx

def write_sites_beds_filtered_any_ag(
    reditools_path: Path,
    contigs: Set[str],
    out_any: Path,
    out_ag: Path,
) -> None:
    has_header, delim, header = detect_header_and_delim(reditools_path)

    kept_any = kept_ag = 0
    dropped_contig = 0

    if has_header and header is not None:
        idx = {h.strip().lower(): i for i, h in enumerate(header)}

        if "region" in idx and "position" in idx:
            i_chr = idx["region"]
            i_pos = idx["position"]
            i_start = i_end = None
        elif "chr" in idx and "start" in idx and "end" in idx:
            i_chr = idx["chr"]; i_start = idx["start"]; i_end = idx["end"]
            i_pos = None
        else:
            raise SystemExit(f"[FATAL] Headered input but cannot find Region/Position or chr/start/end in: {header}")

        i_sub = idx.get("allsubs", None)
        if i_sub is None:
            raise SystemExit(f"[FATAL] Need AllSubs column in headered REDItools table. Columns: {header}")

        with reditools_path.open() as f, out_any.open("w") as wany, out_ag.open("w") as wag:
            first = True
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                if first:
                    first = False
                    continue
                parts = line.rstrip("\n").split(delim)
                try:
                    chrom = sanitize_chrom(parts[i_chr])
                    if chrom not in contigs:
                        dropped_contig += 1
                        continue

                    if i_pos is not None:
                        pos1 = int(float(parts[i_pos]))
                        start = pos1 - 1
                        end = pos1
                    else:
                        start = int(float(parts[i_start]))
                        end = int(float(parts[i_end]))
                    if end <= start:
                        continue
                except Exception:
                    continue

                sub = parts[i_sub].strip() if i_sub < len(parts) else "-"
                edited_any = (sub not in ("-", ".", ""))
                edited_ag = ("AG" in sub.upper())

                if edited_any:
                    wany.write(f"{chrom}\t{start}\t{end}\n")
                    kept_any += 1
                if edited_ag:
                    wag.write(f"{chrom}\t{start}\t{end}\n")
                    kept_ag += 1

        eprint(f"[INFO] SITES contig filter: kept_any {kept_any:,}, kept_ag {kept_ag:,}, dropped {dropped_contig:,} not in .fai")
        return

    # headerless fallback
    coord_mode, sub_idx = infer_site_schema_headerless(reditools_path, delim)
    if sub_idx is None:
        raise SystemExit("[FATAL] Could not infer AllSubs-like column in headerless REDItools file.")

    with reditools_path.open() as f, out_any.open("w") as wany, out_ag.open("w") as wag:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split(delim)
            if len(parts) < 2:
                continue

            chrom = sanitize_chrom(parts[0])
            if chrom not in contigs:
                dropped_contig += 1
                continue

            try:
                if coord_mode == "bed":
                    if len(parts) < 3:
                        continue
                    start = int(parts[1])
                    end = int(parts[2])
                else:
                    pos1 = int(parts[1])
                    start = pos1 - 1
                    end = pos1
                if end <= start:
                    continue
            except Exception:
                continue

            sub = parts[sub_idx].strip() if sub_idx < len(parts) else "-"
            edited_any = (sub not in ("-", ".", ""))
            edited_ag = ("AG" in sub.upper())

            if edited_any:
                wany.write(f"{chrom}\t{start}\t{end}\n")
                kept_any += 1
            if edited_ag:
                wag.write(f"{chrom}\t{start}\t{end}\n")
                kept_ag += 1

    eprint(f"[INFO] SITES contig filter: kept_any {kept_any:,}, kept_ag {kept_ag:,}, dropped {dropped_contig:,} not in .fai")

# --------------------- sorting helpers ---------------------------

def write_contig_rank_tsv(contig_order: List[str], out_path: Path) -> None:
    with out_path.open("w") as w:
        for i, c in enumerate(contig_order, start=1):
            w.write(f"{c}\t{i}\n")

def sort_bed_by_fai_rank(in_bed: Path, out_bed: Path, rank_tsv: Path) -> None:
    cmd = [
        "bash", "-lc",
        (
            "set -euo pipefail; "
            "LC_ALL=C "
            f"awk 'BEGIN{{FS=OFS=\"\\t\"}} "
            f"FNR==NR{{rank[$1]=$2; next}} "
            "{c=$1; if(!(c in rank)) next; "
            "r=rank[c]; "
            "printf \"%s\\t%s\\t%s\\t%s\", r,$1,$2,$3; "
            "for(i=4;i<=NF;i++) printf \"\\t%s\", $i; "
            "printf \"\\n\" }' "
            f"{rank_tsv} {in_bed} "
            "| LC_ALL=C sort -k1,1n -k3,3n -k4,4n "
            "| cut -f2- "
            f"> {out_bed}"
        )
    ]
    run(cmd)

# -------------------- force BED (no header) ----------------------

def parse_force_bed_col13_to_bed_filtered(
    force_bed_in: Path,
    out_bed: Path,
    contigs: Set[str],
    force_col_1based: int = 13,
) -> None:
    """
    Input: BED-like, no header. Assumes:
      col1=chr col2=start col3=end
      col13 (1-based) = dsRNA_force
    Output: BED with columns:
      chr start end dsRNA_force length
    """
    fi = force_col_1based - 1
    kept = dropped = bad = 0
    with force_bed_in.open() as f, out_bed.open("w") as w:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                bad += 1
                continue
            chrom = sanitize_chrom(p[0])
            if chrom not in contigs:
                dropped += 1
                continue
            try:
                start = int(float(p[1]))
                end   = int(float(p[2]))
            except Exception:
                bad += 1
                continue
            if end <= start:
                bad += 1
                continue

            forcev = "NA"
            if len(p) > fi:
                try:
                    forcev = str(float(p[fi]))
                except Exception:
                    forcev = "NA"

            length = end - start
            w.write(f"{chrom}\t{start}\t{end}\t{forcev}\t{length}\n")
            kept += 1

    eprint(f"[INFO] FORCE bed->bed: kept {kept:,}, dropped {dropped:,} not in .fai, bad {bad:,}")

# -------------------- bedtools wrappers ---------------------------

def bedtools_intersect_count_sorted(a_sorted: Path, b_sorted: Path, out_tsv: Path, genome_fai: Path) -> None:
    run([
        "bedtools", "intersect",
        "-sorted", "-g", str(genome_fai),
        "-a", str(a_sorted),
        "-b", str(b_sorted),
        "-c"
    ], stdout_path=out_tsv)

# -------------------- threshold summary ---------------------------

def safe_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")

def summarize_threshold(interval_tsv: Path, threshold: float, out_summary: Path) -> None:
    """
    interval_tsv columns (no header in bedtools outputs we create later):
      chr start end force length n_any n_ag
    We'll write a 2-row summary: FORCE_GE_THRESHOLD and FORCE_LT_THRESHOLD
    """
    hi_len = hi_any = hi_ag = 0.0
    lo_len = lo_any = lo_ag = 0.0

    with interval_tsv.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        # Expect our header; if not, still try by position.
        # We'll use by name if present.
        col = {c: i for i, c in enumerate(header)}
        # required: force,length,n_edited_any,n_edited_AG
        i_force = col.get("dsRNA_force", 3)
        i_len   = col.get("length", 4)
        i_any   = col.get("n_edited_any", 5)
        i_ag    = col.get("n_edited_AG", 6)

        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            forcev = safe_float(p[i_force]) if i_force < len(p) else float("nan")
            length = safe_float(p[i_len]) if i_len < len(p) else float("nan")
            anyc   = safe_float(p[i_any]) if i_any < len(p) else 0.0
            agc    = safe_float(p[i_ag]) if i_ag < len(p) else 0.0
            if not (length and length > 0):
                continue

            if (not math.isnan(forcev)) and (forcev >= threshold):
                hi_len += length
                hi_any += anyc
                hi_ag  += agc
            else:
                lo_len += length
                lo_any += anyc
                lo_ag  += agc

    def per_kb(ct, bp):
        kb = bp / 1000.0
        return (ct / kb) if kb > 0 else float("nan")

    with out_summary.open("w") as w:
        w.write("\t".join([
            "bin", "threshold",
            "total_bp", "total_kb",
            "any", "ag",
            "any_per_kb", "ag_per_kb"
        ]) + "\n")

        w.write("\t".join(map(str, [
            "FORCE_GE_THRESHOLD", threshold,
            hi_len, hi_len/1000.0,
            hi_any, hi_ag,
            per_kb(hi_any, hi_len), per_kb(hi_ag, hi_len)
        ])) + "\n")

        w.write("\t".join(map(str, [
            "FORCE_LT_THRESHOLD", threshold,
            lo_len, lo_len/1000.0,
            lo_any, lo_ag,
            per_kb(lo_any, lo_len), per_kb(lo_ag, lo_len)
        ])) + "\n")

        # ratio row (hi/lo)
        r_any = per_kb(hi_any, hi_len) / per_kb(lo_any, lo_len) if per_kb(lo_any, lo_len) and not math.isnan(per_kb(lo_any, lo_len)) else float("nan")
        r_ag  = per_kb(hi_ag, hi_len)  / per_kb(lo_ag, lo_len)  if per_kb(lo_ag, lo_len)  and not math.isnan(per_kb(lo_ag, lo_len))  else float("nan")

        w.write("\t".join(map(str, [
            "RATIO_HI_OVER_LO", threshold,
            "NA", "NA",
            "NA", "NA",
            r_any, r_ag
        ])) + "\n")

# ------------------------------- main ----------------------------

def main():
    which_or_die("bedtools")
    which_or_die("awk")
    which_or_die("sort")

    ap = argparse.ArgumentParser(description="Count REDItools edits per dsRNA-force interval BED and summarize force>=threshold vs <threshold.")
    ap.add_argument("--reditools_bed", required=True)
    ap.add_argument("--force_bed", required=True, help="Headerless BED. Col13 (1-based) must be dsRNA_force.")
    ap.add_argument("--genome_fai", required=True)
    ap.add_argument("--out_tsv", required=True, help="Per-interval output TSV.")
    ap.add_argument("--out_summary_tsv", required=True, help="2-bin summary TSV (>=thr vs <thr + ratio).")
    ap.add_argument("--tmpdir", default=None)
    ap.add_argument("--force_col", type=int, default=13, help="1-based column index of dsRNA_force in force BED (default 13).")
    ap.add_argument("--threshold", type=float, default=0.15)
    args = ap.parse_args()

    reditools_path = Path(args.reditools_bed)
    force_bed_in   = Path(args.force_bed)
    genome_fai     = Path(args.genome_fai)
    out_tsv        = Path(args.out_tsv)
    out_sum        = Path(args.out_summary_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    out_sum.parent.mkdir(parents=True, exist_ok=True)

    contigs_set, contig_order = load_fai_contigs_and_rank(genome_fai)

    eprint(f"[INFO] REDItools: {reditools_path}")
    eprint(f"[INFO] Force BED: {force_bed_in} (force col={args.force_col})")
    eprint(f"[INFO] Threshold: {args.threshold}")
    eprint(f"[INFO] Genome fai: {genome_fai}")
    eprint(f"[INFO] Out intervals: {out_tsv}")
    eprint(f"[INFO] Out summary:   {out_sum}")

    tmp_parent = args.tmpdir if args.tmpdir else None
    with tempfile.TemporaryDirectory(dir=tmp_parent) as td_str:
        td = Path(td_str)

        rank_tsv = td / "contig.rank.tsv"
        write_contig_rank_tsv(contig_order, rank_tsv)

        # force BED -> filtered minimal bed (chr start end force length)
        force_filt = td / "force.filtered.bed"
        parse_force_bed_col13_to_bed_filtered(force_bed_in, force_filt, contigs_set, force_col_1based=args.force_col)

        # REDItools sites -> ANY + AG
        sites_any = td / "sites.any.bed"
        sites_ag  = td / "sites.ag.bed"
        eprint("[INFO] Extracting REDItools sites -> BED3 (ANY + AG), contig-filtered to fai")
        write_sites_beds_filtered_any_ag(reditools_path, contigs_set, sites_any, sites_ag)

        # sort all
        force_sorted = td / "force.sorted.bed"
        any_sorted   = td / "sites.any.sorted.bed"
        ag_sorted    = td / "sites.ag.sorted.bed"

        eprint("[INFO] Sorting force/sites beds by genome .fai order (GNU sort + rank map)")
        sort_bed_by_fai_rank(force_filt, force_sorted, rank_tsv)
        sort_bed_by_fai_rank(sites_any, any_sorted, rank_tsv)
        sort_bed_by_fai_rank(sites_ag,  ag_sorted,  rank_tsv)

        # counts per interval
        c_any = td / "force.count_any.tsv"
        c_ag  = td / "force.count_ag.tsv"

        eprint("[INFO] bedtools intersect -c (edited_any) [sorted-streaming]")
        bedtools_intersect_count_sorted(force_sorted, any_sorted, c_any, genome_fai)

        eprint("[INFO] bedtools intersect -c (edited_AG) [sorted-streaming]")
        bedtools_intersect_count_sorted(force_sorted, ag_sorted, c_ag, genome_fai)

        # merge into interval output with header
        # force_sorted has: chr start end force length
        # c_any/c_ag add a last column = count
        with force_sorted.open() as ff, c_any.open() as fa, c_ag.open() as fg, out_tsv.open("w") as w:
            w.write("\t".join(["chr","start","end","dsRNA_force","length","n_edited_any","n_edited_AG"]) + "\n")
            for lf, la, lg in zip(ff, fa, fg):
                pf = lf.rstrip("\n").split("\t")
                pa = la.rstrip("\n").split("\t")
                pg = lg.rstrip("\n").split("\t")
                # sanity: pf == pa[:-1] etc, but don't assume perfect; use force + lastcols
                w.write("\t".join(pf[0:5] + [pa[-1], pg[-1]]) + "\n")

        # threshold summary
        summarize_threshold(out_tsv, args.threshold, out_sum)

    eprint("[OK] Done.")

if __name__ == "__main__":
    main()
