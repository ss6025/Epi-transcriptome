#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd

def eprint(*a):
    print(*a, file=sys.stderr)

def is_alu_family(x: str) -> bool:
    s = str(x).strip().lower()
    # force TSV has lowercase family/class columns per your header
    # Accept common Alu labels
    return ("alu" in s)

def has_ir_signal(v) -> bool:
    """
    Use dsRNA_seqA column as the IR indicator.
    Treat non-empty / not '.' / not 'NA' as IR-present.
    """
    if pd.isna(v):
        return False
    s = str(v).strip()
    if s == "" or s == ".":
        return False
    if s.lower() in ("na", "nan", "none"):
        return False
    return True

def write_bed(df, out_path, cols=("chr","start","end")):
    df = df.loc[:, list(cols)].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"], errors="coerce")
    df = df.dropna(subset=["start","end"])
    df["start"] = df["start"].astype(int)
    df["end"]   = df["end"].astype(int)

    # drop invalid intervals
    df = df[(df["start"] >= 0) & (df["end"] > df["start"])]

    # sort
    df = df.sort_values(["chr","start","end"])
    df.to_csv(out_path, sep="\t", header=False, index=False)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force_tsv", required=True,
                    help="squire_seqA_intersect_100%_annot.tsv with columns incl chr,start,end,dsRNA_seqA,class,family,max_mean_dsRNA_force")
    ap.add_argument("--rmsk_bed", required=True,
                    help="mm10_rmsk.bed (BED-ish; may have many columns)")
    ap.add_argument("--out_dir", required=True,
                    help="Output dir for MOR.bed, IR_Alu.bed, nonIR_Alu.bed")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # -----------------------
    # (A) IR_Alu / nonIR_Alu from force TSV
    # -----------------------
    eprint(f"[INFO] Reading force TSV: {args.force_tsv}")
    f = pd.read_csv(args.force_tsv, sep="\t", low_memory=False)

    required = {"chr","start","end","dsRNA_seqA","family"}
    missing = required - set(f.columns)
    if missing:
        raise SystemExit(f"[FATAL] force_tsv missing columns: {sorted(list(missing))}\n"
                         f"Found: {list(f.columns)}")

    # Keep only Alu-family rows
    f_alu = f[f["family"].apply(is_alu_family)].copy()
    eprint(f"[INFO] Force TSV Alu rows: {len(f_alu):,}")

    # IR vs nonIR using dsRNA_seqA presence
    f_alu["is_ir"] = f_alu["dsRNA_seqA"].apply(has_ir_signal)

    ir_bed = os.path.join(args.out_dir, "IR_Alu.bed")
    nonir_bed = os.path.join(args.out_dir, "nonIR_Alu.bed")

    write_bed(f_alu[f_alu["is_ir"]], ir_bed)
    write_bed(f_alu[~f_alu["is_ir"]], nonir_bed)

    eprint(f"[OK] Wrote {ir_bed}")
    eprint(f"[OK] Wrote {nonir_bed}")

    # -----------------------
    # (B) MOR.bed from mm10_rmsk.bed
    # -----------------------
    eprint(f"[INFO] Building MOR.bed from: {args.rmsk_bed}")

    # mm10_rmsk.bed is typically BED-like: chr start end name score strand class family ... (varies)
    # We will scan each line and keep if any field contains "MOR" (case-insensitive).
    mor_out = os.path.join(args.out_dir, "MOR.bed")
    kept = 0
    total = 0

    with open(args.rmsk_bed) as r, open(mor_out, "w") as w:
        for line in r:
            if not line.strip() or line.startswith("#"):
                continue
            total += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            hay = "\t".join(parts).lower()
            if "mor" not in hay:
                continue
            # write first 3 bed columns
            try:
                chrom = parts[0]
                start = int(float(parts[1]))
                end   = int(float(parts[2]))
            except Exception:
                continue
            if end > start >= 0:
                w.write(f"{chrom}\t{start}\t{end}\n")
                kept += 1

    eprint(f"[OK] Wrote {mor_out} (kept {kept:,} / scanned {total:,})")

    # tiny report
    report = os.path.join(args.out_dir, "beds_build_report.txt")
    with open(report, "w") as rep:
        rep.write(f"force_tsv\t{args.force_tsv}\n")
        rep.write(f"rmsk_bed\t{args.rmsk_bed}\n")
        rep.write(f"IR_Alu.bed\t{ir_bed}\n")
        rep.write(f"nonIR_Alu.bed\t{nonir_bed}\n")
        rep.write(f"MOR.bed\t{mor_out}\n")
    eprint(f"[OK] Wrote {report}")

if __name__ == "__main__":
    main()
