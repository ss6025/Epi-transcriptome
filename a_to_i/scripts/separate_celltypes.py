#!/usr/bin/env python3
"""
Separate BAM files by cell type based on CB:Z barcode tags.

1. Reads markers TSV (combined_markers_noNA.tsv) with barcode and mannualAnnot columns
2. Extracts barcodes from BAM files using samtools view
3. Separates reads by cell type into separate BAM files
4. Output files are named {parent}_{celltype}.bam
"""

import argparse
import subprocess
from pathlib import Path
import pandas as pd
from collections import defaultdict
import re


def normalize_barcode(barcode: str) -> str:
    """
    Normalize barcode for matching.
    BAM barcodes: CB:Z:CTTACCGGTCCAGAAG-1
    CSV barcodes: AAACCCACAGTCTACA-1_1
    Match on the base sequence (ACGT only, before any suffix).
    """
    if not barcode:
        return ""
    barcode = barcode.replace("CB:Z:", "")
    match = re.match(r"^([ACGT]+)", barcode)
    if match:
        return match.group(1)
    return barcode


def find_column(df: pd.DataFrame, candidates: list) -> str:
    cols_lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


def load_markers_tsv(markers_file: Path) -> dict:
    print(f"[INFO] Loading markers from: {markers_file}")
    df = pd.read_csv(markers_file, sep='\t')
    if len(df.columns) == 1 and ',' in str(df.columns[0]):
        df = pd.read_csv(markers_file, sep=',')

    barcode_col  = find_column(df, ['barcode', 'Barcode'])
    celltype_col = find_column(df, ['mannualAnnot', 'manualAnnot', 'ManualAnnot'])

    if not barcode_col:
        raise SystemExit(f"[FATAL] Missing barcode column. Available: {list(df.columns)}")
    if not celltype_col:
        raise SystemExit(f"[FATAL] Missing mannualAnnot/manualAnnot column. Available: {list(df.columns)}")

    barcode_map = {}
    for _, row in df.iterrows():
        barcode_full = str(row[barcode_col]).strip()
        cell_type    = str(row[celltype_col]).strip() or 'Unknown'
        normalized   = normalize_barcode(barcode_full)
        if normalized:
            barcode_map[normalized] = {'celltype': cell_type, 'barcode_full': barcode_full}

    print(f"[INFO] Loaded {len(barcode_map)} barcode mappings from {celltype_col}")
    print(f"[INFO] Unique cell types: {df[celltype_col].nunique()}")
    return barcode_map


def extract_barcodes_from_bam(bam_file: Path) -> dict:
    print(f"[INFO] Extracting barcodes from: {bam_file}")
    cmd = ['samtools', 'view', str(bam_file)]
    try:
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        read_barcodes = {}
        total_reads   = 0
        reads_with_cb = 0

        for line in process.stdout:
            total_reads += 1
            if total_reads % 1_000_000 == 0:
                print(f"  Processed {total_reads:,} reads...", end='\r')
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            read_name = fields[0]
            barcode   = None
            for field in fields[11:]:
                if field.startswith('CB:Z:'):
                    barcode = field[5:]
                    break
            if barcode:
                reads_with_cb += 1
                normalized = normalize_barcode(barcode)
                if normalized:
                    read_barcodes[read_name] = normalized

        process.wait()
        if process.returncode != 0:
            stderr = process.stderr.read()
            raise subprocess.CalledProcessError(process.returncode, cmd, stderr=stderr)

        print(f"\n[INFO] Total reads: {total_reads:,}")
        print(f"[INFO] Reads with CB tag: {reads_with_cb:,} ({100*reads_with_cb/max(total_reads,1):.1f}%)")
        print(f"[INFO] Unique barcodes found: {len(set(read_barcodes.values())):,}")
        return read_barcodes

    except subprocess.CalledProcessError as e:
        raise SystemExit(f"[FATAL] Error running samtools: {e.stderr}")
    except FileNotFoundError:
        raise SystemExit("[FATAL] samtools not found.")


def separate_bam_by_celltype(bam_file: Path, markers_file: Path, output_dir: Path):
    barcode_map   = load_markers_tsv(markers_file)
    read_barcodes = extract_barcodes_from_bam(bam_file)

    celltype_reads = defaultdict(set)
    unknown_reads  = set()

    print("[INFO] Grouping reads by cell type...")
    for read_name, nb in read_barcodes.items():
        if nb in barcode_map:
            cell_type_safe = re.sub(r'[^\w\-]', '_', barcode_map[nb]['celltype'])
            celltype_reads[cell_type_safe].add(read_name)
        else:
            unknown_reads.add(read_name)

    print(f"[INFO] Found {len(celltype_reads)} cell types")
    print(f"[INFO] Reads with unknown barcodes: {len(unknown_reads):,}")

    parent_name = bam_file.stem
    parent_name = parent_name.replace('possorted_genome_', '').replace('_bam', '')

    output_dir.mkdir(parents=True, exist_ok=True)
    read_list_dir = output_dir / "read_lists"
    read_list_dir.mkdir(exist_ok=True)

    celltype_files = {}
    for cell_type, read_names in celltype_reads.items():
        read_list_file = read_list_dir / f"{parent_name}_{cell_type}_reads.txt"
        with open(read_list_file, 'w') as f:
            for rn in read_names:
                f.write(f"{rn}\n")
        output_bam = output_dir / f"{parent_name}_{cell_type}.bam"
        celltype_files[cell_type] = {
            'read_list': read_list_file,
            'output_bam': output_bam,
            'read_count': len(read_names),
        }

    print("[INFO] Extracting reads to separate BAM files...")
    for cell_type, info in celltype_files.items():
        print(f"  Processing {cell_type} ({info['read_count']:,} reads)...")
        cmd = ['samtools', 'view', '-b', '-N', str(info['read_list']), str(bam_file)]
        try:
            with open(info['output_bam'], 'wb') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, check=True)
            subprocess.run(['samtools', 'index', str(info['output_bam'])], check=True, stderr=subprocess.PIPE)
            print(f"    Created: {info['output_bam']}")
        except subprocess.CalledProcessError as e:
            print(f"    [ERROR] Failed for {cell_type}: {e.stderr.decode()}")

    summary_file = output_dir / f"{parent_name}_separation_summary.txt"
    with open(summary_file, 'w') as f:
        f.write(f"Parent BAM: {bam_file}\nMarkers: {markers_file}\n\n")
        f.write(f"{'Cell Type':<50} {'Read Count':>15} {'Output BAM':<50}\n")
        f.write(f"{'-'*80}\n")
        for ct, info in sorted(celltype_files.items(), key=lambda x: -x[1]['read_count']):
            f.write(f"{ct:<50} {info['read_count']:>15,} {info['output_bam'].name:<50}\n")
        f.write(f"\nUnknown barcodes: {len(unknown_reads):,}\n")
        f.write(f"Total reads: {len(read_barcodes):,}\n")

    print(f"\n[OK] Done — {len(celltype_files)} cell-type BAMs in {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Separate BAM by cell type via CB:Z tags")
    parser.add_argument('--bam',        required=True,  type=Path)
    parser.add_argument('--markers',    type=Path, help='TSV with barcode + mannualAnnot columns')
    parser.add_argument('--metadata',   type=Path, help='Alias for --markers (backward compat)')
    parser.add_argument('--output_dir', required=True,  type=Path)
    args = parser.parse_args()

    markers_file = args.markers or args.metadata
    if not markers_file:
        raise SystemExit("[FATAL] Must provide --markers or --metadata")
    if not args.bam.exists():
        raise SystemExit(f"[FATAL] BAM not found: {args.bam}")
    if not markers_file.exists():
        raise SystemExit(f"[FATAL] Markers file not found: {markers_file}")

    try:
        subprocess.run(['samtools', '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise SystemExit("[FATAL] samtools not found.")

    separate_bam_by_celltype(args.bam, markers_file, args.output_dir)


if __name__ == "__main__":
    main()
