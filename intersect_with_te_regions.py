#!/usr/bin/env python3
"""
Intersect BED files with TE annotation files to identify TE classes (SINE, LINE, LTR, nonTE)
This script only adds TE class annotation - it does not process dsRNA force data
"""

import argparse
import sys
from pathlib import Path
import subprocess
import tempfile
import os

def sort_bed_file(bed_file, sorted_file):
    """Sort a BED file using system sort command"""
    try:
        subprocess.run(
            ['sort', '-k1,1', '-k2,2n', bed_file],
            stdout=open(sorted_file, 'w'),
            stderr=subprocess.PIPE,
            check=True
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to sort {bed_file}: {e.stderr.decode()}", file=sys.stderr)
        return False

def bedtools_intersect(query_bed, target_bed, output_file):
    """Run bedtools intersect and return count of intersections"""
    try:
        result = subprocess.run(
            ['bedtools', 'intersect', '-a', query_bed, '-b', target_bed, '-wa'],
            stdout=open(output_file, 'w'),
            stderr=subprocess.PIPE,
            check=True
        )
        # Count lines in output
        count = sum(1 for _ in open(output_file)) if os.path.exists(output_file) else 0
        return count
    except subprocess.CalledProcessError as e:
        print(f"[WARN] bedtools intersect failed: {e.stderr.decode()}", file=sys.stderr)
        # Create empty file
        Path(output_file).touch()
        return 0

def read_bed_file(bed_file):
    """Read BED file and return as list of lines (first 10 columns)"""
    sites = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            # Create key from chr, start, end
            key = (fields[0], fields[1], fields[2])
            # Store full line (first 10 columns, pad if needed)
            site_line = '\t'.join(fields[:10]) if len(fields) >= 10 else '\t'.join(fields) + '\t' * (10 - len(fields))
            sites[key] = site_line
    return sites

def process_file(bed_file, sine_bed, line_bed, ltr_bed, out_dir):
    """Process a single BED file and annotate with TE classes"""
    basename = Path(bed_file).stem
    print(f"[INFO] Processing {basename}...")
    
    # Check input file
    input_count = sum(1 for _ in open(bed_file)) if os.path.exists(bed_file) else 0
    print(f"  Input BED file has {input_count} regions")
    
    # Sort input BED file
    sorted_bed = out_dir / f"{basename}.sorted.bed"
    if not sort_bed_file(bed_file, sorted_bed):
        print(f"[ERROR] Failed to sort {bed_file}")
        return False
    
    sorted_count = sum(1 for _ in open(sorted_bed)) if os.path.exists(sorted_bed) else 0
    print(f"  Sorted BED file has {sorted_count} regions")
    
    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory(dir=out_dir, prefix=f"{basename}_") as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Intersect with SINE
        sine_intersect = tmpdir / "sine_intersect.bed"
        print("  Intersecting with SINE.bed...")
        sine_count = bedtools_intersect(sorted_bed, sine_bed, sine_intersect)
        print(f"  Found {sine_count} SINE intersections")
        
        # Intersect with LINE
        line_intersect = tmpdir / "line_intersect.bed"
        print("  Intersecting with LINE.bed...")
        line_count = bedtools_intersect(sorted_bed, line_bed, line_intersect)
        print(f"  Found {line_count} LINE intersections")
        
        # Intersect with LTR
        ltr_intersect = tmpdir / "ltr_intersect.bed"
        print("  Intersecting with LTR.bed...")
        ltr_count = bedtools_intersect(sorted_bed, ltr_bed, ltr_intersect)
        print(f"  Found {ltr_count} LTR intersections")
        
        # Combine with priority (SINE > LINE > LTR)
        print("  Combining intersections with priority (SINE > LINE > LTR)...")
        
        # Read all intersection files
        sine_sites = read_bed_file(sine_intersect) if sine_intersect.exists() else {}
        line_sites = read_bed_file(line_intersect) if line_intersect.exists() else {}
        ltr_sites = read_bed_file(ltr_intersect) if ltr_intersect.exists() else {}
        
        # Assign classes with priority
        site_classes = {}
        site_lines = {}
        
        # SINE has highest priority
        for key, line in sine_sites.items():
            site_classes[key] = "SINE"
            site_lines[key] = line
        
        # LINE (only if not already SINE)
        for key, line in line_sites.items():
            if key not in site_classes:
                site_classes[key] = "LINE"
                site_lines[key] = line
        
        # LTR (only if not already SINE or LINE)
        for key, line in ltr_sites.items():
            if key not in site_classes:
                site_classes[key] = "LTR"
                site_lines[key] = line
        
        # Read all input sites to find non-TE
        all_sites = read_bed_file(sorted_bed)
        
        # Write output
        output_file = out_dir / f"{basename}.te_annotated.bed"
        with open(output_file, 'w') as f:
            # Write TE sites
            for key in sorted(site_classes.keys()):
                f.write(f"{site_lines[key]}\t{site_classes[key]}\n")
            
            # Write non-TE sites
            for key in sorted(all_sites.keys()):
                if key not in site_classes:
                    f.write(f"{all_sites[key]}\tnonTE\n")
        
        # Show breakdown
        total_count = sum(1 for _ in open(output_file)) if os.path.exists(output_file) else 0
        print(f"  Combined TE-annotated file: {output_file}")
        print(f"  Total annotated regions: {total_count}")
        
        if os.path.exists(output_file):
            # Count by class
            sine_final = sum(1 for line in open(output_file) if line.strip().endswith('\tSINE'))
            line_final = sum(1 for line in open(output_file) if line.strip().endswith('\tLINE'))
            ltr_final = sum(1 for line in open(output_file) if line.strip().endswith('\tLTR'))
            nonte_final = sum(1 for line in open(output_file) if line.strip().endswith('\tnonTE'))
            
            print("  TE class breakdown:")
            print(f"    SINE:   {sine_final}")
            print(f"    LINE:   {line_final}")
            print(f"    LTR:    {ltr_final}")
            print(f"    nonTE:  {nonte_final}")
        
        print(f"  [OK] Successfully processed {basename}")
        return True

def main():
    parser = argparse.ArgumentParser(
        description="Intersect BED files with TE annotation files (SINE, LINE, LTR)"
    )
    parser.add_argument("--in_dir", required=True, help="Input directory with BED files")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--te_dir", required=True, help="Directory with TE annotation BED files")
    parser.add_argument("--pattern", default="*.AGTC.cov25.freq0.01.bed", help="File pattern to match")
    parser.add_argument("--file_index", type=int, help="Array job file index (0-based)")
    
    args = parser.parse_args()
    
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    te_dir = Path(args.te_dir)
    
    # Create output directory
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Check TE annotation files
    sine_bed = te_dir / "mm10_SINE.bed.sorted"
    line_bed = te_dir / "mm10_LINE.bed.sorted"
    ltr_bed = te_dir / "mm10_LTR.bed.sorted"
    
    # Sort TE files if needed
    for te_file, sorted_file in [
        (te_dir / "mm10_SINE.bed", sine_bed),
        (te_dir / "mm10_LINE.bed", line_bed),
        (te_dir / "mm10_LTR.bed", ltr_bed),
    ]:
        if not te_file.exists():
            print(f"[ERROR] TE annotation file not found: {te_file}", file=sys.stderr)
            sys.exit(1)
        
        if not sorted_file.exists() or te_file.stat().st_mtime > sorted_file.stat().st_mtime:
            print(f"[INFO] Sorting {te_file.name}...")
            if not sort_bed_file(te_file, sorted_file):
                print(f"[ERROR] Failed to sort {te_file}", file=sys.stderr)
                sys.exit(1)
    
    print("[INFO] TE annotation files ready")
    
    # Find BED files
    bed_files = sorted(in_dir.glob(args.pattern))
    total_files = len(bed_files)
    
    print(f"[INFO] Found {total_files} BED files to process")
    
    # If file_index is provided (array job), process only that file
    if args.file_index is not None:
        if args.file_index < 0 or args.file_index >= total_files:
            print(f"[INFO] File index {args.file_index} is out of range (0-{total_files-1}), skipping")
            sys.exit(0)
        bed_files = [bed_files[args.file_index]]
        print(f"[INFO] Array job mode: Processing file {args.file_index + 1} of {total_files}")
    else:
        print(f"[INFO] Non-array mode: Processing all {total_files} files")
    
    # Process each file
    success_count = 0
    for bed_file in bed_files:
        if process_file(bed_file, sine_bed, line_bed, ltr_bed, out_dir):
            success_count += 1
        print()  # Blank line between files
    
    print(f"[INFO] TE intersection complete. Processed {success_count}/{len(bed_files)} files successfully.")
    print(f"[INFO] Output directory: {out_dir}")

if __name__ == "__main__":
    main()



