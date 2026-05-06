#!/usr/bin/env python3
"""
Run REDItools on all BAM files in a celltypes/ directory.

Processes BAMs from:
  {CELLTYPE_BAMS}/{sample}/celltypes/

Excludes: *_nan.bam by default.
Output: {parent}_{celltype}_RNA_redi.bed, plus a combined all_celltypes BED.
"""
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import re


def find_reditools(reditools_path: Path = None) -> str | None:
    if reditools_path and reditools_path.exists():
        return str(reditools_path)
    default = Path("/data1/greenbab/users/rotembc1/pipeline2.0/REDItools2/src/cineca/reditools.py")
    if default.exists():
        return str(default)
    for cmd in ['reditools.py', 'REDItools2.py']:
        try:
            r = subprocess.run(['which', cmd], capture_output=True, text=True, check=False)
            if r.returncode == 0:
                return r.stdout.strip()
        except Exception:
            pass
    return None


def build_ref_matching_bam(bam_file: Path, ref_file: Path, work_dir: Path,
                            samtools: str = "samtools") -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    ref_fai = Path(str(ref_file) + '.fai')
    if not ref_fai.exists():
        subprocess.run([samtools, 'faidx', str(ref_file)], check=True)

    ref_contigs_file = work_dir / "ref.contigs.txt"
    with open(ref_fai) as f:
        contigs = [line.split('\t')[0] for line in f if line.strip()]
    with open(ref_contigs_file, 'w') as f:
        for c in sorted(set(contigs)):
            f.write(f"{c}\n")

    bam_contigs_file = work_dir / "bam.contigs.txt"
    r = subprocess.run([samtools, 'idxstats', str(bam_file)], capture_output=True, text=True, check=True)
    bam_contigs = [l.split('\t')[0] for l in r.stdout.split('\n') if l.strip() and l.split('\t')[0] != '*']
    with open(bam_contigs_file, 'w') as f:
        for c in sorted(set(bam_contigs)):
            f.write(f"{c}\n")

    keep_contigs_file = work_dir / "keep.contigs.txt"
    r = subprocess.run(['comm', '-12', str(ref_contigs_file), str(bam_contigs_file)],
                       capture_output=True, text=True, check=True)
    common = [l.strip() for l in r.stdout.split('\n') if l.strip()]
    with open(keep_contigs_file, 'w') as f:
        for c in common:
            f.write(f"{c}\n")

    if len(common) < 5:
        print("[WARN] Too few shared contigs — using original BAM")
        return bam_file

    refmatch_bam = work_dir / f"{bam_file.stem}.refmatch.bam"
    awk_script   = work_dir / "filter_contigs.awk"
    with open(awk_script, 'w') as f:
        f.write(f'''BEGIN {{
    while((getline < "{keep_contigs_file}") > 0) {{ ok[$1] = 1 }}
}}
/^@/ {{
    if($1 == "@SQ") {{
        sn = ""
        for(i=1; i<=NF; i++) {{ if($i ~ /^SN:/) {{ sn = substr($i, 4); break }} }}
        if(sn == "" || ok[sn]) print
    }} else print
    next
}}
{{ r=$3; rn=$7; if(!ok[r]) next; if(rn!="=" && rn!="*" && !ok[rn]) next; print }}''')

    cmd = (f'{samtools} view -h "{bam_file}" | awk -f "{awk_script}" '
           f'| {samtools} view -b - | {samtools} sort -@ 4 -o "{refmatch_bam}" -')
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"[WARN] Failed to build REF-matching BAM: {r.stderr}")
        return bam_file
    subprocess.run([samtools, 'index', str(refmatch_bam)], check=True)
    return refmatch_bam


def run_reditools(bam_file: Path, ref_file: Path, output_file: Path,
                  reditools_path: str, threads: int = 8,
                  build_refmatch: bool = True, work_dir: Path = None) -> bool:
    if build_refmatch and work_dir:
        bam_to_use = build_ref_matching_bam(bam_file, ref_file, work_dir)
    else:
        bam_to_use = bam_file

    py_cmd = 'python' if ('cineca' in reditools_path or 'REDItools2' in reditools_path) else 'python3'
    cmd = [py_cmd, reditools_path, '-f', str(bam_to_use), '-r', str(ref_file),
           '-o', str(output_file), '-t', str(threads)]
    print(f"[INFO] {' '.join(cmd)}")

    stderr_file = output_file.parent / f"{output_file.stem}.reditools.stderr.log"
    with open(stderr_file, 'w') as sf:
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=sf, text=True)

    if r.returncode != 0:
        print(f"[ERROR] REDItools failed (exit {r.returncode})")
        return False
    if not output_file.exists() or output_file.stat().st_size == 0:
        print("[WARN] Empty output from REDItools")
        return False
    return True


def output_prefix_and_celltype(bam_file: Path, sample_folder: str) -> tuple[str, str]:
    stem = bam_file.stem
    if sample_folder and stem.startswith(sample_folder + "_"):
        return sample_folder, stem[len(sample_folder) + 1:]
    parts = stem.rsplit('_', 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (stem, "unknown")


def list_bams(celltypes_dir: Path, exclude_nan: bool = True,
              bam_regex: str | None = None) -> list[Path]:
    bams = sorted(f for f in celltypes_dir.glob("*.bam") if not f.name.endswith(".bai"))
    if exclude_nan:
        bams = [f for f in bams if not f.stem.endswith("_nan")]
    if bam_regex:
        cre = re.compile(bam_regex)
        bams = [f for f in bams if cre.search(f.name)]
    return bams


def add_celltype_column(output_file: Path, celltype: str) -> bool:
    import gzip
    is_gz = output_file.suffix == '.gz'
    tmp   = output_file.parent / f"{output_file.stem}.tmp{output_file.suffix}"
    try:
        fin  = gzip.open(output_file, 'rt') if is_gz else open(output_file)
        fout = gzip.open(tmp, 'wt')         if is_gz else open(tmp, 'w')
        hdr  = fin.readline()
        if not hdr or 'CellType' in hdr:
            fin.close(); fout.close(); tmp.unlink(missing_ok=True)
            return True
        fout.write(hdr.rstrip() + "\tCellType\n")
        for line in fin:
            if line.strip():
                fout.write(line.rstrip() + f"\t{celltype}\n")
        fin.close(); fout.close()
        output_file.unlink()
        tmp.rename(output_file)
        return True
    except Exception as e:
        print(f"[ERROR] add_celltype_column: {e}")
        tmp.unlink(missing_ok=True)
        return False


def combine_beds(files: list[Path], combined: Path) -> bool:
    import gzip
    is_gz   = combined.suffix == '.gz'
    written = 0
    try:
        out = gzip.open(combined, 'wt') if is_gz else open(combined, 'w')
        hdr_done = False
        for f in sorted(files):
            if not f.exists():
                continue
            inp = gzip.open(f, 'rt') if f.suffix == '.gz' else open(f)
            hdr = inp.readline()
            if not hdr:
                inp.close(); continue
            if not hdr_done:
                out.write(hdr); hdr_done = True
            for line in inp:
                if line.strip():
                    out.write(line); written += 1
            inp.close()
        out.close()
        print(f"[OK] Combined {written:,} lines -> {combined.name}")
        return True
    except Exception as e:
        print(f"[ERROR] combine_beds: {e}")
        return False


def process_celltypes_dir(celltypes_dir: Path, output_dir: Path, ref_file: Path,
                           reditools_path: str, threads: int = 8,
                           build_refmatch: bool = True, combine: bool = True,
                           exclude_nan: bool = True, bam_regex: str | None = None,
                           combine_suffix: str = "all_celltypes") -> int:
    bams = list_bams(celltypes_dir, exclude_nan=exclude_nan, bam_regex=bam_regex)
    if not bams:
        print(f"[WARN] No BAMs in {celltypes_dir}")
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    sample_name  = celltypes_dir.parent.name
    work_base    = output_dir / "_work"
    success      = 0
    output_files = []

    for bam in bams:
        print(f"\n{'='*70}\n[INFO] {bam.name}\n{'='*70}")
        parent, cell_type = output_prefix_and_celltype(bam, sample_name)
        ct_safe           = re.sub(r'[^\w\-]', '_', cell_type)
        tag               = f"{parent}_{ct_safe}"
        out_bed           = output_dir / f"{tag}_RNA_redi.bed"
        work_dir          = work_base / tag

        if run_reditools(bam, ref_file, out_bed, reditools_path, threads,
                         build_refmatch, work_dir):
            add_celltype_column(out_bed, cell_type)
            output_files.append(out_bed)
            success += 1

    if combine and output_files:
        combined_path = output_dir / f"{sample_name}_{combine_suffix}_RNA_redi.bed"
        combine_beds(output_files, combined_path)

    return success


def main():
    parser = argparse.ArgumentParser(description="Run REDItools on cell-type BAMs")
    parser.add_argument('--base_dir',     type=Path)
    parser.add_argument('--celltypes_dir',type=Path)
    parser.add_argument('--output_dir',   required=True, type=Path)
    parser.add_argument('--ref',          required=True, type=Path)
    parser.add_argument('--reditools_path', type=Path)
    parser.add_argument('--threads',      type=int, default=8)
    parser.add_argument('--no_refmatch',  action='store_true')
    parser.add_argument('--no_combine',   action='store_true')
    parser.add_argument('--include_nan',  action='store_true')
    parser.add_argument('--bam-regex',    type=str, default=None)
    parser.add_argument('--combine-suffix', type=str, default='all_celltypes')
    args = parser.parse_args()

    if not args.base_dir and not args.celltypes_dir:
        raise SystemExit("[FATAL] Must provide --base_dir or --celltypes_dir")
    if args.base_dir and args.celltypes_dir:
        raise SystemExit("[FATAL] Cannot provide both")
    if not args.ref.exists():
        raise SystemExit(f"[FATAL] Reference not found: {args.ref}")

    rt = find_reditools(args.reditools_path)
    if not rt:
        raise SystemExit("[FATAL] REDItools not found")

    kw = dict(ref_file=args.ref, reditools_path=rt, threads=args.threads,
               build_refmatch=not args.no_refmatch, combine=not args.no_combine,
               exclude_nan=not args.include_nan, bam_regex=args.bam_regex,
               combine_suffix=args.combine_suffix)

    if args.celltypes_dir:
        n = process_celltypes_dir(args.celltypes_dir, args.output_dir, **kw)
        print(f"\n[OK] Processed {n} BAM file(s)")
    else:
        dirs = sorted(d for d in args.base_dir.glob("*/celltypes") if d.is_dir())
        if not dirs:
            raise SystemExit(f"[FATAL] No */celltypes dirs in {args.base_dir}")
        total = 0
        for d in dirs:
            total += process_celltypes_dir(d, args.output_dir, **kw)
        print(f"\n[OK] Processed {total} BAM file(s) across {len(dirs)} samples")


if __name__ == "__main__":
    main()
