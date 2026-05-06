#!/usr/bin/env python3
"""
QC for filtered AGTC editing sites (step 4 output: *.AGTC.bed).

Plots per sample:
  - Frequency distribution (all sites)
  - Coverage distribution (AG and TC separately)
  - Coverage distribution per cell type (subplot grid)
  - Summary bar charts: mean freq and mean coverage per cell type

Usage:
  python3 qc_agtc.py --input_dir filtered_agtc/ --out_dir qc_agtc/

  # One sample (for parallel sbatch):
  python3 qc_agtc.py --input_dir filtered_agtc/ --out_dir qc_agtc/ --sample MySample
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd


CHUNKSIZE = 500_000


# ── column name helpers ──────────────────────────────────────────────────────

def find_col(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lc = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lc:
            return lc[cand.lower()]
    return None


# ── histogram helpers ────────────────────────────────────────────────────────

def weighted_quantile(xs: np.ndarray, ys: np.ndarray, q: float) -> float | None:
    if xs.size == 0:
        return None
    cum = np.cumsum(ys)
    if cum[-1] == 0:
        return None
    idx = int(np.searchsorted(cum, q * int(cum[-1]), side="left"))
    return float(xs[min(idx, len(xs) - 1)])


def hist_arrays(hist: dict) -> tuple[np.ndarray, np.ndarray]:
    xs = np.array(sorted(hist), dtype=float)
    ys = np.array([hist[x] for x in xs], dtype=np.int64)
    return xs, ys


# ── individual plot functions ─────────────────────────────────────────────────

def plot_freq_dist(hist: dict, title: str = "Editing frequency distribution") -> plt.Figure | None:
    if not hist:
        return None
    xs, ys = hist_arrays(hist)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.step(xs, ys, where="mid", linewidth=1.6, color="steelblue")
    med = weighted_quantile(xs, ys, 0.5)
    p99 = weighted_quantile(xs, ys, 0.99)
    if med is not None:
        ax.axvline(med, linestyle="--", linewidth=1.2, color="coral",  label=f"median = {med:.3f}")
    if p99 is not None:
        ax.axvline(p99, linestyle=":",  linewidth=1.2, color="purple", label=f"p99    = {p99:.3f}")
    ax.set_xlabel("Frequency")
    ax.set_ylabel("Sites")
    ax.set_title(title)
    ax.set_xlim(0, 1)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_cov_dist(hist: dict, sub: str = "", logy: bool = True,
                  x_max_q: float = 0.999) -> plt.Figure | None:
    if not hist:
        return None
    xs, ys = hist_arrays(hist)
    x_max = weighted_quantile(xs, ys, x_max_q)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.step(xs, ys, where="mid", linewidth=1.6, color="seagreen")
    med = weighted_quantile(xs, ys, 0.5)
    if med is not None:
        ax.axvline(med, linestyle="--", linewidth=1.2, color="coral", label=f"median = {med:.0f}")
    ax.set_xlabel("Coverage")
    ax.set_ylabel("Sites")
    ax.set_title(f"Coverage distribution{' — ' + sub if sub else ''}")
    if x_max:
        ax.set_xlim(0, x_max)
    if logy:
        ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_cov_per_celltype(hists: dict[str, dict], logy: bool = True,
                           x_max_q: float = 0.999) -> plt.Figure | None:
    """Subplot grid: one coverage histogram per cell type."""
    celltypes = sorted(hists)
    if not celltypes:
        return None
    n = len(celltypes)
    n_cols = min(4, n)
    n_rows = (n + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows), squeeze=False)

    # global x_max across all cell types
    all_x: list[float] = []
    for ct in celltypes:
        xs, ys = hist_arrays(hists[ct])
        q = weighted_quantile(xs, ys, x_max_q)
        if q is not None:
            all_x.append(q)
    x_max = max(all_x) if all_x else 100.0
    x_max = max(x_max, 10)

    for i, ct in enumerate(celltypes):
        ax = axes.flat[i]
        xs, ys = hist_arrays(hists[ct])
        if xs.size > 0:
            ax.step(xs, ys, where="mid", linewidth=1.2, color="steelblue")
        ax.set_xlabel("Coverage")
        ax.set_ylabel("Sites")
        ax.set_title(ct, fontsize=9)
        ax.set_xlim(0, x_max)
        if logy:
            ax.set_yscale("log")
        ax.grid(alpha=0.25)

    for j in range(n, n_rows * n_cols):
        axes.flat[j].set_visible(False)

    fig.suptitle("Coverage distribution per cell type (AGTC sites)", y=1.01, fontsize=11)
    fig.tight_layout()
    return fig


def plot_bar(celltypes: list[str], values: list[float], ylabel: str,
             title: str, color: str = "steelblue") -> plt.Figure:
    fig, ax = plt.subplots(figsize=(max(8, len(celltypes) * 0.55), 5))
    ax.bar(celltypes, values, color=color, edgecolor="black", alpha=0.85)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Cell type")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=45)
    plt.setp(ax.get_xticklabels(), ha="right")
    ax.grid(axis="y", alpha=0.3)
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter("{x:,.3f}"))
    fig.tight_layout()
    return fig


def savefig(fig: plt.Figure, path: Path) -> None:
    for ext in (".pdf", ".png"):
        fig.savefig(path.with_suffix(ext), dpi=150, bbox_inches="tight")
    plt.close(fig)


# ── process one AGTC BED file ─────────────────────────────────────────────────

def process_file(path: Path) -> dict:
    """
    Returns aggregated stats dict:
      freq_hist, cov_ag_hist, cov_tc_hist, cov_per_celltype,
      celltype_stats { ct: {n, freq_sum, freq_n, cov_sum, cov_n} }
    """
    header = pd.read_csv(path, sep="\t", nrows=0)
    freq_col     = find_col(header, ["Frequency", "freq", "frequency"])
    cov_col      = find_col(header, ["Coverage", "Coverage-q30", "coverageq30", "cov"])
    subs_col     = find_col(header, ["AllSubs", "allsubs", "subs", "sub"])
    celltype_col = find_col(header, ["CellType", "celltype", "cell_type", "cell"])

    freq_hist:        dict[float, int]             = defaultdict(int)
    cov_ag_hist:      dict[int,   int]             = defaultdict(int)
    cov_tc_hist:      dict[int,   int]             = defaultdict(int)
    cov_per_celltype: dict[str, dict[int, int]]    = defaultdict(lambda: defaultdict(int))
    ct_stats:         dict[str, dict]              = defaultdict(
        lambda: dict(n=0, freq_sum=0.0, freq_n=0, cov_sum=0.0, cov_n=0)
    )

    for chunk in pd.read_csv(path, sep="\t", chunksize=CHUNKSIZE, low_memory=False):
        if freq_col and freq_col in chunk.columns:
            chunk[freq_col] = pd.to_numeric(chunk[freq_col], errors="coerce")
        if cov_col and cov_col in chunk.columns:
            chunk[cov_col] = pd.to_numeric(chunk[cov_col], errors="coerce")

        for row in chunk.itertuples(index=False):
            ct   = str(getattr(row, celltype_col, "Unknown") if celltype_col else "Unknown").strip() or "Unknown"
            sub  = str(getattr(row, subs_col, "")            if subs_col    else "").strip().upper()
            freq = getattr(row, freq_col, None)               if freq_col    else None
            cov  = getattr(row, cov_col,  None)               if cov_col     else None

            s = ct_stats[ct]
            s["n"] += 1

            if freq is not None and pd.notna(freq):
                fv = float(freq)
                if np.isfinite(fv) and 0 <= fv <= 1:
                    freq_hist[round(fv, 3)] += 1
                    s["freq_sum"] += fv
                    s["freq_n"]   += 1

            if cov is not None and pd.notna(cov):
                cv = int(round(float(cov)))
                if cv >= 0:
                    cov_per_celltype[ct][cv] += 1
                    s["cov_sum"] += cv
                    s["cov_n"]   += 1
                    if sub == "AG":
                        cov_ag_hist[cv] += 1
                    elif sub == "TC":
                        cov_tc_hist[cv] += 1

    return dict(
        freq_hist=dict(freq_hist),
        cov_ag_hist=dict(cov_ag_hist),
        cov_tc_hist=dict(cov_tc_hist),
        cov_per_celltype={ct: dict(h) for ct, h in cov_per_celltype.items()},
        ct_stats={ct: dict(s) for ct, s in ct_stats.items()},
    )


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(
        description="QC plots for filtered AGTC editing sites (coverage + frequency)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--input_dir", required=True, type=Path,
                    help="Directory containing *.AGTC.bed files (step 4 output)")
    ap.add_argument("--out_dir",   type=Path, default=None,
                    help="Output directory (default: input_dir/qc_agtc/)")
    ap.add_argument("--sample",    type=str, default=None,
                    help="Process only files whose name starts with this prefix (for parallel sbatch)")
    args = ap.parse_args()

    out_dir = args.out_dir or (args.input_dir / "qc_agtc")
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Output: {out_dir}")

    bed_files = sorted(args.input_dir.glob("*.AGTC.bed"))
    if not bed_files:
        raise SystemExit(f"[FATAL] No *.AGTC.bed files in {args.input_dir}")

    if args.sample:
        bed_files = [f for f in bed_files if f.name.startswith(args.sample)]
        if not bed_files:
            raise SystemExit(f"[FATAL] No files matching sample '{args.sample}'")
        out_dir = out_dir / args.sample
        out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Processing {len(bed_files)} file(s)")

    # Accumulators across all files
    all_freq_hist:        dict[float, int]          = defaultdict(int)
    all_cov_ag:           dict[int, int]            = defaultdict(int)
    all_cov_tc:           dict[int, int]            = defaultdict(int)
    all_cov_per_ct:       dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    all_ct_stats:         dict[str, dict]           = defaultdict(
        lambda: dict(n=0, freq_sum=0.0, freq_n=0, cov_sum=0.0, cov_n=0)
    )

    for bed in bed_files:
        print(f"[INFO]  {bed.name}")
        r = process_file(bed)

        for k, v in r["freq_hist"].items():
            all_freq_hist[k] += v
        for k, v in r["cov_ag_hist"].items():
            all_cov_ag[k] += v
        for k, v in r["cov_tc_hist"].items():
            all_cov_tc[k] += v
        for ct, h in r["cov_per_celltype"].items():
            for k, v in h.items():
                all_cov_per_ct[ct][k] += v
        for ct, s in r["ct_stats"].items():
            g = all_ct_stats[ct]
            g["n"]        += s["n"]
            g["freq_sum"] += s["freq_sum"]
            g["freq_n"]   += s["freq_n"]
            g["cov_sum"]  += s["cov_sum"]
            g["cov_n"]    += s["cov_n"]

    # ── Frequency distribution ────────────────────────────────────────────────
    fig = plot_freq_dist(all_freq_hist)
    if fig:
        savefig(fig, out_dir / "qc_frequency_distribution")
        print("[OK] qc_frequency_distribution.pdf/.png")

    # ── Coverage distribution — AG and TC ────────────────────────────────────
    for sub, hist in [("AG", all_cov_ag), ("TC", all_cov_tc)]:
        fig = plot_cov_dist(hist, sub=sub)
        if fig:
            savefig(fig, out_dir / f"qc_coverage_{sub}")
            print(f"[OK] qc_coverage_{sub}.pdf/.png")

    # ── Coverage per cell type ────────────────────────────────────────────────
    if all_cov_per_ct:
        fig = plot_cov_per_celltype({ct: dict(h) for ct, h in all_cov_per_ct.items()})
        if fig:
            savefig(fig, out_dir / "qc_coverage_per_celltype")
            print("[OK] qc_coverage_per_celltype.pdf/.png")

    # ── Bar charts: mean freq and mean coverage per cell type ─────────────────
    agg_rows = []
    for ct, s in sorted(all_ct_stats.items()):
        agg_rows.append({
            "celltype":  ct,
            "n_sites":   s["n"],
            "mean_freq": s["freq_sum"] / s["freq_n"] if s["freq_n"] > 0 else None,
            "mean_cov":  s["cov_sum"]  / s["cov_n"]  if s["cov_n"]  > 0 else None,
        })
    agg = pd.DataFrame(agg_rows).sort_values("n_sites", ascending=False)
    agg.to_csv(out_dir / "qc_summary_per_celltype.tsv", sep="\t", index=False)
    print("[OK] qc_summary_per_celltype.tsv")

    valid_freq = agg.dropna(subset=["mean_freq"])
    if not valid_freq.empty:
        fig = plot_bar(
            list(valid_freq["celltype"]), list(valid_freq["mean_freq"]),
            ylabel="Mean editing frequency", color="coral",
            title="Mean editing frequency per cell type (AGTC sites)",
        )
        savefig(fig, out_dir / "qc_mean_frequency_per_celltype")
        print("[OK] qc_mean_frequency_per_celltype.pdf/.png")

    valid_cov = agg.dropna(subset=["mean_cov"])
    if not valid_cov.empty:
        fig = plot_bar(
            list(valid_cov["celltype"]), list(valid_cov["mean_cov"]),
            ylabel="Mean coverage", color="seagreen",
            title="Mean coverage per cell type (AGTC sites)",
        )
        savefig(fig, out_dir / "qc_mean_coverage_per_celltype")
        print("[OK] qc_mean_coverage_per_celltype.pdf/.png")

    # ── Print summary ─────────────────────────────────────────────────────────
    total = int(agg["n_sites"].sum())
    print(f"\n{'─'*55}")
    print(f"  Total AGTC sites : {total:,}")
    print(f"  Cell types       : {len(agg)}")
    print(f"  Output           : {out_dir}")
    print(f"{'─'*55}")
    print(agg.to_string(index=False))


if __name__ == "__main__":
    main()
