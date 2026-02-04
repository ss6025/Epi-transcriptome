#!/usr/bin/env python3
import os, math, argparse, re
from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def read_header(path):
    with open(path, "r") as f:
        return f.readline().rstrip("\n").split("\t")


def find_col(hdr, want):
    """
    Find a column by exact match or substring match.
    """
    if want in hdr:
        return want
    w = want.lower()
    for c in hdr:
        if w == c.lower():
            return c
    for c in hdr:
        if w in c.lower():
            return c
    return None


def subs_mask(series_subs: pd.Series, tokens=("AG","TC")):
    """
    True for rows where AllSubs contains any token as a whitespace-delimited item.
    Handles multi-subs like 'GT GA'. Ignores '-', NA.
    """
    tokens = [t.upper() for t in tokens]
    s = series_subs.astype(str).str.upper()
    ok = (~s.isin(["-", "NA", "", "NAN"])) & s.notna()
    pat = r"(^|\s)(" + "|".join(map(re.escape, tokens)) + r")($|\s)"
    return ok & s.str.contains(pat, regex=True)


def weighted_quantile_from_hist(xs, ys, q):
    cum = np.cumsum(ys)
    n = int(cum[-1])
    if n == 0:
        return None
    idx = int(np.searchsorted(cum, q * n, side="left"))
    return float(xs[min(idx, len(xs) - 1)])


def stream_counts_and_covhist_tokenized(
    bed_path,
    col_sub,
    col_freq,
    col_cov,
    min_freq=0.0,
    min_cov=0.0,
    chunksize=500_000,
    cov_round=None,   # None for "noisy" coverage
):
    """
    Tokenize substitutions: "GT GA" -> ["GT","GA"], "AG" -> ["AG"]
    Returns:
      sub_counts: dict[sub] -> count
      cov_hist: dict[sub] -> dict[cov_rounded] -> count
      n_kept: total rows kept after filtering
    """
    sub_counts = defaultdict(int)
    cov_hist = defaultdict(lambda: defaultdict(int))
    n_kept_rows = 0

    usecols = [col_sub, col_freq, col_cov]

    for chunk in pd.read_csv(
        bed_path,
        sep="\t",
        usecols=usecols,
        comment="#",
        chunksize=chunksize,
        low_memory=True,
    ):
        freq = pd.to_numeric(chunk[col_freq], errors="coerce").to_numpy()
        cov  = pd.to_numeric(chunk[col_cov],  errors="coerce").to_numpy()
        subs = chunk[col_sub].astype(str).to_numpy()

        m = np.isfinite(freq) & np.isfinite(cov) & (freq >= min_freq) & (cov >= min_cov)
        subs_m = subs[m]
        cov_m  = cov[m]
        if subs_m.size == 0:
            continue

        # keep coverage as-is unless explicitly rounding
        covv = cov_m if cov_round is None else np.round(cov_m, int(cov_round))

        # Tokenize: "GT GA" -> ["GT","GA"], "AG" -> ["AG"]
        tokens = []
        cov_tokens = []
        for s, c in zip(subs_m, covv):
            s = str(s).strip().upper()
            if s in ("-", "NA", "", "NAN"):
                continue
            for tok in s.split():
                tok = tok.strip().upper()
                if len(tok) == 2 and set(tok) <= set("ACGT"):
                    tokens.append(tok)
                    cov_tokens.append(float(c))

        if not tokens:
            continue

        n_kept_rows += int(subs_m.shape[0])

        tokens = np.asarray(tokens, dtype=object)
        cov_tokens = np.asarray(cov_tokens, dtype=float)

        # counts
        u_subs, u_counts = np.unique(tokens, return_counts=True)
        for s, c in zip(u_subs, u_counts):
            sub_counts[s] += int(c)

        # coverage histogram per substitution
        for s in u_subs:
            mm = (tokens == s)
            vals = cov_tokens[mm]
            u_cov, c_cov = np.unique(vals, return_counts=True)
            h = cov_hist[s]
            for v, c in zip(u_cov, c_cov):
                h[float(v)] += int(c)

    return sub_counts, cov_hist, n_kept_rows


def stream_frequency_hist_filtered(
    bed_path,
    col_freq,
    col_sub,
    tokens=("AG","TC"),
    chunksize=500_000,
    round_dp=3,
    min_val=0.0,
):
    """
    Build frequency histogram for rows matching tokens.
    """
    hist = defaultdict(int)
    usecols = [col_freq, col_sub]

    for chunk in pd.read_csv(
        bed_path,
        sep="\t",
        usecols=usecols,
        comment="#",
        chunksize=chunksize,
        low_memory=True,
    ):
        m = subs_mask(chunk[col_sub], tokens=tokens)
        if not m.any():
            continue

        arr = pd.to_numeric(chunk.loc[m, col_freq], errors="coerce").to_numpy()
        arr = arr[np.isfinite(arr)]
        if min_val is not None:
            arr = arr[arr > float(min_val)]
        if arr.size == 0:
            continue

        if round_dp is not None:
            arr = np.round(arr, int(round_dp))

        vals, counts = np.unique(arr, return_counts=True)
        for v, c in zip(vals, counts):
            hist[float(v)] += int(c)

    return hist


def hist_to_arrays(hist):
    xs = np.array(sorted(hist.keys()), dtype=float)
    ys = np.array([hist[x] for x in xs], dtype=int)
    return xs, ys


def weighted_quantile(xs, ys, q):
    cum = np.cumsum(ys)
    n = int(cum[-1])
    if n == 0:
        return None
    idx = int(np.searchsorted(cum, q * n, side="left"))
    return float(xs[min(idx, len(xs) - 1)])


def plot_mismatches(sub_counts, title=None, top_n=30, exclude=None):
    """
    Plot substitution counts bar chart.
    """
    exclude = {str(x).upper().strip() for x in (exclude or set())}

    # normalize keys + merge duplicates
    cleaned = {}
    for k, v in sub_counts.items():
        ku = str(k).upper().strip()
        if ku in exclude:
            continue
        cleaned[ku] = cleaned.get(ku, 0) + int(v)

    items = sorted(cleaned.items(), key=lambda x: x[1], reverse=True)
    if not items:
        raise ValueError("No substitutions after filtering (nothing to plot).")
    items = items[:top_n]

    labels = [k for k, _ in items]
    counts = np.array([v for _, v in items], dtype=np.int64)

    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(labels)), 4.8))
    ax.bar(labels, counts)
    ax.set_ylabel("Count")
    ax.set_title(title or "Substitution counts")
    ax.tick_params(axis="x", rotation=90)
    ax.grid(axis="y", alpha=0.2)

    # no scientific notation, no offset
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x:,.0f}'))

    fig.tight_layout()
    return fig


def plot_cov_overlay_step(cov_hist, subs=("AG","TC"), title=None, x_max=200, logy=False):
    """
    Plot coverage overlay for multiple substitutions (step plot).
    """
    fig, ax = plt.subplots(figsize=(7.6, 4.8))

    keys = {str(k).upper(): k for k in cov_hist.keys()}

    for sub in subs:
        su = str(sub).upper()
        if su not in keys:
            continue
        h = cov_hist[keys[su]]
        xs = np.array(sorted(h.keys()), dtype=float)
        ys = np.array([h[float(x)] for x in xs], dtype=int)

        # cap x-axis only
        m = xs <= x_max
        xs, ys = xs[m], ys[m]
        if xs.size == 0:
            continue

        ax.step(xs, ys, where="mid", linewidth=1.6, label=su)

    ax.set_xlim(0, x_max)
    ax.set_xlabel("Coverage-q30")
    ax.set_ylabel("Sites")
    ax.set_title(title or f"Coverage-q30 for {' vs '.join(subs)} (x≤{x_max})")
    if logy:
        ax.set_yscale("log")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    return fig


def plot_frequency_histogram(hist, title=None, zoom_q=0.999, logy=False):
    """
    Plot frequency histogram with full and zoomed views.
    """
    xs, ys = hist_to_arrays(hist)
    if xs.size == 0:
        return None, None

    med = weighted_quantile(xs, ys, 0.5)
    p99 = weighted_quantile(xs, ys, 0.99)
    zoom_xmax = weighted_quantile(xs, ys, zoom_q)

    # Full plot
    fig1, ax1 = plt.subplots(figsize=(10, 2.8))
    ax1.step(xs, ys, where="mid", linewidth=1.5)
    if med is not None:
        ax1.axvline(med, linestyle="--", linewidth=1, label="median")
    if p99 is not None:
        ax1.axvline(p99, linestyle=":", linewidth=1, label="p99")
    ax1.set_title(f"{title}: Frequency (full)" if title else "Frequency (full)")
    ax1.set_xlabel("Frequency")
    ax1.set_ylabel("rows")
    if logy:
        ax1.set_yscale("log")
    ax1.grid(alpha=0.25)
    ax1.legend(frameon=False)
    fig1.tight_layout()

    # Zoomed plot
    fig2 = None
    if zoom_xmax is not None:
        m = xs <= zoom_xmax
        fig2, ax2 = plt.subplots(figsize=(10, 2.8))
        ax2.step(xs[m], ys[m], where="mid", linewidth=1.5)
        ax2.set_title(f"{title}: Frequency (≤ q={zoom_q})" if title else f"Frequency (≤ q={zoom_q})")
        ax2.set_xlabel("Frequency")
        ax2.set_ylabel("rows")
        if logy:
            ax2.set_yscale("log")
        ax2.grid(alpha=0.25)
        fig2.tight_layout()

    return fig1, fig2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True, help="Input BED with header")
    ap.add_argument("--outdir", required=True, help="Output directory for plots + TSVs")
    ap.add_argument("--min_freq", type=float, default=0.0, help="Keep rows with Frequency >= this")
    ap.add_argument("--min_cov", type=float, default=0.0, help="Keep rows with Coverage-q30 >= this")
    ap.add_argument("--chunksize", type=int, default=500_000)
    ap.add_argument("--logy", action="store_true", help="log scale on y-axis")
    ap.add_argument("--cov_round", type=int, default=None, help="Round coverage-q30 to N dp (None = no rounding)")
    ap.add_argument("--tokens", nargs="+", default=["AG", "TC"], help="Substitution tokens to analyze (default: AG TC)")
    ap.add_argument("--exclude", nargs="+", default=["TA", "AT", "GA", "GT", "CT", "CA"], help="Substitutions to exclude from counts plot")
    ap.add_argument("--top_n", type=int, default=20, help="Top N substitutions to show in counts plot")
    ap.add_argument("--cov_xmax", type=int, default=200, help="Maximum x-axis for coverage plot")
    ap.add_argument("--freq_zoom_q", type=float, default=0.999, help="Quantile for frequency zoom plot")
    args = ap.parse_args()

    if args.outdir.strip() == "":
        raise SystemExit("[FATAL] --outdir is empty.")
    os.makedirs(args.outdir, exist_ok=True)

    hdr = read_header(args.bed)
    col_sub = find_col(hdr, "AllSubs")
    col_freq = find_col(hdr, "Frequency")
    col_cov = find_col(hdr, "Coverage-q30")

    missing = [("AllSubs", col_sub), ("Frequency", col_freq), ("Coverage-q30", col_cov)]
    missing = [name for name, got in missing if got is None]
    if missing:
        raise SystemExit(f"[FATAL] Missing required columns: {missing}\nHeader: {hdr}")

    sample = os.path.basename(args.bed).rsplit(".", 1)[0]

    # Compute substitution counts and coverage histograms (with tokenization)
    print(f"[INFO] Processing {args.bed}...")
    sub_counts, cov_hist, n_kept = stream_counts_and_covhist_tokenized(
        args.bed,
        col_sub=col_sub,
        col_freq=col_freq,
        col_cov=col_cov,
        min_freq=args.min_freq,
        min_cov=args.min_cov,
        chunksize=args.chunksize,
        cov_round=args.cov_round,
    )

    print(f"[INFO] Kept rows: {n_kept:,}")
    print(f"[INFO] AG count (tokenized): {sub_counts.get('AG', 0):,}")
    print(f"[INFO] TC count (tokenized): {sub_counts.get('TC', 0):,}")

    tag = f"minFreq{args.min_freq}.minCov{args.min_cov}"

    # 1. Substitution counts plot
    out_counts_tsv = os.path.join(args.outdir, f"{sample}.subs_counts.{tag}.tsv")
    pd.DataFrame(sorted(sub_counts.items(), key=lambda x: x[1], reverse=True),
                 columns=["AllSubs", "Count"]).to_csv(out_counts_tsv, sep="\t", index=False)

    try:
        fig_counts = plot_mismatches(
            sub_counts,
            title=f"{sample}: substitution counts (tokenized)",
            exclude=args.exclude,
            top_n=args.top_n
        )
        out_counts_png = os.path.join(args.outdir, f"{sample}.subs_counts.{tag}.png")
        fig_counts.savefig(out_counts_png, dpi=250, bbox_inches="tight")
        plt.close(fig_counts)
        print(f"[OK] Saved: {out_counts_png}")
    except Exception as e:
        print(f"[WARN] Failed to create counts plot: {e}")

    # 2. Coverage overlay plot (AG vs TC)
    tokens_tuple = tuple(args.tokens)
    try:
        fig_cov = plot_cov_overlay_step(
            cov_hist,
            subs=tokens_tuple,
            title=f"{sample}: Coverage-q30 (x≤{args.cov_xmax})",
            x_max=args.cov_xmax,
            logy=args.logy
        )
        out_cov_png = os.path.join(args.outdir, f"{sample}.coverage_overlay_{'_'.join(args.tokens)}.{tag}.png")
        fig_cov.savefig(out_cov_png, dpi=250, bbox_inches="tight")
        plt.close(fig_cov)
        print(f"[OK] Saved: {out_cov_png}")
    except Exception as e:
        print(f"[WARN] Failed to create coverage overlay plot: {e}")

    # 3. Frequency histogram (for AG/TC only)
    try:
        freq_hist = stream_frequency_hist_filtered(
            args.bed,
            col_freq=col_freq,
            col_sub=col_sub,
            tokens=tokens_tuple,
            chunksize=args.chunksize,
            round_dp=3,
            min_val=0.0,
        )
        if freq_hist:
            fig_freq_full, fig_freq_zoom = plot_frequency_histogram(
                freq_hist,
                title=f"{sample}: Frequency for {tokens_tuple} rows",
                zoom_q=args.freq_zoom_q,
                logy=args.logy
            )
            if fig_freq_full:
                out_freq_full = os.path.join(args.outdir, f"{sample}.frequency_full.{tag}.png")
                fig_freq_full.savefig(out_freq_full, dpi=250, bbox_inches="tight")
                plt.close(fig_freq_full)
                print(f"[OK] Saved: {out_freq_full}")
            if fig_freq_zoom:
                out_freq_zoom = os.path.join(args.outdir, f"{sample}.frequency_zoom.{tag}.png")
                fig_freq_zoom.savefig(out_freq_zoom, dpi=250, bbox_inches="tight")
                plt.close(fig_freq_zoom)
                print(f"[OK] Saved: {out_freq_zoom}")
    except Exception as e:
        print(f"[WARN] Failed to create frequency histogram: {e}")

    print("[OK] QC analysis complete")
    print("[OK] wrote:")
    print(" ", out_counts_tsv)
    if 'out_counts_png' in locals():
        print(" ", out_counts_png)
    if 'out_cov_png' in locals():
        print(" ", out_cov_png)


if __name__ == "__main__":
    main()
