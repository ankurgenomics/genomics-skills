#!/usr/bin/env python3
# Author: Ankur Sharma | Agentic AI · Machine Learning · Bioinformatics · Data Science
# GitHub: https://github.com/ankurgenomics | Portfolio: agentic-genomics
"""plot-volcano — publication-quality volcano plot (300 DPI PNG + SVG).

Usage
-----
python plot_volcano.py --input de_results.tsv --title "Treated vs Control"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

COLORS = {"up": "#c0392b", "down": "#2980b9", "ns": "#bdc3c7"}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Publication-quality volcano plot")
    p.add_argument("--input", required=True, help="TSV/CSV with FC + p-value columns")
    p.add_argument("--fc-col",    default="log2FoldChange")
    p.add_argument("--pval-col",  default="padj")
    p.add_argument("--label-col", default="gene")
    p.add_argument("--fc-thresh", type=float, default=1.0)
    p.add_argument("--p-thresh",  type=float, default=0.05)
    p.add_argument("--top-label", type=int,   default=15)
    p.add_argument("--title",     default="Volcano plot")
    p.add_argument("--outdir",    default="results/")
    p.add_argument("--dpi",       type=int,   default=300)
    return p


def _load(path: str, fc_col: str, pval_col: str) -> pd.DataFrame:
    sep = "\t" if path.endswith(".tsv") else ","
    df = pd.read_csv(path, sep=sep)
    missing = [c for c in (fc_col, pval_col) if c not in df.columns]
    if missing:
        raise ValueError(f"Columns not found: {missing}. Available: {df.columns.tolist()}")
    df = df.dropna(subset=[fc_col, pval_col])
    df["_neg_log10_p"] = -np.log10(df[pval_col].clip(1e-300, 1))
    return df


def _classify(df: pd.DataFrame, fc_col: str, pval_col: str,
               fc_thresh: float, p_thresh: float) -> pd.DataFrame:
    up   = (df[fc_col] >=  fc_thresh) & (df[pval_col] <= p_thresh)
    down = (df[fc_col] <= -fc_thresh) & (df[pval_col] <= p_thresh)
    df["significance"] = "ns"
    df.loc[up,   "significance"] = "up"
    df.loc[down, "significance"] = "down"
    return df


def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = _load(args.input, args.fc_col, args.pval_col)
    df = _classify(df, args.fc_col, args.pval_col, args.fc_thresh, args.p_thresh)

    n_up   = (df["significance"] == "up").sum()
    n_down = (df["significance"] == "down").sum()
    print(f"[info] Up: {n_up}  Down: {n_down}  NS: {len(df) - n_up - n_down}")

    fig, ax = plt.subplots(figsize=(9, 7))

    for sig, color in COLORS.items():
        sub = df[df["significance"] == sig]
        ax.scatter(sub[args.fc_col], sub["_neg_log10_p"],
                   c=color, s=8 if sig == "ns" else 18,
                   alpha=0.5 if sig == "ns" else 0.8,
                   label=f"{sig.upper()} (n={len(sub)})", rasterized=True)

    # Threshold lines
    ax.axvline( args.fc_thresh,  color="gray", linestyle="--", linewidth=0.8)
    ax.axvline(-args.fc_thresh,  color="gray", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(args.p_thresh), color="gray", linestyle="--", linewidth=0.8)

    # Top gene labels (non-overlapping)
    sig_df = df[df["significance"] != "ns"].copy()
    top = sig_df.nlargest(args.top_label, "_neg_log10_p")
    if args.label_col in df.columns:
        texts = [
            ax.text(row[args.fc_col], row["_neg_log10_p"], row[args.label_col],
                    fontsize=7, ha="center")
            for _, row in top.iterrows()
        ]
        adjust_text(texts, ax=ax, arrowprops={"arrowstyle": "-", "color": "gray", "lw": 0.5})

    ax.set_xlabel(f"log₂ Fold Change", fontsize=12)
    ax.set_ylabel("-log₁₀ adjusted p-value", fontsize=12)
    ax.set_title(args.title, fontsize=14, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.8)

    # Quadrant counts
    ax.text(0.98, 0.98, f"↑ {n_up}  ↓ {n_down}", transform=ax.transAxes,
            ha="right", va="top", fontsize=10,
            bbox={"boxstyle": "round", "fc": "white", "ec": "lightgray", "alpha": 0.9})

    plt.tight_layout()
    for ext in ("png", "svg"):
        out = outdir / f"volcano.{ext}"
        fig.savefig(out, dpi=args.dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[output] {out}")
    plt.close(fig)

    ann_tsv = outdir / "volcano_annotated.tsv"
    df.drop(columns=["_neg_log10_p"]).to_csv(ann_tsv, sep="\t", index=False)
    print(f"[output] {ann_tsv}")
    print(f"\n✓ plot-volcano complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
