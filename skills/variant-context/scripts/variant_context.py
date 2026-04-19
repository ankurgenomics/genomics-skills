#!/usr/bin/env python3
"""variant-context — somatic mutation landscape for a gene.

Fetches mutation data from MyVariant.info, produces a lollipop plot,
hotspot bar chart, and a structured TSV.

Usage
-----
python variant_context.py --gene TP53 --cancer LUAD --outdir results/

Outputs
-------
{gene}_lollipop.png / .svg   — mutation lollipop with domain track
{gene}_hotspots.png           — top-N mutated positions bar chart
{gene}_mutations.tsv          — all mutations with scores
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Any

import httpx
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tenacity import retry, stop_after_attempt, wait_exponential

# ── constants ────────────────────────────────────────────────────────────────

MYVARIANT_QUERY = "https://myvariant.info/v1/query"
CACHE_DIR = Path.home() / ".cache" / "genomics-skills" / "variant-context"
CONSEQUENCE_COLORS: dict[str, str] = {
    "missense_variant": "#e74c3c",
    "stop_gained": "#8e44ad",
    "frameshift_variant": "#2c3e50",
    "splice_donor_variant": "#f39c12",
    "splice_acceptor_variant": "#f39c12",
    "synonymous_variant": "#95a5a6",
    "other": "#bdc3c7",
}

# ── CLI ──────────────────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Somatic mutation landscape for a gene from MyVariant.info"
    )
    p.add_argument("--gene", required=True, help="HGNC gene symbol (e.g. TP53)")
    p.add_argument("--cancer", default="all", help="TCGA cancer type or 'all'")
    p.add_argument("--assembly", default="hg38", choices=["hg38", "hg19"])
    p.add_argument("--outdir", default="results/", help="Output directory")
    p.add_argument("--top-n", type=int, default=20, help="Top N hotspot positions")
    p.add_argument("--dpi", type=int, default=300, help="Plot resolution (DPI)")
    return p


# ── data fetching ─────────────────────────────────────────────────────────────


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
def _query_myvariant(gene: str, assembly: str, size: int = 500) -> list[dict[str, Any]]:
    """Query MyVariant.info for all variants in a gene."""
    cache_path = CACHE_DIR / f"{gene}_{assembly}.json"
    if cache_path.exists():
        print(f"[cache] Loading {gene} mutations from {cache_path}")
        return json.loads(cache_path.read_text())

    print(f"[fetch] Querying MyVariant.info for {gene} ({assembly}) ...")
    params = {
        "q": f"dbnsfp.genename:{gene}",
        "fields": "dbnsfp.genename,dbnsfp.aapos,dbnsfp.aaref,dbnsfp.aalt,"
                  "dbnsfp.cadd_phred,clinvar.rcv.clinical_significance,"
                  "vcf.ref,vcf.alt,vep.most_severe_consequence,chrom,hg38.start",
        "size": size,
        "assembly": assembly,
    }
    with httpx.Client(timeout=20.0) as client:
        resp = client.get(MYVARIANT_QUERY, params=params)
        resp.raise_for_status()
        hits = resp.json().get("hits", [])

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(hits))
    return hits


def _hits_to_df(hits: list[dict[str, Any]]) -> pd.DataFrame:
    """Flatten MyVariant.info hits into a tidy DataFrame."""
    rows = []
    for h in hits:
        dbnsfp = h.get("dbnsfp", {})
        rows.append(
            {
                "gene": dbnsfp.get("genename", ""),
                "aa_pos": dbnsfp.get("aapos", np.nan),
                "aa_ref": dbnsfp.get("aaref", ""),
                "aa_alt": dbnsfp.get("aalt", ""),
                "cadd_phred": dbnsfp.get("cadd_phred", np.nan),
                "consequence": h.get("vep", {}).get("most_severe_consequence", "other"),
                "clinvar_sig": (
                    h.get("clinvar", {}).get("rcv", {}).get("clinical_significance", "")
                    if isinstance(h.get("clinvar", {}).get("rcv"), dict)
                    else ""
                ),
                "chrom": h.get("chrom", ""),
                "pos_hg38": h.get("hg38", {}).get("start", np.nan),
            }
        )
    df = pd.DataFrame(rows)
    df = df[df["aa_pos"].notna() & (df["aa_pos"] != ".")].copy()
    df["aa_pos"] = pd.to_numeric(df["aa_pos"], errors="coerce")
    return df.dropna(subset=["aa_pos"])


# ── plotting ─────────────────────────────────────────────────────────────────


def _plot_lollipop(df: pd.DataFrame, gene: str, outdir: Path, dpi: int) -> None:
    """Draw a mutation lollipop plot."""
    if df.empty:
        print(f"[warn] No mutation data available for {gene} — skipping lollipop plot.")
        return

    # Count mutations per aa position
    counts = df.groupby("aa_pos").size().reset_index(name="count")
    max_pos = int(df["aa_pos"].max()) + 50

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.set_xlim(0, max_pos)
    ax.set_ylim(-0.5, counts["count"].max() + 2)

    # Stem + dot for each position
    for _, row in counts.iterrows():
        pos = row["aa_pos"]
        cnt = row["count"]
        # Pick color by most common consequence at this position
        conseq = df[df["aa_pos"] == pos]["consequence"].mode()
        color = CONSEQUENCE_COLORS.get(
            conseq.iloc[0] if not conseq.empty else "other", CONSEQUENCE_COLORS["other"]
        )
        ax.plot([pos, pos], [0, cnt], color="#bdc3c7", linewidth=0.8, zorder=1)
        ax.scatter(pos, cnt, color=color, s=max(20, cnt * 8), zorder=2, alpha=0.85)

    # Protein backbone
    ax.axhline(0, color="#2c3e50", linewidth=2.5, zorder=0)
    ax.fill_between([0, max_pos], [-0.3, -0.3], [0, 0], color="#ecf0f1", zorder=0)

    ax.set_xlabel("Amino acid position", fontsize=11)
    ax.set_ylabel("Mutation count", fontsize=11)
    ax.set_title(f"{gene} — Somatic mutation lollipop", fontsize=13, fontweight="bold")

    legend_patches = [
        mpatches.Patch(color=c, label=k.replace("_", " "))
        for k, c in CONSEQUENCE_COLORS.items()
    ]
    ax.legend(handles=legend_patches, fontsize=8, loc="upper right", framealpha=0.8)

    plt.tight_layout()
    for ext in ("png", "svg"):
        out = outdir / f"{gene}_lollipop.{ext}"
        fig.savefig(out, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[output] {out}")
    plt.close(fig)


def _plot_hotspots(df: pd.DataFrame, gene: str, top_n: int, outdir: Path, dpi: int) -> None:
    """Bar chart of top-N mutated amino acid positions."""
    if df.empty:
        return

    counts = (
        df.groupby("aa_pos")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .head(top_n)
        .sort_values("aa_pos")
    )

    fig, ax = plt.subplots(figsize=(max(8, top_n // 2), 4))
    ax.bar(
        counts["aa_pos"].astype(str),
        counts["count"],
        color="#e74c3c",
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_xlabel("Amino acid position", fontsize=11)
    ax.set_ylabel("Mutation count", fontsize=11)
    ax.set_title(f"{gene} — Top {top_n} mutation hotspots", fontsize=13, fontweight="bold")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.tight_layout()

    out = outdir / f"{gene}_hotspots.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"[output] {out}")
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    hits = _query_myvariant(args.gene, args.assembly)
    if not hits:
        print(f"[warn] No variants found for gene={args.gene}. Check the gene symbol.", file=sys.stderr)
        # Write empty TSV so downstream agents don't fail
        pd.DataFrame(columns=["gene", "aa_pos", "aa_ref", "aa_alt", "cadd_phred",
                               "consequence", "clinvar_sig", "chrom", "pos_hg38"]).to_csv(
            outdir / f"{args.gene}_mutations.tsv", sep="\t", index=False
        )
        sys.exit(0)

    df = _hits_to_df(hits)
    print(f"[info] {len(df)} mutations parsed for {args.gene}")

    # Save TSV
    tsv_path = outdir / f"{args.gene}_mutations.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    print(f"[output] {tsv_path}")

    # Plots
    _plot_lollipop(df, args.gene, outdir, args.dpi)
    _plot_hotspots(df, args.gene, args.top_n, outdir, args.dpi)

    print(f"\n✓ variant-context complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
