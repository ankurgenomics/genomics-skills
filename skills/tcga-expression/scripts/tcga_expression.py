#!/usr/bin/env python3
"""tcga-expression — pan-cancer and tumor-vs-normal expression from TCGA via cBioPortal.

Uses the cBioPortal REST API (no auth, no file downloads) to fetch
RNA-seq v2 RSEM expression for a gene across TCGA pan-cancer atlas 2018.

Usage
-----
python tcga_expression.py --gene TP53 --mode pan-cancer --outdir results/
python tcga_expression.py --gene TP53 --mode tumor-vs-normal --cancer LUAD --outdir results/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy import stats

# ── constants ────────────────────────────────────────────────────────────────

CBIO_BASE    = "https://www.cbioportal.org/api"
CACHE_DIR    = Path.home() / ".cache" / "genomics-skills" / "tcga-expression"
TCGA_PALETTE = sns.color_palette("tab20", 20)

TCGA_CANCERS = [
    "acc","blca","brca","cesc","chol","coad","dlbc","esca","gbm","hnsc",
    "kich","kirc","kirp","laml","lgg","lihc","luad","lusc","meso","ov",
    "paad","pcpg","prad","sarc","skcm","stad","tgct","thca","thym",
    "ucec","ucs","uvm",
]


# ── CLI ──────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="TCGA gene expression via cBioPortal")
    p.add_argument("--gene", required=True, help="HGNC gene symbol, e.g. TP53")
    p.add_argument("--mode", default="pan-cancer",
                   choices=["pan-cancer", "tumor-vs-normal"])
    p.add_argument("--cancer", default="all",
                   help="Cancer abbreviation for tumor-vs-normal, e.g. LUAD")
    p.add_argument("--outdir", default="results/")
    p.add_argument("--dpi", type=int, default=300)
    p.add_argument("--top-n", type=int, default=15,
                   help="Max cancer types in pan-cancer plot")
    return p


# ── helpers ───────────────────────────────────────────────────────────────────

def _entrez_id(gene: str) -> int:
    resp = requests.get(f"{CBIO_BASE}/genes/{gene}", timeout=10)
    resp.raise_for_status()
    return resp.json()["entrezGeneId"]


def _fetch_one(cancer_prefix: str, entrez_id: int) -> list[dict]:
    profile = f"{cancer_prefix}_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"
    sample_list = f"{cancer_prefix}_tcga_pan_can_atlas_2018_all"
    try:
        resp = requests.get(
            f"{CBIO_BASE}/molecular-profiles/{profile}/molecular-data",
            params={"entrezGeneId": entrez_id, "sampleListId": sample_list},
            headers={"Accept": "application/json"}, timeout=15,
        )
        if resp.status_code == 404:
            return []
        resp.raise_for_status()
        return resp.json()
    except Exception as e:
        print(f"[warn] {cancer_prefix}: {e}", file=sys.stderr)
        return []


# ── fetcher ───────────────────────────────────────────────────────────────────

def _fetch_expression(gene: str, cancer_filter: str | None = None) -> pd.DataFrame:
    cache_key = f"{gene}_{cancer_filter or 'pan'}.parquet"
    cache_path = CACHE_DIR / cache_key
    if cache_path.exists():
        print(f"[cache] {cache_path}")
        return pd.read_parquet(cache_path)

    print(f"[fetch] cBioPortal RNA-seq for {gene} ...")
    entrez = _entrez_id(gene)

    cancers = ([cancer_filter.lower().replace("tcga-", "")]
               if cancer_filter and cancer_filter.lower() != "all"
               else TCGA_CANCERS)

    rows = []
    for cancer in cancers:
        data = _fetch_one(cancer, entrez)
        project = f"TCGA-{cancer.upper()}"
        for rec in data:
            val = rec.get("value")
            if val is None:
                continue
            rows.append({"project": project, "sample_type": "Primary Tumor",
                         "expression": np.log2(max(float(val), 0) + 1)})
        if data:
            print(f"  {project}: {len(data)} samples")

    df = pd.DataFrame(rows)
    if not df.empty:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache_path, index=False)
        print(f"[info] {len(df)} samples · {df['project'].nunique()} projects")
    return df


# ── plotting ──────────────────────────────────────────────────────────────────

def _plot_pancancer(df: pd.DataFrame, gene: str, top_n: int, outdir: Path, dpi: int) -> None:
    medians = df.groupby("project")["expression"].median().sort_values(ascending=False)
    order = medians.head(top_n).index.tolist()
    df_top = df[df["project"].isin(order)].copy()
    df_top["_hue"] = df_top["project"]
    palette_map = {p: TCGA_PALETTE[i % 20] for i, p in enumerate(order)}

    fig, ax = plt.subplots(figsize=(max(10, len(order)), 5))
    sns.boxplot(data=df_top, x="project", y="expression", hue="_hue", order=order,
                palette=palette_map, legend=False, linewidth=0.8, fliersize=2, ax=ax)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=9)
    ax.set_xlabel("TCGA project", fontsize=11)
    ax.set_ylabel("Expression (log₂ RSEM+1)", fontsize=11)
    ax.set_title(f"{gene} — Pan-cancer expression (TCGA)", fontsize=13, fontweight="bold")
    plt.tight_layout()
    for ext in ("png", "svg"):
        out = outdir / f"{gene}_pancancer_expression.{ext}"
        fig.savefig(out, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[output] {out}")
    plt.close(fig)


def _plot_tumor_vs_normal(df: pd.DataFrame, gene: str, cancer: str,
                           outdir: Path, dpi: int) -> None:
    project = f"TCGA-{cancer.upper().replace('TCGA-', '')}"
    df_c = df[df["project"] == project].copy()
    if df_c.empty:
        print(f"[warn] No data for {project}", file=sys.stderr)
        return
    fig, ax = plt.subplots(figsize=(4, 5))
    sns.boxplot(data=df_c, y="expression", color="#e74c3c", linewidth=1.0, ax=ax)
    sns.stripplot(data=df_c, y="expression", color="black", alpha=0.3, size=3, ax=ax)
    ax.set_ylabel("Expression (log₂ RSEM+1)", fontsize=11)
    ax.set_title(f"{gene} — {project}", fontsize=12, fontweight="bold")
    plt.tight_layout()
    out = outdir / f"{gene}_{project}_expression.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"[output] {out}")
    plt.close(fig)


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cancer = args.cancer if args.mode == "tumor-vs-normal" else None
    df = _fetch_expression(args.gene, cancer_filter=cancer)

    if df.empty:
        print("[error] No expression data retrieved.", file=sys.stderr)
        sys.exit(1)

    tsv = outdir / f"{args.gene}_expression.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"[output] {tsv}")

    if args.mode == "pan-cancer":
        _plot_pancancer(df, args.gene, args.top_n, outdir, args.dpi)
    else:
        _plot_tumor_vs_normal(df, args.gene, args.cancer, outdir, args.dpi)

    print(f"\n✓ tcga-expression complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
