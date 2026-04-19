#!/usr/bin/env python3
"""go-enrichment — GO/KEGG/Reactome enrichment via g:Profiler REST API.

Usage
-----
python go_enrichment.py --genes TP53,BRCA1,ATM --outdir results/
python go_enrichment.py --genes gene_list.txt --sources GO:BP,KEGG --outdir results/
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from tenacity import retry, stop_after_attempt, wait_exponential

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
CACHE_DIR = Path.home() / ".cache" / "genomics-skills" / "go-enrichment"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="GO/KEGG/Reactome enrichment via g:Profiler")
    p.add_argument("--genes", required=True,
                   help="Comma-separated gene symbols or path to text file (one per line)")
    p.add_argument("--organism", default="hsapiens")
    p.add_argument("--sources", default="GO:BP,KEGG,REAC",
                   help="Comma-separated sources: GO:BP, GO:MF, GO:CC, KEGG, REAC, WP")
    p.add_argument("--top-n", type=int, default=20)
    p.add_argument("--outdir", default="results/")
    p.add_argument("--dpi", type=int, default=300)
    return p


def _load_genes(genes_arg: str) -> list[str]:
    p = Path(genes_arg)
    if p.exists():
        return [g.strip() for g in p.read_text().splitlines() if g.strip()]
    return [g.strip() for g in genes_arg.split(",") if g.strip()]


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
def _query_gprofiler(genes: list[str], organism: str, sources: list[str]) -> list[dict]:
    cache_key = hashlib.md5("|".join(sorted(genes) + sources).encode()).hexdigest()
    cache_path = CACHE_DIR / f"{cache_key}.json"
    if cache_path.exists():
        print(f"[cache] Loading enrichment results from cache")
        return json.loads(cache_path.read_text())

    print(f"[fetch] g:Profiler enrichment for {len(genes)} genes ...")
    payload = {
        "organism": organism,
        "query": genes,
        "sources": sources,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
    }
    resp = requests.post(GPROFILER_URL, json=payload, timeout=30)
    resp.raise_for_status()
    results = resp.json().get("result", [])
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(results))
    return results


def _results_to_df(results: list[dict]) -> pd.DataFrame:
    rows = []
    for r in results:
        rows.append({
            "source": r.get("source", ""),
            "term_id": r.get("native", ""),
            "term_name": r.get("name", ""),
            "p_value": r.get("p_value", 1.0),
            "significant": r.get("significant", False),
            "query_size": r.get("query_size", 0),
            "intersection_size": r.get("intersection_size", 0),
            "term_size": r.get("term_size", 0),
            "genes": ",".join(r.get("intersections", {}).get("gene", [])),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["gene_ratio"] = df["intersection_size"] / df["query_size"]
    df["neg_log10_p"] = -np.log10(df["p_value"].clip(1e-300, 1))
    return df.sort_values("p_value")


def _plot_bubble(df: pd.DataFrame, top_n: int, outdir: Path, dpi: int) -> None:
    top = df.head(top_n).copy()
    if top.empty:
        print("[warn] No significant terms to plot")
        return

    source_colors = {"GO:BP": "#3498db", "GO:MF": "#2ecc71", "GO:CC": "#e74c3c",
                     "KEGG": "#f39c12", "REAC": "#9b59b6", "WP": "#1abc9c"}
    colors = top["source"].map(lambda s: source_colors.get(s, "#95a5a6"))

    fig, ax = plt.subplots(figsize=(10, max(5, top_n // 2)))
    sc = ax.scatter(
        top["gene_ratio"],
        range(len(top)),
        s=top["intersection_size"] * 15,
        c=top["neg_log10_p"],
        cmap="RdYlBu_r",
        alpha=0.85,
        edgecolors="white",
        linewidths=0.5,
    )
    plt.colorbar(sc, ax=ax, label="-log₁₀(adj. p-value)", shrink=0.6)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(
        [f"[{r['source']}] {r['term_name'][:55]}" for _, r in top.iterrows()],
        fontsize=8,
    )
    ax.set_xlabel("Gene ratio", fontsize=11)
    ax.set_title(f"Pathway enrichment — top {top_n} terms", fontsize=12, fontweight="bold")
    plt.tight_layout()
    for ext in ("png", "svg"):
        out = outdir / f"enrichment_bubble.{ext}"
        fig.savefig(out, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[output] {out}")
    plt.close(fig)


def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = _load_genes(args.genes)
    sources = [s.strip() for s in args.sources.split(",")]
    print(f"[info] {len(genes)} genes, sources: {sources}")

    results = _query_gprofiler(genes, args.organism, sources)
    df = _results_to_df(results)

    tsv = outdir / "enrichment_results.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"[output] {tsv} ({len(df)} terms)")

    _plot_bubble(df[df["significant"]], args.top_n, outdir, args.dpi)
    print(f"\n✓ go-enrichment complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
