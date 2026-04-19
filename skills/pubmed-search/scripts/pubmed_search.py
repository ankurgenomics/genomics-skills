#!/usr/bin/env python3
# Author: Ankur Sharma | Agentic AI · Machine Learning · Bioinformatics · Data Science
# GitHub: https://github.com/ankurgenomics | Portfolio: agentic-genomics
"""pubmed-search — literature retrieval and trend analysis via NCBI E-utilities.

Usage
-----
python pubmed_search.py --query "TP53 R175H pathogenicity" --max-results 30 --outdir results/

Environment
-----------
NCBI_EMAIL   — required by NCBI policy
NCBI_API_KEY — optional; raises rate limit from 3 → 10 req/s
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import urllib.parse
from datetime import datetime
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from dotenv import load_dotenv
from tenacity import retry, stop_after_attempt, wait_exponential

load_dotenv()

CACHE_DIR   = Path.home() / ".cache" / "genomics-skills" / "pubmed-search"
ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
EFETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EPMC_URL    = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
NCBI_EMAIL  = os.getenv("NCBI_EMAIL", "user@example.com")
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="PubMed literature search + trend analysis")
    p.add_argument("--query", required=True)
    p.add_argument("--max-results", type=int, default=50)
    p.add_argument("--date-from", default="2018")
    p.add_argument("--date-to", default=str(datetime.now().year))
    p.add_argument("--outdir", default="results/")
    p.add_argument("--dpi", type=int, default=300)
    return p


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=4))
def _esearch(query: str, max_results: int, date_from: str, date_to: str) -> list[str]:
    """Return a list of PMIDs."""
    params = {
        "db": "pubmed", "term": query,
        "retmax": max_results, "retmode": "json",
        "datetype": "pdat", "mindate": date_from, "maxdate": date_to,
        "email": NCBI_EMAIL,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    resp = requests.get(ESEARCH_URL, params=params, timeout=15)
    resp.raise_for_status()
    return resp.json()["esearchresult"].get("idlist", [])


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=4))
def _efetch(pmids: list[str]) -> list[dict]:
    """Fetch abstracts + metadata for a list of PMIDs."""
    if not pmids:
        return []
    params = {
        "db": "pubmed", "id": ",".join(pmids),
        "retmode": "xml", "email": NCBI_EMAIL,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    resp = requests.get(EFETCH_URL, params=params, timeout=30)
    resp.raise_for_status()

    root = ET.fromstring(resp.text)
    records = []
    for article in root.findall(".//PubmedArticle"):
        try:
            pmid  = article.findtext(".//PMID", "")
            title = article.findtext(".//ArticleTitle", "")
            abstract = " ".join(
                t.text or "" for t in article.findall(".//AbstractText")
            )
            year = article.findtext(".//PubDate/Year") or \
                   article.findtext(".//PubDate/MedlineDate", "")[:4]
            journal = article.findtext(".//Journal/Title", "")
            authors = ", ".join(
                f"{a.findtext('LastName', '')} {a.findtext('Initials', '')}".strip()
                for a in article.findall(".//Author")[:3]
            )
            records.append({"pmid": pmid, "title": title, "abstract": abstract,
                             "year": year, "journal": journal, "authors": authors,
                             "citations": 0})
        except Exception:
            continue
    return records


def _fetch_citations(records: list[dict]) -> list[dict]:
    """Add citation counts from Europe PMC."""
    pmids = [r["pmid"] for r in records if r["pmid"]]
    if not pmids:
        return records
    try:
        query = " OR ".join(f"EXT_ID:{p}" for p in pmids[:20])
        resp = requests.get(EPMC_URL,
                            params={"query": query, "format": "json", "pageSize": 20},
                            timeout=15)
        resp.raise_for_status()
        cite_map = {
            r["id"]: r.get("citedByCount", 0)
            for r in resp.json().get("resultList", {}).get("result", [])
        }
        for rec in records:
            rec["citations"] = cite_map.get(rec["pmid"], 0)
    except Exception as e:
        print(f"[warn] Europe PMC citation enrichment failed: {e}", file=sys.stderr)
    return records


def _plot_trend(df: pd.DataFrame, outdir: Path, dpi: int) -> None:
    year_counts = df["year"].value_counts().sort_index()
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(year_counts.index.astype(str), year_counts.values,
           color="#3498db", edgecolor="white")
    ax.set_xlabel("Year", fontsize=11)
    ax.set_ylabel("Publications", fontsize=11)
    ax.set_title("PubMed publications per year", fontsize=12, fontweight="bold")
    plt.xticks(rotation=45)
    plt.tight_layout()
    out = outdir / "pubmed_trend.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"[output] {out}")
    plt.close(fig)


def _write_digest(df: pd.DataFrame, query: str, outdir: Path) -> None:
    top = df.head(5)
    lines = [f"# PubMed digest: `{query}`\n",
             f"*{len(df)} papers retrieved · Generated {datetime.now().date()}*\n"]
    for i, (_, row) in enumerate(top.iterrows(), 1):
        lines += [
            f"\n## {i}. {row['title']}",
            f"**Authors:** {row['authors']}  ",
            f"**Journal:** {row['journal']} ({row['year']})  ",
            f"**PMID:** [{row['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{row['pmid']}/)  ",
            f"**Citations:** {row['citations']}\n",
            f"> {row['abstract'][:500]}{'...' if len(str(row['abstract'])) > 500 else ''}\n",
        ]
    digest_path = outdir / "pubmed_digest.md"
    digest_path.write_text("\n".join(lines))
    print(f"[output] {digest_path}")


def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cache_key = hashlib.md5(f"{args.query}{args.max_results}{args.date_from}{args.date_to}".encode()).hexdigest()
    cache_path = CACHE_DIR / f"{cache_key}.json"

    if cache_path.exists():
        print("[cache] Loading results from cache")
        records = json.loads(cache_path.read_text())
    else:
        print(f"[fetch] PubMed search: {args.query!r} ...")
        pmids = _esearch(args.query, args.max_results, args.date_from, args.date_to)
        if not pmids:
            print("[warn] No results found for query", file=sys.stderr)
            sys.exit(0)
        print(f"[info] {len(pmids)} PMIDs found. Fetching metadata ...")
        records = _efetch(pmids)
        records = _fetch_citations(records)
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(json.dumps(records))

    df = pd.DataFrame(records).sort_values("citations", ascending=False)
    tsv = outdir / "pubmed_results.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"[output] {tsv} ({len(df)} papers)")

    _plot_trend(df, outdir, args.dpi)
    _write_digest(df, args.query, outdir)
    print(f"\n✓ pubmed-search complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
