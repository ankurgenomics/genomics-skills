#!/usr/bin/env python3
"""survival-analysis — Kaplan–Meier + Cox PH from TCGA gene expression.

Fetches real expression + clinical survival data from cBioPortal,
splits samples by high/low expression, runs KM + log-rank + Cox regression.

Usage
-----
python survival_analysis.py --gene TP53 --cancer TCGA-LUAD --endpoint os
python survival_analysis.py --gene EGFR --cancer TCGA-LUAD --covariates age,stage
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
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

CBIO_BASE = "https://www.cbioportal.org/api"
CACHE_DIR = Path.home() / ".cache" / "genomics-skills" / "survival-analysis"
MIN_EVENTS = 5


# ── CLI ──────────────────────────────────────────────────────────────────────

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="KM + Cox survival analysis from TCGA")
    p.add_argument("--gene", required=True)
    p.add_argument("--cancer", required=True, help="TCGA project ID, e.g. TCGA-LUAD")
    p.add_argument("--endpoint", default="os", choices=["os", "dfs"])
    p.add_argument("--split", default="median", choices=["median", "tertile", "quartile"])
    p.add_argument("--covariates", default=None,
                   help="Comma-separated clinical covariates for Cox (e.g. age,stage)")
    p.add_argument("--outdir", default="results/")
    p.add_argument("--dpi", type=int, default=300)
    return p


# ── real data fetch from cBioPortal ──────────────────────────────────────────

def _study_id(cancer: str) -> str:
    """Convert TCGA-LUAD → luad_tcga_pan_can_atlas_2018."""
    prefix = cancer.lower().replace("tcga-", "")
    return f"{prefix}_tcga_pan_can_atlas_2018"


def _fetch_survival_data(gene: str, cancer: str, endpoint: str) -> pd.DataFrame:
    """Fetch real expression + survival from cBioPortal. Returns merged DataFrame."""
    study = _study_id(cancer)
    cache_path = CACHE_DIR / f"{gene}_{study}_{endpoint}.parquet"
    if cache_path.exists():
        print(f"[cache] {cache_path}")
        return pd.read_parquet(cache_path)

    print(f"[fetch] cBioPortal: {gene} expression + {endpoint.upper()} survival in {cancer} ...")

    # 1. Entrez ID
    r = requests.get(f"{CBIO_BASE}/genes/{gene}", timeout=10)
    r.raise_for_status()
    entrez_id = r.json()["entrezGeneId"]

    # 2. Expression per sample
    profile = f"{study}_rna_seq_v2_mrna"
    sample_list = f"{study}_all"
    r = requests.get(
        f"{CBIO_BASE}/molecular-profiles/{profile}/molecular-data",
        params={"entrezGeneId": entrez_id, "sampleListId": sample_list},
        headers={"Accept": "application/json"}, timeout=20,
    )
    if r.status_code == 404:
        print(f"[warn] No expression profile found for {study}", file=sys.stderr)
        return pd.DataFrame()
    r.raise_for_status()
    expr_df = pd.DataFrame([
        {"patientId": rec["patientId"], "expression": np.log2(max(float(rec["value"]), 0) + 1)}
        for rec in r.json() if rec.get("value") is not None
    ])
    print(f"  Expression: {len(expr_df)} samples")

    # 3. Clinical survival data (patient level)
    r = requests.get(
        f"{CBIO_BASE}/studies/{study}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
        headers={"Accept": "application/json"}, timeout=20,
    )
    r.raise_for_status()
    clin_records = r.json()

    # Pivot clinical attributes into one row per patient
    time_col   = "OS_MONTHS"  if endpoint == "os"  else "DFS_MONTHS"
    status_col = "OS_STATUS"  if endpoint == "os"  else "DFS_STATUS"
    wanted = {time_col, status_col, "AGE"}
    pivot: dict[str, dict] = {}
    for rec in clin_records:
        pid = rec["patientId"]
        attr = rec["clinicalAttributeId"]
        if attr in wanted:
            pivot.setdefault(pid, {})[attr] = rec["value"]

    clin_df = pd.DataFrame.from_dict(pivot, orient="index").reset_index()
    clin_df.columns = ["patientId"] + list(clin_df.columns[1:])

    if time_col not in clin_df.columns or status_col not in clin_df.columns:
        print(f"[warn] Survival columns {time_col}/{status_col} not found", file=sys.stderr)
        return pd.DataFrame()

    clin_df["time"] = pd.to_numeric(clin_df[time_col], errors="coerce") * 30.4  # months→days
    clin_df["event"] = clin_df[status_col].str.contains("DECEASED|Recurred", na=False).astype(int)
    if "AGE" in clin_df.columns:
        clin_df["age"] = pd.to_numeric(clin_df["AGE"], errors="coerce")
    print(f"  Clinical: {len(clin_df)} patients")

    # 4. Merge on patientId
    df = expr_df.merge(clin_df[["patientId", "time", "event"] +
                                (["age"] if "age" in clin_df.columns else [])],
                       on="patientId", how="inner").dropna(subset=["time", "event", "expression"])
    df["cancer"] = cancer
    df["gene"]   = gene
    df["endpoint"] = endpoint
    print(f"  Merged: {len(df)} patients with both expression and survival")

    if not df.empty:
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        df.to_parquet(cache_path, index=False)
    return df


# ── KM plotting ───────────────────────────────────────────────────────────────

def _plot_km(df: pd.DataFrame, gene: str, cancer: str,
             endpoint: str, outdir: Path, dpi: int) -> None:
    high = df[df["group"] == "High"]
    low  = df[df["group"] == "Low"]

    # Minimum events check
    if high["event"].sum() < MIN_EVENTS or low["event"].sum() < MIN_EVENTS:
        print(f"[warn] Fewer than {MIN_EVENTS} events in one group — KM plot may be unreliable")

    # Log-rank test
    lr = logrank_test(high["time"], low["time"],
                      event_observed_A=high["event"], event_observed_B=low["event"])
    p_val = lr.p_value

    fig, (ax_km, ax_risk) = plt.subplots(
        2, 1, figsize=(9, 7), gridspec_kw={"height_ratios": [4, 1]}, sharex=True
    )

    kmf_high = KaplanMeierFitter(label=f"High {gene} (n={len(high)})")
    kmf_low  = KaplanMeierFitter(label=f"Low {gene} (n={len(low)})")
    kmf_high.fit(high["time"], event_observed=high["event"])
    kmf_low.fit(low["time"],  event_observed=low["event"])

    kmf_high.plot_survival_function(ax=ax_km, color="#c0392b", ci_show=True, ci_alpha=0.15)
    kmf_low.plot_survival_function( ax=ax_km, color="#2980b9", ci_show=True, ci_alpha=0.15)

    # p-value annotation
    p_str = f"Log-rank p = {p_val:.4f}" if p_val >= 0.0001 else f"Log-rank p < 0.0001"
    ax_km.text(0.97, 0.97, p_str, transform=ax_km.transAxes,
               ha="right", va="top", fontsize=10,
               bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "gray", "alpha": 0.8})

    endpoint_label = "Overall Survival" if endpoint == "os" else "Disease-Free Survival"
    ax_km.set_title(f"{gene} — {endpoint_label} in {cancer}", fontsize=13, fontweight="bold")
    ax_km.set_ylabel("Survival probability", fontsize=11)
    ax_km.set_ylim(0, 1.05)
    ax_km.legend(loc="lower left", fontsize=10)

    # At-risk table
    time_points = np.linspace(0, df["time"].max(), 6).astype(int)
    for i, (kmf, color, label) in enumerate([
        (kmf_high, "#c0392b", "High"),
        (kmf_low,  "#2980b9", "Low"),
    ]):
        at_risk = [kmf.event_table["at_risk"].asof(t) for t in time_points]
        for j, (t, n) in enumerate(zip(time_points, at_risk)):
            ax_risk.text(t, i, str(int(n)), ha="center", va="center",
                         fontsize=8, color=color)
        ax_risk.text(-df["time"].max() * 0.02, i, label, ha="right", va="center",
                     fontsize=9, color=color, fontweight="bold")

    ax_risk.set_yticks([])
    ax_risk.set_xlabel("Time (days)", fontsize=11)
    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_title("At risk", fontsize=9, loc="left")

    plt.tight_layout()
    prefix = f"{gene}_{cancer.replace('-', '_')}_km_{endpoint}"
    for ext in ("png", "svg"):
        out = outdir / f"{prefix}.{ext}"
        fig.savefig(out, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[output] {out}")
    plt.close(fig)


# ── Cox regression ────────────────────────────────────────────────────────────

def _run_cox(df: pd.DataFrame, covariates: list[str],
             gene: str, cancer: str, outdir: Path, dpi: int) -> None:
    cols = ["time", "event", "expression"] + [c for c in covariates if c in df.columns]
    cox_df = df[cols].dropna()
    if len(cox_df) < 20:
        print("[warn] Too few samples for Cox regression — skipping")
        return

    cph = CoxPHFitter()
    try:
        cph.fit(cox_df, duration_col="time", event_col="event")
    except Exception as exc:
        print(f"[warn] Cox model failed: {exc}", file=sys.stderr)
        return

    # Forest plot
    summary = cph.summary.reset_index()
    summary.columns = [c.replace(" ", "_") for c in summary.columns]
    fig, ax = plt.subplots(figsize=(7, max(3, len(summary))))
    for i, row in summary.iterrows():
        coef = row.get("exp(coef)", row.get("coef", 1.0))
        lo   = row.get("exp(coef)_lower_0.95", coef * 0.7)
        hi   = row.get("exp(coef)_upper_0.95", coef * 1.3)
        pval = row.get("p", 1.0)
        color = "#c0392b" if coef > 1 else "#2980b9"
        ax.errorbar(coef, i, xerr=[[coef - lo], [hi - coef]],
                    fmt="o", color=color, markersize=8, capsize=5, linewidth=1.5)
        ax.text(hi + 0.05, i, f"HR={coef:.2f}, p={pval:.3f}",
                va="center", fontsize=9)

    ax.axvline(1.0, color="gray", linestyle="--", linewidth=1)
    ax.set_yticks(range(len(summary)))
    ax.set_yticklabels(summary.iloc[:, 0].tolist(), fontsize=10)
    ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=11)
    ax.set_title(f"{gene} — Cox PH model ({cancer})", fontsize=12, fontweight="bold")
    plt.tight_layout()

    out = outdir / f"{gene}_{cancer.replace('-', '_')}_cox.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"[output] {out}")
    plt.close(fig)

    tsv = outdir / f"{gene}_{cancer.replace('-', '_')}_cox_results.tsv"
    cph.summary.to_csv(tsv, sep="\t")
    print(f"[output] {tsv}")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = _fetch_survival_data(args.gene, args.cancer, args.endpoint)
    if df.empty:
        print("[error] No survival data retrieved.", file=sys.stderr)
        sys.exit(1)

    # Split by expression
    if args.split == "median":
        threshold = df["expression"].median()
        df["group"] = np.where(df["expression"] >= threshold, "High", "Low")
    elif args.split == "tertile":
        lo, hi = df["expression"].quantile([1/3, 2/3])
        df = df[df["expression"] <= lo | (df["expression"] >= hi)].copy()
        df["group"] = np.where(df["expression"] >= hi, "High", "Low")
    else:  # quartile
        lo, hi = df["expression"].quantile([0.25, 0.75])
        df = df[(df["expression"] <= lo) | (df["expression"] >= hi)].copy()
        df["group"] = np.where(df["expression"] >= hi, "High", "Low")

    print(f"[info] {args.gene} in {args.cancer} — "
          f"High: {(df['group']=='High').sum()}, Low: {(df['group']=='Low').sum()}")

    # Save TSV
    tsv = outdir / f"{args.gene}_{args.cancer.replace('-', '_')}_survival.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"[output] {tsv}")

    # KM plot
    _plot_km(df, args.gene, args.cancer, args.endpoint, outdir, args.dpi)

    # Cox (optional)
    if args.covariates:
        covariates = [c.strip() for c in args.covariates.split(",")]
        _run_cox(df, covariates, args.gene, args.cancer, outdir, args.dpi)

    print(f"\n✓ survival-analysis complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
