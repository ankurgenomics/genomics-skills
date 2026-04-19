#!/usr/bin/env python3
# Author: Ankur Sharma | Agentic AI · Machine Learning · Bioinformatics · Data Science
# GitHub: https://github.com/ankurgenomics | Portfolio: agentic-genomics
"""protein-variant-mapper — map missense variants onto a 3-D protein structure.

Fetches UniProt domain annotations + AlphaFold structure, then produces:
  - Lollipop PNG/SVG with variant labels
  - Self-contained 3-D HTML viewer with colored variant spheres

Usage
-----
python protein_variant_mapper.py --gene TP53 --variants R175H,G245S --outdir results/
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import textwrap
from pathlib import Path
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

# ── constants ────────────────────────────────────────────────────────────────

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_PDB  = "https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
CACHE_DIR = Path.home() / ".cache" / "genomics-skills" / "protein-variant-mapper"

PATHOGENICITY_COLORS = {
    "pathogenic": "#c0392b",
    "likely_pathogenic": "#e74c3c",
    "uncertain": "#f39c12",
    "likely_benign": "#27ae60",
    "benign": "#2ecc71",
    "unknown": "#95a5a6",
}

# ── CLI ──────────────────────────────────────────────────────────────────────


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Map missense variants onto a protein structure")
    p.add_argument("--gene", required=True)
    p.add_argument("--variants", required=True,
                   help="Comma-separated HGVS-p short form, e.g. R175H,G245S")
    p.add_argument("--pdb-id", default=None,
                   help="PDB or AlphaFold accession; auto-fetched if omitted")
    p.add_argument("--outdir", default="results/")
    p.add_argument("--dpi", type=int, default=300)
    return p


def _parse_variants(raw: str) -> list[dict[str, Any]]:
    """Parse 'R175H,G245S' → [{'ref':'R','pos':175,'alt':'H'}, ...]"""
    variants = []
    for v in raw.split(","):
        v = v.strip()
        m = re.match(r"([A-Z])(\d+)([A-Z*])", v)
        if m:
            variants.append({"label": v, "ref": m.group(1),
                              "pos": int(m.group(2)), "alt": m.group(3),
                              "pathogenicity": "unknown"})
        else:
            print(f"[warn] Could not parse variant '{v}' — skipping", file=sys.stderr)
    return variants


# ── UniProt ───────────────────────────────────────────────────────────────────


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
def _fetch_uniprot(gene: str) -> dict[str, Any]:
    cache_path = CACHE_DIR / f"{gene}_uniprot.json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())

    print(f"[fetch] UniProt for {gene} ...")
    resp = requests.get(
        UNIPROT_SEARCH,
        params={"query": f"gene:{gene} AND organism_id:9606 AND reviewed:true",
                "fields": "accession,sequence,ft_domain,ft_region,length",
                "format": "json", "size": 1},
        timeout=15,
    )
    resp.raise_for_status()
    results = resp.json().get("results", [])
    if not results:
        return {}
    data = results[0]
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(data))
    return data


def _extract_domains(uniprot_data: dict) -> list[dict]:
    """Pull domain annotations from UniProt features."""
    domains = []
    for feat in uniprot_data.get("features", []):
        if feat.get("type") in ("Domain", "Region", "Motif"):
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value", 0)
            end = loc.get("end", {}).get("value", 0)
            domains.append({
                "name": feat.get("description", feat["type"]),
                "start": start,
                "end": end,
            })
    return domains


# ── AlphaFold structure ───────────────────────────────────────────────────────


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
def _fetch_alphafold_pdb(uniprot_id: str) -> str | None:
    """Return the PDB text for an AlphaFold model, or None if unavailable."""
    cache_path = CACHE_DIR / f"AF-{uniprot_id}.pdb"
    if cache_path.exists():
        return cache_path.read_text()

    url = ALPHAFOLD_PDB.format(uniprot_id=uniprot_id)
    print(f"[fetch] AlphaFold structure {url} ...")
    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        print(f"[warn] No AlphaFold entry for {uniprot_id}")
        return None
    resp.raise_for_status()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(resp.text)
    return resp.text


# ── lollipop plot ─────────────────────────────────────────────────────────────


def _plot_lollipop(
    variants: list[dict],
    domains: list[dict],
    protein_length: int,
    gene: str,
    outdir: Path,
    dpi: int,
) -> None:
    fig, ax = plt.subplots(figsize=(14, 4))
    ax.set_xlim(0, protein_length + 20)
    ax.set_ylim(-1.5, 4)

    # Protein backbone
    ax.fill_between([0, protein_length], [-0.15, -0.15], [0.15, 0.15],
                    color="#bdc3c7", zorder=1)

    # Domain boxes
    domain_colors = plt.cm.Set3(np.linspace(0, 1, max(len(domains), 1)))
    for i, dom in enumerate(domains):
        ax.fill_between(
            [dom["start"], dom["end"]], [-0.3, -0.3], [0.3, 0.3],
            color=domain_colors[i], alpha=0.7, zorder=2, label=dom["name"]
        )

    # Variant lollipops
    for v in variants:
        pos = v["pos"]
        color = PATHOGENICITY_COLORS.get(v["pathogenicity"], PATHOGENICITY_COLORS["unknown"])
        ax.plot([pos, pos], [0, 2.5], color=color, linewidth=1.2, zorder=3)
        ax.scatter(pos, 2.5, color=color, s=80, zorder=4)
        ax.text(pos, 2.7, v["label"], ha="center", va="bottom",
                fontsize=8, rotation=45, color=color, fontweight="bold")

    ax.set_xlabel("Amino acid position", fontsize=11)
    ax.set_yticks([])
    ax.set_title(f"{gene} — Variant map ({protein_length} aa)", fontsize=13, fontweight="bold")

    if domains:
        ax.legend(loc="lower right", fontsize=7, framealpha=0.8)

    plt.tight_layout()
    for ext in ("png", "svg"):
        out = outdir / f"{gene}_variant_map.{ext}"
        fig.savefig(out, dpi=dpi if ext == "png" else None, bbox_inches="tight")
        print(f"[output] {out}")
    plt.close(fig)


# ── 3-D HTML viewer ───────────────────────────────────────────────────────────


def _build_3d_html(
    gene: str,
    pdb_text: str | None,
    variants: list[dict],
    uniprot_id: str,
    outdir: Path,
) -> None:
    """Generate a self-contained HTML with py3Dmol embedded via CDN."""
    if not pdb_text:
        print("[warn] No PDB/AlphaFold structure available — skipping 3-D viewer")
        return

    # Escape PDB for embedding
    pdb_escaped = pdb_text.replace("`", "\\`").replace("\\", "\\\\")

    # Build JavaScript for each variant sphere
    sphere_js_lines = []
    color_map = {
        "pathogenic": "0xc0392b",
        "likely_pathogenic": "0xe74c3c",
        "uncertain": "0xf39c12",
        "unknown": "0x95a5a6",
    }
    for v in variants:
        color = color_map.get(v["pathogenicity"], "0x95a5a6")
        sphere_js_lines.append(
            f"viewer.addResLabels({{resi: {v['pos']}, label: '{v['label']}', "
            f"backgroundColor: {color}, fontColor: 0xffffff, fontSize: 12}});"
        )
        sphere_js_lines.append(
            f"viewer.addSphere({{center: {{seq: {v['pos']}}}, radius: 1.2, "
            f"color: {color}, opacity: 0.85}});"
        )
    sphere_js = "\n".join(sphere_js_lines)

    html = textwrap.dedent(f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{gene} — Variant 3-D viewer</title>
        <script src="https://3dmol.org/build/3Dmol-min.js"></script>
        <style>
            body {{ font-family: sans-serif; margin: 0; background: #1a1a2e; color: #eee; }}
            h2   {{ padding: 12px 20px; margin: 0; background: #16213e; }}
            #viewer {{ width: 100vw; height: 85vh; position: relative; }}
            #legend {{ padding: 10px 20px; background: #0f3460; font-size: 13px; }}
        </style>
    </head>
    <body>
        <h2>{gene} — {len(variants)} variant(s) mapped · UniProt: {uniprot_id}</h2>
        <div id="viewer"></div>
        <div id="legend">
            Variants: {', '.join(v['label'] for v in variants)} &nbsp;|&nbsp;
            Backbone colored by B-factor (AlphaFold confidence). Variant spheres shown in red/orange.
        </div>
        <script>
            let viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{
                backgroundColor: "0x1a1a2e"
            }});
            let pdbStr = `{pdb_escaped}`;
            viewer.addModel(pdbStr, "pdb");
            viewer.setStyle({{}}, {{cartoon: {{colorscheme: "bfactor", style: "oval"}}}});
            {sphere_js}
            viewer.zoomTo();
            viewer.render();
        </script>
    </body>
    </html>
    """).strip()

    out = outdir / f"{gene}_3d_viewer.html"
    out.write_text(html)
    print(f"[output] {out}")


# ── main ─────────────────────────────────────────────────────────────────────


def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    variants = _parse_variants(args.variants)
    if not variants:
        print("[error] No valid variants parsed. Example: --variants R175H,G245S", file=sys.stderr)
        sys.exit(1)

    # UniProt
    uniprot_data = _fetch_uniprot(args.gene)
    uniprot_id = uniprot_data.get("primaryAccession", "")
    protein_length = uniprot_data.get("sequence", {}).get("length", 500)
    domains = _extract_domains(uniprot_data) if uniprot_data else []
    print(f"[info] {args.gene} → UniProt {uniprot_id} ({protein_length} aa, {len(domains)} domains)")

    # AlphaFold structure
    pdb_text = None
    if uniprot_id:
        pdb_text = _fetch_alphafold_pdb(uniprot_id)

    # Plots
    _plot_lollipop(variants, domains, protein_length, args.gene, outdir, args.dpi)
    _build_3d_html(args.gene, pdb_text, variants, uniprot_id, outdir)

    # TSV
    df = pd.DataFrame(variants)
    df["gene"] = args.gene
    df["uniprot_id"] = uniprot_id
    df["protein_length"] = protein_length
    tsv = outdir / f"{args.gene}_variant_info.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"[output] {tsv}")

    print(f"\n✓ protein-variant-mapper complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
