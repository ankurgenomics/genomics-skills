#!/usr/bin/env python3
# Author: Ankur Sharma | Agentic AI · Machine Learning · Bioinformatics · Data Science
# GitHub: https://github.com/ankurgenomics | Portfolio: agentic-genomics
"""protein-structure-viewer — 3-D HTML viewer + pocket detection + B-factor plot.

Fetches AlphaFold or PDB structure, generates a self-contained HTML viewer,
a per-residue B-factor/pLDDT line plot, and a simple pocket TSV.

Usage
-----
python protein_structure_viewer.py --gene EGFR --outdir results/
python protein_structure_viewer.py --pdb-id 2GS6 --outdir results/
"""
from __future__ import annotations

import argparse
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

UNIPROT_SEARCH  = "https://rest.uniprot.org/uniprotkb/search"
ALPHAFOLD_PDB   = "https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
PDB_RCSB        = "https://files.rcsb.org/download/{pdb_id}.pdb"
CACHE_DIR = Path.home() / ".cache" / "genomics-skills" / "protein-structure-viewer"


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="3-D protein structure viewer + pocket detection")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene",   help="HGNC gene symbol (AlphaFold auto-fetched)")
    group.add_argument("--pdb-id", help="PDB 4-letter accession")
    p.add_argument("--outdir", default="results/")
    p.add_argument("--dpi",    type=int, default=300)
    return p


# ── structure fetch ───────────────────────────────────────────────────────────

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
def _uniprot_id(gene: str) -> str:
    resp = requests.get(UNIPROT_SEARCH,
                        params={"query": f"gene:{gene} AND organism_id:9606 AND reviewed:true",
                                "fields": "accession", "format": "json", "size": 1},
                        timeout=15)
    resp.raise_for_status()
    results = resp.json().get("results", [])
    return results[0]["primaryAccession"] if results else ""


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=1, max=8))
def _fetch_pdb(url: str, cache_path: Path) -> str | None:
    if cache_path.exists():
        return cache_path.read_text()
    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(resp.text)
    return resp.text


def _get_pdb_text(gene: str | None, pdb_id: str | None) -> tuple[str, str, str | None]:
    """Returns (label, source_desc, pdb_text_or_None)."""
    if pdb_id:
        url = PDB_RCSB.format(pdb_id=pdb_id)
        cache = CACHE_DIR / f"pdb_{pdb_id}.pdb"
        label = pdb_id
        source = f"PDB {pdb_id}"
    else:
        uid = _uniprot_id(gene)
        if not uid:
            print(f"[warn] No UniProt entry for {gene}", file=sys.stderr)
            return gene, "AlphaFold (not found)", None
        url = ALPHAFOLD_PDB.format(uid=uid)
        cache = CACHE_DIR / f"AF-{uid}.pdb"
        label = gene
        source = f"AlphaFold · UniProt {uid}"

    print(f"[fetch] {source} ...")
    pdb_text = _fetch_pdb(url, cache)
    return label, source, pdb_text


# ── PDB parsing (no biopython needed for basic ops) ──────────────────────────

def _parse_ca_atoms(pdb_text: str) -> pd.DataFrame:
    """Extract Cα ATOM records: residue number, x, y, z, b_factor."""
    rows = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                rows.append({
                    "resnum": int(line[22:26].strip()),
                    "resname": line[17:20].strip(),
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "b_factor": float(line[60:66]),
                })
            except ValueError as e:
                print(f"[warn] Skipping malformed PDB ATOM record: {e}", file=sys.stderr)
    return pd.DataFrame(rows)


# ── pocket detection (simple geometric cavity search on Cα grid) ──────────────

def _detect_pockets(ca_df: pd.DataFrame, grid_spacing: float = 2.0,
                    min_pocket_size: int = 8) -> pd.DataFrame:
    """
    Very simple cavity detector: voxelizes Cα positions, finds grid points that
    are surrounded by protein atoms on multiple sides. Returns pocket centroids.
    """
    if ca_df.empty:
        return pd.DataFrame(columns=["pocket_id", "centroid_x", "centroid_y",
                                     "centroid_z", "residues", "avg_bfactor"])
    coords = ca_df[["x", "y", "z"]].values
    pockets = []
    # Sample grid around protein
    mins = coords.min(axis=0) - 5
    maxs = coords.max(axis=0) + 5
    gx = np.arange(mins[0], maxs[0], grid_spacing)
    gy = np.arange(mins[1], maxs[1], grid_spacing)
    gz = np.arange(mins[2], maxs[2], grid_spacing)

    for pocket_id, (cx, cy, cz) in enumerate(
        [(x, y, z) for x in gx[::3] for y in gy[::3] for z in gz[::3]]
    ):
        dists = np.sqrt(((coords - np.array([cx, cy, cz])) ** 2).sum(axis=1))
        nearby_mask = dists < 8.0
        if nearby_mask.sum() >= min_pocket_size:
            nearby = ca_df[nearby_mask]
            pockets.append({
                "pocket_id": len(pockets) + 1,
                "centroid_x": round(cx, 2),
                "centroid_y": round(cy, 2),
                "centroid_z": round(cz, 2),
                "residues": len(nearby),
                "avg_bfactor": round(nearby["b_factor"].mean(), 2),
            })
        if len(pockets) >= 10:
            break

    return pd.DataFrame(pockets)


# ── B-factor plot ─────────────────────────────────────────────────────────────

def _plot_bfactor(ca_df: pd.DataFrame, label: str, outdir: Path, dpi: int) -> None:
    if ca_df.empty:
        return
    fig, ax = plt.subplots(figsize=(12, 3))
    ax.plot(ca_df["resnum"], ca_df["b_factor"], linewidth=0.8, color="#3498db")
    ax.fill_between(ca_df["resnum"], 0, ca_df["b_factor"], alpha=0.25, color="#3498db")
    ax.axhline(70, color="#2ecc71", linestyle="--", linewidth=0.8, label="pLDDT > 70 (confident)")
    ax.set_xlabel("Residue number", fontsize=11)
    ax.set_ylabel("B-factor / pLDDT", fontsize=11)
    ax.set_title(f"{label} — per-residue B-factor / AlphaFold pLDDT", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = outdir / f"{label}_bfactor.png"
    fig.savefig(out, dpi=dpi, bbox_inches="tight")
    print(f"[output] {out}")
    plt.close(fig)


# ── 3-D HTML viewer ───────────────────────────────────────────────────────────

def _build_viewer_html(label: str, source: str, pdb_text: str, outdir: Path) -> None:
    pdb_escaped = pdb_text.replace("`", "\\`").replace("\\", "\\\\")
    html = textwrap.dedent(f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>{label} — 3-D Structure Viewer</title>
        <script src="https://3dmol.org/build/3Dmol-min.js"></script>
        <style>
            body  {{ margin: 0; background: #0d1117; color: #c9d1d9; font-family: sans-serif; }}
            h2    {{ padding: 10px 20px; margin: 0; background: #161b22;
                     border-bottom: 1px solid #30363d; }}
            #view {{ width: 100vw; height: 90vh; position: relative; }}
            #info {{ padding: 8px 20px; background: #161b22; font-size: 12px;
                     border-top: 1px solid #30363d; }}
        </style>
    </head>
    <body>
        <h2>{label} · {source}</h2>
        <div id="view"></div>
        <div id="info">
            Backbone colored by B-factor / pLDDT (blue=high confidence, red=low).
            Scroll to zoom · Click+drag to rotate · Right-click to pan.
        </div>
        <script>
        let viewer = $3Dmol.createViewer(document.getElementById("view"), {{
            backgroundColor: "0x0d1117"
        }});
        let pdb = `{pdb_escaped}`;
        viewer.addModel(pdb, "pdb");
        viewer.setStyle({{}}, {{cartoon: {{colorscheme: "bfactor", style: "oval"}}}});
        viewer.addSurface($3Dmol.SurfaceType.SAS, {{
            opacity: 0.15, colorscheme: "whiteCarbon"
        }});
        viewer.zoomTo();
        viewer.render();
        </script>
    </body>
    </html>
    """).strip()
    out = outdir / f"{label}_3d_viewer.html"
    out.write_text(html)
    print(f"[output] {out}")


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    args = _build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    label, source, pdb_text = _get_pdb_text(args.gene, args.pdb_id)

    if not pdb_text:
        print(f"[error] Could not fetch structure for {label}.", file=sys.stderr)
        sys.exit(1)

    ca_df = _parse_ca_atoms(pdb_text)
    print(f"[info] {label}: {len(ca_df)} Cα atoms parsed")

    _plot_bfactor(ca_df, label, outdir, args.dpi)

    pockets = _detect_pockets(ca_df)
    pocket_tsv = outdir / f"{label}_pockets.tsv"
    pockets.to_csv(pocket_tsv, sep="\t", index=False)
    print(f"[output] {pocket_tsv} ({len(pockets)} pocket candidates)")

    _build_viewer_html(label, source, pdb_text, outdir)

    print(f"\n✓ protein-structure-viewer complete. Outputs in {outdir}/")


if __name__ == "__main__":
    main()
