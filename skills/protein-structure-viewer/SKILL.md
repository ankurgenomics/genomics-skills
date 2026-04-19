# SKILL: protein-structure-viewer

## Purpose
Fetch a protein structure (PDB or AlphaFold) and produce an interactive 3-D
HTML viewer, a pocket detection TSV, and a B-factor / confidence plot — no
local PyMOL or ChimeraX installation required.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--gene` | str | ✅* | — | HGNC gene symbol; UniProt + AlphaFold ID auto-resolved |
| `--pdb-id` | str | ✅* | — | PDB 4-letter code (overrides `--gene` structure lookup) |
| `--outdir` | str | ❌ | `results/` | Output directory |
| `--dpi` | int | ❌ | `300` | B-factor plot resolution |

\* Provide `--gene` OR `--pdb-id` (at least one required).

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `{gene}_3d_viewer.html` | HTML | Self-contained 3-D viewer (py3Dmol via CDN) |
| `{gene}_bfactor.png` | PNG | B-factor / AlphaFold pLDDT per-residue line plot |
| `{gene}_pockets.tsv` | TSV | Pocket candidates: centroid residues, volume estimate, avg B-factor |

## Execution policy

- AlphaFold DB or PDB RCSB used for structure download (free, no auth)
- Pocket detection: geometric cavity search on Cα coordinates (no external binary needed)
- Self-contained HTML uses py3Dmol CDN — opens offline once downloaded
- Caches PDB files to `~/.cache/genomics-skills/protein-structure-viewer/`

## Trigger phrases (for agent routing)

- "show 3D structure of {gene}"
- "visualize {gene} protein structure"
- "open structure viewer for {gene}"
- "protein structure for {pdb_id}"
- "detect pockets in {gene}"
- "AlphaFold structure {gene}"

## Example invocation

```bash
python skills/protein-structure-viewer/scripts/protein_structure_viewer.py \
    --gene EGFR --outdir results/

python skills/protein-structure-viewer/scripts/protein_structure_viewer.py \
    --pdb-id 2GS6 --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.
