# SKILL: protein-variant-mapper

## Purpose
Given a gene symbol and a list of amino acid variants (e.g. `R175H`, `G245S`),
fetch the UniProt entry and AlphaFold / PDB structure, then produce:
- A lollipop map PNG with each variant labeled at its position and colored by
  predicted or known pathogenicity
- An interactive 3-D HTML viewer with all variants displayed as labeled spheres
  on the protein backbone

Designed as a direct downstream step from `agentic-genomics` variant interpretation:
take the top-ranked missense variant(s), drop them onto the protein structure for
structural context.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--gene` | str | ✅ | — | HGNC gene symbol (e.g. `TP53`) |
| `--variants` | str | ✅ | — | Comma-separated variant list in HGVS-p short form (e.g. `R175H,G245S`) |
| `--pdb-id` | str | ❌ | auto | PDB or AlphaFold ID to use; auto-fetched from UniProt if omitted |
| `--outdir` | str | ❌ | `results/` | Directory where outputs are written |
| `--dpi` | int | ❌ | `300` | Lollipop plot resolution |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `{gene}_variant_map.png` | PNG (300 DPI) | Lollipop with domain track + variant labels |
| `{gene}_variant_map.svg` | SVG | Vector lollipop |
| `{gene}_3d_viewer.html` | HTML | Interactive py3Dmol viewer with labeled variant spheres |
| `{gene}_variant_info.tsv` | TSV | Variant list with position, UniProt annotation, AlphaFold model source |

## Execution policy

- UniProt REST API used for protein length and domain annotations (free, no auth)
- AlphaFold DB used for structure (`https://alphafold.ebi.ac.uk/`) when no PDB ID given
- Caches structure files to `~/.cache/genomics-skills/protein-variant-mapper/`
- 3-D HTML viewer is self-contained (no internet required to open)
- Fails gracefully if AlphaFold has no entry for the gene (falls back to annotation-only lollipop)

## Trigger phrases (for agent routing)

- "map variant {X} onto the structure of {gene}"
- "show {gene} {variant} on 3-D structure"
- "where does {variant} land on {gene} protein"
- "visualize {gene} missense variants structurally"
- "protein variant map for {gene}"
- "lollipop with variants {variant list} for {gene}"

## Example invocation

```bash
python skills/protein-variant-mapper/scripts/protein_variant_mapper.py \
    --gene TP53 --variants R175H,G245S,R248W --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.
