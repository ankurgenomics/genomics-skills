# SKILL: variant-context

## Purpose
Given a gene symbol and optional cancer type, fetch somatic mutation data from
MyVariant.info / COSMIC / cBioPortal and produce:
- A mutation lollipop plot (PNG + SVG) with domain annotations
- A somatic hotspot bar chart showing the top mutated positions
- A TSV of all mutations with consequence, frequency, and functional scores

This skill provides the somatic mutation landscape for a candidate gene —
complementing single-variant ACMG classification with population-level context.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--gene` | str | ✅ | — | HGNC gene symbol (e.g. `TP53`, `BRCA1`) |
| `--cancer` | str | ❌ | `all` | TCGA cancer type abbreviation (e.g. `LUAD`, `BRCA`) or `all` |
| `--assembly` | str | ❌ | `hg38` | Genome assembly (`hg38` or `hg19`) |
| `--outdir` | str | ❌ | `results/` | Directory where outputs are written |
| `--top-n` | int | ❌ | `20` | Number of top hotspot positions to highlight |
| `--dpi` | int | ❌ | `300` | Output plot resolution |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `{gene}_lollipop.png` | PNG (300 DPI) | Mutation lollipop with protein domain track |
| `{gene}_lollipop.svg` | SVG | Vector version of lollipop for publication |
| `{gene}_hotspots.png` | PNG | Bar chart of top N mutated positions |
| `{gene}_mutations.tsv` | TSV | All mutations: position, consequence, frequency, CADD, ClinVar |

## Execution policy

- Fetches from MyVariant.info REST API (free, no auth required)
- Caches responses to `~/.cache/genomics-skills/variant-context/`
- Retries up to 3× with exponential backoff on transient HTTP errors
- Fails gracefully if gene has no somatic mutation data (empty TSV + warning)
- Runtime: ~10–30 seconds depending on gene size and cache state

## Trigger phrases (for agent routing)

- "show mutation landscape for {gene}"
- "lollipop plot for {gene}"
- "somatic hotspots in {gene}"
- "mutation frequency in {gene} in {cancer}"
- "what positions are mutated in {gene}"
- "COSMIC / cBioPortal mutation map for {gene}"

## Example invocation

```bash
python skills/variant-context/scripts/variant_context.py \
    --gene TP53 --cancer LUAD --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.
