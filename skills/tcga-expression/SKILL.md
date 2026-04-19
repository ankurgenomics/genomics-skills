# SKILL: tcga-expression

## Purpose
Query the **cBioPortal REST API** to retrieve real TCGA Pan-Cancer Atlas 2018
RNA-seq v2 RSEM expression for a gene (9,479 samples · 31 cancer types), then produce:
- Pan-cancer expression box plots across all TCGA cohorts (PNG + SVG)
- Tumor-vs-normal comparison for a specific cohort (if requested)
- Expression TSV for downstream statistical analysis or survival modeling

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--gene` | str | ✅ | — | HGNC gene symbol (e.g. `TP53`) |
| `--mode` | str | ❌ | `pan-cancer` | `pan-cancer` or `tumor-vs-normal` |
| `--cancer` | str | ❌ | all | TCGA project ID for `tumor-vs-normal` mode (e.g. `TCGA-LUAD`) |
| `--outdir` | str | ❌ | `results/` | Directory where outputs are written |
| `--dpi` | int | ❌ | `300` | Plot resolution |
| `--top-n` | int | ❌ | `15` | Number of cancer types to show in pan-cancer plot |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `{gene}_pancancer_expression.png` | PNG (300 DPI) | Box plots across TCGA cohorts |
| `{gene}_pancancer_expression.svg` | SVG | Vector version |
| `{gene}_tumor_vs_normal.png` | PNG | Tumor vs normal comparison (mode=tumor-vs-normal) |
| `{gene}_expression.tsv` | TSV | Expression values per sample with project + sample_type columns |

## Execution policy

- Uses GDC REST API (free, no auth required for open-access tier)
- Caches JSON responses to `~/.cache/genomics-skills/tcga-expression/`
- Queries gene by Ensembl ID (auto-resolved via MyGene.info)
- Retries 3× with exponential backoff; warns and exits gracefully on failure
- Note: GDC open-access data returns log2(FPKM+1) values

## Trigger phrases (for agent routing)

- "TCGA expression for {gene}"
- "pan-cancer expression of {gene}"
- "is {gene} upregulated in {cancer}"
- "expression across cancer types for {gene}"
- "tumor vs normal expression {gene} {cancer}"
- "how is {gene} expressed in TCGA"

## Example invocation

```bash
python skills/tcga-expression/scripts/tcga_expression.py \
    --gene TP53 --mode pan-cancer --outdir results/

python skills/tcga-expression/scripts/tcga_expression.py \
    --gene TP53 --mode tumor-vs-normal --cancer TCGA-LUAD --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.

---
**Author:** Ankur Sharma — Agentic AI · Machine Learning · Bioinformatics · Data Science  
**GitHub:** [ankurgenomics](https://github.com/ankurgenomics)  
**Project:** [genomics-skills](https://github.com/ankurgenomics/genomics-skills)
