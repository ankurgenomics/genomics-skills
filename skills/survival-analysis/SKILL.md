# SKILL: survival-analysis

## Purpose
Given a gene symbol and TCGA cohort, split samples into high/low expression
groups and produce Kaplan–Meier overall survival (OS) and disease-free survival
(DFS) curves, with log-rank test p-values and at-risk tables.

Optionally runs multivariate Cox proportional hazards regression with clinical
covariates (age, stage, grade) when clinical data is available.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--gene` | str | ✅ | — | HGNC gene symbol |
| `--cancer` | str | ✅ | — | TCGA project ID (e.g. `TCGA-LUAD`) |
| `--endpoint` | str | ❌ | `os` | Survival endpoint: `os` (overall) or `dfs` (disease-free) |
| `--split` | str | ❌ | `median` | How to split: `median`, `tertile`, or `quartile` |
| `--covariates` | str | ❌ | — | Comma-separated clinical covariates for Cox (e.g. `age,stage`) |
| `--outdir` | str | ❌ | `results/` | Output directory |
| `--dpi` | int | ❌ | `300` | Plot resolution |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `{gene}_{cancer}_km_{endpoint}.png` | PNG (300 DPI) | KM curves with at-risk table + log-rank p |
| `{gene}_{cancer}_km_{endpoint}.svg` | SVG | Vector version |
| `{gene}_{cancer}_cox.png` | PNG | Forest plot (Cox HR) — only if `--covariates` given |
| `{gene}_{cancer}_survival.tsv` | TSV | Sample-level: group label, time, event |
| `{gene}_{cancer}_cox_results.tsv` | TSV | Cox hazard ratios, CIs, p-values |

## Execution policy

- Uses `lifelines` for KM estimation, log-rank test, and Cox PH model
- Clinical / survival data pulled from GDC open-access tier
- Caches data to `~/.cache/genomics-skills/survival-analysis/`
- Minimum 10 events per group required for KM curve; warns and skips if not met
- Cox model only run if `--covariates` flag is provided

## Trigger phrases (for agent routing)

- "survival analysis for {gene} in {cancer}"
- "KM curve for {gene} {cancer}"
- "is {gene} expression prognostic in {cancer}"
- "overall survival {gene} TCGA {cancer}"
- "disease-free survival {gene}"
- "Cox regression {gene} {cancer} with age and stage"
- "does {gene} predict survival"

## Example invocation

```bash
python skills/survival-analysis/scripts/survival_analysis.py \
    --gene TP53 --cancer TCGA-LUAD --endpoint os --split median --outdir results/

python skills/survival-analysis/scripts/survival_analysis.py \
    --gene EGFR --cancer TCGA-LUAD --endpoint os \
    --covariates age,stage --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.

---
**Author:** Ankur Sharma — Agentic AI · Machine Learning · Bioinformatics · Data Science  
**GitHub:** [ankurgenomics](https://github.com/ankurgenomics)  
**Project:** [genomics-skills](https://github.com/ankurgenomics/genomics-skills)
