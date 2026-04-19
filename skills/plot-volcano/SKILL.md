# SKILL: plot-volcano

## Purpose
Given a differential expression or association result table containing fold-change
and p-value columns, generate a publication-quality volcano plot at 300 DPI.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--input` | str | ✅ | — | Path to TSV/CSV with FC and p-value columns |
| `--fc-col` | str | ❌ | `log2FoldChange` | Column name for log2 fold-change |
| `--pval-col` | str | ❌ | `padj` | Column name for adjusted p-value |
| `--label-col` | str | ❌ | `gene` | Column name for gene/feature labels |
| `--fc-thresh` | float | ❌ | `1.0` | |log2FC| threshold |
| `--p-thresh` | float | ❌ | `0.05` | Adjusted p-value threshold |
| `--top-label` | int | ❌ | `15` | Number of top genes to label |
| `--title` | str | ❌ | `Volcano plot` | Plot title |
| `--outdir` | str | ❌ | `results/` | Output directory |
| `--dpi` | int | ❌ | `300` | Plot resolution |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `volcano.png` | PNG (300 DPI) | Publication-quality volcano plot |
| `volcano.svg` | SVG | Vector version |
| `volcano_annotated.tsv` | TSV | Input table with `significance` category column added |

## Execution policy

- Pure Python: matplotlib + adjustText (no R/ggplot2)
- adjustText used for non-overlapping gene labels
- Color scheme: up = red, down = blue, NS = grey
- Dashed lines at FC and p-value thresholds

## Trigger phrases (for agent routing)

- "volcano plot for {dataset}"
- "plot DE results as volcano"
- "show significant genes as volcano"
- "make a volcano from my DE table"

## Example invocation

```bash
python skills/plot-volcano/scripts/plot_volcano.py \
    --input de_results.tsv --fc-col log2FoldChange --pval-col padj \
    --title "Treated vs Control" --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.
