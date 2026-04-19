# SKILL: go-enrichment

## Purpose
Given a gene list, run GO (Biological Process / Molecular Function / Cellular Component),
KEGG, and Reactome pathway enrichment using the g:Profiler REST API, then produce:
- An enrichment results TSV sorted by adjusted p-value
- A bubble plot with top significant terms (x = gene ratio, y = term, size = gene count,
  color = -log10 adjusted p-value)

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--genes` | str | ✅ | — | Comma-separated gene symbols OR path to a text file (one gene per line) |
| `--organism` | str | ❌ | `hsapiens` | g:Profiler organism code |
| `--sources` | str | ❌ | `GO:BP,KEGG,REAC` | Comma-separated enrichment sources |
| `--top-n` | int | ❌ | `20` | Number of top terms to show in bubble plot |
| `--outdir` | str | ❌ | `results/` | Output directory |
| `--dpi` | int | ❌ | `300` | Plot resolution |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `enrichment_results.tsv` | TSV | All significant terms with p-values, gene counts, gene lists |
| `enrichment_bubble.png` | PNG (300 DPI) | Bubble plot of top terms |
| `enrichment_bubble.svg` | SVG | Vector bubble plot |

## Execution policy

- Uses g:Profiler REST API (`https://biit.cs.ut.ee/gprofiler/api/gost/profile/`) — free, no auth
- Input can be symbols, Ensembl IDs, or Entrez IDs (g:Profiler handles all)
- FDR correction applied (Benjamini–Hochberg via g:Profiler)
- Minimum gene set size: 5; maximum: 500 (g:Profiler defaults)
- Caches API responses to `~/.cache/genomics-skills/go-enrichment/`

## Trigger phrases (for agent routing)

- "GO enrichment for gene list {genes}"
- "pathway analysis for {gene list}"
- "KEGG enrichment {genes}"
- "what pathways are enriched in {gene list}"
- "Reactome enrichment {genes}"
- "functional enrichment for my DE genes"

## Example invocation

```bash
python skills/go-enrichment/scripts/go_enrichment.py \
    --genes TP53,BRCA1,ATM,CHEK2,MDM2 --outdir results/

python skills/go-enrichment/scripts/go_enrichment.py \
    --genes my_gene_list.txt --sources GO:BP,KEGG --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.

---
**Author:** Ankur Sharma — Agentic AI · Machine Learning · Bioinformatics · Data Science  
**GitHub:** [ankurgenomics](https://github.com/ankurgenomics)  
**Project:** [genomics-skills](https://github.com/ankurgenomics/genomics-skills)
