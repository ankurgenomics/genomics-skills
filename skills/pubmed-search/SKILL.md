# SKILL: pubmed-search

## Purpose
Search PubMed via NCBI E-utilities for literature on a gene, variant, or topic.
Produces a ranked TSV of papers (with citation count from Europe PMC), a
structured markdown digest of the top 5 papers, and a publication trend chart.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--query` | str | ✅ | — | PubMed query string (e.g. `TP53 variant pathogenicity`) |
| `--max-results` | int | ❌ | `50` | Maximum papers to fetch |
| `--date-from` | str | ❌ | `2018` | Year range start |
| `--date-to` | str | ❌ | current | Year range end |
| `--outdir` | str | ❌ | `results/` | Output directory |
| `--dpi` | int | ❌ | `300` | Plot DPI |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `pubmed_results.tsv` | TSV | PMID, title, authors, year, journal, abstract, citation count |
| `pubmed_digest.md` | Markdown | Structured digest of top 5 papers |
| `pubmed_trend.png` | PNG | Publications per year bar chart |

## Execution policy

- Uses NCBI E-utilities (free; set `NCBI_EMAIL` env var; `NCBI_API_KEY` raises rate limit)
- Fetches citation counts from Europe PMC REST API
- Caches query results for 24 h in `~/.cache/genomics-skills/pubmed-search/`
- Rate-limited: 3 req/s without API key, 10 req/s with key

## Trigger phrases (for agent routing)

- "search PubMed for {query}"
- "recent papers on {gene} {variant}"
- "literature for {gene}"
- "find papers about {topic}"
- "publication trend for {gene}"

## Example invocation

```bash
python skills/pubmed-search/scripts/pubmed_search.py \
    --query "TP53 R175H pathogenicity" --max-results 30 --outdir results/
```

## Dependencies

See `requirements.txt` in this directory.
