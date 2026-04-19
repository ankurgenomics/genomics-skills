# Integration Guide: genomics-skills ↔ agentic-genomics

## Overview

`agentic-genomics` identifies top candidate variants from a proband VCF using
the `GenomicsCopilot` LangGraph agent. `genomics-skills` provides downstream
genomic context for each candidate gene.

The integration point is an `OmicsContextAgent` that fires after variant
prioritization and calls genomics-skills as LangGraph tools.

---

## Architecture

```
proband.vcf + HPO terms
        │
        ▼
  GenomicsCopilot (agentic-genomics)
  → ACMG-lite evidence chains
  → Ranked candidate variants
        │
        ▼  top N candidate genes extracted from rankings
  OmicsContextAgent (calls genomics-skills)
  → variant-context         somatic hotspot landscape
  → tcga-expression         pan-cancer expression context
  → survival-analysis       is this gene prognostic?
  → protein-variant-mapper  patient variant on 3-D structure
  → go-enrichment           pathway context
  → pubmed-search           recent literature
        │
        ▼
  Consolidated markdown report + figures
```

---

## Wiring genomics-skills into agentic-genomics as LangGraph tools

Install both packages:

```bash
pip install -e /path/to/agentic-genomics
pip install -e /path/to/genomics-skills[all]
```

Example tool wrapper (add to `src/agentic_genomics/agents/omics_context/tools/`):

```python
from langchain_core.tools import tool
from genomics_skills.runner import run_skill

@tool
def tcga_expression(gene: str, mode: str = "pan-cancer", cancer: str = "all") -> str:
    """
    Fetch TCGA gene expression for a candidate gene.
    Returns path to output TSV and PNG plot.
    """
    import tempfile, os
    outdir = os.path.join(tempfile.mkdtemp(), "tcga_expression")
    rc = run_skill("tcga-expression", [
        "--gene", gene, "--mode", mode, "--cancer", cancer,
        "--outdir", outdir,
    ])
    if rc != 0:
        return f"tcga-expression failed for {gene}"
    return f"Outputs written to {outdir}/"

@tool
def survival_analysis(gene: str, cancer: str, endpoint: str = "os") -> str:
    """Run Kaplan-Meier survival analysis for a gene in a TCGA cohort."""
    import tempfile, os
    outdir = os.path.join(tempfile.mkdtemp(), "survival")
    rc = run_skill("survival-analysis", [
        "--gene", gene, "--cancer", cancer, "--endpoint", endpoint,
        "--outdir", outdir,
    ])
    return f"Outputs written to {outdir}/" if rc == 0 else f"survival-analysis failed for {gene}"
```

---

## Calling skills directly (no agent)

Each skill is a standalone Python script — call it from any pipeline:

```bash
# After agentic-genomics identifies TP53 as top candidate:
python skills/variant-context/scripts/variant_context.py \
    --gene TP53 --cancer LUAD --outdir reports/TP53/

python skills/tcga-expression/scripts/tcga_expression.py \
    --gene TP53 --mode pan-cancer --outdir reports/TP53/

python skills/survival-analysis/scripts/survival_analysis.py \
    --gene TP53 --cancer TCGA-LUAD --endpoint os --outdir reports/TP53/

python skills/protein-variant-mapper/scripts/protein_variant_mapper.py \
    --gene TP53 --variants R175H --outdir reports/TP53/
```

---

## Using with bioinfor-claw / OpenClaw

Register the skills directory in `~/.openclaw/openclaw.json`:

```json
{
  "skills": {
    "load": {
      "extraDirs": ["/path/to/genomics-skills/skills"]
    }
  }
}
```

All 8 skills will be auto-discovered on next OpenClaw restart.

---

## Output compatibility

All skill outputs use stable, predictable paths:
- `{gene}_{skill_type}.tsv` — for downstream pandas reads
- `{gene}_{skill_type}.png` — for report embedding
- `{gene}_{skill_type}.html` — for interactive viewers

This lets you chain skills without glue code:

```python
import pandas as pd

# Read survival output from survival-analysis skill
surv = pd.read_csv("reports/TP53/TP53_TCGA_LUAD_survival.tsv", sep="\t")

# Feed high-expression gene list into go-enrichment skill
high_genes = surv[surv["group"] == "High"]["gene_id"].tolist()
# ... write to file, then pass to go-enrichment --genes
```
