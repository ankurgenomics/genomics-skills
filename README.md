# genomics-skills

> **Built by [Ankur Sharma](https://github.com/ankurgenomics) — Agentic AI · Machine Learning · Bioinformatics · Data Science**  
> A production-quality, agent-friendly genomics skill library demonstrating end-to-end ML/AI engineering for life sciences.

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue?style=flat-square)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green?style=flat-square)](LICENSE)
[![Skills: 8](https://img.shields.io/badge/skills-8-0b8a7a?style=flat-square)](skills/)
[![Tests: 14 passed](https://img.shields.io/badge/tests-14%20passed-brightgreen?style=flat-square)](tests/)
[![Agentic AI](https://img.shields.io/badge/Agentic-AI-blueviolet?style=flat-square)]()
[![LLM Routing](https://img.shields.io/badge/LLM-Claude%20Haiku-orange?style=flat-square)]()

---

## About the Author

**Ankur Sharma** — AI Engineer · ML Scientist · Bioinformatician  
🔗 [GitHub: @ankurgenomics](https://github.com/ankurgenomics) · [agentic-genomics](https://github.com/ankurgenomics/agentic-genomics)

### Skills demonstrated in this repo

| Domain | Technologies |
|---|---|
| **Agentic AI** | LLM-backed skill routing (Claude Haiku), SKILL.md agent contracts, tool-calling patterns |
| **Machine Learning** | Survival analysis (Cox regression), statistical modelling, feature engineering from omics data |
| **Data Science** | Pandas, NumPy, SciPy, Seaborn/Matplotlib, Parquet caching, TSV/JSON pipelines |
| **Bioinformatics** | TCGA, cBioPortal REST API, Kaplan-Meier, GO/KEGG enrichment, PDB protein structure, PubMed |
| **Software Engineering** | Modular Python package, argparse CLI, pytest (14 tests), pyproject.toml, CI-ready |
| **MLOps / Data Eng.** | REST API integration, Parquet caching, reproducible outputs, gitignore hygiene |
| **NLP / LLMs** | Natural-language skill routing via Anthropic API, prompt engineering |
| **DevOps** | Git, GitHub, virtual environments, pip-installable package with optional dependency groups |

---

## Why this repo exists

My flagship project [**agentic-genomics**](https://github.com/ankurgenomics/agentic-genomics) is a full
**LangGraph agentic pipeline** that takes a patient VCF + HPO phenotype terms and returns ranked candidate
variants with ACMG-lite evidence chains.

This repo is the **downstream skill layer** — 8 deterministic, agent-callable tools that answer follow-up
genomics questions using real public data:

> *"What does expression of this gene look like across TCGA cancer types?"*  
> *"Is high expression prognostic in LUAD?"*  
> *"Where does this missense land on the 3-D protein structure?"*  
> *"What do the last 20 papers say about this variant?"*

---

## Example outputs

### Pan-cancer expression — TP53 across 31 TCGA projects (9,479 real patient samples)

![TP53 pan-cancer expression](examples/TP53_pancancer_expression.png)

*Data: TCGA Pan-Cancer Atlas 2018 · RNA-seq v2 RSEM via cBioPortal REST API · log₂(RSEM+1)*

---

### Kaplan-Meier survival — TP53 expression in TCGA-LUAD (501 real patients)

![TP53 LUAD Kaplan-Meier](examples/TP53_TCGA_LUAD_km_os.png)

*High vs. low expression split at median · Overall survival · log-rank p-value computed*

---

## Skill catalog

| # | Skill | Agentic use case | ML / DS techniques | Key outputs |
|---|-------|------------------|--------------------|-------------|
| 1 | [`tcga-expression`](skills/tcga-expression/) | "Show expression across cancers" | Data fetching, log-transform, boxplots | TSV + PNG/SVG, 9,479 samples |
| 2 | [`survival-analysis`](skills/survival-analysis/) | "Is this gene prognostic?" | Cox PH regression, KM estimator, log-rank test | KM curves, survival TSV |
| 3 | [`variant-context`](skills/variant-context/) | "Annotate this variant" | REST API integration, lollipop viz | Lollipop PNG, hotspot TSV |
| 4 | [`protein-variant-mapper`](skills/protein-variant-mapper/) | "Show variant on 3-D structure" | PDB parsing, 3-D structure rendering | Interactive HTML, lollipop PNG |
| 5 | [`protein-structure-viewer`](skills/protein-structure-viewer/) | "View the protein structure" | PDB file parsing, pocket detection | HTML viewer, B-factor plot |
| 6 | [`go-enrichment`](skills/go-enrichment/) | "What pathways is this gene in?" | Gene set enrichment, REST API | GO/KEGG TSV, bubble plot |
| 7 | [`pubmed-search`](skills/pubmed-search/) | "Recent papers on this gene?" | NLP digest, citation trends | TSV, markdown digest, timeline plot |
| 8 | [`plot-volcano`](skills/plot-volcano/) | "Make a publication figure" | DE result visualization, FDR thresholds | 300 DPI PNG + SVG |

---

## Architecture

```
genomics-skills/
│
├── skills/                          ← one folder per skill (agent-discoverable)
│   ├── tcga-expression/
│   │   ├── SKILL.md                 ← agent-readable contract (inputs · outputs · triggers)
│   │   ├── scripts/tcga_expression.py   ← pure Python, cBioPortal REST API
│   │   └── requirements.txt
│   └── ... (7 more skills)
│
├── src/genomics_skills/
│   ├── cli.py      ← `genomics-skill run / list / info / suggest` (Typer CLI)
│   └── runner.py   ← skill discovery engine + LLM routing (Claude Haiku)
│
├── examples/       ← real output images from test runs
├── tests/          ← pytest suite (14 tests: discovery + smoke)
└── pyproject.toml  ← pip-installable package with optional dependency groups
```

**Key engineering decisions:**
- **Parquet caching** at `~/.cache/genomics-skills/` — repeat queries are instant
- **LLM routing** via Claude Haiku: `genomics-skill suggest "natural language query"` maps to the right skill
- **No auth, no downloads** — all data fetched via public REST APIs (cBioPortal, MyVariant.info, NCBI)
- **Stable output formats** (TSV, PNG/SVG, JSON) for downstream chaining without glue code

---

## Quick start

```bash
git clone https://github.com/ankurgenomics/genomics-skills.git
cd genomics-skills
python -m venv .venv && source .venv/bin/activate
pip install -e ".[all]"

# Run a skill directly
python skills/tcga-expression/scripts/tcga_expression.py \
    --gene TP53 --mode pan-cancer --outdir results/

# Or via CLI
genomics-skill run tcga-expression --gene TP53 --mode pan-cancer
genomics-skill list
genomics-skill suggest "show me survival data for BRCA1 in breast cancer"
```

---

## CLI reference

```
genomics-skill list                          # list all 8 skills
genomics-skill info <skill-name>             # show SKILL.md contract
genomics-skill run  <skill-name> [args...]   # execute a skill
genomics-skill suggest "<natural language>"  # AI routing via Claude Haiku
```

### More example runs

```bash
# Survival analysis — KRAS in colorectal cancer
genomics-skill run survival-analysis --gene KRAS --cancer TCGA-COAD --endpoint os --outdir results/

# GO enrichment for a gene list
genomics-skill run go-enrichment --genes TP53,BRCA1,EGFR,KRAS --outdir results/

# PubMed search
genomics-skill run pubmed-search --query "BRCA1 homologous recombination PARP inhibitor" --max-results 25 --outdir results/

# AI-powered skill routing
genomics-skill suggest "what does EGFR expression look like in lung cancer?"
genomics-skill suggest "run survival analysis for MYC in glioblastoma"
```

---

## Data sources

| Skill | Source | Scale |
|---|---|---|
| `tcga-expression` | [cBioPortal REST API](https://www.cbioportal.org/api) | 9,479 samples · 31 TCGA cancer types |
| `survival-analysis` | [cBioPortal REST API](https://www.cbioportal.org/api) | Real OS/DFS clinical data |
| `variant-context` | [MyVariant.info](https://myvariant.info) | ClinVar + COSMIC annotations |
| `protein-variant-mapper` | [UniProt](https://www.uniprot.org) + [AlphaFold](https://alphafold.ebi.ac.uk) | AF2 + PDB structures |
| `protein-structure-viewer` | [RCSB PDB](https://www.rcsb.org) | Any PDB entry |
| `pubmed-search` | [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25497/) + Europe PMC | Real-time PubMed |
| `go-enrichment` | [g:Profiler REST API](https://biit.cs.ut.ee/gprofiler/) | GO · KEGG · Reactome |
| `plot-volcano` | Local DE results TSV | User-provided |

---

## How this fits into the agentic pipeline

```
proband.vcf + HPO terms
        │
        ▼
  agentic-genomics  ← LangGraph agent (see: github.com/ankurgenomics/agentic-genomics)
  → variant ranking + ACMG-lite evidence
        │
        ▼  top candidate genes
  genomics-skills  ← this repo
  → tcga-expression     → dysregulation across cancer types
  → survival-analysis   → prognostic value
  → variant-context     → somatic hotspot landscape
  → protein-variant-mapper → 3-D structural context
  → go-enrichment       → pathway membership
  → pubmed-search       → recent literature
        │
        ▼
  Consolidated research report (markdown + figures)
```

---

## Running tests

```bash
pytest tests/ -v   # → 14 passed in ~6s
```

Covers: skill auto-discovery, SKILL.md contract validation, and `--help` smoke tests for all 8 skills.

---

## Roadmap

| Skill | Status |
|-------|--------|
| `tcga-expression` | 🟢 Live — 9,479 samples · 31 cancers · cBioPortal |
| `survival-analysis` | 🟢 Live — KM + Cox from real TCGA clinical data |
| `variant-context` | 🟢 MVP |
| `protein-variant-mapper` | 🟢 MVP |
| `protein-structure-viewer` | 🟢 MVP |
| `go-enrichment` | 🟢 MVP |
| `pubmed-search` | 🟢 MVP |
| `plot-volcano` | 🟢 MVP |
| `crispr-sgrna-design` | 🔵 Planned |
| `coexpression-network` | 🔵 Planned — STRING PPI |
| `depmap-essentiality` | 🔵 Planned — DepMap |
| `single-cell-marker` | 🔵 Planned — scRNA-seq |

---

## Author

**Ankur Sharma**  
AI Engineer · Machine Learning Scientist · Bioinformatician · Data Scientist

📧 learning.ankur.ai@gmail.com  
🔗 [GitHub: @ankurgenomics](https://github.com/ankurgenomics)  
🧬 [agentic-genomics](https://github.com/ankurgenomics/agentic-genomics) — LangGraph variant interpretation pipeline  
🛠 [genomics-skills](https://github.com/ankurgenomics/genomics-skills) — this repo

**Core competencies demonstrated here:**  
`Agentic AI` · `LLM Integration` · `REST API Engineering` · `Machine Learning` · `Survival Analysis` · `Data Science` · `Python` · `Bioinformatics` · `TCGA / Genomics` · `Software Engineering` · `pytest` · `Modular Architecture`

---

## License

MIT — see [LICENSE](LICENSE).

> ⚠️ Research demonstration only. Not for clinical use. Not a medical device.  
> Data sourced from public APIs — see each source's terms of use.
