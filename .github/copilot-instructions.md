# Copilot / AI coding agent instructions for genomics-skills

These brief, repo-specific instructions help an automated coding agent be immediately productive.

1. Big picture
- This repo is a modular "skill" library: one folder per skill under `skills/`. Each skill ships a `SKILL.md`, a `requirements.txt`, and a single pure-Python script under `scripts/` that implements the contract.
- The package entrypoints are `src/genomics_skills/cli.py` (CLI) and `src/genomics_skills/runner.py` (skill discovery + execution engine). The runner auto-discovers skills by scanning `skills/*/SKILL.md`.
- Design principle: skills compute (deterministic), agents decide. Do not add agent logic into skill scripts.

2. How skills are structured (examples)
- Example skill: `skills/tcga-expression/`
  - Contract: `skills/tcga-expression/SKILL.md` (inputs, outputs, trigger phrases). Agents rely on the contract — preserve it when editing.
  - Implementation: `skills/tcga-expression/scripts/tcga_expression.py` (uses `argparse` flags and writes results to `--outdir`).
  - Dependencies: `skills/tcga-expression/requirements.txt`.
- When authoring or modifying a skill, update these three places: `SKILL.md`, the script (keep `argparse` flags stable), `requirements.txt`.

3. Developer workflows (commands you may need to run)
- Create virtualenv and install everything: `python -m venv .venv && source .venv/bin/activate && pip install -e ".[all]"`
- Run a skill directly (script):
  `python skills/tcga-expression/scripts/tcga_expression.py --gene TP53 --mode pan-cancer --outdir results/`
- Run via CLI: `genomics-skill run tcga-expression --gene TP53 --mode pan-cancer`
- List skills: `genomics-skill list`
- Tests: run `pytest` from the repo root. Relevant tests: `tests/test_runner.py`, `tests/test_skills_smoke.py`.

4. Project-specific conventions & patterns
- One purpose per skill: each script focuses on a single, narrowly-scoped operation. Avoid adding multi-responsibility behavior.
- SKILL.md is the canonical contract. Agents parse it. Do not remove or radically change fields unless you update dependent agent wiring (see `docs/integration-guide.md`).
- Scripts expose `argparse` CLI flags and write outputs to an `--outdir` or explicit output path. Search `skills/*/scripts/*.py` to find flag naming patterns.
- Output formats are stable and machine-readable: TSV/CSV, JSON manifests, PNG/SVG. Prefer these formats for new outputs unless a SKILL.md says otherwise.
- Pure Python stack only — numpy / pandas /matplotlib /lifelines /biopython. Avoid introducing large non-Python toolchains.

5. Integration points & external dependencies
- Skills call public APIs/databases: TCGA/GDC, MyVariant.info, UniProt, NCBI (PubMed). Sensitive env vars: `NCBI_EMAIL` and optionally `ANTHROPIC_API_KEY` used by `pubmed-search` and LLM-backed runner. See `.env.example` referenced in `README.md`.
- The repo is meant to be consumed by `agentic-genomics` and other agent runners that scan `SKILL.md`. Preserve backward-compatible SKILL.md changes.

6. What to check before making edits
- If changing a skill's CLI flags, update `SKILL.md` accordingly and run `tests/test_skills_smoke.py`.
- If adding dependencies, add them to the skill's `requirements.txt` and ensure `pip install -e ".[all]"` still completes.
- If touching `src/genomics_skills/runner.py`, run `tests/test_runner.py` — the runner is the discovery + execution integration point.

7. Examples of actionable tasks for an AI agent
- Add a new skill: create `skills/<name>/SKILL.md`, `skills/<name>/scripts/<name>.py` with `argparse`, `requirements.txt`. Follow `docs/skill-authoring-guide.md`.
- Fix a SKILL.md input mismatch: update `SKILL.md` to match the script's `argparse` flags or vice versa; run smoke tests.
- Update outputs to include a JSON manifest: add file writing in script and declare file in `SKILL.md` outputs.

8. Where to look for more context
- `README.md` (root): repo architecture, quick-start commands, examples.
- `docs/skill-authoring-guide.md` and `docs/integration-guide.md`: authoring checklist and agent wiring.
- `src/genomics_skills/cli.py` and `src/genomics_skills/runner.py`: discovery and CLI behaviour.
- Representative skills under `skills/` (e.g. `tcga-expression`, `pubmed-search`, `variant-context`).

If anything in these instructions is unclear or you want more detail on any section (tests, a specific skill, or the runner), tell me which part to expand and I will iterate.
