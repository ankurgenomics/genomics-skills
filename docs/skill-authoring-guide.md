# Skill Authoring Guide

## Overview

Every skill in `genomics-skills` is a self-contained directory under `skills/`.
The structure is:

```
skills/
└── your-skill-name/
    ├── SKILL.md          ← agent-readable contract (REQUIRED)
    ├── scripts/
    │   └── your_skill.py ← pure-Python implementation (REQUIRED)
    └── requirements.txt  ← per-skill pip dependencies (REQUIRED)
```

Naming: use lowercase kebab-case for the directory (`my-skill`), snake_case for
the Python file (`my_skill.py`).

---

## SKILL.md template

```markdown
# SKILL: <skill-name>

## Purpose
One-paragraph description of what this skill does, what inputs it accepts,
and what it produces. Write for both humans and agents.

## Inputs

| Flag | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `--flag` | str | ✅ | — | Description |

## Outputs

| File | Format | Description |
|------|--------|-------------|
| `output.tsv` | TSV | Description |

## Execution policy

- Data sources used (URLs, APIs)
- Caching behaviour
- Rate limits / auth requirements
- Failure modes

## Trigger phrases (for agent routing)

- "natural language phrase that would invoke this skill"
- "another phrase"

## Example invocation

\```bash
python skills/your-skill/scripts/your_skill.py --flag value
\```

## Dependencies

See `requirements.txt` in this directory.
```

---

## Python script conventions

1. **argparse interface** — every flag documented in SKILL.md must correspond to
   an `add_argument()` call. The `--help` output is what agents call
   `list_skill_scripts` to inspect.

2. **Explicit outputs** — print `[output] <path>` for every file written so the
   runner can collect them.

3. **Graceful failure** — if data is unavailable, write an empty TSV with correct
   headers and exit with code 0 with a `[warn]` message. Only exit non-zero on
   programming errors.

4. **Caching** — use `~/.cache/genomics-skills/<skill-name>/` for API responses.
   Key by a hash of the relevant inputs.

5. **Output formats** — prefer TSV over CSV, PNG+SVG over PDF, JSON manifests
   over custom formats. Name files predictably: `{gene}_{skill}_{type}.{ext}`.

6. **No R dependency** — pure Python stack only: numpy, pandas, matplotlib,
   seaborn, scipy, lifelines, biopython, requests, httpx.

---

## Checklist before opening a PR

- [ ] `SKILL.md` complete (purpose, inputs table, outputs table, trigger phrases, example)
- [ ] `scripts/your_skill.py` has `--help` that matches SKILL.md inputs
- [ ] `requirements.txt` lists all non-stdlib dependencies
- [ ] `[output] <path>` printed for every file written
- [ ] Tested locally: `python skills/your-skill/scripts/your_skill.py --help` exits 0
- [ ] `genomics-skill list` shows your new skill
- [ ] At least one smoke test added to `tests/test_skills_smoke.py`

---
*Authored by [Ankur Sharma](https://github.com/ankurgenomics) — Agentic AI · Machine Learning · Bioinformatics · Data Science*
