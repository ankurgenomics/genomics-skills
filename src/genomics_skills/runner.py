"""Skill discovery and execution engine.
# Author: Ankur Sharma | Agentic AI · Machine Learning · Bioinformatics · Data Science
# GitHub: https://github.com/ankurgenomics | Portfolio: agentic-genomics

Scans the skills/ directory for SKILL.md files and resolves the matching
Python script to run. Includes an optional LLM-backed skill router that
uses ANTHROPIC_API_KEY to match a natural-language query to the right skill.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

# Root of the repo (two levels up from src/genomics_skills/)
_REPO_ROOT = Path(__file__).parent.parent.parent
_SKILLS_DIR = _REPO_ROOT / "skills"


def list_skills() -> list[dict]:
    """Return metadata for every discovered skill."""
    skills = []
    for skill_dir in sorted(_SKILLS_DIR.iterdir()):
        skill_md = skill_dir / "SKILL.md"
        scripts_dir = skill_dir / "scripts"
        if not skill_md.exists():
            continue
        # Extract Purpose line from SKILL.md
        purpose = ""
        for line in skill_md.read_text().splitlines():
            if line.startswith("## Purpose"):
                continue
            if line.strip() and purpose == "":
                purpose = line.strip()
            if purpose:
                break
        scripts = sorted(scripts_dir.glob("*.py")) if scripts_dir.exists() else []
        skills.append(
            {
                "name": skill_dir.name,
                "path": skill_dir,
                "skill_md": skill_md,
                "scripts": scripts,
                "purpose": purpose,
            }
        )
    return skills


def find_skill(name: str) -> dict | None:
    """Find a skill by exact name."""
    for s in list_skills():
        if s["name"] == name:
            return s
    return None


def run_skill(name: str, extra_args: list[str]) -> int:
    """Locate skill script and execute it with extra_args."""
    skill = find_skill(name)
    if not skill:
        print(f"[error] Skill '{name}' not found. Run `genomics-skill list` to see available skills.",
              file=sys.stderr)
        return 1

    if not skill["scripts"]:
        print(f"[error] No scripts found in {skill['path'] / 'scripts'}/", file=sys.stderr)
        return 1

    script = skill["scripts"][0]
    cmd = [sys.executable, str(script)] + extra_args
    print(f"[run] {' '.join(cmd)}")
    result = subprocess.run(cmd)
    return result.returncode


def suggest_skill(query: str) -> dict | None:
    """Use Claude (Anthropic) to match a natural-language query to the best skill.

    Requires ANTHROPIC_API_KEY in the environment.
    Returns the matched skill dict, or None if no match / API unavailable.

    Example
    -------
    >>> skill = suggest_skill("show me TP53 expression across cancer types")
    >>> print(skill["name"])  # tcga-expression
    """
    import urllib.request

    api_key = os.environ.get("ANTHROPIC_API_KEY", "")
    if not api_key:
        print("[llm] ANTHROPIC_API_KEY not set — skipping LLM routing.", file=sys.stderr)
        return None

    skills = list_skills()
    skill_list = "\n".join(
        f'- {s["name"]}: {s["purpose"]}' for s in skills
    )

    prompt = (
        f"You are a genomics skill router. Given a user query, return ONLY the skill name "
        f"(exact, no explanation) that best matches from this list:\n\n{skill_list}\n\n"
        f"If nothing matches, return: none\n\nUser query: {query}\n\nSkill name:"
    )

    payload = json.dumps({
        "model": "claude-3-5-haiku-20241022",
        "max_tokens": 32,
        "messages": [{"role": "user", "content": prompt}],
    }).encode()

    req = urllib.request.Request(
        "https://api.anthropic.com/v1/messages",
        data=payload,
        headers={
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            body = json.loads(resp.read())
        skill_name = body["content"][0]["text"].strip().lower()
        if skill_name == "none":
            print(f"[llm] No matching skill found for: {query!r}", file=sys.stderr)
            return None
        matched = find_skill(skill_name)
        if matched:
            print(f"[llm] Matched skill: {matched['name']}")
        else:
            print(f"[llm] LLM returned unknown skill name: {skill_name!r}", file=sys.stderr)
        return matched
    except Exception as exc:
        print(f"[llm] API call failed: {exc}", file=sys.stderr)
        return None
