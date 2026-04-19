"""Tests for the skill runner — discovery and metadata."""
from __future__ import annotations

from pathlib import Path

import pytest

from genomics_skills import runner

EXPECTED_SKILLS = {
    "variant-context",
    "protein-variant-mapper",
    "protein-structure-viewer",
    "tcga-expression",
    "survival-analysis",
    "go-enrichment",
    "pubmed-search",
    "plot-volcano",
}


def test_list_skills_finds_all():
    skills = runner.list_skills()
    names = {s["name"] for s in skills}
    assert EXPECTED_SKILLS.issubset(names), f"Missing: {EXPECTED_SKILLS - names}"


def test_each_skill_has_skill_md():
    for skill in runner.list_skills():
        assert skill["skill_md"].exists(), f"SKILL.md missing for {skill['name']}"


def test_each_skill_has_script():
    for skill in runner.list_skills():
        assert skill["scripts"], f"No scripts/ found for {skill['name']}"
        assert skill["scripts"][0].exists()


def test_find_skill_returns_correct():
    skill = runner.find_skill("survival-analysis")
    assert skill is not None
    assert skill["name"] == "survival-analysis"


def test_find_skill_unknown_returns_none():
    assert runner.find_skill("does-not-exist") is None


def test_skill_purposes_nonempty():
    for skill in runner.list_skills():
        assert skill["purpose"], f"Empty purpose for {skill['name']}"
