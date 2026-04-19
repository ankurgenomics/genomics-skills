"""Smoke tests for individual skill scripts — check they parse args and exit cleanly."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

SKILLS_DIR = Path(__file__).parent.parent / "skills"


def _script(skill_name: str) -> Path:
    scripts = list((SKILLS_DIR / skill_name / "scripts").glob("*.py"))
    assert scripts, f"No script found for {skill_name}"
    return scripts[0]


def test_plot_volcano_help():
    """plot-volcano should print help without error."""
    result = subprocess.run(
        [sys.executable, str(_script("plot-volcano")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--input" in result.stdout


def test_variant_context_help():
    result = subprocess.run(
        [sys.executable, str(_script("variant-context")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--gene" in result.stdout


def test_tcga_expression_help():
    result = subprocess.run(
        [sys.executable, str(_script("tcga-expression")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--gene" in result.stdout


def test_survival_analysis_help():
    result = subprocess.run(
        [sys.executable, str(_script("survival-analysis")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--gene" in result.stdout


def test_go_enrichment_help():
    result = subprocess.run(
        [sys.executable, str(_script("go-enrichment")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--genes" in result.stdout


def test_pubmed_search_help():
    result = subprocess.run(
        [sys.executable, str(_script("pubmed-search")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--query" in result.stdout


def test_protein_structure_viewer_help():
    result = subprocess.run(
        [sys.executable, str(_script("protein-structure-viewer")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--gene" in result.stdout


def test_protein_variant_mapper_help():
    result = subprocess.run(
        [sys.executable, str(_script("protein-variant-mapper")), "--help"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0
    assert "--gene" in result.stdout
