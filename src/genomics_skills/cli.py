"""genomics-skill CLI — list and run genomics skills from the terminal.

Usage
-----
genomics-skill list
genomics-skill run variant-context --gene TP53 --cancer LUAD
genomics-skill info variant-context
"""
from __future__ import annotations

import sys
import typer
from rich.console import Console
from rich.table import Table

from genomics_skills import runner

cli = typer.Typer(
    name="genomics-skill",
    help="Modular genomics skill library — list, inspect, and run skills.",
    add_completion=False,
)
app = cli  # entrypoint alias (pyproject.toml: genomics-skill = "genomics_skills.cli:app")
console = Console()


@cli.command("list")
def cmd_list():
    """List all available skills with their purpose."""
    skills = runner.list_skills()
    table = Table(title="genomics-skills", show_lines=True)
    table.add_column("Skill", style="bold cyan", no_wrap=True)
    table.add_column("Purpose", style="white")
    table.add_column("Scripts", style="dim")
    for s in skills:
        table.add_row(
            s["name"],
            s["purpose"][:80],
            str(len(s["scripts"])),
        )
    console.print(table)
    console.print(
        f"\n[dim]Run a skill:[/dim]  genomics-skill run <name> [args]\n"
        f"[dim]Read docs:[/dim]   genomics-skill info <name>"
    )


@cli.command("info")
def cmd_info(name: str = typer.Argument(..., help="Skill name")):
    """Print the SKILL.md for a skill."""
    skill = runner.find_skill(name)
    if not skill:
        console.print(f"[red]Skill '{name}' not found.[/red] Run `genomics-skill list`.")
        raise typer.Exit(1)
    console.print(skill["skill_md"].read_text())


@cli.command("run", context_settings={"allow_extra_args": True, "ignore_unknown_options": True})
def cmd_run(
    name: str = typer.Argument(..., help="Skill name"),
    extra: list[str] = typer.Argument(None),
):
    """Run a skill, passing all remaining arguments to its script."""
    rc = runner.run_skill(name, extra or [])
    raise typer.Exit(rc)


@cli.command("suggest")
def cmd_suggest(
    query: str = typer.Argument(..., help="Natural-language description of what you want to do"),
):
    """Use AI to find the right skill for a natural-language query (requires ANTHROPIC_API_KEY)."""
    skill = runner.suggest_skill(query)
    if not skill:
        console.print("[red]No matching skill found.[/red]")
        raise typer.Exit(1)
    console.print(f"\n[bold cyan]Suggested skill:[/bold cyan] {skill['name']}")
    console.print(f"[dim]{skill['purpose']}[/dim]")
    console.print(f"\nRun it with:  [bold]genomics-skill run {skill['name']} --help[/bold]")


def main():
    cli()


if __name__ == "__main__":
    main()
