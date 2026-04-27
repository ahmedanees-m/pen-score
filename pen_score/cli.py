"""CLI for pen-score."""

from __future__ import annotations

from typing import Any

import click
from rich.console import Console

console = Console()


@click.group()
@click.version_option()
def main() -> None:
    """pen-score: Multi-axis scoring framework for programmable genome editors."""


@main.command("score-editor")
@click.argument("accession")
@click.option(
    "--use-case",
    default="human_therapeutic_aav_insertion",
    show_default=True,
    help="Use-case profile key.",
)
@click.option("--output", "-o", default=None, help="Output Parquet path.")
@click.option(
    "--exclude-axes",
    default=None,
    help=(
        "Comma-separated axis names to exclude from the composite PenScore.  "
        "Weights are renormalised over remaining axes.  "
        "Example: --exclude-axes S_Mature,S_Immuno"
    ),
)
def score_editor(
    accession: str, use_case: str, output: str | None, exclude_axes: str | None
) -> None:
    """Score a single editor by UniProt accession."""
    from pen_score.api import Scorer

    exclude = [ax.strip() for ax in exclude_axes.split(",")] if exclude_axes else None
    scorer = Scorer.load()
    result = scorer.score_editor(accession=accession, use_case=use_case, exclude_axes=exclude)
    console.print(result.model_dump())
    if output:
        import pandas as pd

        pd.DataFrame([result.model_dump()]).to_parquet(output, index=False)
        console.print(f"[green]Saved to {output}[/green]")


@main.command("score-all-editors")
@click.option("--output", "-o", default="scorecard.parquet", show_default=True)
@click.option(
    "--use-case",
    default="human_therapeutic_aav_insertion",
    show_default=True,
)
def score_all_editors(output: str, use_case: str) -> None:
    """Score all editors in the curated universe and write the public scorecard."""
    from pen_score.api import Scorer

    scorer = Scorer.load()
    scorecard = scorer.get_scorecard(use_case=use_case)
    scorecard.to_parquet(output, index=False)
    console.print(f"[green]Scorecard written to {output} ({len(scorecard)} editors)[/green]")


@main.command("select")
@click.option("--use-case", required=True, help="Use-case profile key.")
@click.option("--top-k", default=5, show_default=True, help="Number of editors to return.")
@click.option(
    "--require-dsb-free",
    is_flag=True,
    default=False,
    help="Restrict candidates to DSB-free editors (S_DSB >= 0.85).",
)
@click.option("--output", "-o", default=None, help="Save results to Parquet path.")
def select_editor(use_case: str, top_k: int, require_dsb_free: bool, output: str | None) -> None:
    """Return top-k ranked editors for a given use case, with axis reasoning."""
    from rich import box
    from rich.table import Table

    from pen_score.api import Scorer

    scorer = Scorer.load()
    ranked = scorer.select_editor(use_case=use_case, top_k=top_k, require_dsb_free=require_dsb_free)

    _AXES = ["S_DSB", "S_Spec", "S_Cargo", "S_Deliv", "S_Immuno", "S_Prog", "S_Mature", "S_Energy"]

    def fmt(v: Any) -> str:
        if v is None:
            return "[dim]--[/dim]"
        if v >= 0.75:
            return f"[green]{v:.4f}[/green]"
        if v >= 0.40:
            return f"[yellow]{v:.4f}[/yellow]"
        return f"[red]{v:.4f}[/red]"

    def reasoning(row: Any) -> str:
        strengths = [
            ax.replace("S_", "") for ax in _AXES if row.get(ax) is not None and row[ax] >= 0.75
        ]
        weaknesses = [
            ax.replace("S_", "") for ax in _AXES if row.get(ax) is None or row[ax] <= 0.35
        ]
        parts: list[str] = []
        if strengths:
            parts.append(f"[green]Strong[/green]: {', '.join(strengths)}")
        if weaknesses:
            parts.append(f"[red]Weak/missing[/red]: {', '.join(weaknesses)}")
        return " | ".join(parts) if parts else "All axes within range"

    dsb_note = " [DSB-free only]" if require_dsb_free else ""
    table = Table(
        box=box.ROUNDED,
        title=f"Top-{top_k} editors -- {use_case}{dsb_note}",
        show_lines=False,
    )
    table.add_column("Rank", style="dim", width=5, justify="center")
    table.add_column("Editor", style="bold white", min_width=14)
    table.add_column("PenScore", justify="right", style="cyan")
    table.add_column("DSB", justify="right")
    table.add_column("Spec", justify="right")
    table.add_column("Cargo", justify="right")
    table.add_column("Deliv", justify="right")
    table.add_column("Immuno", justify="right")
    table.add_column("Prog", justify="right")
    table.add_column("Mature", justify="right")
    table.add_column("Energy", justify="right")
    table.add_column("Reasoning", no_wrap=False, max_width=52)

    for rank, (_, row) in enumerate(ranked.iterrows(), 1):
        table.add_row(
            str(rank),
            str(row.get("editor_id", "")),
            fmt(row.get("PenScore")),
            fmt(row.get("S_DSB")),
            fmt(row.get("S_Spec")),
            fmt(row.get("S_Cargo")),
            fmt(row.get("S_Deliv")),
            fmt(row.get("S_Immuno")),
            fmt(row.get("S_Prog")),
            fmt(row.get("S_Mature")),
            fmt(row.get("S_Energy")),
            reasoning(row),
        )

    console.print(table)
    if output:
        ranked.to_parquet(output, index=False)
        console.print(f"[green]Saved to {output}[/green]")


if __name__ == "__main__":
    main()
