"""Analyze simulated or future OVA measurements and create figures."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str((Path.cwd() / ".matplotlib_cache").resolve())
)

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import f_oneway


REQUIRED_COLUMNS = {
    "sample_id",
    "branch",
    "treatment",
    "saccharide",
    "replicate",
    "free_amino_groups_pct",
    "ige_binding_pct",
    "particle_size_nm",
    "surface_hydrophobicity_au",
}
BRANCH_ORDER = ["Non-enzymatic", "TG-assisted", "Tyr/CA-assisted"]
SACCHARIDE_ORDER = ["Control", "Lactose", "Mannotriose", "Dextran"]
PALETTE = {
    "Control": "#8C8C8C",
    "Lactose": "#4C78A8",
    "Mannotriose": "#F2A541",
    "Dextran": "#59A14F",
}


def validate_data(data: pd.DataFrame) -> None:
    """Raise a helpful error if the input schema is incomplete."""
    missing = REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    duplicated = data["sample_id"].duplicated()
    if duplicated.any():
        duplicates = data.loc[duplicated, "sample_id"].tolist()
        raise ValueError(f"Duplicate sample_id values: {duplicates}")


def compute_outcomes(data: pd.DataFrame) -> pd.DataFrame:
    """Calculate control-relative modification and IgE reduction."""
    validate_data(data)
    result = data.copy()
    controls = (
        result.loc[result["saccharide"] == "Control"]
        .groupby("branch")[["free_amino_groups_pct", "ige_binding_pct"]]
        .mean()
        .rename(
            columns={
                "free_amino_groups_pct": "control_free_amino_groups_pct",
                "ige_binding_pct": "control_ige_binding_pct",
            }
        )
    )
    missing_controls = set(result["branch"]) - set(controls.index)
    if missing_controls:
        raise ValueError(
            "Every branch needs a Control group. Missing controls for: "
            f"{sorted(missing_controls)}"
        )
    result = result.join(controls, on="branch")
    result["modification_extent_pct"] = (
        (
            result["control_free_amino_groups_pct"]
            - result["free_amino_groups_pct"]
        )
        / result["control_free_amino_groups_pct"]
        * 100
    )
    result["ige_reduction_vs_control_pct"] = (
        (result["control_ige_binding_pct"] - result["ige_binding_pct"])
        / result["control_ige_binding_pct"]
        * 100
    )
    return result


def summarize_treatments(data: pd.DataFrame) -> pd.DataFrame:
    """Create a tidy mean, SD, and count table."""
    metrics = [
        "free_amino_groups_pct",
        "modification_extent_pct",
        "ige_binding_pct",
        "ige_reduction_vs_control_pct",
        "particle_size_nm",
        "surface_hydrophobicity_au",
    ]
    summary = (
        data.groupby(["branch", "saccharide"])[metrics]
        .agg(["mean", "std", "count"])
        .round(3)
    )
    summary.columns = ["_".join(column) for column in summary.columns]
    return summary.reset_index()


def run_anova(data: pd.DataFrame) -> pd.DataFrame:
    """Run exploratory one-way ANOVA within each proposed branch."""
    records = []
    for branch, branch_data in data.groupby("branch"):
        groups = [
            group["modification_extent_pct"].dropna().to_numpy()
            for _, group in branch_data.groupby("saccharide")
            if len(group) >= 2
        ]
        if len(groups) < 2:
            continue
        statistic, p_value = f_oneway(*groups)
        records.append(
            {
                "branch": branch,
                "outcome": "modification_extent_pct",
                "f_statistic": statistic,
                "p_value": p_value,
                "interpretation": (
                    "Simulated-data demonstration only; not a biological result."
                ),
            }
        )
    return pd.DataFrame(records)


def save_bar_plot(
    data: pd.DataFrame,
    y: str,
    ylabel: str,
    title: str,
    output_path: Path,
) -> None:
    """Save a branch-by-saccharide mean plot with SD error bars."""
    sns.set_theme(style="whitegrid", context="notebook", font_scale=1.05)
    figure, axis = plt.subplots(figsize=(10, 5.8))
    sns.barplot(
        data=data,
        x="branch",
        y=y,
        hue="saccharide",
        order=BRANCH_ORDER,
        hue_order=SACCHARIDE_ORDER,
        palette=PALETTE,
        errorbar="sd",
        capsize=0.08,
        ax=axis,
    )
    axis.set_title(title, weight="bold", fontsize=16, pad=12)
    axis.set_xlabel("")
    axis.set_ylabel(ylabel, fontsize=11)
    axis.tick_params(axis="both", labelsize=10)
    axis.legend(
        title="Saccharide",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        fontsize=9,
        title_fontsize=10,
    )
    sns.despine(figure)
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def save_scatter_plot(data: pd.DataFrame, output_path: Path) -> None:
    """Plot modification extent against simulated IgE binding."""
    sns.set_theme(style="whitegrid", context="notebook", font_scale=1.05)
    figure, axis = plt.subplots(figsize=(9, 5.8))
    sns.scatterplot(
        data=data,
        x="modification_extent_pct",
        y="ige_binding_pct",
        hue="saccharide",
        style="branch",
        hue_order=SACCHARIDE_ORDER,
        palette=PALETTE,
        s=85,
        alpha=0.85,
        ax=axis,
    )
    axis.set_title(
        "Modification Extent vs. IgE Binding (Simulated Data)",
        weight="bold",
        fontsize=16,
        pad=12,
    )
    axis.set_xlabel("Modification extent vs. branch control (%)", fontsize=11)
    axis.set_ylabel("IgE binding (% of native OVA)", fontsize=11)
    axis.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=True,
        fontsize=9,
    )
    sns.despine(figure)
    figure.tight_layout()
    figure.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the OVA analysis prototype.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/example_simulated_ova_data.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    analyzed = compute_outcomes(pd.read_csv(args.input))
    table_dir = args.output_dir / "tables"
    figure_dir = args.output_dir / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    analyzed.to_csv(table_dir / "analyzed_measurements.csv", index=False)
    summarize_treatments(analyzed).to_csv(
        table_dir / "treatment_summary.csv", index=False
    )
    run_anova(analyzed).to_csv(table_dir / "anova_results.csv", index=False)
    save_bar_plot(
        analyzed,
        "modification_extent_pct",
        "Modification extent vs. branch control (%)",
        "OVA Modification Extent by Saccharide (Simulated Data)",
        figure_dir / "modification_extent.png",
    )
    save_bar_plot(
        analyzed,
        "ige_binding_pct",
        "IgE binding (% of native OVA)",
        "OVA IgE Binding by Saccharide (Simulated Data)",
        figure_dir / "ige_binding.png",
    )
    save_scatter_plot(analyzed, figure_dir / "modification_vs_ige.png")
    print(f"Analyzed {len(analyzed)} rows; outputs written to {args.output_dir}")
    print("Simulated-data demonstration only; no biological conclusion.")


if __name__ == "__main__":
    main()
