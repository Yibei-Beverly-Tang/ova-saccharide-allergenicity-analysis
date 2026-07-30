"""Generate reproducible simulated data for the OVA analysis prototype.

These values are not experimental measurements and must not be reported as
biological findings.
"""

from pathlib import Path

import numpy as np
import pandas as pd


OUTPUT_PATH = Path("data/example_simulated_ova_data.csv")
RANDOM_SEED = 2026
REPLICATES = 5

DESIGN = {
    "Non-enzymatic": {
        "Control": (96, 98, 110, 48),
        "Lactose": (64, 72, 145, 40),
        "Mannotriose": (71, 79, 175, 38),
        "Dextran": (82, 75, 280, 34),
    },
    "TG-assisted": {
        "Control": (88, 94, 165, 44),
        "Lactose": (52, 61, 195, 35),
        "Mannotriose": (60, 69, 225, 33),
        "Dextran": (68, 65, 345, 29),
    },
    "Tyr/CA-assisted": {
        "Control": (92, 96, 145, 46),
        "Lactose": (70, 69, 180, 38),
        "Dextran": (80, 67, 315, 31),
    },
}

STANDARD_DEVIATIONS = {
    "free_amino_groups_pct": 3.0,
    "ige_binding_pct": 4.0,
    "particle_size_nm": 12.0,
    "surface_hydrophobicity_au": 3.0,
}


def build_dataset() -> pd.DataFrame:
    """Return a tidy, reproducible simulated dataset."""
    rng = np.random.default_rng(RANDOM_SEED)
    records = []
    metric_names = list(STANDARD_DEVIATIONS)

    for branch, treatments in DESIGN.items():
        for saccharide, means in treatments.items():
            treatment = (
                f"{branch} control"
                if saccharide == "Control"
                else f"{branch} {saccharide}"
            )
            for replicate in range(1, REPLICATES + 1):
                record = {
                    "sample_id": (
                        f"{branch[:3].upper()}-{saccharide[:3].upper()}-R{replicate}"
                    ),
                    "branch": branch,
                    "treatment": treatment,
                    "saccharide": saccharide,
                    "replicate": replicate,
                }
                for metric, mean in zip(metric_names, means):
                    value = rng.normal(mean, STANDARD_DEVIATIONS[metric])
                    record[metric] = round(max(value, 0), 2)
                records.append(record)

    return pd.DataFrame.from_records(records)


def main() -> None:
    data = build_dataset()
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUTPUT_PATH, index=False)
    print(f"Wrote {len(data)} simulated rows to {OUTPUT_PATH}")
    print("These values are simulated and are not experimental results.")


if __name__ == "__main__":
    main()
