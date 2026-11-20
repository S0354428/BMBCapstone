"""
load_bradford.py  --  Figure 2, step 1.

Loads the BSA standards and the Plate 2 sample plate and maps every well to
its sample, condition, dilution and replicate. No fitting or plotting here.

Run:     python load_bradford.py
Import:  from load_bradford import load_standards, load_plate
"""

import csv
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# Paths  (Src/ and Data/ are siblings inside BMBCap/)
# --------------------------------------------------------------------------
DATA_DIR = Path(__file__).resolve().parent.parent / "Data"
STANDARDS_FILE = DATA_DIR / "bradford_exploratory_selected_triplicates.csv"
PLATE_FILE = DATA_DIR / "Bradford__Plate_2_-_Sheet1_.csv"

# --------------------------------------------------------------------------
# PLATE LAYOUT  -- edit here if any of this is wrong
# --------------------------------------------------------------------------
# Rows are replicates, columns are solutions (same convention as the standards).
CONDITION_ROWS = {"Induced": "ABC", "Uninduced": "DEF"}

# Each group takes three adjacent columns, left to right = undiluted, 1:2, 1:10
GROUP_START_COL = {"Crude": 1, "E1": 4, "E2": 7, "FT": 10}
DILUTIONS = [("undiluted", 1), ("1:2", 2), ("1:10", 10)]

# Wash wells (our reading of the handwritten sketch): three wells across
# columns 1-3 of rows G and H, one dilution only. Confirm before relying on them.
WASH_WELLS = {"W1 Induced": "G", "W1 Uninduced": "H"}

ROW_LETTERS = "ABCDEFGH"


# --------------------------------------------------------------------------
# Standards
# --------------------------------------------------------------------------
def load_standards():
    """
    Return a DataFrame with columns:
        BSA_mg_per_mL, rep1, rep2, rep3, mean
    Only the 11 data rows are read; the summary table below them is ignored.
    """
    raw = pd.read_csv(STANDARDS_FILE, nrows=11)
    df = pd.DataFrame({
        "BSA_mg_per_mL": raw["BSA_mg_per_mL"].astype(float),
        "rep1": raw["replicate_1_A595"].astype(float),
        "rep2": raw["replicate_2_A595"].astype(float),
        "rep3": raw["replicate_3_A595"].astype(float),
    })
    df["mean"] = df[["rep1", "rep2", "rep3"]].mean(axis=1)
    return df


# --------------------------------------------------------------------------
# Plate
# --------------------------------------------------------------------------
def read_plate_grid(path=PLATE_FILE):
    """Find the 8 x 12 block of absorbances in the reader export."""
    with open(path, newline="") as f:
        rows = list(csv.reader(f))

    header_idx = next(
        i for i, r in enumerate(rows)
        if len(r) > 13 and [c.strip() for c in r[2:14]] == [str(n) for n in range(1, 13)]
    )
    grid = {}
    for r in rows[header_idx + 1: header_idx + 9]:
        letter = r[1].strip()
        grid[letter] = [float(v) for v in r[2:14]]
    return pd.DataFrame(grid, index=range(1, 13)).T  # rows A-H, columns 1-12


def load_plate():
    """
    Return a tidy DataFrame, one row per measured well:
        well, condition, group, dilution, dilution_factor, replicate, A595
    Wells not in the layout (empty wells) are left out.
    """
    grid = read_plate_grid()
    records = []

    for condition, letters in CONDITION_ROWS.items():
        for group, start in GROUP_START_COL.items():
            for k, (dil_name, dil_factor) in enumerate(DILUTIONS):
                col = start + k
                for rep, letter in enumerate(letters, start=1):
                    records.append({
                        "well": f"{letter}{col}",
                        "condition": condition,
                        "group": group,
                        "dilution": dil_name,
                        "dilution_factor": dil_factor,
                        "replicate": rep,
                        "A595": grid.loc[letter, col],
                    })

    for name, letter in WASH_WELLS.items():
        condition = name.split()[1]
        for rep, col in enumerate((1, 2, 3), start=1):
            records.append({
                "well": f"{letter}{col}",
                "condition": condition,
                "group": "W1",
                "dilution": "undiluted",
                "dilution_factor": 1,
                "replicate": rep,
                "A595": grid.loc[letter, col],
            })

    return pd.DataFrame(records)


# --------------------------------------------------------------------------
# Sanity check
# --------------------------------------------------------------------------
def sanity_check():
    std = load_standards()
    plate = load_plate()
    grid = read_plate_grid()

    print("=== Standards ===")
    print(std.round(3).to_string(index=False))
    print(f"\nBlank (0 mg/mL mean): {std['mean'].iloc[0]:.3f}")

    print("\n=== Plate 2 raw grid (A595) ===")
    print(grid.round(3).to_string())

    print("\n=== Layout check: wells -> sample ===")
    print("(each line = one group at one dilution; wells listed in replicate order)")
    summary = (plate.groupby(["condition", "group", "dilution"], sort=False)["well"]
               .apply(lambda s: ", ".join(s)))
    print(summary.to_string())

    used = set(plate["well"])
    unused = [f"{r}{c}" for r in ROW_LETTERS for c in range(1, 13) if f"{r}{c}" not in used]
    print(f"\nWells used: {len(used)}   Wells not assigned: {len(unused)}")
    print("Unassigned wells read (min / max):",
          f"{min(grid.loc[w[0], int(w[1:])] for w in unused):.3f} /",
          f"{max(grid.loc[w[0], int(w[1:])] for w in unused):.3f}")

    print("\n=== Mean A595 by condition / group / dilution ===")
    means = (plate.groupby(["condition", "group", "dilution"], sort=False)["A595"]
             .agg(["mean", "std"]).round(3))
    print(means.to_string())


if __name__ == "__main__":
    sanity_check()
