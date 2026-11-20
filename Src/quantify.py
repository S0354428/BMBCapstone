"""
quantify.py  --  Figure 2, step 3.

Converts Plate 2 absorbances to protein concentration (mg/mL) using the
standard curve from fit_standards.py, then picks one dilution per sample.
No plotting here.

Per well:      C_original = ((A595 - blank) - b) / m  *  dilution_factor
Per dilution:  absorbance is averaged over the 3 replicate wells and given a status:
    below detection : corrected A <= intercept + 1 SE   (curve can't resolve it)
    linear          : above that and <= LINEAR_LIMIT_A  (~0.5 mg/mL standard)
    high            : above the linear part, still inside the standards' range
    out of range    : above the highest standard

Rule for the reported value (one per condition x group):
    1. among 'linear' dilutions, use the one with the strongest signal
    2. otherwise, among 'high' dilutions, use the MOST dilute one (least saturated);
       reported as an approximate value
    3. otherwise: below detection (reported as 0)
Every well is kept; nothing is dropped as an outlier. Error = SD of the 3
replicate concentrations at the chosen dilution.

Run:     python quantify.py
Import:  from quantify import quantify
Writes:  ../Results/concentrations_all_dilutions.csv, ../Results/concentrations_reported.csv
"""

from pathlib import Path

import numpy as np
import pandas as pd

from fit_standards import fit_standards
from load_bradford import load_plate

RESULTS_DIR = Path(__file__).resolve().parent.parent / "Results"

# Corrected A595 above which the standards flatten out (~0.5 mg/mL standard).
LINEAR_LIMIT_A = 0.90

STATUS_NOTE = {
    "linear": "",
    "high": "approximate: above the linear range of the standard curve",
    "below detection": "below detection",
    "out of range": "above the highest standard",
}


def status_of(corrected_a, fit):
    floor = fit["intercept"] + fit["intercept_se"]
    if corrected_a <= floor:
        return "below detection"
    if corrected_a <= LINEAR_LIMIT_A:
        return "linear"
    if corrected_a <= fit["y"].max():
        return "high"
    return "out of range"


def quantify():
    """Return (all_dilutions_df, reported_df)."""
    fit = fit_standards()
    plate = load_plate()
    m, b, blank = fit["slope"], fit["intercept"], fit["blank"]

    plate["corrected_A"] = plate["A595"] - blank
    plate["conc_mg_mL"] = (plate["corrected_A"] - b) / m * plate["dilution_factor"]

    rows = []
    keys = ["condition", "group", "dilution", "dilution_factor"]
    for (cond, grp, dil, fac), g in plate.groupby(keys, sort=False):
        mean_corr = g["corrected_A"].mean()
        rows.append({
            "condition": cond, "group": grp, "dilution": dil, "dilution_factor": fac,
            "mean_A595": g["A595"].mean(),
            "mean_corrected_A": mean_corr,
            "status": status_of(mean_corr, fit),
            "conc_mean_mg_mL": g["conc_mg_mL"].mean(),
            "conc_sd_mg_mL": g["conc_mg_mL"].std(ddof=1),
        })
    all_df = pd.DataFrame(rows)

    reported = []
    for (cond, grp), g in all_df.groupby(["condition", "group"], sort=False):
        linear = g[g["status"] == "linear"]
        high = g[g["status"] == "high"]
        if len(linear):
            pick = linear.loc[linear["mean_corrected_A"].idxmax()]
            conc, sd, status = pick["conc_mean_mg_mL"], pick["conc_sd_mg_mL"], "linear"
        elif len(high):
            pick = high.loc[high["dilution_factor"].idxmax()]
            conc, sd, status = pick["conc_mean_mg_mL"], pick["conc_sd_mg_mL"], "high"
        else:
            pick = g.loc[g["dilution_factor"].idxmin()]
            conc, sd, status = 0.0, 0.0, "below detection"
        reported.append({
            "condition": cond, "group": grp,
            "dilution_used": pick["dilution"] if status != "below detection" else "n/a",
            "conc_mg_mL": conc, "sd_mg_mL": sd,
            "status": status, "note": STATUS_NOTE[status],
        })
    rep_df = pd.DataFrame(reported)

    RESULTS_DIR.mkdir(exist_ok=True)
    all_df.to_csv(RESULTS_DIR / "concentrations_all_dilutions.csv", index=False)
    rep_df.to_csv(RESULTS_DIR / "concentrations_reported.csv", index=False)
    return all_df, rep_df


def sanity_check():
    fit = fit_standards()
    all_df, rep_df = quantify()
    floor = fit["intercept"] + fit["intercept_se"]

    print(f"Curve: A = {fit['slope']:.4f}*C + {fit['intercept']:.4f}   blank = {fit['blank']:.3f}")
    print(f"Detection floor (corrected A): {floor:.3f}   "
          f"linear limit: {LINEAR_LIMIT_A:.2f}   highest standard: {fit['y'].max():.3f}\n")

    print("=== Every sample at every dilution ===")
    show = all_df.copy()
    for c in ("mean_A595", "mean_corrected_A", "conc_mean_mg_mL", "conc_sd_mg_mL"):
        show[c] = show[c].round(3)
    print(show.drop(columns="dilution_factor").to_string(index=False))

    print("\n=== Reported values (one per sample) ===")
    out = rep_df.copy()
    out["conc_mg_mL"] = out["conc_mg_mL"].round(2)
    out["sd_mg_mL"] = out["sd_mg_mL"].round(2)
    print(out.to_string(index=False))
    print(f"\nSaved CSVs to {RESULTS_DIR}")


if __name__ == "__main__":
    sanity_check()
