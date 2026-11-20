"""
fit_standards.py  --  Figure 2, step 2.

Blank-corrects the BSA standards and fits the standard curve.
No sample quantification and no plotting here.

Method (follows the lab handout):
    1. average the three replicates for each standard and the blank
    2. subtract the mean blank (0 mg/mL standard) from every standard
    3. fit corrected A595 (y) against [BSA] in mg/mL (x), straight line,
       free intercept:  A = m * C + b

Run:     python fit_standards.py
Import:  from fit_standards import fit_standards
"""

import numpy as np
from scipy import stats

from load_bradford import load_standards

R2_CUTOFF = 0.9800  # handout: R2 below this means repeat the assay


def fit_standards():
    """
    Return a dict:
        x            BSA concentrations (mg/mL)
        y            blank-corrected mean A595
        blank        mean A595 of the 0 mg/mL standard
        slope, intercept, r2, slope_se, intercept_se
        r2_all_reps  R2 if all 33 individual readings are fit instead (reference)
    """
    std = load_standards()
    x = std["BSA_mg_per_mL"].to_numpy()
    blank = std["mean"].iloc[0]
    y = std["mean"].to_numpy() - blank

    fit = stats.linregress(x, y)

    # reference only: same line type using every individual reading
    reps = std[["rep1", "rep2", "rep3"]].to_numpy() - blank
    fit_all = stats.linregress(np.repeat(x, 3), reps.flatten())

    return {
        "x": x,
        "y": y,
        "blank": blank,
        "slope": fit.slope,
        "intercept": fit.intercept,
        "r2": fit.rvalue ** 2,
        "slope_se": fit.stderr,
        "intercept_se": fit.intercept_stderr,
        "r2_all_reps": fit_all.rvalue ** 2,
    }


def equation_text(f, digits=3):
    """e.g.  y = 1.340x + 0.202"""
    sign = "+" if f["intercept"] >= 0 else "-"
    return f"y = {f['slope']:.{digits}f}x {sign} {abs(f['intercept']):.{digits}f}"


def sanity_check():
    f = fit_standards()

    print("=== Blank-corrected standards ===")
    print(f"Blank (mean A595 at 0 mg/mL): {f['blank']:.3f}")
    print(f"{'BSA (mg/mL)':>12} {'corrected A595':>15} {'fitted':>8} {'residual':>9}")
    for xi, yi in zip(f["x"], f["y"]):
        fitted = f["slope"] * xi + f["intercept"]
        print(f"{xi:>12.1f} {yi:>15.3f} {fitted:>8.3f} {yi - fitted:>9.3f}")

    print("\n=== Linear fit (corrected A595 vs mg/mL, free intercept) ===")
    print(f"Equation : {equation_text(f)}   (A = m*C + b)")
    print(f"Slope    : {f['slope']:.4f} +/- {f['slope_se']:.4f}")
    print(f"Intercept: {f['intercept']:.4f} +/- {f['intercept_se']:.4f}")
    print(f"R2 (11 means)       : {f['r2']:.4f}")
    print(f"R2 (all 33 readings): {f['r2_all_reps']:.4f}   (reference)")
    verdict = "meets" if f["r2"] >= R2_CUTOFF else "is below"
    print(f"R2 {verdict} the handout cutoff of {R2_CUTOFF:.4f}")
    print(f"\nHighest corrected standard: {f['y'].max():.3f}  "
          "(sample wells above this are outside the calibration range)")


if __name__ == "__main__":
    sanity_check()
