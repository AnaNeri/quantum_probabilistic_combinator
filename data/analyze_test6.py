#!/usr/bin/env python3
"""Compare Test 6's exact vs. simulated (Monte Carlo) break-even results.

Test 6 answers the same question test2 raises (is bit-flip correction worth
it, or is it better to just build the circuit from better gates?) two ways:
  - test6exact:  closed-form density matrices (quantumChoiceExact, no
                 sampling), break-even found by bisection
  - test6sim:    Monte Carlo trials (quantumChoice/gateWithBFInQubit), same
                 mechanism as test2/runTest2Sweep, break-even found by
                 linear interpolation over a decreaseProb grid
This script reads both sets of CSVs (produced by `tests test6 <id>`,
`tests test6exact <id>`, or `tests test6sim <n> <id>`) and reports how the
simulated break-even estimate compares to the exact one, per p.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DATA_DIR = Path(__file__).parent


def plot_comparison(idOfTest, curves, sim_points, breakevens):
    ps = sorted(sim_points["prob_error"].unique()) or sorted(curves["prob_error"].unique())[:2]
    fig, axes = plt.subplots(1, len(ps), figsize=(5.5 * len(ps), 4.5), sharey=True)
    if len(ps) == 1:
        axes = [axes]
    for axis, p in zip(axes, ps):
        curve = curves[curves["prob_error"] == p].sort_values("decrease_prob")
        axis.plot(curve["decrease_prob"], curve["1_corr"], color="tab:blue", label="correction (exact)")
        axis.plot(curve["decrease_prob"], curve["1_reduced_only"], color="tab:orange", label="better_gates (exact)")

        points = sim_points[sim_points["prob_error"] == p].sort_values("decrease_prob")
        if not points.empty:
            axis.scatter(points["decrease_prob"], points["1_corr"], color="tab:blue", marker="o",
                         edgecolor="black", zorder=3, label="correction (sim)")
            axis.scatter(points["decrease_prob"], points["1_reduced_only"], color="tab:orange", marker="s",
                         edgecolor="black", zorder=3, label="better_gates (sim)")

        match = breakevens[breakevens["prob_error"] == p]
        if not match.empty and pd.notna(match["breakeven_decrease_prob_exact"].iloc[0]):
            d_star = match["breakeven_decrease_prob_exact"].iloc[0]
            axis.axvline(d_star, color="gray", linestyle="--", label=f"break-even (exact) d={d_star:.3f}")

        axis.set_title(f"p = {p:g}")
        axis.set_xlabel("decreaseProb")
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=7)
    axes[0].set_ylabel("error rate ('1' population)")
    fig.suptitle("Test 6: correction vs. better-gates-only, exact curve vs. simulated points")
    fig.tight_layout()
    plot_path = DATA_DIR / f"test6_comparison_{idOfTest}.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    return plot_path


def main():
    idOfTest = sys.argv[1] if len(sys.argv) > 1 else "*"

    exact_files = sorted(DATA_DIR.glob(f"out_test6_breakeven_exact_{idOfTest}.csv"))
    sim_files = sorted(DATA_DIR.glob(f"out_test6_breakeven_sim_{idOfTest}.csv"))
    curve_files = sorted(DATA_DIR.glob(f"out_test6_exact_{idOfTest}.csv"))
    sim_point_files = sorted(DATA_DIR.glob(f"out_test6_sim_{idOfTest}.csv"))
    if not exact_files:
        raise SystemExit(f"No out_test6_breakeven_exact_{idOfTest}.csv found. Run `tests test6exact <id>` first.")
    if not sim_files:
        raise SystemExit(f"No out_test6_breakeven_sim_{idOfTest}.csv found. Run `tests test6sim <n> <id>` first.")

    exact = pd.concat([pd.read_csv(p) for p in exact_files], ignore_index=True)
    sim = pd.concat([pd.read_csv(p) for p in sim_files], ignore_index=True)

    merged = pd.merge(sim, exact, on="prob_error", how="inner", suffixes=("_sim", "_exact"))
    merged["abs_diff"] = (merged["breakeven_decrease_prob_sim"] - merged["breakeven_decrease_prob_exact"]).abs()

    print(merged.to_string(index=False))
    print()
    print("Note: test6sim only sweeps a coarse decreaseProb grid ([0.0, 0.5, 1.0]) with a")
    print("small trial count, so it can only report a break-even point if the sign of")
    print("(1_corr - 1_reduced_only) actually changes between two adjacent grid points;")
    print("test6exact bisects to ~2^-60 precision over the full [0,1] range regardless of p.")

    if curve_files and sim_point_files:
        curves = pd.concat([pd.read_csv(p) for p in curve_files], ignore_index=True)
        sim_points = pd.concat([pd.read_csv(p) for p in sim_point_files], ignore_index=True)
        plot_path = plot_comparison(idOfTest, curves, sim_points, merged)
        print(f"\nGraphic saved to {plot_path}")
    else:
        print("\nSkipping plot: need both out_test6_exact_*.csv and out_test6_sim_*.csv to draw it.")


if __name__ == "__main__":
    main()
