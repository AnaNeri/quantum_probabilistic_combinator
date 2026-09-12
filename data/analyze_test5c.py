#!/usr/bin/env python3
"""Aggregate Test 5c raw runs and save a comparison plot."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


DATA_DIR = Path(__file__).parent
RAW_PATTERN = "test5c_raw_Combined_1110_*.csv"


def main():
    raw_files = sorted(DATA_DIR.glob(RAW_PATTERN))
    if not raw_files:
        raise SystemExit("No combined Test 5c |1110> raw CSV files found. Run the Haskell test first.")

    frames = []
    for path in raw_files:
        frame = pd.read_csv(path)
        frame = frame.rename(
            columns={
                "target_qubit": "low_error_qubit",
                "target_probability": "low_error_probability",
                "fidelity_qubit_0": "fidelity_target_qubit_0",
            }
        )
        frames.append(frame)
    raw = pd.concat(frames, ignore_index=True)
    grouped = (
        raw.groupby(["noise_model", "probability", "low_error_qubit", "low_error_probability"])
        .agg(
            avg_fidelity=("fidelity_target_qubit_0", "mean"),
            std_fidelity=("fidelity_target_qubit_0", "std"),
            count=("fidelity_target_qubit_0", "size"),
        )
        .reset_index()
        .sort_values(["noise_model", "probability", "low_error_qubit"])
    )
    grouped["std_fidelity"] = grouped["std_fidelity"].fillna(0.0)
    clean_path = DATA_DIR / "test5c_clean.csv"
    grouped.to_csv(clean_path, index=False)

    fig, axis = plt.subplots(figsize=(8, 5))
    for low_error_qubit, target_frame in grouped.groupby("low_error_qubit"):
        target_frame = target_frame.sort_values("probability")
        axis.errorbar(
            target_frame["probability"],
            target_frame["avg_fidelity"],
            yerr=target_frame["std_fidelity"],
            marker="o",
            capsize=4,
            linewidth=2,
            label=f"low-error qubit {low_error_qubit}",
        )
    axis.set_xlabel("Other-gate fault probability")
    axis.set_ylabel("Average fidelity of target qubit 0")
    axis.set_xscale("log")
    axis.set_ylim(0, 1.05)
    axis.grid(True, alpha=0.3)
    axis.legend()
    fig.suptitle("Test 5c: combined reset and Z noise")
    fig.tight_layout()
    plot_path = DATA_DIR / "test5c_comparison.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")

    print(f"Clean summary saved to {clean_path}")
    print(f"Graphic saved to {plot_path}")
    print(grouped.to_string(index=False))


if __name__ == "__main__":
    main()
