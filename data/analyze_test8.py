#!/usr/bin/env python3
"""
Display script for Test 8 results (out_test8_spam_multi_*.csv).

What Test 8 tests
------------------
SPAM (State Preparation And Measurement) noise: a bit-flip (X) error applied
once, during state preparation only, before any circuit gates run. No noise
is applied on any other step of the circuit; measurement noise is not
modeled (the quantamorphism defers measurement).

The circuit is a single, self-contained 3-qubit Toffoli (CCX): qubits 0 and
1 are the controls, qubit 2 is the target. Two CCX applications (steps=2)
are run after the (possibly corrupted) state preparation.

Each of the 3 qubits gets its own independent bit-flip probability at state
preparation: two qubits get probability p, and one "protected" qubit gets
the lower probability p*0.1. Two cases are compared:
  - testA_target_protected:  target (qubit 2) gets p*0.1; qubits 0,1 get p
  - testB_control0_protected: control 0 (qubit 0) gets p*0.1; qubits 1,2 get p

Parameters used
---------------
  n_qubits    = 3
  steps       = 2   (CCX applied twice, noiseless after state prep)
  trials      = 1000 (Monte Carlo IO samples)
  probability sweep p = [0.0005, 0.005, 0.05]

Columns reported per (case, p)
-------------------------------
  qubit_probs               : per-qubit (index, probability) list used at state prep
  fidelity_exact_vs_ideal   : exact (Dist, branch-enumerated) state vs the noiseless ideal
  fidelity_mc_vs_ideal      : Monte Carlo (IO, sampled) state vs the noiseless ideal
  fidelity_exact_vs_mc      : exact state vs the Monte Carlo estimate (sampler check)
"""

import glob
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DATA_DIR = Path(__file__).parent
PATTERN = str(DATA_DIR / "out_test8_spam_multi_*.csv")


def load_latest():
    files = sorted(glob.glob(PATTERN))
    if not files:
        raise SystemExit("No out_test8_spam_multi_*.csv files found. Run `cabal run exe:tests -- test8` first.")
    latest = files[-1]
    return parse_csv(latest), latest


def parse_csv(path):
    """Parse the Test 8 CSV by hand: the `qubit_probs` field contains
    unquoted embedded commas (e.g. "[(0,5.0e-4),(1,5.0e-4),(2,5.0e-5)]"),
    which breaks a plain pandas.read_csv. The other 6 trailing fields
    (probability,steps,trials,fidelity_exact_vs_ideal,fidelity_mc_vs_ideal,
    fidelity_exact_vs_mc) and the leading `case` field are always single
    tokens, so anything in between is reassembled as qubit_probs.
    """
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    header = ["case", "qubit_probs", "probability", "steps", "trials",
              "fidelity_exact_vs_ideal", "fidelity_mc_vs_ideal", "fidelity_exact_vs_mc"]
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        tokens = line.split(",")
        case = tokens[0]
        tail = tokens[-6:]
        qubit_probs = ",".join(tokens[1:-6])
        rows.append([case, qubit_probs] + tail)
    df = pd.DataFrame(rows, columns=header)
    for col in ["probability", "steps", "trials", "fidelity_exact_vs_ideal", "fidelity_mc_vs_ideal", "fidelity_exact_vs_mc"]:
        df[col] = df[col].astype(float)
    df["steps"] = df["steps"].astype(int)
    df["trials"] = df["trials"].astype(int)
    return df


def main():
    df, source_file = load_latest()

    print(__doc__)
    print(f"Source file: {os.path.basename(source_file)}\n")

    display_cols = [
        "case",
        "qubit_probs",
        "probability",
        "steps",
        "trials",
        "fidelity_exact_vs_ideal",
        "fidelity_mc_vs_ideal",
        "fidelity_exact_vs_mc",
    ]
    table = df[display_cols].copy()
    for col in ["fidelity_exact_vs_ideal", "fidelity_mc_vs_ideal", "fidelity_exact_vs_mc"]:
        table[col] = table[col].map(lambda v: f"{v:.6f}")

    print(table.to_string(index=False))

    fig, axis = plt.subplots(figsize=(8, 5))
    for case_name, case_frame in df.groupby("case"):
        case_frame = case_frame.sort_values("probability")
        axis.plot(
            case_frame["probability"],
            case_frame["fidelity_exact_vs_ideal"],
            marker="o",
            linewidth=2,
            label=f"{case_name} (exact)",
        )
        axis.plot(
            case_frame["probability"],
            case_frame["fidelity_mc_vs_ideal"],
            marker="x",
            linestyle="--",
            linewidth=1.5,
            label=f"{case_name} (mc)",
        )
    axis.set_xlabel("State-prep bit-flip probability p")
    axis.set_ylabel("Fidelity vs ideal (noiseless)")
    axis.set_xscale("log")
    axis.set_ylim(0, 1.05)
    axis.grid(True, alpha=0.3)
    axis.legend()
    fig.suptitle("Test 8: SPAM noise, per-qubit probabilities p,p,p*0.1")
    fig.tight_layout()
    plot_path = DATA_DIR / "test8_comparison.png"
    fig.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"\nGraphic saved to {plot_path}")


if __name__ == "__main__":
    main()
