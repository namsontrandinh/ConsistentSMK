#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_rkcl_atkc.py

Generates 3 figures from an RKCL/ATKC(/MultiStream) result CSV
(e.g. facebook_all_algos.csv):

  1. f(S) vs Budget B                         -> fig1_value_vs_B.png
  2. Consistency C vs Budget B (vs bound 2q)  -> fig2_consistency_vs_B.png
  3. Internal repair recourse: actual vs bound -> fig3_recourse_vs_B.png

All chart text (titles, axis labels, legends) is in English, for direct
use in the CSONET 2026 camera-ready. Code comments stay in Vietnamese
(they never appear in the figures).

Each figure has one subplot per distinct delta value in the data (the
theoretical bounds differ by delta, so mixing them on one axis would
be misleading).

Usage:
    python3 plot_rkcl_atkc.py facebook_all_algos.csv [output_dir]
"""

import sys
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # no display needed, just save files
import matplotlib.pyplot as plt


def load_and_unify(csv_path: str) -> pd.DataFrame:
    """Read the CSV and build unified 'delta', 'q', 'N',
    'successful_exchanges', 'theoretical_recourse' columns shared by
    RKCL-IM and ATKC (each algorithm stores these under its own
    prefixed columns: atkc_* / rkcl_*). MultiStream rows (if present)
    carry no delta of their own; they are tagged with the delta of the
    matching cost dataset so fig1 can overlay them fairly."""
    df = pd.read_csv(csv_path)

    df["is_rkcl"] = df["algo"].str.upper().str.startswith("RKCL")
    df["is_atkc"] = df["algo"].str.upper() == "ATKC"
    df["is_multistream"] = df["algo"].str.upper() == "MULTISTREAM"

    df["delta"] = df["rkcl_delta_used"].where(df["is_rkcl"], df["atkc_delta"])
    # MultiStream has no delta parameter; tag it by which cost dataset
    # it ran on, matching the RKCL/ATKC delta that used the same
    # dataset (fb_cost1020 <-> delta=0.1, fb_cost4860 <-> delta=0.33).
    df.loc[df["is_multistream"] & df["graph"].str.contains("1020"), "delta"] = 0.1
    df.loc[df["is_multistream"] & df["graph"].str.contains("4860"), "delta"] = 0.33

    df["q"] = df["rkcl_q"].where(df["is_rkcl"], df["atkc_q"])
    df["N"] = df["rkcl_N"].where(df["is_rkcl"], df["atkc_N"])
    df["successful_exchanges"] = df["rkcl_successful_exchanges"].where(
        df["is_rkcl"], df["atkc_successful_exchanges"]
    )

    # Theoretical internal repair recourse bound (Theorem 1 / Theorem 2):
    #   RKCL: N*(q+1)
    #   ATKC: (J+1)*N*(q+1)
    rkcl_bound = df["N"] * (df["q"] + 1)
    atkc_bound = (df["atkc_J"] + 1) * df["N"] * (df["q"] + 1)
    df["theoretical_recourse"] = rkcl_bound.where(df["is_rkcl"], atkc_bound)

    # Short legend label
    df["algo_label"] = df["algo"].replace({"RKCL-IM": "RKCL"})

    return df


def plot_value_vs_B(df: pd.DataFrame, out_dir: str):
    deltas = sorted(df["delta"].unique())
    fig, axes = plt.subplots(1, len(deltas), figsize=(6 * len(deltas), 5), sharey=True)
    if len(deltas) == 1:
        axes = [axes]

    for ax, delta in zip(axes, deltas):
        sub = df[df["delta"] == delta]
        for algo in sub["algo_label"].unique():
            s = sub[sub["algo_label"] == algo].sort_values("B")
            ax.plot(s["B"], s["f_value"], marker="o", label=algo)
        ax.set_title(f"δ = {delta}")
        ax.set_xlabel("Budget B")
        ax.grid(alpha=0.3)
        ax.legend()

    axes[0].set_ylabel("f(S) — IC spread")
    fig.suptitle("Solution value f(S) vs Budget B")
    fig.tight_layout()
    path = os.path.join(out_dir, "fig1_value_vs_B.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_consistency_vs_B(df: pd.DataFrame, out_dir: str):
    deltas = sorted(df["delta"].unique())
    fig, axes = plt.subplots(1, len(deltas), figsize=(6 * len(deltas), 5), sharey=False)
    if len(deltas) == 1:
        axes = [axes]

    for ax, delta in zip(axes, deltas):
        sub = df[df["delta"] == delta]
        q = sub["q"].iloc[0]
        bound_2q = 2 * q

        for algo in sub["algo_label"].unique():
            s = sub[sub["algo_label"] == algo].sort_values("B")
            if algo.lower() == "multistream" or s["C"].isna().all():
                continue  # MultiStream has no consistency data (not measured)
            ax.plot(s["B"], s["C"], marker="o", label=f"{algo} (measured)")

        ax.axhline(bound_2q, color="red", linestyle="--", label=f"Theoretical bound 2q={bound_2q}")
        ax.set_title(f"δ = {delta} (q={q})")
        ax.set_xlabel("Budget B")
        ax.set_ylabel("C (max |Yt△Yt-1|)")
        ax.grid(alpha=0.3)
        ax.legend()

    fig.suptitle("Consistency: measured vs theoretical bound 2q")
    fig.tight_layout()
    path = os.path.join(out_dir, "fig2_consistency_vs_B.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def plot_recourse_vs_B(df: pd.DataFrame, out_dir: str):
    deltas = sorted(df["delta"].unique())
    fig, axes = plt.subplots(1, len(deltas), figsize=(6 * len(deltas), 5), sharey=False)
    if len(deltas) == 1:
        axes = [axes]

    for ax, delta in zip(axes, deltas):
        sub = df[df["delta"] == delta]

        for algo in sub["algo_label"].unique():
            s = sub[sub["algo_label"] == algo].sort_values("B")
            if algo.lower() == "multistream" or s["successful_exchanges"].isna().all():
                continue  # MultiStream has no recourse data (not measured)
            ax.plot(s["B"], s["successful_exchanges"], marker="o",
                     label=f"{algo} (measured)")
            ax.plot(s["B"], s["theoretical_recourse"], marker="x", linestyle="--",
                     label=f"{algo} (theoretical bound)")

        ax.set_title(f"δ = {delta}")
        ax.set_xlabel("Budget B")
        ax.set_ylabel("Number of successful exchanges")
        # Linear scale (not log): RKCL can legitimately have 0
        # successful exchanges, which a log scale cannot display.
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("Internal repair recourse: measured vs theoretical bound")
    fig.tight_layout()
    path = os.path.join(out_dir, "fig3_recourse_vs_B.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plot_rkcl_atkc.py <csv_path> [output_dir]")
        sys.exit(1)

    csv_path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) >= 3 else "."
    os.makedirs(out_dir, exist_ok=True)

    df = load_and_unify(csv_path)

    plot_value_vs_B(df, out_dir)
    plot_consistency_vs_B(df, out_dir)
    plot_recourse_vs_B(df, out_dir)

    print("\nDone. 3 PNG files created in:", os.path.abspath(out_dir))


if __name__ == "__main__":
    main()
