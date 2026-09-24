#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_rkcl_atkc.py

Ve 3 nhom bieu do tu file CSV ket qua RKCL/ATKC (vd facebook_rkcl_atkc_final.csv):

  1. f(S) theo Budget B          -> fig1_value_vs_B.png
  2. Consistency C theo Budget B (so voi can ly thuyet 2q) -> fig2_consistency_vs_B.png
  3. Repair recourse thuc te so voi can ly thuyet           -> fig3_recourse_vs_B.png

Moi bieu do ve rieng cho tung gia tri delta (vi cac duong delta khac nhau
co thang do/can ly thuyet khac nhau, gop chung se roi).

Usage:
    python3 plot_rkcl_atkc.py facebook_rkcl_atkc_final.csv [output_dir]
"""

import sys
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # khong can man hinh, chi luu file
import matplotlib.pyplot as plt


def load_and_unify(csv_path: str) -> pd.DataFrame:
    """Doc CSV va tao cac cot 'delta', 'q', 'N', 'successful_exchanges',
    'theoretical_recourse' dung chung cho ca RKCL-IM va ATKC, vi hai
    thuat toan luu cac gia tri nay o cac cot rieng (atkc_* / rkcl_*)."""
    df = pd.read_csv(csv_path)

    df["is_rkcl"] = df["algo"].str.upper().str.startswith("RKCL")
    df["is_atkc"] = df["algo"].str.upper() == "ATKC"

    df["delta"] = df["rkcl_delta_used"].where(df["is_rkcl"], df["atkc_delta"])
    df["q"] = df["rkcl_q"].where(df["is_rkcl"], df["atkc_q"])
    df["N"] = df["rkcl_N"].where(df["is_rkcl"], df["atkc_N"])
    df["successful_exchanges"] = df["rkcl_successful_exchanges"].where(
        df["is_rkcl"], df["atkc_successful_exchanges"]
    )

    # Can ly thuyet cho internal repair recourse (Theorem 1 / Theorem 2):
    #   RKCL: N*(q+1)
    #   ATKC: (J+1)*N*(q+1)
    rkcl_bound = df["N"] * (df["q"] + 1)
    atkc_bound = (df["atkc_J"] + 1) * df["N"] * (df["q"] + 1)
    df["theoretical_recourse"] = rkcl_bound.where(df["is_rkcl"], atkc_bound)

    # Nhan gon cho legend
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
    fig.suptitle("Giá trị nghiệm f(S) theo Budget B")
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
            ax.plot(s["B"], s["C"], marker="o", label=f"{algo} (C thực tế)")

        ax.axhline(bound_2q, color="red", linestyle="--", label=f"Cận lý thuyết 2q={bound_2q}")
        ax.set_title(f"δ = {delta} (q={q})")
        ax.set_xlabel("Budget B")
        ax.set_ylabel("C (max |Yt△Yt-1|)")
        ax.grid(alpha=0.3)
        ax.legend()

    fig.suptitle("Consistency thực tế so với cận lý thuyết 2q")
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
            ax.plot(s["B"], s["successful_exchanges"], marker="o",
                     label=f"{algo} (thực tế)")
            ax.plot(s["B"], s["theoretical_recourse"], marker="x", linestyle="--",
                     label=f"{algo} (cận lý thuyết)")

        ax.set_title(f"δ = {delta}")
        ax.set_xlabel("Budget B")
        ax.set_ylabel("Số successful exchanges")
        ax.set_yscale("log")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("Internal repair recourse: thực tế vs cận lý thuyết (thang log)")
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

    print("\nDone. 3 file PNG da duoc tao trong:", os.path.abspath(out_dir))


if __name__ == "__main__":
    main()
