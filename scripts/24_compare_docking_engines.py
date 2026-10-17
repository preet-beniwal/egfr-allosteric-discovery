#!/usr/bin/env python3
"""
Compare Vina and AutoDock4 docking results on the same 516-ligand library.

Inputs:
    data/results/                Vina output (existing)
    data/results_ad4/scores.tsv  AD4 output (from batch)

Outputs:
    data/results_ad4/roc_comparison.png
    data/results_ad4/engine_comparison.tsv

Metrics: ROC-AUC and enrichment factor at 1% (EF1).
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve

# ---- ground truth ----
KNOWN_ACTIVES = {
    "DDC4002",
    "EA1001",
    "EAI001",
    "EAI045",
    "EAI045_analog_1",
    "EAI045_analog_2",
    "EAI045_analog_3",
    "EAI045_analog_4",
    "EAI045_analog_5",
    "EAI045_analog_6",
    "EAI045_analog_7",
    "EAI045_analog_8",
    "EAI045_analog_9",
    "EAI045_analog_10",
    "JBJ-02-112-05",
    "JBJ-04-125-02",
}
def is_active(name):
    base = os.path.basename(name).replace(".pdbqt", "").replace("_out", "")
    return int(base in KNOWN_ACTIVES)

# ---- Vina scores ----
def load_vina_scores():
    """Read the Vina summary CSV written by 06_dock_all.py."""
    df = pd.read_csv("data/results/docking_scores.csv")
    # Normalize column names in case they differ
    df.columns = [c.strip().lower() for c in df.columns]
    if "ligand" not in df.columns or "affinity_kcal_mol" not in df.columns:
        raise ValueError(
            f"Unexpected Vina CSV columns: {df.columns.tolist()}"
        )
    df = df.rename(columns={"affinity_kcal_mol": "vina_score"})
    df["vina_score"] = pd.to_numeric(df["vina_score"], errors="coerce")
    return df.dropna(subset=["vina_score"])

# ---- AD4 scores ----
def load_ad4_scores():
    df = pd.read_csv("data/results_ad4/scores.tsv", sep="\t", header=None,
                     names=["ligand", "ad4_score"])
    df["ad4_score"] = pd.to_numeric(df["ad4_score"], errors="coerce")
    return df.dropna(subset=["ad4_score"])

def compute_roc(scores, higher_is_better):
    """Vina: more negative = better → lower_is_better.
       AD4:  more negative = better → lower_is_better (same direction)."""
    y_true = scores["active"].values
    raw = scores["score"].values
    # flip so higher score = better, which sklearn expects
    y_score = -raw if not higher_is_better else raw
    fpr, tpr, _ = roc_curve(y_true, y_score)
    auc = roc_auc_score(y_true, y_score)
    # EF1%
    k = max(1, int(0.01 * len(scores)))
    top_k = scores.sort_values("score").head(k)  # most negative
    ef1 = (top_k["active"].sum() / k) / (y_true.sum() / len(y_true))
    return auc, ef1, fpr, tpr

def main():
    print("[Comparison] Vina vs AutoDock4 on 516-ligand library\n")

    vina = load_vina_scores()
    ad4  = load_ad4_scores()
    print(f"  Vina ligands: {len(vina)}")
    print(f"  AD4  ligands: {len(ad4)}")

    vina["active"] = vina["ligand"].apply(is_active)
    ad4["active"]  = ad4["ligand"].apply(is_active)

    vina = vina.rename(columns={"vina_score": "score"})
    ad4  = ad4.rename(columns={"ad4_score": "score"})

    vina_auc, vina_ef1, vina_fpr, vina_tpr = compute_roc(vina, higher_is_better=False)
    ad4_auc,  ad4_ef1,  ad4_fpr,  ad4_tpr  = compute_roc(ad4,  higher_is_better=False)

    print(f"\n  Vina  ROC-AUC = {vina_auc:.3f}   EF1% = {vina_ef1:.2f}")
    print(f"  AD4   ROC-AUC = {ad4_auc:.3f}   EF1% = {ad4_ef1:.2f}")

    # ---- figure ----
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(vina_fpr, vina_tpr, color="#1f77b4", lw=2,
            label=f"Vina (AUC = {vina_auc:.3f})")
    ax.plot(ad4_fpr,  ad4_tpr,  color="#d62728", lw=2,
            label=f"AutoDock4 (AUC = {ad4_auc:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1)
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("Docking engine comparison on 516-ligand library")
    ax.legend(loc="lower right")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("data/results_ad4/roc_comparison.png", dpi=300)
    print("\n  wrote data/results_ad4/roc_comparison.png")

    # ---- summary table ----
    summary = pd.DataFrame([
        {"engine": "Vina 1.2.0",       "roc_auc": round(vina_auc, 3),
         "ef1_percent": round(vina_ef1, 2), "n_ligands": len(vina)},
        {"engine": "AutoDock4 (smina)", "roc_auc": round(ad4_auc, 3),
         "ef1_percent": round(ad4_ef1, 2), "n_ligands": len(ad4)},
    ])
    summary.to_csv("data/results_ad4/engine_comparison.tsv", sep="\t", index=False)
    print("  wrote data/results_ad4/engine_comparison.tsv")
    print("\n", summary.to_string(index=False))

if __name__ == "__main__":
    main()
