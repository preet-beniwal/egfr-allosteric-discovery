"""
Phase 2 - Step 2.5
Analyze docking results:
  1. Rank all 516 ligands by affinity
  2. Split into actives vs decoys
  3. Compute ROC-AUC + enrichment factor
  4. Save the top 10 hits for Phase 3 (MD)
"""

import matplotlib
matplotlib.use('Agg')  # headless mode for WSL
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score

# ---------- LOAD ----------
df = pd.read_csv('data/results/docking_scores.csv')
df = df.dropna(subset=['affinity_kcal_mol'])
df = df.sort_values('affinity_kcal_mol').reset_index(drop=True)
print(f"Total ligands docked: {len(df)}")

# ---------- LABELS ----------
df['is_active'] = ~df['ligand'].str.startswith('decoy_')
n_actives = int(df['is_active'].sum())
n_decoys  = int((~df['is_active']).sum())
print(f"Actives: {n_actives}")
print(f"Decoys:  {n_decoys}")

# ---------- DISTRIBUTION ----------
print("\n=== Score Statistics ===")
print(f"Actives  mean: {df[df.is_active].affinity_kcal_mol.mean():.2f} kcal/mol")
print(f"Decoys   mean: {df[~df.is_active].affinity_kcal_mol.mean():.2f} kcal/mol")
print(f"Actives  best: {df[df.is_active].affinity_kcal_mol.min():.2f} kcal/mol")
print(f"Decoys   best: {df[~df.is_active].affinity_kcal_mol.min():.2f} kcal/mol")

# ---------- ROC-AUC ----------
y_true  = df['is_active'].astype(int).values
y_score = -df['affinity_kcal_mol'].values  # higher = better binder

auc_val = roc_auc_score(y_true, y_score)
print(f"\n=== ROC-AUC: {auc_val:.3f} ===")
if auc_val >= 0.85:
    print("PUBLICATION QUALITY")
elif auc_val >= 0.75:
    print("EXCELLENT discrimination")
elif auc_val >= 0.65:
    print("ACCEPTABLE discrimination")
else:
    print("WEAK - needs investigation")

# ---------- PLOT ----------
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

axes[0].hist(df[df.is_active].affinity_kcal_mol,  bins=20,
             alpha=0.7, label='Actives', color='green')
axes[0].hist(df[~df.is_active].affinity_kcal_mol, bins=20,
             alpha=0.5, label='Decoys', color='red')
axes[0].set_xlabel('Binding affinity (kcal/mol)')
axes[0].set_ylabel('Count')
axes[0].set_title('Score Distribution: Actives vs Decoys')
axes[0].legend()

fpr, tpr, _ = roc_curve(y_true, y_score)
axes[1].plot(fpr, tpr, color='navy', lw=2, label=f'ROC (AUC = {auc_val:.3f})')
axes[1].plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
axes[1].set_xlabel('False Positive Rate')
axes[1].set_ylabel('True Positive Rate')
axes[1].set_title('ROC Curve - Docking Benchmark')
axes[1].legend()

plt.tight_layout()
plt.savefig('data/results/roc_analysis.png', dpi=150)
print("Plot saved: data/results/roc_analysis.png")

# ---------- ENRICHMENT FACTOR (Top 1%) ----------
top_1pct = max(1, int(len(df) * 0.01))
top_hits = df.head(top_1pct)
ef1 = (top_hits.is_active.sum() / top_1pct) / (n_actives / len(df))
print(f"\n=== EF1% (Top 1% enrichment): {ef1:.2f} ===")
print(f"Top {top_1pct} ligands contained {top_hits.is_active.sum()} actives")

# ---------- TOP 10 ----------
top10 = df.head(10)[['ligand', 'affinity_kcal_mol']]
top10.to_csv('data/results/top10_hits.csv', index=False)
print(f"\n=== TOP 10 HITS ===")
print(top10.to_string(index=False))
print("\nSaved -> data/results/top10_hits.csv")
