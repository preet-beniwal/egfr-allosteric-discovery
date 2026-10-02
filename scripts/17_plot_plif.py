"""Phase 3.7 — PLIF heatmap comparing actives vs decoys."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Load mapped results
df = pd.read_csv("data/md/plip/results/plip_interactions_EGFR.csv")

# Focus on key allosteric residues
KEY = ['LYS745', 'PHE856', 'LEU858', 'LEU718',
       'PHE723', 'ASP855', 'MET790', 'MET793']

df_key = df[df['Residue_EGFR'].isin(KEY)]

# Create binary matrix (any interaction per residue-ligand)
matrix = df_key.groupby(['Residue_EGFR', 'Ligand']).size().unstack(fill_value=0)

# Order ligands: actives first
LIG_ORDER = ['EAI045', 'EAI001', 'decoy_0086', 'decoy_0454']
matrix = matrix.reindex(columns=[l for l in LIG_ORDER if l in matrix.columns])
matrix = matrix.reindex(KEY)

fig, ax = plt.subplots(figsize=(9, 5))
im = ax.imshow(matrix.values, cmap='YlOrRd', aspect='auto', vmin=0, vmax=2)

ax.set_xticks(range(len(matrix.columns)))
ax.set_xticklabels(matrix.columns, rotation=30, ha='right', fontsize=11)
ax.set_yticks(range(len(matrix.index)))
ax.set_yticklabels(matrix.index, fontsize=11)

# Annotate cells
for i in range(len(matrix.index)):
    for j in range(len(matrix.columns)):
        val = matrix.values[i, j]
        if val > 0:
            ax.text(j, i, str(val), ha='center', va='center',
                    fontweight='bold', color='white' if val > 1 else 'black')

# Separator between actives and decoys
ax.axvline(1.5, color='black', linewidth=2)

plt.colorbar(im, label='Number of interactions')
plt.title('PLIF: EGFR Allosteric Residue Contacts\n(Actives vs Docking False Positives)',
          fontsize=12, fontweight='bold')
plt.tight_layout()
plt.savefig("data/md/plip/results/plif_heatmap.png", dpi=150)
print("✅ Heatmap saved: data/md/plip/results/plif_heatmap.png")
