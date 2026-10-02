"""Phase 3.6 — Plot MM/GBSA binding energies."""
import matplotlib.pyplot as plt
import numpy as np

data = {
    'EAI001':  (-47.60, 'active'),
    'EAI045':  (-46.57, 'active'),
    'decoy_0086': (-42.22, 'decoy'),
    'decoy_0454': (-38.49, 'decoy'),
}

names = list(data.keys())
values = [data[n][0] for n in names]
colors = ['#2e7d32' if data[n][1] == 'active' else '#c62828' for n in names]

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(names, values, color=colors, edgecolor='black', linewidth=1.2)

for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, val - 1.5,
            f'{val:.1f}', ha='center', va='top',
            fontweight='bold', fontsize=11, color='white')

ax.set_ylabel('ΔG_bind (kcal/mol)', fontsize=12)
ax.set_title('MM/GBSA Binding Free Energies', fontsize=13, fontweight='bold')
ax.axhline(0, color='black', linewidth=0.5)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#2e7d32', label='Known actives'),
                   Patch(facecolor='#c62828', label='Docking false positives')]
ax.legend(handles=legend_elements, loc='lower right')

plt.tight_layout()
plt.savefig('data/md/mmpbsa/mmpbsa_ranking.png', dpi=150)
print("Plot saved: data/md/mmpbsa/mmpbsa_ranking.png")
