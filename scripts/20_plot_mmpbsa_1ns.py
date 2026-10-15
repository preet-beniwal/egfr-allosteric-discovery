"""Updated MM/GBSA figure with 1 ns (100 frames) results + error bars."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 1 ns results
data = {
    'EAI001':  (-49.55, 0.81, 'active'),
    'EAI045':  (-46.91, 1.63, 'active'),
    'decoy_0454': (-39.21, 1.72, 'decoy'),
    'decoy_0086': (-35.51, 1.94, 'decoy'),
}

names = list(data.keys())
values = [data[n][0] for n in names]
errors = [data[n][1] for n in names]
colors = ['#2e7d32' if data[n][2] == 'active' else '#c62828' for n in names]

fig, ax = plt.subplots(figsize=(9, 5.5))
bars = ax.bar(names, values, yerr=errors, capsize=8,
              color=colors, edgecolor='black', linewidth=1.2,
              error_kw={'elinewidth': 2, 'capthick': 2})

for bar, val, err in zip(bars, values, errors):
    ax.text(bar.get_x() + bar.get_width()/2, val - err - 1.5,
            f'{val:.2f}\n± {err:.2f}',
            ha='center', va='top',
            fontweight='bold', fontsize=10, color='white')

ax.set_ylabel('ΔG_bind (kcal/mol)', fontsize=12)
ax.set_title('MM/GBSA Binding Free Energies (1 ns, 100 frames)',
             fontsize=13, fontweight='bold')
ax.axhline(0, color='black', linewidth=0.5)

# Actives vs decoys separator
ax.axvline(1.5, color='black', linewidth=1, linestyle='--', alpha=0.5)

from matplotlib.patches import Patch
legend_elements = [Patch(facecolor='#2e7d32', label='Known actives'),
                   Patch(facecolor='#c62828', label='Docking false positives')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=10)

plt.tight_layout()
plt.savefig('data/md/mmpbsa_1ns/mmpbsa_ranking_1ns.png', dpi=150)
print("✅ Figure saved: data/md/mmpbsa_1ns/mmpbsa_ranking_1ns.png")
