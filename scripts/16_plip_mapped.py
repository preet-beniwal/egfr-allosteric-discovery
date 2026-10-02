"""
Phase 3.7b — Map PLIP residue numbers back to EGFR canonical numbering.
"""
import pandas as pd

# Offset determined from pdb2gmx renumbering
# Original first residue in 5D41 chain A = 698
OFFSET = 697

# Load PLIP matrix
df = pd.read_csv("data/md/plip/results/plip_interactions.csv")

# Extract residue number from "LYS48" style strings
import re
def map_residue(res_str):
    m = re.match(r'([A-Z]+)(\d+)', res_str)
    if not m:
        return res_str, 0
    name, num = m.group(1), int(m.group(2))
    new_num = num + OFFSET
    return f"{name}{new_num}", new_num

df[['Residue_EGFR', 'EGFR_Number']] = df['Residue'].apply(
    lambda x: pd.Series(map_residue(x))
)

# Save mapped version
df.to_csv("data/md/plip/results/plip_interactions_EGFR.csv", index=False)

print("=== Residue mapping ===")
print("  Original → EGFR canonical")
for _, row in df.drop_duplicates('Residue').iterrows():
    print(f"  {row['Residue']:10s} → {row['Residue_EGFR']:10s}")

# Highlight key EGFR residues
KEY_RESIDUES = ['LYS745', 'LEU718', 'MET790', 'MET793', 'ASP855', 'PHE856',
                'PHE723', 'THR854', 'LEU858', 'VAL726']

print(f"\n{'='*60}")
print(f"  KEY EGFR RESIDUE INTERACTIONS")
print(f"{'='*60}")

key_df = df[df['Residue_EGFR'].isin(KEY_RESIDUES)]

# Pivot: ligand × (residue, interaction)
if not key_df.empty:
    matrix = key_df.groupby(['Residue_EGFR', 'Interaction', 'Ligand']).size().unstack(fill_value=0)
    print(matrix.to_string())
    matrix.to_csv("data/md/plip/results/plip_key_residues.csv")
else:
    print("  No interactions with key residues found")

# Summary per ligand
print(f"\n{'='*60}")
print(f"  PER-LIGAND INTERACTION SUMMARY")
print(f"{'='*60}")
for lig in df['Ligand'].unique():
    sub = df[df['Ligand'] == lig]
    print(f"\n  {lig}:")
    print(f"    Total interactions: {len(sub)}")
    print(f"    H-bonds:            {len(sub[sub['Interaction']=='HBond'])}")
    print(f"    Hydrophobic:        {len(sub[sub['Interaction']=='Hydrophobic'])}")
    print(f"    Unique residues:    {sub['Residue_EGFR'].nunique()}")
