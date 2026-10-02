"""
Phase 3.7 — Protein-Ligand Interaction Fingerprints (PLIF)
Fixed: use plf.Molecule.from_mda() instead of MDAnalysis direct conversion.
"""
import MDAnalysis as mda
import prolif as plf
import pandas as pd
import numpy as np
import os

LIGANDS = ['EAI045', 'EAI001', 'decoy_0086', 'decoy_0454']

INTERACTION_TYPES = [
    'HBAcceptor', 'HBDonor', 'Hydrophobic', 'PiStacking',
    'CationPi', 'PiCation', 'VdWContact', 'MetalAcceptor'
]

all_results = {}

for LIG in LIGANDS:
    print(f"\n{'='*60}")
    print(f"  Analyzing: {LIG}")
    print(f"{'='*60}")
    
    plif_dir = f"data/md/plif/{LIG}"
    
    # Load dry trajectory
    u = mda.Universe(f"{plif_dir}/npt_dry.pdb", f"{plif_dir}/npt_dry.xtc")
    print(f"  Atoms: {len(u.atoms)}")
    print(f"  Frames: {len(u.trajectory)}")
    
    # Select protein and ligand
    protein = u.select_atoms("protein")
    ligand  = u.select_atoms("resname UNL")
    
    print(f"  Protein atoms: {len(protein)}")
    print(f"  Ligand atoms: {len(ligand)}")
    
    if len(ligand) == 0:
        print(f"  ❌ No ligand found")
        continue
    
    # Use ProLIF's own converter (handles H detection correctly)
    try:
        lig_mol = plf.Molecule.from_mda(ligand)
        prot_mol = plf.Molecule.from_mda(protein)
        print(f"  ✅ Converted to RDKit: lig={lig_mol.n_atoms}, prot={prot_mol.n_atoms}")
    except Exception as e:
        print(f"  ❌ ProLIF conversion failed: {e}")
        continue
    
    # Run ProLIF fingerprint
    fp = plf.Fingerprint(INTERACTION_TYPES)
    
    try:
        fp.run(u.trajectory, lig_mol, prot_mol, progress=False)
    except Exception as e:
        print(f"  ❌ ProLIF run failed: {e}")
        continue
    
    # Aggregate results
    df = fp.to_dataframe()
    
    if df.empty or df.shape[1] == 0:
        print(f"  ⚠️  No interactions detected")
        all_results[LIG] = pd.Series(dtype=float)
        continue
    
    # Occupancy = fraction of frames each interaction appears in
    occupancy = df.sum(axis=0) / len(df)
    occupancy = occupancy[occupancy > 0.05]  # >5% occupancy
    occupancy = occupancy.sort_values(ascending=False)
    
    all_results[LIG] = occupancy
    
    print(f"  Interactions detected: {len(occupancy)}")
    print(f"  Top interactions:")
    for (residue, interaction), occ in occupancy.head(10).items():
        try:
            res_str = f"{residue[1]}{residue[2]}"
        except:
            res_str = str(residue)
        print(f"    {res_str:20s} {interaction:20s} {occ:.1%}")

# --- Save combined results ---
os.makedirs("data/md/plif/results", exist_ok=True)

# Build a comparison dataframe
rows = []
for LIG, occ in all_results.items():
    for (residue, interaction), frac in occ.items():
        try:
            res_str = f"{residue[1]}{residue[2]}"
        except:
            res_str = str(residue)
        rows.append({
            'Ligand': LIG,
            'Residue': res_str,
            'Interaction': interaction,
            'Occupancy': frac
        })

results_df = pd.DataFrame(rows)
results_df.to_csv("data/md/plif/results/plif_full.csv", index=False)

# Pivot: interactions × ligands
if not results_df.empty:
    pivot = results_df.pivot_table(
        index=['Residue', 'Interaction'],
        columns='Ligand',
        values='Occupancy',
        fill_value=0
    )
    pivot.to_csv("data/md/plif/results/plif_matrix.csv")
    
    print(f"\n{'='*60}")
    print(f"  PLIF ANALYSIS COMPLETE")
    print(f"{'='*60}")
    print(f"  Total interactions tracked: {len(pivot)}")
    print(f"  Results: data/md/plif/results/plif_matrix.csv")
else:
    print("\n⚠️  No interactions found across any system")
