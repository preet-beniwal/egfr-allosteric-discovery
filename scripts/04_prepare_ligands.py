"""
Phase 1 - Step 1.5d
Convert all SMILES (actives + decoys) into 3D PDBQT ligands 
ready for AutoDock Vina docking.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
import os, subprocess

os.makedirs('data/processed/ligands', exist_ok=True)

# Combine actives + decoys
with open('data/raw/known_actives.smi', 'r') as f:
    actives = f.readlines()
with open('data/raw/decoys.smi', 'r') as f:
    decoys = f.readlines()

all_ligands = actives + decoys
success, failed = 0, 0

for line in all_ligands:
    parts = line.strip().split('\t')
    if len(parts) < 2:
        continue
    smi, name = parts[0], parts[1]

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        failed += 1
        continue

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        AllChem.EmbedMultipleConfs(mol, numConfs=5, params=params)
        AllChem.MMFFOptimizeMoleculeConfs(mol)
    except Exception:
        failed += 1
        continue

    pdb_file = f'data/processed/ligands/{name}.pdb'
    pdbqt_file = f'data/processed/ligands/{name}.pdbqt'

    Chem.MolToPDBFile(mol, pdb_file)

    # Convert to PDBQT with OpenBabel (adds hydrogens + assigns Gasteiger charges)
    subprocess.run(
        ['obabel', pdb_file, '-O', pdbqt_file, '-h', '--partialcharge', 'gasteiger'],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )

    if os.path.exists(pdbqt_file) and os.path.getsize(pdbqt_file) > 0:
        success += 1
    else:
        failed += 1

print(f"Successfully prepared: {success} ligands")
print(f"Failed: {failed} ligands")
print(f"Output directory: data/processed/ligands/")
