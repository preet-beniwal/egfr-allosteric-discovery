"""
Phase 1 - Step 1.5c
Generate 500 drug-like decoy molecules using RDKit fragment assembly.
These serve as inactives for ROC-AUC benchmarking.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
import random

random.seed(42)  # reproducibility

# Common drug-like fragments and linkers
fragments = [
    'c1ccccc1', 'c1ccncc1', 'C1CCNCC1', 'C1CCOCC1',
    'c1ccsc1', 'c1cc[nH]c1', 'C(=O)N', 'C(=O)O',
    'S(=O)(=O)N', 'F', 'Cl', 'C#N', 'OC',
    'C(F)(F)F', 'C', 'O', 'N'
]

decoys = set()
while len(decoys) < 500:
    n_frags = random.randint(2, 4)
    smi = ''.join(random.choices(fragments, k=n_frags))
    if random.random() > 0.5:
        smi += '(' + random.choice(fragments) + ')'
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None and 200 < Descriptors.MolWt(mol) < 500:
            decoys.add(Chem.MolToSmiles(mol))
    except Exception:
        pass

with open('data/raw/decoys.smi', 'w') as f:
    for i, smi in enumerate(decoys):
        f.write(f"{smi}\tdecoy_{i:04d}\n")

print(f"Generated {len(decoys)} decoys -> data/raw/decoys.smi")
