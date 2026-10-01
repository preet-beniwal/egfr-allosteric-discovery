"""
Repair missing atoms in the EGFR structure using pdbfixer.
This replaces the manual ChimeraX/PyMOL approach with an automated one.
"""
from pdbfixer import PDBFixer
from openmm.app import PDBFile

INPUT  = 'data/processed/egfr_clean.pdb'
OUTPUT = 'data/processed/egfr_clean_fixed.pdb'

print(f"Loading: {INPUT}")
fixer = PDBFixer(filename=INPUT)

# Detect problems
print("Finding missing residues...")
fixer.findMissingResidues()
print(f"  Missing residues: {dict(fixer.missingResidues)}")

print("Finding missing atoms...")
fixer.findMissingAtoms()
print(f"  Missing atoms: {dict(fixer.missingAtoms)}")
print(f"  Missing terminals: {dict(fixer.missingTerminals)}")

# Repair
print("Adding missing atoms...")
fixer.addMissingAtoms()

print("Adding hydrogens at pH 7.4...")
fixer.addMissingHydrogens(pH=7.4)

# Save
with open(OUTPUT, 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

print(f"\n✅ Repaired structure saved: {OUTPUT}")
