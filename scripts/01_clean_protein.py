"""
Phase 1 - Step 1.4
Clean EGFR T790M structure (PDB: 5D41) for allosteric docking.

KEEPS:
  - Protein (EGFR kinase domain)
  - 57N (EAI045 - allosteric inhibitor, reference ligand)

REMOVES:
  - HOH (water)
  - ANP (ATP analog)
  - MG  (magnesium ion)
"""

from pymol import cmd

INPUT       = "data/raw/5D41.pdb"
OUT_PROTEIN = "data/processed/egfr_clean.pdb"
OUT_LIGAND  = "data/processed/eai045_reference.pdb"

# Load raw structure
cmd.load(INPUT, "egfr")
print(f"Loaded {INPUT}: {cmd.count_atoms('egfr')} atoms")

# ---- 1. Extract EAI045 (57N) BEFORE removing it ----
cmd.extract("lig", "resn 57N")
cmd.save(OUT_LIGAND, "lig")
print(f"Saved reference ligand -> {OUT_LIGAND} ({cmd.count_atoms('lig')} atoms)")

# ---- 2. Remove everything that is NOT protein ----
cmd.remove("not polymer")
print(f"After removing all HETATM: {cmd.count_atoms('egfr')} atoms")

# ---- 3. Remove alternate conformations (keep altloc A only) ----
cmd.remove("not alt ''+A")

# ---- 4. Add hydrogens (approximate; PDB2PQR refinement comes later) ----
cmd.h_add("all")

# ---- 5. Save cleaned protein ----
cmd.save(OUT_PROTEIN, "egfr")
print(f"Saved cleaned protein -> {OUT_PROTEIN} ({cmd.count_atoms('egfr')} atoms)")

cmd.reinitialize()
print("DONE.")
