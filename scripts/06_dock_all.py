"""
Phase 2 - Step 2.4
High-throughput docking of all 516 ligands against EGFR T790M
using AutoDock Vina Python API.

Run with:  nohup python scripts/06_dock_all.py > logs/docking.log 2>&1 &
"""
import os
import csv
import time
import numpy as np
from vina import Vina
from glob import glob

# ---------- CONFIGURATION ----------
RECEPTOR     = 'data/processed/egfr_receptor.pdbqt'
LIGAND_DIR   = 'data/processed/ligands'
OUTPUT_CSV   = 'data/results/docking_scores.csv'
OUTPUT_POSES = 'data/results/poses'

center = np.loadtxt('data/processed/pocket_center.txt')
BOX_SIZE        = [25, 25, 25]
EXHAUSTIVENESS  = 4
N_POSES         = 3

# ---------- SETUP ----------
os.makedirs(OUTPUT_POSES, exist_ok=True)

ligand_files = sorted(glob(f'{LIGAND_DIR}/*.pdbqt'))
print(f"Found {len(ligand_files)} ligand files", flush=True)
print(f"Receptor: {RECEPTOR}", flush=True)
print(f"Pocket center: {center.tolist()}", flush=True)
print(f"Exhaustiveness: {EXHAUSTIVENESS}", flush=True)
print("=" * 60, flush=True)

# Initialize Vina ONCE — this is a huge speedup
v = Vina(sf_name='vina', cpu=0, verbosity=0)
v.set_receptor(RECEPTOR)
v.compute_vina_maps(center=center.tolist(), box_size=BOX_SIZE)
print("Vina maps computed. Starting docking loop...\n", flush=True)

# ---------- DOCKING LOOP ----------
results = []
t_start = time.time()

for i, lig_path in enumerate(ligand_files, 1):
    name = os.path.basename(lig_path).replace('.pdbqt', '')
    try:
        v.set_ligand_from_file(lig_path)
        v.dock(exhaustiveness=EXHAUSTIVENESS, n_poses=N_POSES)
        energies = v.energies(n_poses=N_POSES)
        best_score = float(energies[0][0])

        v.write_poses(f'{OUTPUT_POSES}/{name}_out.pdbqt',
                      n_poses=N_POSES, overwrite=True)
        results.append({'ligand': name,
                        'affinity_kcal_mol': round(best_score, 3)})
        print(f"[{i:3d}/{len(ligand_files)}] {name:25s} -> {best_score:7.2f} kcal/mol", flush=True)

    except Exception as e:
        print(f"[{i:3d}/{len(ligand_files)}] {name:25s} -> ERROR: {e}", flush=True)
        results.append({'ligand': name, 'affinity_kcal_mol': None})

# ---------- SAVE ----------
with open(OUTPUT_CSV, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['ligand', 'affinity_kcal_mol'])
    writer.writeheader()
    writer.writerows(results)

elapsed = (time.time() - t_start) / 60
print(f"\n{'=' * 60}", flush=True)
print(f"✅ Docking complete in {elapsed:.1f} minutes", flush=True)
print(f"Results saved: {OUTPUT_CSV}", flush=True)
print(f"Poses saved:   {OUTPUT_POSES}/", flush=True)
