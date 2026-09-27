"""
Phase 2 - Step 2.2b
Fix ligand PDBQT files: keep only the FIRST MODEL of each multi-model file.
Vina requires single-model ligand files.
"""
import os
from glob import glob

LIGAND_DIR = 'data/processed/ligands'
fixed = 0
skipped = 0
errors = []

for path in sorted(glob(f'{LIGAND_DIR}/*.pdbqt')):
    try:
        with open(path, 'r') as f:
            lines = f.readlines()

        if not any(l.startswith('MODEL') for l in lines):
            skipped += 1
            continue

        out_lines = []
        in_first = False
        for line in lines:
            if line.startswith('MODEL'):
                if in_first:
                    break
                in_first = True
                continue
            if line.startswith('ENDMDL'):
                break
            out_lines.append(line)

        with open(path, 'w') as f:
            f.writelines(out_lines)
        fixed += 1

    except Exception as e:
        errors.append((path, str(e)))

print(f"Fixed (multi-model -> single): {fixed}")
print(f"Already single-model (skipped): {skipped}")
print(f"Errors: {len(errors)}")
if errors:
    for p, e in errors[:5]:
        print(f"  {p}: {e}")
