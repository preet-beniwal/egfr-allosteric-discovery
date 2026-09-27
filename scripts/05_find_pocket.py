"""
Phase 2 - Step 2.1
Calculate the centroid of the EAI045 reference ligand.
This defines the allosteric pocket center for docking.
"""
import numpy as np

coords = []
with open('data/processed/eai045_reference.pdb', 'r') as f:
    for line in f:
        if line.startswith(('ATOM', 'HETATM')):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])

coords = np.array(coords)
center = coords.mean(axis=0)

print(f"Reference ligand atoms: {len(coords)}")
print(f"Pocket center (x, y, z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
print(f"Ligand bounding box size (A): {coords.max(0) - coords.min(0)}")

np.savetxt('data/processed/pocket_center.txt', center, fmt='%.3f')
print("Saved -> data/processed/pocket_center.txt")
