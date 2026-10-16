#!/usr/bin/env python3
"""
Parse block-averaged MM/GBSA results from FINAL_RESULTS.dat.

The DAT file contains a "Delta (Complex - Receptor - Ligand):" section
whose ΔTOTAL row holds the block-average binding free energy in the
first numeric column. This is what we extract.
"""

import os
import pandas as pd
import numpy as np

BLOCKS = "data/md/mmpbsa_blocks"
LIGANDS = ["EAI045", "EAI001", "decoy_0086", "decoy_0454"]


def parse_dat(path):
    """Return the Delta TOTAL (kcal/mol) from a gmx_MMPBSA DAT file."""
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    in_delta = False
    for line in lines:
        if line.startswith("Delta (Complex"):
            in_delta = True
            continue
        if in_delta and "TOTAL" in line:
            parts = line.split()
            # parts[0] is ΔTOTAL (may be garbled if encoding fails);
            # parts[1] is the average binding energy
            try:
                return float(parts[1])
            except (IndexError, ValueError):
                continue
    return None


rows = []
for lig in LIGANDS:
    block_dg = []
    for b in range(4):
        dat = f"{BLOCKS}/{lig}_block{b}/FINAL_RESULTS.dat"
        if not os.path.exists(dat):
            print(f"  missing {dat}")
            block_dg.append(None)
            continue
        dg = parse_dat(dat)
        block_dg.append(round(dg, 2) if dg is not None else None)

    valid = [v for v in block_dg if v is not None]
    if len(valid) == 4:
        rows.append({
            "ligand": lig,
            "block0": block_dg[0], "block1": block_dg[1],
            "block2": block_dg[2], "block3": block_dg[3],
            "mean": round(np.mean(valid), 2),
            "std":  round(np.std(valid, ddof=1), 2),
            "range": round(max(valid) - min(valid), 2),
            "category": "Active" if lig in ("EAI045", "EAI001") else "Decoy",
        })

df = pd.DataFrame(rows)
out = f"{BLOCKS}/block_averaging_summary.tsv"
df.to_csv(out, sep="\t", index=False)
print(df.to_string(index=False))
print(f"\nwrote {out}")

