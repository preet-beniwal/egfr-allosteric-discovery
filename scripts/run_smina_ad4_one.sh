#!/usr/bin/env bash
# Per-ligand smina AD4 docking wrapper for parallel batch execution.
# Args: <ligand_pdbqt> <output_dir> <cx> <cy> <cz>

LIG="$1"; OUTDIR="$2"; CX="$3"; CY="$4"; CZ="$5"
NAME=$(basename "$LIG" .pdbqt)

smina \
  --receptor data/processed/egfr_receptor.pdbqt \
  --ligand "$LIG" \
  --center_x "$CX" --center_y "$CY" --center_z "$CZ" \
  --size_x 25 --size_y 25 --size_z 25 \
  --scoring ad4_scoring \
  --exhaustiveness 8 \
  --num_modes 1 \
  --out "$OUTDIR/${NAME}_out.pdbqt" \
  > "$OUTDIR/${NAME}.log" 2>&1

AFF=$(grep -m1 "minimizedAffinity" "$OUTDIR/${NAME}_out.pdbqt" 2>/dev/null | awk '{print $3}')
[ -z "$AFF" ] && AFF="NaN"
echo -e "${NAME}\t${AFF}"
