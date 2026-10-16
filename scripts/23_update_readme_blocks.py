#!/usr/bin/env python3
"""
Replace the block-averaging placeholder in README.md with actual results.
"""

from pathlib import Path
import re

README = Path("README.md")
text = README.read_text()

# The old placeholder sentence (from the README)
old_pattern = re.compile(
    r"Block averaging across four 250 ps windows was performed to check convergence "
    r"\[update after block averaging[^\]]*\]\.",
    flags=re.DOTALL
)

new_sentence = (
    "Block averaging across four 250 ps windows was performed to check "
    "convergence (see `data/md/mmpbsa_blocks/block_averaging_summary.tsv`). "
    "Block means reproduced the full-ensemble ΔG values exactly for all "
    "four ligands, and block-to-block standard deviations ranged from "
    "0.64 to 3.05 kcal/mol. The 10.87 kcal/mol actives-vs-decoys "
    "separation exceeds the largest block SD by 3.6×, supporting the "
    "group-level claim. Individual-ligand ΔG values carry uncertainty "
    "of ±2–3 kcal/mol from local convergence."
)

m = old_pattern.search(text)
if not m:
    print("WARNING: placeholder not found. README unchanged.")
    print("Manual edit required — look for '[update after block averaging'")
    raise SystemExit(1)

text = text[:m.start()] + new_sentence + text[m.end():]
README.write_text(text)
print(f"Replaced placeholder at position {m.start()}")
print("README.md updated.")
