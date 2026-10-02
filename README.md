# EGFR Allosteric Inhibitor Discovery

**A multi-scale computational pipeline for allosteric inhibitor discovery against drug-resistant EGFR T790M.**

---

## Motivation

Third-generation EGFR inhibitors like osimertinib work well against T790M-mutant lung cancer until the C797S mutation emerges. At that point, every ATP-competitive inhibitor fails. Allosteric inhibition — targeting the pocket between the αC-helix and the ATP site — offers a route around this problem because the allosteric binding mode does not depend on the ATP cleft geometry that C797S disrupts.

But virtual screening alone cannot solve this. Three specific problems:

1. Allosteric pockets change shape when ligands bind. Rigid docking misses this.
2. Empirical scoring functions cannot reliably distinguish actives from geometrically compatible decoys.
3. Classical force fields ignore the electronic effects — charge transfer, polarization — that determine real binding affinity.

This project builds a pipeline that addresses all three. It starts from the crystal structure of T790M EGFR bound to the allosteric inhibitor EAI045 (PDB 5D41) and carries a curated compound library through docking, MD, MM/GBSA, and interaction fingerprinting. The goal is not to invent a new inhibitor. The goal is to demonstrate that a well-designed multi-scale pipeline can separate true binders from docking artifacts — a result that would generalise to other allosteric targets.

---

## Key Findings

Three results, in the order they matter for the project.

### 1. Docking validation: ROC-AUC 0.994

The AutoDock Vina protocol separated 16 literature-validated allosteric EGFR inhibitors from 500 drug-like decoys with ROC-AUC = 0.994. For allosteric site docking, where published benchmarks typically land between 0.60 and 0.75, this is very high. It means the pocket definition, box placement, and scoring setup are correct before we spend compute on MD.

### 2. MM/GBSA demoted both docking false positives

Both known actives and two false positives sat in the top 5 docking hits. MM/GBSA on equilibrated MD ensembles separated them cleanly:

| Rank | Ligand | MM/GBSA ΔG_bind | Category |
|------|--------|-----------------|----------|
| 1 | EAI001 | −47.60 kcal/mol | Active |
| 2 | EAI045 | −46.57 kcal/mol | Active (crystal ligand) |
| 3 | decoy_0086 | −42.22 kcal/mol | Docking false positive |
| 4 | decoy_0454 | −38.49 kcal/mol | Docking false positive |

Actives average −47.1 kcal/mol. Decoys average −40.4 kcal/mol. The 6.7 kcal/mol gap is the whole point — Vina could not separate them, but MM/GBSA on the MD ensemble could.

### 3. Interaction fingerprinting identified the allosteric signature

PLIP interaction analysis, remapped to EGFR canonical numbering, showed that both actives engage the same three-residue signature and both false positives do not:

- **LYS745** (catalytic lysine) — hydrogen bond + hydrophobic contact
- **PHE856** (DFG motif phenylalanine) — hydrogen bond + hydrophobic contact
- **LEU858** (allosteric pocket) — hydrophobic contact

This matches the experimental binding mode of EAI045 in the 5D41 crystal structure (Jia et al., *Nature* 2016). The pipeline rediscovered the known binding mode independently, from docking geometry alone. That is the validation that matters.

---

## Pipeline Overview

Four phases, in order.

| Phase | What | Key tools |
|-------|------|-----------|
| **1** | Target preparation | PyMOL, PDBFixer, RDKit |
| **2** | High-throughput docking | AutoDock Vina 1.2.0 |
| **3** | MD + MM/GBSA | GROMACS 2026.3, gmx_MMPBSA |
| **4** | Interaction analysis | PLIP, custom renumbering |

Phase 1 cleans PDB 5D41, generates the ligand library (16 actives from literature, 500 decoys from RDKit), and produces 3D conformers for all 516 compounds.

Phase 2 docks every ligand against the allosteric pocket and computes ROC-AUC and enrichment factor at 1%. The result is a ranked list of hits and a validated protocol.

Phase 3 takes the top 4 hits through 100 ps NVT and 100 ps NPT equilibration in GROMACS with AMBER14SB and TIP3P. MM/GBSA binding energies are computed on the equilibrated ensemble with gmx_MMPBSA.

Phase 4 runs PLIP on the MD-snapshotted protein–ligand complexes, extracts hydrogen bonds, hydrophobic contacts, and π-stacking, then remaps residue numbering from the PDB structure to EGFR canonical numbering.

---

## Repository Structure

    egfr-allosteric-discovery/
    |-- README.md
    |-- LICENSE
    |-- environment_egfr.yml
    |-- scripts/
    |   |-- 01_clean_protein.py
    |   |-- 02_build_actives.py
    |   |-- 03_generate_decoys.py
    |   |-- 04_prepare_ligands.py
    |   |-- 05_find_pocket.py
    |   |-- 06_dock_all.py
    |   |-- 07_fix_multimodel.py
    |   |-- 08_roc_analysis.py
    |   |-- 13_plot_mmpbsa.py
    |   |-- 14_plif_analysis.py
    |   |-- 15_parse_plip.py
    |   |-- 16_plip_mapped.py
    |   +-- 17_plot_plif.py
    |-- data/
    |   |-- raw/               # PDB 5D41, curated SMILES
    |   |-- processed/         # Cleaned receptor, prepared ligands
    |   |-- results/           # Docking results + ROC plot
    |   +-- md/                # MD workflow
    |       |-- mdp/           # GROMACS parameter files
    |       |-- topologies/    # Ligand GAFF2 topologies
    |       |-- systems/       # Per-ligand MD directories
    |       |-- mmpbsa/        # MM/GBSA results
    |       +-- plip/          # Interaction fingerprints + figures
    +-- logs/

---

## Reproducibility

Two conda environments, both declared as YAML files.

    # Main analysis environment
    conda env create -f environment_egfr.yml
    conda activate egfr

    # MD environment (separate — GROMACS and gmx_MMPBSA pin different CUDA/deps)
    conda create -n md -c conda-forge python=3.11 gromacs ambertools gmx_MMPBSA
    conda activate md

    # Full pipeline
    python scripts/01_clean_protein.py
    python scripts/02_build_actives.py
    python scripts/03_generate_decoys.py
    python scripts/04_prepare_ligands.py
    python scripts/05_find_pocket.py
    python scripts/06_dock_all.py             # ~30 min, CPU-parallel
    python scripts/08_roc_analysis.py         # produces AUC + EF1%

    # MD + MM/GBSA (per ligand)
    # run from data/md/systems/<ligand>/
    gmx grompp -f ../../mdp/minim.mdp -c complex.gro ...
    # ... see per-system README in data/md/systems/ for full commands

    # PLIP interaction fingerprints
    plip -f data/md/plip/structures/EAI045_dry.pdb -o data/md/plip/results/EAI045 -x
    python scripts/15_parse_plip.py
    python scripts/16_plip_mapped.py
    python scripts/17_plot_plif.py

The docking script runs in about 30 minutes on 4 CPU threads. MD equilibration for the four selected systems took roughly 12 hours overnight. MM/GBSA on 100-ps windows takes 20–40 min per system.

---

## Software

- **Structure preparation**: PyMOL, PDBFixer, RDKit
- **Docking**: AutoDock Vina 1.2.0
- **MD engine**: GROMACS 2026.3 (AMBER14SB force field, TIP3P water)
- **MM/GBSA**: gmx_MMPBSA
- **Ligand parameterization**: ACPYPE (GAFF2)
- **Interaction analysis**: PLIP 2.3
- **Analysis + plotting**: Python 3.11, MDAnalysis, NumPy, pandas, matplotlib

---

## Limitations

Several things this pipeline does **not** do, and they matter for interpreting the results.

- **The MD windows are short.** 100 ps NVT + 100 ps NPT is enough to relax the complex from a docking pose, but not enough to sample ligand binding–unbinding events or slow loop motions. The MM/GBSA energies are converged to the extent that the equilibrated ensemble is representative — for well-ordered allosteric pockets this holds, but it is a real assumption.

- **MM/GBSA ignores explicit solvent entropy.** The separation between actives and decoys is clean (6.7 kcal/mol), but absolute ΔG values carry systematic error. The ranking is what matters, not the numbers.

- **Docking box was defined around EAI045's pose** from the crystal structure. This biases the protocol toward EAI045-like binding modes. A ligand that engages a different sub-pocket of the allosteric site would score poorly despite being a genuine binder.

- **The decoy set is RDKit-generated.** This gives chemically reasonable molecules but does not sample the distribution of real high-throughput screening actives. A decoy set drawn from true experimental negatives would be a stronger benchmark.

- **No experimental validation.** Everything is computational. The MM/GBSA separation between actives and decoys is a model result, not a measured one. Biochemical assays (SPR, ITC, or cell-based TKI-resistance reversal) would be required to claim any of these compounds as real allosteric binders.

- **Single crystal structure.** PDB 5D41 is one conformation of one mutant (T790M/V948R). Allosteric pockets are known to be structurally plastic, and a different starting structure could change the results.

---

## Data Availability

All code is in this repository. The receptor structure (5D41) is public from RCSB. Ligand SMILES are curated in `scripts/02_build_actives.py` — the actives come from published literature, the decoys are generated deterministically from a fixed random seed.

MD trajectories (.xtc files, ~2 GB each) are not committed — they are regenerated from the pipeline. Intermediate GROMACS files (.edr, .log, .tpr) are also excluded. All final analysis outputs are in `data/md/mmpbsa/` and `data/md/plip/`.

---

## Citation

    @software{beniwal2026egfr,
      title  = {Multi-Scale Computational Pipeline for Allosteric
                Inhibitor Discovery Against EGFR T790M},
      author = {Beniwal, Preet},
      year   = {2026},
      url    = {https://github.com/preet-beniwal/egfr-allosteric-discovery},
      note   = {Preprint in preparation}
    }

Key methods references:

- Eberhardt et al. (2021) AutoDock Vina 1.2.0 — *J. Chem. Inf. Model.*
- Abraham et al. (2015) GROMACS — *SoftwareX*
- Valdés-Tresanco et al. (2021) gmx_MMPBSA — *J. Chem. Theory Comput.*
- Adasme et al. (2021) PLIP 2021 — *Nucleic Acids Res.*
- Jia et al. (2016) EAI045 crystal structure — *Nature*

---

## Author

**Preet Beniwal** — Independent bioinformatics portfolio project.
GitHub: [@preet-beniwal](https://github.com/preet-beniwal)

---

## License

MIT — see [LICENSE](LICENSE).

---

## Disclaimer

This is computational work. Every compound identified here is a hypothesis, not a drug. Experimental validation — biochemical assays, cell-based TKI-resistance reversal, and eventually animal models — would be required before any therapeutic claim could be made. The pipeline is provided as a research methodology template.
