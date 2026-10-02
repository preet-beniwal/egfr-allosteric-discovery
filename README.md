# EGFR Allosteric Inhibitor Discovery

**A docking and interaction-fingerprinting pipeline for allosteric inhibitor discovery against drug-resistant EGFR T790M, with preliminary MM/GBSA rescoring.**

---

## Motivation

Third-generation EGFR inhibitors like osimertinib work against the T790M gatekeeper mutation until the C797S mutation emerges. At that point every approved ATP-competitive inhibitor fails, because the resistance mechanism directly disrupts the ATP-site binding geometry. Allosteric inhibition — targeting the pocket between the αC-helix and the ATP cleft — offers a route around this problem because allosteric binding does not depend on the ATP site.

The computational difficulty is that allosteric pockets are more plastic than ATP sites, and docking scoring functions are worse at discriminating true binders from geometrically compatible decoys in flexible pockets. This project tests whether a docking pipeline can be validated against known biology, and whether interaction fingerprinting can recover the experimentally observed binding mode.

The starting structure is the crystal structure of T790M/V948R EGFR bound to the allosteric inhibitor EAI045 (PDB 5D41). The compound library is 16 literature-validated allosteric EGFR inhibitors plus 500 drug-like decoys generated with RDKit. The pipeline docks all 516 compounds against the allosteric pocket, benchmarks the protocol against the actives/decoys split, and analyses the interaction geometry of the top hits.

---

## Key Findings

### 1. Docking validation: ROC-AUC 0.994

AutoDock Vina separated 16 literature-validated allosteric EGFR inhibitors from 500 RDKit decoys with ROC-AUC 0.994.

For allosteric site docking this is very high. Published benchmarks for allosteric pockets typically land between 0.60 and 0.75, because the shallower, more flexible geometry makes scoring harder. A value of 0.994 means the pocket definition, box placement, and scoring parameters are correct before any MD is attempted. This is the strongest result in the project because it does not depend on anything downstream.

### 2. Interaction fingerprinting recovered the experimental binding mode

The top 4 docking hits (2 known actives, 2 false positives) were simulated in MD and analysed with PLIP. Residue numbers were remapped to EGFR canonical numbering.

Both known actives engage the same three-residue interaction signature:

| Residue | Interaction type | Notes |
|---------|-----------------|-------|
| LYS745 | H-bond + hydrophobic | Catalytic lysine, αC-helix |
| PHE856 | H-bond + hydrophobic | DFG motif phenylalanine |
| LEU858 | Hydrophobic | Allosteric pocket |

Both docking false positives bind adjacent but distinct residues and do not engage this signature.

This LYS745/PHE856/LEU858 pattern is the **experimentally observed binding mode of EAI045** in the 5D41 co-crystal structure (Jia et al., *Nature* 2016). Recovering it from docking poses alone — without any prior knowledge of the interaction geometry — is a genuine validation that the docking protocol produces biologically meaningful poses, not just numerically favourable ones.

### 3. Preliminary MM/GBSA rescoring demoted both false positives

This result is preliminary and is reported with that caveat. See the Limitations section for the reason.

MM/GBSA was computed on a 100 ps equilibrated ensemble for each of the 4 simulated systems:

| Rank | Ligand | MM/GBSA ΔG_bind | Category |
|------|--------|-----------------|----------|
| 1 | EAI001 | −47.60 kcal/mol | Known active |
| 2 | EAI045 | −46.57 kcal/mol | Known active (crystal ligand) |
| 3 | decoy_0086 | −42.22 kcal/mol | Docking false positive |
| 4 | decoy_0454 | −38.49 kcal/mol | Docking false positive |

Mean active: −47.1 kcal/mol. Mean decoy: −40.4 kcal/mol. Gap: 6.7 kcal/mol.

Both known actives ranked above both docking false positives. The direction of the separation is consistent with the docking result and the interaction fingerprint. But the sampling is far too short for MM/GBSA convergence (see Limitations), so this should be treated as a directional check rather than a validated binding energy calculation.

---

## Pipeline Overview

Four phases.

| Phase | What | Tools |
|-------|------|-------|
| **1** | Target and ligand preparation | PyMOL, PDBFixer, RDKit |
| **2** | High-throughput docking | AutoDock Vina 1.2.0 |
| **3** | Short MD + MM/GBSA rescoring | GROMACS 2026.3, gmx_MMPBSA |
| **4** | Interaction fingerprinting | PLIP + custom renumbering |

**Phase 1** cleans PDB 5D41, removes waters and the co-crystallised ligand, generates 3D conformers for 16 literature actives and 500 RDKit decoys (516 total).

**Phase 2** docks every ligand against the allosteric pocket. The box is centred on the EAI045 binding site. Scoring is AutoDock Vina default. The protocol is benchmarked with ROC-AUC and enrichment factor at 1%.

**Phase 3** takes the top 4 docking hits through minimisation, 100 ps NVT, and 100 ps NPT in GROMACS with AMBER14SB and TIP3P. Ligand parameters are generated with ACPYPE (GAFF2). MM/GBSA is computed on the equilibrated ensemble with gmx_MMPBSA.

**Phase 4** runs PLIP on MD-snapshotted complexes and remaps residue numbering from the PDB structure numbering to EGFR canonical numbering, so that the interaction fingerprint is directly comparable to the published crystal structure.

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
    |   +-- md/
    |       |-- mdp/           # GROMACS parameter files
    |       |-- topologies/    # GAFF2 ligand topologies (ACPYPE)
    |       |-- systems/       # One directory per simulated ligand
    |       |-- mmpbsa/        # MM/GBSA results
    |       +-- plip/          # Interaction fingerprints + figures
    +-- logs/

---

## Reproducibility

Two conda environments. They are kept separate because GROMACS and gmx_MMPBSA pin different CUDA and Python dependencies.

    # Main analysis environment
    conda env create -f environment_egfr.yml
    conda activate egfr

    # MD environment
    conda create -n md -c conda-forge python=3.11 gromacs ambertools gmx_MMPBSA
    conda activate md

Full pipeline:

    python scripts/01_clean_protein.py
    python scripts/02_build_actives.py
    python scripts/03_generate_decoys.py
    python scripts/04_prepare_ligands.py
    python scripts/05_find_pocket.py
    python scripts/06_dock_all.py             # ~30 min on 4 CPU threads
    python scripts/08_roc_analysis.py         # produces ROC-AUC + EF1%

    # MD per ligand (see data/md/systems/<ligand>/ for exact gmx commands)
    # MM/GBSA after MD (gmx_MMPBSA)
    python scripts/13_plot_mmpbsa.py

    # PLIP interaction fingerprints
    plip -f data/md/plip/structures/EAI045_dry.pdb \
         -o data/md/plip/results/EAI045 -x
    python scripts/15_parse_plip.py
    python scripts/16_plip_mapped.py
    python scripts/17_plot_plif.py

Docking takes about 30 minutes. MD equilibration for the four selected systems took roughly 12 hours. MM/GBSA on 100 ps ensembles takes 20–40 minutes per system.

---

## Software

- Structure preparation: PyMOL, PDBFixer, RDKit
- Docking: AutoDock Vina 1.2.0
- MD engine: GROMACS 2026.3 (AMBER14SB, TIP3P)
- MM/GBSA: gmx_MMPBSA
- Ligand parameterisation: ACPYPE (GAFF2)
- Interaction analysis: PLIP 2.3
- Analysis and plotting: Python 3.11, MDAnalysis, NumPy, pandas, matplotlib

---

## Limitations

This is the section that matters most for interpreting the results, because the pipeline has one real weakness.

**The MD sampling is far too short for MM/GBSA convergence.** Each of the four simulated complexes was run for 100 ps NVT + 100 ps NPT. Published MM/GBSA protocols typically use 10–50 ns of production MD per complex, with block averaging across independent windows to demonstrate convergence. At 100 ps, the ensemble is equilibrated enough to relax side chains and water around the docking pose, but not enough to sample slow motions such as αC-helix breathing or DFG-motif flipping. The estimated error on a 100 ps MM/GBSA calculation is easily 3–8 kcal/mol — comparable to the 6.7 kcal/mol separation reported here. The MM/GBSA result should be read as a directional filter, not a validated binding energy.

**Other limitations:**

- The docking box was defined around the EAI045 pose from the crystal structure. Ligands that engage a different sub-pocket of the allosteric site would score poorly even if they are genuine binders. This is a common bias in allosteric docking and is not corrected here.
- The decoy set is RDKit-generated. It is chemically reasonable but does not sample the distribution of real experimental screening actives. A benchmark against confirmed experimental negatives would be stronger.
- Only one crystal structure was used. Allosteric pockets are plastic, and a different starting conformation could change the ranking.
- No experimental validation. Everything here is computational. The MM/GBSA separation is a model result, not a measured one.
- Absolute MM/GBSA values carry systematic error from the implicit solvent model and the neglect of explicit entropy. Only the relative ranking is meaningful, and only if the sampling is sufficient — which for this project it is not.

---

## Future Work

Three extensions would strengthen this project into a publishable study:

1. **Extended MD.** Rerun the four systems for 10 ns each with block-averaged MM/GBSA across three independent windows. This is the minimum needed to make the MM/GBSA result defensible. Estimated compute: 40–80 CPU-hours.
2. **Free energy perturbation.** For the two known actives, run FEP or TI against a thermodynamic cycle to establish a converged reference. If FEP confirms the MM/GBSA ranking, the short-MD result gains credibility.
3. **Experimental validation.** The pipeline output is a hypothesis. Biochemical assays (SPR or ITC against the isolated EGFR kinase domain) or cell-based assays (TKI-resistance reversal in C797S-mutant lines) would be required to claim any of these compounds as real allosteric binders.

The docking protocol and interaction fingerprint are already solid. The extension work is what turns the preliminary MM/GBSA hint into a validated result.

---

## Data Availability

All code is in this repository. The receptor structure (5D41) is public from RCSB. Ligand SMILES are in `scripts/02_build_actives.py` — actives come from published literature, decoys are generated deterministically with a fixed random seed.

MD trajectories (.xtc, ~2 GB each) are not committed but are regenerated from the pipeline. Intermediate GROMACS files (.edr, .log, .tpr) are also excluded. All final outputs are in `data/md/mmpbsa/` and `data/md/plip/`.

---

## Citation

    @software{beniwal2026egfr,
      title  = {Docking and Interaction Fingerprinting Pipeline for
                Allosteric Inhibitor Discovery Against EGFR T790M},
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
