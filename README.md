# EGFR Allosteric Inhibitor Discovery

A multi-scale computational pipeline integrating deep learning docking,
molecular dynamics, and quantum mechanical refinement to identify
allosteric inhibitors of the drug-resistant EGFR L858R/T790M/C797S mutant.

## Pipeline
- Phase 1: Target preparation & ligand library
- Phase 2: High-throughput docking + ROC-AUC benchmarking
- Phase 3: 100 ns MD + MM/GBSA
- Phase 4: DFT refinement (HOMO-LUMO, charge transfer)
- Phase 5: Reproducibility & documentation

## Environment
conda env create -f environment_egfr.yml
