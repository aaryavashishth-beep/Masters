# TDP-43 RRM1/RRM2 Docking Project

## Research work till date
TDP-43 is an RNA-binding protein whose aggregation and loss of normal nuclear function are implicated in amyotrophic lateral sclerosis (ALS) and frontotemporal dementia (FTD), making its RNA-recognition motifs (RRM1 and RRM2) a candidate target for small-molecule intervention. This project aims to identify fragment-sized compounds capable of engaging the TDP-43 RRM1/RRM2 RNA-binding interface through structure-based virtual screening. An ensemble docking strategy was used across five structurally diverse NMR conformers (PDB 4BS2) to account for the flexibility of the binding pocket, with a docking search space validated directly against literature-anchored pocket residues (Lukavsky et al., 2013; Qin et al., 2019) rather than relying on a single static structure. A GPU-accelerated docking pipeline (UniDock-Pro) was used to screen an 11,267-compound ChemDiv fragment library against the corrected binding-site grid, with hits ranked by predicted binding affinity and ligand efficiency, and cross-validated against known pocket residues to distinguish genuine binding-site engagement from scoring-function artifacts. The resulting hit list and residue-contact hotspot data are intended to inform future structure-activity and pharmacophore development for TDP-43-targeted small molecules.

## Structure
- `phase1/` — NMR ensemble structure prep, pocket diversity analysis, prepared receptor PDBQTs
- `conformation_docking/` — confirmation docking pipeline (20-ligand test set), grid verification and correction
- `phase2/` — full ChemDiv fragment library screen (11,267 compounds x 5 receptor conformers)
- `reports/` — write-ups: pocket diversity analysis, docking strategy, grid correction technical note

## Key results
- Corrected docking grid: center (24.64, 15.50, -11.43), size (36, 28, 25) A
- Literature pocket-residue contact rate: 86% -> 97% after grid correction (confirmation set)
- Phase 2 top hits ranked by score and ligand efficiency in `phase2/phase2_top50_overall.csv`

## Note on repo contents
Raw SDF libraries, individual ligand PDBQT files, and full docking pose outputs
are excluded (see `.gitignore`) as these are large generated datasets
(tens of thousands of files). Scripts here can regenerate them from the
original library files, which are stored on the lab GPU workstation (atlas).

## 
