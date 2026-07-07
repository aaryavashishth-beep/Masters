# TDP-43 RRM1/RRM2 Docking Project

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
