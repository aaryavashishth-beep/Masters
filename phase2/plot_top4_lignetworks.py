#!/usr/bin/env python3
"""
================================================================================
 TOP 4 LIGNETWORK DIAGRAMS — interactive 2D interaction diagrams
================================================================================

Reuses the exact same, already-validated pipeline from batch_prolif_analysis.py
(pose extraction, coordinate transplant onto correctly-bonded library
molecules, receptor parsing) to regenerate fingerprints for the top 4 hits
that already succeeded in the batch run, then renders each as a standalone
interactive HTML network diagram via ProLIF's LigNetwork.

Usage:
    source prolif_env/bin/activate
    python3 plot_top4_lignetworks.py \
        --per-hit-csv prolif_batch/per_hit_interactions.csv \
        --receptor-dir receptors \
        --results-dir docking_results \
        --library-sdf Download_SDF_Fragments_Library_11269.sdf \
        --out-dir prolif_batch/lignetworks

Outputs one HTML file per hit, e.g.:
    prolif_batch/lignetworks/1_model_02_ligand_3805.html
Open any of these directly in a browser -- fully interactive, no server needed.
"""

import argparse
import csv
from pathlib import Path

import prolif as plf
from rdkit import Chem

from batch_prolif_analysis import (
    extract_pose,
    obabel_pdbqt_to_sdf,
    obabel_pdbqt_to_pdb,
    transplant_coords,
    load_library_molecules,
)
import re
import tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-hit-csv", default="prolif_batch/per_hit_interactions.csv")
    ap.add_argument("--receptor-dir", default="receptors")
    ap.add_argument("--results-dir", default="docking_results")
    ap.add_argument("--library-sdf", default="Download_SDF_Fragments_Library_11269.sdf")
    ap.add_argument("--out-dir", default="prolif_batch/lignetworks")
    ap.add_argument("--top-n", type=int, default=4)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.per_hit_csv) as f:
        rows = list(csv.DictReader(f))

    # Only rows that actually succeeded (have real interactions, not "NONE")
    successful = [r for r in rows if r["interactions"] != "NONE"]
    successful.sort(key=lambda r: float(r["score_kcal_mol"]))
    top_n = successful[: args.top_n]

    if not top_n:
        raise SystemExit("No successful hits found in per-hit CSV -- rerun batch_prolif_analysis.py first")

    print(f"[*] Generating LigNetwork diagrams for top {len(top_n)} hits:")
    for r in top_n:
        print(f"    {r['receptor_model']} {r['ligand']}  score={r['score_kcal_mol']}")

    def ligand_index(ligand_name):
        n = int(re.search(r"ligand_(\d+)", ligand_name).group(1))
        return n - 1

    needed_indices = {ligand_index(r["ligand"]) for r in top_n}
    library_mols = load_library_molecules(args.library_sdf, needed_indices)

    receptor_cache = {}
    for model in sorted({r["receptor_model"] for r in top_n}):
        receptor_pdbqt = Path(args.receptor_dir) / f"{model}.pdbqt"
        with tempfile.TemporaryDirectory() as tmpdir:
            receptor_pdb = Path(tmpdir) / "receptor.pdb"
            obabel_pdbqt_to_pdb(receptor_pdbqt, receptor_pdb)
            protein_rdkit = Chem.MolFromPDBFile(str(receptor_pdb), removeHs=False)
        receptor_cache[model] = plf.Molecule(protein_rdkit)

    fp_engine = plf.Fingerprint(["HBDonor", "HBAcceptor", "PiStacking",
                                  "Hydrophobic", "PiCation", "CationPi", "Anionic", "Cationic"])

    for rank, r in enumerate(top_n, start=1):
        model = r["receptor_model"]
        ligand = r["ligand"]
        score = float(r["score_kcal_mol"])
        idx = ligand_index(ligand)
        orig_mol = library_mols[idx]

        pose_pdbqt = Path(args.results_dir) / model / f"{ligand}.pdbqt"

        with tempfile.TemporaryDirectory() as tmpdir:
            extracted_pose = Path(tmpdir) / "pose.pdbqt"
            extract_pose(pose_pdbqt, score, extracted_pose)

            pose_sdf = Path(tmpdir) / "pose.sdf"
            obabel_pdbqt_to_sdf(extracted_pose, pose_sdf)

            mol_with_coords = transplant_coords(orig_mol, pose_sdf)

        lig_mol = plf.Molecule(mol_with_coords)
        prot_mol = receptor_cache[model]

        fp_engine.run_from_iterable([lig_mol], prot_mol)

        net = fp_engine.plot_lignetwork(lig_mol, kind="frame", frame=0)
        out_path = out_dir / f"{rank}_{model}_{ligand.replace('_out', '')}.html"
        net.save(str(out_path))
        print(f"[{rank}/{len(top_n)}] saved -> {out_path}")

    print(f"\n[+] Done. Open any .html file in {out_dir} directly in a browser.")


if __name__ == "__main__":
    main()
