#!/usr/bin/env python3
"""
================================================================================
 BATCH PROLIF INTERACTION FINGERPRINT — Top 50 Phase 2 hits
================================================================================

Validated pipeline (confirmed working manually before this script was written):
  1. Extract the best-scoring MODEL block from the docked pose PDBQT
  2. Convert it to SDF via obabel (coordinates only -- do NOT trust its
     guessed bond orders, PDBQT/naive-SDF round-trips can invent impossible
     valences)
  3. Load the ORIGINAL, correctly-bonded molecule from the source ChemDiv
     SDF library (validated by RDKit during prep_ligands.py) and transplant
     the docked pose's 3D coordinates onto it -- this gives real chemistry
     with real docked geometry
  4. Parse the receptor with RDKit's PDB parser (template-based, not
     distance-guessed bonds)
  5. Run ProLIF (RDKit/MDAnalysis-based -- no OpenBabel-Python dependency,
     sidesteps the venv/plugin issues PLIP's CLI hit) for HBDonor,
     HBAcceptor, PiStacking, and Hydrophobic interactions
  6. Tally a consensus interaction fingerprint across all 50 hits

Requirements (already confirmed installed in prolif_env):
    pip install prolif

Usage:
    source prolif_env/bin/activate
    python3 batch_prolif_analysis.py \
        --top50-csv phase2_top50_overall.csv \
        --receptor-dir receptors \
        --results-dir docking_results \
        --library-sdf Download_SDF_Fragments_Library_11269.sdf \
        --out-dir prolif_batch
"""

import argparse
import csv
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import prolif as plf
from rdkit import Chem

OBABEL_BIN = "/usr/local/bin/obabel"


def extract_pose(pdbqt_path, target_score, out_path):
    """Pull the MODEL block whose score matches target_score most closely."""
    import subprocess
    best_diff, best_lines = None, None
    cur_lines, cur_score = [], None
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                cur_lines, cur_score = [], None
            elif "RESULT:" in line:
                m = re.search(r"RESULT:\s*(-?\d+\.?\d*)", line)
                if m:
                    cur_score = float(m.group(1))
                cur_lines.append(line)
            elif line.startswith("ENDMDL"):
                if cur_score is not None:
                    diff = abs(cur_score - target_score)
                    if best_diff is None or diff < best_diff:
                        best_diff, best_lines = diff, list(cur_lines)
            else:
                cur_lines.append(line)
    if best_lines is None:
        return False
    with open(out_path, "w") as f:
        f.writelines(best_lines)
    return True


def obabel_pdbqt_to_sdf(pdbqt_path, sdf_path):
    import subprocess
    subprocess.run([OBABEL_BIN, str(pdbqt_path), "-O", str(sdf_path)],
                    check=True, capture_output=True)


def obabel_pdbqt_to_pdb(pdbqt_path, pdb_path):
    import subprocess
    subprocess.run([OBABEL_BIN, str(pdbqt_path), "-O", str(pdb_path)],
                    check=True, capture_output=True)


def transplant_coords(orig_mol, pose_sdf_path):
    """
    Return a copy of orig_mol with coordinates from pose_sdf_path.

    Matches on HEAVY ATOMS ONLY: PDBQT is a united-atom format (only polar
    H are explicit; nonpolar H are merged into their parent atom), while
    the source library SDF is heavy-atom-only. Comparing raw atom counts
    between the two is unreliable -- any molecule with a polar H (-OH,
    -NH, etc.) mismatches by exactly that many atoms. Instead: strip both
    down to heavy atoms, match on that, transplant those coordinates, then
    let RDKit regenerate reasonable hydrogen positions from the resulting
    3D heavy-atom geometry.
    """
    pose_mol = Chem.MolFromMolFile(str(pose_sdf_path), sanitize=False)
    if pose_mol is None:
        return None

    orig_heavy = Chem.RemoveHs(orig_mol, sanitize=False)
    pose_heavy = Chem.RemoveHs(pose_mol, sanitize=False)

    if orig_heavy.GetNumAtoms() != pose_heavy.GetNumAtoms():
        return None

    mol = Chem.Mol(orig_heavy)
    conf = pose_heavy.GetConformer()
    new_conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        new_conf.SetAtomPosition(i, conf.GetAtomPosition(i))
    mol.RemoveAllConformers()
    mol.AddConformer(new_conf)

    try:
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol, addCoords=True)
    except Exception:
        return None

    return mol


def load_library_molecules(sdf_path, needed_indices):
    """Single pass through the (large) library SDF, keeping only what's needed."""
    needed = set(needed_indices)
    found = {}
    suppl = Chem.SDMolSupplier(str(sdf_path))
    for idx, mol in enumerate(suppl):
        if idx in needed:
            found[idx] = mol
            if len(found) == len(needed):
                break
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top50-csv", default="phase2_top50_overall.csv")
    ap.add_argument("--receptor-dir", default="receptors")
    ap.add_argument("--results-dir", default="docking_results")
    ap.add_argument("--library-sdf", default="Download_SDF_Fragments_Library_11269.sdf")
    ap.add_argument("--out-dir", default="prolif_batch")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    with open(args.top50_csv) as f:
        hits = list(csv.DictReader(f))

    def ligand_index(ligand_name):
        # "ligand_3805_out" -> 3804 (0-indexed position in the source SDF)
        n = int(re.search(r"ligand_(\d+)", ligand_name).group(1))
        return n - 1

    print("[*] Loading required molecules from library SDF (single pass)...")
    needed_indices = {ligand_index(h["ligand"]) for h in hits}
    library_mols = load_library_molecules(args.library_sdf, needed_indices)
    print(f"    Loaded {len(library_mols)}/{len(needed_indices)} unique ligands")

    print("[*] Parsing receptor structures (once per model)...")
    receptor_cache = {}
    for model in sorted({h["receptor_model"] for h in hits}):
        receptor_pdbqt = Path(args.receptor_dir) / f"{model}.pdbqt"
        with tempfile.TemporaryDirectory() as tmpdir:
            receptor_pdb = Path(tmpdir) / "receptor.pdb"
            obabel_pdbqt_to_pdb(receptor_pdbqt, receptor_pdb)
            protein_rdkit = Chem.MolFromPDBFile(str(receptor_pdb), removeHs=False)
        if protein_rdkit is None:
            sys.exit(f"Failed to parse receptor {model} -- cannot continue")
        receptor_cache[model] = plf.Molecule(protein_rdkit)
        print(f"    {model}: {protein_rdkit.GetNumAtoms()} atoms")

    fp_engine = plf.Fingerprint(["HBDonor", "HBAcceptor", "PiStacking",
                                  "Hydrophobic", "PiCation", "CationPi", "Anionic", "Cationic"])

    per_hit_rows = []
    hotspot = defaultdict(lambda: defaultdict(int))
    failures = []

    for i, hit in enumerate(hits, start=1):
        model = hit["receptor_model"]
        ligand = hit["ligand"]
        score = float(hit["score_kcal_mol"])
        tag = f"{model}_{ligand}"
        print(f"[{i}/{len(hits)}] {tag} (score {score})")

        idx = ligand_index(ligand)
        orig_mol = library_mols.get(idx)
        if orig_mol is None:
            failures.append((tag, "not found in library SDF"))
            continue

        pose_pdbqt = Path(args.results_dir) / model / f"{ligand}.pdbqt"
        if not pose_pdbqt.exists():
            failures.append((tag, "pose file missing"))
            continue

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                extracted_pose = Path(tmpdir) / "pose.pdbqt"
                if not extract_pose(pose_pdbqt, score, extracted_pose):
                    failures.append((tag, "could not extract matching pose"))
                    continue

                pose_sdf = Path(tmpdir) / "pose.sdf"
                obabel_pdbqt_to_sdf(extracted_pose, pose_sdf)

                mol_with_coords = transplant_coords(orig_mol, pose_sdf)
                if mol_with_coords is None:
                    failures.append((tag, "atom count mismatch, coordinate transplant failed"))
                    continue

                lig_mol = plf.Molecule(mol_with_coords)
                prot_mol = receptor_cache[model]

                fp_engine.run_from_iterable([lig_mol], prot_mol)
                df = fp_engine.to_dataframe()

        except Exception as e:
            failures.append((tag, f"unexpected error: {e}"))
            continue

        interactions_found = []
        if not df.empty:
            for col in df.columns:
                _, residue, itype = col
                if df[col].iloc[0]:
                    interactions_found.append((residue, itype))
                    hotspot[residue][itype] += 1

        per_hit_rows.append({
            "receptor_model": model,
            "ligand": ligand,
            "score_kcal_mol": score,
            "interactions": ";".join(f"{r}:{t}" for r, t in interactions_found) or "NONE",
            "n_interactions": len(interactions_found),
        })

    with open(out_dir / "per_hit_interactions.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["receptor_model", "ligand", "score_kcal_mol",
                                            "interactions", "n_interactions"])
        w.writeheader()
        w.writerows(per_hit_rows)

    all_types = sorted({t for res in hotspot for t in hotspot[res]})
    with open(out_dir / "consensus_hotspot_by_type.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["residue"] + all_types + ["total"])
        for res in sorted(hotspot, key=lambda r: -sum(hotspot[r].values())):
            row = [res] + [hotspot[res].get(t, 0) for t in all_types]
            row.append(sum(hotspot[res].values()))
            w.writerow(row)

    with open(out_dir / "failures.log", "w") as f:
        for tag, reason in failures:
            f.write(f"{tag}\t{reason}\n")

    print(f"\n[+] Done. {len(per_hit_rows)}/{len(hits)} hits successfully analyzed.")
    if failures:
        print(f"    {len(failures)} failures logged in {out_dir/'failures.log'}")
    print(f"    Per-hit interactions: {out_dir/'per_hit_interactions.csv'}")
    print(f"    Consensus hotspot map: {out_dir/'consensus_hotspot_by_type.csv'}")

    print("\nTop 10 consensus hotspot residues:")
    with open(out_dir / "consensus_hotspot_by_type.csv") as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in list(reader)[:10]:
            nonzero = [f"{h}={v}" for h, v in zip(header[1:-1], row[1:-1]) if v != "0"]
            print(f"  {row[0]:12s} total={row[-1]:>3s}  " + "  ".join(nonzero))


if __name__ == "__main__":
    main()
