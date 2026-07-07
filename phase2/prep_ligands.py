#!/usr/bin/env python3
"""
Ligand preparation for Phase 2 ChemDiv screening.

Converts each molecule individually (RDKit split -> per-molecule OpenBabel
call) rather than one giant batch OpenBabel run. This means a single bad
or slow-converting molecule gets logged and skipped instead of silently
killing the entire conversion partway through (as happened with the
batch approach: it stopped at exactly 457/11269 with no error message).

Usage:
    python3 prep_ligands.py

Outputs:
    prepped_ligands/ligand_<n>.pdbqt   -- one file per successfully converted molecule
    ligand_index.txt                   -- list of paths, for run_docking.py
    prep_failures.log                  -- which molecule indices failed and why
"""

import subprocess
import sys
import tempfile
from pathlib import Path
from rdkit import Chem

INPUT_SDF = Path("Download_SDF_Fragments_Library_11269.sdf")
OUTPUT_DIR = Path("prepped_ligands")
INDEX_FILE = Path("ligand_index.txt")
FAIL_LOG = Path("prep_failures.log")
TIMEOUT_SECONDS = 20  # generous ceiling per molecule; normal conversion takes <1s


def prepare_ligands():
    if not INPUT_SDF.exists():
        sys.exit(f"Error: Could not find {INPUT_SDF}")

    OUTPUT_DIR.mkdir(exist_ok=True)

    print(f"[*] Reading molecules from {INPUT_SDF} ...")
    suppl = Chem.SDMolSupplier(str(INPUT_SDF))

    n_total = 0
    n_rdkit_invalid = 0
    n_converted = 0
    n_failed = 0
    failures = []

    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, mol in enumerate(suppl, start=1):
            n_total += 1
            if idx % 500 == 0:
                print(f"    ... processed {idx} molecules "
                      f"({n_converted} ok, {n_failed} failed, {n_rdkit_invalid} rdkit-invalid)")

            if mol is None:
                n_rdkit_invalid += 1
                failures.append((idx, "rdkit_parse_failed"))
                continue

            tmp_sdf = Path(tmpdir) / f"mol_{idx}.sdf"
            out_pdbqt = OUTPUT_DIR / f"ligand_{idx}.pdbqt"

            writer = Chem.SDWriter(str(tmp_sdf))
            writer.write(mol)
            writer.close()

            cmd = ["obabel", str(tmp_sdf), "-O", str(out_pdbqt),
                   "-h", "--partialcharge", "gasteiger"]
            try:
                result = subprocess.run(cmd, capture_output=True, text=True,
                                          timeout=TIMEOUT_SECONDS)
                if result.returncode == 0 and out_pdbqt.exists():
                    # Exit code 0 and file existing isn't enough -- OpenBabel
                    # can write a structurally empty PDBQT (0 ATOM records)
                    # and still return success. Verify real atom content.
                    has_atoms = False
                    with open(out_pdbqt) as pf:
                        for line in pf:
                            if line.startswith("ATOM") or line.startswith("HETATM"):
                                has_atoms = True
                                break
                    if has_atoms:
                        n_converted += 1
                    else:
                        n_failed += 1
                        out_pdbqt.unlink(missing_ok=True)
                        failures.append((idx, "empty_pdbqt_no_atoms"))
                else:
                    n_failed += 1
                    failures.append((idx, f"obabel_error: {result.stderr.strip()[:200]}"))
            except subprocess.TimeoutExpired:
                n_failed += 1
                failures.append((idx, f"timeout_after_{TIMEOUT_SECONDS}s"))

    with open(FAIL_LOG, "w") as f:
        for idx, reason in failures:
            f.write(f"{idx}\t{reason}\n")

    with open(INDEX_FILE, "w") as f:
        for p in sorted(OUTPUT_DIR.glob("*.pdbqt"),
                          key=lambda p: int(p.stem.split("_")[1])):
            f.write(f"{p.resolve()}\n")

    print(f"\n[+] Done. {n_total} molecules read from input.")
    print(f"    Converted successfully: {n_converted}")
    print(f"    RDKit-invalid (skipped): {n_rdkit_invalid}")
    print(f"    OpenBabel failed/timeout: {n_failed}")
    print(f"    Index written to: {INDEX_FILE}")
    if failures:
        print(f"    Failure details in: {FAIL_LOG}")


if __name__ == "__main__":
    prepare_ligands()
