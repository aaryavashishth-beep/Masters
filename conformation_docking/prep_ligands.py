#!/usr/bin/env python3
import os
import sys
import shutil
import subprocess
from pathlib import Path

# --- CONFIGURATION ---
# Set this to your massive DrugBank file
INPUT_SDF = Path("drugbank_all_3d_structures.sdf/3D structures.sdf") 
OUTPUT_DIR = Path("prepped_ligands")
INDEX_FILE = Path("ligand_index.txt")
MAX_LIGANDS = 20  # Keep this as 20 for your test run

def prepare_ligands():
    if not INPUT_SDF.exists():
        print(f"Error: Could not find {INPUT_SDF}")
        sys.exit(1)

    OUTPUT_DIR.mkdir(exist_ok=True)
    
    print(f"[*] Preparing {MAX_LIGANDS} ligands for UniDock-Pro...")
    
    # We use a temporary small file for the test
    test_sdf = Path("drugbank_test_20.sdf")
    
    # Extract first N molecules if not already done
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(str(INPUT_SDF))
    writer = Chem.SDWriter(str(test_sdf))
    count = 0
    for mol in suppl:
        if mol and count < MAX_LIGANDS:
            writer.write(mol)
            count += 1
        elif count >= MAX_LIGANDS:
            break
    writer.close()

    # Convert to PDBQT using OpenBabel
    # UniDock-Pro requires PDBQT or SDF with charges
    cmd = f"obabel {test_sdf} -O {OUTPUT_DIR}/ligand_.pdbqt -m -h --partialcharge gasteiger"
    subprocess.run(cmd, shell=True, check=True)

    # Generate the index file for UniDock-Pro
    with open(INDEX_FILE, "w") as f:
        for p in sorted(OUTPUT_DIR.glob("*.pdbqt")):
            f.write(f"{p.resolve()}\n")
            
    print(f"[+] Prep complete. {count} ligands ready.")

if __name__ == "__main__":
    prepare_ligands()
