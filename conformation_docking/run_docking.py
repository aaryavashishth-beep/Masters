#!/usr/bin/env python3
import os
from pathlib import Path

UNIDOCK_BIN = "/home/s2831761/bin/udp"
RECEPTOR_DIR = Path("phase1/receptors")
LIGAND_INDEX = Path("ligand_index.txt")
OUTPUT_DIR = Path("docking_results")

def generate_sge():
    OUTPUT_DIR.mkdir(exist_ok=True)

    # CORRECTED grid — computed directly from Ca coordinates of all 15
    # literature pocket residues across all 5 receptor models, +5A padding.
    # Replaces the RNA-centroid-derived box, which excluded RRM1 on the Y axis.
    center_x, center_y, center_z = 24.64, 15.50, -11.43
    size_x, size_y, size_z = 36, 28, 25

    script_content = f"""#!/bin/bash
for model in model_01 model_02 model_05 model_11 model_19; do
    echo "Docking against ${{model}}..."
    {UNIDOCK_BIN} \\
        --receptor {RECEPTOR_DIR}/${{model}}.pdbqt \\
        --ligand_index {LIGAND_INDEX.resolve()} \\
        --center_x {center_x} --center_y {center_y} --center_z {center_z} \\
        --size_x {size_x} --size_y {size_y} --size_z {size_z} \\
        --search_mode balance \\
        --dir {OUTPUT_DIR.resolve()}/${{model}}
done
"""
    with open("submit_docking.sge", "w") as f:
        f.write(script_content)
    print("[+] Created submit_docking.sge with corrected grid")

if __name__ == "__main__":
    generate_sge()
