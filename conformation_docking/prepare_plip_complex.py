#!/usr/bin/env python3
"""
================================================================================
 PREPARE PLIP COMPLEX — merge a receptor + one docked ligand pose into a
 single complex.pdb that PLIP (local install or web server) can read.
================================================================================

Requires: OpenBabel's `obabel` on PATH.
  Check with: which obabel
  On atlas, ADFRsuite already provides this (confirmed earlier in this
  project — it's on your PATH already).

Usage:
    python3 prepare_plip_complex.py \
        --receptor phase1/receptors/model_01.pdbqt \
        --ligand-poses docking_results/model_01/ligand_12_out.pdbqt \
        --pose-rank 1 \
        --out complex_model01_ligand12.pdb

--pose-rank 1 = best-scoring pose (default). Use 2, 3... for other modes.
"""

import argparse
import os
import re
import subprocess
import sys
import tempfile


def extract_pose(pdbqt_path, rank, out_path):
    """Pull the rank-th best-scoring MODEL block out of a multi-pose PDBQT."""
    poses = []
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
                    poses.append((cur_score, cur_lines))
            else:
                cur_lines.append(line)

    if not poses:
        sys.exit(f"No scored poses found in {pdbqt_path} "
                  f"(expected MODEL/ENDMDL blocks with a RESULT: score line)")

    poses.sort(key=lambda p: p[0])  # most negative (best) score first
    if rank > len(poses):
        sys.exit(f"--pose-rank {rank} requested but only {len(poses)} poses present")

    score, lines = poses[rank - 1]
    with open(out_path, "w") as f:
        f.writelines(lines)
    print(f"Extracted pose rank {rank} (score {score} kcal/mol) -> {out_path}")
    return score


def obabel_to_pdb(in_path, out_path, resname=None):
    subprocess.run(["obabel", in_path, "-O", out_path],
                    check=True, capture_output=True)
    if resname:
        # Tag the ligand residue name so PLIP clearly IDs it as HETATM/LIG
        with open(out_path) as f:
            lines = f.readlines()
        with open(out_path, "w") as f:
            for line in lines:
                if line.startswith(("ATOM", "HETATM")):
                    line = "HETATM" + line[6:17] + f"{resname:<3}" + line[20:]
                f.write(line)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True,
                     help="Receptor .pdbqt (e.g. phase1/receptors/model_01.pdbqt)")
    ap.add_argument("--ligand-poses", required=True,
                     help="Docked ligand pose .pdbqt with multiple MODEL blocks")
    ap.add_argument("--pose-rank", type=int, default=1,
                     help="1 = best-scoring pose (default)")
    ap.add_argument("--out", required=True, help="Output complex PDB path")
    ap.add_argument("--ligand-resname", default="LIG")
    args = ap.parse_args()

    tmpdir = tempfile.mkdtemp()
    pose_pdbqt = os.path.join(tmpdir, "pose.pdbqt")
    extract_pose(args.ligand_poses, args.pose_rank, pose_pdbqt)

    receptor_pdb = os.path.join(tmpdir, "receptor.pdb")
    ligand_pdb = os.path.join(tmpdir, "ligand.pdb")
    obabel_to_pdb(args.receptor, receptor_pdb)
    obabel_to_pdb(pose_pdbqt, ligand_pdb, resname=args.ligand_resname)

    with open(args.out, "w") as out:
        with open(receptor_pdb) as f:
            for line in f:
                if line.startswith(("ATOM", "TER")):
                    out.write(line)
        with open(ligand_pdb) as f:
            for line in f:
                if line.startswith("HETATM"):
                    out.write(line)
        out.write("END\n")

    print(f"\nComplex written to: {args.out}")
    print("\nOption A — run PLIP locally on atlas:")
    print(f"    pip install plip --break-system-packages   # one-time")
    print(f"    plip -f {args.out} -o plip_report/ --pymol")
    print("\nOption B — no install, use the web server:")
    print(f"    scp s2831761@atlas.bch.ed.ac.uk:{os.path.abspath(args.out)} ~/Downloads/")
    print("    then upload the file at https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index")


if __name__ == "__main__":
    main()
