#!/usr/bin/env python3
"""
================================================================================
 CONTACT-RESIDUE ANALYSIS — TDP-43 RRM1/RRM2 Docking Poses
 Checks which receptor residues each docked ligand pose actually contacts,
 and flags whether those contacts match the literature-defined pocket.
================================================================================

Usage:
    python3 analyze_contacts.py \
        --receptor-dir ~/conformation_docking/phase1/receptors \
        --results-dir  ~/conformation_docking/docking_results \
        --cutoff 4.0 \
        --top-n 1

What it does:
  1. Parses each receptor PDBQT (model_01.pdbqt, model_02.pdbqt, ...) to get
     per-atom coordinates tagged with residue name/number.
  2. Parses every ligand pose PDBQT produced by udp under docking_results/,
     reading the Vina-style "REMARK VINA RESULT" score for each MODEL block.
  3. For the top-N poses (by score) of each ligand, finds every receptor
     residue with at least one heavy atom within --cutoff Å of any ligand
     atom.
  4. Cross-references contacted residues against the literature pocket
     lists (RRM1: Lukavsky et al. 2013; RRM2: Qin et al. 2019) and reports
     hit/miss per pose.
  5. Writes a CSV summary and prints a short console report.

No third-party dependencies — pure stdlib.
"""

import argparse
import csv
import glob
import math
import os
import re
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Literature-defined pocket residues (see Phase1_Pocket_Diversity_Analysis /
# TDP43_Docking_Strategy_Report for citations)
# ---------------------------------------------------------------------------
RRM1_POCKET = {113: "Trp113", 147: "Phe147", 149: "Phe149", 151: "Arg151",
               152: "Asn152", 181: "Lys181", 183: "Asp183"}
RRM2_POCKET = {245: "Gly245", 246: "Glu246", 247: "Asp247", 248: "Lys248",
               256: "His256", 257: "Ile257", 258: "Ser258", 259: "Asn259"}
LIT_POCKET = {**RRM1_POCKET, **RRM2_POCKET}


def parse_receptor(pdbqt_path):
    """Return dict: resnum -> (resname, [(x,y,z), ...] heavy-atom coords)."""
    residues = defaultdict(lambda: [None, []])
    with open(pdbqt_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            residues[resnum][0] = resname
            residues[resnum][1].append((x, y, z))
    return residues


def parse_ligand_poses(pdbqt_path):
    """
    Yield (model_index, score, [(x,y,z), ...]) for each MODEL block in a
    Vina/UniDock-style multi-pose PDBQT output file.
    """
    model_idx = None
    score = None
    coords = []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                model_idx = int(line.split()[1])
                score = None
                coords = []
            elif line.startswith("REMARK VINA RESULT") or "RESULT:" in line:
                m = re.search(r"RESULT:\s*(-?\d+\.?\d*)", line)
                if m:
                    score = float(m.group(1))
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    pass
            elif line.startswith("ENDMDL"):
                if model_idx is not None and score is not None:
                    yield model_idx, score, coords
                model_idx = None

    # Fallback: single-pose file with no MODEL/ENDMDL wrapper
    if model_idx is None and coords == [] :
        pass


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def contacted_residues(ligand_coords, receptor_residues, cutoff):
    hits = {}
    for resnum, (resname, atoms) in receptor_residues.items():
        for latom in ligand_coords:
            close = False
            for ratom in atoms:
                if dist(latom, ratom) <= cutoff:
                    close = True
                    break
            if close:
                hits[resnum] = resname
                break
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--cutoff", type=float, default=4.0,
                     help="Contact distance cutoff in Angstroms (default 4.0)")
    ap.add_argument("--top-n", type=int, default=1,
                     help="Number of top-scoring poses per ligand to check (default 1)")
    ap.add_argument("--out-csv", default="contact_summary.csv")
    args = ap.parse_args()

    receptor_files = sorted(glob.glob(os.path.join(args.receptor_dir, "*.pdbqt")))
    if not receptor_files:
        sys.exit(f"No receptor .pdbqt files found in {args.receptor_dir}")

    rows = []
    for rec_path in receptor_files:
        model_name = os.path.splitext(os.path.basename(rec_path))[0]
        receptor_residues = parse_receptor(rec_path)
        model_result_dir = os.path.join(args.results_dir, model_name)
        if not os.path.isdir(model_result_dir):
            print(f"[skip] no results dir for {model_name}: {model_result_dir}")
            continue

        pose_files = sorted(glob.glob(os.path.join(model_result_dir, "*.pdbqt")))
        for pose_path in pose_files:
            ligand_name = os.path.splitext(os.path.basename(pose_path))[0]
            poses = sorted(parse_ligand_poses(pose_path), key=lambda t: t[1])[:args.top_n]
            for rank, (midx, score, coords) in enumerate(poses, start=1):
                if not coords:
                    continue
                hits = contacted_residues(coords, receptor_residues, args.cutoff)
                # Use the literature label (e.g. "Gly245"), not the raw PDBQT
                # resname (e.g. "GLY"), so it always carries a residue number.
                lit_hits = {n: LIT_POCKET[n] for n in hits if n in LIT_POCKET}
                rrm1_hit = any(n in RRM1_POCKET for n in hits)
                rrm2_hit = any(n in RRM2_POCKET for n in hits)
                domain = "RRM1" if rrm1_hit and not rrm2_hit else \
                         "RRM2" if rrm2_hit and not rrm1_hit else \
                         "RRM1+RRM2" if rrm1_hit and rrm2_hit else "neither"

                rows.append({
                    "receptor_model": model_name,
                    "ligand": ligand_name,
                    "pose_rank": rank,
                    "score_kcal_mol": score,
                    "n_contacted_residues": len(hits),
                    "literature_pocket_contacts": ";".join(
                        sorted(lit_hits.values(), key=lambda s: int(re.sub(r"\D", "", s)))
                    ) or "NONE",
                    "domain_matched": domain,
                    "all_contacted_residues": ";".join(
                        f"{nm}{n}" for n, nm in sorted(hits.items())
                    ),
                })

    if not rows:
        sys.exit("No poses parsed — check --receptor-dir / --results-dir paths and "
                 "that pose PDBQTs contain REMARK VINA RESULT lines.")

    with open(args.out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # ---- Aggregated hotspot frequency (residue -> how often it's contacted)
    freq = defaultdict(int)
    for r in rows:
        if r["literature_pocket_contacts"] == "NONE":
            continue
        for label in r["literature_pocket_contacts"].split(";"):
            freq[label] += 1

    freq_csv = os.path.splitext(args.out_csv)[0] + "_residue_frequency.csv"
    with open(freq_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["residue", "resnum", "domain", "n_poses_contacted",
                          "pct_of_all_poses"])
        for resnum, label in sorted(LIT_POCKET.items()):
            n = freq.get(label, 0)
            domain = "RRM1" if resnum in RRM1_POCKET else "RRM2"
            writer.writerow([label, resnum, domain, n,
                              round(100 * n / len(rows), 1)])
    print(f"Per-residue hotspot frequency written to: {freq_csv}")

    # Console summary
    n_total = len(rows)
    n_lit = sum(1 for r in rows if r["literature_pocket_contacts"] != "NONE")
    n_neither = sum(1 for r in rows if r["domain_matched"] == "neither")
    print(f"\nAnalyzed {n_total} top-{args.top_n} poses across "
          f"{len(receptor_files)} receptor models.")
    print(f"Poses touching >=1 literature pocket residue: {n_lit}/{n_total} "
          f"({100*n_lit/n_total:.1f}%)")
    print(f"Poses touching NEITHER RRM1 nor RRM2 pocket list: {n_neither}/{n_total}")
    print(f"\nFull results written to: {args.out_csv}")
    print("\nSample rows:")
    for r in rows[:5]:
        print(f"  {r['receptor_model']:10s} {r['ligand']:15s} "
              f"score={r['score_kcal_mol']:>7.2f}  "
              f"domain={r['domain_matched']:10s}  "
              f"lit_contacts={r['literature_pocket_contacts']}")


if __name__ == "__main__":
    main()
