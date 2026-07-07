#!/usr/bin/env python3
"""
================================================================================
 PHASE 2 RESULTS SUMMARY — TDP-43 ChemDiv Fragment Screen (11,267 ligands x 5
 receptor models = 56,335 top-1 poses)
================================================================================

Parses every docked pose once, then produces three ranked views:

  1. phase2_top50_overall.csv     -- best 50 scores across all models/ligands
  2. phase2_top_per_model.csv     -- best N (default 20) scores within each
                                     of the 5 receptor models
  3. phase2_literature_hits.csv   -- every pose that contacts >=1 literature
                                     pocket residue (RRM1/RRM2), sorted by score

  phase2_all_results.csv          -- the full underlying table (56,335 rows),
                                     in case you need to re-slice it later

Usage:
    python3 summarize_phase2.py \
        --receptor-dir receptors \
        --results-dir docking_results \
        --cutoff 4.0 \
        --top-n-per-model 20 \
        --top-n-overall 50

Expect this to take a while (tens of thousands of poses, pure-Python distance
checks) -- run it in tmux, not in the foreground:
    tmux new -s summarize -d 'python3 summarize_phase2.py > summarize.log 2>&1'
"""

import argparse
import csv
import glob
import math
import os
import re
import sys
import time
from collections import defaultdict

RRM1_POCKET = {113: "Trp113", 147: "Phe147", 149: "Phe149", 151: "Arg151",
               152: "Asn152", 181: "Lys181", 183: "Asp183"}
RRM2_POCKET = {245: "Gly245", 246: "Glu246", 247: "Asp247", 248: "Lys248",
               256: "His256", 257: "Ile257", 258: "Ser258", 259: "Asn259"}
LIT_POCKET = {**RRM1_POCKET, **RRM2_POCKET}


def parse_receptor(pdbqt_path):
    residues = defaultdict(lambda: [None, []])
    with open(pdbqt_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                resname = line[17:20].strip()
                resnum = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            residues[resnum][0] = resname
            residues[resnum][1].append((x, y, z))
    return residues


def best_pose(pdbqt_path):
    """Return (score, [(x,y,z), ...]) for the best-scoring MODEL block."""
    best_score, best_coords = None, None
    cur_score, cur_coords = None, []
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL"):
                cur_score, cur_coords = None, []
            elif "RESULT:" in line:
                m = re.search(r"RESULT:\s*(-?\d+\.?\d*)", line)
                if m:
                    cur_score = float(m.group(1))
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                    cur_coords.append((x, y, z))
                except ValueError:
                    pass
            elif line.startswith("ENDMDL"):
                if cur_score is not None and (best_score is None or cur_score < best_score):
                    best_score, best_coords = cur_score, cur_coords
    return best_score, best_coords or []


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def contacted_residues(ligand_coords, receptor_residues, cutoff):
    hits = {}
    for resnum, (resname, atoms) in receptor_residues.items():
        for latom in ligand_coords:
            if any(dist(latom, ratom) <= cutoff for ratom in atoms):
                hits[resnum] = resname
                break
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor-dir", default="receptors")
    ap.add_argument("--results-dir", default="docking_results")
    ap.add_argument("--cutoff", type=float, default=4.0)
    ap.add_argument("--top-n-per-model", type=int, default=20)
    ap.add_argument("--top-n-overall", type=int, default=50)
    args = ap.parse_args()

    receptor_files = sorted(glob.glob(os.path.join(args.receptor_dir, "*.pdbqt")))
    if not receptor_files:
        sys.exit(f"No receptor .pdbqt files found in {args.receptor_dir}")

    all_rows = []
    t0 = time.time()

    for rec_path in receptor_files:
        model_name = os.path.splitext(os.path.basename(rec_path))[0]
        receptor_residues = parse_receptor(rec_path)
        model_dir = os.path.join(args.results_dir, model_name)
        if not os.path.isdir(model_dir):
            print(f"[skip] no results dir for {model_name}")
            continue

        pose_files = sorted(glob.glob(os.path.join(model_dir, "*.pdbqt")))
        print(f"[*] {model_name}: {len(pose_files)} ligand pose files")

        for i, pose_path in enumerate(pose_files, start=1):
            if i % 1000 == 0:
                elapsed = time.time() - t0
                print(f"    ... {model_name}: {i}/{len(pose_files)} "
                      f"({elapsed/60:.1f} min elapsed)")

            ligand_name = os.path.splitext(os.path.basename(pose_path))[0]
            score, coords = best_pose(pose_path)
            if score is None or not coords:
                continue

            hits = contacted_residues(coords, receptor_residues, args.cutoff)
            lit_hits = {n: LIT_POCKET[n] for n in hits if n in LIT_POCKET}
            rrm1_hit = any(n in RRM1_POCKET for n in hits)
            rrm2_hit = any(n in RRM2_POCKET for n in hits)
            domain = "RRM1+RRM2" if rrm1_hit and rrm2_hit else \
                     "RRM1" if rrm1_hit else "RRM2" if rrm2_hit else "neither"

            all_rows.append({
                "receptor_model": model_name,
                "ligand": ligand_name,
                "score_kcal_mol": score,
                "n_contacted_residues": len(hits),
                "literature_pocket_contacts": ";".join(
                    sorted(lit_hits.values(), key=lambda s: int(re.sub(r"\D", "", s)))
                ) or "NONE",
                "domain_matched": domain,
            })

    if not all_rows:
        sys.exit("No poses parsed -- check --receptor-dir / --results-dir")

    fieldnames = list(all_rows[0].keys())

    # 1. Full underlying table
    with open("phase2_all_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(all_rows)

    # 2. Top N overall by score (most negative first)
    top_overall = sorted(all_rows, key=lambda r: r["score_kcal_mol"])[:args.top_n_overall]
    with open("phase2_top50_overall.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(top_overall)

    # 3. Top N per model
    by_model = defaultdict(list)
    for r in all_rows:
        by_model[r["receptor_model"]].append(r)
    with open("phase2_top_per_model.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for model_name in sorted(by_model):
            top_n = sorted(by_model[model_name], key=lambda r: r["score_kcal_mol"])[:args.top_n_per_model]
            w.writerows(top_n)

    # 4. Every pose contacting a literature residue, sorted by score
    lit_hits_rows = [r for r in all_rows if r["literature_pocket_contacts"] != "NONE"]
    lit_hits_rows.sort(key=lambda r: r["score_kcal_mol"])
    with open("phase2_literature_hits.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(lit_hits_rows)

    elapsed = time.time() - t0
    print(f"\n[+] Done in {elapsed/60:.1f} minutes.")
    print(f"    Total poses analyzed: {len(all_rows)}")
    print(f"    Poses contacting >=1 literature residue: {len(lit_hits_rows)} "
          f"({100*len(lit_hits_rows)/len(all_rows):.1f}%)")
    print(f"\n    phase2_all_results.csv       -- full table ({len(all_rows)} rows)")
    print(f"    phase2_top50_overall.csv     -- best {args.top_n_overall} scores overall")
    print(f"    phase2_top_per_model.csv     -- best {args.top_n_per_model} per model (5 models)")
    print(f"    phase2_literature_hits.csv   -- {len(lit_hits_rows)} literature-contacting poses, ranked")

    print("\n    Top 10 overall:")
    for r in top_overall[:10]:
        print(f"      {r['receptor_model']:10s} {r['ligand']:18s} "
              f"score={r['score_kcal_mol']:>7.2f}  domain={r['domain_matched']:10s}  "
              f"lit={r['literature_pocket_contacts']}")


if __name__ == "__main__":
    main()
