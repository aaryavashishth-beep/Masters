#!/usr/bin/env python3
"""
Streaming Rule-of-Three analysis for a large SDF library (e.g. Enamine
in-stock screening_collection_3d.sdf on Atlas/Apollo).

Run THIS on the cluster where the 16GB file actually lives — do not try
to transfer the file elsewhere. Requires RDKit:
    module load rdkit        # or: pip install rdkit --user
    python3 analyze_screening_collection.py screening_collection_3d.sdf

Streams record-by-record via RDKit's ForwardSDMolSupplier, so memory use
stays flat regardless of file size. Prints a compact summary you can paste
back into chat — no need to share the underlying file.
"""

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def main(path):
    n = 0
    n_parsed = 0
    n_ro3 = 0
    mw_sum = mw_lt300 = 0
    logp_sum = logp_ok = 0
    hbd_sum = hbd_ok = 0
    hba_sum = hba_ok = 0
    rot_sum = rot_ok = 0
    psa_sum = psa_ok = 0

    suppl = Chem.ForwardSDMolSupplier(path, sanitize=True, removeHs=False)
    for mol in suppl:
        n += 1
        if mol is None:
            continue
        try:
            mw   = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd  = Lipinski.NumHDonors(mol)
            hba  = Lipinski.NumHAcceptors(mol)
            rot  = Descriptors.NumRotatableBonds(mol)
            psa  = rdMolDescriptors.CalcTPSA(mol)
        except Exception:
            continue

        n_parsed += 1
        mw_sum += mw; logp_sum += logp; hbd_sum += hbd
        hba_sum += hba; rot_sum += rot; psa_sum += psa

        if mw <= 300:  mw_lt300 += 1
        if logp <= 3:  logp_ok += 1
        if hbd <= 3:   hbd_ok += 1
        if hba <= 3:   hba_ok += 1
        if rot <= 3:   rot_ok += 1
        if psa <= 60:  psa_ok += 1

        ro3 = (mw <= 300 and logp <= 3 and hbd <= 3 and hba <= 3
               and rot <= 3 and psa <= 60)
        if ro3:
            n_ro3 += 1

        if n % 200000 == 0:
            print(f"...processed {n:,} records", file=sys.stderr)

    print("\n=== Rule-of-Three Summary ===")
    print(f"Total records read      : {n:,}")
    print(f"Successfully parsed     : {n_parsed:,}")
    if n_parsed == 0:
        print("No molecules parsed — check file path / RDKit install.")
        return
    print(f"MW mean                 : {mw_sum/n_parsed:.1f}   (<=300: {100*mw_lt300/n_parsed:.1f}%)")
    print(f"logP mean               : {logp_sum/n_parsed:.2f}   (<=3: {100*logp_ok/n_parsed:.1f}%)")
    print(f"HBD mean                : {hbd_sum/n_parsed:.2f}   (<=3: {100*hbd_ok/n_parsed:.1f}%)")
    print(f"HBA mean                : {hba_sum/n_parsed:.2f}   (<=3: {100*hba_ok/n_parsed:.1f}%)")
    print(f"RotBonds mean           : {rot_sum/n_parsed:.2f}   (<=3: {100*rot_ok/n_parsed:.1f}%)")
    print(f"PSA mean                : {psa_sum/n_parsed:.1f}   (<=60: {100*psa_ok/n_parsed:.1f}%)")
    print(f"Full Ro3 pass (all 6)   : {n_ro3:,}  ({100*n_ro3/n_parsed:.1f}%)")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 analyze_screening_collection.py <path_to.sdf>")
        sys.exit(1)
    main(sys.argv[1])
