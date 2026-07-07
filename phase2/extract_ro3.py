#!/usr/bin/env python3
"""
Extract only Rule-of-Three compliant molecules from a large SDF into a new,
much smaller SDF. Run on Atlas, next to the original file.

Usage:
    python3 extract_ro3.py screening_collection_3d.sdf ro3_compliant.sdf
"""

import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors

def main(in_path, out_path):
    n = 0
    n_written = 0
    writer = Chem.SDWriter(out_path)
    suppl = Chem.ForwardSDMolSupplier(in_path, sanitize=True, removeHs=False)

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

        if (mw <= 300 and logp <= 3 and hbd <= 3 and hba <= 3
                and rot <= 3 and psa <= 60):
            writer.write(mol)
            n_written += 1

        if n % 200000 == 0:
            print(f"...processed {n:,} | kept {n_written:,}", file=sys.stderr)

    writer.close()
    print(f"\nDone. Read {n:,} records, wrote {n_written:,} Ro3-compliant "
          f"molecules to {out_path}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 extract_ro3.py <input.sdf> <output.sdf>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
