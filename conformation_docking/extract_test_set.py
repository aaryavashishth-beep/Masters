#!/usr/bin/env python3
from rdkit import Chem
from rdkit.Chem import Descriptors

# Update this path to where your unzipped DrugBank file actually lives
input_sdf = "drugbank_all_3d_structures.sdf/3D structures.sdf" 
output_sdf = "drugbank_test_20.sdf"

supplier = Chem.SDMolSupplier(input_sdf)
writer = Chem.SDWriter(output_sdf)

count = 0
print("Extracting 20 small compounds from DrugBank dataset...")

for mol in supplier:
    if mol is None:
        continue
        
    try:
        mw = Descriptors.MolWt(mol)
        # Filter for true small fragment properties (MW <= 250)
        if mw <= 250.0:
            writer.write(mol)
            count += 1
            print(f"  [{count}/20] Retained Mol ID: {mol.GetProp('_Name')} (MW: {mw:.2f})")
    except Exception:
        continue
        
    if count >= 20:
        break

writer.close()
print(f"\nSuccess! Created {output_sdf} containing exactly 20 small molecules.")
