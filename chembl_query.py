from chembl_webresource_client.new_client import new_client
from rdkit import Chem
import pandas as pd

def run_query(smarts, targets):
    query_mol = Chem.MolFromSmarts(smarts)
    activity = new_client.activity

    results = []

    for target_name, target_id in targets.items():
        acts = activity.filter(
            target_chembl_id=target_id,
            standard_type__in=["IC50", "EC50"],
            assay_type="B"
        )

        for a in acts:
            smiles = a.get("canonical_smiles")
            if not smiles:
                continue

            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.HasSubstructMatch(query_mol):
                results.append({
                    "target": target_name,
                    "molecule_chembl_id": a["molecule_chembl_id"],
                    "smiles": smiles,
                    "standard_type": a["standard_type"],
                    "standard_value": a["standard_value"],
                    "standard_units": a["standard_units"],
                })

    return pd.DataFrame(results)