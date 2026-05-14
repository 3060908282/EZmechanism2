#!/usr/bin/env python3
"""Get molecular information for a SMILES string.

Usage:
    python3 mol_info.py <smiles> [--explicit-h]

Output: JSON to stdout

When --explicit-h is passed, returns explicit-h SMILES representation
for each molecule component in a dot-separated SMILES string.
"""

import json
import sys

from rdkit import RDLogger, Chem
from rdkit.Chem import rdMolDescriptors, Descriptors
RDLogger.DisableLog('rdApp.*')

from shared import validate_smiles


def explicit_h_smiles(smiles: str) -> dict:
    """Convert SMILES to explicit-H SMILES representation."""
    parts = [s.strip() for s in smiles.split('.') if s.strip()]
    explicit_parts = []
    part_details = []
    for part in parts:
        mol = Chem.MolFromSmiles(part)
        if mol is None:
            explicit_parts.append(f'[INVALID:{part}]')
            part_details.append({'original': part, 'explicit_h': f'[INVALID:{part}]'})
            continue
        mol_h = Chem.AddHs(mol)
        smi_h = Chem.MolToSmiles(mol_h, isomericSmiles=True)
        explicit_parts.append(smi_h)
        part_details.append({'original': part, 'explicit_h': smi_h})

    return {
        'smiles': smiles,
        'explicit_h_smiles': '.'.join(explicit_parts),
        'parts': part_details,
    }


def main() -> None:
    if len(sys.argv) < 2:
        print(json.dumps({'error': 'Usage: python3 mol_info.py <smiles> [--explicit-h]'}))
        sys.exit(1)

    smiles = sys.argv[1]
    explicit_h = '--explicit-h' in sys.argv

    validate_smiles(smiles)

    if explicit_h:
        result = explicit_h_smiles(smiles)
        print(json.dumps(result))
        return

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(json.dumps({'valid': False, 'error': 'Invalid SMILES'}))
        sys.exit(1)

    info = {
        'valid': True,
        'smiles': smiles,
        'canonical_smiles': Chem.MolToSmiles(mol),
        'formula': rdMolDescriptors.CalcMolFormula(mol),
        'molecular_weight': round(Descriptors.MolWt(mol), 4),
        'num_atoms': mol.GetNumAtoms(),
        'num_heavy_atoms': mol.GetNumHeavyAtoms(),
        'num_rotatable_bonds': Descriptors.NumRotatableBonds(mol),
        'num_h_bond_donors': Descriptors.NumHDonors(mol),
        'num_h_bond_acceptors': Descriptors.NumHAcceptors(mol),
        'tpsa': round(Descriptors.TPSA(mol), 2),
        'logp': round(Descriptors.MolLogP(mol), 4),
        'ring_count': Descriptors.RingCount(mol),
        'aromatic_rings': Descriptors.NumAromaticRings(mol),
    }
    print(json.dumps(info))


if __name__ == '__main__':
    main()
