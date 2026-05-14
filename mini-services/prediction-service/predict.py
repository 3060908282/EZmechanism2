#!/usr/bin/env python3
"""Predict reaction products from SMILES using M-CSA rules.

Usage:
    python3 predict.py <smiles> [max_steps]

    max_steps: integer 1-10 (defaults to 3)

Output: JSON to stdout
"""

import json
import sys
import time

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

from shared import validate_smiles, simplify_smarts, load_rules_from_xlsx, RULES_FILE


def main() -> None:
    if len(sys.argv) < 2:
        print(json.dumps({'error': 'Usage: python3 predict.py <smiles> [max_steps]'}))
        sys.exit(1)

    smiles = sys.argv[1]
    validate_smiles(smiles)

    max_steps = 3
    if len(sys.argv) >= 3:
        try:
            max_steps = int(sys.argv[2])
            if max_steps < 1 or max_steps > 10:
                print(json.dumps({'error': 'max_steps must be an integer between 1 and 10'}))
                sys.exit(1)
        except ValueError:
            print(json.dumps({'error': 'max_steps must be an integer between 1 and 10'}))
            sys.exit(1)

    # Load rules
    rules = load_rules_from_xlsx(RULES_FILE)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(json.dumps({'error': 'Invalid SMILES'}))
        sys.exit(1)

    t_start = time.time()
    steps = []
    current = [smiles]
    visited = {smiles}

    for step_num in range(1, max_steps + 1):
        products = []
        step_rules = []
        seen = set()
        for sub in current:
            for rule in rules:
                matched = False
                score = 0.0
                simplified = simplify_smarts(rule['reactant_smarts'])
                for comp in simplified.split('.'):
                    comp = comp.strip()
                    if not comp or len(comp) < 3 or len(comp) > 300:
                        continue
                    try:
                        pat = Chem.MolFromSmarts(comp)
                        if pat and mol.HasSubstructMatch(pat):
                            matched, score = True, 0.8
                            break
                    except Exception:
                        pass
                if matched and score >= 0.5:
                    try:
                        rxn = AllChem.ReactionFromSmarts(rule['reaction_smarts'])
                        if rxn:
                            results = rxn.RunReactants([mol], maxProducts=5)
                            for ptuple in results:
                                for p in ptuple:
                                    try:
                                        Chem.SanitizeMol(p)
                                        smi = Chem.MolToSmiles(p)
                                        if smi not in seen and smi not in visited:
                                            seen.add(smi)
                                            products.append(smi)
                                            step_rules.append({
                                                'rule_id': rule['rule_id'],
                                                'mcsa_id': rule.get('mcsa_id'),
                                                'match_score': score,
                                                'category': 'Enzyme Catalysis',
                                                'substrate': sub,
                                                'product': smi,
                                            })
                                    except Exception:
                                        pass
                    except Exception:
                        pass
        if not products:
            break
        steps.append({
            'step': step_num,
            'substrates': list(current),
            'products': products,
            'rules_applied': step_rules,
            'rules_count': len(step_rules),
        })
        for p in products:
            visited.add(p)
        current = products

    elapsed = time.time() - t_start
    result = {
        'smiles': smiles,
        'steps': steps,
        'max_steps': max_steps,
        'total_steps': len(steps),
        'total_rules_available': len(rules),
        'elapsed_seconds': round(elapsed, 3),
        'unique_products': len(visited) - 1,
    }
    print(json.dumps(result))


if __name__ == '__main__':
    main()
