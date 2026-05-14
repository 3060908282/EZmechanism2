#!/usr/bin/env python3
"""Match SMILES against M-CSA reaction rules.

Usage:
    python3 match_rules.py <smiles>

Output: JSON to stdout
"""

import json
import sys
import time

from rdkit import RDLogger, Chem
RDLogger.DisableLog('rdApp.*')

from shared import validate_smiles, simplify_smarts, load_rules_from_xlsx, RULES_FILE


def main() -> None:
    if len(sys.argv) < 2:
        print(json.dumps({'error': 'Usage: python3 match_rules.py <smiles>'}))
        sys.exit(1)

    smiles = sys.argv[1]
    validate_smiles(smiles)

    # Load rules
    rules = load_rules_from_xlsx(RULES_FILE)

    # Match
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(json.dumps({'error': 'Invalid SMILES'}))
        sys.exit(1)

    t0 = time.time()
    matches = []
    for i, rule in enumerate(rules):
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
        if not matched:
            for comp in rule['reactant_smarts'].split('.'):
                comp = comp.strip()
                if not comp or len(comp) > 300:
                    continue
                try:
                    pat = Chem.MolFromSmarts(comp)
                    if pat and mol.HasSubstructMatch(pat):
                        matched, score = True, 1.0
                        break
                except Exception:
                    pass
        if matched:
            matches.append({
                'rule_id': rule['rule_id'],
                'rule_name': rule.get('rule_name', ''),
                'mcsa_id': rule.get('mcsa_id'),
                'mechanism_id': rule.get('mechanism_id'),
                'step_id': rule.get('step_id'),
                'reaction_smarts': rule['reaction_smarts'],
                'reactant_smarts': rule['reactant_smarts'],
                'match_score': score,
                'category': rule.get('category', ''),
                'source': rule.get('source', 'mcsa'),
                'products': [],
            })

    elapsed = time.time() - t0
    result = {
        'smiles': smiles,
        'matches': matches,
        'total_rules': len(rules),
        'matches_count': len(matches),
        'elapsed_seconds': round(elapsed, 3)
    }
    print(json.dumps(result))


if __name__ == '__main__':
    main()
