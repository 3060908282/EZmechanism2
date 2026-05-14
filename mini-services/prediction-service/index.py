"""
M-CSA Enzyme Mechanism Prediction Service

A Flask-based microservice that loads reaction rules from an Excel file,
matches substrate SMILES against rule reactant SMARTS patterns, and
predicts enzyme reaction mechanism products using RDKit.
"""

import os
import sys
import json
import time
import re
import logging
import traceback
from typing import List, Dict, Optional, Any, Tuple

from flask import Flask, request, jsonify
from flask_cors import CORS

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Disable RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS, rdMolDescriptors, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

# Shared functions (factored out from index.py)
from shared import (
    simplify_smarts,
    get_builtin_rules,
    load_rules_from_db,
    load_rules_from_xlsx,
    RULES_FILE,
    RULES_FILE_TEST,
)

app = Flask(__name__)
CORS(app)

# ---- Global state (RULES_FILE / RULES_FILE_TEST imported from shared.py) ----
rules_cache: List[Dict[str, Any]] = []
rules_cache_2: List[Dict[str, Any]] = []
builtin_rules: List[Dict[str, Any]] = []
rules_loaded = False
START_TIME = time.time()
LAST_LOAD_TIME: Optional[float] = None


# (get_builtin_rules, load_rules_from_xlsx, simplify_smarts now imported from shared.py)


def load_rules() -> List[Dict[str, Any]]:
    """Load reaction rules from Excel files and combine with built-in rules.
    
    Rules are cached in memory after first load for fast subsequent access.
    """
    global rules_cache, rules_cache_2, builtin_rules, rules_loaded, LAST_LOAD_TIME

    if rules_loaded:
        return rules_cache + rules_cache_2 + builtin_rules

    t0 = time.time()

    # Try SQLite first (much faster: 0.3s vs 12s for 50k rules)
    try:
        rules_cache = load_rules_from_db(max_rules=10000)
        logger.info(f"  Loaded {len(rules_cache)} rules from SQLite")
    except Exception as e:
        logger.warning(f"  SQLite load failed: {e}, falling back to xlsx")
        rules_cache = []

    # Fallback to xlsx
    if not rules_cache:
        try:
            rules_cache = load_rules_from_xlsx(RULES_FILE, max_rules=10000, source_tag='mcsa')
        except Exception as e:
            logger.error(f"Error loading primary rules file: {e}")
            rules_cache = []

    # Load supplementary test rules file
    try:
        if os.path.exists(RULES_FILE_TEST):
            existing_ids = {r['rule_id'] for r in rules_cache if r.get('rule_id')}
            test_rules = load_rules_from_xlsx(RULES_FILE_TEST, source_tag='mcsa-test')
            rules_cache_2 = [r for r in test_rules if r.get('rule_id') not in existing_ids]
            logger.info(f"  After dedup: {len(rules_cache_2)} unique rules from test file")
        else:
            rules_cache_2 = []
    except Exception as e:
        logger.error(f"Error loading supplementary rules file: {e}")
        rules_cache_2 = []

    builtin_rules = get_builtin_rules()
    rules_loaded = True
    LAST_LOAD_TIME = time.time()

    total = len(rules_cache) + len(rules_cache_2) + len(builtin_rules)
    logger.info(f"Total rules: {len(rules_cache)} (mcsa) + {len(rules_cache_2)} (mcsa-test) + {len(builtin_rules)} (builtin) = {total} in {time.time() - t0:.2f}s")
    return rules_cache + rules_cache_2 + builtin_rules


def match_smiles_against_smarts(smiles: str, reactant_smarts: str) -> Tuple[bool, float]:
    """Check if a substrate SMILES matches a reactant SMARTS pattern.
    
    Returns (matched: bool, score: float) where score indicates match quality.
    
    Optimized for M-CSA rules: uses aggressive simplification + component-level matching.
    M-CSA rules describe multi-molecule enzymatic reactions (substrate + cofactors + residues),
    so we match against individual components rather than the full pattern.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (False, 0.0)

        # 1. Try simplified individual components (best strategy for M-CSA rules)
        simplified = simplify_smarts(reactant_smarts)
        for comp in simplified.split('.'):
            comp = comp.strip()
            if not comp or len(comp) < 3 or len(comp) > 300:
                continue
            try:
                pat = Chem.MolFromSmarts(comp)
                if pat is not None and mol.HasSubstructMatch(pat):
                    return (True, 0.8)
            except:
                pass

        # 2. Try original SMARTS components
        for comp in reactant_smarts.split('.'):
            comp = comp.strip()
            if not comp or len(comp) > 300:
                continue
            try:
                pat = Chem.MolFromSmarts(comp)
                if pat is not None and mol.HasSubstructMatch(pat):
                    return (True, 1.0)
            except:
                pass

        # 3. Try full simplified pattern
        if simplified != reactant_smarts:
            try:
                simp_pattern = Chem.MolFromSmarts(simplified)
                if simp_pattern is not None and mol.HasSubstructMatch(simp_pattern):
                    return (True, 0.6)
            except:
                pass

        # 4. Try full original pattern
        try:
            pattern = Chem.MolFromSmarts(reactant_smarts)
            if pattern is not None and mol.HasSubstructMatch(pattern):
                return (True, 1.0)
        except:
            pass

        return (False, 0.0)

    except Exception:
        return (False, 0.0)


def apply_rule(smiles: str, reaction_smarts: str) -> List[str]:
    """Apply a reaction SMARTS rule to a substrate and return predicted product SMILES."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return []

        # Check how many reactants the rule expects
        reactant_part = reaction_smarts.split('>>')[0]
        num_reactants = reactant_part.count('.') + 1

        results = []

        try:
            rxn = AllChem.ReactionFromSmarts(reaction_smarts)
            if rxn is None:
                return []
        except:
            return []

        if num_reactants == 1:
            results = rxn.RunReactants([mol], maxProducts=10)
        else:
            # Multi-reactant rule
            components = [c.strip() for c in reactant_part.split('.')]

            # Find which component the substrate matches best
            substrate_idx = -1
            best_atoms = 0
            for i, comp in enumerate(components):
                try:
                    pat = Chem.MolFromSmarts(comp)
                    if pat is not None:
                        match = mol.HasSubstructMatch(pat)
                        if match and pat.GetNumAtoms() > best_atoms:
                            best_atoms = pat.GetNumAtoms()
                            substrate_idx = i
                except:
                    pass

            # Build reactant list
            reactants = []
            for i in range(num_reactants):
                if i == substrate_idx or substrate_idx < 0:
                    reactants.append(mol)
                else:
                    comp_smarts = components[i]
                    gen_mol = None
                    try:
                        gen_mol = Chem.MolFromSmiles(comp_smarts)
                    except:
                        pass
                    if gen_mol is None:
                        try:
                            gen_mol = Chem.MolFromSmarts(comp_smarts)
                        except:
                            pass
                    reactants.append(gen_mol if gen_mol else mol)

            if len(reactants) == num_reactants:
                try:
                    results = rxn.RunReactants(reactants, maxProducts=10)
                except:
                    results = []

        if not results:
            return []

        product_smiles_list = []
        seen = set()

        for product_tuple in results:
            for p in product_tuple:
                try:
                    Chem.SanitizeMol(p)
                    smi = Chem.MolToSmiles(p)
                    if smi not in seen:
                        seen.add(smi)
                        product_smiles_list.append(smi)
                except:
                    continue

        return product_smiles_list

    except Exception as e:
        logger.debug(f"Error applying rule: {e}")
        return []


@app.route('/api/health', methods=['GET'])
def health():
    """Enhanced health check endpoint with uptime, rule counts, and timestamps."""
    rules = load_rules()
    uptime_seconds = round(time.time() - START_TIME, 1)
    uptime_str = f"{int(uptime_seconds // 3600)}h {int((uptime_seconds % 3600) // 60)}m {int(uptime_seconds % 60)}s"

    result = {
        'status': 'ok',
        'service': 'M-CSA Mechanism Predictor',
        'uptime_seconds': uptime_seconds,
        'uptime': uptime_str,
        'rules': {
            'mcsa_rules_natmet': len(rules_cache),
            'mcsa_rules_test': len(rules_cache_2),
            'builtin_rules': len(builtin_rules),
            'total_rules': len(rules),
        },
        'last_load_timestamp': LAST_LOAD_TIME,
        'last_load_iso': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(LAST_LOAD_TIME)) if LAST_LOAD_TIME else None,
        'start_time': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(START_TIME)),
    }
    return jsonify(result)


@app.route('/api/match-rules', methods=['POST'])
def match_rules():
    """Match substrate SMILES against reaction rules.
    
    Body: { "smiles": "CC(=O)O" }
    Response: { "matches": [...], "total_rules": N, "smiles": "..." }
    """
    rules = load_rules()

    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()

        if not smiles:
            return jsonify({'error': 'SMILES string is required'}), 400

        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': f'Invalid SMILES: {smiles}'}), 400

        t0 = time.time()
        matches = []

        for rule in rules:
            reactant_smarts = rule.get('reactant_smarts', '')
            if not reactant_smarts:
                # For built-in rules, parse from reaction_smarts
                rxn = rule.get('reaction_smarts', '')
                if '>>' in rxn:
                    reactant_smarts = rxn.split('>>')[0]
                else:
                    continue

            matched, score = match_smiles_against_smarts(smiles, reactant_smarts)
            if matched:
                # Try to apply the rule
                products = apply_rule(smiles, rule['reaction_smarts'])

                match_entry = {
                    'rule_id': rule['rule_id'],
                    'rule_name': rule.get('rule_name', f'Rule {rule["rule_id"]}'),
                    'mcsa_id': rule.get('mcsa_id'),
                    'step_id': rule.get('step_id'),
                    'mechanism_id': rule.get('mechanism_id'),
                    'reaction_smarts': rule['reaction_smarts'],
                    'reactant_smarts': reactant_smarts,
                    'product_smarts': rule.get('product_smarts', ''),
                    'products': products,
                    'match_score': round(score, 3),
                    'category': rule.get('category', 'Unknown'),
                    'enzyme': rule.get('enzyme', ''),
                    'source': rule.get('source', 'mcsa'),
                }
                matches.append(match_entry)

        # Sort by match score
        matches.sort(key=lambda x: x['match_score'], reverse=True)

        elapsed = time.time() - t0
        logger.info(f"Matched {smiles}: {len(matches)}/{len(rules)} rules in {elapsed:.3f}s")

        return jsonify({
            'smiles': smiles,
            'matches': matches,
            'total_rules': len(rules),
            'matches_count': len(matches),
            'elapsed_seconds': round(elapsed, 3)
        })

    except Exception as e:
        logger.error(f"Error in match-rules: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/predict', methods=['POST'])
def predict():
    """Predict multi-step reaction mechanism.
    
    Body: { "smiles": "CC(=O)O", "max_steps": 5 }
    Response: { "steps": [...], "total_time": ..., "rules_checked": ... }
    """
    rules = load_rules()

    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()
        max_steps = min(data.get('max_steps', 5), 10)

        if not smiles:
            return jsonify({'error': 'SMILES string is required'}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': f'Invalid SMILES: {smiles}'}), 400

        t_start = time.time()
        steps = []
        current_smiles_list = [smiles]
        total_rules_checked = 0
        visited = {smiles}

        for step_num in range(1, max_steps + 1):
            all_products = []
            step_rules = []
            seen_products = set()

            for substrate in current_smiles_list:
                for rule in rules:
                    total_rules_checked += 1
                    reactant_smarts = rule.get('reactant_smarts', '')
                    if not reactant_smarts:
                        rxn = rule.get('reaction_smarts', '')
                        if '>>' in rxn:
                            reactant_smarts = rxn.split('>>')[0]
                        else:
                            continue

                    matched, score = match_smiles_against_smarts(substrate, reactant_smarts)
                    if matched and score >= 0.5:
                        products = apply_rule(substrate, rule['reaction_smarts'])
                        for p in products:
                            if p not in seen_products and p not in visited:
                                seen_products.add(p)
                                all_products.append(p)
                                step_rules.append({
                                    'rule_id': rule['rule_id'],
                                    'rule_name': rule.get('rule_name', f'Rule {rule["rule_id"]}'),
                                    'mcsa_id': rule.get('mcsa_id'),
                                    'reaction_smarts': rule['reaction_smarts'],
                                    'substrate': substrate,
                                    'product': p,
                                    'match_score': round(score, 3),
                                    'category': rule.get('category', ''),
                                    'enzyme': rule.get('enzyme', ''),
                                })

            if not all_products:
                break

            steps.append({
                'step': step_num,
                'substrates': list(current_smiles_list),
                'products': all_products,
                'rules_applied': step_rules,
                'rules_count': len(step_rules),
            })

            for p in all_products:
                visited.add(p)
            current_smiles_list = all_products

        elapsed = time.time() - t_start
        logger.info(f"Prediction for {smiles}: {len(steps)} steps, {total_rules_checked} rules checked in {elapsed:.3f}s")

        return jsonify({
            'smiles': smiles,
            'steps': steps,
            'max_steps': max_steps,
            'total_steps': len(steps),
            'total_rules_checked': total_rules_checked,
            'total_rules_available': len(rules),
            'elapsed_seconds': round(elapsed, 3),
            'unique_products': len(visited) - 1
        })

    except Exception as e:
        logger.error(f"Error in predict: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/mol-svg', methods=['POST'])
def mol_svg():
    """Generate a 2D molecule SVG from a SMILES string.

    Body: { "smiles": "CCO" }
    Query: ?numbered=true  (include atom numbering)
    Response: { "smiles": "...", "svg": "...", "width": N, "height": N }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()
        numbered = request.args.get('numbered', 'false').lower() == 'true'

        if not smiles:
            return jsonify({'error': 'SMILES string is required'}), 400

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': f'Invalid SMILES: {smiles}'}), 400

        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)

        img_width = 300
        img_height = 300

        drawer = rdMolDraw2D.MolDraw2DSVG(img_width, img_height)

        if numbered:
            opts = drawer.drawOptions()
            opts.addAtomIndices = True

        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg_text = drawer.GetDrawingText()

        return jsonify({
            'smiles': smiles,
            'canonical_smiles': Chem.MolToSmiles(mol),
            'svg': svg_text,
            'width': img_width,
            'height': img_height,
            'numbered': numbered,
        })

    except Exception as e:
        logger.error(f"Error in mol-svg: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/reaction-svg', methods=['POST'])
def reaction_svg():
    """Render a reaction as a visual SVG diagram (reactant → product).

    Body: { "reaction_smarts": "CCO>>CC=O", "smiles": "CCO" }
    Response: { "svg": "...", "reactants": [...], "products": [...] }
    """
    try:
        data = request.get_json()
        reaction_smarts = data.get('reaction_smarts', '').strip()
        smiles = data.get('smiles', '').strip()

        if not reaction_smarts:
            return jsonify({'error': 'reaction_smarts is required'}), 400

        if '>>' not in reaction_smarts:
            return jsonify({'error': 'reaction_smarts must contain >> to separate reactants and products'}), 400

        parts = reaction_smarts.split('>>', 1)
        reactant_smarts = parts[0].strip()
        product_smarts = parts[1].strip()

        # Try to build molecules from SMARTS/SMILES
        def parse_mol_part(s):
            """Parse a dot-separated SMARTS/SMILES string into a list of RDKit mol objects."""
            mols = []
            for comp in s.split('.'):
                comp = comp.strip()
                if not comp:
                    continue
                mol = None
                # Try SMILES first (more likely to produce valid molecules)
                try:
                    mol = Chem.MolFromSmiles(comp)
                except:
                    pass
                if mol is None:
                    try:
                        mol = Chem.MolFromSmarts(comp)
                    except:
                        pass
                if mol is not None:
                    AllChem.Compute2DCoords(mol)
                    mols.append(mol)
            return mols

        reactant_mols = parse_mol_part(reactant_smarts)
        product_mols = parse_mol_part(product_smarts)

        # If SMARTS didn't produce molecules but we have a substrate SMILES, use it as reactant
        if not reactant_mols and smiles:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                AllChem.Compute2DCoords(mol)
                reactant_mols = [mol]

        if not reactant_mols and not product_mols:
            return jsonify({'error': 'Could not parse any valid molecules from the reaction SMARTS'}), 400

        # Generate SVG for each side
        def mols_to_svg(mols, label):
            if not mols:
                return ''
            # Combine mols into a single molecule for drawing
            combined = mols[0]
            for m in mols[1:]:
                combined = Chem.CombineMols(combined, m)

            w, h = 400, 200
            drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
            drawer.DrawMolecule(combined)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            return svg

        reactant_svgs = mols_to_svg(reactant_mols, 'Reactants')
        product_svgs = mols_to_svg(product_mols, 'Products')

        # Compose the reaction diagram SVG
        # Calculate total width based on components
        total_width = 900
        arrow_width = 80
        reactant_width = 350 if reactant_mols else 0
        product_width = 350 if product_mols else 0
        total_width = max(total_width, reactant_width + arrow_width + product_width + 60)
        height = 250

        # Build composed SVG
        composed_svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="{total_width}" height="{height}" viewBox="0 0 {total_width} {height}">
  <style>.mol-img {{ image-rendering: optimizeQuality; }}</style>
  <rect width="100%" height="100%" fill="white" rx="8"/>
  <text x="15" y="20" font-family="Arial, sans-serif" font-size="12" fill="#666">Reaction Diagram</text>
'''

        x_offset = 20
        if reactant_mols:
            # Extract the inner SVG content
            inner = reactant_svgs
            # Replace the outer svg tag
            inner = re.sub(r'<svg[^>]*>', '', inner)
            inner = inner.replace('</svg>', '')
            composed_svg += f'  <g transform="translate({x_offset}, 25)" class="mol-img">{inner}</g>\n'
            composed_svg += f'  <text x="{x_offset + 10}" y="{height - 10}" font-family="Arial, sans-serif" font-size="11" fill="#888">Reactant(s)</text>\n'
            x_offset += reactant_width

        # Arrow
        arrow_x = x_offset + 10
        arrow_end = arrow_x + arrow_width
        composed_svg += f'  <line x1="{arrow_x}" y1="{height // 2}" x2="{arrow_end}" y2="{height // 2}" stroke="#333" stroke-width="2" marker-end="url(#arrowhead)"/>\n'
        composed_svg += '  <defs><marker id="arrowhead" markerWidth="10" markerHeight="7" refX="10" refY="3.5" orient="auto"><polygon points="0 0, 10 3.5, 0 7" fill="#333"/></marker></defs>\n'
        x_offset = arrow_end + 10

        if product_mols:
            inner = product_svgs
            inner = re.sub(r'<svg[^>]*>', '', inner)
            inner = inner.replace('</svg>', '')
            composed_svg += f'  <g transform="translate({x_offset}, 25)" class="mol-img">{inner}</g>\n'
            composed_svg += f'  <text x="{x_offset + 10}" y="{height - 10}" font-family="Arial, sans-serif" font-size="11" fill="#888">Product(s)</text>\n'

        composed_svg += '</svg>'

        return jsonify({
            'svg': composed_svg,
            'reactants': [Chem.MolToSmiles(m) if m is not None else None for m in reactant_mols],
            'products': [Chem.MolToSmiles(m) if m is not None else None for m in product_mols],
            'reaction_smarts': reaction_smarts,
        })

    except Exception as e:
        logger.error(f"Error in reaction-svg: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/batch-predict', methods=['POST'])
def batch_predict():
    """Run prediction for multiple SMILES in a single request.

    Body: { "smiles_list": ["CCO", "CC(=O)O"], "max_steps": 3 }
    Response: { "results": [...], "total_time": N }
    """
    rules = load_rules()

    try:
        data = request.get_json()
        smiles_list = data.get('smiles_list', [])
        max_steps = min(data.get('max_steps', 3), 10)

        if not isinstance(smiles_list, list) or len(smiles_list) == 0:
            return jsonify({'error': 'smiles_list must be a non-empty array'}), 400

        if len(smiles_list) > 50:
            return jsonify({'error': 'smiles_list can contain at most 50 entries'}), 400

        t_start = time.time()
        results = []

        for idx, smiles in enumerate(smiles_list):
            smiles = smiles.strip() if isinstance(smiles, str) else ''
            entry = {
                'index': idx,
                'smiles': smiles,
                'status': 'ok',
                'steps': [],
                'total_steps': 0,
                'unique_products': 0,
                'elapsed_seconds': 0,
            }

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    entry['status'] = 'invalid_smiles'
                    entry['error'] = f'Invalid SMILES: {smiles}'
                    results.append(entry)
                    continue

                entry_t0 = time.time()
                steps = []
                current_smiles_list = [smiles]
                visited = {smiles}

                for step_num in range(1, max_steps + 1):
                    all_products = []
                    step_rules = []
                    seen_products = set()

                    for substrate in current_smiles_list:
                        for rule in rules:
                            reactant_smarts = rule.get('reactant_smarts', '')
                            if not reactant_smarts:
                                rxn = rule.get('reaction_smarts', '')
                                if '>>' in rxn:
                                    reactant_smarts = rxn.split('>>')[0]
                                else:
                                    continue

                            matched, score = match_smiles_against_smarts(substrate, reactant_smarts)
                            if matched and score >= 0.5:
                                products = apply_rule(substrate, rule['reaction_smarts'])
                                for p in products:
                                    if p not in seen_products and p not in visited:
                                        seen_products.add(p)
                                        all_products.append(p)
                                        step_rules.append({
                                            'rule_id': rule['rule_id'],
                                            'rule_name': rule.get('rule_name', f'Rule {rule["rule_id"]}'),
                                            'substrate': substrate,
                                            'product': p,
                                            'match_score': round(score, 3),
                                        })

                    if not all_products:
                        break

                    steps.append({
                        'step': step_num,
                        'products': all_products,
                        'rules_applied': len(step_rules),
                    })

                    for p in all_products:
                        visited.add(p)
                    current_smiles_list = all_products

                entry['steps'] = steps
                entry['total_steps'] = len(steps)
                entry['unique_products'] = len(visited) - 1
                entry['elapsed_seconds'] = round(time.time() - entry_t0, 3)

            except Exception as ex:
                entry['status'] = 'error'
                entry['error'] = str(ex)

            results.append(entry)

        total_time = round(time.time() - t_start, 3)
        logger.info(f"Batch predict: {len(smiles_list)} inputs, {total_time:.3f}s total")

        return jsonify({
            'results': results,
            'total_count': len(results),
            'max_steps': max_steps,
            'total_rules_available': len(rules),
            'total_time': total_time,
        })

    except Exception as e:
        logger.error(f"Error in batch-predict: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


@app.route('/api/validate', methods=['POST'])
def validate_smiles():
    """Validate a SMILES string."""
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()

        if not smiles:
            return jsonify({'valid': False, 'error': 'Empty SMILES'})

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'valid': False, 'error': f'Invalid SMILES: {smiles}'})

        formula = rdMolDescriptors.CalcMolFormula(mol)
        canonical = Chem.MolToSmiles(mol)
        mw = Descriptors.MolWt(mol)

        return jsonify({
            'valid': True,
            'smiles': smiles,
            'canonical_smiles': canonical,
            'formula': formula,
            'num_atoms': mol.GetNumAtoms(),
            'num_heavy_atoms': mol.GetNumHeavyAtoms(),
            'molecular_weight': round(mw, 4)
        })

    except Exception as e:
        return jsonify({'valid': False, 'error': str(e)}), 500


@app.route('/api/mol-info', methods=['POST'])
def mol_info():
    """Get detailed molecule information."""
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': f'Invalid SMILES: {smiles}'}), 400

        info = {
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

        return jsonify(info)

    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/smiles-explicit-h', methods=['POST'])
def smiles_explicit_h():
    """Convert SMILES to explicit-H SMILES representation.

    Body: { "smiles": "CCO.CCO" }
    Response: { "explicit_h_smiles": "[H]C([H])([H])C([H])([H])O[H].[H]C([H])([H])C([H])([H])O[H]" }
    """
    try:
        data = request.get_json()
        smiles = data.get('smiles', '').strip()
        if not smiles:
            return jsonify({'error': 'No SMILES provided'}), 400

        # Handle multi-molecule SMILES (dot-separated)
        parts = [s.strip() for s in smiles.split('.') if s.strip()]
        explicit_parts = []
        for part in parts:
            mol = Chem.MolFromSmiles(part)
            if mol is None:
                explicit_parts.append(f'[INVALID:{part}]')
                continue
            mol_h = Chem.AddHs(mol)
            smi_h = Chem.MolToSmiles(mol_h, isomericSmiles=True)
            explicit_parts.append(smi_h)

        return jsonify({
            'smiles': smiles,
            'explicit_h_smiles': '.'.join(explicit_parts),
            'parts': [{'original': p, 'explicit_h': ep} for p, ep in zip(parts, explicit_parts)],
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/api/mechanism-search', methods=['POST'])
def mechanism_search():
    """Bidirectional BFS + Dijkstra mechanism search (proxied from mechanism_search.py).

    Body: {
        "reactants": "CCO",
        "products": "CC(=O)O",
        "max_configs": 300,
        "max_rules": 0,        // 0 = load all rules (default, uses lazy parsing)
        "max_bond_length": 6.0
    }
    Response: { "paths": [...], "graph": {...}, "stats": {...} }
    """
    import subprocess
    import tempfile

    try:
        data = request.get_json()
        reactants = data.get('reactants', '').strip()
        products = data.get('products', '').strip()  # Can be empty for forward-only mode

        if not reactants:
            return jsonify({'error': 'Reactants SMILES are required'}), 400

        # Validate SMILES
        if Chem.MolFromSmiles(reactants) is None:
            return jsonify({'error': f'Invalid reactants SMILES: {reactants}'}), 400
        if products and Chem.MolFromSmiles(products) is None:
            return jsonify({'error': f'Invalid products SMILES: {products}'}), 400

        max_configs = min(max(int(data.get('max_configs', 300)), 50), 2000)
        # max_rules controls how many M-CSA rules to load from SQLite/xlsx.
        # Default 0 = load all 51K rules (uses batch processing to avoid OOM).
        # Set to e.g. 10000 for faster search with smaller rule set.
        max_rules = int(data.get('max_rules', 0))
        if max_rules < 0:
            max_rules = 0  # 0 = all
        max_bond_length = min(max(float(data.get('max_bond_length', 6.0)), 1.0), 20.0)
        mcsa_id = int(data.get('mcsa_id', 0))
        if mcsa_id < 0:
            mcsa_id = 0  # 0 = use all rules

        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'mechanism_search.py')

        cmd = [
            sys.executable, script_path,
            '--reactants', reactants,
            '--products', products,
            '--max-configs', str(max_configs),
            '--max-rules', str(max_rules),
            '--max-bond-length', str(max_bond_length),
        ]

        if mcsa_id > 0:
            cmd.extend(['--mcsa-id', str(mcsa_id)])

        # Add optional arguments
        residues_data = data.get('residues')
        if residues_data and isinstance(residues_data, list):
            cmd.extend(['--residues', json.dumps(residues_data)])

        pdb_id_data = data.get('pdb_id')
        if pdb_id_data:
            cmd.extend(['--pdb-id', str(pdb_id_data)])

        pdb_text_data = data.get('pdb_text')
        pdb_temp_file = None
        if pdb_text_data and isinstance(pdb_text_data, str) and len(pdb_text_data) > 100:
            # Write PDB text to temp file to avoid E2BIG (argument list too long)
            pdb_temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
            pdb_temp_file.write(pdb_text_data)
            pdb_temp_file.close()
            cmd.extend(['--pdb-text-file', pdb_temp_file.name])

        logger.info(f"Running mechanism search: {reactants} -> {products}")
        t0 = time.time()

        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                cwd=os.path.dirname(os.path.abspath(__file__)),
            )

            elapsed = time.time() - t0

            # Log subprocess output for debugging matching details
            if result.stdout:
                for line in result.stdout.strip().split('\n')[:50]:
                    logger.info(f"[search-stdout] {line}")
            if result.stderr:
                for line in result.stderr.strip().split('\n')[:30]:
                    logger.warning(f"[search-stderr] {line}")

            if result.returncode != 0:
                error_msg = result.stderr.strip() if result.stderr else 'Unknown error'
                logger.error(f"Mechanism search failed: {error_msg}")
                return jsonify({'error': error_msg}), 500

            try:
                data = json.loads(result.stdout.strip())
                data['flask_elapsed'] = round(elapsed, 3)
                return jsonify(data)
            except json.JSONDecodeError as e:
                logger.error(f"Failed to parse mechanism search output: {e}")
                return jsonify({'error': f'Invalid output from search script: {e}'}), 500
        finally:
            # Clean up temp PDB file
            if pdb_temp_file:
                try:
                    os.unlink(pdb_temp_file.name)
                except Exception:
                    pass

    except subprocess.TimeoutExpired:
        logger.error("Mechanism search timed out after 600s")
        return jsonify({'error': 'Mechanism search timed out after 600 seconds'}), 504
    except Exception as e:
        logger.error(f"Error in mechanism-search: {e}")
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500


# ---- PDB Handler Endpoints (no RDKit required) ----
from pdb_handler import fetch_pdb_info, parse_pdb_text, get_active_site_prediction
import urllib.request

@app.route('/api/pdb-fetch', methods=['POST'])
def pdb_fetch():
    """Fetch PDB structure info from RCSB by PDB ID.
    Body: { "pdb_id": "1EJG" }
    """
    try:
        data = request.get_json()
        pdb_id = data.get('pdb_id', '').strip().upper()
        if not pdb_id or len(pdb_id) != 4:
            return jsonify({'error': 'Valid 4-character PDB ID is required'}), 400
        logger.info(f"PDB fetch: {pdb_id}")
        info = fetch_pdb_info(pdb_id)
        return jsonify(info)
    except Exception as e:
        logger.error(f"PDB fetch error: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/pdb-parse', methods=['POST'])
def pdb_parse():
    """Parse PDB format text.
    Body: { "pdb_text": "HEADER ..." }
    """
    try:
        data = request.get_json()
        pdb_text = data.get('pdb_text', '').strip()
        if not pdb_text or len(pdb_text) < 100:
            return jsonify({'error': 'Valid PDB text required (min 100 chars)'}), 400
        info = parse_pdb_text(pdb_text)
        return jsonify(info)
    except Exception as e:
        logger.error(f"PDB parse error: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/pdb-raw-fetch', methods=['POST'])
def pdb_raw_fetch():
    """Fetch raw PDB text from RCSB (for 3D viewer).
    Body: { "pdb_id": "1EJG" }
    """
    try:
        data = request.get_json()
        pdb_id = data.get('pdb_id', '').strip().upper()
        if not pdb_id or len(pdb_id) != 4:
            return jsonify({'error': 'Valid 4-character PDB ID is required'}), 400
        url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
        req = urllib.request.Request(url, headers={"User-Agent": "EZmechanism/1.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            if resp.status != 200:
                return jsonify({'error': f'RCSB returned status {resp.status}'}), 500
            pdb_text = resp.read().decode('utf-8', errors='replace')
        return jsonify({'pdb_text': pdb_text})
    except Exception as e:
        logger.error(f"PDB raw fetch error: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/pdb-active-site', methods=['POST'])
def pdb_active_site():
    """Predict catalytic active site residues.
    Body: { "pdb_id": "1EJG" }
    """
    try:
        data = request.get_json()
        pdb_id = data.get('pdb_id', '').strip().upper()
        if not pdb_id or len(pdb_id) != 4:
            return jsonify({'error': 'Valid 4-character PDB ID is required'}), 400
        logger.info(f"PDB active-site: {pdb_id}")
        info = fetch_pdb_info(pdb_id)
        result = get_active_site_prediction(info)
        return jsonify(result)
    except Exception as e:
        logger.error(f"PDB active-site error: {e}")
        return jsonify({'error': str(e)}), 500





if __name__ == '__main__':
    load_rules()
    logger.info("Starting M-CSA Mechanism Prediction Service on port 3003...")
    app.run(host='0.0.0.0', port=3003, debug=False, threaded=True)
