#!/usr/bin/env python3
"""
shared.py — Prediction Service Shared Module

Consolidates all duplicated code from match_rules.py, predict.py,
mol_info.py, mol_svg.py, mechanism_search.py, and index.py into
a single source of truth.

Provides:
  - SMILES_RE: validation regex
  - RULES_FILE: default rules Excel path
  - validate_smiles(): CLI-safe SMILES validation (exits on failure)
  - is_valid_smiles(): returns bool (non-exiting version)
  - simplify_smarts(): SMARTS pattern simplification
  - get_builtin_rules(): 10 built-in biochemical reaction rules (merged schema)
  - load_rules_from_xlsx(): load M-CSA rules from Excel
  - load_rules(): unified rule loading (builtin + M-CSA) with pre-parsing
"""

import json
import logging
import os
import re
import sys
from typing import Any, Dict, List, Optional, Tuple

import openpyxl

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RULES_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rules_nat_met.xlsx')
RULES_FILE_TEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rules_test.xlsx')
SMILES_RE = re.compile(r'^[A-Za-z0-9@+\-\[\]()\\/{},%=#.:;!?*]+$')


# ---------------------------------------------------------------------------
# SMILES Validation
# ---------------------------------------------------------------------------
def is_valid_smiles(smiles: str) -> bool:
    """Check whether *smiles* contains only allowed SMILES characters.

    Returns True if valid, False otherwise.  Does NOT call sys.exit().
    """
    return bool(smiles and SMILES_RE.match(smiles))


def validate_smiles(smiles: str) -> None:
    """Validate SMILES and exit with JSON error if invalid.

    Designed for CLI scripts that output JSON to stdout.
    """
    if not is_valid_smiles(smiles):
        print(json.dumps({'error': 'Invalid SMILES characters'}))
        sys.exit(1)


# ---------------------------------------------------------------------------
# SMARTS Simplification
# ---------------------------------------------------------------------------
def simplify_smarts(smarts: str) -> str:
    """Simplify a SMARTS pattern for flexible substructure matching.

    Removes atom mapping (:N), explicit zero charges (+0/-0), and isotope numbers
    that appear INSIDE bracket atom expressions (e.g. [13C] → [C]).
    Preserves ring-closure digits outside brackets and all non-zero charges
    (e.g. [O-], [N+], [NH3+]) because they carry critical chemical meaning.
    Collapses multiple dots into a single dot and strips leading/trailing dots.

    BUG FIX: The old regex ``\\d+(?=[A-Z#\\[\\(])`` destroyed ring-closure digits
    (e.g. [CH]1[CH2]…[CH]1 → [CH][CH2]…[CH]1, breaking cyclohexane patterns).
    Isotope digits only appear INSIDE ``[…]`` brackets, so we restrict removal
    to that context.
    """
    s = re.sub(r':\d+', '', smarts)           # remove :N mapping numbers
    s = re.sub(r'[+-]0(?![0-9])', '', s)      # remove explicit +0/-0 charge only
    # Remove isotope numbers ONLY inside bracket expressions [13C] → [C], [2H] → [H]
    # Pattern: open bracket, optional +/-, digits, then element symbol (uppercase+lowercase)
    s = re.sub(r'(\[)[+-]?\d+(?=[A-Z])', r'\1', s)  # isotope inside brackets
    s = re.sub(r'[Rr]\d*', '', s)             # remove ring membership
    s = re.sub(r'\.+', '.', s)                # collapse multiple dots
    return s.strip('.')


# ---------------------------------------------------------------------------
# Built-in Biochemical Reaction Rules
# ---------------------------------------------------------------------------
def get_builtin_rules() -> List[Dict[str, Any]]:
    """Return 10 built-in simple biochemical reaction rules.

    Merged schema with ALL fields needed by both mechanism_search.py and
    index.py: reactant_smarts, product_smarts, category, enzyme, ec_number,
    description, rule_name, rule_id.
    """
    return [
        {
            'rule_id': 'BUILTIN-001',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[CH2:1][CH2:2][OH:3]>>[CH:1]=[O:4].[CH3:2]',
            'reactant_smarts': '[CH2:1][CH2:2][OH:3]',
            'product_smarts': '[CH:1]=[O:4].[CH3:2]',
            'rule_name': 'Alcohol Oxidation (Ethanol → Acetaldehyde)',
            'category': 'Oxidation',
            'enzyme': 'Alcohol Dehydrogenase',
            'ec_number': 'EC 1.1.1.1',
            'description': 'Oxidation of primary alcohol to aldehyde using NAD+ as cofactor',
        },
        {
            'rule_id': 'BUILTIN-002',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[CH3:1][CH:2]=[O:3]>>[CH3:1][C:2]([OH:4])=[O:3]',
            'reactant_smarts': '[CH3:1][CH:2]=[O:3]',
            'product_smarts': '[CH3:1][C:2]([OH:4])=[O:3]',
            'rule_name': 'Aldehyde Oxidation (Acetaldehyde → Acetic Acid)',
            'category': 'Oxidation',
            'enzyme': 'Aldehyde Dehydrogenase',
            'ec_number': 'EC 1.2.1.3',
            'description': 'Oxidation of aldehyde to carboxylic acid',
        },
        {
            'rule_id': 'BUILTIN-003',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[C:1]([O:2][C:3])(=[O:4])[C:5]>>[C:1]([OH:2])(=[O:4])[C:5].[O:6][C:3]',
            'reactant_smarts': '[C:1]([O:2][C:3])(=[O:4])[C:5]',
            'product_smarts': '[C:1]([OH:2])(=[O:4])[C:5].[O:6][C:3]',
            'rule_name': 'Ester Hydrolysis',
            'category': 'Hydrolysis',
            'enzyme': 'Esterase',
            'ec_number': 'EC 3.1.1.1',
            'description': 'Hydrolysis of ester bonds',
        },
        {
            'rule_id': 'BUILTIN-004',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[OH:1][C:2]>>[O:1][P:3](=O)([O-])[O-].[C:2]',
            'reactant_smarts': '[OH:1][C:2]',
            'product_smarts': '[O:1][P:3](=O)([O-])[O-].[C:2]',
            'rule_name': 'Phosphorylation (ATP Transfer)',
            'category': 'Transferase',
            'enzyme': 'Kinase',
            'ec_number': 'EC 2.7.1.1',
            'description': 'Transfer of phosphate group from ATP to substrate',
        },
        {
            'rule_id': 'BUILTIN-005',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[CH2:1][CH:2][OH:3]>>[CH:1]=[CH:2].[OH:3]',
            'reactant_smarts': '[CH2:1][CH:2][OH:3]',
            'product_smarts': '[CH:1]=[CH:2].[OH:3]',
            'rule_name': 'Dehydrogenation (General)',
            'category': 'Oxidation',
            'enzyme': 'Dehydrogenase',
            'ec_number': 'EC 1.1.-.-',
            'description': 'General dehydrogenation of alcohol to alkene',
        },
        {
            'rule_id': 'BUILTIN-006',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[C:1]([COOH:2])>>[C:1].[C](=O)=O',
            'reactant_smarts': '[C:1]([COOH:2])',
            'product_smarts': '[C:1].[C](=O)=O',
            'rule_name': 'Decarboxylation',
            'category': 'Lyase',
            'enzyme': 'Decarboxylase',
            'ec_number': 'EC 4.1.1.-',
            'description': 'Removal of carboxyl group as CO2',
        },
        {
            'rule_id': 'BUILTIN-007',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[CH:1][NH2:2]>>[CH:1]=[O:3].[NH3]',
            'reactant_smarts': '[CH:1][NH2:2]',
            'product_smarts': '[CH:1]=[O:3].[NH3]',
            'rule_name': 'Amine Deamination',
            'category': 'Oxidation',
            'enzyme': 'Amine Oxidase',
            'ec_number': 'EC 1.4.3.4',
            'description': 'Oxidative deamination of primary amine',
        },
        {
            'rule_id': 'BUILTIN-008',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[C:1][O-:2]>>[C:1][OH:2]',
            'reactant_smarts': '[C:1][O-:2]',
            'product_smarts': '[C:1][OH:2]',
            'rule_name': 'Acid-Base Protonation',
            'category': 'Isomerase',
            'enzyme': 'General Acid-Base Catalyst',
            'ec_number': 'EC 5.-.-.-',
            'description': 'Protonation of alkoxide to form alcohol',
        },
        {
            'rule_id': 'BUILTIN-009',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[C:1]([OH:2])(=[O:3])[C:4].[NH2:5][C:6]>>[C:1]([NH:5][C:6])(=[O:3])[C:4]',
            'reactant_smarts': '[C:1]([OH:2])(=[O:3])[C:4].[NH2:5][C:6]',
            'product_smarts': '[C:1]([NH:5][C:6])(=[O:3])[C:4]',
            'rule_name': 'Amide Bond Formation',
            'category': 'Ligase',
            'enzyme': 'Ligase / Synthetase',
            'ec_number': 'EC 6.3.-.-',
            'description': 'Formation of amide bond between carboxylic acid and amine (ATP-dependent)',
        },
        {
            'rule_id': 'BUILTIN-010',
            'mcsa_id': None,
            'mechanism_id': None,
            'step_id': None,
            'is_reversed': None,
            'reaction_smarts': '[CH3:1][CH2:2][OH:3]>>[CH3:1][CH:2]=[O:4]',
            'reactant_smarts': '[CH3:1][CH2:2][OH:3]',
            'product_smarts': '[CH3:1][CH:2]=[O:4]',
            'rule_name': 'Alcohol to Aldehyde (single product)',
            'category': 'Oxidation',
            'enzyme': 'Alcohol Dehydrogenase',
            'ec_number': 'EC 1.1.1.1',
            'description': 'Oxidation of primary alcohol to aldehyde (single-step)',
        },
    ]


# ---------------------------------------------------------------------------
# Rule Loading from SQLite (preferred) or Excel (fallback)
# ---------------------------------------------------------------------------

# Path to the SQLite database (project db/ directory)
_DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'db', 'custom.db')


def load_rules_from_db(
    max_rules: int = 0,
    source_filter: str = '',
    complete_only: bool = True,
) -> List[Dict[str, Any]]:
    """Load M-CSA reaction rules from SQLite database.

    Much faster than xlsx parsing (~0.5s vs ~12s for 50k rules).
    Falls back gracefully if database is not available.

    Args:
        max_rules: Maximum rules to load (0 = all).
        source_filter: Filter by sourceTag ('mcsa', 'mcsa-test', or '' for all).
        complete_only: Only load rules where ruleCompleteInStep=True (default: True).
            Incomplete rules describe partial step transformations and produce
            excessive noise during search. Complete rules (6K vs 51K total) are
            the high-quality subset that fully describe mechanistic steps.

    Returns:
        List of rule dicts in the same format as load_rules_from_xlsx.
    """
    import sqlite3

    db_path = os.path.abspath(_DB_PATH)
    if not os.path.exists(db_path):
        return []

    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row

    try:
        query = "SELECT * FROM ReactionRule"
        params: tuple = ()
        conditions = []
        if complete_only:
            conditions.append("ruleCompleteInStep = 1")
        if source_filter:
            conditions.append("sourceTag = ?")
            params = params + (source_filter,)
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        if max_rules > 0:
            query += " LIMIT ?"
            params = params + (max_rules,)

        rows = conn.execute(query, params).fetchall()

        loaded: List[Dict[str, Any]] = []
        for row in rows:
            rule: Dict[str, Any] = {
                'rule_id': row['ruleId'] or '',
                'mcsa_id': row['mcsaId'] or None,
                'mechanism_id': row['mechanismId'] or None,
                'step_id': row['stepId'] or None,
                'is_reversed': bool(row['isReversed']),
                'reaction_smarts': row['reactionSmarts'] or '',
                'reactant_smarts': row['reactantSmarts'] or '',
                'product_smarts': row['productSmarts'] or '',
                'rule_arrows': row['ruleArrows'] or '',
                'source': row['sourceTag'] or 'mcsa',
                'rule_name': f"M-CSA Step {row['stepId'] or '?'} "
                             f"(Mech {row['mechanismId'] or '?'})",
                'category': 'Enzyme Catalysis',
            }
            loaded.append(rule)

        return loaded
    finally:
        conn.close()


def load_rules_from_xlsx(
    filepath: str,
    max_rules: int = 0,
    source_tag: str = 'mcsa',
    complete_only: bool = True,
) -> List[Dict[str, Any]]:
    """Load M-CSA reaction rules from an Excel file.

    Args:
        filepath: Path to the rules Excel file.
        max_rules: Maximum number of M-CSA rules to load (0 = all).
        source_tag: Tag identifying the source ('mcsa', 'mcsa-test', etc.).
        complete_only: Only load rules where rule_complete_in_step=True (default: True).
            Incomplete rules describe partial step transformations and produce
            excessive noise during search. Complete rules (6K vs 51K total) are
            the high-quality subset that fully describe mechanistic steps.

    Returns:
        List of rule dicts with reaction_smarts, reactant_smarts, product_smarts,
        rule_id, mcsa_id, mechanism_id, step_id, is_reversed, rule_name, source.
    """
    if not os.path.exists(filepath):
        return []

    wb = openpyxl.load_workbook(filepath, read_only=True, data_only=True)
    ws = wb.active

    rows_iter = ws.iter_rows(values_only=True)
    try:
        headers = next(rows_iter)
    except StopIteration:
        wb.close()
        return []

    # Build column index mapping from headers
    col_map: Dict[str, int] = {}
    for i, h in enumerate(headers):
        if h:
            col_map[h.strip().lower()] = i

    smarts_col_idx = col_map.get('reaction smarts', col_map.get('reaction_smarts', 6))
    rule_id_idx = col_map.get('rule_id', 3)
    mcsa_id_idx = col_map.get('mcsa_id', 0)
    mechanism_id_idx = col_map.get('mechanism_id', 1)
    step_id_idx = col_map.get('step_id', 2)
    is_reversed_idx = col_map.get('is_reversed', 4)
    complete_idx = col_map.get('rule_complete_in_step', 5)
    radical_idx = col_map.get('radical_in_step', 7)
    arrows_idx = col_map.get('rule_arrows', 8)

    skipped_incomplete = 0
    loaded: List[Dict[str, Any]] = []
    for row in rows_iter:
        if max_rules > 0 and len(loaded) >= max_rules:
            break
        if not row or not any(row):
            continue

        # Filter: only complete rules (rule_complete_in_step=True)
        if complete_only:
            complete_val = row[complete_idx] if complete_idx < len(row) else None
            is_complete = bool(complete_val) if complete_val is not None else False
            if not is_complete:
                skipped_incomplete += 1
                continue

        reaction_smarts_raw = row[smarts_col_idx] if smarts_col_idx < len(row) else None
        if not reaction_smarts_raw:
            continue

        smarts_str = str(reaction_smarts_raw).strip()
        if not smarts_str or '>>' not in smarts_str:
            continue

        parts = smarts_str.split('>>', 1)
        reactant_smarts = parts[0].strip()
        product_smarts = parts[1].strip()

        if not reactant_smarts or not product_smarts:
            continue

        rule: Dict[str, Any] = {
            'rule_id': str(row[rule_id_idx]) if rule_id_idx < len(row) and row[rule_id_idx] is not None else '',
            'mcsa_id': row[mcsa_id_idx] if mcsa_id_idx < len(row) else None,
            'mechanism_id': row[mechanism_id_idx] if mechanism_id_idx < len(row) else None,
            'step_id': row[step_id_idx] if step_id_idx < len(row) else None,
            'is_reversed': row[is_reversed_idx] if is_reversed_idx < len(row) else None,
            'rule_complete_in_step': bool(row[complete_idx]) if complete_idx < len(row) and row[complete_idx] else True,
            'reaction_smarts': smarts_str,
            'reactant_smarts': reactant_smarts,
            'product_smarts': product_smarts,
            'source': source_tag,
            'rule_name': f"M-CSA Step {row[step_id_idx] if step_id_idx < len(row) and row[step_id_idx] else '?'} "
                         f"(Mech {row[mechanism_id_idx] if mechanism_id_idx < len(row) and row[mechanism_id_idx] else '?'})",
            'category': 'Enzyme Catalysis',
        }
        loaded.append(rule)

    wb.close()
    if skipped_incomplete > 0:
        logger.info("load_rules_from_xlsx: skipped %d incomplete rules, loaded %d complete rules",
                     skipped_incomplete, len(loaded))
    return loaded


def load_rules(
    max_rules: int = 0,
    rules_file: Optional[str] = None,
    preparse: bool = False,
    complete_only: bool = True,
) -> List[Dict[str, Any]]:
    """Load all reaction rules: built-in + M-CSA from SQLite (preferred) or Excel (fallback).

    Args:
        max_rules: Maximum M-CSA rules to load (0 = all).
        rules_file: Path to rules file (fallback if DB not available).
        preparse: Whether to pre-parse simplified SMARTS components for faster matching.
        complete_only: Only load rules where rule_complete_in_step=True (default: True).
            Incomplete rules (~45K out of 51K) describe partial step transformations
            and produce excessive noise during search. Setting True reduces rules from
            ~51K to ~6K, dramatically improving search speed and result quality.

    Returns:
        List of rule dicts.  Built-in rules come first, then M-CSA rules.
    """
    from rdkit.Chem import AllChem

    builtin = get_builtin_rules()

    # Reserve room for builtin rules when max_rules is set
    mcsa_limit = (max_rules - len(builtin)) if max_rules > 0 else 0

    # Try SQLite first (fast: ~0.5s for 6k rules)
    mcsa = load_rules_from_db(max_rules=mcsa_limit, complete_only=complete_only)

    # Fallback to xlsx if database is empty or unavailable
    if not mcsa:
        filepath = rules_file or RULES_FILE
        mcsa = load_rules_from_xlsx(filepath, max_rules=mcsa_limit, source_tag='mcsa',
                                     complete_only=complete_only)

    rules = builtin + mcsa

    logger.info("load_rules: %d builtin + %d M-CSA (complete_only=%s) = %d total",
                len(builtin), len(mcsa), complete_only, len(rules))

    if preparse:
        _preparse_rules(rules)

    return rules


def _preparse_rules(rules: List[Dict[str, Any]], min_atoms: int = 3) -> None:
    """Pre-parse SMARTS components for each rule — Upload source code approach.

    Reproduces the SMARTS preprocessing from Upload models.py PredictionRun.apply_rules:
    1. Split reaction_smarts on '>>' into reactant/product sides
    2. Split each side on '.' into multi-molecule fragments
    3. Parse each fragment with Chem.MolFromSmarts() into mol objects
    4. Analyze atom maps: rule_am -> part_idx, rule_am -> atom, rule_am -> part_size
    5. Compute bond changes between reactant and product sides
    6. Compute strict reaction center atoms (involved in bond changes)
    """
    from rdkit import Chem

    for rule in rules:
        r_smarts = rule.get('reactant_smarts', '')
        p_smarts = rule.get('product_smarts', '')
        full_smarts = rule.get('reaction_smarts', '')

        # ----------------------------------------------------------------
        # Upload approach (models.py line 946-974): Parse SMARTS parts with
        # atom mapping analysis and bond change computation
        # ----------------------------------------------------------------
        # rule_parts: list of Mol objects for each reactant fragment
        rule['_rule_parts'] = []
        # reactant_smarts_rule_parts: list of original SMARTS strings
        rule['_reactant_smarts_rule_parts'] = []
        # rule_am_to_part_id: {atom_map_num: part_index}
        rule['_rule_am_to_part_id'] = {}
        # rule_am_to_atom: ({am: atom_reactant_side}, {am: atom_product_side})
        rule['_rule_am_to_atom'] = ({}, {})
        # rule_am_to_part_size: {atom_map_num: num_heavy_atoms_in_part}
        rule['_rule_am_to_part_size'] = {}

        reactant_bonds_ams = set()
        for part_idx, part_str in enumerate(r_smarts.split('.')):
            part_str = part_str.strip()
            rule['_reactant_smarts_rule_parts'].append(part_str)
            part_mol = Chem.MolFromSmarts(part_str)
            if part_mol is not None:
                rule['_rule_parts'].append(part_mol)
                for atom in part_mol.GetAtoms():
                    am = atom.GetAtomMapNum()
                    if am:
                        rule['_rule_am_to_part_id'][am] = part_idx
                        rule['_rule_am_to_atom'][0][am] = atom
                        rule['_rule_am_to_part_size'][am] = part_mol.GetNumHeavyAtoms()
                for bond in part_mol.GetBonds():
                    reactant_bonds_ams.add(
                        frozenset((bond.GetBeginAtom().GetAtomMapNum(),
                                   bond.GetEndAtom().GetAtomMapNum())))

        # product_rule_parts: list of Mol objects for each product fragment
        rule['_product_rule_parts'] = []
        rule['_product_smarts_rule_parts'] = []
        product_bonds_ams = set()
        for part_idx, part_str in enumerate(p_smarts.split('.')):
            part_str = part_str.strip()
            rule['_product_smarts_rule_parts'].append(part_str)
            part_mol = Chem.MolFromSmarts(part_str)
            if part_mol is not None:
                rule['_product_rule_parts'].append(part_mol)
                for atom in part_mol.GetAtoms():
                    am = atom.GetAtomMapNum()
                    if am:
                        rule['_rule_am_to_atom'][1][am] = atom
                for bond in part_mol.GetBonds():
                    product_bonds_ams.add(
                        frozenset((bond.GetBeginAtom().GetAtomMapNum(),
                                   bond.GetEndAtom().GetAtomMapNum())))

        # Upload approach (models.py line 972-974): Compute bond changes
        #   new_bonds = product_bonds - reactant_bonds (formed)
        #   broken_bonds = reactant_bonds - product_bonds (cleaved)
        rule['_bond_changes_am'] = (
            product_bonds_ams - reactant_bonds_ams,   # new bonds
            reactant_bonds_ams - product_bonds_ams,    # broken bonds
        )

        # Upload approach (models.py line 974): Strict reaction center atoms
        # = all atom map numbers that participate in bond changes
        rule['_strict_ams'] = {
            am for bonds in rule['_bond_changes_am'] for bond in bonds for am in bond
        }

        # Heavy-atom bond changes (replace H with connected heavy atom)
        new_bonds_heavy = set()
        for bond_ams in rule['_bond_changes_am'][0]:
            atoms_am = []
            for atom_am in bond_ams:
                atom = rule['_rule_am_to_atom'][0].get(atom_am)
                if atom and atom.GetSymbol() == "H":
                    try:
                        reactant_bond = [
                            b for b in rule['_bond_changes_am'][1]
                            for a in b if a == atom_am
                        ][0]
                        atom_am = [a for a in reactant_bond if a != atom_am][0]
                    except (IndexError, KeyError):
                        pass
                atoms_am.append(atom_am)
            if len(atoms_am) == 2:
                new_bonds_heavy.add(frozenset(atoms_am))
        rule['_bond_changes_am_heavy'] = (new_bonds_heavy, rule['_bond_changes_am'][1])

        # _num_reactants: used by mechanism_search.py for template setup
        rule['_num_reactants'] = r_smarts.count('.') + 1


def merge_rule_parts(rule_parts, repeated_ams_groups):
    """Merge multiple rule parts that match the same reactant molecule.

    Ported from Upload models.py line 1689-1709.
    When the same reactant matches multiple rule parts, we need to combine
    those parts into a single SMARTS pattern for RunReactants.

    Args:
        rule_parts: List of RDKit Mol objects (rule parts to merge)
        repeated_ams_groups: List of lists of atom map numbers that map to
                            the same reactant atom

    Returns:
        Combined RWMol with duplicate atoms merged and bonds transferred
    """
    from rdkit import Chem

    combined = _combine_mols_list(rule_parts)
    new_rule_part = Chem.RWMol(combined)
    am_to_atom = {atom.GetAtomMapNum(): atom for atom in new_rule_part.GetAtoms()}
    atoms_to_delete = set()

    for rule_ams in repeated_ams_groups:
        rule_ams = list(rule_ams)
        if not rule_ams:
            continue
        atom_to_keep = am_to_atom.get(rule_ams.pop(0))
        if atom_to_keep is None:
            continue
        for rule_am in rule_ams:
            atom_to_delete = am_to_atom.get(rule_am)
            if atom_to_delete is None:
                continue
            atoms_to_delete.add(atom_to_delete.GetIdx())
            bonds = atom_to_delete.GetBonds()
            # If kept atom is dummy, reset to carbon
            if bonds and atom_to_keep.GetAtomicNum() == 0:
                atom_to_keep.SetAtomicNum(6)
                atom_to_keep.SetQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(6))
            for bond in bonds:
                other_atom = bond.GetOtherAtom(atom_to_delete)
                if not new_rule_part.GetBondBetweenAtoms(
                        atom_to_keep.GetIdx(), other_atom.GetIdx()):
                    new_rule_part.AddBond(
                        atom_to_keep.GetIdx(), other_atom.GetIdx(), bond.GetBondType())

    for atom_idx in sorted(atoms_to_delete, reverse=True):
        new_rule_part.RemoveAtom(atom_idx)

    return new_rule_part


def _combine_mols_list(mols):
    """Combine multiple RDKit Mol objects into a single Mol.

    Equivalent to Upload common.utils.combine_mols_list.
    """
    from rdkit import Chem
    combined = Chem.RWMol()
    atom_offset = 0
    for mol in mols:
        for atom in mol.GetAtoms():
            new_atom = Chem.Atom(atom)
            combined.AddAtom(new_atom)
        for bond in mol.GetBonds():
            combined.AddBond(
                bond.GetBeginAtomIdx() + atom_offset,
                bond.GetEndAtomIdx() + atom_offset,
                bond.GetBondType())
        atom_offset += mol.GetNumAtoms()
    return combined.GetMol()


