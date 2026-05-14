#!/usr/bin/env python3
"""
EzMechanism FindMCS-based Molecule Matching Utilities

Ported from the original EzMechanism paper (Nature Methods 2023).
These functions use RDKit's Maximum Common Substructure (MCS) algorithm
for element-type-aware molecule matching, which is more flexible than
simple substructure matching.

Key features:
- FindMCS with CompareIsotopes: Matches molecules by encoding element types in isotopes
- Special molecule handling: Protons (H+), water (H₂O), hydronium (H₃O⁺)
- One-to-many matching: find_n_matching_mols with isotope masking for recursive pairing
- Atom mapping: Returns both idx_mapping and atom_map_num_mapping
"""

import logging
from typing import Dict, Iterable, Tuple, List, Optional
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Geometry import Point2D, Point3D

logger = logging.getLogger("script")


# ============================================================================
# Special molecule detection (from original predict_mechanism_util.py)
# ============================================================================

def mol_is_proton(mol: Chem.Mol) -> bool:
    """Check if mol is a proton (single H atom)."""
    return mol.GetNumAtoms() == 1 and mol.GetAtoms()[0].GetSymbol() == "H"


def first_proton(mol: Chem.Mol) -> Optional[Chem.Atom]:
    """Return the first proton (H atom) in the molecule."""
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            return atom
    return None


def mol_is_water_species(mol: Chem.Mol) -> bool:
    """Check if mol is a water molecule, hydroxide or hydronium."""
    return (mol.GetNumHeavyAtoms() == 1 and
            any([atom.GetSymbol() == "O" for atom in mol.GetAtoms()]))


def mol_is_hydronium(mol: Chem.Mol) -> bool:
    """Check if mol is a hydronium ion (H₃O⁺)."""
    return mol.GetNumHeavyAtoms() == 1 and any(
              [atom.GetSymbol() == "O" and
               atom.GetFormalCharge() == 1 for atom in mol.GetAtoms()])


def mol_is_water_or_proton(mol: Chem.Mol) -> bool:
    """Check if mol is a water molecule or a proton."""
    return mol_is_proton(mol) or mol_is_water_species(mol)


# ============================================================================
# Core matching functions (from original predict_mechanism_util.py)
# ============================================================================

def prepare_mols_for_comparison(mol_list_a: Iterable[Chem.Mol],
                                  mol_list_b: Iterable[Chem.Mol]) -> Dict[Chem.Mol, Chem.Mol]:
    """Prepare molecules for MCS-based comparison.

    Creates RWMol copies, sanitizes them, then:
    - Sets isotope = atomicNum (for CompareIsotopes matching)
    - Sets atomMapNum = 0 (to eliminate atom map interference)

    This is the KEY technique from the EzMechanism paper for handling
    SMILES vs SMARTS with atom map numbers.
    """
    comp = {mol: Chem.RWMol(mol) for group in [mol_list_a, mol_list_b] for mol in group}
    
    errors = []
    for mol in comp.values():
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            errors.append(str(e))
    
    for mol in comp.values():
        for atom in mol.GetAtoms():
            atom.SetIsotope(atom.GetAtomicNum())
            atom.SetAtomMapNum(0)
    
    return comp


def find_matching_mols(mol_list_a: Iterable[Chem.Mol],
                        mol_list_b: Iterable[Chem.Mol])\
        -> Dict[Chem.Mol, Tuple[Chem.Mol, Dict[int, int], Dict[int, int]]]:
    """
    Find the matching mols in two list of mols

    Returns a Dict where mols of the first list are keys and the values are a tuple with the mol
    in the second list and two mappings. The first mapping is based on the the atom idx,
    the second is based on the the atom map number

    :param mol_list_a: a list of rdkit Mol objects
    :param mol_list_b: a list of rdkit Mol objects
    :return: a Dict with the matches and atom mappings

    [对齐 Upload] comp 只做 SanitizeMol，不设置 isotope/atomMapNum。
    仅 find_n_matching_mols 才设置 isotope/atomMapNum（用于 CompareIsotopes）。
    """
    # comparison is based on an aromatized version of the molecules
    comp = {mol: Chem.RWMol(mol) for group in [mol_list_a, mol_list_b] for mol in group}
    [Chem.SanitizeMol(mol) for mol in comp.values()]

    # first calculate all vs all similarity
    similarities = []
    for mol_a in mol_list_a:
        for mol_b in mol_list_b:
            # explicit comparison for single atoms since mcs only calculated for molecules
            if mol_a.GetNumAtoms() == 1 and mol_b.GetNumAtoms() == 1:
                substructure = mol_a
            # map protons to hydronium
            elif mol_is_proton(mol_a) and mol_is_hydronium(mol_b) or\
                    mol_is_proton(mol_b) and mol_is_hydronium(mol_a):
                substructure = mol_a if mol_is_proton(mol_a) else mol_b
            else:
                # timeout is necessary because the algorithm gets stuck sometimes
                mcs = Chem.rdFMCS.FindMCS([comp[mol_a], comp[mol_b]], timeout=2)
                substructure = Chem.MolFromSmarts(mcs.smartsString)
            if substructure.GetNumAtoms():
                similarities.append((mol_a, mol_b, substructure))

    # then pair molecules starting with the most similar ones
    similarities.sort(key=lambda x: -x[2].GetNumAtoms())
    pairs = {}
    no = 0
    while similarities:
        mol_a, mol_b, substructure = similarities[0]
        logger.debug("Paired {} with {}. {}, {}".format(Chem.MolToSmiles(mol_a), Chem.MolToSmiles(mol_b),
                                                        substructure.GetNumAtoms(), Chem.MolToSmiles(substructure)))
        no += 1
        # 如果两个分子都是质子，则直接将它们的原子索引和映射编号进行一对一映射。
        if mol_is_proton(mol_a) and mol_is_proton(mol_b):
            idx_mapping = {mol_a.GetAtoms()[0].GetIdx(): mol_b.GetAtoms()[0].GetIdx()}
            atom_num_map_mapping = {mol_a.GetAtoms()[0].GetIdx(): mol_b.GetAtoms()[0].GetIdx()}
        # 如果其中一个分子是质子，另一个是水合氢离子，则找到各自的质子原子并建立映射关系。
        elif mol_is_proton(mol_a) and mol_is_hydronium(mol_b) or\
                mol_is_proton(mol_b) and mol_is_hydronium(mol_a):
            first_proton_a, first_proton_b = first_proton(mol_a), first_proton(mol_b)
            idx_mapping = {first_proton_a.GetIdx(): first_proton_b.GetIdx()}
            atom_num_map_mapping = {first_proton_a.GetAtomMapNum(): first_proton_b.GetAtomMapNum()}
        # 对于其他分子对，获取它们与子结构匹配的原子顺序，生成原子索引映射和原子映射编号映射。
        else:
            order_a = comp[mol_a].GetSubstructMatch(substructure)
            order_b = comp[mol_b].GetSubstructMatch(substructure)
            idx_mapping = {order_a[no]: order_b[no] for no in range(len(order_a))}
            atom_num_map_mapping = {mol_a.GetAtomWithIdx(idx_a).GetAtomMapNum():
                                    mol_b.GetAtomWithIdx(idx_b).GetAtomMapNum()
                                    for idx_a, idx_b in idx_mapping.items()}
        pairs[mol_a] = (mol_b, idx_mapping, atom_num_map_mapping)
        # 从 similarities 列表中移除已经匹配的分子对，避免重复处理。
        similarities = [sim for sim in similarities if sim[0] is not mol_a and sim[1] is not mol_b]

    return pairs


def find_n_matching_mols(mol_list_a: Iterable[Chem.Mol],
                         mol_list_b: Iterable[Chem.Mol],
                         mol_maps: Optional[Dict] = None,
                         comp: Optional[Dict] = None) -> Dict:
    """Find matching mols with one-to-many support using isotope masking.

    Unlike find_matching_mols (one-to-one), this allows multiple molecules in
    list_b to match parts of a single molecule in list_a.

    Uses isotope masking (1000/2000) after each match to prevent the same atoms
    from being matched again in recursive calls.

    Returns Dict[(mol_a, mol_b)] = (idx_mapping, atom_num_map_mapping)
    """
    if comp is None:
        comp = prepare_mols_for_comparison(mol_list_a, mol_list_b)

    similarities = []
    mol_maps = {} if mol_maps is None else mol_maps
    for mol_a in mol_list_a:
        for mol_b in mol_list_b:
            try:
                mcs = Chem.rdFMCS.FindMCS(
                    [comp[mol_a], comp[mol_b]], timeout=2,
                    atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes
                )
                if mcs.queryMol is not None:
                    similarities.append((mol_a, mol_b, mcs.queryMol))
            except Exception:
                continue

    if similarities:
        similarities.sort(key=lambda x: -x[2].GetNumAtoms())
        mol_a, mol_b, query_mol = similarities[0]

        try:
            order_a = comp[mol_a].GetSubstructMatch(query_mol)
            order_b = comp[mol_b].GetSubstructMatch(query_mol)
            idx_mapping = {order_a[no]: order_b[no] for no in range(len(order_a))}
            atom_num_map_mapping = {
                mol_a.GetAtomWithIdx(idx_a).GetAtomMapNum():
                mol_b.GetAtomWithIdx(idx_b).GetAtomMapNum()
                for idx_a, idx_b in idx_mapping.items()
            }
        except Exception:
            return mol_maps

        if (mol_a, mol_b) not in mol_maps:
            mol_maps[(mol_a, mol_b)] = (idx_mapping, atom_num_map_mapping)

        # Mask matched atoms with different isotope values to prevent re-matching
        for idx_a, idx_b in idx_mapping.items():
            try:
                comp[mol_a].GetAtomWithIdx(idx_a).SetIsotope(1000)
                comp[mol_b].GetAtomWithIdx(idx_b).SetIsotope(2000)
            except Exception:
                pass

        # Recursive call to find additional matches
        return find_n_matching_mols(mol_list_a, mol_list_b, mol_maps, comp)
    else:
        return mol_maps


def compare_mols(mol_a: Chem.Mol, mol_b: Chem.Mol,
                 atom_compare=Chem.rdFMCS.AtomCompare.CompareElements,
                 bond_compare=Chem.rdFMCS.BondCompare.CompareOrder) -> float:
    """Return the percentage of similar atoms between two molecules."""
    try:
        mcs = Chem.rdFMCS.FindMCS([mol_a, mol_b],
                                 atomCompare=atom_compare,
                                 bondCompare=bond_compare, timeout=2)
        max_atoms = max(mol_a.GetNumAtoms(), mol_b.GetNumAtoms())
        if max_atoms == 0:
            return 0.0
        return mcs.numAtoms / max_atoms
    except Exception:
        return 0.0


# ============================================================================
# Convenience functions for SMILES-based matching
# ============================================================================

def smiles_to_mol(smiles: str, sanitize: bool = True, add_hs: bool = False) -> Optional[Chem.Mol]:
    """Convert a SMILES string to an RDKit Mol object."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if sanitize:
            Chem.SanitizeMol(mol)
        if add_hs:
            mol = Chem.AddHs(mol)
        return mol
    except Exception:
        return None


def find_matching_mols_from_smiles(smiles_list_a: List[str],
                                    smiles_list_b: List[str]) -> Dict:
    """Find matching molecules between two SMILES lists.

    Wrapper around find_matching_mols that converts SMILES to Mol objects first.
    """
    mols_a = [m for s in smiles_list_a if (m := smiles_to_mol(s)) is not None]
    mols_b = [m for s in smiles_list_b if (m := smiles_to_mol(s)) is not None]
    if not mols_a or not mols_b:
        return {}
    return find_matching_mols(mols_a, mols_b)


def match_smiles_with_smarts(user_smiles: str, rule_smarts: str) -> Tuple[bool, float, int]:
    """Match a user's SMILES against a rule's SMARTS pattern using FindMCS.

    This is the KEY function for the EzMechanism search algorithm.
    It handles the SMILES vs SMARTS difference by:
    1. Parsing both into Mol objects
    2. Setting isotope=atomicNum and atomMapNum=0 (normalize)
    3. Running FindMCS with CompareIsotopes

    Returns: (matched: bool, score: float, num_matched_atoms: int)
    """
    try:
        # IMPORTANT: simplify M-CSA SMARTS before parsing
        # M-CSA SMARTS uses :N atom map format that can cause parse errors
        # when combined with charges (e.g., [O:2-] fails, [O:2;-1] works)
        from shared import simplify_smarts
        simplified = simplify_smarts(rule_smarts)

        # Parse simplified SMARTS pattern (rule side)
        rule_mol = Chem.MolFromSmarts(simplified)
        if rule_mol is None:
            # Fallback: try original (may work for some formats)
            rule_mol = Chem.MolFromSmarts(rule_smarts)
        if rule_mol is None:
            return (False, 0.0, 0)

        # Parse user SMILES
        user_mol = Chem.MolFromSmiles(user_smiles)
        if user_mol is None:
            return (False, 0.0, 0)

        # Create RWMol copies for comparison
        rule_copy = Chem.RWMol(rule_mol)
        user_copy = Chem.RWMol(user_mol)

        # Sanitize
        try:
            Chem.SanitizeMol(rule_copy)
        except Exception:
            pass
        try:
            Chem.SanitizeMol(user_copy)
        except Exception:
            pass

        # Normalize: set isotope = atomicNum, clear atom map numbers
        for mol in [rule_copy, user_copy]:
            for atom in mol.GetAtoms():
                atom.SetIsotope(atom.GetAtomicNum())
                atom.SetAtomMapNum(0)

        # Run FindMCS with CompareIsotopes
        mcs = Chem.rdFMCS.FindMCS([rule_copy, user_copy], timeout=2,
                                    atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes)

        if mcs.queryMol is None:
            return (False, 0.0, 0)

        num_matched = mcs.numAtoms

        # Score based on coverage ratio
        min_atoms = min(rule_copy.GetNumAtoms(), user_copy.GetNumAtoms())
        max_atoms = max(rule_copy.GetNumAtoms(), user_copy.GetNumAtoms())
        if max_atoms == 0:
            return (False, 0.0, 0)

        # Coverage score: how much of the rule pattern is matched by the user molecule
        coverage = num_matched / min_atoms if min_atoms > 0 else 0

        # Quality score: how much of the user molecule is used
        quality = num_matched / max_atoms

        # Combined score (higher = better match)
        score = 0.4 * coverage + 0.6 * quality

        return (True, round(score, 3), num_matched)
    except Exception:
        return (False, 0.0, 0)


def match_smiles_with_smarts_components(user_smiles: str,
                                         rule_reactant_smarts: str) -> Tuple[bool, float, int, int]:
    """Match user SMILES against each component of a multi-molecule SMARTS pattern.

    For M-CSA rules, the reactant SMARTS often has multiple components separated by '.'
    (e.g., "substrate.cofactor.residue"). This function checks the user SMILES
    against each component and returns the best match.

    Returns: (matched: bool, best_score: float, num_matched_atoms: int, matched_component_idx: int)
    """
    # Validate user SMILES first
    user_mol = Chem.MolFromSmiles(user_smiles)
    if user_mol is None:
        logger.warning("match_smiles_with_smarts_components: invalid user SMILES: %s", user_smiles[:50])
        return (False, 0.0, 0, -1)

    if not rule_reactant_smarts:
        logger.warning("match_smiles_with_smarts_components: empty rule_reactant_smarts")
        return (False, 0.0, 0, -1)

    components = [c.strip() for c in rule_reactant_smarts.split('.') if c.strip()]
    if not components:
        return (False, 0.0, 0, -1)

    best_score = 0.0
    best_atoms = 0
    best_idx = -1

    for idx, comp_smarts in enumerate(components):
        # Skip very long components (likely residue SMARTS that are too specific)
        if len(comp_smarts) > 200:
            continue

        # Each component is independently error-handled so one bad component
        # doesn't kill the entire function (BUG FIX: was previously a single
        # try/except wrapping the ENTIRE function body)
        try:
            # DO NOT simplify here — match_smiles_with_smarts already simplifies
            # internally. Double simplification corrupts ring closure digits
            # (e.g., C1CCCCC1 → CCCCCC1 after 1st simplify → CCCCCC after 2nd)
            matched, score, num_atoms = match_smiles_with_smarts(user_smiles, comp_smarts)
            if matched and num_atoms > best_atoms:
                best_score = score
                best_atoms = num_atoms
                best_idx = idx
                logger.debug("  comp[%d] matched: score=%.3f, atoms=%d, smarts=%s",
                             idx, score, num_atoms, comp_smarts[:80])
        except Exception as e:
            logger.debug("  comp[%d] failed: %s, smarts=%s",
                         idx, str(e), comp_smarts[:80])
            continue

    if best_idx >= 0:
        return (True, best_score, best_atoms, best_idx)
    return (False, 0.0, 0, -1)


# ============================================================================
# Quick similarity check (lightweight alternative to full MCS matching)
# ============================================================================

def quick_similarity_check(smiles_a: str, smiles_b: str) -> float:
    """Quick similarity check using element-level MCS comparison.

    Returns a score between 0 and 1 (1 = identical).
    Much faster than full FindMCS matching.
    """
    try:
        mol_a = Chem.MolFromSmiles(smiles_a)
        mol_b = Chem.MolFromSmiles(smiles_b)
        if mol_a is None or mol_b is None:
            return 0.0
        return compare_mols(mol_a, mol_b)
    except Exception:
        return 0.0
