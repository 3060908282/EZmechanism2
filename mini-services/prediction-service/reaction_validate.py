#!/usr/bin/env python3
"""
Validate reaction atom balance and detect bond changes (cleaved/formed bonds).

Usage:
    # Atom balance only
    python reaction_validate.py --reactants "CCO.O" --products "CC=O.O"

    # Full mapping check with bond change detection
    python reaction_validate.py --check-mapping --reactants "CC(=O)O.CCO" --products "CC(=O)OCC.O"

Output (JSON):
    {
        "valid": true,
        "reactant_atoms": {...},
        "product_atoms": {...},
        "imbalanced": [],
        "bond_changes": {
            "cleaved": [{"elements": "C-O", "bond_type": 1.0, ...}],
            "formed": [{"elements": "C-O", "bond_type": 1.0, ...}],
            "changed": [{"elements": "C-C", "from_type": 1.0, "to_type": 2.0, ...}],
            "mapping_coverage": 1.0,
            "mapped_atoms": 17,
            "total_atoms": 17
        },
        ...
    }
"""

import sys
import json
import argparse


def count_atoms_with_h(mol) -> dict:
    """Count all atoms including implicit hydrogens."""
    atoms = {}
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol()
        atoms[sym] = atoms.get(sym, 0) + 1
        total_h = atom.GetTotalNumHs()
        if total_h > 0:
            atoms["H"] = atoms.get("H", 0) + total_h
    return atoms


def get_formula(mol) -> str:
    """Get molecular formula string from RDKit mol."""
    try:
        from rdkit.Chem import rdMolDescriptors
        return rdMolDescriptors.CalcMolFormula(mol)
    except Exception:
        return "?"


def parse_multi_smiles(smiles_str: str) -> list:
    """Split dot-separated SMILES into individual molecules."""
    parts = smiles_str.strip().split(".")
    results = []
    for part in parts:
        part = part.strip()
        if not part:
            continue
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(part)
            if mol is None:
                results.append({
                    "smiles": part, "formula": "Invalid",
                    "atom_count": 0, "valid": False
                })
                continue
            atoms = count_atoms_with_h(mol)
            results.append({
                "smiles": part, "formula": get_formula(mol),
                "atom_count": sum(atoms.values()), "valid": True,
                "atoms": atoms, "mol": mol
            })
        except Exception:
            results.append({
                "smiles": part, "formula": "Error",
                "atom_count": 0, "valid": False
            })
    return results


def _build_atom_mapping_greedy(r_combined, p_combined):
    """Build full atom mapping using greedy assignment.
    
    Groups atoms by element type, then within each group sorts by
    a connectivity signature (degree + neighbor element multiset).
    Matches atoms with identical signatures first, then falls back
    to degree-only matching for remaining unmatched atoms.
    """
    from collections import Counter
    
    def connectivity_sig(mol, idx):
        atom = mol.GetAtomWithIdx(idx)
        nbr_elems = sorted([n.GetSymbol() for n in atom.GetNeighbors()])
        return (atom.GetDegree(), tuple(nbr_elems))
    
    def connectivity_sig_extended(mol, idx):
        atom = mol.GetAtomWithIdx(idx)
        nbr_counts = Counter(n.GetSymbol() for n in atom.GetNeighbors())
        nbr_key = tuple(sorted(nbr_counts.items()))
        return (atom.GetDegree(), atom.GetIsAromatic(), nbr_key)
    
    # Group atoms by element
    r_by_elem = {}
    p_by_elem = {}
    for atom in r_combined.GetAtoms():
        r_by_elem.setdefault(atom.GetSymbol(), []).append(atom.GetIdx())
    for atom in p_combined.GetAtoms():
        p_by_elem.setdefault(atom.GetSymbol(), []).append(atom.GetIdx())
    
    mapping = {}  # r_idx -> p_idx
    used_p = set()
    
    for elem in r_by_elem:
        if elem not in p_by_elem:
            continue
        r_list = r_by_elem[elem]
        p_list = list(p_by_elem[elem])
        
        # Phase 1: Match by extended signature (degree + aromaticity + neighbor counts)
        r_sigs_ext = [(connectivity_sig_extended(r_combined, i), i) for i in r_list]
        p_sigs_ext = [(connectivity_sig_extended(p_combined, i), i) for i in p_list]
        
        # Build lookup: sig -> list of (p_idx) for exact matching
        p_by_sig = {}
        for sig, p_idx in p_sigs_ext:
            p_by_sig.setdefault(sig, []).append(p_idx)
        
        unmatched_r = []
        for sig, r_idx in r_sigs_ext:
            if sig in p_by_sig:
                candidates = [p_idx for p_idx in p_by_sig[sig] if p_idx not in used_p]
                if candidates:
                    chosen = candidates[0]
                    mapping[r_idx] = chosen
                    used_p.add(chosen)
                    continue
            unmatched_r.append(r_idx)
        
        # Phase 2: Match remaining by degree only
        if unmatched_r:
            p_remaining = [i for i in p_list if i not in used_p]
            p_by_degree = {}
            for p_idx in p_remaining:
                deg = p_combined.GetAtomWithIdx(p_idx).GetDegree()
                p_by_degree.setdefault(deg, []).append(p_idx)
            
            for r_idx in unmatched_r:
                deg = r_combined.GetAtomWithIdx(r_idx).GetDegree()
                if deg in p_by_degree:
                    candidates = [p_idx for p_idx in p_by_degree[deg] if p_idx not in used_p]
                    if candidates:
                        chosen = candidates[0]
                        mapping[r_idx] = chosen
                        used_p.add(chosen)
                        p_by_degree[deg].remove(chosen)
    
    return mapping


def analyze_bond_changes(r_mols_data: list, p_mols_data: list) -> dict:
    """Detect cleaved, formed, and changed bonds between reactants and products.
    
    Uses greedy atom mapping (element + connectivity signature) to map ALL
    atoms between combined reactant and product molecules, then compares bonds.
    """
    from rdkit import Chem

    # Extract valid RDKit molecules
    r_mols = [m["mol"] for m in r_mols_data if m.get("valid") and "mol" in m]
    p_mols = [m["mol"] for m in p_mols_data if m.get("valid") and "mol" in m]

    if not r_mols or not p_mols:
        return {
            "cleaved": [], "formed": [], "changed": [],
            "mapped_atoms": 0, "total_atoms": 0,
            "mapping_coverage": 0, "unmapped_reactants": 0, "unmapped_products": 0
        }

    # Combine all fragments into single molecules
    r_combined = r_mols[0]
    for m in r_mols[1:]:
        r_combined = Chem.CombineMols(r_combined, m)
    p_combined = p_mols[0]
    for m in p_mols[1:]:
        p_combined = Chem.CombineMols(p_combined, m)

    total_r = r_combined.GetNumAtoms()
    total_p = p_combined.GetNumAtoms()

    # Build atom mapping using greedy algorithm
    r_to_p = _build_atom_mapping_greedy(r_combined, p_combined)
    p_to_r = {v: k for k, v in r_to_p.items()}
    
    mapped_count = len(r_to_p)

    # Compare bonds on mapped atoms
    cleaved = []
    formed = []
    changed = []

    # Check all bonds in reactants
    for bond in r_combined.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in r_to_p and j in r_to_p:
            pi, pj = r_to_p[i], r_to_p[j]
            p_bond = p_combined.GetBondBetweenAtoms(pi, pj)
            r_elem = f"{r_combined.GetAtomWithIdx(i).GetSymbol()}-{r_combined.GetAtomWithIdx(j).GetSymbol()}"
            r_btype = round(bond.GetBondTypeAsDouble(), 1)
            if p_bond is None:
                # Bond exists in reactants but not in products -> cleaved
                cleaved.append({
                    "elements": r_elem,
                    "bond_type": r_btype,
                    "bond_label": _bond_type_label(r_btype),
                    "r_atom_indices": [i, j],
                    "p_atom_indices": [pi, pj],
                })
            else:
                p_btype = round(p_bond.GetBondTypeAsDouble(), 1)
                if r_btype != p_btype:
                    # Bond type changed
                    changed.append({
                        "elements": r_elem,
                        "from_type": r_btype,
                        "to_type": p_btype,
                        "from_label": _bond_type_label(r_btype),
                        "to_label": _bond_type_label(p_btype),
                        "r_atom_indices": [i, j],
                        "p_atom_indices": [pi, pj],
                    })

    # Check all bonds in products (for newly formed bonds)
    for bond in p_combined.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in p_to_r and j in p_to_r:
            ri, rj = p_to_r[i], p_to_r[j]
            r_bond = r_combined.GetBondBetweenAtoms(ri, rj)
            if r_bond is None:
                # Bond exists in products but not in reactants -> formed
                p_elem = f"{p_combined.GetAtomWithIdx(i).GetSymbol()}-{p_combined.GetAtomWithIdx(j).GetSymbol()}"
                p_btype = round(bond.GetBondTypeAsDouble(), 1)
                formed.append({
                    "elements": p_elem,
                    "bond_type": p_btype,
                    "bond_label": _bond_type_label(p_btype),
                    "r_atom_indices": [ri, rj],
                    "p_atom_indices": [i, j],
                })

    unmapped_r = total_r - mapped_count
    unmapped_p = total_p - mapped_count
    coverage = mapped_count / max(total_r, total_p) if max(total_r, total_p) > 0 else 0

    return {
        "cleaved": cleaved,
        "formed": formed,
        "changed": changed,
        "mapped_atoms": mapped_count,
        "total_atoms": max(total_r, total_p),
        "mapping_coverage": round(coverage, 2),
        "unmapped_reactants": unmapped_r,
        "unmapped_products": unmapped_p,
        "summary": _build_bond_summary(cleaved, formed, changed),
    }


def _bond_type_label(btype: float) -> str:
    """Convert bond type number to human-readable label."""
    labels = {
        1.0: "Single",
        1.5: "Aromatic",
        2.0: "Double",
        3.0: "Triple",
        0.0: "Zero",
        12.0: "Dative/Coordinate",
    }
    return labels.get(btype, f"Order {btype}")


def _build_bond_summary(cleaved, formed, changed) -> str:
    """Build a human-readable summary of bond changes."""
    parts = []
    if cleaved:
        parts.append(f"{len(cleaved)} bond(s) cleaved")
    if formed:
        parts.append(f"{len(formed)} bond(s) formed")
    if changed:
        parts.append(f"{len(changed)} bond(s) changed order")
    if not parts:
        return "No bond changes detected — reactants and products have identical bonding"
    return "; ".join(parts)


def validate_reaction(reactants: str, products: str, check_mapping: bool = False) -> dict:
    """Compare atom counts between reactants and products."""
    reactant_mols = parse_multi_smiles(reactants)
    product_mols = parse_multi_smiles(products)

    # Sum up atom counts
    reactant_atoms = {}
    reactant_formulas = []
    for m in reactant_mols:
        if m.get("valid") and "atoms" in m:
            for elem, count in m["atoms"].items():
                reactant_atoms[elem] = reactant_atoms.get(elem, 0) + count
        reactant_formulas.append(m["formula"])

    product_atoms = {}
    product_formulas = []
    for m in product_mols:
        if m.get("valid") and "atoms" in m:
            for elem, count in m["atoms"].items():
                product_atoms[elem] = product_atoms.get(elem, 0) + count
        product_formulas.append(m["formula"])

    # Find imbalanced elements
    all_elements = sorted(set(list(reactant_atoms.keys()) + list(product_atoms.keys())))
    imbalanced = []
    for elem in all_elements:
        r_count = reactant_atoms.get(elem, 0)
        p_count = product_atoms.get(elem, 0)
        if r_count != p_count:
            diff = p_count - r_count
            imbalanced.append({
                "element": elem,
                "reactant_count": r_count,
                "product_count": p_count,
                "difference": diff,
                "side": "products" if diff > 0 else "reactants",
                "detail": f"{elem}: reactants={r_count}, products={p_count} (missing {abs(diff)} on {'products' if diff > 0 else 'reactants'} side)"
            })

    total_r = sum(reactant_atoms.values())
    total_p = sum(product_atoms.values())

    result = {
        "valid": len(imbalanced) == 0,
        "reactant_atoms": dict(sorted(reactant_atoms.items())),
        "product_atoms": dict(sorted(product_atoms.items())),
        "imbalanced": imbalanced,
        "reactant_molecules": [{k: v for k, v in m.items() if k != "mol"} for m in reactant_mols],
        "product_molecules": [{k: v for k, v in m.items() if k != "mol"} for m in product_mols],
        "total_reactant_atoms": total_r,
        "total_product_atoms": total_p,
        "reactant_formula": ".".join(reactant_formulas),
        "product_formula": ".".join(product_formulas),
    }

    # Add bond change analysis if requested and atoms are balanced
    if check_mapping:
        if imbalanced:
            result["bond_changes"] = {
                "cleaved": [], "formed": [], "changed": [],
                "mapped_atoms": 0, "total_atoms": 0,
                "mapping_coverage": 0,
                "error": "Atom counts are not balanced — fix atom imbalance before checking bond changes"
            }
        else:
            result["bond_changes"] = analyze_bond_changes(reactant_mols, product_mols)

    return result


def main():
    parser = argparse.ArgumentParser(description="Validate reaction atom balance and check bond changes")
    parser.add_argument("reactants", nargs="?", help="Reactants SMILES (dot-separated)")
    parser.add_argument("products", nargs="?", help="Products SMILES (dot-separated)")
    parser.add_argument("--reactants", dest="reactants_flag", help="Reactants SMILES (flag form)")
    parser.add_argument("--products", dest="products_flag", help="Products SMILES (flag form)")
    parser.add_argument("--check-mapping", action="store_true", help="Also detect bond changes (cleaved/formed bonds)")
    args = parser.parse_args()

    reactants = args.reactants_flag or args.reactants or ""
    products = args.products_flag or args.products or ""

    if not reactants and not products:
        print(json.dumps({"error": "Both reactants and products SMILES are required"}))
        sys.exit(1)

    result = validate_reaction(reactants, products, check_mapping=args.check_mapping)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
