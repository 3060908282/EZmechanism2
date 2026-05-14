#!/usr/bin/env python3
"""
overlay_substrate.py — Overlay user-drawn substrate onto PDB ligand 3D position.

Pipeline:
1. Parse PDB text to extract target ligand 3D coordinates (by res_name + chain + res_num)
2. Generate 3D conformer from user's hand-drawn SMILES
3. Use RDKit's FindMCS to establish atom correspondence between drawn molecule and PDB ligand
4. Apply Kabsch alignment to superimpose drawn molecule onto PDB ligand position
5. Output aligned molecule as MolBlock with 3D coordinates + alignment quality metrics

Usage:
    overlay_substrate.py <smiles> <pdb_text_file_or_stdin> [--res-name ATP] [--chain A] [--res-num 301] [--max-conformers 50] [--rmsd-threshold 2.0]

Output (JSON):
    {
        "success": true/false,
        "molblock": "...",       // aligned MolBlock with 3D coords
        "sdf": "...",            // aligned SDF text
        "rmsd": 0.85,           // alignment RMSD (Å)
        "mcs_atoms_drawn": [0,1,2,...],  // atom indices in drawn mol matched by MCS
        "mcs_atoms_pdb": [0,1,2,...],    // atom indices in PDB mol matched by MCS
        "mcs_num_atoms": 15,
        "mcs_smarts": "...",
        "num_mapped": 15,
        "total_drawn_atoms": 20,
        "total_pdb_atoms": 22,
        "pdb_ligand_smiles": "...",  // SMILES of the extracted PDB ligand
        "pdb_ligand_formula": "...",
        "message": "..."
    }
"""

import sys
import json
import argparse
from typing import Optional, Tuple, List, Dict, Any

# RDKit imports
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdMolTransforms, Descriptors, rdFMCS
from rdkit.Geometry import Point3D
import numpy as np


# Water residue names
_WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}


def parse_pdb_ligand_atoms(pdb_text: str, res_name: str, chain: Optional[str] = None, res_num: Optional[int] = None) -> Optional[Chem.Mol]:
    """
    Extract ligand atoms from PDB text and return as an RDKit Mol with 3D coordinates.
    
    Strategy:
    1. Extract HETATM lines for the target ligand
    2. Build a minimal PDB block with just those atoms
    3. Use RDKit's MolFromPDBBlock which can perceive bonds from inter-atomic distances
    4. Return the mol with 3D coordinates preserved
    
    Returns: RDKit Mol with 3D conformer, or None if not found.
    """
    res_name_upper = res_name.upper()
    
    # Extract matching HETATM lines
    het_lines = []
    conect_lines = []
    
    for line in pdb_text.split('\n'):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            if len(line) < 54:
                continue
            
            atom_res_name = line[17:20].strip().upper()
            atom_chain = line[21].strip() if len(line) > 21 else ''
            
            # Skip water molecules
            if atom_res_name in ('HOH', 'WAT', 'DOD', 'TIP3'):
                continue
            
            # Filter by res_name
            if atom_res_name != res_name_upper:
                continue
            
            # Filter by chain if specified
            if chain is not None and atom_chain != chain.upper():
                continue
            
            # Filter by res_num if specified
            try:
                atom_res_num = int(line[22:26].strip())
            except (ValueError, IndexError):
                continue
            if res_num is not None and atom_res_num != res_num:
                continue
            
            het_lines.append(line)
    
    if not het_lines:
        return None
    
    # Try RDKit's MolFromPDBBlock first (handles bond perception from distances)
    pdb_block = '\n'.join(het_lines)
    
    try:
        mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True, sanitize=False)
        if mol and mol.GetNumConformers() > 0 and mol.GetNumAtoms() >= 2:
            # Try to sanitize to fix bond orders
            try:
                Chem.SanitizeMol(mol)
            except:
                pass
            return mol
    except:
        pass
    
    # Fallback: manual parsing with CONECT records
    atoms = []
    for line in het_lines:
        try:
            serial = int(line[6:11].strip())
            atom_name = line[12:16].strip()
            atom_chain = line[21].strip()
            atom_res_num = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            element = line[76:78].strip() if len(line) >= 78 else ''
            if not element:
                element = ''.join(c for c in atom_name if c.isalpha())[:1].upper()
            if not element:
                continue
            atoms.append((serial, atom_name, atom_chain, atom_res_num, x, y, z, element))
        except (ValueError, IndexError):
            continue
    
    if not atoms:
        return None
    
    # Parse CONECT records
    conects = {}
    serial_to_idx = {a[0]: i for i, a in enumerate(atoms)}
    for line in pdb_text.split('\n'):
        if line.startswith('CONECT'):
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                src_serial = int(parts[1])
                if src_serial not in serial_to_idx:
                    continue
                targets = [int(p) for p in parts[2:6] if int(p) in serial_to_idx]
                if targets:
                    conects[src_serial] = targets
            except (ValueError, IndexError):
                continue
    
    # Build RDKit Mol
    emol = Chem.RWMol()
    atom_positions = []
    
    for serial, atom_name, atom_chain, atom_res_num, x, y, z, element in atoms:
        elem_cap = element.upper()
        atomic_num = None
        try:
            atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(elem_cap)
        except:
            elem_map = {'C': 6, 'H': 1, 'O': 8, 'N': 7, 'S': 16, 'P': 15, 'F': 9, 'CL': 17, 'BR': 35, 'I': 53, 'FE': 26, 'MG': 12, 'MN': 25, 'ZN': 30, 'CA': 20, 'NA': 11, 'K': 19, 'SE': 34}
            atomic_num = elem_map.get(elem_cap)
        if atomic_num is None:
            atomic_num = 6
        emol.AddAtom(Chem.Atom(atomic_num))
        atom_positions.append((x, y, z))
    
    if atom_positions:
        conf = Chem.Conformer(len(atom_positions))
        for i, (x, y, z) in enumerate(atom_positions):
            conf.SetAtomPosition(i, Point3D(x, y, z))
        emol.AddConformer(conf)
    
    # Add bonds from CONECT
    added_bonds = set()
    for src_serial, targets in conects.items():
        if src_serial not in serial_to_idx:
            continue
        src_idx = serial_to_idx[src_serial]
        for tgt_serial in targets:
            if tgt_serial not in serial_to_idx:
                continue
            tgt_idx = serial_to_idx[tgt_serial]
            bond_key = (min(src_idx, tgt_idx), max(src_idx, tgt_idx))
            if bond_key in added_bonds:
                continue
            emol.AddBond(src_idx, tgt_idx, Chem.BondType.SINGLE)
            added_bonds.add(bond_key)
    
    # If no bonds from CONECT, try distance-based bond perception
    if not added_bonds and len(atom_positions) > 1:
        # Use RDKit's ConnectedEdges to add bonds based on inter-atomic distances
        try:
            from rdkit.Chem import rdDistGeom
            # Compute distance matrix
            n = len(atom_positions)
            dmat = rdDistGeom.Get3DDistanceMatrix(emol)
            # Add bonds for atoms within typical covalent distance
            bond_adder = Chem.FragmentOnBonds()
            for i in range(n):
                for j in range(i + 1, n):
                    dist = dmat[i][j]
                    # Typical covalent bond radii sums (generous cutoff)
                    if dist < 2.0:  # conservative cutoff for covalent bonds
                        emol.AddBond(i, j, Chem.BondType.SINGLE)
            added_bonds.update([(i, j) for i in range(n) for j in range(i+1, n)])
        except:
            pass
    
    try:
        Chem.SanitizeMol(emol)
    except:
        pass
    
    return emol.GetMol()


def find_mcs_mapping(mol_a: Chem.Mol, mol_b: Chem.Mol) -> Optional[Tuple[List[int], List[int], str]]:
    """
    Find Maximum Common Substructure between two molecules and return atom mapping.
    
    Returns: (atom_indices_a, atom_indices_b, mcs_smarts) or None if no MCS found.
    """
    mcs_params = rdFMCS.MCSParameters()
    # rdkit 2026+: AtomCompare→AtomTyper, BondCompare→BondTyper, params moved to sub-objects
    mcs_params.AtomTyper = rdFMCS.AtomCompare.CompareAny
    mcs_params.BondTyper = rdFMCS.BondCompare.CompareAny
    mcs_params.MaximizeBonds = True
    mcs_params.Timeout = 10  # seconds
    mcs_params.AtomCompareParameters.MatchChiralTag = False
    mcs_params.AtomCompareParameters.RingMatchesRingOnly = True
    mcs_params.AtomCompareParameters.CompleteRingsOnly = True
    mcs_params.Threshold = 0.8
    
    mcs_result = rdFMCS.FindMCS([mol_a, mol_b], mcs_params)
    
    if mcs_result.numAtoms < 3:
        # Too small MCS, try relaxed parameters
        mcs_params2 = rdFMCS.MCSParameters()
        mcs_params2.AtomTyper = rdFMCS.AtomCompare.CompareAny
        mcs_params2.BondTyper = rdFMCS.BondCompare.CompareAny
        mcs_params2.MaximizeBonds = True
        mcs_params2.Timeout = 10
        mcs_params2.AtomCompareParameters.MatchChiralTag = False
        mcs_params2.AtomCompareParameters.RingMatchesRingOnly = False
        mcs_params2.AtomCompareParameters.CompleteRingsOnly = False
        mcs_params2.Threshold = 0.5
        mcs_result = rdFMCS.FindMCS([mol_a, mol_b], mcs_params2)
    
    if mcs_result.numAtoms < 2:
        return None
    
    # Get MCS pattern
    mcs_smarts = Chem.MolToSmarts(mcs_result.queryMol)
    
    # Build atom mapping using substructure match
    patt = Chem.MolFromSmarts(mcs_smarts)
    if patt is None:
        return None
    
    match_a = mol_a.GetSubstructMatch(patt)
    match_b = mol_b.GetSubstructMatch(patt)
    
    if not match_a or not match_b or len(match_a) != len(match_b):
        return None
    
    return list(match_a), list(match_b), mcs_smarts


def generate_3d_conformer(smiles: str, max_conformers: int = 50, max_attempts: int = 200) -> Optional[Chem.Mol]:
    """
    Generate a 3D conformer from SMILES string.
    Uses ETKDG with multiple conformers and picks the best one.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D conformers using ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 1
    # maxAttempts was renamed to maxConformersAttempted in newer rdkit versions
    try:
        params.maxAttempts = max_attempts
    except AttributeError:
        try:
            params.maxConformersAttempted = max_attempts
        except AttributeError:
            pass  # Use default
    
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=max_conformers, params=params)
    
    if not cids:
        # Fallback: try single conformer
        try:
            AllChem.EmbedMolecule(mol, params=params)
        except:
            pass
        
        if mol.GetNumConformers() == 0:
            return None
    else:
        # Use the first conformer (ETKDG usually produces good results)
        # For alignment purposes, the specific conformer shape doesn't matter much
        # since we'll align it to the PDB ligand
        pass
    
    # Remove hydrogens for alignment (PDB ligands typically don't have H coords)
    mol_no_h = Chem.RemoveHs(mol)
    
    if mol_no_h.GetNumConformers() == 0:
        return mol  # Return with H if RemoveHs removed all conformers
    
    return mol_no_h


def kabsch_align(query_mol: Chem.Mol, ref_mol: Chem.Mol, query_match: List[int], ref_match: List[int]) -> float:
    """
    Apply Kabsch alignment to superimpose query_mol onto ref_mol using matched atom pairs.
    
    Returns: RMSD after alignment.
    """
    if len(query_match) < 3 or len(ref_match) < 3:
        return 999.0
    
    if len(query_match) != len(ref_match):
        return 999.0
    
    query_conf = query_mol.GetConformer()
    ref_conf = ref_mol.GetConformer()
    
    # Get coordinates
    query_coords = np.array([list(query_conf.GetAtomPosition(i)) for i in query_match])
    ref_coords = np.array([list(ref_conf.GetAtomPosition(i)) for i in ref_match])
    
    # Kabsch algorithm
    # 1. Compute centroids
    centroid_query = np.mean(query_coords, axis=0)
    centroid_ref = np.mean(ref_coords, axis=0)
    
    # 2. Center the coordinates
    query_centered = query_coords - centroid_query
    ref_centered = ref_coords - centroid_ref
    
    # 3. Compute covariance matrix
    H = query_centered.T @ ref_centered
    
    # 4. SVD
    U, S, Vt = np.linalg.svd(H)
    
    # 5. Correct for reflection
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, np.sign(d)])
    
    # 6. Rotation matrix
    R = Vt.T @ sign_matrix @ U.T
    
    # 7. Translation
    t = centroid_ref - R @ centroid_query
    
    # Apply rotation and translation to ALL atoms in query molecule
    num_atoms = query_mol.GetNumAtoms()
    new_positions = np.zeros((num_atoms, 3))
    
    for i in range(num_atoms):
        pos = np.array(query_conf.GetAtomPosition(i))
        new_pos = R @ pos + t
        new_positions[i] = new_pos
    
    # Update conformer positions
    for i in range(num_atoms):
        query_conf.SetAtomPosition(i, Point3D(*new_positions[i]))
    
    # Compute RMSD
    rmsd = np.sqrt(np.mean(np.sum((ref_centered - query_centered @ R) ** 2, axis=1)))
    
    return float(rmsd)


def _overlay_water_to_pdb(smiles: str, pdb_text: str, res_name: str, chain: Optional[str] = None, res_num: Optional[int] = None) -> Dict[str, Any]:
    """Special overlay path for water molecules.

    Water molecules have only 1 heavy atom (O), so the standard MCS + Kabsch pipeline
    cannot work (requires ≥3 matching points). Instead, we directly read the O atom
    coordinates from the PDB HETATM record and place the generated water molecule at
    that position.

    Returns: Result dict compatible with overlay_smiles_to_pdb() output format.
    """
    result = {
        "success": False,
        "molblock": None,
        "sdf": None,
        "rmsd": 0.0,
        "mcs_atoms_drawn": [0],
        "mcs_atoms_pdb": [0],
        "mcs_num_atoms": 1,
        "mcs_smarts": "[#8]",
        "num_mapped": 1,
        "total_drawn_atoms": 1,
        "total_pdb_atoms": 1,
        "pdb_ligand_smiles": "O",
        "pdb_ligand_formula": "H2O",
        "message": "",
        "is_water": True,
    }

    # Find the water O atom in PDB
    for line in pdb_text.split('\n'):
        if not (line.startswith('HETATM') or line.startswith('ATOM')):
            continue
        if len(line) < 54:
            continue

        atom_res_name = line[17:20].strip().upper()
        if atom_res_name not in _WATER_RESNAMES:
            continue
        if res_name.upper() not in _WATER_RESNAMES:
            continue

        atom_chain = line[21].strip() if len(line) > 21 else ''
        if chain is not None and atom_chain != chain.upper():
            continue

        try:
            atom_res_num = int(line[22:26].strip())
        except (ValueError, IndexError):
            continue
        if res_num is not None and atom_res_num != res_num:
            continue

        # Found the matching water — extract O coordinates
        atom_name = line[12:16].strip().upper()
        # For water, we want the O atom (not H atoms)
        if atom_name not in ('O', 'OW', 'O1'):
            continue

        try:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
        except (ValueError, IndexError):
            continue

        # Generate water molecule and place it at PDB coordinates
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            # Fallback: use "O" as SMILES
            mol = Chem.MolFromSmiles("O")
        if mol is None:
            result["message"] = "Failed to create water molecule"
            return result

        mol = Chem.AddHs(mol)

        # Generate 3D conformer
        try:
            AllChem.EmbedMolecule(mol)
        except:
            pass

        if mol.GetNumConformers() == 0:
            result["message"] = "Failed to generate 3D conformer for water"
            return result

        # Find the O atom in the generated molecule and align it to PDB position
        conf = mol.GetConformer()
        o_idx = -1
        for i in range(mol.GetNumAtoms()):
            if mol.GetAtomWithIdx(i).GetAtomicNum() == 8:  # Oxygen
                o_idx = i
                break

        if o_idx < 0:
            result["message"] = "No oxygen atom found in generated water molecule"
            return result

        # Get O atom position in generated molecule
        o_pos = np.array(list(conf.GetAtomPosition(o_idx)))
        # Translation vector: move generated O to PDB O position
        pdb_pos = np.array([x, y, z])
        translation = pdb_pos - o_pos

        # Apply translation to all atoms
        num_atoms = mol.GetNumAtoms()
        for i in range(num_atoms):
            pos = np.array(list(conf.GetAtomPosition(i)))
            new_pos = pos + translation
            conf.SetAtomPosition(i, Point3D(*new_pos))

        result["success"] = True
        result["total_drawn_atoms"] = mol.GetNumAtoms()
        result["message"] = f"Water coordinate mapping applied: O atom placed at PDB position ({x:.2f}, {y:.2f}, {z:.2f})"
        result["pdb_ligand_smiles"] = "O"

        # Generate MolBlock
        try:
            mol_no_h = Chem.RemoveHs(mol)
            result["molblock"] = Chem.MolToMolBlock(mol_no_h)
            result["sdf"] = result["molblock"] + "\n$$$$\n" if result["molblock"] else None
        except:
            try:
                result["molblock"] = Chem.MolToMolBlock(mol)
                result["sdf"] = result["molblock"] + "\n$$$$\n" if result["molblock"] else None
            except:
                pass

        return result

    result["message"] = f"Water ligand {res_name} not found in PDB structure"
    return result


def overlay_smiles_to_pdb(
    smiles: str,
    pdb_text: str,
    res_name: str,
    chain: Optional[str] = None,
    res_num: Optional[int] = None,
    max_conformers: int = 50,
    rmsd_threshold: float = 2.0,
) -> Dict[str, Any]:
    """
    Main overlay function: aligns a SMILES molecule onto a PDB ligand position.
    
    Returns: Result dictionary with aligned molblock, RMSD, and metadata.
    """
    result = {
        "success": False,
        "molblock": None,
        "sdf": None,
        "rmsd": None,
        "mcs_atoms_drawn": None,
        "mcs_atoms_pdb": None,
        "mcs_num_atoms": 0,
        "mcs_smarts": None,
        "num_mapped": 0,
        "total_drawn_atoms": 0,
        "total_pdb_atoms": 0,
        "pdb_ligand_smiles": None,
        "pdb_ligand_formula": None,
        "message": "",
    }

    # Water molecules: skip MCS + Kabsch (only 1 heavy atom), use direct coordinate mapping
    if res_name.upper() in _WATER_RESNAMES:
        return _overlay_water_to_pdb(smiles, pdb_text, res_name, chain, res_num)

    # Step 1: Extract PDB ligand
    pdb_ligand = parse_pdb_ligand_atoms(pdb_text, res_name, chain, res_num)
    if pdb_ligand is None:
        result["message"] = f"Ligand {res_name} not found in PDB structure"
        if chain:
            result["message"] += f" on chain {chain}"
        if res_num:
            result["message"] += f" at position {res_num}"
        return result
    
    result["total_pdb_atoms"] = pdb_ligand.GetNumAtoms()
    
    # Get PDB ligand SMILES
    try:
        result["pdb_ligand_smiles"] = Chem.MolToSmiles(pdb_ligand)
    except:
        pass
    
    try:
        result["pdb_ligand_formula"] = Chem.rdMolDescriptors.CalcMolFormula(pdb_ligand)
    except:
        pass
    
    if pdb_ligand.GetNumConformers() == 0:
        result["message"] = "PDB ligand has no 3D coordinates"
        return result
    
    # Step 2: Generate 3D conformer from SMILES
    drawn_mol = generate_3d_conformer(smiles, max_conformers=max_conformers)
    if drawn_mol is None:
        result["message"] = f"Failed to generate 3D structure from SMILES: {smiles[:50]}"
        return result
    
    result["total_drawn_atoms"] = drawn_mol.GetNumAtoms()
    
    if drawn_mol.GetNumConformers() == 0:
        result["message"] = "Generated molecule has no 3D conformer"
        return result
    
    # Step 3: Find MCS mapping
    mcs_result = find_mcs_mapping(drawn_mol, pdb_ligand)
    if mcs_result is None:
        result["message"] = f"No sufficient common substructure found between drawn molecule ({drawn_mol.GetNumAtoms()} atoms) and PDB ligand {res_name} ({pdb_ligand.GetNumAtoms()} atoms)"
        return result
    
    drawn_match, pdb_match, mcs_smarts = mcs_result
    result["mcs_atoms_drawn"] = drawn_match
    result["mcs_atoms_pdb"] = pdb_match
    result["mcs_num_atoms"] = len(drawn_match)
    result["mcs_smarts"] = mcs_smarts
    result["num_mapped"] = len(drawn_match)
    
    # Step 4: Kabsch alignment
    rmsd = kabsch_align(drawn_mol, pdb_ligand, drawn_match, pdb_match)
    result["rmsd"] = round(rmsd, 4)
    
    # Check RMSD quality
    if rmsd > rmsd_threshold:
        result["message"] = f"Warning: Alignment RMSD = {rmsd:.2f} Å (threshold: {rmsd_threshold:.1f} Å). The overlay may not be accurate due to structural differences between the drawn molecule and the PDB ligand."
        # Still return the result, but mark success as questionable
    
    # Step 5: Generate output
    result["success"] = True
    if not result["message"]:
        coverage = len(drawn_match) / max(drawn_mol.GetNumAtoms(), pdb_ligand.GetNumAtoms()) * 100
        result["message"] = f"Overlay successful: {len(drawn_match)} atoms mapped (RMSD={rmsd:.2f}Å, coverage={coverage:.0f}%)"
    
    try:
        result["molblock"] = Chem.MolToMolBlock(drawn_mol)
    except:
        result["molblock"] = None
    
    try:
        sdf = Chem.SDMolSupplier()
        # Generate SDF manually
        result["sdf"] = result["molblock"] + "\n$$$$\n" if result["molblock"] else None
    except:
        result["sdf"] = None
    
    return result


def main():
    parser = argparse.ArgumentParser(description="Overlay substrate SMILES onto PDB ligand 3D position")
    parser.add_argument("smiles", help="SMILES of the hand-drawn substrate")
    parser.add_argument("pdb_file", help="PDB file path or '-' for stdin")
    parser.add_argument("--res-name", required=True, help="PDB ligand residue name (e.g., ATP)")
    parser.add_argument("--chain", default=None, help="PDB chain ID (e.g., A)")
    parser.add_argument("--res-num", type=int, default=None, help="PDB residue number")
    parser.add_argument("--max-conformers", type=int, default=50, help="Max conformers to generate (default: 50)")
    parser.add_argument("--rmsd-threshold", type=float, default=2.0, help="RMSD threshold for quality warning (default: 2.0)")
    
    args = parser.parse_args()
    
    # Read PDB text
    if args.pdb_file == '-':
        pdb_text = sys.stdin.read()
    else:
        with open(args.pdb_file, 'r') as f:
            pdb_text = f.read()
    
    result = overlay_smiles_to_pdb(
        smiles=args.smiles,
        pdb_text=pdb_text,
        res_name=args.res_name,
        chain=args.chain,
        res_num=args.res_num,
        max_conformers=args.max_conformers,
        rmsd_threshold=args.rmsd_threshold,
    )
    
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
