#!/usr/bin/env python3
"""Ligand comparison module: generate SMILES from PDB ligand name, compute MCS, produce highlight SVGs.

Usage:
    python3 ligand_compare.py compare <smiles_a> <smiles_b> [--width 350] [--height 300]
    python3 ligand_compare.py pdb_smiles <res_name>
    python3 ligand_compare.py similarity <smiles_a> <smiles_b>
    python3 ligand_compare.py extract_smiles <pdb_text> <res_name> [chain]

Output: JSON to stdout
"""

import json
import sys
import urllib.request
import urllib.error
from typing import Dict, List, Optional, Tuple
from collections import defaultdict
from math import sqrt

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdFMCS, rdmolfiles, rdmolops
from rdkit.Chem.Draw import MolDraw2DSVG

RDLogger.DisableLog('rdApp.*')


# Common PDB ligand res_name → SMILES mapping (small dictionary of well-known ligands)
_PDB_LIGAND_SMILES: Dict[str, str] = {
    "ATP": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O)N",
    "ADP": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O)N",
    "AMP": "c1nc(c2c(n1)n(cn2)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N",
    "NAD": "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@H]2OC(n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2OP(=O)(O)O)[C@@H](O)[C@H]1O",
    "NAP": "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(=O)([O-])OP(=O)(O)OC[C@H]2OC(n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",
    "FAD": "Cc1cc2nc3c(=O)[nH]c(=O)nc3n(C[C@H](O)[C@H](O)[C@H](O)COP(=O)(O)OP(=O)(O)OC[C@H]3OC(n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3OP(=O)(O)O)c2cc1C",
    "GTP": "c1nc2c(n1)n(c(=O)[nH]c2=N)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OP(=O)(O)O)O)O",
    "GDP": "c1nc2c(n1)n(c(=O)[nH]c2=N)[C@@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)O)O)O",
    "GMP": "c1nc2c(n1)n(c(=O)[nH]c2=N)[C@@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O",
    "HEME": "CC1=C(CCC2=C1C=C(C=C2)CC3=[N]4C=C5C=C(C=C5[N]=4)[N-]3)C",
    "HEM": "CC1=C(CCC2=C1C=C(C=C2)CC3=[N]4C=C5C=C(C=C5[N]=4)[N-]3)C",
    "FE": "[Fe]",
    "ZN": "[Zn]",
    "MG": "[Mg]",
    "MN": "[Mn]",
    "CA": "[Ca]",
    "CU": "[Cu]",
    "CO": "[Co]",
    "NI": "[Ni]",
    "NA": "[Na]",
    "CL": "[Cl-]",
    "K": "[K]",
    "SO4": "O=S(=O)([O-])[O-]",
    "PO4": "O=P(=O)([O-])[O-]",
    "ACE": "CC(=O)",
    "NH2": "N",
    "TRP": "CC(C)(C)C1=CC=C(C=C1)[C@H]2C(=O)NC(=O)[C@@H]2N",
    "AP5": "CC(C)(C)CC(CC(C)(C)C(=O)O)P(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O",
    "EOH": "CCO",
    "NAG": "CC(N)C(O)C1OC(O)C(O)C(O)C1O",
    "MAN": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "GAL": "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",
    "FUC": "CC(O)C(O)C(O)C1OC(O)C(O)C1O",
    "BMA": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "GLC": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "SIA": "CC(=O)OC(C)(C)C(=O)O",
    "NDG": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "PAR": "C(=O)C(=O)N[C@@H](CC(=O)O)C(=O)O",
    "ARA": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "XYP": "C1OC(CO)C(O)C(O)C1O",
    "RIB": "OC[C@H]1OC(O)[C@H](O)[C@@H]1O",
    "DMS": "CS(=O)C",
    "EDO": "CCO",
    "PEG": "OCCO",
    "GOL": "OCC(O)CO",
    "IPA": "CC(O)C",
    "BU2": "CCCC",
    "PGE": r"OC(=O)CC[C@@H]1C[C@H]2[C@@H]3C(=O)/C(=C\C4=CC(=O)CC[C@@]4(C)[C@H]3CC[C@]12C)C(=O)O",
    "LDA": "CC(C)C1=CC=CC=C1C(C)C",
    "MLI": "CC(C)C1=CC=CC=C1C(C)C",
    "CIT": "CC(O)(CC(=O)[O-])C(=O)[O-]",
    "FMT": "C=O",
    "FOR": "C(=O)[O-]",
    "ACT": "CC(=O)C",
    "SAM": "C[S+](CC[C@@H](N)C(=O)O)C",
    "SAH": "C[S+](CC[C@@H](N)C(=O)O)CC1OC(n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",
    "BIOT": "CC(C)(C)[C@H](N)C(=O)N[C@@H]1CSSCC[C@@H](NC(=O)CC2CCCCC2)C1",
    "THA": "CCCCCCCCCCCCCCO",
    "MYR": "CCCCCCCCCCCCCC(=O)O",
    "PLM": "CCCCCCCCCCCCCCCCCCCCCC(=O)O",
    "STE": "CCCCCCCCCCCCCCCC(=O)O",
    "OLE": r"CCCCCCCC/C=C\CCCCCCCC(=O)O",
    "CHD": "OC[C@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",
    "UDP": "C1OC(COP(=O)(O)OP(=O)(O)OC[C@H]2OC(n3cnc4c(O)ncnc43)[C@H](O)[C@@H]2O)C(O)C(O)C1O",
    "CMP": "C1OC(COP(=O)(O)O)C(O)C1O[n]2cnc3c(N)ncnc32",
    "UMP": "C1OC(COP(=O)(O)O)C(O)C1Oc2nc(C)c3c(O)ncnc23",
    "TMP": "C1OC(COP(=O)(O)O)C(O)C1OC2=NC3=C(N)N=C(NC3=N2)N",
    # Additional common ligands
    "PO3": "O=P(=O)([O-])[O-]",
    "CO3": "C(=O)([O-])[O-]",
    "IOD": "[I-]",
    "BR": "[Br-]",
    "F": "[F-]",
    "IMD": "C1=CN=C[NH]1",
    "CYS": "N[C@@H](CS)C(=O)O",
    "SER": "N[C@@H](CO)C(=O)O",
    "THR": "C[C@@H](O)[C@H](N)C(=O)O",
    "LYS": "NCCCC[C@@H](N)C(=O)O",
    "ARG": "N[C@@H](CCCNC(N)=N)C(=O)O",
    "HIS": "N[C@@H](Cc1cnc[nH]1)C(=O)O",
    "ASP": "N[C@@H](CC(=O)O)C(=O)O",
    "GLU": "N[C@@H](CCC(=O)O)C(=O)O",
    "ASN": "N[C@@H](CC(=O)N)C(=O)O",
    "GLN": "N[C@@H](CCC(=O)N)C(=O)O",
    "PHE": "N[C@@H](Cc1ccccc1)C(=O)O",
    "TYR": "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "TRP": "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",
    "MET": "N[C@@H](CCSC)C(=O)O",
    "PRO": "OC(=O)[C@@H]1CCCN1",
    "VAL": "N[C@@H](C(C)C)C(=O)O",
    "LEU": "N[C@@H](CC(C)C)C(=O)O",
    "ILE": "N[C@@H]([C@H](C)CC)C(=O)O",
    "ALA": "N[C@@H](C)C(=O)O",
    "GLY": "NCC(=O)O",
    "MES": "CC(S(=O)(=O)O)C",
    "HEPES": "CC(S(=O)(=O)O)C(CO)NCCO",
    "TRS": "N[C@@H](CC(=O)O)CO",
    "BME": "CSCCO",
    "DTT": "C(CS)CS",
    "AZI": "N=[N+]=[N-]",
    "CYT": "C1OC(COP(=O)(O)O)C(O)C1O[n]1cnc2c(N)ncnc21",
    "URA": "C1OC(COP(=O)(O)O)C(O)C1Oc2nc(C)cc(=O)[nH]2",
    "GUA": "C1OC(COP(=O)(O)O)C(O)C1Oc2nc3c(N)ncnc3n2",
    "ADE": "C1OC(COP(=O)(O)O)C(O)C1Oc2ncnc3c2N=CN=C3N",
    "INH": "NN=C1C=CC=CN1",
    "IPA": "CC(O)C",
    "TAR": "O[C@@H](C(=O)O)[C@H](O)C(=O)O",
    "SUC": "C(C(=O)O)C(=O)O",
    "MLA": "CC(C)C1=CC=CC=C1C(C)C",
    "PLP": "CC1=CN(C=C2C(=NC(=N)N)C(=O)N2C1)CO",
    "THM": "CC1CN1",
    "AMU": "CC(=O)NC1CCC(NC(=O)C)CC1",
    "DRN": "NC1=CC=C(C=C1)C(=O)O",
    "ETH": "CC",
    "BTB": "CC(C)(C)C(=O)O",
    "CLF": "Cl",
    "EDO": "CCO",
    "ENO": "CC(N)=O",
    "FUC": "CC(O)C(O)C(O)C1OC(O)C(O)C1O",
    "FUL": "C1=CC=CC=C1",
    "GAL": "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",
    "MTE": "CSCC(C(=O)O)N",
    "NO3": "O=[N+]([O-])[O-]",
    "OCG": "C1OC(n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",
    "P6G": "C1OC(COP(=O)(O)OP(=O)(O)OC[C@H]2OC(n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)C(O)C(O)C1O",
}

# Water residue names to skip
_WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}


def get_ligand_smiles(res_name: str) -> Optional[str]:
    """Get SMILES for a PDB ligand residue name.

    Uses a built-in dictionary for common ligands.
    Returns None if not found.
    """
    return _PDB_LIGAND_SMILES.get(res_name.upper())


def _lookup_pdbechem(res_name: str, timeout_sec: int = 5) -> Optional[str]:
    """Look up SMILES from PDBeChem API (https://www.ebi.ac.uk/pdbe/chem APIs).

    PDBeChem maintains a database of all chemical components found in PDB structures.
    This is a reliable fallback for ligands not in our built-in dictionary.

    Args:
        res_name: 3-letter PDB residue name
        timeout_sec: request timeout in seconds

    Returns:
        SMILES string or None
    """
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/compound/summary/{res_name}"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout_sec) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            # Navigate: {res_name: [{... compound ...}]}
            compounds = data.get(res_name, [])
            if compounds:
                compound = compounds[0]

                # SMILES field can be a list of dicts: [{"program": ..., "name": "CCO"}, ...]
                smiles_field = compound.get("smiles")
                if isinstance(smiles_field, list) and len(smiles_field) > 0:
                    # Prefer CACTVS or OpenEye, otherwise take the first one
                    for entry in smiles_field:
                        prog = entry.get("program", "").upper()
                        if "CACTVS" in prog:
                            return entry.get("name")
                    # Fallback to first entry
                    return smiles_field[0].get("name")
                elif isinstance(smiles_field, str) and len(smiles_field) > 0:
                    return smiles_field

                # Try InChI -> SMILES conversion
                inchi = compound.get("inchi") or compound.get("standard_inchi")
                if inchi and isinstance(inchi, str):
                    mol = Chem.MolFromInchi(inchi)
                    if mol:
                        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        pass
    return None


def get_ligand_smiles_with_fallback(res_name: str) -> dict:
    """Get SMILES for a ligand with dictionary + PDBeChem API fallback.

    Returns dict with:
        - res_name: the residue name
        - smiles: SMILES string or None
        - source: "dictionary" | "pdbechem" | "water_default" | None
    """
    # 0. Water molecules: return "O" directly — no API call needed
    if res_name.upper() in _WATER_RESNAMES:
        return {"res_name": res_name, "smiles": "O", "source": "water_default"}

    # 1. Try built-in dictionary
    dict_smiles = get_ligand_smiles(res_name)
    if dict_smiles:
        return {"res_name": res_name, "smiles": dict_smiles, "source": "dictionary"}

    # 2. Try PDBeChem API
    pdbe_smiles = _lookup_pdbechem(res_name)
    if pdbe_smiles:
        return {"res_name": res_name, "smiles": pdbe_smiles, "source": "pdbechem"}

    return {"res_name": res_name, "smiles": None, "source": None}


def extract_smiles_from_pdb(pdb_text: str, res_name: str, chain: Optional[str] = None, res_num: Optional[int] = None) -> Optional[str]:
    """Extract SMILES for a specific ligand from PDB structure coordinates.

    Uses RDKit's MolFromPDBBlock to parse the ligand atoms and infer bonds.
    Falls back to CONECT records if available.

    Args:
        pdb_text: Full PDB file text
        res_name: Ligand residue name (3-letter code)
        chain: Optional chain filter
        res_num: Optional residue number filter (critical when multiple instances
                 of the same ligand exist on the same chain)

    Returns:
        Canonical SMILES or None if extraction fails
    """
    if not pdb_text or not res_name:
        return None

    res_name = res_name.upper().strip()
    chain = chain.upper().strip() if chain else None

    # Extract HETATM records for the target ligand
    hetatm_lines = []
    conect_map: Dict[str, List[str]] = defaultdict(list)

    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec == "HETATM":
            rn = line[17:20].strip().upper()
            cid = line[21].strip() if len(line) > 21 else " "
            if rn == res_name and (chain is None or cid == chain):
                # Filter by res_num if specified — prevents merging atoms
                # from multiple instances of the same ligand on one chain
                if res_num is not None:
                    try:
                        line_res_num = int(line[22:26].strip())
                    except (ValueError, IndexError):
                        continue
                    if line_res_num != res_num:
                        continue
                hetatm_lines.append(line)
        elif rec == "CONECT":
            # Parse CONECT records for bond information
            parts = line.split()
            if len(parts) >= 3:
                atom1 = parts[1]
                for atom2 in parts[2:]:
                    conect_map[atom1].append(atom2)

    if not hetatm_lines:
        return None

    # Build a PDB block for just this ligand
    pdb_block_lines = []

    # Renumber atoms sequentially for RDKit
    atom_serial_map: Dict[str, int] = {}
    for idx, line in enumerate(hetatm_lines, start=1):
        orig_serial = line[6:11].strip()
        atom_serial_map[orig_serial] = idx
        new_line = line[:6] + f"{idx:>5}" + line[11:]
        pdb_block_lines.append(new_line)

    # Add CONECT records with renumbered atoms
    seen_conects = set()
    for line in pdb_text.splitlines():
        if line[:6].strip() == "CONECT":
            parts = line.split()
            if len(parts) >= 3:
                orig_a1 = parts[1]
                if orig_a1 in atom_serial_map:
                    for orig_a2 in parts[2:]:
                        if orig_a2 in atom_serial_map:
                            a1 = atom_serial_map[orig_a1]
                            a2 = atom_serial_map[orig_a2]
                            key = (min(a1, a2), max(a1, a2))
                            if key not in seen_conects:
                                seen_conects.add(key)
                                pdb_block_lines.append(f"CONECT{a1:>5}{a2:>5}")

    pdb_block = "\n".join(pdb_block_lines) + "\nEND\n"

    # Try to parse with RDKit
    mol = Chem.MolFromPDBBlock(pdb_block, removeHs=False, sanitize=False)
    if mol is None:
        return None

    # Try to sanitize - if this fails, the molecule is too broken
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        try:
            # Try with less strict sanitization
            Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        except Exception:
            # Last resort: try kekulization only
            try:
                rdmolops.Kekulize(mol, clearAromaticFlags=True)
            except Exception:
                return None

    # Remove hydrogens for cleaner SMILES
    mol = Chem.RemoveHs(mol, sanitize=False)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    # Get canonical SMILES
    try:
        smiles = Chem.MolToSmiles(mol, canonical=True)
        if smiles and len(smiles) >= 1:
            return smiles
    except Exception:
        pass

    return None


def compute_ligand_distances(
    ligands: List[dict],
    active_site_residues: List[dict],
    chain_residues_atoms: dict,
) -> List[float]:
    """Compute minimum distance from each ligand to the nearest selected residue.

    For each ligand, finds the closest atom to any selected residue's CA atom
    and returns the Euclidean distance.

    Args:
        ligands: List of ligand dicts with x, y, z fields
        active_site_residues: List of selected residues with chain, res_num
        chain_residues_atoms: dict of chain -> residue_key -> atom list

    Returns:
        List of distances (same order as ligands), 999.0 if no residues
    """
    if not active_site_residues or not ligands:
        return [999.0] * len(ligands)

    # Collect CA atom positions for selected residues
    residue_positions: List[Tuple[float, float, float]] = []
    for res in active_site_residues:
        cid = res.get("chain", "")
        rnum = res.get("res_num", 0)
        chain_data = chain_residues_atoms.get(cid, {})
        for rkey, atoms in chain_data.items():
            if atoms and atoms[0].get("res_num") == rnum:
                for atom in atoms:
                    if atom.get("atom_name") == "CA":
                        pos = (atom.get("x", 0), atom.get("y", 0), atom.get("z", 0))
                        if all(p != 0 for p in pos):
                            residue_positions.append(pos)
                        break

    if not residue_positions:
        return [999.0] * len(ligands)

    distances = []
    for lig in ligands:
        lx, ly, lz = lig.get("x", 0), lig.get("y", 0), lig.get("z", 0)
        if lx == 0 and ly == 0 and lz == 0:
            distances.append(999.0)
            continue

        min_dist = float("inf")
        for rx, ry, rz in residue_positions:
            d = sqrt((lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2)
            if d < min_dist:
                min_dist = d
        distances.append(round(min_dist, 2))

    return distances


def _find_mcs(smiles_a: str, smiles_b: str, timeout: int = 2) -> dict:
    """Core MCS computation shared by similarity and compare endpoints.

    This is the single source of truth for MCS parameters (atomCompare, bondCompare, timeout).
    Both compute_similarity_only and compute_mcs_highlight delegate to this function,
    ensuring algorithmic consistency and unified timeout behavior.

    Args:
        smiles_a: First molecule SMILES
        smiles_b: Second molecule SMILES
        timeout: MCS search timeout in seconds (default 2, matching original compare_mols)

    Returns:
        Dict with keys:
            - mcs_num_atoms: int (0 if failed)
            - mcs_num_bonds: int (0 if failed)
            - mcs_smarts: str ("" if failed)
            - similarity: float (0.0 if failed)
            - mol_a: RDKit Mol or None
            - mol_b: RDKit Mol or None
    """
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    if mol_a is None or mol_b is None:
        return {
            "mcs_num_atoms": 0, "mcs_num_bonds": 0, "mcs_smarts": "",
            "similarity": 0.0, "mol_a": mol_a, "mol_b": mol_b,
        }
    mx = max(mol_a.GetNumAtoms(), mol_b.GetNumAtoms())
    try:
        mcs = rdFMCS.FindMCS(
            [mol_a, mol_b],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            timeout=timeout,
        )
        return {
            "mcs_num_atoms": mcs.numAtoms,
            "mcs_num_bonds": mcs.numBonds,
            "mcs_smarts": mcs.smartsString,
            "similarity": mcs.numAtoms / mx if mx > 0 else 0.0,
            "mol_a": mol_a, "mol_b": mol_b,
        }
    except Exception:
        return {
            "mcs_num_atoms": 0, "mcs_num_bonds": 0, "mcs_smarts": "",
            "similarity": 0.0, "mol_a": mol_a, "mol_b": mol_b,
        }


def compute_mcs_highlight(
    smiles_a: str,
    smiles_b: str,
    width: int = 350,
    height: int = 300,
) -> dict:
    """Compute MCS between two molecules and return SVGs with green highlights.

    Delegates MCS computation to _find_mcs() for consistency with compute_similarity_only.
    """
    result = _find_mcs(smiles_a, smiles_b, timeout=2)
    mol_a = result["mol_a"]
    mol_b = result["mol_b"]

    if mol_a is None:
        return {"error": f"Invalid SMILES A: {smiles_a}"}
    if mol_b is None:
        return {"error": f"Invalid SMILES B: {smiles_b}"}

    mcs_smarts = result["mcs_smarts"]
    mcs_num_atoms = result["mcs_num_atoms"]
    mcs_num_bonds = result["mcs_num_bonds"]
    similarity = result["similarity"]

    AllChem.Compute2DCoords(mol_a)
    AllChem.Compute2DCoords(mol_b)

    match_a = []
    match_b = []
    if mcs_num_atoms > 0:
        mcs_mol = Chem.MolFromSmarts(mcs_smarts)
        if mcs_mol:
            try:
                match_a = list(mol_a.GetSubstructMatch(mcs_mol))
            except Exception:
                match_a = []
            try:
                match_b = list(mol_b.GetSubstructMatch(mcs_mol))
            except Exception:
                match_b = []

    GREEN = (0.0, 0.75, 0.2, 0.5)

    def draw_highlighted(mol, match_atoms, w, h):
        drawer = MolDraw2DSVG(w, h)
        drawer.SetFontSize(0.75)
        highlight_colors = {a: GREEN for a in match_atoms}
        if match_atoms:
            drawer.DrawMolecule(mol, highlightAtoms=match_atoms, highlightAtomColors=highlight_colors)
        else:
            drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    svg_a = draw_highlighted(mol_a, match_a, width, height)
    svg_b = draw_highlighted(mol_b, match_b, width, height)

    return {
        "mcs_smarts": mcs_smarts,
        "mcs_num_atoms": mcs_num_atoms,
        "mcs_num_bonds": mcs_num_bonds,
        "similarity": round(similarity, 4),
        "svg_a": svg_a,
        "svg_b": svg_b,
        "smiles_a": smiles_a,
        "smiles_b": smiles_b,
        "num_atoms_a": mol_a.GetNumAtoms(),
        "num_atoms_b": mol_b.GetNumAtoms(),
        "match_a": match_a,
        "match_b": match_b,
    }


def compute_similarity_only(smiles_a: str, smiles_b: str) -> dict:
    """Quick similarity score between two SMILES using MCS atom ratio.

    Delegates to _find_mcs() to ensure algorithmic consistency with compute_mcs_highlight
    and unified timeout (2s, matching original compare_mols).

    Returns:
        Dict with keys: similarity, mcs_num_atoms, mcs_num_bonds, mcs_smarts, num_atoms_a, num_atoms_b
    """
    result = _find_mcs(smiles_a, smiles_b, timeout=2)
    return {
        "similarity": round(result["similarity"], 4),
        "mcs_num_atoms": result["mcs_num_atoms"],
        "mcs_num_bonds": result["mcs_num_bonds"],
        "mcs_smarts": result["mcs_smarts"],
        "num_atoms_a": result["mol_a"].GetNumAtoms() if result["mol_a"] else 0,
        "num_atoms_b": result["mol_b"].GetNumAtoms() if result["mol_b"] else 0,
    }


def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: ligand_compare.py compare|pdb_smiles|similarity|extract_smiles ..."}))
        sys.exit(1)

    action = sys.argv[1]

    if action == "compare":
        if len(sys.argv) < 4:
            print(json.dumps({"error": "Usage: ligand_compare.py compare <smiles_a> <smiles_b> [--width 350] [--height 300]"}))
            sys.exit(1)
        smi_a = sys.argv[2]
        smi_b = sys.argv[3]
        w, h = 350, 300
        for i, arg in enumerate(sys.argv):
            if arg == "--width" and i + 1 < len(sys.argv):
                w = int(sys.argv[i + 1])
            if arg == "--height" and i + 1 < len(sys.argv):
                h = int(sys.argv[i + 1])
        result = compute_mcs_highlight(smi_a, smi_b, w, h)
        print(json.dumps(result))

    elif action == "pdb_smiles":
        if len(sys.argv) < 3:
            print(json.dumps({"error": "Usage: ligand_compare.py pdb_smiles <res_name>"}))
            sys.exit(1)
        res_name = sys.argv[2].strip().upper()
        result = get_ligand_smiles_with_fallback(res_name)
        if result["smiles"]:
            print(json.dumps(result))
        else:
            print(json.dumps({"res_name": res_name, "smiles": None, "source": None, "error": f"No SMILES found for {res_name}"}))

    elif action == "similarity":
        if len(sys.argv) < 4:
            print(json.dumps({"error": "Usage: ligand_compare.py similarity <smiles_a> <smiles_b>"}))
            sys.exit(1)
        result = compute_similarity_only(sys.argv[2], sys.argv[3])
        print(json.dumps(result))

    elif action == "extract_smiles":
        """Extract SMILES from PDB structure for a specific ligand.
        
        Usage: ligand_compare.py extract_smiles <res_name> [chain] [res_num]
        
        Reads PDB text from stdin or from file. First checks built-in dictionary,
        then tries to extract SMILES from the actual PDB coordinates using RDKit.
        """
        if len(sys.argv) < 3:
            print(json.dumps({"error": "Usage: ligand_compare.py extract_smiles <res_name> [chain] [res_num]  (reads PDB from stdin)"}))
            sys.exit(1)

        res_name = sys.argv[2].strip().upper()
        chain = sys.argv[3].strip() if len(sys.argv) > 3 else None
        res_num = int(sys.argv[4].strip()) if len(sys.argv) > 4 else None

        # First try the built-in dictionary
        dict_smiles = get_ligand_smiles(res_name)

        # Read PDB text from stdin
        pdb_text = sys.stdin.read() if not sys.stdin.isatty() else ""

        if not pdb_text:
            # No PDB text provided, just return dictionary result
            if dict_smiles:
                print(json.dumps({"res_name": res_name, "smiles": dict_smiles, "source": "dictionary"}))
            else:
                print(json.dumps({"res_name": res_name, "smiles": None, "error": f"No SMILES for {res_name} (no PDB text provided)"}))
            sys.exit(0)

        if dict_smiles:
            print(json.dumps({"res_name": res_name, "smiles": dict_smiles, "source": "dictionary"}))
            sys.exit(0)

        # Try PDBeChem API
        pdbe_smiles = _lookup_pdbechem(res_name)
        if pdbe_smiles:
            print(json.dumps({"res_name": res_name, "smiles": pdbe_smiles, "source": "pdbechem"}))
            sys.exit(0)

        # Try extracting from PDB structure
        extracted = extract_smiles_from_pdb(pdb_text, res_name, chain, res_num)
        if extracted:
            print(json.dumps({"res_name": res_name, "smiles": extracted, "source": "pdb_structure"}))
        else:
            print(json.dumps({"res_name": res_name, "smiles": None, "source": None, "error": f"Could not find SMILES for {res_name} (tried dictionary, PDBeChem API, and PDB structure extraction)"}))

    else:
        print(json.dumps({"error": f"Unknown action: {action}"}))
        sys.exit(1)


if __name__ == "__main__":
    main()
