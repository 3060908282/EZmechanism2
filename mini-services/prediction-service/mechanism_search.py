#!/usr/bin/env python3
"""
M-CSA Bidirectional Mechanism Search Engine

*** This module is the ACTIVE mechanism search implementation. ***
*** Based on Upload source code (upload_models.py PredictionRun). ***

This file contains:
- bidirectional_search(): Main entry point for mechanism search
- PredictionRunState: State management (equivalent to Upload's PredictionRun)
- Upload-aligned scoring: (max(4, bl)-3)^2, bonds_away_from_rc discount
- frozenset(mols) config dedup (same as Upload source code)

Based on the EzMechanism paper (Nature Methods 2023).

Usage:
    python3 mechanism_search.py --reactants "CCO" --products "CC(=O)O" --max-configs 300 --max-rules 5000

Output: JSON to stdout
"""

import argparse
import gc
import itertools
import json
import logging
import math
import sys
import time
from collections import defaultdict
from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple

logger = logging.getLogger("mechanism_search")

# Ensure logging is configured when run as CLI (Flask subprocess won't inherit basicConfig)
if not logger.handlers and not logging.getLogger().handlers:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')


# ---------------------------------------------------------------------------
# Async progress file support
# ---------------------------------------------------------------------------
def write_progress(progress_file: str, state: str, explored_nodes: int = 0,
                   total_nodes: int = 0, current_iteration: int = 0,
                   elapsed_seconds: float = 0.0, error: str = '',
                   result_file: str = '') -> None:
    """Write progress status to a JSON file for async polling.

    Called from the search loop after each iteration so that
    mechanism-search-status can read the file and report progress to the
    frontend without waiting for the full search to complete.

    Args:
        progress_file: Path to the progress JSON file
        state: One of "RUNNING", "DONE", "ERROR"
        explored_nodes: Number of explored configurations
        total_nodes: Total configurations in the search graph
        current_iteration: Current search iteration number
        elapsed_seconds: Time elapsed since search started
        error: Error message (only when state="ERROR")
        result_file: Path to the result JSON file (only when state="DONE")
    """
    from datetime import datetime, timezone
    payload = {
        "state": state,
        "explored_nodes": explored_nodes,
        "total_nodes": total_nodes,
        "current_iteration": current_iteration,
        "elapsed_seconds": round(elapsed_seconds, 1),
        "updated_at": datetime.now(timezone.utc).isoformat(),
    }
    if error:
        payload["error"] = error
    if result_file:
        payload["result_file"] = result_file
    try:
        with open(progress_file, 'w') as f:
            json.dump(payload, f)
    except Exception as e:
        logger.warning("Failed to write progress file %s: %s", progress_file, e)

# Disable RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from rdkit import Chem
from rdkit.Chem import AllChem

import networkx as nx

# From shared module (replaces local duplicates)
from shared import (
    load_rules as shared_load_rules,
    merge_rule_parts,
    RULES_FILE,
)

# EzMechanism paper's FindMCS-based matching (ported from original codebase)
from predict_mechanism_util import (
    mol_is_water_or_proton,
)

# Reaction result processing (Upload reaction_result.py)
from reaction_result import ReactionResult, optimize_mol_2d, detect_mol_stereo


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_MAX_CONFIGS = 300
DEFAULT_MAX_RULES = 0  # 0 = load all rules

# Batch size for rule processing — prevents OOM when processing 51K+ rules.
# Each batch parses SMARTS into RDKit mol objects, then frees them after use.
RULE_BATCH_SIZE = 2000
# Run gc.collect() every N batches (not every batch — gc.collect is ~0.5s each)
GC_INTERVAL = 10

# 3D geometry / bond length constraints
MAX_BOND_LENGTH = 6.0          # Angstroms: filter out reactions forming bonds > this
FAVORABLE_BOND_LENGTH = 3.0   # Angstroms: reference for favorable geometry

# Residue side chain SMILES for common catalytic residues.
#
# Design principle: "2-bond range" from reaction center (aligned with original
# EzMechanism's side_chain_*.mol files). Each fragment contains the catalytic
# functional group plus 1-2 bonds of context, WITHOUT backbone atoms (N-Cα-C=O)
# or connection point markers (*).
#
# Rationale for NO * connection point:
#   - M-CSA rules use fully-specified SMARTS (e.g., [#6+0:1]-[#8-:2]) with
#     ZERO wildcards ([*], [R], [#0] not found in any of 6083 complete rules).
#   - In RDKit, * (dummy atom, atomicNum=0) CANNOT match [#6+0] (carbon) in
#     rule patterns. Adding * reduces SER matching by 98.5% (465→7 rules)!
#   - The original EzMechanism uses ChemAxon Marvin where R-atoms match C in
#     substructure search, but RDKit's * behaves fundamentally differently.
#
# Rationale for minimal fragments:
#   - Full amino acid SMILES (e.g., N[C@@H](CO)C(=O)O for SER) include the
#     backbone (N-Cα-C=O) which creates many equivalent H atoms (NH₂, Cα-H,
#     carboxyl OH). Each equivalent H multiplies RunReactants products.
#   - Minimal fragments reduce equivalent H atoms by 20-81%, dramatically
#     reducing product explosion in the search graph.
#   - The C/H equivalence principle (M-CSA paper): at 2-bond distance from
#     reaction center, C and H are interchangeable in rules, so explicit H
#     atoms at that distance are unnecessary. SMILES implicit H handles this.
#
# Experimental validation (6083 complete rules tested):
#   SER: CO=465 rules vs full_aa=1742 vs *CO=7 (98.5% loss with *)
#   CYS: CS=261 rules vs full_aa=1568 vs *CS=5  (98% loss with *)
#   LYS: CCN=510 rules vs CCCCN=513 (minimal is fine)
#   ASP: CC(=O)O=950 rules vs full_aa=1468
#   HIS: c1cnc[nH]1=5 rules (all variants same, very specific)
RESIDUE_SIDECHAIN_SMILES: Dict[str, str] = {
    'SER': 'CO',                     # -CH₂-OH (hydroxyl + 1-bond C context)
    'CYS': 'CS',                     # -CH₂-SH (thiol + 1-bond C context)
    'LYS': 'CCN',                    # -CH₂-CH₂-NH₂ (amine + 2-bond context)
    'ARG': 'CNC(=N)N',               # -CH₂-NH-C(=NH)-NH₂ (guanidinium core)
    'HIS': 'c1cnc[nH]1',             # Imidazole ring (catalytic N-H)
    'ASP': 'CC(=O)O',                # -CH₂-COOH (carboxylate + 1-bond context)
    'GLU': 'CC(=O)O',                # Same fragment as ASP (2-bond range equivalent)
    'ASN': 'CC(=O)N',                # -CH₂-CONH₂ (amide + 1-bond context)
    'GLN': 'CCC(=O)N',              # -CH₂-CH₂-CONH₂ (amide + 2-bond context)
    'TYR': 'c1ccc(O)cc1',            # Phenol ring (aromatic OH)
    'THR': 'C(C)O',                  # -CH(OH)-CH₃ (secondary alcohol)
    'MET': 'CCSC',                   # -CH₂-CH₂-S-CH₃ (thioether)
    'TRP': 'c1c[nH]c2ccccc12',      # Indole ring (catalytic N-H)
    'PHE': 'c1ccccc1',              # Phenyl ring (aromatic context)
    'LEU': 'CC(C)C',                # -CH₂-CH(CH₃)₂ (isobutyl)
    'ILE': 'C(C)CC',                # -CH(CH₃)-CH₂-CH₃ (sec-butyl)
    'VAL': 'C(C)C',                 # -CH(CH₃)₂ (isopropyl)
    'PRO': 'C1CCCN1',              # Pyrrolidine ring
    'ALA': 'C',                     # -CH₃ (methyl)
    'GLY': '',                       # No side chain (H only) — skip or use main chain
}

# Main chain fragment SMILES for catalytic backbone groups.
# Used when a residue's 'part' is 'main_chain' instead of 'side_chain'.
# Front-end currently only supports side_chain / main_chain toggle.
MAINCHAIN_FRAGMENT_SMILES: Dict[str, str] = {
    'main_chain': 'NCC(=O)',               # N-Cα-C=O (generic backbone)
    'main_chain_amine': 'NC',              # N-H amine fragment
    'main_chain_carbonyl': 'CC(=O)',       # C=O carbonyl fragment
    'main_chain_n_terminus': 'NCC(=O)O',   # NH₃⁺-Cα-COOH
    'main_chain_c_terminus': 'NCC(=O)O',   # NH₂-Cα-COOH
}


def get_residue_smiles(residues: Optional[List[Dict[str, Any]]]) -> List[str]:
    """Generate residue SMILES from residue info.

    Upload includes catalytic residue side chains as molecules in the initial
    configuration so M-CSA rules (which describe residue+substrate chemistry)
    can substructure-match against them.

    Uses minimal side chain fragments (2-bond range from reaction center)
    instead of full amino acid SMILES to reduce equivalent H atom explosion.

    Args:
        residues: List of dicts with:
            - 'res_name' or 'name' (3-letter code, e.g. 'SER')
            - 'part' (optional: 'side_chain' (default) or 'main_chain')

    Returns:
        List of SMILES strings for catalytic residues.
    """
    if not residues:
        return []
    result = []
    for res in residues:
        name = res.get('res_name', '') or res.get('name', '')
        name = name.upper()
        part = (res.get('part', 'side_chain') or 'side_chain').lower()

        if part == 'main_chain':
            # Main chain catalytic group (backbone N-H, C=O, etc.)
            smi = MAINCHAIN_FRAGMENT_SMILES.get('main_chain', '')
        elif name in RESIDUE_SIDECHAIN_SMILES:
            smi = RESIDUE_SIDECHAIN_SMILES[name]
        else:
            smi = ''

        if smi:  # Skip empty strings (e.g., GLY has no side chain)
            result.append(smi)
    return result

# Cache sizing
# 3D Coordinate & PDB Parsing
# ---------------------------------------------------------------------------
def parse_pdb_coords(pdb_text: str) -> Dict[str, Tuple[float, float, float]]:
    """Parse PDB format text and return atom symbol -> (x, y, z) coordinates.

    Handles ATOM and HETATM records. Returns coordinates keyed by
    element symbol (with occurrence index for duplicates, e.g. C, C1, C2).
    """
    coords: Dict[str, Tuple[float, float, float]] = {}
    element_counts: Dict[str, int] = {}

    for line in pdb_text.strip().split('\n'):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip() if len(line) >= 78 else line[12:14].strip()
            element = element.capitalize()
            # Deduplicate by appending index
            if element in element_counts:
                element_counts[element] += 1
                key = f"{element}{element_counts[element]}"
            else:
                element_counts[element] = 1
                key = element
            coords[key] = (x, y, z)
        except (ValueError, IndexError):
            continue

    return coords


def generate_3d_coords(smiles: str) -> Tuple[Optional[Chem.Mol], List[Tuple[float, float, float]]]:
    """Generate 3D coordinates for a molecule using RDKit.

    Returns (mol_with_conformers, list_of_(x,y,z) tuples).
    Returns (None, []) if SMILES is invalid.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (None, [])
        mol = Chem.AddHs(mol)
        # Try embedding with random coordinates, fallback to distance geometry
        res = AllChem.EmbedMolecule(mol, randomSeed=42, maxAttempts=100)
        if res == -1:
            # Fallback: use distance geometry explicitly
            AllChem.Compute2DCoords(mol)
            conf = mol.GetConformer()
        else:
            # Optimize with MMFF
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                pass
            conf = mol.GetConformer()

        coords = []
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords.append((pos.x, pos.y, pos.z))

        return (mol, coords)
    except Exception:
        return (None, [])


def mol_to_atom_coords(mol: Chem.Mol) -> Dict[int, Tuple[float, float, float]]:
    """Extract atom index -> (x, y, z) from an RDKit mol with a conformer.
    Returns empty dict if no conformer is present.
    """
    coords: Dict[int, Tuple[float, float, float]] = {}
    try:
        if mol.GetNumConformers() == 0:
            return coords
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            coords[i] = (pos.x, pos.y, pos.z)
    except Exception:
        pass
    return coords


def _coords_cache_key(smi: str) -> str:
    """Generate a canonical key for coords_cache lookup.

    Strips atom map numbers and canonicalizes SMILES to ensure consistent
    matching regardless of whether the SMILES comes from frontend input,
    RXN parsing (with AM numbers), or an AddHs'd mol.
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return smi
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol)


def _build_pdb_mol_with_conformer(
    smi: str,
    atom_coords: Dict[int, Tuple[float, float, float]],
) -> Optional[Chem.Mol]:
    """Build a Mol object from SMILES with PDB 3D coordinates as a conformer.

    Used to store PDB coordinates in coords_cache as a Mol object (instead of
    a bare dict), which enables proper MCS-based coordinate injection in
    _build_am_to_3d_coord (matching the source code pattern).
    """
    from rdkit.Geometry.rdGeometry import Point3D
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        conf = Chem.Conformer(mol.GetNumAtoms())
        n_set = 0
        for atom_idx, (x, y, z) in atom_coords.items():
            if atom_idx < mol.GetNumAtoms():
                conf.SetAtomPosition(atom_idx, Point3D(x, y, z))
                n_set += 1
        mol.AddConformer(conf, assignId=True)
        if n_set == 0:
            return None
        return mol
    except Exception:
        return None


def load_pdb_coords_for_smiles(
    pdb_text: str,
    smiles: str,
    res_name: Optional[str] = None,
    chain_id: Optional[str] = None,
) -> Optional[Chem.Mol]:
    """Load 3D coordinates from PDB text and return a Mol with PDB coords as conformer.

    Attempts to create an RDKit Mol from the PDB text, find the best match
    to the target SMILES via MCS, and return a target Mol with PDB coordinates
    injected via the MCS atom correspondence.

    Args:
        pdb_text: PDB format text
        smiles: Target SMILES string to match against
        res_name: Optional residue name filter
        chain_id: Optional chain ID filter

    Returns:
        RDKit Mol with 3D conformer from PDB coordinates, or None if mapping fails.
    """
    from rdkit.Geometry.rdGeometry import Point3D
    try:
        # Parse PDB text into RDKit mol
        pdb_mol = Chem.MolFromPDBBlock(pdb_text, removeHs=True)
        if pdb_mol is None:
            pdb_mol = Chem.MolFromPDBBlock(pdb_text, removeHs=False)
        if pdb_mol is None:
            return None

        if pdb_mol.GetNumConformers() == 0:
            return None

        # Filter by chain if specified
        if chain_id:
            filtered_atoms = []
            for atom in pdb_mol.GetAtoms():
                info = atom.GetPDBResidueInfo()
                if info and info.GetChainId().strip() == chain_id.strip():
                    filtered_atoms.append(atom.GetIdx())
            if not filtered_atoms:
                return None

        # Parse target SMILES
        target_mol = Chem.MolFromSmiles(smiles)
        if target_mol is None:
            return None

        # Find MCS between PDB mol and target SMILES to establish atom correspondence
        try:
            pdb_copy = Chem.RWMol(pdb_mol)
            target_copy = Chem.RWMol(target_mol)
            try:
                Chem.SanitizeMol(pdb_copy)
                Chem.SanitizeMol(target_copy)
            except Exception:
                pass

            # Normalize for comparison: set isotope=atomicNum, atomMapNum=0
            for mol in [pdb_copy, target_copy]:
                for atom in mol.GetAtoms():
                    atom.SetIsotope(atom.GetAtomicNum())
                    atom.SetAtomMapNum(0)

            mcs = Chem.rdFMCS.FindMCS(
                [pdb_copy, target_copy],
                atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes,
                timeout=3,
            )

            if mcs.queryMol is None or mcs.numAtoms < 3:
                return None

            mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
            if mcs_mol is None:
                return None

            pdb_match = pdb_copy.GetSubstructMatch(mcs_mol)
            target_match = target_copy.GetSubstructMatch(mcs_mol)

            if not pdb_match or not target_match:
                return None

            # Inject PDB coords into target_mol via MCS correspondence
            pdb_conf = pdb_mol.GetConformer(0)
            target_conf = Chem.Conformer(target_mol.GetNumAtoms())
            n_set = 0
            for pdb_idx, target_idx in zip(pdb_match, target_match):
                pos = pdb_conf.GetAtomPosition(pdb_idx)
                target_conf.SetAtomPosition(target_idx, pos)
                n_set += 1

            if n_set < 3:
                return None

            target_mol.AddConformer(target_conf, assignId=True)
            return target_mol
        except Exception:
            return None
    except Exception:
        return None


def build_coords_cache_from_pdb(
    pdb_text: str,
    molecules: List[str],
) -> Dict[str, Optional[Chem.Mol]]:
    """Build a coords_cache from PDB 3D coordinates for a list of molecules.

    Zero-tolerance policy: only uses real PDB coordinates. No fallback to
    generated coordinates. Molecules without PDB mapping get None, which
    causes _build_am_to_3d_coord to skip them and the search to abort
    if they are core ligands.

    Args:
        pdb_text: PDB format text with 3D coordinates
        molecules: List of SMILES strings

    Returns:
        Dict mapping canonical SMILES to RDKit Mol objects with 3D conformers,
        or None for molecules where PDB coordinate extraction failed.
    """
    coords_cache = {}

    for smi in molecules:
        key = _coords_cache_key(smi)
        if key in coords_cache:
            continue

        # Try PDB-based mapping (no fallback to generated coordinates)
        pdb_mol = load_pdb_coords_for_smiles(pdb_text, smi)
        if pdb_mol is not None and pdb_mol.GetNumConformers() > 0:
            coords_cache[key] = pdb_mol
        else:
            # Zero tolerance: no mapping → None
            logger.warning("No PDB coordinates for %s, set to None (zero tolerance)", key[:40])
            coords_cache[key] = None

    return coords_cache


# ---------------------------------------------------------------------------
# PDB residue fragment extraction
# ---------------------------------------------------------------------------
_MAIN_CHAIN_HEAVY_ATOMS = frozenset({'N', 'CA', 'C', 'O'})


def extract_pdb_residue_fragment(
    pdb_text: str,
    res_name: str,
    res_num: int,
    chain: str,
    part: str,
) -> Optional[Chem.Mol]:
    """Extract a 3D fragment of a specific residue from PDB text.

    Mirrors Upload's structure.get_mol_from_pdb_residue():
      - Filters ATOM records by res_name/res_num/chain
      - Filters atoms by part ("side_chain" or "main_chain")
      - Uses Chem.MolFromPDBBlock for automatic bond inference

    Args:
        pdb_text: Full PDB format text
        res_name: 3-letter residue name (e.g., "SER")
        res_num: Residue sequence number
        chain: Chain identifier
        part: "side_chain" (CA+CB+sidechain) or "main_chain" (N,CA,C,O)

    Returns:
        RDKit Mol with 3D conformer (PDB coordinates) and inferred bonds,
        or None if extraction fails.
    """
    res_name = res_name.upper()
    chain = chain.strip()
    part = part.strip().lower()

    # Step 1: Filter ATOM records for target residue
    residue_lines = []
    for line in pdb_text.strip().split('\n'):
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            continue
        try:
            rec_name = line[17:20].strip().upper()
            rec_chain = line[21].strip()
            rec_seq = int(line[22:26])
            atom_name = line[12:16].strip().upper()

            if rec_name != res_name or rec_chain != chain or rec_seq != res_num:
                continue

            # Exclude hydrogen atoms
            element = line[76:78].strip().upper() if len(line) >= 78 else ''
            if element == 'H':
                continue

            # Filter by part
            if part == 'main_chain' and atom_name not in _MAIN_CHAIN_HEAVY_ATOMS:
                continue
            # side_chain: keep all non-H atoms (CA, CB, sidechain)

            residue_lines.append(line)
        except (ValueError, IndexError):
            continue

    if len(residue_lines) < 1:
        return None

    # Step 2: Build 3D Mol from filtered PDB block (RDKit auto-infers bonds)
    pdb_block = '\n'.join(residue_lines) + '\nEND\n'
    frag_mol = Chem.MolFromPDBBlock(pdb_block, removeHs=True)
    if frag_mol is None or frag_mol.GetNumConformers() == 0:
        return None

    return frag_mol


def _align_fragment_to_mol(
    frag_mol: Chem.Mol,
    target_mol: Chem.Mol,
) -> Optional[Chem.Mol]:
    """Align a PDB fragment (extracted ligand) to a target SMILES molecule via MCS.

    Returns a copy of target_mol with PDB 3D coordinates injected as a conformer,
    or None if alignment fails.
    """
    from rdkit.Geometry.rdGeometry import Point3D
    try:
        # Build isotope-labeled versions for MCS CompareIsotopes
        frag_copy = Chem.RWMol(frag_mol)
        target_copy = Chem.RWMol(target_mol)
        for mol in [frag_copy, target_copy]:
            for atom in mol.GetAtoms():
                atom.SetIsotope(atom.GetAtomicNum())
                atom.SetAtomMapNum(0)
        try:
            Chem.SanitizeMol(frag_copy)
            Chem.SanitizeMol(target_copy)
        except Exception:
            pass

        mcs = Chem.rdFMCS.FindMCS(
            [frag_copy.GetMol(), target_copy.GetMol()],
            atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes,
            timeout=3,
        )
        if mcs.queryMol is None or mcs.numAtoms < 3:
            return None

        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        if mcs_mol is None:
            return None

        frag_match = frag_copy.GetMol().GetSubstructMatch(mcs_mol)
        target_match = target_copy.GetMol().GetSubstructMatch(mcs_mol)
        if not frag_match or not target_match:
            return None

        # Inject PDB fragment coords into target_mol via MCS correspondence
        frag_conf = frag_mol.GetConformer(0)
        result_mol = Chem.RWMol(target_mol)
        result_conf = Chem.Conformer(result_mol.GetNumAtoms())
        n_set = 0
        for frag_idx, target_idx in zip(frag_match, target_match):
            pos = frag_conf.GetAtomPosition(frag_idx)
            result_conf.SetAtomPosition(target_idx, pos)
            n_set += 1

        if n_set < 3:
            return None

        result_mol.AddConformer(result_conf, assignId=True)
        return result_mol.GetMol()
    except Exception:
        return None


def _euclidean_distance(
    p1: Tuple[float, float, float],
    p2: Tuple[float, float, float],
) -> float:
    """Compute Euclidean distance between two 3D points."""
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2)


def _build_am_to_3d_coord(
    run_state: 'PredictionRunState',
    coords_cache: Optional[Dict[str, Chem.Mol]],
    reactants_smiles: str,
    residues: Optional[List[Dict[str, Any]]] = None,
    pdb_text: Optional[str] = None,
    residue_smiles_list: Optional[List[str]] = None,
    ligand_mappings: Optional[List[Dict[str, Any]]] = None,
) -> None:
    """Build am_to_3d_coord from PDB coordinates on reactant molecules.

    Processing order (precise → general, earlier steps take priority):
      ① Ligand mapping path: per-mapping PDB fragment extraction → MCS → MMFF
      ② Residue path: per-residue PDB fragment extraction → MCS → MMFF
      ③ General coords_cache path: whole-molecule MCS → MMFF
      ④ Build am_to_3d_coord from last conformers

    A shared _assigned_mol_ids set tracks which mols have been assigned
    coordinates. Once a mol is assigned by an earlier step, later steps
    MUST skip it to prevent overwriting precise coordinates with fuzzy ones.

    Aligns with Upload upload_models.py add_3d_conformer_to_inital_mols_from_pdb()
    (L1392-1437).
    """
    from rdkit.Geometry.rdGeometry import Point3D

    # Shared set: tracks mol indices that have already been assigned PDB coordinates.
    # All three processing steps (ligand → residue → general) share this set
    # to ensure 1:1 mapping and prevent coordinate overwrite.
    _assigned_mol_ids: Set[int] = set()

    def _get_mol_canonical_key(mol: Chem.Mol) -> Optional[str]:
        """Get canonical SMILES key for a mol (strip AM + isotope + H)."""
        try:
            mol_copy = Chem.RWMol(mol)
            for atom in mol_copy.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetIsotope(0)
            mol_no_h = Chem.RemoveHs(mol_copy.GetMol())
            return Chem.MolToSmiles(mol_no_h)
        except Exception:
            return None

    def _find_next_unassigned_mol(target_key: str) -> Optional[Tuple[int, Chem.Mol]]:
        """Find the next unassigned mol in initial_mols matching the canonical key.

        Returns (index, mol) tuple, or None if not found.
        Marks the found index as assigned in _assigned_mol_ids.
        """
        for idx, mol in enumerate(run_state.initial_mols):
            if idx in _assigned_mol_ids:
                continue
            mol_key = _get_mol_canonical_key(mol)
            if mol_key is not None and mol_key == target_key:
                _assigned_mol_ids.add(idx)
                return (idx, mol)
        return None

    def _mcs_align_and_mmff(target_mol: Chem.Mol, frag_mol: Chem.Mol,
                             min_mcs_atoms: int = 1) -> bool:
        """MCS-align a PDB fragment to a target mol, then MMFF-optimize.

        Adds the optimized conformer to target_mol if successful.
        Returns True if successful, False otherwise.

        Args:
            target_mol: The mol in initial_mols to receive PDB coordinates
            frag_mol: PDB fragment mol with 3D conformer
            min_mcs_atoms: Minimum number of MCS-matched atoms required
        """
        try:
            # Build tmp_mol with isotope=atomicNum for MCS CompareIsotopes
            tmp_mol = Chem.Mol(target_mol)
            for atom in tmp_mol.GetAtoms():
                if atom.GetSymbol() in {"*", }:
                    atom.SetAtomicNum(6)
                atom.SetIsotope(atom.GetAtomicNum())
            Chem.rdmolops.SanitizeMol(tmp_mol)
            AllChem.EmbedMolecule(tmp_mol)

            # Build frag_mol_for_mcs with isotope=atomicNum
            frag_mol_for_mcs = Chem.RWMol(frag_mol)
            for atom in frag_mol_for_mcs.GetAtoms():
                atom.SetIsotope(atom.GetAtomicNum())
            frag_mol_for_mcs = frag_mol_for_mcs.GetMol()

            if tmp_mol.GetNumAtoms() == 1:
                # Single atom: directly copy PDB coordinate
                conf0 = tmp_mol.GetConformer()
                pdb_conf = frag_mol.GetConformer()
                conf0.SetAtomPosition(0, pdb_conf.GetAtomPosition(0))
            else:
                # MCS matching
                mcs = Chem.rdFMCS.FindMCS(
                    (tmp_mol, frag_mol_for_mcs), timeout=2,
                    atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes,
                )
                if mcs.queryMol is None or mcs.numAtoms < min_mcs_atoms:
                    return False

                mcs2frag = frag_mol_for_mcs.GetSubstructMatch(mcs.queryMol)
                mcs2tmp = tmp_mol.GetSubstructMatch(mcs.queryMol)
                mol_to_frag_idx = dict(zip(mcs2tmp, mcs2frag))

                if not mol_to_frag_idx:
                    return False

                # MMFF force field
                mp = AllChem.MMFFGetMoleculeProperties(tmp_mol)
                if mp is None:
                    return False
                ff = AllChem.MMFFGetMoleculeForceField(tmp_mol, mp, confId=0)

                # Copy PDB fragment coords for MCS-matched atoms + constraints
                tmp_conformer = tmp_mol.GetConformer()
                frag_conformer = frag_mol.GetConformer()
                n_fixed = 0
                for atom in tmp_mol.GetAtoms():
                    atom_idx = atom.GetIdx()
                    if atom_idx in mol_to_frag_idx:
                        frag_coord = frag_conformer.GetAtomPosition(mol_to_frag_idx[atom_idx])
                        tmp_conformer.SetAtomPosition(atom_idx, frag_coord)
                        ff.MMFFAddPositionConstraint(atom_idx, .1, 100)
                        n_fixed += 1

                if n_fixed == 0:
                    return False

                # Constrained MMFF minimization
                ff.Minimize(maxIts=2000)

            # Add optimized conformer to original mol
            target_mol.AddConformer(tmp_mol.GetConformer(), assignId=True)
            return True

        except Exception:
            return False

    # ================================================================
    # ① Ligand mapping path: per-mapping PDB fragment extraction
    # Each mapping entry gets its own independent PDB coordinates.
    # This correctly handles multiple molecules with the same SMILES
    # (e.g., multiple water molecules at different PDB sites).
    # ================================================================
    if ligand_mappings and pdb_text:
        logger.info("Processing %d ligand mappings (consume-one-by-one mode)",
                     len(ligand_mappings))

        for mapping in ligand_mappings:
            smi = mapping.get('smiles', '').strip()
            res_name = (mapping.get('res_name', '') or '').upper()
            res_num = int(mapping.get('res_num', 0))
            chain = (mapping.get('chain', '') or '').strip()
            ligand_type = mapping.get('ligand_type', 'substrate')

            if not smi:
                continue

            target_key = _coords_cache_key(smi)
            result = _find_next_unassigned_mol(target_key)
            if result is None:
                logger.warning("Ligand mapping: no unassigned mol for %s (mapping: %s %s %d)",
                               target_key[:40], res_name, chain, res_num)
                continue

            idx, target_mol = result

            # Extract PDB residue fragment
            frag_mol = extract_pdb_residue_fragment(
                pdb_text, res_name, res_num, chain, ''
            )
            if frag_mol is None or frag_mol.GetNumConformers() == 0:
                logger.warning("Ligand mapping: PDB fragment extraction failed for %s %s %d",
                               res_name, chain, res_num)
                # Un-assign this mol so other steps can try
                _assigned_mol_ids.discard(idx)
                continue

            # MCS align + MMFF optimize
            success = _mcs_align_and_mmff(target_mol, frag_mol, min_mcs_atoms=1)
            if success:
                logger.info("  Ligand mapping: %s -> PDB %s %s %d (type=%s) ✓",
                            target_key[:40], res_name, chain, res_num, ligand_type)
            else:
                logger.warning("Ligand mapping: MCS/MMFF failed for %s (mapping: %s %s %d)",
                               target_key[:40], res_name, chain, res_num)
                # Un-assign this mol so other steps can try
                _assigned_mol_ids.discard(idx)

    # ================================================================
    # ② Residue-specific path: per-residue PDB fragment extraction
    # Mirrors Upload but uses per-residue PDB fragments.
    # ================================================================
    if residues and pdb_text and residue_smiles_list:
        _res_smi_keys = {smi: _coords_cache_key(smi) for smi in residue_smiles_list}

        for res_info in residues:
            res_name = (res_info.get('res_name', '') or '').upper()
            res_num = res_info.get('res_num', 0)
            chain = (res_info.get('chain', '') or '').strip()
            part = (res_info.get('part', 'side_chain') or 'side_chain').strip().lower()

            if res_name not in RESIDUE_SIDECHAIN_SMILES:
                continue

            res_smi = RESIDUE_SIDECHAIN_SMILES[res_name]
            target_key = _res_smi_keys.get(res_smi)
            if target_key is None:
                continue

            result = _find_next_unassigned_mol(target_key)
            if result is None:
                logger.debug("No unassigned mol for residue %s %d %s", res_name, res_num, chain)
                continue

            idx, target_mol = result

            # Extract PDB residue fragment
            frag_mol = extract_pdb_residue_fragment(pdb_text, res_name, res_num, chain, part)
            if frag_mol is None or frag_mol.GetNumConformers() == 0:
                logger.debug("PDB fragment extraction failed for %s %d %s (%s)",
                             res_name, res_num, chain, part)
                _assigned_mol_ids.discard(idx)
                continue

            # MCS alignment + MMFF constrained optimization
            success = _mcs_align_and_mmff(target_mol, frag_mol, min_mcs_atoms=2)
            if success:
                logger.debug("Residue 3D aligned: %s %d %s (%s) ✓",
                             res_name, res_num, chain, part)
            else:
                logger.debug("Residue MCS/MMFF failed for %s %d %s",
                             res_name, res_num, chain)
                _assigned_mol_ids.discard(idx)

    # ================================================================
    # ③ General coords_cache path: whole-molecule MCS
    # Handles mols without ligand/residue mappings that still have
    # coordinates in coords_cache. MUST skip already-assigned mols.
    # ================================================================
    if coords_cache:
        for idx, mol in enumerate(run_state.initial_mols):
            # CRITICAL: Skip mols already assigned by ligand/residue processing.
            # Without this, the general path would overwrite precise PDB coordinates
            # with fuzzy MCS results, and downstream am_to_3d_coord would pick
            # the wrong (last) conformer.
            if idx in _assigned_mol_ids:
                continue

            lookup_key = _get_mol_canonical_key(mol)
            if lookup_key is None:
                continue

            pdb_mol = coords_cache.get(lookup_key)
            if pdb_mol is None or pdb_mol.GetNumConformers() == 0:
                continue

            # Upload-style: MCS align + MMFF
            try:
                tmp_mol = Chem.Mol(mol)
                for atom in tmp_mol.GetAtoms():
                    if atom.GetSymbol() in {"*", }:
                        atom.SetAtomicNum(6)
                    atom.SetIsotope(atom.GetAtomicNum())
                Chem.rdmolops.SanitizeMol(tmp_mol)
                AllChem.EmbedMolecule(tmp_mol)

                pdb_mol_for_mcs = Chem.RWMol(pdb_mol)
                for atom in pdb_mol_for_mcs.GetAtoms():
                    atom.SetIsotope(atom.GetAtomicNum())
                pdb_mol_for_mcs = pdb_mol_for_mcs.GetMol()

                if tmp_mol.GetNumAtoms() == 1:
                    conf0 = tmp_mol.GetConformer()
                    pdb_conf = pdb_mol.GetConformer()
                    conf0.SetAtomPosition(0, pdb_conf.GetAtomPosition(0))
                else:
                    mcs = Chem.rdFMCS.FindMCS(
                        (tmp_mol, pdb_mol_for_mcs), timeout=2,
                        atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes,
                    )
                    if mcs.queryMol is None or mcs.numAtoms < 1:
                        continue

                    mcs2pdb = pdb_mol_for_mcs.GetSubstructMatch(mcs.queryMol)
                    mcs2tmp = tmp_mol.GetSubstructMatch(mcs.queryMol)
                    mol_to_pdb_idx = dict(zip(mcs2tmp, mcs2pdb))

                    if not mol_to_pdb_idx:
                        continue

                    mp = AllChem.MMFFGetMoleculeProperties(tmp_mol)
                    if mp is None:
                        continue

                    ff = AllChem.MMFFGetMoleculeForceField(tmp_mol, mp, confId=0)

                    pdb_conformer = pdb_mol.GetConformer()
                    tmp_conformer = tmp_mol.GetConformer()
                    n_fixed = 0
                    for atom in tmp_mol.GetAtoms():
                        atom_idx = atom.GetIdx()
                        if atom_idx in mol_to_pdb_idx:
                            pdb_coord = pdb_conformer.GetAtomPosition(mol_to_pdb_idx[atom_idx])
                            tmp_conformer.SetAtomPosition(atom_idx, pdb_coord)
                            ff.MMFFAddPositionConstraint(atom_idx, .1, 100)
                            n_fixed += 1

                    if n_fixed == 0:
                        continue

                    ff.Minimize(maxIts=2000)

                mol.AddConformer(tmp_mol.GetConformer(), assignId=True)

            except Exception as e:
                logger.debug("General coords_cache: MMFF failed for mol %s: %s", lookup_key[:30], e)
                continue

    # ================================================================
    # ④ Build am_to_3d_coord from mols with PDB 3D conformer (Upload L816-825)
    # Uses the last conformer (most accurate PDB-aligned 3D).
    # ================================================================
    am_to_3d_coord: Dict[int, Any] = {}
    for mol in run_state.initial_mols:
        if mol.GetNumConformers() < 1:
            continue  # No conformer at all
        try:
            # Use the last conformer (PDB-aligned 3D), which is the most accurate
            conf_id = mol.GetNumConformers() - 1
            conformer_3d = mol.GetConformer(conf_id)
            for atom in mol.GetAtoms():
                if atom.GetSymbol() != "H":
                    am = atom.GetAtomMapNum()
                    if am != 0:
                        am_to_3d_coord[am] = conformer_3d.GetAtomPosition(atom.GetIdx())
        except Exception:
            continue

    run_state.am_to_3d_coord = am_to_3d_coord
# ---------------------------------------------------------------------------
# Canonical SMILES cache
# ---------------------------------------------------------------------------
_canonical_cache: Dict[str, Optional[str]] = {}
_equivalent_cache: Dict[str, Optional[str]] = {}  # Separate cache for equivalent_smiles


def canonical_smiles(smiles: str) -> Optional[str]:
    """Return canonical SMILES for a SMILES string, with caching."""
    if smiles in _canonical_cache:
        return _canonical_cache[smiles]
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            _canonical_cache[smiles] = None
            return None
        result = Chem.MolToSmiles(mol)
        _canonical_cache[smiles] = result
        return result
    except Exception:
        _canonical_cache[smiles] = None
        return None


def make_mol_equivalent(mol: Chem.Mol) -> Chem.Mol:
    """Normalise a molecule so that connectivity-equivalent copies yield the same SMILES.
    
    Removes explicit hydrogens, strips atom-mapping numbers, sanitises.
    This is the paper's dedup technique - molecules with different H representations
    but same connectivity become identical.
    """
    try:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
            atom.SetIsotope(0)  # FIX: 重置 isotope，避免同位素标注导致等效分子产生不同 SMILES
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        Chem.SanitizeMol(mol_no_h)
        return mol_no_h
    except Exception:
        return Chem.Mol()


def equivalent_smiles(smiles: str) -> Optional[str]:
    """Return SMILES after make_mol_equivalent normalization (with caching).
    
    Uses a SEPARATE cache from canonical_smiles to avoid cache collision.
    """
    if smiles in _equivalent_cache:
        return _equivalent_cache[smiles]
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            _equivalent_cache[smiles] = None
            return None
        equiv = make_mol_equivalent(mol)
        result = Chem.MolToSmiles(equiv)
        _equivalent_cache[smiles] = result
        return result
    except Exception:
        _equivalent_cache[smiles] = None
        return None


# ===========================================================================
# Atom Mapping Verification & Correction
# ===========================================================================
def _count_bond_changes(
    r_mols: List[Chem.Mol],
    p_mols: List[Chem.Mol],
) -> Tuple[int, Set[tuple], Set[tuple]]:
    """Count bond changes between reactants and products using current AM numbers.

    Returns (num_changed_bonds, bonds_only_in_reactants, bonds_only_in_products).
    """
    r_bond_set: Set[tuple] = set()
    for mol in r_mols:
        for bond in mol.GetBonds():
            am1 = bond.GetBeginAtom().GetAtomMapNum()
            am2 = bond.GetEndAtom().GetAtomMapNum()
            if am1 != 0 and am2 != 0:
                r_bond_set.add((frozenset((am1, am2)), bond.GetBondTypeAsDouble()))

    p_bond_set: Set[tuple] = set()
    for mol in p_mols:
        for bond in mol.GetBonds():
            am1 = bond.GetBeginAtom().GetAtomMapNum()
            am2 = bond.GetEndAtom().GetAtomMapNum()
            if am1 != 0 and am2 != 0:
                p_bond_set.add((frozenset((am1, am2)), bond.GetBondTypeAsDouble()))

    only_r = r_bond_set - p_bond_set
    only_p = p_bond_set - r_bond_set
    return len(only_r) + len(only_p), only_r, only_p


def remap_reaction_by_mcs(
    reactants_smiles: str,
    products_smiles: str,
) -> Tuple[str, str]:
    """Re-assign atom map numbers using MCS-based optimal mapping.

    Indigo's automap can produce suboptimal mappings that inflate the reaction
    center (e.g., swapping a carboxyl =O with a water molecule). This function
    uses Maximum Common Substructure (MCS) to find the mapping that minimises
    the number of bond changes, which is chemically more correct.

    Algorithm:
      1. Parse reactant/product SMILES into individual molecules
      2. For each (reactant, product) pair, compute MCS
      3. Assign AM numbers based on MCS atom correspondence
      4. Handle unmatched atoms (consumed/produced) separately
      5. Verify the new mapping has fewer bond changes

    Returns (corrected_reactants_smiles, corrected_products_smiles).
    If correction fails or doesn't improve, returns the original SMILES.
    """
    try:
        from rdkit.Chem import rdFMCS

        # Parse molecules
        r_smi_list = [s.strip() for s in reactants_smiles.split('.') if s.strip()]
        p_smi_list = [s.strip() for s in products_smiles.split('.') if s.strip()]

        r_mols_orig = [Chem.MolFromSmiles(s) for s in r_smi_list]
        p_mols_orig = [Chem.MolFromSmiles(s) for s in p_smi_list]

        r_mols_orig = [m for m in r_mols_orig if m is not None]
        p_mols_orig = [m for m in p_mols_orig if m is not None]

        if not r_mols_orig or not p_mols_orig:
            return reactants_smiles, products_smiles

        # Count bond changes with current mapping
        current_changes, _, _ = _count_bond_changes(r_mols_orig, p_mols_orig)

        # Only re-map if there are bond changes AND some atoms are mapped
        has_am = any(atom.GetAtomMapNum() != 0 for mol in r_mols_orig for atom in mol.GetAtoms())
        if not has_am or current_changes == 0:
            return reactants_smiles, products_smiles

        # ---- Step 1: Build a bipartite matching using MCS ----
        # For each (reactant_mol, product_mol) pair, compute MCS and count
        # matched atoms. Use this to find the best assignment.

        # Strip AM for MCS computation
        def strip_am(mol: Chem.Mol) -> Chem.Mol:
            """Create a copy with AM and isotope stripped."""
            copy = Chem.RWMol(mol)
            for atom in copy.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetIsotope(0)
            return Chem.RemoveHs(copy.GetMol())

        # Compute MCS scores for all (reactant, product) pairs
        mcs_scores: Dict[Tuple[int, int], int] = {}  # (r_idx, p_idx) -> num_mcs_atoms
        mcs_matches: Dict[Tuple[int, int], Tuple[Tuple, Tuple]] = {}  # (r_idx, p_idx) -> (r_match, p_match)

        for ri, r_mol in enumerate(r_mols_orig):
            r_stripped = strip_am(r_mol)
            for pi, p_mol in enumerate(p_mols_orig):
                p_stripped = strip_am(p_mol)
                try:
                    mcs = rdFMCS.FindMCS(
                        [r_stripped, p_stripped],
                        atomCompare=rdFMCS.AtomCompare.CompareElements,
                        bondCompare=rdFMCS.BondCompare.CompareOrder,
                        timeout=5,
                    )
                    if mcs.numAtoms > 0:
                        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                        if mcs_mol:
                            r_match = r_stripped.GetSubstructMatch(mcs_mol)
                            p_match = p_stripped.GetSubstructMatch(mcs_mol)
                            mcs_scores[(ri, pi)] = mcs.numAtoms
                            mcs_matches[(ri, pi)] = (r_match, p_match)
                except Exception:
                    pass

        # ---- Step 2: Greedy assignment: assign each reactant to best product ----
        # Use a greedy approach: sort all (r, p) pairs by MCS score descending,
        # then assign each reactant to its best unassigned product.
        assigned_r: Set[int] = set()
        assigned_p: Set[int] = set()
        assignments: List[Tuple[int, int]] = []  # (r_idx, p_idx)

        for (ri, pi), score in sorted(mcs_scores.items(), key=lambda x: -x[1]):
            if ri in assigned_r or pi in assigned_p:
                continue
            # Only assign if MCS covers a significant portion of the smaller mol
            r_stripped = strip_am(r_mols_orig[ri])
            p_stripped = strip_am(p_mols_orig[pi])
            min_atoms = min(r_stripped.GetNumAtoms(), p_stripped.GetNumAtoms())
            if min_atoms > 0 and score >= min_atoms * 0.5:
                assigned_r.add(ri)
                assigned_p.add(pi)
                assignments.append((ri, pi))

        # ---- Step 3: Assign atom map numbers based on MCS matching ----
        # Create working copies
        r_mols_work = [Chem.RWMol(m) for m in r_mols_orig]
        p_mols_work = [Chem.RWMol(m) for m in p_mols_orig]

        # Clear all existing AM numbers
        for mol in r_mols_work + p_mols_work:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
                atom.SetIsotope(0)

        am_counter = 1

        # For each assigned (r, p) pair, map atoms via MCS
        for ri, pi in assignments:
            r_match, p_match = mcs_matches[(ri, pi)]
            r_stripped = strip_am(r_mols_orig[ri])
            p_stripped = strip_am(p_mols_orig[pi])

            # Get matches on working copies
            r_work_stripped = strip_am(r_mols_work[ri].GetMol())
            p_work_stripped = strip_am(p_mols_work[pi].GetMol())

            try:
                mcs_result = rdFMCS.FindMCS(
                    [r_work_stripped, p_work_stripped],
                    atomCompare=rdFMCS.AtomCompare.CompareElements,
                    bondCompare=rdFMCS.BondCompare.CompareOrder,
                    timeout=5,
                )
                mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString) if mcs_result.numAtoms > 0 else None
            except Exception:
                mcs_mol = None

            if mcs_mol is None:
                continue

            r_work_match = r_mols_work[ri].GetMol().GetSubstructMatch(mcs_mol)
            p_work_match = p_mols_work[pi].GetMol().GetSubstructMatch(mcs_mol)

            if not r_work_match or not p_work_match:
                # Fallback: try with RemoveHs
                r_no_h = Chem.RemoveHs(r_mols_work[ri].GetMol())
                p_no_h = Chem.RemoveHs(p_mols_work[pi].GetMol())
                r_work_match = r_no_h.GetSubstructMatch(mcs_mol)
                p_work_match = p_no_h.GetSubstructMatch(mcs_mol)

            if not r_work_match or not p_work_match:
                continue

            # Assign AM numbers for matched atoms
            r_matched_set = set(r_work_match)
            p_matched_set = set(p_work_match)
            for r_idx, p_idx in zip(r_work_match, p_work_match):
                if r_idx < r_mols_work[ri].GetNumAtoms():
                    r_mols_work[ri].GetAtomWithIdx(r_idx).SetAtomMapNum(am_counter)
                if p_idx < p_mols_work[pi].GetNumAtoms():
                    p_mols_work[pi].GetAtomWithIdx(p_idx).SetAtomMapNum(am_counter)
                am_counter += 1

            # Match unmatched atoms by element between paired reactant/product
            # This handles cases like: water O consumed → product OH produced
            # Both should share the same AM number for 3D coordinate mapping
            r_unmatched = [(atom.GetIdx(), atom.GetSymbol())
                           for atom in r_mols_work[ri].GetAtoms()
                           if atom.GetIdx() not in r_matched_set and atom.GetSymbol() != "H"
                           and atom.GetAtomMapNum() == 0]
            p_unmatched = [(atom.GetIdx(), atom.GetSymbol())
                           for atom in p_mols_work[pi].GetAtoms()
                           if atom.GetIdx() not in p_matched_set and atom.GetSymbol() != "H"
                           and atom.GetAtomMapNum() == 0]

            # Greedy element matching: pair up same-element atoms
            p_used: Set[int] = set()
            for r_idx, r_elem in r_unmatched:
                for p_local_idx, (p_idx, p_elem) in enumerate(p_unmatched):
                    if p_local_idx in p_used:
                        continue
                    if r_elem == p_elem:
                        r_mols_work[ri].GetAtomWithIdx(r_idx).SetAtomMapNum(am_counter)
                        p_mols_work[pi].GetAtomWithIdx(p_idx).SetAtomMapNum(am_counter)
                        p_used.add(p_local_idx)
                        am_counter += 1
                        break

            # Assign unique AM numbers for remaining unmatched reactant atoms
            for r_idx, r_elem in r_unmatched:
                atom = r_mols_work[ri].GetAtomWithIdx(r_idx)
                if atom.GetAtomMapNum() == 0:
                    atom.SetAtomMapNum(am_counter)
                    am_counter += 1

            # Assign unique AM numbers for remaining unmatched product atoms
            for p_local_idx, (p_idx, p_elem) in enumerate(p_unmatched):
                if p_local_idx in p_used:
                    continue
                atom = p_mols_work[pi].GetAtomWithIdx(p_idx)
                if atom.GetAtomMapNum() == 0:
                    atom.SetAtomMapNum(am_counter)
                    am_counter += 1

        # Assign AM for unassigned molecules (cross-mapping by element)
        # Try to match unassigned reactant mols with unassigned product mols
        # by element composition. This handles cases like water → water.
        unassigned_r_indices = [ri for ri in range(len(r_mols_work)) if ri not in assigned_r]
        unassigned_p_indices = [pi for pi in range(len(p_mols_work)) if pi not in assigned_p]

        # Build element signatures for matching
        def _elem_signature(mol: Chem.Mol) -> str:
            """Element composition string for quick comparison."""
            from collections import Counter
            counts = Counter(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() != "H")
            return ''.join(f'{elem}{count}' for elem, count in sorted(counts.items()))

        # Greedy matching of unassigned mols by element signature
        p_assigned_now: Set[int] = set()
        for ri in unassigned_r_indices:
            r_sig = _elem_signature(r_mols_work[ri].GetMol())
            best_pi = None
            for pi in unassigned_p_indices:
                if pi in p_assigned_now:
                    continue
                p_sig = _elem_signature(p_mols_work[pi].GetMol())
                if r_sig == p_sig:
                    best_pi = pi
                    break
            if best_pi is not None:
                # Match atoms by element between these two mols
                p_assigned_now.add(best_pi)
                r_atoms = [(a.GetIdx(), a.GetSymbol()) for a in r_mols_work[ri].GetAtoms() if a.GetSymbol() != "H"]
                p_atoms = [(a.GetIdx(), a.GetSymbol()) for a in p_mols_work[best_pi].GetAtoms() if a.GetSymbol() != "H"]

                p_atom_used: Set[int] = set()
                for r_idx, r_elem in r_atoms:
                    for p_local_idx, (p_idx, p_elem) in enumerate(p_atoms):
                        if p_local_idx in p_atom_used:
                            continue
                        if r_elem == p_elem:
                            r_mols_work[ri].GetAtomWithIdx(r_idx).SetAtomMapNum(am_counter)
                            p_mols_work[best_pi].GetAtomWithIdx(p_idx).SetAtomMapNum(am_counter)
                            p_atom_used.add(p_local_idx)
                            am_counter += 1
                            break
                    else:
                        # No matching element found — assign unique AM
                        r_mols_work[ri].GetAtomWithIdx(r_idx).SetAtomMapNum(am_counter)
                        am_counter += 1
                # Assign remaining unmatched product atoms
                for p_local_idx, (p_idx, p_elem) in enumerate(p_atoms):
                    if p_local_idx not in p_atom_used:
                        p_mols_work[best_pi].GetAtomWithIdx(p_idx).SetAtomMapNum(am_counter)
                        am_counter += 1
            else:
                # No matching product mol — assign unique AMs
                for atom in r_mols_work[ri].GetAtoms():
                    if atom.GetSymbol() != "H" and atom.GetAtomMapNum() == 0:
                        atom.SetAtomMapNum(am_counter)
                        am_counter += 1

        # Assign remaining unassigned product mols
        for pi in unassigned_p_indices:
            if pi in p_assigned_now:
                continue
            for atom in p_mols_work[pi].GetAtoms():
                if atom.GetSymbol() != "H" and atom.GetAtomMapNum() == 0:
                    atom.SetAtomMapNum(am_counter)
                    am_counter += 1

        # ---- Step 3.5: Fix orphan atoms (atom conservation) ----
        # After all assignments, atoms that only appear in reactants ("lost") or
        # only in products ("new") violate atom conservation. This happens when:
        #   - A reactant molecule is consumed (e.g., water providing OH to product)
        #     but its atoms appear in a different product molecule than expected
        #   - The greedy MCS pairing assigns the reactant mol to a different product,
        #     leaving the consumed atoms orphaned
        #
        # Fix: match lost reactant atoms with new product atoms by element,
        # replacing the "new" AM number with the "lost" AM number.
        # Example: water O:17 consumed → product OH gets new AM:15
        #          → Replace AM:15 with AM:17 (the O came from water 17)
        r_am_info: Dict[int, Tuple[str, int, int]] = {}  # am -> (element, mol_idx, atom_idx)
        p_am_info: Dict[int, Tuple[str, int, int]] = {}  # am -> (element, mol_idx, atom_idx)

        for ri, mol in enumerate(r_mols_work):
            for atom in mol.GetAtoms():
                am = atom.GetAtomMapNum()
                if am > 0:
                    r_am_info[am] = (atom.GetSymbol(), ri, atom.GetIdx())

        for pi, mol in enumerate(p_mols_work):
            for atom in mol.GetAtoms():
                am = atom.GetAtomMapNum()
                if am > 0:
                    p_am_info[am] = (atom.GetSymbol(), pi, atom.GetIdx())

        lost_ams = set(r_am_info.keys()) - set(p_am_info.keys())
        new_ams = set(p_am_info.keys()) - set(r_am_info.keys())

        if lost_ams and new_ams:
            new_used: Set[int] = set()
            merge_count = 0
            for lost_am in sorted(lost_ams):
                lost_elem = r_am_info[lost_am][0]
                for new_am in sorted(new_ams):
                    if new_am in new_used:
                        continue
                    if p_am_info[new_am][0] == lost_elem:
                        # Replace the "new" AM with the "lost" AM in the product
                        pi = p_am_info[new_am][1]
                        atom_idx = p_am_info[new_am][2]
                        p_mols_work[pi].GetAtomWithIdx(atom_idx).SetAtomMapNum(lost_am)
                        new_used.add(new_am)
                        merge_count += 1
                        break
            if merge_count > 0:
                logger.info(
                    "remap_reaction_by_mcs Step 3.5: merged %d orphan atom(s) "
                    "(lost=%s, new=%s)",
                    merge_count, sorted(lost_ams), sorted(new_ams),
                )

        # ---- Step 4: Convert to SMILES and verify ----
        r_new_smiles_list = []
        for mol in r_mols_work:
            try:
                smi = Chem.MolToSmiles(mol.GetMol(), isomericSmiles=True)
                if smi:
                    r_new_smiles_list.append(smi)
            except Exception:
                pass

        p_new_smiles_list = []
        for mol in p_mols_work:
            try:
                smi = Chem.MolToSmiles(mol.GetMol(), isomericSmiles=True)
                if smi:
                    p_new_smiles_list.append(smi)
            except Exception:
                pass

        if not r_new_smiles_list or not p_new_smiles_list:
            return reactants_smiles, products_smiles

        # Count bond changes with new mapping
        r_new_mols = [Chem.MolFromSmiles(s) for s in r_new_smiles_list]
        p_new_mols = [Chem.MolFromSmiles(s) for s in p_new_smiles_list]
        r_new_mols = [m for m in r_new_mols if m is not None]
        p_new_mols = [m for m in p_new_mols if m is not None]

        new_changes, _, _ = _count_bond_changes(r_new_mols, p_new_mols)

        # Check atom conservation for both mappings
        def _count_orphan_atoms(r_mols_list, p_mols_list) -> int:
            """Count atoms that violate conservation (appear in only one side)."""
            r_ams = set(a.GetAtomMapNum() for m in r_mols_list for a in m.GetAtoms() if a.GetAtomMapNum() > 0)
            p_ams = set(a.GetAtomMapNum() for m in p_mols_list for a in m.GetAtoms() if a.GetAtomMapNum() > 0)
            return len(r_ams - p_ams) + len(p_ams - r_ams)

        original_orphans = _count_orphan_atoms(r_mols_orig, p_mols_orig)
        new_orphans = _count_orphan_atoms(r_new_mols, p_new_mols)

        # Use new mapping if:
        # 1. It reduces bond changes, OR
        # 2. Bond changes are equal AND it improves atom conservation
        use_new = False
        reason = ""
        if new_changes < current_changes:
            use_new = True
            reason = f"fewer bond changes ({current_changes} -> {new_changes})"
        elif new_changes == current_changes and new_orphans < original_orphans:
            use_new = True
            reason = f"better atom conservation (orphans: {original_orphans} -> {new_orphans})"
        elif new_changes == current_changes and new_orphans == 0 and original_orphans == 0:
            # Both are equivalent — keep original for stability
            pass

        if use_new:
            logger.info("remap_reaction_by_mcs: using MCS mapping — %s", reason)
            return '.'.join(r_new_smiles_list), '.'.join(p_new_smiles_list)
        else:
            logger.info(
                "remap_reaction_by_mcs: keeping original mapping "
                "(bond changes: %d->%d, orphans: %d->%d)",
                current_changes, new_changes, original_orphans, new_orphans,
            )
            return reactants_smiles, products_smiles

    except Exception as e:
        logger.warning("remap_reaction_by_mcs failed: %s — keeping original mapping", e)
        return reactants_smiles, products_smiles


# ===========================================================================
# Reaction Center Identification
# ===========================================================================
def identify_reaction_center_atom_maps(
    reactants_smiles: str,
    products_smiles: str,
) -> Set[int]:
    """Identify reaction center atom map numbers from overall reaction.

    Equivalent to Upload's Prediction.step_rc_ams (upload_models.py L700-706):
    Returns a set of atom map numbers that participate in bond formation/cleavage,
    using the bond-level symmetric difference method.

    In the Upload paper, step_rc_ams comes from:
      1. reaction.get_reaction_centres_atoms() → rxn_reaction_centres_am (Django DB)
      2. rxn_am_to_scheme_am maps them to scheme-level atom map numbers
      3. Filtered to exclude water/proton molecules

    In our system (no Django DB), we compute the RC atom map numbers directly
    from the atom-mapped SMILES (from Ketcher RXN or user input):
      1. Parse reactants/products SMILES → Mol objects
      2. Bond-level symmetric difference → changed bonds
      3. Return the atom map numbers of changed bond atoms
      4. Filter out water/proton molecules (same as paper)

    Returns:
        Set of atom map numbers (NOT atomic numbers) involved in the reaction center.
    """
    try:
        from predict_mechanism_util import mol_is_water_or_proton

        r_mols = [Chem.MolFromSmiles(s.strip()) for s in reactants_smiles.split('.') if s.strip()]
        p_mols = [Chem.MolFromSmiles(s.strip()) for s in products_smiles.split('.') if s.strip()]
        if not r_mols or not p_mols:
            return set()
        r_mols = [m for m in r_mols if m is not None]
        p_mols = [m for m in p_mols if m is not None]
        if not r_mols or not p_mols:
            return set()

        # ---- Bond-level symmetric difference (Upload's method) ----
        # Collect bond sets: (frozenset(atom_map1, atom_map2), bond_type)
        r_has_am = any(atom.GetAtomMapNum() != 0 for mol in r_mols for atom in mol.GetAtoms())
        p_has_am = any(atom.GetAtomMapNum() != 0 for mol in p_mols for atom in mol.GetAtoms())

        if r_has_am and p_has_am:
            # --- Atom-mapped path: return atom map numbers (like paper's step_rc_ams) ---
            r_bond_set: Set[tuple] = set()
            for mol in r_mols:
                for bond in mol.GetBonds():
                    am1 = bond.GetBeginAtom().GetAtomMapNum()
                    am2 = bond.GetEndAtom().GetAtomMapNum()
                    if am1 != 0 and am2 != 0:
                        r_bond_set.add((frozenset((am1, am2)), bond.GetBondTypeAsDouble()))

            p_bond_set: Set[tuple] = set()
            for mol in p_mols:
                for bond in mol.GetBonds():
                    am1 = bond.GetBeginAtom().GetAtomMapNum()
                    am2 = bond.GetEndAtom().GetAtomMapNum()
                    if am1 != 0 and am2 != 0:
                        p_bond_set.add((frozenset((am1, am2)), bond.GetBondTypeAsDouble()))

            # Symmetric difference gives changed bonds
            changed_bonds = r_bond_set ^ p_bond_set
            rc_atom_maps: Set[int] = set()
            for bond_info in changed_bonds:
                rc_atom_maps.update(bond_info[0])  # frozenset of (am1, am2)

            # Filter out water/proton molecules (same as paper's step_rc_ams L706)
            filtered_rc: Set[int] = set()
            for am in rc_atom_maps:
                # Find the molecule containing this atom map
                for mol in r_mols + p_mols:
                    for atom in mol.GetAtoms():
                        if atom.GetAtomMapNum() == am:
                            if not mol_is_water_or_proton(mol):
                                filtered_rc.add(am)
                            break
                    else:
                        continue
                    break
            return filtered_rc

        # --- Fallback (no atom maps): return atomic numbers for reaction center elements ---
        r_counter: Dict[int, int] = {}
        for mol in r_mols:
            for atom in mol.GetAtoms():
                r_counter[atom.GetAtomicNum()] = r_counter.get(atom.GetAtomicNum(), 0) + 1

        p_counter: Dict[int, int] = {}
        for mol in p_mols:
            for atom in mol.GetAtoms():
                p_counter[atom.GetAtomicNum()] = p_counter.get(atom.GetAtomicNum(), 0) + 1

        changed: Set[int] = set()
        all_elems = set(r_counter.keys()) | set(p_counter.keys())
        for elem in all_elems:
            if abs(r_counter.get(elem, 0) - p_counter.get(elem, 0)) > 0:
                changed.add(elem)
        return changed
    except Exception:
        return set()


# ===========================================================================
# Upload Source Code-Based Matching Logic
# ===========================================================================
# Implemented from Upload models.py PredictionRun.apply_rules():
# - PredictionRunState: State management class (equivalent to Upload's PredictionRun)
# - _lazy_parse_rule(): Lazy per-rule SMARTS parsing (avoids OOM with 51k rules)
# - _get_templates_for_duplicated_reactants(): Duplicate reactant handling (Upload line 1027-1082)
# - _get_reaction_data_for_products(): ReactionResult + scoring (Upload line 1116-1148)
# - _apply_rule_to_configuration(): Per-config rule application (Upload line 1151-1190)
# - _get_next_node_to_explore(): Dijkstra exploration (Upload line 1192-1219)
# - bidirectional_search(): Main entry -- SMILES->Mol parsing, graph setup, Dijkstra search
#
# Key Upload differences from our old code:
# 1. Uses GetSubstructMatches (PLURAL) -- returns ALL matches per mol x part
# 2. Stores matches as {rule_am: mol_atom_idx} dicts for atom mapping
# 3. Uses nx.MultiGraph with integer node IDs
# 4. Uses nx.single_source_dijkstra_path_length for exploration ordering
# 5. Handles duplicate reactants with merge_rule_parts + MolToSmarts dynamic rewriting
# 6. ReactionResult computes max_bond_length from 3D coords with rotation correction
# ===========================================================================


# ---------------------------------------------------------------------------
# Step 1: SMILES -> Mol Conversion & Equivalence Caching
# (Upload models.py: add_and_get_equivalent_mols, set_and_get_configuration_id)
# ---------------------------------------------------------------------------
class PredictionRunState:
    """Holds all state for a mechanism prediction run.

    Equivalent to Upload's PredictionRun class, but without Django dependencies.
    Manages:
    - smiles_to_mol: equivalence-normalized Mol cache
    - mols_2_conf_id: configuration dedup (frozenset of mols -> int ID)
    - G: NetworkX MultiGraph for the configuration landscape
    - mol_to_part_matches: precomputed substructure match cache
    - rule_to_bond_changes_am_heavy: bond change analysis per rule
    """

    def __init__(self, max_nodes: int = DEFAULT_MAX_CONFIGS):
        self.smiles_to_mol: Dict[str, Chem.Mol] = {}
        self.mols_2_conf_id: Dict[FrozenSet, int] = {}
        self.G: nx.MultiGraph = nx.MultiGraph()
        self.mol_to_part_matches: Dict[Chem.Mol, Dict[Chem.Mol, List[Dict[int, int]]]] = {}
        self.rule_reactants_to_products: Dict[str, Dict] = {}
        self.max_nodes = max_nodes
        self.initial_mols: List[Chem.Mol] = []
        self.initial_product_mols: List[Chem.Mol] = []
        self.to_explore: Set[int] = set()
        self.explored: Set[int] = set()
        self.explore_from_r: bool = True  # alternating flag for bidirectional
        self.total_expanded: int = 0
        # Part C: Overall reaction center atom map numbers (equivalent to Upload's step_rc_ams)
        # Used by get_bonds_away_from_rxn_rc for score × 0.50 discount
        self.overall_rc_atom_maps: Set[int] = set()

        # Part D: Atom map number → 3D coordinate mapping (Upload upload_models.py L816-825)
        # Built from PDB coordinates of reactant molecules. When available, bond length
        # is computed using reactant PDB coords (Upload's _quick_parse method).
        # When not available, falls back to product 3D coords (compute_max_bond_length_from_3d).
        self.am_to_3d_coord: Dict[int, Any] = {}  # {atom_map_num: Point3D}

        # NOTE: mol_to_h_mol cache has been REMOVED.
        # After Fix #1 (graph nodes store AddHs mol), all mols in the graph already
        # have explicit H. No separate AddHs cache is needed anymore.
        # Previously, the graph stored implicit-H mol but RunReactants needed explicit-H
        # mol (mol_h), causing set(combo) ≠ graph node mols → set difference failure.

    def add_and_get_equivalent_mols(self, *mols: Chem.Mol) -> List[Chem.Mol]:
        """Normalize molecules via make_mol_equivalent, dedup by SMILES.

        Upload models.py line 796-806.
        Uses canonical SMILES of the equivalent (de-H'd, de-mapped) form as cache key.
        """
        equivalent_mols = []
        for mol in mols:
            try:
                smiles = Chem.MolToSmiles(make_mol_equivalent(Chem.Mol(mol)))
            except Exception:
                smiles = Chem.MolToSmiles(mol)
            if smiles not in self.smiles_to_mol:
                self.smiles_to_mol[smiles] = mol
                for atom in mol.GetAtoms():
                    atom.SetProp("smiles", smiles)
            equivalent_mols.append(self.smiles_to_mol[smiles])
        return equivalent_mols

    def set_and_get_configuration_id(self, mols: Tuple[Chem.Mol, ...], **kwargs) -> int:
        """Check if configuration is unique; add to graph if new. Returns config ID.

        Upload models.py line 808-814.
        """
        fs = frozenset(mols)
        configuration_id = self.mols_2_conf_id.setdefault(fs, len(self.mols_2_conf_id) + 1)
        if configuration_id not in self.G.nodes:
            self.G.add_node(configuration_id, **kwargs)
            self.G.nodes[configuration_id]["mols"] = set(mols)
        return configuration_id


# ---------------------------------------------------------------------------
# Step 2: mol_to_part_matches Precomputation
# (Upload models.py line 1003-1013)
# ---------------------------------------------------------------------------
def _lazy_parse_rule(rule: Dict[str, Any]) -> None:
    """Lazily parse a rule's SMARTS parts on first encounter.

    Instead of preparsing all 51k rules upfront (which causes OOM),
    parse only the rules that the search loop actually processes.
    This mirrors the Upload source code approach of _preparse_rules
    but applied lazily per-rule.
    """
    if '_rule_parts' in rule:
        return  # Already parsed

    r_smarts = rule.get('reactant_smarts', '')
    p_smarts = rule.get('product_smarts', '')

    # Lazy ReactionFromSmarts — None until RunReactants needs it
    rule['_reactant_rxn'] = None

    # Parse reactant SMARTS parts (essential for substructure matching)
    rule['_rule_parts'] = []
    rule['_reactant_smarts_rule_parts'] = []
    rule['_rule_am_to_part_id'] = {}
    rule['_rule_am_to_atom'] = ({}, {})
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

    # Product rule parts — only parse if we have reactant parts (needed for
    # duplicated reactants handling and bond change computation)
    rule['_product_rule_parts'] = []
    rule['_product_smarts_rule_parts'] = []
    product_bonds_ams = set()
    if rule['_rule_parts']:
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

    # Bond changes and strict reaction centers
    rule['_bond_changes_am'] = (
        product_bonds_ams - reactant_bonds_ams,
        reactant_bonds_ams - product_bonds_ams,
    )
    rule['_strict_ams'] = {
        am for bonds in rule['_bond_changes_am'] for bond in bonds for am in bond
    }

    # Heavy-atom bond changes
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

    rule['_num_reactants'] = r_smarts.count('.') + 1


def _release_rule_rdkit_objects(rule: Dict[str, Any]) -> None:
    """Free RDKit mol objects from a rule to reclaim memory.

    Keeps the raw dict fields (reactant_smarts, product_smarts, etc.)
    but removes parsed RDKit objects (_rule_parts, _reactant_rxn, etc.).
    The rule can be re-parsed by _lazy_parse_rule() on next encounter.
    """
    for key in (
        '_rule_parts', '_reactant_smarts_rule_parts',
        '_product_rule_parts', '_product_smarts_rule_parts',
        '_reactant_rxn',
        '_rule_am_to_part_id', '_rule_am_to_atom',
        '_rule_am_to_part_size', '_bond_changes_am',
        '_strict_ams', '_bond_changes_am_heavy', '_num_reactants',
    ):
        rule.pop(key, None)


# ---------------------------------------------------------------------------
# Step 2.1: SMARTS string pre-filtering (the key OOM + speed optimization)
# ---------------------------------------------------------------------------
# Common organic elements — used for pre-filtering rules against molecules.
# We include '*' and '#' because SMARTS uses these for wildcard/any-atom matches.
# '[' and ']' are structural markers in SMARTS, not element symbols.
_ORGANIC_ELEMENTS = frozenset('CNOSPFIHBrcl')

# Map SMARTS element notations (inside brackets) to canonical element symbols.
# E.g., "[nH]" → n (aromatic nitrogen), "[se]" → Se (selenium)
_SMARTS_ELEMENT_MAP = {
    'c': 'C', 'n': 'N', 'o': 'O', 's': 'S', 'p': 'P',  # aromatic
    'b': 'B', 'f': 'F', 'i': 'I',                        # 2-letter: Br, Cl, Se handled separately
}


def _extract_mol_elements(mol: Chem.Mol) -> Set[str]:
    """Extract the set of element symbols present in a molecule (uppercase).

    Used for fast string-level pre-filtering: only rules whose SMARTS
    contain at least one of these elements can possibly match.
    """
    return {atom.GetSymbol() for atom in mol.GetAtoms()}


def _smarts_has_element(smarts: str, element: str) -> bool:
    """Check if a SMARTS string references a specific element.

    Handles both explicit symbols and SMARTS bracket notation:
    - Element 'C': matches 'C' outside brackets, '[C', '[c' (aromatic C)
    - Element 'Br': matches 'Br' outside brackets, '[Br'
    - Element 'Cl': matches 'Cl' outside brackets, '[Cl'
    """
    if not smarts:
        return False
    # Check for '#' (any atom) or '*' (wildcard) — these match anything
    if '#' in smarts or '*' in smarts:
        return True

    # Two-letter elements (check first to avoid false positives)
    if len(element) == 2:
        # Br: check for 'Br' or '[Br'
        return element in smarts or f'[{element}' in smarts

    # One-letter elements
    e_upper = element.upper()
    e_lower = element.lower()

    # Inside bracket notation: [C, [c, [C@], [#6], etc.
    if f'[{e_upper}' in smarts:
        return True
    # Aromatic form inside brackets
    if e_lower != e_upper and f'[{e_lower}' in smarts:
        return True

    # Outside brackets — find bare element symbol
    # Simple approach: check if the uppercase element appears outside brackets
    # A more accurate approach would parse bracket nesting, but for speed
    # we use a heuristic that's good enough for M-CSA SMARTS patterns.
    # We scan the SMARTS string outside brackets:
    i = 0
    in_bracket = False
    n = len(smarts)
    while i < n:
        ch = smarts[i]
        if ch == '[':
            in_bracket = True
            i += 1
            continue
        if ch == ']':
            in_bracket = False
            i += 1
            continue
        if not in_bracket:
            # Check for element symbol
            if ch == e_upper:
                return True
            if ch == e_lower:
                return True
            # Check for 2-letter elements starting at this position
            if i + 1 < n and smarts[i:i+2] == element:
                return True
        i += 1

    return False


def _extract_smarts_elements(smarts: str) -> FrozenSet[str]:
    """Extract element symbols referenced in a SMARTS string.

    Scans the SMARTS pattern and collects element symbols that appear
    outside of brackets (e.g., 'C' in CC(=O)O) and inside brackets
    (e.g., 'C' in [C@H], 'N' in [#7]).

    IMPORTANT: [#6] means carbon (atomic number 6), NOT a wildcard.
    Only bare '#' (without digits following) and '*' are wildcards.
    We parse [#N] notation into actual element symbols.

    Returns a frozenset of uppercase element symbols (or '*' for true wildcards).
    """
    if not smarts:
        return frozenset()

    # Atomic number → element symbol mapping (most common in M-CSA)
    _ATOMIC_NUMBERS = {
        1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O',
        9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P',
        16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 24: 'Cr', 25: 'Mn',
        26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 33: 'As', 34: 'Se',
        35: 'Br', 46: 'Pd', 47: 'Ag', 53: 'I', 55: 'Cs', 56: 'Ba', 78: 'Pt',
        79: 'Au', 80: 'Hg', 82: 'Pb',
    }

    elements = set()
    i = 0
    n = len(smarts)
    while i < n:
        ch = smarts[i]
        if ch == '[':
            # Inside bracket — scan for element symbol
            i += 1
            # Skip optional chirality @, @@
            while i < n and smarts[i] == '@':
                i += 1
            # Skip optional isotope digits
            while i < n and smarts[i].isdigit():
                i += 1
            # Skip optional charge +/- and digits
            while i < n and smarts[i] in '+-':
                i += 1
                while i < n and smarts[i].isdigit():
                    i += 1
            # Now we should be at the element symbol
            if i < n and smarts[i] == '#':
                # Atomic number notation: [#6] = carbon, [#7] = nitrogen, etc.
                i += 1  # skip '#'
                num_str = ''
                while i < n and smarts[i].isdigit():
                    num_str += smarts[i]
                    i += 1
                if num_str:
                    atomic_num = int(num_str)
                    elem = _ATOMIC_NUMBERS.get(atomic_num)
                    if elem:
                        elements.add(elem)
            elif i < n and smarts[i] == '*':
                # True wildcard inside brackets: [*]
                elements.add('*')
                i += 1
            elif i < n and smarts[i].isalpha():
                # Read the element symbol (1-2 chars)
                elem = smarts[i].upper()
                i += 1
                if i < n and smarts[i].islower():
                    elem += smarts[i]  # e.g., 'Br', 'Cl'
                    i += 1
                elements.add(elem)
            # Skip to closing bracket
            while i < n and smarts[i] != ']':
                i += 1
            if i < n:
                i += 1  # skip ']'
        elif ch == '(' or ch == ')' or ch == '.' or ch == '=' or ch == ':' or ch == ';':
            i += 1
        elif ch == '#':
            # Bare # outside brackets — check if followed by digits (ring closure)
            # In SMILES, #C means triple bond to C; but in SMARTS it's unusual
            # Skip it (likely part of ring notation or syntax we don't need)
            i += 1
        elif ch == '*':
            # True wildcard outside brackets
            elements.add('*')
            i += 1
        elif ch == '-' or ch == '+' or ch == '/':
            i += 1
        elif ch.isdigit():
            # Ring closure digit
            i += 1
            if i < n and smarts[i].isdigit():
                i += 1  # two-digit ring closure
        elif ch.isalpha():
            # Element symbol outside brackets
            elem = ch.upper()
            i += 1
            if i < n and smarts[i].islower():
                elem += smarts[i]
                i += 1
            elements.add(elem)
        else:
            i += 1

    return frozenset(elements)


def _build_element_index(
    rules: List[Dict[str, Any]],
) -> Tuple[
    Dict[str, List[int]],
    List[FrozenSet[str]],
]:
    """Build element inverted index for fast rule pre-filtering.

    Returns:
        element_to_indices: element symbol → list of rule indices
        rule_elements: per-rule frozenset of element symbols
    """
    element_to_indices: Dict[str, List[int]] = {}
    rule_elements: List[FrozenSet[str]] = []

    for idx, rule in enumerate(rules):
        smarts = rule.get('reactant_smarts', '')
        elems = _extract_smarts_elements(smarts)
        rule_elements.append(elems)
        for elem in elems:
            if elem not in element_to_indices:
                element_to_indices[elem] = []
            element_to_indices[elem].append(idx)

    return element_to_indices, rule_elements


def _prefilter_rules_by_elements(
    rules: List[Dict[str, Any]],
    mol_elements: Set[str],
    element_to_indices: Dict[str, List[int]],
    rule_elements: List[FrozenSet[str]],
) -> List[Dict[str, Any]]:
    """Fast pre-filter: find rules that share at least one element with the molecule.

    Uses the pre-built inverted index for O(num_elements) lookup instead of
    scanning all rules. Returns the filtered rule list.

    Key optimization: for molecules with rare elements (S, P, F, etc.),
    this dramatically reduces the rule set. For C/H/O/N-only molecules,
    we still get significant reduction because we require rules to match
    at least 2 of the molecule's element types (most rules only involve 1-2).
    """
    if not mol_elements:
        return rules

    # Collect all rule indices that match ANY of the molecule's elements
    candidate_indices: Set[int] = set()

    # Check wildcard rules first (contain '*' in element set)
    wildcard_indices = element_to_indices.get('*', [])
    candidate_indices.update(wildcard_indices)

    for elem in mol_elements:
        indices = element_to_indices.get(elem, [])
        candidate_indices.update(indices)

    # If we still have too many candidates (e.g., C/H/O/N molecule matching
    # most rules), apply a stricter filter: rule must share at least 2 elements
    if len(candidate_indices) > len(rules) * 0.5 and len(mol_elements) >= 2:
        # More than 50% of rules match — apply stricter filter
        # Each rule must share at least min(2, len(mol_elements)) elements
        min_overlap = min(2, len(mol_elements))
        filtered = []
        for idx in candidate_indices:
            overlap = len(rule_elements[idx] & mol_elements)
            if overlap >= min_overlap:
                filtered.append(rules[idx])
        return filtered

    # Return filtered rules (maintain original order)
    return [rules[idx] for idx in sorted(candidate_indices)]


def _compute_mol_to_part_match(
    run_state: PredictionRunState,
    mol: Chem.Mol,
    rule_part: Chem.Mol,
) -> None:
    """Compute substructure matches for a single (mol, rule_part) pair.

    Upload models.py line 1003-1013:
        matches = mol.GetSubstructMatches(rule_part)
        for match in matches:
            self.mol_to_part_matches[mol][rule_part].append(
                {rule_part.GetAtomWithIdx(rule_idx).GetAtomMapNum(): r_idx
                 for rule_idx, r_idx in enumerate(match)})
    """
    if mol in run_state.mol_to_part_matches and rule_part in run_state.mol_to_part_matches[mol]:
        return  # Already computed
    if mol not in run_state.mol_to_part_matches:
        run_state.mol_to_part_matches[mol] = {}

    run_state.mol_to_part_matches[mol][rule_part] = []
    try:
        matches = mol.GetSubstructMatches(rule_part)
        for match in matches:
            # FIX: 对齐原始 models.py (L1012-1013)，包含 atomMapNum==0 的匹配。
            # 过滤 atomMapNum==0 会丢失匹配信息，影响后续 _get_templates_for_duplicated_reactants。
            match_dict = {
                rule_part.GetAtomWithIdx(rule_idx).GetAtomMapNum(): r_idx
                for rule_idx, r_idx in enumerate(match)
            }
            if match_dict:
                run_state.mol_to_part_matches[mol][rule_part].append(match_dict)
    except Exception:
        pass


# ---------------------------------------------------------------------------
# Step 3: get_templates_for_duplicated_reactants
# (Upload models.py line 1027-1082)
# ---------------------------------------------------------------------------
def _get_templates_for_duplicated_reactants(
    run_state: PredictionRunState,
    reactants: Tuple[Chem.Mol, ...],
    rule: Dict[str, Any],
) -> List[Tuple]:
    """Handle duplicate reactants (same mol matching multiple rule parts).

    When a single reactant molecule matches multiple rule parts, we need to
    merge those parts into a combined SMARTS pattern for RunReactants.

    Upload models.py line 1027-1082. Three cases:
    1. Same strict reaction center -> ignore
    2. Completely different reaction centers -> add parentheses
    3. Different but overlapping -> merge_rule_parts

    Returns list of (rxn_template, template_reactants) tuples.
    """
    reactant_to_rule_parts_idx = defaultdict(list)
    for rule_part_idx, reactant in enumerate(reactants):
        reactant_to_rule_parts_idx[reactant].append(rule_part_idx)

    leftside_new_rule_parts = []
    smarts_to_repeated_am_groups = {}

    rule_parts = rule.get('_rule_parts', [])
    reactant_smarts_parts = rule.get('_reactant_smarts_rule_parts', [])
    product_rule_parts = rule.get('_product_rule_parts', [])
    strict_ams = rule.get('_strict_ams', set())

    for reactant, rule_parts_idx in reactant_to_rule_parts_idx.items():
        parts = [rule_parts[r_idx] for r_idx in rule_parts_idx if r_idx < len(rule_parts)]

        if len(parts) == 1:
            r_idx = rule_parts_idx[0]
            if r_idx < len(reactant_smarts_parts):
                leftside_new_rule_parts.append([reactant_smarts_parts[r_idx]])
            else:
                leftside_new_rule_parts.append([])
        else:
            leftside_new_rule_parts.append([])
            # Fix #1: reactants are already AddHs (graph nodes store AddHs mol),
            # so query mol_to_part_matches directly — no _get_mol_h conversion.
            rule_part_matches = [
                run_state.mol_to_part_matches.get(reactant, {}).get(rule_part, [])
                for rule_part in parts
            ]
            if not any(rule_part_matches):
                continue

            for matches_combination in itertools.product(*rule_part_matches):
                if not any(matches_combination):
                    continue
                reactant_idx_to_rule_ams = defaultdict(set)
                for match in matches_combination:
                    for rule_am, reactant_atom_idx in match.items():
                        reactant_idx_to_rule_ams[reactant_atom_idx].add(rule_am)
                repeated_am_groups = [
                    list(ams) for ams in reactant_idx_to_rule_ams.values() if len(ams) > 1
                ]
                repeated_ams = {am for am_group in repeated_am_groups for am in am_group}

                if strict_ams & repeated_ams:
                    continue  # Case 1: same strict reaction center -> skip

                # Case 2 or 3: merge rule parts
                try:
                    new_rule_part = merge_rule_parts(parts, repeated_am_groups)
                    smarts = Chem.MolToSmarts(new_rule_part).replace("&", "")
                    smarts_to_repeated_am_groups[smarts] = repeated_am_groups
                    leftside_new_rule_parts[-1].append(smarts)
                except Exception:
                    pass

    # Build reaction templates from merged parts
    rxn_templates = []
    if not any(leftside_new_rule_parts):
        return rxn_templates

    try:
        for rule_parts_smarts in itertools.product(*leftside_new_rule_parts):
            if not rule_parts_smarts:
                continue
            leftside_smarts = ".".join([
                f"({rp})" if "." in rp else rp for rp in rule_parts_smarts
            ])
            repeated_ams = []
            for smarts in rule_parts_smarts:
                if (repeated := smarts_to_repeated_am_groups.get(smarts)):
                    repeated_ams.extend(repeated)

            if product_rule_parts:
                try:
                    p_rule_parts = merge_rule_parts(product_rule_parts, repeated_ams)
                    rigthside_smarts = Chem.MolToSmarts(p_rule_parts).replace("&", "")
                    reaction_smarts = f"{leftside_smarts}>>{rigthside_smarts}"
                    template = Chem.AllChem.ReactionFromSmarts(reaction_smarts)
                    template.smarts = reaction_smarts
                    rxn_templates.append((template, list(reactant_to_rule_parts_idx.keys())))
                except Exception:
                    pass
    except Exception:
        pass

    return rxn_templates


# ---------------------------------------------------------------------------
# Step 4: Product Validation & Scoring via ReactionResult
# (Upload models.py line 1116-1148)
# ---------------------------------------------------------------------------
def _compute_bond_length_from_products(
    products: List[Chem.Mol],
    rule: Dict[str, Any],
    rule_am_to_product_am: Dict[int, int],
) -> float:
    """Compute max bond length of new bonds using product 3D coordinates.

    Replicates Upload's approach:
    1. Get new bond atom map pairs from rule['_bond_changes_am_heavy']
    2. Map rule atom maps to product atom maps (via old_mapno)
    3. Generate 3D coords for products (with sanitization)
    4. Compute Euclidean distance between bonded atoms

    Important: RDKit RunReactants stores rule atom map numbers in the
    'old_mapno' property. Our ReactionResult._quick_parse reads this
    and sets atom.SetAtomMapNum(rule_am), so product atoms have am > 0.
    """
    from rdkit.Chem import AllChem

    new_bonds_heavy = rule.get('_bond_changes_am_heavy', (set(), set()))[0]
    if not new_bonds_heavy:
        return 0.0

    # Collect which atom maps we need (only new bond endpoints)
    needed_ams: Set[int] = set()
    for bond_ams in new_bonds_heavy:
        for am in bond_ams:
            p_am = rule_am_to_product_am.get(am)
            if p_am is not None:
                needed_ams.add(p_am)

    if not needed_ams:
        return 0.0

    # Build product atom map -> 3D position lookup
    am_to_pos: Dict[int, Tuple[float, float, float]] = {}

    for product in products:
        # Check if this product contains any needed atom maps
        product_has_needed = False
        for atom in product.GetAtoms():
            if atom.GetAtomMapNum() in needed_ams:
                product_has_needed = True
                break
        if not product_has_needed:
            continue

        # Generate 3D coordinates if needed
        if product.GetNumConformers() == 0:
            try:
                mol_copy = Chem.Mol(product)
                Chem.SanitizeMol(mol_copy)
                mol_h = Chem.AddHs(mol_copy)
                # Use EmbedMolecule only (skip MMFFOptimize to avoid
                # RingInfo crash on RunReactants products)
                res = AllChem.EmbedMolecule(mol_h, randomSeed=42, maxAttempts=50)
                if res == -1:
                    # Embedding failed, try without Hs
                    mol_h = mol_copy
                    res = AllChem.EmbedMolecule(mol_h, randomSeed=42, maxAttempts=20)
                if res != -1 and mol_h.GetNumConformers() > 0:
                    # Try MMFF optimization (may fail, non-critical)
                    try:
                        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=100)
                    except Exception:
                        pass
                    conf = mol_h.GetConformer(0)
                    for atom in mol_h.GetAtoms():
                        am = atom.GetAtomMapNum()
                        if am != 0 and atom.GetAtomicNum() > 1 and am in needed_ams:
                            pos = conf.GetAtomPosition(atom.GetIdx())
                            am_to_pos[am] = (pos.x, pos.y, pos.z)
            except Exception:
                pass
        else:
            conf = product.GetConformer(0)
            for atom in product.GetAtoms():
                am = atom.GetAtomMapNum()
                if am != 0 and am in needed_ams:
                    pos = conf.GetAtomPosition(atom.GetIdx())
                    am_to_pos[am] = (pos.x, pos.y, pos.z)

    if not am_to_pos:
        return 0.0

    max_bl = 0.0
    computed_count = 0
    for bond_ams in new_bonds_heavy:
        atoms_ams = sorted(bond_ams)
        am1, am2 = atoms_ams

        # Map rule atom maps to product atom maps
        p_am1 = rule_am_to_product_am.get(am1)
        p_am2 = rule_am_to_product_am.get(am2)

        if p_am1 is None or p_am2 is None:
            continue

        pos1 = am_to_pos.get(p_am1)
        pos2 = am_to_pos.get(p_am2)
        if pos1 is None or pos2 is None:
            continue

        dist = math.sqrt((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2)
        max_bl = max(max_bl, dist)
        computed_count += 1

    logger.debug(
        "Bond length: computed %d/%d bonds, max=%.3fÅ",
        computed_count, len(new_bonds_heavy), max_bl,
    )
    return max_bl


def _get_reaction_data_for_products(
    template_reactants: Tuple[Chem.Mol, ...],
    products: Tuple[Chem.Mol, ...],
    rule: Dict[str, Any],
    run_state: PredictionRunState,
    max_bond_length: float = MAX_BOND_LENGTH,
    coords_cache: Optional[Dict[str, Chem.Mol]] = None,
) -> Optional[Tuple[FrozenSet, Dict[str, Any]]]:
    """Check product validity and compute reaction data using ReactionResult.

    对齐 Upload upload_models.py L1116-1148 的 get_reaction_data_for_products():
      Upload 搜索路径仅做 3 件事（不调用 process_full_pipeline）:
        1. ReactionResult.__init__ → _quick_parse (设置 max_bond_length, rcs_ams)
        2. add_and_get_equivalent_mols (去重 + AtomValenceException 检查)
        3. score = (max(4, max_bond_length) - 3)^2
        4. bonds_away_from_rc discount: if 1 → score *= 0.50

    Upload 的 process_full_pipeline / sanitize / minimise_2d_3d / get_strict_rc /
    has_bond_to_carbon_residue 等函数在搜索期间均未调用（仅在结果可视化时使用）。
    """
    try:
        new_bonds_heavy = rule.get('_bond_changes_am_heavy', (set(), set()))[0]
        rr = ReactionResult(
            reactants=template_reactants,
            products=list(products),
            rule_smarts=rule.get('reaction_smarts', ''),
            am_to_3d_coord=run_state.am_to_3d_coord,
            new_bonds_from_rule=new_bonds_heavy if new_bonds_heavy else None,
        )
    except (KeyError, Exception):
        return None

    # Upload L1135: max_bond_length 由 _quick_parse 设置
    score_bl = rr.max_bond_length

    # NOTE: No hard bond length filter here.
    # Upload uses scoring alone to prune: score = (max(4, bond_length) - 3)^2
    # pushes implausible intermediates (bond_length >> 3Å → score 25-64) to
    # the back of Dijkstra's queue. A hard filter (score_bl > max_bond_length)
    # was previously used when we lacked reliable PDB coordinates, but with the
    # zero-tolerance coordinate policy now in place, all scored bonds come from
    # real PDB coordinates, so the Upload scoring-based pruning is sufficient.

    # Upload L1125-1133: add_and_get_equivalent_mols + AtomValenceException 检查
    # Upload 原始代码:
    #   for product in reaction_result.products:
    #       try:
    #           checked_products.add(self.add_and_get_equivalent_mols(product)[0])
    #           product.UpdatePropertyCache()
    #       except Chem.rdchem.AtomValenceException as e:
    #           print("problem with rule: ", rule.id, e)
    #           ...
    #           return
    checked_products = set()
    for product in rr.products:
        try:
            equiv = run_state.add_and_get_equivalent_mols(product)
            checked_products.add(equiv[0])
            product.UpdatePropertyCache()
        except Chem.rdchem.AtomValenceException as e:
            logger.debug("AtomValenceException for rule %s: %s", rule.get('rule_id', '?'), e)
            return None
        except Exception:
            return None

    # Upload L1137: score = (max(4, max_bond_length) - 3)^2
    # With zero-tolerance coordinate policy, PDB mode results with missing
    # coordinates are already discarded by KeyError in _quick_parse.
    # score_bl=0 now means: either no PDB data (non-PDB mode) or no new bonds
    # to score — both valid cases using the standard formula.
    score = (max(4, score_bl) - 3) ** 2

    # Upload L1139-1141: bonds_away_from_rc discount
    # get_bonds_away_from_rxn_rc 检查:
    #   self.reaction_centres_am & self.predict_mechanism.prediction.step_rc_ams
    # 对齐: run_state.overall_rc_atom_maps 等效于 Upload 的 step_rc_ams
    try:
        bonds_away = rr.get_bonds_away_from_rxn_rc(run_state.overall_rc_atom_maps)
        if bonds_away == 1:
            score *= 0.50
    except Exception:
        pass

    # Upload L1143-1148: 构建 reaction_data
    reaction_data = {
        "max_bond_length": round(score_bl, 2),
        "score": round(score, 2),
        "rcs_ams": list(rr.reaction_centres_am),
    }

    return frozenset(checked_products), reaction_data


# ---------------------------------------------------------------------------
# Step 5: Apply Rule to Configuration
# (Upload models.py line 1151-1190)
# ---------------------------------------------------------------------------
def _apply_rule_to_configuration(
    run_state: PredictionRunState,
    node_id: int,
    rule: Dict[str, Any],
    max_bond_length: float = MAX_BOND_LENGTH,
    coords_cache: Optional[Dict[str, Chem.Mol]] = None,
    _diag: Optional[Dict] = None,
) -> None:
    """Apply a single rule to a single configuration node.

    Upload models.py line 1151-1190:
    1. For each rule part, find which mols in config match
    2. Build reactant_lists via itertools.product
    3. For each combination, call _get_products_for_reactants_and_rule
    4. Create new config node and add edge to graph
    """
    rule_parts = rule.get('_rule_parts', [])
    if not rule_parts:
        return

    reactants_lists = []
    any_part_match = False
    for part in rule_parts:
        mol_matches = []
        for mol in run_state.G.nodes[node_id]["mols"]:
            # Fix #1: Graph nodes now store AddHs mol (with explicit H), same as
            # Upload's Ketcher/PDB molecules. No _get_mol_h conversion needed.
            # This ensures set(combo) uses the SAME Python objects as graph node mols,
            # so the set difference (mols - combo) correctly removes reactants.
            matches = run_state.mol_to_part_matches.get(mol, {}).get(part, [])
            if matches:
                mol_matches.append(mol)
        if mol_matches:
            any_part_match = True
        reactants_lists.append(mol_matches)

    if _diag is not None:
        if any_part_match:
            _diag["any_match"] += 1
        if all(reactants_lists):  # all parts have at least one match
            _diag["all_match"] += 1

    for combo in itertools.product(*reactants_lists):
        products_possibilities = _get_products_for_reactants_and_rule(
            run_state, combo, rule, max_bond_length,
            coords_cache=coords_cache,
            _diag=_diag,
        )
        if products_possibilities is None:
            if _diag is not None:
                _diag["rr_fail"] += 1
            continue

        if _diag is not None and products_possibilities:
            _diag["rr_ok"] += 1

        for products, reaction_data in products_possibilities.items():
            configuration_mols = (run_state.G.nodes[node_id]["mols"] - set(combo)) | products
            new_node_id = run_state.set_and_get_configuration_id(
                tuple(configuration_mols),
                source=f"rule_{rule.get('rule_id', '?')}",
            )

            # Avoid reverse edges (A->B and B->A with same rule)
            # Upload upload_models.py L1169-1179:
            #   ignore = False
            #   for _, d in self.G[node_id][new_node_id].items():
            #       if rule.id == d["rule_id"] and new_bl > d["max_bond_length"]:
            #           ignore = True  # new edge is WORSE → skip
            #   if ignore: continue
            if new_node_id in run_state.G[node_id]:
                if new_node_id == run_state.G[node_id][new_node_id][0].get("start"):
                    continue
                ignore = False
                for _, d in run_state.G[node_id][new_node_id].items():
                    if rule.get('rule_id') == d.get("rule_id") and \
                       reaction_data["max_bond_length"] > d.get("max_bond_length", 0):
                        ignore = True  # new edge is WORSE than existing → skip
                if ignore:
                    continue

            edge_no = run_state.G.add_edge(node_id, new_node_id)
            edge = run_state.G.edges[node_id, new_node_id, edge_no]
            edge["start"] = node_id
            edge["target"] = new_node_id
            edge["edge_no"] = edge_no
            edge["max_bond_length"] = reaction_data["max_bond_length"]
            edge["score"] = reaction_data["score"]
            edge["rule_id"] = rule.get('rule_id', '')
            edge["rule_name"] = rule.get('rule_name', '')
            edge["reaction_smarts"] = rule.get('reaction_smarts', '')
            edge["rcs_ams"] = reaction_data["rcs_ams"]


def _get_products_for_reactants_and_rule(
    run_state: PredictionRunState,
    reactants: Tuple[Chem.Mol, ...],
    rule: Dict[str, Any],
    max_bond_length: float = MAX_BOND_LENGTH,
    coords_cache: Optional[Dict[str, Chem.Mol]] = None,
    _diag: Optional[Dict] = None,
) -> Optional[Dict[Tuple, Dict[str, Any]]]:
    """Get possible products for given reactants and rule (with caching).

    Upload models.py line 1084-1114.
    """
    rule_id = rule.get('rule_id', '')
    # FIX: 使用 tuple 而非 frozenset 作为缓存键。
    # frozenset 会丢失顺序和重复信息：(mol_A, mol_A) → frozenset({mol_A})
    # 这导致同一分子匹配规则不同部分时缓存结果错误。
    # 原始 models.py 使用 tuple（L1086），保持顺序和重复。
    reactants_key = tuple(reactants)

    if rule_id not in run_state.rule_reactants_to_products:
        run_state.rule_reactants_to_products[rule_id] = {}

    if reactants_key in run_state.rule_reactants_to_products[rule_id]:
        return run_state.rule_reactants_to_products[rule_id][reactants_key]

    # Handle duplicate reactants
    if len(reactants) != len(set(reactants)):
        rxn_templates = _get_templates_for_duplicated_reactants(run_state, reactants, rule)
    else:
        rxn_obj = rule.get('_reactant_rxn')
        # Lazy creation of ReactionFromSmarts (Upload source code approach)
        # Created only when a rule actually matches molecules in current config,
        # avoiding OOM from creating 51k Reaction objects upfront.
        if rxn_obj is None:
            full_smarts = rule.get('reaction_smarts', '')
            if full_smarts and '>>' in full_smarts:
                try:
                    from rdkit.Chem import AllChem
                    rxn_obj = AllChem.ReactionFromSmarts(full_smarts)
                    rxn_obj.smarts = full_smarts
                    rule['_reactant_rxn'] = rxn_obj
                except Exception:
                    pass
        if rxn_obj is None:
            run_state.rule_reactants_to_products[rule_id][reactants_key] = {}
            return {}
        rxn_templates = [(rxn_obj, list(reactants))]

    checked_products_possibilities = {}
    rule_parts = rule.get('_rule_parts', [])

    for rxn_template, template_reactants in rxn_templates:
        # Determine maxProducts (Upload models.py L1094-1097)
        # Fix #1: reactants are already AddHs (graph nodes store AddHs mol),
        # so query mol_to_part_matches directly — no _get_mol_h conversion.
        max_products = 1
        for i, reactant in enumerate(reactants):
            if i < len(rule_parts):
                part = rule_parts[i]
                n_matches = len(run_state.mol_to_part_matches.get(reactant, {}).get(part, []))
                if n_matches > 1:
                    max_products = 100  # Reduced from 1000 to prevent iteration stall
                    # Upload source uses 1000 but has more efficient scoring;
                    # 100 is sufficient for mechanism discovery while keeping
                    # per-iteration time reasonable (~10s vs 10+ min with 1000)

        try:
            possible_prods = rxn_template.RunReactants(
                list(template_reactants), maxProducts=max_products
            )
        except Exception:
            continue

        # Dedup products with multiple matches
        if max_products != 1:
            try:
                possible_prods = list({
                    tuple({Chem.MolToSmiles(p) for p in ps}): ps
                    for ps in possible_prods
                }.values())
            except Exception:
                pass

        # [对齐 Upload] 不再调用 filter_products 做预过滤
        # Upload 的 get_products_for_reactants_and_rule (upload_models.py L1084-1114)
        # 仅做 SMILES 去重，然后直接把所有产物传给 get_reaction_data_for_products。
        # 过滤逻辑（sanitize、残基碳）已在 _get_reaction_data_for_products 中处理。
        #
        # [已注释掉] 原 filter_products 调用（创建重复 ReactionResult + 重复 pipeline）：
        # try:
        #     filtered_prods = filter_products(
        #         list(possible_prods),
        #         list(template_reactants),
        #         max_bond_length=max_bond_length,
        #         remove_water_proton=True,
        #     )
        #     if filtered_prods:
        #         possible_prods = [tuple(filtered_prods)]
        #     else:
        #         possible_prods = []
        # except Exception as e:
        #     logger.debug("filter_products failed for rule %s: %s", rule_id, e)

        for products in possible_prods:
            result = _get_reaction_data_for_products(
                template_reactants, products, rule, run_state, max_bond_length,
                coords_cache=coords_cache,
            )
            if result is not None:
                checked_products, reaction_data = result
                checked_products_possibilities[frozenset(checked_products)] = reaction_data
                if _diag is not None:
                    _diag["rd_ok"] += 1
            else:
                if _diag is not None:
                    _diag["rd_fail"] += 1

    run_state.rule_reactants_to_products[rule_id][reactants_key] = checked_products_possibilities
    return checked_products_possibilities if checked_products_possibilities else None


# ---------------------------------------------------------------------------
# Step 6: Dijkstra Exploration (Upload models.py line 1192-1219)
# ---------------------------------------------------------------------------
def _get_next_node_to_explore(
    run_state: PredictionRunState,
) -> Optional[int]:
    """Get next node to explore using Dijkstra path length from both sides.

    Upload models.py line 1192-1219:
    Uses nx.single_source_dijkstra_path_length from BOTH reactant (node 1)
    and product (node 2) sides, alternating between them.

    Returns:
        next_node_id, or None if exhausted.
    """
    if not run_state.to_explore:
        if len(run_state.explored) >= run_state.max_nodes:
            return None

        # Alternating bidirectional exploration (Upload approach)
        n_in_a_row = 1 if len(run_state.explored) < 100 else 10

        start_nodes = (1, 2) if run_state.explore_from_r else (2, 1)
        run_state.explore_from_r = not run_state.explore_from_r

        for start_node in start_nodes:
            if start_node not in run_state.G.nodes:
                continue
            try:
                dijkstra_lengths = nx.single_source_dijkstra_path_length(
                    run_state.G, start_node, weight="score"
                )
                for k in dijkstra_lengths:
                    if k not in run_state.explored:
                        run_state.to_explore.add(k)
                        if len(run_state.to_explore) >= n_in_a_row:
                            break
                else:
                    continue
                break
            except (nx.NetworkXNoPath, nx.NodeNotFound):
                continue

    if run_state.to_explore:
        next_node = run_state.to_explore.pop()
        run_state.explored.add(next_node)
        if next_node in run_state.G.nodes:
            run_state.G.nodes[next_node]["explored"] = True
        return next_node

    return None


# ---------------------------------------------------------------------------
# Main: bidirectional_search - Upload-style Implementation
# ---------------------------------------------------------------------------
def _get_step_arrow_info(reaction_smarts: str) -> Dict[str, Any]:
    """Get step arrow inference info."""
    try:
        from arrow_transform import get_arrow_info_for_step as _get_arrow
        return _get_arrow(reaction_smarts)
    except Exception:
        return {"arrows": [], "arrow_svg": ""}


def bidirectional_search(
    reactants_smiles: str,
    products_smiles: str,
    rules: List[Dict[str, Any]],
    max_configs: int = DEFAULT_MAX_CONFIGS,
    max_depth: int = 10,
    coords_cache: Optional[Dict[str, Chem.Mol]] = None,
    max_bond_length: float = MAX_BOND_LENGTH,
    residues: Optional[List[Dict[str, Any]]] = None,
    pdb_id: str = '',
    pdb_text: str = '',
    progress_file: str = '',
    result_file: str = '',
    mcsa_id: int = 0,
    ligand_mappings: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Upload-style bidirectional mechanism search.

    Adapts Upload models.py PredictionRun.apply_rules() for our system.

    Upload's approach (models.py line 904-1025):
    1. Setup: Parse SMILES -> Mol, add residues, init MultiGraph
    2. Precompute: rule parts, bond changes, strict reaction centers
    3. Search loop:
       a. Get next node via Dijkstra (nx.single_source_dijkstra_path_length)
       b. Precompute mol.GetSubstructMatches for all mols x rule parts
       c. Apply each rule to current config
       d. RunReactants -> ReactionResult -> score -> new config
    4. Export results
    """
    t_start = time.time()
    search_mode = "bidirectional"
    forward_only = not products_smiles.strip()

    # Filter rules by mcsa_id if specified (Own Rules mode)
    if mcsa_id > 0:
        original_count = len(rules)
        rules = [r for r in rules if r.get('mcsa_id') == mcsa_id]
        logger.info("mcsa_id=%d filter: %d → %d rules (Own Rules mode)",
                    mcsa_id, original_count, len(rules))
        if not rules:
            logger.warning("No rules found for mcsa_id=%d", mcsa_id)

    # Build element inverted index for fast rule pre-filtering.
    # Element filter is SAFE: substructure matching requires element presence,
    # so it cannot exclude valid rules.
    # Size filter has been REMOVED (Defect 3 fix): SMARTS [*:1] wildcards
    # are not counted by :N patterns, causing systematic underestimation.
    element_to_indices, rule_elements = _build_element_index(rules)

    if forward_only:
        search_mode = "forward-only"
        logger.info("Forward-only mode (no products specified)")

    # The rule database already contains both forward and reverse entries for each
    # reaction, making the rule pool bidirectionally complete. Substructure matching
    # naturally filters directionality: a rule A+B→C+D matches reactant-side mols,
    # while its counterpart C+D→A+B matches product-side mols. No direction branch
    # is needed in the search loop — matching the original models.py approach.

    # ==================================================================
    # STEP 1: Parse SMILES -> Mol objects (cached)
    # Upload models.py PredictionRun.setup() line 861-900
    # ==================================================================
    run_state = PredictionRunState(max_nodes=max_configs)

    # Parse reactant SMILES -> Mol objects
    reactant_mol_list: List[Chem.Mol] = []
    for smi in reactants_smiles.split('.'):
        smi = smi.strip()
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            logger.warning("Invalid reactant SMILES: %s", smi)
            continue
        reactant_mol_list.append(mol)

    if not reactant_mol_list:
        return _empty_result(reactants_smiles, products_smiles, rules)

    # Add residue side chains as molecules (Upload scheme_mols_reactants includes residues)
    residue_smiles_list = get_residue_smiles(residues)
    res_mol_objects: List[Chem.Mol] = []  # Keep refs for AM assignment & product copy
    for res_smi in residue_smiles_list:
        res_mol = Chem.MolFromSmiles(res_smi)
        if res_mol is not None:
            reactant_mol_list.append(res_mol)
            res_mol_objects.append(res_mol)

    # Assign unique atom map numbers to residue heavy atoms.
    # Upload: residues share the same canvas as substrates, so every atom has
    # a unique AM.  Since users don't draw residues in our system, we assign
    # automatically starting from max(substrate AM) + 1.
    if res_mol_objects:
        max_am = 0
        for mol in reactant_mol_list:
            if mol not in res_mol_objects:
                for atom in mol.GetAtoms():
                    am = atom.GetAtomMapNum()
                    if am > max_am:
                        max_am = am
        for res_mol in res_mol_objects:
            am_counter = max_am + 1
            for atom in res_mol.GetAtoms():
                if atom.GetSymbol() == "H":
                    atom.SetAtomMapNum(0)
                else:
                    atom.SetAtomMapNum(am_counter)
                    am_counter += 1
            max_am = am_counter - 1

    # CRITICAL: Set isotope = atom_map_num on all reactant atoms BEFORE caching.
    # Upload upload_models.py L244-246 (scheme_mols_reactants):
    #   if atom.GetSymbol() == "H": atom.SetAtomMapNum(0)
    #   atom.SetIsotope(atom.GetAtomMapNum())
    # RDKit RunReactants preserves the isotope of matched atoms. After RunReactants,
    # _quick_parse uses atom.GetIsotope() to recover the original reactant AM
    # (which RDKit stores in old_mapno = template AM). Without this, the mapping
    # rule_am → product_am is broken and ALL downstream logic fails:
    #   - get_strict_reaction_centres returns wrong symmetric difference
    #   - get_bonds_away_from_rxn_rc uses wrong AMs for intersection
    #   - Upload path bond length cannot find coords in am_to_3d_coord
    for mol in reactant_mol_list:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                atom.SetAtomMapNum(0)
            atom.SetIsotope(atom.GetAtomMapNum())

    # CRITICAL FIX #1: AddHs BEFORE add_and_get_equivalent_mols.
    # Upload's molecules come from Ketcher/PDB with explicit H already present.
    # Our molecules come from MolFromSmiles (implicit H). We must AddHs now so
    # that graph nodes store the SAME mol object used for substructure matching
    # and RunReactants. This ensures set(combo) uses the same Python objects as
    # graph node mols, so the set difference (mols - combo) correctly removes
    # reactants — fixing the 0-paths bug where configurations kept all molecules
    # and never converged.
    #
    # Order matters: AddHs must come AFTER isotope=AM is set, because AddHs
    # adds new H atoms that inherit isotope=0 from the implicit-H state, and
    # we want H atoms to have isotope=0 (since their AM is 0).
    for i, mol in enumerate(reactant_mol_list):
        reactant_mol_list[i] = Chem.AddHs(mol)
    # Update res_mol_objects refs (they are the same list tail)
    res_mol_objects = reactant_mol_list[len(reactant_mol_list) - len(res_mol_objects):]

    # Equivalence-normalize substrates only; residues bypass dedup
    # (same SMILES at different PDB positions need separate mol instances).
    substrate_mol_list = reactant_mol_list[:len(reactant_mol_list) - len(res_mol_objects)]
    run_state.initial_mols = run_state.add_and_get_equivalent_mols(*substrate_mol_list)
    for res_mol in res_mol_objects:
        run_state.initial_mols.append(res_mol)
    r_node_id = run_state.set_and_get_configuration_id(
        tuple(run_state.initial_mols), is_initial_reactant=True
    )

    # --- Part D: Build am_to_3d_coord from PDB coordinates ---
    # Upload upload_models.py L816-825:
    #   am_to_3d_coord = {atom.GetAtomMapNum(): conformer_3d.GetAtomPosition(atom.GetIdx())}
    #   for each mol in initial_mols, using conformer 1 (PDB 3D).
    #
    # When PDB coordinates are available, inject them as conformers on
    # reactant mols and build am_to_3d_coord for Upload-style bond length computation.
    # Trigger on either coords_cache having real coordinates OR ligand_mappings being provided
    # (ligand_mappings are processed by _build_am_to_3d_coord's per-mapping "consume" mode).
    has_pdb_data = (
        (coords_cache and any(v is not None and v.GetNumConformers() > 0 for v in coords_cache.values()))
        or (ligand_mappings and pdb_text)
    )
    if has_pdb_data:
        _build_am_to_3d_coord(run_state, coords_cache, reactants_smiles,
                              residues=residues, pdb_text=pdb_text,
                              residue_smiles_list=residue_smiles_list,
                              ligand_mappings=ligand_mappings)
        if run_state.am_to_3d_coord:
            logger.info("am_to_3d_coord: %d atoms with PDB 3D coordinates (Upload path enabled)",
                        len(run_state.am_to_3d_coord))

            # --- Coordinate completeness check (zero tolerance for ALL ligands) ---
            # Aligned with original EzMechanism: in PDB mode, ALL molecules (including
            # water) must have complete PDB coordinates. The original relies on Django
            # DB to guarantee mapping completeness; our standalone version enforces this
            # explicitly at startup. Without water coordinates, product-side edges
            # involving water atoms are silently discarded by _quick_parse (KeyError),
            # causing bidirectional search imbalance (product frontier starved).
            for mol in run_state.initial_mols:
                is_water = mol_is_water_or_proton(mol)

                # ALL molecules: must have ALL non-H atoms in am_to_3d_coord
                missing_ams = []
                for atom in mol.GetAtoms():
                    if atom.GetAtomicNum() > 1 and atom.GetAtomMapNum() > 0:
                        if atom.GetAtomMapNum() not in run_state.am_to_3d_coord:
                            missing_ams.append(atom.GetAtomMapNum())
                if missing_ams:
                    mol_smi = Chem.MolToSmiles(Chem.RemoveHs(mol)) if mol.GetNumAtoms() > 0 else '?'
                    mol_type = "Water molecule" if is_water else "Core ligand"
                    logger.error("COORDINATE INTEGRITY CHECK FAILED: %s %s has %d atoms "
                                 "missing PDB coordinates (AMs: %s). Cannot start search in PDB mode. "
                                 "%s",
                                 mol_type, mol_smi[:50], len(missing_ams), str(missing_ams[:10]),
                                 "Please map all water molecules to specific HOH residues in the PDB structure."
                                 if is_water else "Ensure all ligands have proper PDB residue mappings.")
                    # Return empty result — cannot search without complete coordinates
                    empty_result = {
                        "status": "error",
                        "message": (f"{mol_type} {mol_smi[:30]} has {len(missing_ams)} atoms "
                                    f"missing PDB coordinates. "
                                    + ("Please map all water molecules to specific HOH residues in the PDB structure."
                                       if is_water else
                                       "Ensure all ligands have proper PDB residue mappings.")),
                        "search_path": [],
                        "total_configs_explored": 0,
                        "pdb_id": pdb_id,
                        "mcsa_id": mcsa_id,
                    }
                    if progress_file:
                        try:
                            with open(progress_file, 'w') as f:
                                json.dump({"status": "error", "progress": 0, "message": empty_result["message"]}, f)
                        except Exception:
                            pass
                    if result_file:
                        try:
                            with open(result_file, 'w') as f:
                                json.dump(empty_result, f)
                        except Exception:
                            pass
                    return empty_result

            logger.info("Coordinate completeness check passed: all ligands (including water) have PDB coordinates")
        else:
            logger.info("am_to_3d_coord: no PDB coordinates available (fallback path)")

    # Parse product SMILES -> Mol objects
    # Separate substrate mols from residue mols (same-SMILES-different-position
    # residues must bypass equivalence dedup, same as reactant side).
    product_substrate_list: List[Chem.Mol] = []
    if products_smiles.strip():
        for smi in products_smiles.split('.'):
            smi = smi.strip()
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                logger.warning("Invalid product SMILES: %s", smi)
                continue
            product_substrate_list.append(mol)

    # Deep copy residue mols for product side (preserves AM + conformers)
    product_residue_list = [Chem.Mol(res_mol) for res_mol in res_mol_objects]

    # Set isotope = atom_map_num on product substrates (for backward search)
    for mol in product_substrate_list:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                atom.SetAtomMapNum(0)
            atom.SetIsotope(atom.GetAtomMapNum())

    # CRITICAL FIX #1 (product side): AddHs BEFORE add_and_get_equivalent_mols.
    # Same reason as reactant side — product graph nodes must store AddHs mol
    # so that backward search substructure matching and set difference work.
    for i, mol in enumerate(product_substrate_list):
        product_substrate_list[i] = Chem.AddHs(mol)

    # Residue mols already have isotope set and AddHs from reactant side;
    # deep copy preserves both.

    # Equivalence-normalize product substrates only; residues bypass dedup
    p_node_id = None
    if product_substrate_list or product_residue_list:
        initial_prods = run_state.add_and_get_equivalent_mols(*product_substrate_list)
        for res_mol in product_residue_list:
            initial_prods.append(res_mol)
        run_state.initial_product_mols = initial_prods

        # Upload setup() L892-899: 复制反应物 3D 坐标到产物分子 conformer 1
        # 作用: 让产物分子也有 PDB 参考坐标，后续产物可视化时可以直接使用
        if run_state.am_to_3d_coord:
            for mol in run_state.initial_product_mols:
                if len(mol.GetConformers()) == 1:
                    mol.AddConformer(mol.GetConformer(), assignId=True)
                    conformer_3d = mol.GetConformer(1)
                    conformer_3d.Set3D(True)
                    for atom in mol.GetAtoms():
                        if atom.GetSymbol() != "H":
                            am = atom.GetAtomMapNum()
                            if am != 0 and am in run_state.am_to_3d_coord:
                                try:
                                    conformer_3d.SetAtomPosition(
                                        atom.GetIdx(),
                                        run_state.am_to_3d_coord[am],
                                    )
                                except Exception:
                                    pass

        p_node_id = run_state.set_and_get_configuration_id(
            tuple(run_state.initial_product_mols), is_initial_product=True
        )

    logger.info(
        "Setup: %d reactant mols, %d product mols, %d residue mols, %d rules",
        len(reactant_mol_list), len(product_substrate_list) + len(product_residue_list),
        len(residue_smiles_list), len(rules),
    )

    # Log initial molecule SMILES for diagnostics
    for i, mol in enumerate(run_state.initial_mols):
        try:
            smi = Chem.MolToSmiles(mol)
            n_atoms = mol.GetNumAtoms()
            logger.info("  initial_mol[%d]: %s  (%d atoms)", i, smi[:120], n_atoms)
        except Exception:
            pass
    if run_state.initial_product_mols:
        for i, mol in enumerate(run_state.initial_product_mols):
            try:
                smi = Chem.MolToSmiles(mol)
                n_atoms = mol.GetNumAtoms()
                logger.info("  product_mol[%d]: %s  (%d atoms)", i, smi[:120], n_atoms)
            except Exception:
                pass

    # ==================================================================
    # STEP 1.5: Identify overall reaction center (Part C)
    # Upload's equivalent: Prediction.step_rc_ams (upload_models.py L700-706)
    #
    # Uses atom-mapped SMILES from Ketcher RXN (or user input) to compute
    # reaction center atom map numbers via bond-level symmetric difference.
    # This is our equivalent of the paper's step_rc_ams which comes from
    # Django DB (admin-annotated from M-CSA).
    #
    # These atom map numbers are stored in run_state and used later in
    # _get_reaction_data_for_products() for the score × 0.50 discount
    # (upload_models.py L1139-1141).
    # ==================================================================
    if products_smiles.strip():
        run_state.overall_rc_atom_maps = identify_reaction_center_atom_maps(
            reactants_smiles, products_smiles
        )
        if run_state.overall_rc_atom_maps:
            logger.info(
                "Reaction center: %d atom map numbers identified (step_rc_ams equivalent)",
                len(run_state.overall_rc_atom_maps),
            )
        else:
            logger.info("Reaction center: no atom-mapped SMILES available, RC scoring disabled")

    # ==================================================================
    # STEP 2: Lazy rule parsing (no upfront preparse — avoids OOM)
    #
    # Instead of calling _lightweight_preparse_rules on all 51k rules
    # (which creates MolFromSmarts for every rule and needs 4.2GB+),
    # we lazily parse each rule's SMARTS only when the search loop
    # first encounters it via _lazy_parse_rule().
    # ==================================================================
    t_preparse = time.time()
    preparse_time = 0.0  # Will be measured inline as rules are parsed
    logger.info(
        "Rules: %d loaded, element index: %d elements, lazy parsing enabled",
        len(rules), len(element_to_indices),
    )

    # ==================================================================
    # STEP 2.5: Initialize search
    # Upload models.py apply_rules() line 996-998
    # ==================================================================
    if p_node_id is not None:
        run_state.to_explore = {r_node_id, p_node_id}
    else:
        run_state.to_explore = {r_node_id}
    run_state.explored = set()

    # ==================================================================
    # STEP 3: Main search loop (Upload-style Dijkstra)
    # Upload models.py apply_rules() line 999-1020
    # ==================================================================
    iteration = 0
    meeting_points_found = 0

    # Write initial progress for async mode
    if progress_file:
        write_progress(progress_file, "RUNNING",
                       explored_nodes=0, total_nodes=run_state.G.number_of_nodes(),
                       current_iteration=0, elapsed_seconds=0.0)

    while True:
        next_node_id = _get_next_node_to_explore(run_state)
        if next_node_id is None:
            break
        iteration += 1
        if next_node_id not in run_state.G.nodes:
            continue

        node_mols = run_state.G.nodes[next_node_id].get("mols", set())
        logger.info(
            "Iter %d: node %d (%d mols, explored=%d, nodes=%d)",
            iteration, next_node_id,
            len(node_mols),
            len(run_state.explored), run_state.G.number_of_nodes(),
        )

        # Upload source code approach (apply_rules line 999-1020):
        # For each rule, lazily parse its rule_parts if not already done,
        # then compute mol-to-part matches inline.
        #
        # IMPORTANT: M-CSA rule SMARTS reference explicit H atoms (e.g., [#1+0:2]).
        # The graph stores de-H'd molecules (via make_mol_equivalent), so we
        # must AddHs before substructure matching.
        nodes_before = run_state.G.number_of_nodes()
        edges_before = run_state.G.number_of_edges()
        diag = {"any_match": 0, "all_match": 0, "rr_ok": 0, "rr_fail": 0,
                "rd_ok": 0, "rd_fail": 0}

        # ---- PRE-FILTER RULES (element types only) ----
        # Element filter: rules must share elements with molecule (inverted index, fast)
        # This is SAFE: substructure matching requires element presence.
        # The rule pool is bidirectionally complete — no direction branch needed.
        t_filter_start = time.time()
        all_mol_elements: Set[str] = set()
        for mol in node_mols:
            all_mol_elements.update(_extract_mol_elements(mol))

        # Element-based pre-filter (uses inverted index) → set of indices
        candidate_indices: Set[int] = set()
        # Wildcard rules
        candidate_indices.update(element_to_indices.get('*', []))
        for elem in all_mol_elements:
            candidate_indices.update(element_to_indices.get(elem, []))
        # If too many matches (C/H/O/N molecule), require >=2 element overlap
        if len(candidate_indices) > len(rules) * 0.5 and len(all_mol_elements) >= 2:
            min_overlap = min(2, len(all_mol_elements))
            filtered_idx: Set[int] = set()
            for idx in candidate_indices:
                if len(rule_elements[idx] & all_mol_elements) >= min_overlap:
                    filtered_idx.add(idx)
            candidate_indices = filtered_idx

        # Convert indices to rule list (sorted for deterministic order)
        iter_rules = [rules[idx] for idx in sorted(candidate_indices)]
        filter_time = time.time() - t_filter_start

        # Use batch processing when the FILTERED rule count is large.
        # With pre-filtering, this is typically 1-5K rules, well within memory limits.
        # Batch mode is only needed as a safety net for edge cases.
        use_batches = len(iter_rules) > 8000

        if use_batches:
            # Fix #1: node_mols are already AddHs, use them directly as cache keys
            current_mol_keys = list(node_mols)

        # Process filtered rules (in batches if needed for memory safety)
        total_filtered = len(iter_rules)
        step = RULE_BATCH_SIZE if use_batches else total_filtered
        batch_idx = 0
        for batch_start in range(0, total_filtered, step):
            batch_end = min(batch_start + step, total_filtered)

            for rule in iter_rules[batch_start:batch_end]:
                # Lazy parse rule parts on first encounter
                if '_rule_parts' not in rule:
                    _lazy_parse_rule(rule)

                rule_parts = rule.get('_rule_parts', [])
                if not rule_parts:
                    continue

                # Compute mol-to-part matches for this rule's parts
                # Fix #1: node_mols are already AddHs, use them directly
                for mol in node_mols:
                    for part in rule_parts:
                        _compute_mol_to_part_match(run_state, mol, part)

                # Apply the rule to this configuration
                _apply_rule_to_configuration(
                    run_state, next_node_id, rule, max_bond_length,
                    coords_cache=coords_cache,
                    _diag=diag,
                )

            # Batch cleanup: free RDKit objects to prevent OOM
            if use_batches:
                for rule in iter_rules[batch_start:batch_end]:
                    _release_rule_rdkit_objects(rule)

                # Clear stale cache entries (O(num_mols) — just pop entire inner dicts)
                for mol_key in current_mol_keys:
                    run_state.mol_to_part_matches.pop(mol_key, None)

                # Periodic gc.collect — every GC_INTERVAL batches to amortize cost
                batch_idx += 1
                if batch_idx % GC_INTERVAL == 0:
                    gc.collect()

        nodes_after = run_state.G.number_of_nodes()
        edges_after = run_state.G.number_of_edges()
        logger.info(
            "Iter %d DIAG: filtered=%d/%d rules (elems=%s, %.1fms) "
            "any_match=%d all_match=%d "
            "RunReactants_ok=%d RunReactants_fail=%d reaction_data_ok=%d reaction_data_fail=%d "
            "new_nodes=%d new_edges=%d",
            iteration, total_filtered, len(rules),
            sorted(all_mol_elements),
            filter_time * 1000,
            diag["any_match"], diag["all_match"],
            diag["rr_ok"], diag["rr_fail"],
            diag["rd_ok"], diag["rd_fail"],
            nodes_after - nodes_before, edges_after - edges_before,
        )

        run_state.total_expanded += 1

        # Write progress after each iteration for async polling
        if progress_file:
            write_progress(progress_file, "RUNNING",
                           explored_nodes=len(run_state.explored),
                           total_nodes=run_state.G.number_of_nodes(),
                           current_iteration=iteration,
                           elapsed_seconds=time.time() - t_start)

        # Check for meeting points (bidirectional only) — logging only, no early termination
        # Upload source code runs until _get_next_node_to_explore() returns None or max_nodes reached
        if p_node_id is not None:
            if nx.has_path(run_state.G, r_node_id, p_node_id):
                try:
                    path_len = nx.shortest_path_length(
                        run_state.G, r_node_id, p_node_id
                    )
                    if path_len <= max_depth:
                        meeting_points_found += 1
                        logger.info(
                            "Meeting point! path_len=%d (iter=%d)",
                            path_len, iteration,
                        )
                except nx.NetworkXNoPath:
                    pass

    search_time = time.time() - t_start
    logger.info(
        "Done: %d iters, %d nodes, %d edges, %.2fs",
        iteration, run_state.G.number_of_nodes(),
        run_state.G.number_of_edges(), search_time,
    )

    # ==================================================================
    # STEP 4: Build results
    # ==================================================================
    result = _build_upload_result(
        run_state=run_state,
        reactants_smiles=reactants_smiles,
        products_smiles=products_smiles,
        rules=rules,
        r_node_id=r_node_id,
        p_node_id=p_node_id,
        search_mode=search_mode,
        search_time=search_time,
        meeting_points_found=meeting_points_found,
        max_bond_length=max_bond_length,
        original_reactant_mols=reactant_mol_list,
        original_product_mols=product_substrate_list + product_residue_list if (product_substrate_list or product_residue_list) else None,
    )
    result["preparse_time"] = round(preparse_time, 3)

    # Write final progress + result file for async mode
    if progress_file:
        if result_file:
            try:
                with open(result_file, 'w') as f:
                    json.dump(result, f)
            except Exception as e:
                logger.warning("Failed to write result file %s: %s", result_file, e)
        write_progress(progress_file, "DONE",
                       explored_nodes=len(run_state.explored),
                       total_nodes=run_state.G.number_of_nodes(),
                       current_iteration=iteration,
                       elapsed_seconds=search_time,
                       result_file=result_file)

    return result


def _random_label_suffix(length: int = 3) -> str:
    """Generate a random short letter string for config labels."""
    import random, string
    return ''.join(random.choices(string.ascii_lowercase, k=length))


def _build_upload_result(
    run_state: PredictionRunState,
    reactants_smiles: str,
    products_smiles: str,
    rules: List[Dict[str, Any]],
    r_node_id: int,
    p_node_id: Optional[int],
    search_mode: str,
    search_time: float,
    meeting_points_found: int,
    max_bond_length: float,
    original_reactant_mols: Optional[List[Chem.Mol]] = None,
    original_product_mols: Optional[List[Chem.Mol]] = None,
) -> Dict[str, Any]:
    """Build result JSON from Upload-style MultiGraph."""
    nodes_list = []
    edges_list = []

    # Arrow info cache: reaction_smarts -> {"arrows": [...], "arrow_svg": "..."}
    # Avoids recomputing arrow info for the same rule across edges and steps
    _arrow_info_cache: Dict[str, Dict[str, Any]] = {}

    # Pre-compute SMILES for original input molecules (before dedup)
    # so reactant/product nodes show ALL input molecules, not just unique ones.
    original_reactant_smiles = sorted(
        Chem.MolToSmiles(m) for m in original_reactant_mols if m is not None
    ) if original_reactant_mols else []
    original_product_smiles = sorted(
        Chem.MolToSmiles(m) for m in original_product_mols if m is not None
    ) if original_product_mols else []

    # =====================================================================
    # Pre-compute global Dijkstra distances (2 calls instead of per-edge)
    # Used for: node distance labeling, edge direction, step direction.
    # Complexity: O(V log V + E) instead of O(E × V log V).
    # =====================================================================
    _INF = float('inf')
    dist_from_r: Dict[int, float] = {}
    dist_from_p: Dict[int, float] = {}
    try:
        dist_from_r = nx.single_source_dijkstra_path_length(
            run_state.G, r_node_id, weight="score"
        )
    except Exception:
        pass
    if p_node_id is not None:
        try:
            dist_from_p = nx.single_source_dijkstra_path_length(
                run_state.G, p_node_id, weight="score"
            )
        except Exception:
            pass

    # Node distances from reactant for labeling (unweighted hop count)
    node_distances: Dict[int, int] = {}
    try:
        lengths = nx.single_source_shortest_path_length(run_state.G, r_node_id)
        node_distances = dict(lengths)
    except Exception:
        pass

    # =====================================================================
    # Path extraction: single optimal paths via Dijkstra
    # =====================================================================
    # Instead of enumerating all simple paths (exponential, times out on
    # large graphs), we compute two deterministic paths:
    #
    # 1. Score-optimal path: dijkstra_path(weight="score")
    #    - Minimises total score sum → best 3D-geometry mechanism
    #    - O(V log V + E), completes in < 1s even on 11K-node graphs
    #
    # 2. Step-minimal path: shortest_path() (unweighted)
    #    - Minimises edge count → fewest reaction steps
    #    - O(V + E), sub-second
    #
    # If both paths are identical, only one is kept.  This matches the
    # original EzMechanism behaviour (single best path) while providing
    # an optional comparison when the two criteria disagree.
    # =====================================================================
    all_paths: List[List[int]] = []   # list of node-ID lists
    path_labels: List[str] = []       # "score_optimal" or "step_minimal"
    rule_path_count: Dict[str, int] = {}  # rule_id -> count across paths (for edge coloring)

    if p_node_id is not None and p_node_id in dist_from_r:
        # 1) Score-optimal path (weighted Dijkstra)
        try:
            best_path = nx.dijkstra_path(
                run_state.G, r_node_id, p_node_id, weight="score"
            )
            all_paths.append(best_path)
            path_labels.append("score_optimal")
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            pass

        # 2) Step-minimal path (unweighted BFS) — only if different from #1
        try:
            shortest_path = nx.shortest_path(run_state.G, r_node_id, p_node_id)
            if not all_paths or shortest_path != all_paths[0]:
                all_paths.append(shortest_path)
                path_labels.append("step_minimal")
        except (nx.NetworkXNoPath, nx.NodeNotFound):
            pass

        # Collect path node IDs + count rules across paths
        path_node_ids: Set[int] = set()
        for path_keys in all_paths:
            path_node_ids.update(path_keys)
            for i in range(len(path_keys) - 1):
                u, v = path_keys[i], path_keys[i + 1]
                edge_data = run_state.G.get_edge_data(u, v)
                if edge_data:
                    edge_info = edge_data.get(0, edge_data) if isinstance(edge_data, dict) else edge_data[0]
                    rid = str(edge_info.get("rule_id", ""))
                    if rid:
                        rule_path_count[rid] = rule_path_count.get(rid, 0) + 1

        logger.info(
            "Path extraction: %d path(s) found [%s]",
            len(all_paths),
            ", ".join(path_labels),
        )
    else:
        path_node_ids: Set[int] = set()
        if r_node_id is not None:
            path_node_ids.add(r_node_id)
        if p_node_id is not None:
            path_node_ids.add(p_node_id)

    # ---------------------------------------------------------------------------
    # Phase A: Optimize 2D coordinates + detect stereochemistry for path nodes
    # ---------------------------------------------------------------------------
    # The search produces intermediate Mol objects with raw 2D coordinates from
    # RunReactants. These are often distorted (overlapping atoms, bad angles).
    # optimize_mol_2d() re-computes clean 2D layouts via Compute2DCoords.
    # detect_mol_stereo() assigns chiral tags from 3D conformers when available.
    #
    # Safety guards:
    # - optimize_mol_2d: only called when NumConformers == 1 (won't touch 3D)
    # - detect_mol_stereo: only called when a true 3D conformer exists (confId=1)
    # - Both return bool; failures are logged but never crash the result build
    # - Only path nodes are optimized (non-path intermediates are left as-is)
    # ---------------------------------------------------------------------------
    _opt_2d_ok = 0
    _opt_2d_skip = 0
    _stereo_ok = 0
    _stereo_skip = 0
    for node_id in path_node_ids:
        mols_set = run_state.G.nodes.get(node_id, {}).get("mols", set())
        for mol in mols_set:
            if mol is None:
                continue
            # 2D optimization: safe for NumConformers == 1 molecules only
            if optimize_mol_2d(mol):
                _opt_2d_ok += 1
            else:
                _opt_2d_skip += 1
            # Stereo detection: only succeeds for molecules with 3D conformer
            if detect_mol_stereo(mol):
                _stereo_ok += 1
            else:
                _stereo_skip += 1
    if _opt_2d_ok > 0 or _stereo_ok > 0:
        logger.debug(
            "_build_upload_result Phase A: 2d_opt=%d/%d, stereo=%d/%d",
            _opt_2d_ok, _opt_2d_ok + _opt_2d_skip,
            _stereo_ok, _stereo_ok + _stereo_skip,
        )

    # Export all nodes
    for node_id, data in run_state.G.nodes(data=True):
        mols_set = data.get("mols", set())
        mol_smiles = []
        for mol in mols_set:
            try:
                mol_smiles.append(Chem.MolToSmiles(mol))
            except Exception:
                pass

        # For reactant/product nodes, use the original input molecules
        # (preserves duplicates like 2×SER, 2×LYS, 2×water)
        if node_id == r_node_id and original_reactant_smiles:
            display_mols = original_reactant_smiles
        elif node_id == p_node_id and original_product_smiles:
            display_mols = original_product_smiles
        else:
            display_mols = sorted(mol_smiles)

        # Determine node source category
        if node_id == r_node_id:
            node_source = "reactant"
            display_label = "R"
        elif node_id == p_node_id:
            node_source = "product"
            display_label = "P"
        elif node_id in path_node_ids:
            # Config is on at least one found path = exploration config
            node_source = "exploration"
            dist = node_distances.get(node_id, 0)
            display_label = f"{dist}_{_random_label_suffix()}"
        else:
            # Generated but not on any path
            node_source = "intermediate"
            dist = node_distances.get(node_id, 0)
            display_label = f"{dist}_{_random_label_suffix()}"

        nodes_list.append({
            "id": f"node_{node_id}",
            "label": display_label,
            "smiles_label": ", ".join(display_mols[:3]) + (
                f" (+{len(display_mols)-3} more)" if len(display_mols) > 3 else ""
            ),
            "molecules": display_mols,
            "source": node_source,
            "depth": node_distances.get(node_id, data.get("depth", -1)),
        })

    # Export all edges
    for u, v, data in run_state.G.edges(data=True):
        # Determine edge direction: forward if u closer to reactant, backward if closer to product
        # Uses pre-computed dist_from_r / dist_from_p dictionaries (O(1) lookup)
        direction = "forward"
        if p_node_id is not None and u == p_node_id:
            direction = "reverse"
        elif p_node_id is not None:
            dist_u_fwd = dist_from_r.get(u, _INF) if u != r_node_id else 0
            dist_v_fwd = dist_from_r.get(v, _INF) if v != r_node_id else 0
            dist_u_bwd = dist_from_p.get(u, _INF) if u != p_node_id else 0
            dist_v_bwd = dist_from_p.get(v, _INF) if v != p_node_id else 0
            if dist_u_bwd < dist_u_fwd or dist_v_bwd < dist_v_fwd:
                direction = "reverse"

        # Determine rule_status for edge coloring:
        # - "cross_mechanism": rule appears in multiple paths (orange)
        # - "database_only": rule appears in only one path, simulating db-only (red)
        # - "normal": default (gray)
        rid = str(data.get("rule_id", ""))
        rule_count = rule_path_count.get(rid, 0)
        if rule_count > 1:
            rule_status = "cross_mechanism"
        elif rule_count == 1:
            rule_status = "database_only"
        else:
            rule_status = "normal"

        # Arrow info for this edge (cached by reaction_smarts)
        arrow_info = _arrow_info_cache.get(data.get("reaction_smarts", ""))
        if arrow_info is None:
            arrow_info = _get_step_arrow_info(data.get("reaction_smarts", ""))
            _arrow_info_cache[data.get("reaction_smarts", "")] = arrow_info

        edges_list.append({
            "id": f"edge_{len(edges_list)}",
            "source": f"node_{u}",
            "target": f"node_{v}",
            "rule_id": rid,
            "rule_name": data.get("rule_name", ""),
            "reaction_smarts": data.get("reaction_smarts", ""),
            "direction": direction,
            "label": data.get("rule_name", ""),
            "score": round(data.get("score", 0), 3),
            "max_bond_length": round(data.get("max_bond_length", 0), 3),
            "rule_status": rule_status,
            "arrows": arrow_info.get("arrows", []),
            "arrow_svg": arrow_info.get("arrow_svg", ""),
        })

    # Build path output from extracted paths (1 or 2 paths max)
    paths = []
    if p_node_id is not None and all_paths:
        try:
            for pidx, path_keys in enumerate(all_paths):
                steps = []
                total_score = 0.0
                for i in range(len(path_keys) - 1):
                    from_key, to_key = path_keys[i], path_keys[i + 1]
                    edge_data = run_state.G.get_edge_data(from_key, to_key)
                    if edge_data:
                        edge_info = (
                            edge_data.get(0, edge_data)
                            if isinstance(edge_data, dict)
                            else edge_data[0]
                        )
                        # Build from_config / to_config from graph node data
                        # IMPORTANT: G.nodes stores Mol objects — must convert to SMILES for JSON
                        # For reactant/product nodes, use original input (preserves duplicates)
                        from_node_data = run_state.G.nodes.get(from_key, {})
                        to_node_data = run_state.G.nodes.get(to_key, {})
                        if from_key == r_node_id and original_reactant_smiles:
                            from_mols = list(original_reactant_smiles)
                        else:
                            from_mols = sorted(
                                Chem.MolToSmiles(mol) for mol in from_node_data.get("mols", set())
                                if mol is not None
                            )
                        if to_key == p_node_id and original_product_smiles:
                            to_mols = list(original_product_smiles)
                        else:
                            to_mols = sorted(
                                Chem.MolToSmiles(mol) for mol in to_node_data.get("mols", set())
                                if mol is not None
                            )
                        # Determine step direction based on from_node position
                        # Uses pre-computed dist_from_r / dist_from_p (O(1) lookup)
                        step_direction = "forward"
                        if p_node_id is not None:
                            d_r = dist_from_r.get(from_key, _INF) if from_key != r_node_id else 0
                            d_p = dist_from_p.get(from_key, _INF) if from_key != p_node_id else 0
                            if d_p < d_r:
                                step_direction = "reverse"
                        # Arrow info for this step (cached by reaction_smarts)
                        step_smarts = edge_info.get("reaction_smarts", "")
                        step_arrow_info = _arrow_info_cache.get(step_smarts)
                        if step_arrow_info is None:
                            step_arrow_info = _get_step_arrow_info(step_smarts)
                            _arrow_info_cache[step_smarts] = step_arrow_info

                        steps.append({
                            "step_num": i + 1,
                            "from_config": {
                                "molecules": from_mols,
                                "source": from_node_data.get("source", "unknown"),
                            },
                            "to_config": {
                                "molecules": to_mols,
                                "source": to_node_data.get("source", "unknown"),
                            },
                            "rule_id": edge_info.get("rule_id", ""),
                            "rule_name": edge_info.get("rule_name", ""),
                            "reaction_smarts": edge_info.get("reaction_smarts", ""),
                            "direction": step_direction,
                            "score": round(edge_info.get("score", 0), 3),
                            "max_bond_length": round(edge_info.get("max_bond_length", 0), 3),
                            "arrows": step_arrow_info.get("arrows", []),
                            "arrow_svg": step_arrow_info.get("arrow_svg", ""),
                        })
                        total_score += edge_info.get("score", 0)

                # Determine path_type and path_label for this path
                path_type = path_labels[pidx] if pidx < len(path_labels) else "score_optimal"
                if path_type == "score_optimal":
                    path_label = "Score-Optimal Path"
                elif path_type == "step_minimal":
                    path_label = "Shortest-Step Path"
                else:
                    path_label = f"Path {pidx + 1}"

                paths.append({
                    "path_id": pidx,
                    "path_type": path_type,
                    "path_label": path_label,
                    "steps": steps,
                    "total_score": round(total_score, 3),
                    "total_steps": len(steps),
                    "num_steps": len(steps),
                })
        except Exception as e:
            logger.error("Failed to build path steps: %s", e, exc_info=True)
    elif search_mode == "forward-only":
        # Forward-only: each edge from reactant node is a single-step path
        if r_node_id in run_state.G:
            for fi, neighbor in enumerate(run_state.G.neighbors(r_node_id)):
                edge_data = run_state.G.get_edge_data(r_node_id, neighbor)
                if edge_data:
                    edge_info = (
                        edge_data.get(0, edge_data)
                        if isinstance(edge_data, dict)
                        else edge_data[0]
                    )
                    # Build from_config / to_config for forward-only single-step path
                    # IMPORTANT: G.nodes stores Mol objects — must convert to SMILES for JSON
                    # For reactant node, use original input (preserves duplicates)
                    from_node_data = run_state.G.nodes.get(r_node_id, {})
                    to_node_data = run_state.G.nodes.get(neighbor, {})
                    if original_reactant_smiles:
                        from_mols = list(original_reactant_smiles)
                    else:
                        from_mols = sorted(
                            Chem.MolToSmiles(mol) for mol in from_node_data.get("mols", set())
                            if mol is not None
                        )
                    to_mols = sorted(
                        Chem.MolToSmiles(mol) for mol in to_node_data.get("mols", set())
                        if mol is not None
                    )
                    # Arrow info for forward-only step (cached by reaction_smarts)
                    fo_smarts = edge_info.get("reaction_smarts", "")
                    fo_arrow_info = _arrow_info_cache.get(fo_smarts)
                    if fo_arrow_info is None:
                        fo_arrow_info = _get_step_arrow_info(fo_smarts)
                        _arrow_info_cache[fo_smarts] = fo_arrow_info

                    paths.append({
                        "path_id": fi,
                        "steps": [{
                            "step_num": 1,
                            "from_config": {
                                "molecules": from_mols,
                                "source": from_node_data.get("source", "unknown"),
                            },
                            "to_config": {
                                "molecules": to_mols,
                                "source": to_node_data.get("source", "unknown"),
                            },
                            "rule_id": edge_info.get("rule_id", ""),
                            "rule_name": edge_info.get("rule_name", ""),
                            "reaction_smarts": edge_info.get("reaction_smarts", ""),
                            "direction": "forward",
                            "score": round(edge_info.get("score", 0), 3),
                            "max_bond_length": round(edge_info.get("max_bond_length", 0), 3),
                            "arrows": fo_arrow_info.get("arrows", []),
                            "arrow_svg": fo_arrow_info.get("arrow_svg", ""),
                        }],
                        "total_score": round(edge_info.get("score", 0), 3),
                        "total_steps": 1,
                        "num_steps": 1,
                    })

    return {
        "reactants": reactants_smiles,
        "products": products_smiles,
        "paths": paths,
        "graph": {"nodes": nodes_list, "edges": edges_list},
        "search_time": round(search_time, 3),
        "total_time": round(search_time, 3),
        "stats": {
            "search_mode": search_mode,
            "forward_explored": len(run_state.explored),
            "backward_explored": (
                len(run_state.explored) if search_mode == "bidirectional" else 0
            ),
            "forward_configs": run_state.G.number_of_nodes(),
            "backward_configs": (
                run_state.G.number_of_nodes() if search_mode == "bidirectional" else 0
            ),
            "meeting_points": meeting_points_found,
            "total_rules": len(rules),
            "paths_found": len(paths),
            "max_depth_reached": max(
                (len(p.get("steps", [])) for p in paths), default=0
            ),
            "total_nodes": run_state.G.number_of_nodes(),
            "total_edges": run_state.G.number_of_edges(),
            "total_expanded": run_state.total_expanded,
            "search_time": round(search_time, 3),
        },
    }


def _empty_result(
    reactants_smiles: str,
    products_smiles: str,
    rules: List[Dict[str, Any]],
) -> Dict[str, Any]:
    return {
        "reactants": reactants_smiles,
        "products": products_smiles,
        "paths": [],
        "graph": {"nodes": [], "edges": []},
        "stats": {
            "search_mode": "error",
            "forward_explored": 0, "backward_explored": 0,
            "forward_configs": 0, "backward_configs": 0,
            "meeting_points": 0, "total_rules": len(rules),
            "paths_found": 0, "max_depth_reached": 0,
            "error": "No valid search could be performed",
        },
    }




# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="M-CSA Bidirectional Mechanism Search Engine",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 mechanism_search.py --reactants "CCO" --products "CC(=O)O"
  python3 mechanism_search.py --reactants "CCO" --products ""  (forward-only mode)
  python3 mechanism_search.py --reactants "CCO" --products "CC(=O)O" --max-configs 300 --max-rules 5000
        """,
    )
    parser.add_argument('--reactants', '-r', required=True,
                        help='Reactant SMILES (dot-separated for multiple)')
    parser.add_argument('--products', '-p', required=False, default='',
                        help='Product SMILES (dot-separated for multiple). Leave empty for forward-only search.')
    parser.add_argument('--max-configs', '-c', type=int, default=DEFAULT_MAX_CONFIGS,
                        help=f'Max configurations per side (default: {DEFAULT_MAX_CONFIGS})')
    parser.add_argument('--max-rules', type=int, default=DEFAULT_MAX_RULES,
                        help=f'Max M-CSA rules to load, 0=all (default: {DEFAULT_MAX_RULES}). '
                             'Note: frontend now uses --max-bond-length instead of --max-rules.')
    parser.add_argument('--max-depth', '-d', type=int, default=10,
                        help='Maximum search depth per side (default: 10)')
    parser.add_argument('--rules-file', type=str, default=RULES_FILE,
                        help=f'Path to rules Excel file (default: {RULES_FILE})')
    parser.add_argument('--reactant-pdb', type=str, default='',
                        help='PDB file text or path for reactant 3D coordinates')
    parser.add_argument('--product-pdb', type=str, default='',
                        help='PDB file text or path for product 3D coordinates')
    parser.add_argument('--max-bond-length', type=float, default=MAX_BOND_LENGTH,
                        help=f'Max bond length threshold in Angstroms (default: {MAX_BOND_LENGTH})')
    parser.add_argument('--residues', type=str, default='',
                        help='Catalytic residues as JSON string, e.g. '
                             '"[{"res_name":"SER","res_num":195,"chain":"A","part":"side_chain"}]"')
    parser.add_argument('--pdb-id', type=str, default='',
                        help='PDB identifier for the protein structure (e.g. "1HXR")')
    parser.add_argument('--pdb-text-file', type=str, default='',
                        help='Path to file containing PDB format text (avoids E2BIG with large PDB)')
    parser.add_argument('--rxn-text-file', type=str, default='',
                        help='Path to file containing RXN format reaction (atom-mapped, from Ketcher)')
    parser.add_argument('--ligand-mappings-file', type=str, default='',
                        help='Path to JSON file with ligand mappings [{smiles, res_name, chain, res_num}, ...]')
    parser.add_argument('--progress-file', type=str, default='',
                        help='Path to write progress JSON for async polling. '
                             'When set, search writes progress after each iteration. '
                             'Without this flag, operates in sync mode (JSON to stdout).')
    parser.add_argument('--result-file', type=str, default='',
                        help='Path to write final result JSON (used with --progress-file). '
                             'In sync mode (no --progress-file), result goes to stdout.')
    parser.add_argument('--mcsa-id', type=int, default=0,
                        help='Filter rules by M-CSA enzyme ID for "Own Rules" mode '
                             '(0 = use all rules, default)')
    parser.add_argument('--verify-mapping', action='store_true', default=False,
                        help='Verify atom mapping quality and output correction info as JSON')

    args = parser.parse_args()

    # ---- Quick mode: verify atom mapping only ----
    if args.verify_mapping:
        if not args.reactants or not args.products:
            print(json.dumps({"error": "Both --reactants and --products are required for verification"}))
            sys.exit(1)
        r, p = args.reactants, args.products
        corrected_r, corrected_p = remap_reaction_by_mcs(r, p)
        r_mols = [Chem.MolFromSmiles(s.strip()) for s in r.split('.') if s.strip()]
        p_mols = [Chem.MolFromSmiles(s.strip()) for s in p.split('.') if s.strip()]
        r_mols = [m for m in r_mols if m is not None]
        p_mols = [m for m in p_mols if m is not None]
        orig_changes, _, _ = _count_bond_changes(r_mols, p_mols)

        cr_mols = [Chem.MolFromSmiles(s.strip()) for s in corrected_r.split('.') if s.strip()]
        cp_mols = [Chem.MolFromSmiles(s.strip()) for s in corrected_p.split('.') if s.strip()]
        cr_mols = [m for m in cr_mols if m is not None]
        cp_mols = [m for m in cp_mols if m is not None]
        new_changes, only_r, only_p = _count_bond_changes(cr_mols, cp_mols)

        corrected = corrected_r != r or corrected_p != p
        result = {
            "original_bond_changes": orig_changes,
            "corrected_bond_changes": new_changes,
            "corrected": corrected,
            "original_reactants": r,
            "original_products": p,
            "corrected_reactants": corrected_r,
            "corrected_products": corrected_p,
            "bonds_broken": [{"atoms": sorted(bond_key), "type": bond_type} for bond_key, bond_type in only_r],
            "bonds_formed": [{"atoms": sorted(bond_key), "type": bond_type} for bond_key, bond_type in only_p],
        }
        print(json.dumps(result, indent=2))
        sys.exit(0)

    # Read PDB text from file if --pdb-text-file is provided
    # (avoids E2BIG error when PDB text is passed as CLI argument)
    pdb_text_from_file = ''
    if args.pdb_text_file:
        try:
            with open(args.pdb_text_file, 'r') as f:
                pdb_text_from_file = f.read()
            logger.info(f"Read PDB text from file: {len(pdb_text_from_file)} chars")
        except Exception as e:
            logger.warning(f"Failed to read PDB text file: {e}")

    # Read ligand_mappings from JSON file if --ligand-mappings-file is provided
    # Each mapping: {smiles, res_name, chain, res_num} for precise PDB coord extraction
    ligand_mappings: List[Dict[str, Any]] = []
    if args.ligand_mappings_file:
        try:
            with open(args.ligand_mappings_file, 'r') as f:
                ligand_mappings = json.load(f)
            logger.info(f"Read {len(ligand_mappings)} ligand mappings from file")
        except Exception as e:
            logger.warning(f"Failed to read ligand mappings file: {e}")

    # Read RXN text from file if --rxn-text-file is provided
    # RXN file from Ketcher may contain atom mapping (if user clicked Auto Map)
    rxn_text_from_file = ''
    if args.rxn_text_file:
        try:
            with open(args.rxn_text_file, 'r') as f:
                rxn_text_from_file = f.read()
            logger.info(f"Read RXN text from file: {len(rxn_text_from_file)} chars")
        except Exception as e:
            logger.warning(f"Failed to read RXN text file: {e}")

    # If RXN text is provided, parse it to extract reactant and product SMILES
    # (with atom map numbers preserved) — overrides CLI --reactants/--products
    if rxn_text_from_file:
        try:
            from rdkit.Chem import rdChemReactions as rxn_parser

            # Try RDKit's built-in parser first
            rxn = rxn_parser.ReactionFromRxnBlock(rxn_text_from_file)
            if rxn is not None:
                r_mols = rxn.GetReactants()
                p_mols = rxn.GetProducts()
            else:
                raise ValueError("ReactionFromRxnBlock returned None")

        except Exception as rxn_err:
            # Fallback: manual MOL block parsing (handles Ketcher-specific formats)
            logger.info(f"RDKit RXN parser failed ({rxn_err}), trying manual MOL block extraction")
            try:
                # Parse header to get reactant/product counts
                lines = rxn_text_from_file.strip().split('\n')
                n_reactants = 0
                n_products = 0
                for i, line in enumerate(lines):
                    stripped = line.strip()
                    if stripped.startswith('$RXN'):
                        for j in range(i + 1, min(i + 6, len(lines))):
                            parts = lines[j].split()
                            if len(parts) >= 2:
                                try:
                                    n_reactants = int(parts[0])
                                    n_products = int(parts[1])
                                    break
                                except ValueError:
                                    continue
                        break

                # Split into MOL blocks
                mol_sections = rxn_text_from_file.split('$MOL')[1:]
                r_mols = []
                p_mols = []
                for idx, section in enumerate(mol_sections):
                    mol = Chem.MolFromMolBlock('$MOL' + section)
                    if mol is not None:
                        if idx < n_reactants:
                            r_mols.append(mol)
                        else:
                            p_mols.append(mol)

                if not r_mols:
                    raise ValueError(f"Manual parsing found 0 reactant mols from {len(mol_sections)} MOL blocks")
            except Exception as manual_err:
                logger.warning(f"Manual RXN parsing also failed: {manual_err} — falling back to SMILES input")
                r_mols = []
                p_mols = []

        # Convert parsed mols to SMILES with atom map numbers preserved
        if r_mols:
            r_smiles_list = []
            for mol in r_mols:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True)
                if smi:
                    r_smiles_list.append(smi)
            p_smiles_list = []
            for mol in p_mols:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True)
                if smi:
                    p_smiles_list.append(smi)
            if r_smiles_list:
                args = argparse.Namespace(**vars(args))
                args.reactants = '.'.join(r_smiles_list)
                args.products = '.'.join(p_smiles_list)
                logger.info(f"RXN parsed: {len(r_smiles_list)} reactants, {len(p_smiles_list)} products")
                for i, smi in enumerate(r_smiles_list[:3]):
                    logger.info(f"  Reactant[{i}]: {smi[:100]}")
                for i, smi in enumerate(p_smiles_list[:3]):
                    logger.info(f"  Product[{i}]: {smi[:100]}")

    # ---- Atom mapping verification & correction ----
    # Indigo's automap can produce suboptimal mappings that inflate the reaction
    # center (e.g., swapping a carboxyl =O with a water molecule). This step
    # verifies the mapping and corrects it using MCS if needed.
    if args.products.strip() and args.reactants.strip():
        corrected_r, corrected_p = remap_reaction_by_mcs(args.reactants, args.products)
        if corrected_r != args.reactants or corrected_p != args.products:
            logger.info("Atom mapping corrected by MCS re-mapping")
            logger.info("  Original reactants: %s", args.reactants[:120])
            logger.info("  Corrected reactants: %s", corrected_r[:120])
            logger.info("  Original products: %s", args.products[:120])
            logger.info("  Corrected products: %s", corrected_p[:120])
            args = argparse.Namespace(**vars(args))
            args.reactants = corrected_r
            args.products = corrected_p

    t_total_start = time.time()

    t_load_start = time.time()
    # Always use lazy parsing — eager preparse of all rules causes OOM with 15K+
    # rules and even 12K rules takes ~25s. The lazy parser (_lazy_parse_rule) in the
    # search loop handles SMARTS parsing on demand, which is much faster overall.
    rules = shared_load_rules(max_rules=args.max_rules, rules_file=args.rules_file, preparse=False)
    load_time = time.time() - t_load_start
    logger.info("Loaded %d rules in %.3fs from SQLite/xlsx", len(rules), load_time)

    # Build coordinates cache from PDB inputs (optional)
    # coords_cache now ONLY stores coordinates for molecules WITHOUT ligand_mappings.
    # Molecules with ligand_mappings are handled by _build_am_to_3d_coord's
    # "consume-one-by-one" mode, which correctly assigns independent PDB coordinates
    # to each mapping entry (including multiple waters at different sites).
    coords_cache: Optional[Dict[str, Chem.Mol]] = None
    if pdb_text_from_file or args.reactant_pdb or args.product_pdb:
        coords_cache = {}
        reactant_mols = [s.strip() for s in args.reactants.split('.') if s.strip()]
        product_mols = [s.strip() for s in args.products.split('.') if s.strip()]
        all_mols_smiles = reactant_mols + product_mols

        # Build set of SMILES that have ligand_mappings — these will be handled
        # by _build_am_to_3d_coord's per-mapping processing, not by coords_cache.
        mapped_smiles_keys: Set[str] = set()
        for m in ligand_mappings:
            smi = m.get('smiles', '').strip()
            if smi:
                mapped_smiles_keys.add(_coords_cache_key(smi))

        if pdb_text_from_file:
            logger.info("Building coords_cache from PDB text (%d chars, %d molecules, %d ligand mappings)",
                        len(pdb_text_from_file), len(all_mols_smiles), len(ligand_mappings))
            logger.info("  %d SMILES have ligand mappings (handled by _build_am_to_3d_coord)",
                        len(mapped_smiles_keys))

            for smi in all_mols_smiles:
                key = _coords_cache_key(smi)
                if key in coords_cache:
                    continue

                if key in mapped_smiles_keys:
                    # This molecule has ligand mappings — skip coords_cache.
                    # _build_am_to_3d_coord will handle it with per-mapping PDB extraction.
                    continue

                # No ligand mapping — zero tolerance, no fallback
                logger.warning("No ligand mapping for %s, coordinates set to None (zero tolerance)", key[:40])
                coords_cache[key] = None
        else:
            # No PDB text available — not in PDB mode, no coords_cache needed
            logger.info("No PDB text available, skipping coords_cache build (non-PDB mode)")
            coords_cache = None

        if coords_cache is not None:
            logger.info("coords_cache: %d/%d molecules have 3D coordinates",
                         sum(1 for v in coords_cache.values() if v is not None and v.GetNumConformers() > 0),
                         len(all_mols_smiles))

    # Parse residues JSON (optional)
    residues: Optional[List[Dict[str, Any]]] = None
    if args.residues:
        try:
            residues = json.loads(args.residues)
            if not isinstance(residues, list):
                residues = None
        except (json.JSONDecodeError, TypeError):
            residues = None

    result = None
    try:
        result = bidirectional_search(
            reactants_smiles=args.reactants,
            products_smiles=args.products,
            rules=rules,
            max_configs=args.max_configs,
            max_depth=args.max_depth,
            coords_cache=coords_cache,
            max_bond_length=args.max_bond_length,
            residues=residues,
            pdb_id=args.pdb_id,
            pdb_text=pdb_text_from_file,
            progress_file=args.progress_file,
            result_file=args.result_file,
            mcsa_id=args.mcsa_id,
            ligand_mappings=ligand_mappings if ligand_mappings else None,
        )
    except Exception as e:
        # In async mode, write ERROR state to progress file
        if args.progress_file:
            write_progress(args.progress_file, "ERROR",
                           error=str(e),
                           elapsed_seconds=time.time() - t_total_start)
            logger.error("Search failed: %s", e, exc_info=True)
            # Don't re-raise — let the process exit cleanly so status endpoint
            # can report the error. Exit with non-zero code for diagnostics.
            sys.exit(1)
        raise

    total_time = time.time() - t_total_start
    result["load_time"] = round(load_time, 3)

    # Preserve backward-compatible output format
    stats = result.get("stats", {})
    search_time = stats.pop("search_time", None) if "search_time" in stats else None
    if search_time is None:
        search_time = round(total_time - load_time, 3)
    result["search_time"] = search_time
    result["total_time"] = round(total_time, 3)

    # In async mode (--progress-file), result is already written to result_file
    # by bidirectional_search(). Only print to stdout in sync mode.
    if not args.progress_file:
        print(json.dumps(result))


if __name__ == '__main__':
    main()
