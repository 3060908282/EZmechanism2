#!/usr/bin/env python3
"""MCS-based Molecule Merging Module."""

import logging
from typing import Dict, List, Optional, Set, Tuple
from rdkit import Chem
from rdkit.Chem import rdFMCS

logger = logging.getLogger(__name__)


def _mol_is_proton(mol):
    return mol.GetNumAtoms() == 1 and mol.GetAtoms()[0].GetSymbol() == "H"


def _mol_is_water_species(mol):
    return (mol.GetNumHeavyAtoms() == 1
            and any(a.GetSymbol() == "O" for a in mol.GetAtoms()))


def _mol_is_hydronium(mol):
    return (mol.GetNumHeavyAtoms() == 1
            and any(a.GetSymbol() == "O" and a.GetFormalCharge() == 1
                    for a in mol.GetAtoms()))


def _first_proton_idx(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            return atom.GetIdx()
    return -1


def merge_molecules_by_mcs(mol_a, mol_b, merge_atoms_a=None, merge_atoms_b=None):
    try:
        rw_a = Chem.RWMol(mol_a)
        rw_b = Chem.RWMol(mol_b)
        maa = merge_atoms_a or set()
        mab = merge_atoms_b or set()
        for atom in rw_a.GetAtoms():
            iso = atom.GetAtomicNum()
            if atom.GetIdx() in maa:
                iso += 100
            ch = atom.GetFormalCharge()
            if ch:
                iso += ((10 + ch) % 10) * 1000
            atom.SetIsotope(iso)
        for atom in rw_b.GetAtoms():
            iso = atom.GetAtomicNum()
            if atom.GetIdx() in mab:
                iso += 100
            ch = atom.GetFormalCharge()
            if ch:
                iso += ((10 + ch) % 10) * 1000
            atom.SetIsotope(iso)
        mcs = Chem.rdFMCS.FindMCS(
            [rw_a, rw_b],
            atomCompare=Chem.rdFMCS.AtomCompare.CompareIsotopes,
            timeout=2)
        if mcs.numAtoms == 0:
            return None
        mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
        if mcs_mol is None:
            return None
        ma = rw_a.GetSubstructMatch(mcs_mol)
        mb = rw_b.GetSubstructMatch(mcs_mol)
        if not ma or not mb:
            return None
        b2a = {ib: ia for ia, ib in zip(ma, mb)}
        matched_b = set(mb)
        offset = mol_a.GetNumAtoms()
        comb = Chem.RWMol(Chem.CombineMols(mol_a, mol_b))
        for bond in rw_b.GetBonds():
            bg, ed = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bm, em = bg in b2a, ed in b2a
            if bm and em:
                continue
            elif bm and not em:
                ai, bi = b2a[bg], ed + offset
                if not comb.GetBondBetweenAtoms(ai, bi):
                    comb.AddBond(ai, bi, bond.GetBondType())
            elif not bm and em:
                ai, bi = b2a[ed], bg + offset
                if not comb.GetBondBetweenAtoms(ai, bi):
                    comb.AddBond(ai, bi, bond.GetBondType())
            else:
                nb, ne = bg + offset, ed + offset
                if not comb.GetBondBetweenAtoms(nb, ne):
                    comb.AddBond(nb, ne, bond.GetBondType())
        for idx in sorted([b + offset for b in matched_b], reverse=True):
            comb.RemoveAtom(idx)
        for atom in comb.GetAtoms():
            atom.SetIsotope(0)
        comb.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(comb)
        except Exception:
            pass
        return comb.GetMol()
    except Exception:
        return None


def find_matching_molecules(mol_list_a, mol_list_b):
    if not mol_list_a or not mol_list_b:
        return {}
    comp = {}
    for mol in mol_list_a + mol_list_b:
        rw = Chem.RWMol(mol)
        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass
        comp[id(mol)] = rw
    sims = []
    for ma in mol_list_a:
        for mb in mol_list_b:
            if ma.GetNumAtoms() == 1 and mb.GetNumAtoms() == 1:
                sub = ma
            elif (_mol_is_proton(ma) and _mol_is_hydronium(mb)) or                     (_mol_is_proton(mb) and _mol_is_hydronium(ma)):
                sub = ma if _mol_is_proton(ma) else mb
            else:
                try:
                    mcs = Chem.rdFMCS.FindMCS(
                        [comp[id(ma)], comp[id(mb)]], timeout=2)
                    sub = Chem.MolFromSmarts(mcs.smartsString)
                    if sub is None:
                        continue
                except Exception:
                    continue
            if sub and sub.GetNumAtoms() > 0:
                sims.append((ma, mb, sub))
    if not sims:
        return {}
    sims.sort(key=lambda x: -x[2].GetNumAtoms())
    pairs = {}
    ua, ub = set(), set()
    for ma, mb, sub in sims:
        if id(ma) in ua or id(mb) in ub:
            continue
        if _mol_is_proton(ma) and _mol_is_proton(mb):
            im = {ma.GetAtoms()[0].GetIdx(): mb.GetAtoms()[0].GetIdx()}
            am = dict(im)
        elif (_mol_is_proton(ma) and _mol_is_hydronium(mb)) or                 (_mol_is_proton(mb) and _mol_is_hydronium(ma)):
            pia = _first_proton_idx(ma)
            pib = _first_proton_idx(mb)
            if pia < 0 or pib < 0:
                continue
            im = {pia: pib}
            am = {ma.GetAtomWithIdx(pia).GetAtomMapNum():
                  mb.GetAtomWithIdx(pib).GetAtomMapNum()}
        else:
            oa = comp[id(ma)].GetSubstructMatch(sub)
            ob = comp[id(mb)].GetSubstructMatch(sub)
            if not oa or not ob or len(oa) != len(ob):
                continue
            im = {oa[i]: ob[i] for i in range(len(oa))}
            am = {}
            for ia, ib in im.items():
                ama = ma.GetAtomWithIdx(ia).GetAtomMapNum()
                amb = mb.GetAtomWithIdx(ib).GetAtomMapNum()
                if ama != 0 and amb != 0:
                    am[ama] = amb
        pairs[ma] = (mb, im, am)
        ua.add(id(ma))
        ub.add(id(mb))
    return pairs


def compare_molecules(mol_a, mol_b):
    if mol_a is None or mol_b is None:
        return 0.0
    mx = max(mol_a.GetNumAtoms(), mol_b.GetNumAtoms())
    if mx == 0:
        return 1.0
    try:
        mcs = Chem.rdFMCS.FindMCS(
            [mol_a, mol_b],
            atomCompare=Chem.rdFMCS.AtomCompare.CompareElements,
            bondCompare=Chem.rdFMCS.BondCompare.CompareOrder,
            timeout=2)
        return mcs.numAtoms / mx if mcs.numAtoms else 0.0
    except Exception:
        return 0.0


def merge_multi_component_configs(cma, cmb):
    if not cma or not cmb:
        r = list(cma) if cma else list(cmb)
        return [Chem.Mol(m) for m in r]
    mp = find_matching_molecules(cma, cmb)
    res = []
    mb_ids = set()
    for ma, (mb, im, am) in mp.items():
        merged = merge_molecules_by_mcs(ma, mb)
        res.append(merged if merged else Chem.Mol(ma))
        mb_ids.add(id(mb))
    for ma in cma:
        if ma not in mp:
            res.append(Chem.Mol(ma))
    for mb in cmb:
        if id(mb) not in mb_ids:
            res.append(Chem.Mol(mb))
    return res


def merge_configs_from_smiles(sa, sb):
    ma_list = [Chem.MolFromSmiles(s) for s in sa]
    mb_list = [Chem.MolFromSmiles(s) for s in sb]
    ma_list = [m for m in ma_list if m]
    mb_list = [m for m in mb_list if m]
    if not ma_list or not mb_list:
        return ([s for s in sa + sb if s], 0.0)
    ts, c = 0.0, 0
    for a in ma_list:
        for b in mb_list:
            ts += compare_molecules(a, b)
            c += 1
    avg = ts / c if c else 0.0
    merged = merge_multi_component_configs(ma_list, mb_list)
    smis = [Chem.MolToSmiles(m) for m in merged if Chem.MolToSmiles(m)]
    return (smis, round(avg, 4))
