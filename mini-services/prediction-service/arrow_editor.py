#!/usr/bin/env python3
"""
arrow_editor.py — 电子推动箭头编辑器模块

本模块提供了交互式电子推动箭头（Curly Arrow）编辑 API，支持：
- 分子原子索引查询（get_molecule_with_atom_indices）
- 箭头编辑应用（apply_arrow_edit）
- 反应箭头推断（infer_arrows_from_reaction）

依赖：rdkit, shared.py, arrow_transform.py（仅限 RDKit，无 Django 框架依赖）
"""

import argparse
import json
import math
import re
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

# 禁用 RDKit 警告信息
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdFMCS
from rdkit.Chem.rdchem import BondType

from shared import is_valid_smiles

# 导入 arrow_transform 模块的核心函数
from arrow_transform import (
    apply_arrows_to_smiles,
    validate_arrow_sequence,
    generate_reaction_scheme_svg,
    simulate_reaction_mechanism,
    parse_reaction_to_arrows,
)


# ---------------------------------------------------------------------------
# 分子原子信息查询
# ---------------------------------------------------------------------------

def get_molecule_with_atom_indices(smiles: str) -> Dict[str, Any]:
    """获取分子的原子索引、键信息及 2D 结构 SVG。

    将 SMILES 解析为 RDKit Mol 对象，遍历所有原子和键，
    提取原子索引、元素符号、原子序数、芳香性、度数、形式电荷等信息，
    以及键的连接关系和键级。

    Parameters
    ----------
    smiles : str
        输入分子的 SMILES 字符串。

    Returns
    -------
    Dict[str, Any]
        分子信息字典，包含：
        - smiles: 原始 SMILES
        - canonical_smiles: 规范 SMILES
        - num_atoms: 重原子数
        - atoms: 原子信息列表（index, symbol, atomic_num, is_aromatic, degree, formal_charge）
        - bonds: 键信息列表（source, target, order, type）
        - svg: 2D 分子结构 SVG 字符串

    示例
    -----
    >>> info = get_molecule_with_atom_indices("CCO")
    >>> info["num_atoms"]
    9
    """
    # 验证 SMILES
    if not is_valid_smiles(smiles):
        return {
            "smiles": smiles,
            "error": "Invalid SMILES characters",
        }

    # 解析分子
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "smiles": smiles,
            "error": "RDKit 无法解析此 SMILES",
        }

    # 规范 SMILES
    canonical = Chem.MolToSmiles(mol, canonical=True)

    # 提取原子信息
    atoms: List[Dict[str, Any]] = []
    for atom in mol.GetAtoms():
        atoms.append({
            "index": atom.GetIdx(),
            "symbol": atom.GetSymbol(),
            "atomic_num": atom.GetAtomicNum(),
            "is_aromatic": atom.GetIsAromatic(),
            "degree": atom.GetDegree(),
            "formal_charge": atom.GetFormalCharge(),
        })

    # 提取键信息
    bonds: List[Dict[str, Any]] = []
    for bond in mol.GetBonds():
        bond_type_str = _bond_type_to_string(bond.GetBondType())
        bonds.append({
            "source": bond.GetBeginAtomIdx(),
            "target": bond.GetEndAtomIdx(),
            "order": bond.GetBondTypeAsDouble(),
            "type": bond_type_str,
        })

    # 生成 2D 坐标和 SVG
    _ensure_2d_coords(mol)
    svg = _mol_to_svg_string(mol, width=350, height=300)

    return {
        "smiles": smiles,
        "canonical_smiles": canonical,
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "atoms": atoms,
        "bonds": bonds,
        "svg": svg,
    }


# ---------------------------------------------------------------------------
# 箭头编辑应用
# ---------------------------------------------------------------------------

def apply_arrow_edit(
    smiles: str,
    arrows: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """将箭头编辑应用到分子上并返回产物信息。

    接受一个 SMILES 和一个箭头编辑列表，对每个箭头编辑执行化学变换，
    生成产物并报告键变化详情。

    Parameters
    ----------
    smiles : str
        输入分子的 SMILES 字符串。
    arrows : List[Dict[str, Any]]
        箭头编辑列表，每个编辑为字典格式：
        ``{"source_atoms": [int, ...], "target_atoms": [int, ...], "electrons": int}``

    Returns
    -------
    Dict[str, Any]
        编辑结果，包含：
        - input_smiles: 输入 SMILES
        - product_smiles: 产物 SMILES（非规范）
        - product_canonical: 产物规范 SMILES
        - arrows_applied: 成功应用的箭头数量
        - success: 是否成功
        - validation: 验证结果（valid, warnings, charge_imbalance, electron_count）
        - bond_changes: 键变化列表（type, atoms, symbols）
        - product_svg: 产物 2D 结构 SVG
        - arrow_svg: 反应机理 SVG

    示例
    -----
    >>> result = apply_arrow_edit("CCO", [{"source_atoms": [1], "target_atoms": [2], "electrons": 2}])
    """
    # 验证输入
    if not is_valid_smiles(smiles):
        return {
            "input_smiles": smiles,
            "product_smiles": "",
            "product_canonical": "",
            "arrows_applied": 0,
            "success": False,
            "validation": {
                "valid": False,
                "warnings": ["Invalid SMILES characters"],
                "charge_imbalance": 0,
                "electron_count": 0,
            },
            "bond_changes": [],
            "product_svg": "",
            "arrow_svg": "",
            "error": "Invalid SMILES characters",
        }

    # 解析分子
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "input_smiles": smiles,
            "product_smiles": "",
            "product_canonical": "",
            "arrows_applied": 0,
            "success": False,
            "validation": {
                "valid": False,
                "warnings": ["RDKit cannot parse SMILES"],
                "charge_imbalance": 0,
                "electron_count": 0,
            },
            "bond_changes": [],
            "product_svg": "",
            "arrow_svg": "",
            "error": "RDKit cannot parse SMILES",
        }

    # 验证箭头序列
    validation = validate_arrow_sequence(mol, arrows)
    electron_count = sum(a.get("electrons", 2) for a in arrows)

    # 使用 simulate_reaction_mechanism 获取详细的键变化信息
    sim_result = simulate_reaction_mechanism(smiles, arrows)

    if not sim_result.get("success", False):
        return {
            "input_smiles": smiles,
            "product_smiles": sim_result.get("product_smiles", ""),
            "product_canonical": "",
            "arrows_applied": 0,
            "success": False,
            "validation": {
                "valid": validation.get("valid", False),
                "warnings": validation.get("warnings", []),
                "charge_imbalance": validation.get("charge_imbalance", 0),
                "electron_count": electron_count,
            },
            "bond_changes": [],
            "product_svg": "",
            "arrow_svg": "",
            "error": sim_result.get("error", "Arrow application failed"),
        }

    product_smiles = sim_result.get("product_smiles", "")

    # 计算规范 SMILES
    product_canonical = ""
    try:
        prod_mol = Chem.MolFromSmiles(product_smiles)
        if prod_mol is not None:
            product_canonical = Chem.MolToSmiles(prod_mol, canonical=True)
    except Exception:
        pass

    # 提取键变化
    bond_changes: List[Dict[str, Any]] = []
    all_step_changes = sim_result.get("bond_changes_per_step", [])
    for step_changes in all_step_changes:
        for change in step_changes:
            # 简化键变化输出格式
            simplified_change = {
                "type": change.get("type", ""),
                "atoms": change.get("atoms", []),
                "symbols": change.get("symbols", []),
            }
            bond_changes.append(simplified_change)

    # 生成产物 SVG
    product_svg = ""
    try:
        prod_mol = Chem.MolFromSmiles(product_smiles)
        if prod_mol is not None:
            _ensure_2d_coords(prod_mol)
            product_svg = _mol_to_svg_string(prod_mol, width=350, height=300)
    except Exception:
        pass

    # 生成反应机理 SVG
    arrow_svg = ""
    try:
        arrow_svg = generate_reaction_scheme_svg(smiles, arrows, product_smiles)
    except Exception:
        pass

    return {
        "input_smiles": smiles,
        "product_smiles": product_smiles,
        "product_canonical": product_canonical,
        "arrows_applied": len(arrows),
        "success": True,
        "validation": {
            "valid": validation.get("valid", False),
            "warnings": validation.get("warnings", []),
            "charge_imbalance": validation.get("charge_imbalance", 0),
            "electron_count": electron_count,
        },
        "bond_changes": bond_changes,
        "product_svg": product_svg,
        "arrow_svg": arrow_svg,
    }


# ---------------------------------------------------------------------------
# 反应箭头推断
# ---------------------------------------------------------------------------

def infer_arrows_from_reaction(
    reactant_smiles: str,
    product_smiles: str,
) -> Dict[str, Any]:
    """从反应物到产物推断电子推动箭头序列。

    使用 RDKit 的 MCS（最大公共子结构）算法查找反应物和产物之间的
    原子对应关系，然后比较键集合差异推断出键断裂和键形成事件，
    最终生成电子推动箭头序列。

    Parameters
    ----------
    reactant_smiles : str
        反应物 SMILES 字符串（可以是多组分，用 . 分隔）。
    product_smiles : str
        产物 SMILES 字符串（可以是多组分，用 . 分隔）。

    Returns
    -------
    Dict[str, Any]
        推断结果，包含：
        - reactant_smiles: 反应物 SMILES
        - product_smiles: 产物 SMILES
        - arrows: 推断的箭头列表
        - bond_changes: 键变化列表
        - reaction_svg: 反应机理 SVG
        - confidence: 推断置信度（0-1）

    示例
    -----
    >>> result = infer_arrows_from_reaction("CCO", "CC=O")
    """
    # 验证输入
    if not is_valid_smiles(reactant_smiles) or not is_valid_smiles(product_smiles):
        return {
            "reactant_smiles": reactant_smiles,
            "product_smiles": product_smiles,
            "arrows": [],
            "bond_changes": [],
            "reaction_svg": "",
            "confidence": 0.0,
            "error": "Invalid SMILES characters",
        }

    # 解析反应物和产物分子
    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
    product_mol = Chem.MolFromSmiles(product_smiles)

    if reactant_mol is None or product_mol is None:
        return {
            "reactant_smiles": reactant_smiles,
            "product_smiles": product_smiles,
            "arrows": [],
            "bond_changes": [],
            "reaction_svg": "",
            "confidence": 0.0,
            "error": "RDKit cannot parse one or both SMILES",
        }

    # 使用 MCS 查找原子映射
    mapping = _find_atom_mapping_mcs(reactant_mol, product_mol)

    if not mapping:
        # 尝试不要求完整匹配的 MCS
        mapping = _find_atom_mapping_mcs_relaxed(reactant_mol, product_mol)

    if not mapping:
        # 如果仍然没有映射，使用原子索引直接对应（可能不准确）
        mapping = _fallback_atom_mapping(reactant_mol, product_mol)

    # 构建产物侧的键集合（使用映射后的反应物索引）
    product_bonds = _get_mapped_bond_set(product_mol, mapping)

    # 构建反应物的键集合
    reactant_bonds = _get_bond_set(reactant_mol)

    # 比较键集合差异
    reactant_pairs = {(a, b) for a, b, _ in reactant_bonds}
    product_pairs = {(a, b) for a, b, _ in product_bonds}
    reactant_order_map = {(a, b): o for a, b, o in reactant_bonds}
    product_order_map = {(a, b): o for a, b, o in product_bonds}

    arrows: List[Dict[str, Any]] = []
    bond_changes: List[Dict[str, Any]] = []

    # 断裂的键
    broken = reactant_pairs - product_pairs
    for a, b in broken:
        sym_a = reactant_mol.GetAtomWithIdx(a).GetSymbol() if a < reactant_mol.GetNumAtoms() else '?'
        sym_b = reactant_mol.GetAtomWithIdx(b).GetSymbol() if b < reactant_mol.GetNumAtoms() else '?'
        arrows.append({
            "source_atoms": [a, b],
            "target_atoms": [a],
            "electrons": 2,
            "type": "bond_break",
        })
        bond_changes.append({
            "type": "broken",
            "atoms": [a, b],
            "symbols": [sym_a, sym_b],
            "old_order": reactant_order_map.get((a, b), 1),
            "new_order": 0,
        })

    # 新键
    formed = product_pairs - reactant_pairs
    for a, b in formed:
        sym_a = reactant_mol.GetAtomWithIdx(a).GetSymbol() if a < reactant_mol.GetNumAtoms() else '?'
        sym_b = reactant_mol.GetAtomWithIdx(b).GetSymbol() if b < reactant_mol.GetNumAtoms() else '?'
        arrows.append({
            "source_atoms": [a],
            "target_atoms": [b],
            "electrons": 2,
            "type": "bond_form",
        })
        bond_changes.append({
            "type": "formed",
            "atoms": [a, b],
            "symbols": [sym_a, sym_b],
            "old_order": 0,
            "new_order": product_order_map.get((a, b), 1),
        })

    # 键级变化
    common = reactant_pairs & product_pairs
    for a, b in common:
        old_order = reactant_order_map.get((a, b), 1)
        new_order = product_order_map.get((a, b), 1)
        if abs(old_order - new_order) > 0.01:
            sym_a = reactant_mol.GetAtomWithIdx(a).GetSymbol() if a < reactant_mol.GetNumAtoms() else '?'
            sym_b = reactant_mol.GetAtomWithIdx(b).GetSymbol() if b < reactant_mol.GetNumAtoms() else '?'
            arrows.append({
                "source_atoms": [a],
                "target_atoms": [b],
                "electrons": 2,
                "type": "order_change",
            })
            bond_changes.append({
                "type": "order_changed",
                "atoms": [a, b],
                "symbols": [sym_a, sym_b],
                "old_order": old_order,
                "new_order": new_order,
            })

    # 计算置信度
    confidence = _compute_confidence(
        reactant_mol, product_mol, mapping, arrows
    )

    # 生成反应 SVG
    reaction_svg = ""
    try:
        if arrows:
            reaction_svg = generate_reaction_scheme_svg(
                reactant_smiles, arrows, product_smiles
            )
    except Exception:
        pass

    return {
        "reactant_smiles": reactant_smiles,
        "product_smiles": product_smiles,
        "arrows": arrows,
        "bond_changes": bond_changes,
        "reaction_svg": reaction_svg,
        "confidence": round(confidence, 2),
    }


# ---------------------------------------------------------------------------
# 辅助函数
# ---------------------------------------------------------------------------

def _bond_type_to_string(bond_type: BondType) -> str:
    """将 RDKit BondType 转换为字符串表示。

    Parameters
    ----------
    bond_type : BondType
        RDKit 键类型。

    Returns
    -------
    str
        键类型字符串（SINGLE, DOUBLE, TRIPLE, AROMATIC 等）。
    """
    type_map = {
        BondType.SINGLE: "SINGLE",
        BondType.DOUBLE: "DOUBLE",
        BondType.TRIPLE: "TRIPLE",
        BondType.AROMATIC: "AROMATIC",
        BondType.UNSPECIFIED: "UNSPECIFIED",
    }
    return type_map.get(bond_type, str(bond_type))


def _ensure_2d_coords(mol: Chem.Mol) -> bool:
    """确保分子具有 2D 坐标。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子。

    Returns
    -------
    bool
        是否成功生成 2D 坐标。
    """
    try:
        if mol.GetNumConformers() == 0:
            AllChem.Compute2DCoords(mol)
        return mol.GetNumConformers() > 0
    except Exception:
        return False


def _mol_to_svg_string(mol: Chem.Mol, width: int = 300, height: int = 250) -> str:
    """将 RDKit Mol 对象转换为 SVG 字符串。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子。
    width : int
        SVG 宽度。
    height : int
        SVG 高度。

    Returns
    -------
    str
        SVG 字符串。
    """
    try:
        drawer = Draw.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception:
        return ''


def _get_bond_set(mol: Chem.Mol) -> Set[Tuple[int, int, float]]:
    """获取分子的键集合。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子。

    Returns
    -------
    Set[Tuple[int, int, float]]
        键集合，每项为 (原子1索引, 原子2索引, 键级)。
    """
    bonds: Set[Tuple[int, int, float]] = set()
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()
        bonds.add((min(a1, a2), max(a1, a2), order))
    return bonds


def _find_atom_mapping_mcs(
    reactant_mol: Chem.Mol,
    product_mol: Chem.Mol,
) -> Dict[int, int]:
    """使用 MCS 查找反应物到产物的原子映射。

    Parameters
    ----------
    reactant_mol : Chem.Mol
        反应物分子。
    product_mol : Chem.Mol
        产物分子。

    Returns
    -------
    Dict[int, int]
        反应物原子索引 → 产物原子索引的映射字典。
        空字典表示未找到有效映射。
    """
    try:
        mcs_result = rdFMCS.FindMCS(
            [reactant_mol, product_mol],
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            matchValences=False,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            timeout=10,
        )

        if mcs_result.numAtoms == 0 or mcs_result.smartsString is None:
            return {}

        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is None:
            return {}

        # 获取反应物中的 MCS 匹配
        reactant_match = reactant_mol.GetSubstructMatch(mcs_mol)
        product_match = product_mol.GetSubstructMatch(mcs_mol)

        if not reactant_match or not product_match:
            return {}

        # 构建映射
        mapping: Dict[int, int] = {}
        for r_idx, p_idx in zip(reactant_match, product_match):
            mapping[r_idx] = p_idx

        return mapping

    except Exception:
        return {}


def _find_atom_mapping_mcs_relaxed(
    reactant_mol: Chem.Mol,
    product_mol: Chem.Mol,
) -> Dict[int, int]:
    """使用宽松的 MCS 匹配策略查找原子映射。

    与 _find_atom_mapping_mcs 不同，这里尝试逐组分匹配
    （对于多组分分子更有效）。

    Parameters
    ----------
    reactant_mol : Chem.Mol
        反应物分子。
    product_mol : Chem.Mol
        产物分子。

    Returns
    -------
    Dict[int, int]
        反应物原子索引 → 产物原子索引的映射字典。
    """
    try:
        # 尝试更宽松的匹配参数
        mcs_result = rdFMCS.FindMCS(
            [reactant_mol, product_mol],
            atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            matchValences=False,
            ringMatchesRingOnly=False,
            completeRingsOnly=False,
            maximizeBonds=False,
            timeout=10,
            threshold=0.5,
        )

        if mcs_result.numAtoms == 0 or mcs_result.smartsString is None:
            return {}

        mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
        if mcs_mol is None:
            return {}

        reactant_match = reactant_mol.GetSubstructMatch(mcs_mol)
        product_match = product_mol.GetSubstructMatch(mcs_mol)

        if not reactant_match or not product_match:
            return {}

        mapping: Dict[int, int] = {}
        for r_idx, p_idx in zip(reactant_match, product_match):
            mapping[r_idx] = p_idx

        return mapping

    except Exception:
        return {}


def _fallback_atom_mapping(
    reactant_mol: Chem.Mol,
    product_mol: Chem.Mol,
) -> Dict[int, int]:
    """基于原子元素和邻域相似度的回退原子映射方法。

    当 MCS 无法找到映射时，使用贪心策略根据原子的元素类型和
    邻接关系建立映射。

    Parameters
    ----------
    reactant_mol : Chem.Mol
        反应物分子。
    product_mol : Chem.Mol
        产物分子。

    Returns
    -------
    Dict[int, int]
        反应物原子索引 → 产物原子索引的映射字典。
    """
    mapping: Dict[int, int] = {}

    # 获取两个分子的原子签名（元素+度数）
    def atom_signature(mol: Chem.Mol, idx: int) -> Tuple[int, int, str]:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = sorted(atom.GetNeighbors(), key=lambda a: a.GetIdx())
        neighbor_elems = tuple(sorted(n.GetSymbol() for n in neighbors))
        return (atom.GetAtomicNum(), atom.GetDegree(), neighbor_elems)

    # 构建反应物和产物的签名字典
    reactant_sigs: Dict[Tuple, List[int]] = {}
    for atom in reactant_mol.GetAtoms():
        sig = atom_signature(reactant_mol, atom.GetIdx())
        reactant_sigs.setdefault(sig, []).append(atom.GetIdx())

    product_sigs: Dict[Tuple, List[int]] = {}
    for atom in product_mol.GetAtoms():
        sig = atom_signature(product_mol, atom.GetIdx())
        product_sigs.setdefault(sig, []).append(atom.GetIdx())

    # 按签名匹配（贪心）
    used_r: set = set()
    used_p: set = set()
    for sig in reactant_sigs:
        if sig in product_sigs:
            r_list = [i for i in reactant_sigs[sig] if i not in used_r]
            p_list = [i for i in product_sigs[sig] if i not in used_p]
            for r, p in zip(r_list, p_list):
                mapping[r] = p
                used_r.add(r)
                used_p.add(p)

    return mapping


def _get_mapped_bond_set(
    product_mol: Chem.Mol,
    mapping: Dict[int, int],
) -> Set[Tuple[int, int, float]]:
    """获取产物分子的键集合，使用反应物索引表示。

    将产物中的原子索引通过映射转换为反应物侧的索引。

    Parameters
    ----------
    product_mol : Chem.Mol
        产物分子。
    mapping : Dict[int, int]
        反应物原子索引 → 产物原子索引的映射。
        需要反转以便从产物索引查反应物索引。

    Returns
    -------
    Set[Tuple[int, int, float]]
        映射后的键集合（使用反应物索引）。
    """
    # 反转映射：产物索引 → 反应物索引
    reverse_mapping: Dict[int, int] = {}
    for r_idx, p_idx in mapping.items():
        reverse_mapping[p_idx] = r_idx

    bonds: Set[Tuple[int, int, float]] = set()
    for bond in product_mol.GetBonds():
        p_a1 = bond.GetBeginAtomIdx()
        p_a2 = bond.GetEndAtomIdx()
        order = bond.GetBondTypeAsDouble()

        # 尝试映射到反应物索引
        r_a1 = reverse_mapping.get(p_a1)
        r_a2 = reverse_mapping.get(p_a2)

        if r_a1 is not None and r_a2 is not None:
            bonds.add((min(r_a1, r_a2), max(r_a1, r_a2), order))

    return bonds


def _compute_confidence(
    reactant_mol: Chem.Mol,
    product_mol: Chem.Mol,
    mapping: Dict[int, int],
    arrows: List[Dict[str, Any]],
) -> float:
    """计算箭头推断的置信度。

    置信度基于以下因素：
    1. 原子映射覆盖率（映射原子数 / 反应物原子总数）
    2. 产物侧覆盖率（映射产物原子数 / 产物原子总数）
    3. 键变化的合理性（较少的异常键变化 → 更高置信度）

    Parameters
    ----------
    reactant_mol : Chem.Mol
        反应物分子。
    product_mol : Chem.Mol
        产物分子。
    mapping : Dict[int, int]
        原子映射。
    arrows : List[Dict[str, Any]]
        推断的箭头列表。

    Returns
    -------
    float
        置信度（0.0-1.0）。
    """
    if not mapping:
        return 0.0

    # 原子映射覆盖率
    r_coverage = len(mapping) / max(reactant_mol.GetNumAtoms(), 1)
    p_coverage = len(mapping) / max(product_mol.GetNumAtoms(), 1)
    coverage_score = min(r_coverage, p_coverage)

    # 键变化合理性
    total_arrows = len(arrows)
    if total_arrows == 0:
        # 没有键变化 → 可能是相同分子
        r_smi = Chem.MolToSmiles(reactant_mol, canonical=True)
        p_smi = Chem.MolToSmiles(product_mol, canonical=True)
        if r_smi == p_smi:
            return 1.0
        arrow_score = 0.3  # 可能有未检测到的变化
    elif total_arrows <= 3:
        arrow_score = 0.9
    elif total_arrows <= 6:
        arrow_score = 0.7
    else:
        arrow_score = 0.5

    # 综合置信度
    confidence = coverage_score * 0.6 + arrow_score * 0.4

    return min(max(confidence, 0.0), 1.0)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    """命令行入口：提供三个子命令用于交互式箭头编辑。

    子命令：
        atoms <smiles>                  — 查询分子原子信息
        apply <smiles> <arrows_json>    — 应用箭头编辑
        infer <reactant_smiles> <product_smiles> — 推断反应箭头

    示例：
        python3 arrow_editor.py atoms "CCO"
        python3 arrow_editor.py apply "CCO" '[{"source_atoms":[1],"target_atoms":[2],"electrons":2}]'
        python3 arrow_editor.py infer "CCO" "CC=O"
    """
    parser = argparse.ArgumentParser(
        description="电子推动箭头编辑器 — 查询分子原子信息、应用箭头编辑、推断反应箭头",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
子命令：
  atoms <smiles>                              查询分子原子索引和键信息
  apply <smiles> <arrows_json>                应用箭头编辑到分子
  infer <reactant_smiles> <product_smiles>    从反应推断电子推动箭头

示例：
  python3 arrow_editor.py atoms "CCO"
  python3 arrow_editor.py apply "CCO" '[{"source_atoms":[1],"target_atoms":[2],"electrons":2}]'
  python3 arrow_editor.py infer "CCO" "CC=O"
        """,
    )
    subparsers = parser.add_subparsers(dest='command', help='子命令')

    # atoms 子命令
    atoms_parser = subparsers.add_parser('atoms', help='查询分子原子索引和键信息')
    atoms_parser.add_argument('smiles', help='SMILES 字符串')

    # apply 子命令
    apply_parser = subparsers.add_parser('apply', help='应用箭头编辑到分子')
    apply_parser.add_argument('smiles', help='输入分子 SMILES')
    apply_parser.add_argument('arrows_json', help='箭头编辑 JSON 字符串')

    # infer 子命令
    infer_parser = subparsers.add_parser('infer', help='从反应推断电子推动箭头')
    infer_parser.add_argument('reactant', help='反应物 SMILES')
    infer_parser.add_argument('product', help='产物 SMILES')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    if args.command == 'atoms':
        result = get_molecule_with_atom_indices(args.smiles)
        print(json.dumps(result, ensure_ascii=False))

    elif args.command == 'apply':
        # 解析箭头 JSON
        try:
            arrows = json.loads(args.arrows_json)
            if not isinstance(arrows, list):
                print(json.dumps({'error': 'arrows_json 必须是一个数组'}))
                sys.exit(1)
        except json.JSONDecodeError as e:
            print(json.dumps({'error': f'JSON 解析失败: {str(e)}'}))
            sys.exit(1)

        result = apply_arrow_edit(args.smiles, arrows)
        print(json.dumps(result, ensure_ascii=False))

    elif args.command == 'infer':
        result = infer_arrows_from_reaction(args.reactant, args.product)
        print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
