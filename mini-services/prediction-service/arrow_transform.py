#!/usr/bin/env python3
"""
电子推动箭头（Curly Arrow）反应变换模块

本模块提供了基于 EzMechanism 论文的高层次电子推动箭头变换接口，
用于多步箭头序列模拟、反应机理分析以及反应 SVG 可视化。

主要功能：
- 多步箭头序列模拟（simulate_reaction_mechanism）
- 反应 SMARTS 到箭头序列的推断（parse_reaction_to_arrows）
- SMILES + 箭头快捷变换（apply_arrows_to_smiles）
- 箭头序列化学合理性验证（validate_arrow_sequence）
- 反应机理 SVG 图生成（generate_reaction_scheme_svg）

依赖：rdkit（仅限 RDKit，无 Django 框架依赖）
"""

from typing import Any, Dict, FrozenSet, List, Optional, Set, Tuple
import math
import io

# 禁用 RDKit 警告信息
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdchem import BondType

# 导入论文中的电子推动箭头变换函数
from common.utils import get_mol_product_from_arrow


# ===========================================================================
# 辅助函数
# ===========================================================================

def _apply_single_arrow_safe(
    mol: Chem.Mol,
    source_map_nums: FrozenSet[int],
    target_map_nums: FrozenSet[int],
    electrons: int,
) -> Chem.Mol:
    """安全地应用单个电子推动箭头，自动过滤自键（source == target）对。

    底层 ``get_mol_product_from_arrow`` 对 source × target 做笛卡尔积，
    当 source 和 target 有重叠原子时会产生自键错误。
    本函数先将重叠原子拆分为仅 source 和仅 target 子集，
    再分别调用底层函数，最后合并电荷效果。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子（需有原子映射编号）。
    source_map_nums : FrozenSet[int]
        源原子映射编号集合。
    target_map_nums : FrozenSet[int]
        目标原子映射编号集合。
    electrons : int
        电子数（1 或 2）。

    Returns
    -------
    Chem.Mol
        变换后的产物分子。
    """
    # 分离仅 source、仅 target 和重叠原子
    only_source = source_map_nums - target_map_nums
    only_target = target_map_nums - source_map_nums
    overlap = source_map_nums & target_map_nums

    # 构建映射查找表
    am_to_idx: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        amn = atom.GetAtomMapNum()
        if amn != 0:
            am_to_idx[amn] = atom.GetIdx()

    rwmol = Chem.RWMol(mol)

    # 处理仅 source × 仅 target 的有效键操作
    effective_source = only_source | overlap
    effective_target = only_target | overlap

    # 仅对不构成自键的原子对执行键操作
    processed_pairs: Set[Tuple[int, int]] = set()
    for s_am in effective_source:
        for t_am in effective_target:
            if s_am == t_am:
                continue  # 跳过自键
            pair = (min(s_am, t_am), max(s_am, t_am))
            if pair in processed_pairs:
                continue
            processed_pairs.add(pair)

            s_idx = am_to_idx.get(s_am)
            t_idx = am_to_idx.get(t_am)
            if s_idx is None or t_idx is None:
                continue

            bond = rwmol.GetBondBetweenAtoms(s_idx, t_idx)
            if bond is not None:
                current_order = bond.GetBondTypeAsDouble()
                new_order = current_order + electrons
                if new_order <= 0:
                    rwmol.RemoveBond(s_idx, t_idx)
                elif new_order >= 3:
                    bond.SetBondType(BondType.TRIPLE)
                elif new_order >= 2:
                    bond.SetBondType(BondType.DOUBLE)
                elif new_order >= 1.5:
                    bond.SetBondType(BondType.AROMATIC)
                else:
                    bond.SetBondType(BondType.SINGLE)
            else:
                rwmol.AddBond(s_idx, t_idx, BondType.SINGLE)

            # 电荷调整（仅对 2-电子箭头生效）
            if electrons == 2:
                s_atom = rwmol.GetAtomWithIdx(s_idx)
                t_atom = rwmol.GetAtomWithIdx(t_idx)
                s_atom.SetFormalCharge(s_atom.GetFormalCharge() + electrons)
                t_atom.SetFormalCharge(t_atom.GetFormalCharge() - electrons)

    # 对重叠原子（既是 source 又是 target）：调整电荷但不操作键
    # 重叠意味着电子从该原子"出发"又"回到"自身 → 净电荷不变
    # 但如果 source 有多个原子，重叠原子的电荷需要 +electrons
    # 如果 target 有多个原子，重叠原子的电荷需要 -electrons
    # 净效果：+electrons（作为 source 的一部分）- electrons（作为 target 的一部分）= 0

    try:
        Chem.SanitizeMol(rwmol)
    except Exception:
        pass

    return rwmol.GetMol()


def _smiles_to_mol(smiles: str, add_hs: bool = False) -> Optional[Chem.Mol]:
    """将 SMILES 字符串转换为 RDKit Mol 对象。

    Parameters
    ----------
    smiles : str
        输入的 SMILES 字符串。
    add_hs : bool
        是否添加显式氢原子（默认 False）。

    Returns
    -------
    Optional[Chem.Mol]
        转换后的 Mol 对象，解析失败时返回 None。
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        if add_hs:
            mol = Chem.AddHs(mol)
        return mol
    except Exception:
        return None


def _mol_to_smiles(mol: Chem.Mol, canonical: bool = True) -> str:
    """将 RDKit Mol 对象转换为 SMILES 字符串。

    Parameters
    ----------
    mol : Chem.Mol
        输入的 Mol 对象。
    canonical : bool
        是否使用规范 SMILES（默认 True）。

    Returns
    -------
    str
        SMILES 字符串，失败时返回空字符串。
    """
    try:
        smi = Chem.MolToSmiles(mol, canonical=canonical)
        return smi if smi else ''
    except Exception:
        return ''


def _assign_atom_map_nums(mol: Chem.Mol) -> Chem.Mol:
    """为分子中的每个重原子按索引顺序分配原子映射编号。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子。

    Returns
    -------
    Chem.Mol
        带有原子映射编号的分子副本。
    """
    rw = Chem.RWMol(mol)
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() > 0:
            atom.SetAtomMapNum(atom.GetIdx() + 1)  # 1-based map nums
    try:
        Chem.SanitizeMol(rw)
    except Exception:
        pass
    return rw.GetMol()


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


def _get_bond_set(mol: Chem.Mol) -> Set[Tuple[int, int, float]]:
    """获取分子的键集合。

    每个键表示为 (begin_atom_idx, end_atom_idx, bond_order) 的元组。

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


def _get_atom_count_dict(mol: Chem.Mol) -> Dict[int, int]:
    """获取分子中各元素的原子计数。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子。

    Returns
    -------
    Dict[int, int]
        原子序数到计数的映射。
    """
    counter: Dict[int, int] = {}
    for atom in mol.GetAtoms():
        counter[atom.GetAtomicNum()] = counter.get(atom.GetAtomicNum(), 0) + 1
    return counter


# ===========================================================================
# 核心函数
# ===========================================================================

def simulate_reaction_mechanism(
    reactant_smiles: str,
    arrow_sequence: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """模拟多步电子推动箭头反应机理。

    将反应物 SMILES 和一个箭头描述序列作为输入，依次应用每个箭头变换，
    生成各步中间体和最终产物。

    Parameters
    ----------
    reactant_smiles : str
        反应物的 SMILES 字符串。
    arrow_sequence : List[Dict[str, Any]]
        箭头描述序列，每个箭头为字典格式：
        ``{"source_atoms": [int, ...], "target_atoms": [int, ...], "electrons": int}``
        其中 ``source_atoms`` 和 ``target_atoms`` 中的数字为 0-based 原子索引，
        ``electrons`` 为电子数（通常为 1 或 2）。

    Returns
    -------
    Dict[str, Any]
        结果字典，包含以下字段：
        - ``product_smiles`` (str): 最终产物的 SMILES 字符串
        - ``intermediates`` (List[str]): 每步箭头应用后的中间体 SMILES 列表
        - ``bond_changes_per_step`` (List[List[Dict]]): 每步的键变化详情
        - ``success`` (bool): 是否成功完成所有变换
        - ``error`` (Optional[str]): 错误信息（无错误时为 None）

    示例
    -----
    >>> result = simulate_reaction_mechanism('CCO', [
    ...     {"source_atoms": [0, 1], "target_atoms": [1], "electrons": 2}
    ... ])
    """
    if not reactant_smiles:
        return {
            "product_smiles": "",
            "intermediates": [],
            "bond_changes_per_step": [],
            "success": False,
            "error": "空反应物 SMILES",
        }

    # 解析反应物分子
    mol = _smiles_to_mol(reactant_smiles)
    if mol is None:
        return {
            "product_smiles": "",
            "intermediates": [],
            "bond_changes_per_step": [],
            "success": False,
            "error": f"无法解析 SMILES: {reactant_smiles}",
        }

    # 分配原子映射编号（1-based，对应 0-based 索引 + 1）
    mol = _assign_atom_map_nums(mol)

    intermediates: List[str] = []
    bond_changes_per_step: List[List[Dict[str, Any]]] = []
    current_mol = mol

    try:
        for step_idx, arrow in enumerate(arrow_sequence):
            source_atoms = arrow.get("source_atoms", [])
            target_atoms = arrow.get("target_atoms", [])
            electrons = arrow.get("electrons", 2)

            if not source_atoms or not target_atoms:
                return {
                    "product_smiles": _mol_to_smiles(current_mol),
                    "intermediates": intermediates,
                    "bond_changes_per_step": bond_changes_per_step,
                    "success": False,
                    "error": f"第 {step_idx + 1} 步箭头缺少源或目标原子",
                }

            # 记录当前键集合
            prev_bonds = _get_bond_set(current_mol)

            # 将 0-based 索引转换为 1-based 原子映射编号
            source_map_nums = frozenset(idx + 1 for idx in source_atoms)
            target_map_nums = frozenset(idx + 1 for idx in target_atoms)

            # 使用安全的箭头变换函数（自动处理自键）
            current_mol = _apply_single_arrow_safe(
                current_mol, source_map_nums, target_map_nums, electrons
            )

            # 记录变换后的键集合
            post_bonds = _get_bond_set(current_mol)

            # 计算键变化
            step_changes = _compute_bond_changes(prev_bonds, post_bonds, current_mol)
            bond_changes_per_step.append(step_changes)

            # 记录中间体 SMILES
            intermediate_smi = _mol_to_smiles(current_mol)
            intermediates.append(intermediate_smi)

        return {
            "product_smiles": _mol_to_smiles(current_mol),
            "intermediates": intermediates,
            "bond_changes_per_step": bond_changes_per_step,
            "success": True,
            "error": None,
        }

    except Exception as e:
        return {
            "product_smiles": _mol_to_smiles(current_mol),
            "intermediates": intermediates,
            "bond_changes_per_step": bond_changes_per_step,
            "success": False,
            "error": f"第 {step_idx + 1} 步变换失败: {str(e)}",
        }


def _compute_bond_changes(
    prev_bonds: Set[Tuple[int, int, float]],
    post_bonds: Set[Tuple[int, int, float]],
    mol: Chem.Mol,
) -> List[Dict[str, Any]]:
    """计算两个键集合之间的变化详情。

    Parameters
    ----------
    prev_bonds : Set[Tuple[int, int, float]]
        变换前的键集合。
    post_bonds : Set[Tuple[int, int, float]]
        变换后的键集合。
    mol : Chem.Mol
        当前分子（用于获取原子符号）。

    Returns
    -------
    List[Dict[str, Any]]
        键变化列表，每项包含 type, atoms, old_order, new_order 等字段。
    """
    changes: List[Dict[str, Any]] = []

    # 构建仅含键对的集合（不含键级）以便比较
    prev_pairs = {(a, b) for a, b, _ in prev_bonds}
    post_pairs = {(a, b) for a, b, _ in post_bonds}

    # 查找键的键级映射
    prev_order_map = {(a, b): o for a, b, o in prev_bonds}
    post_order_map = {(a, b): o for a, b, o in post_bonds}

    # 新键（产物中有但反应物中没有的）
    formed = post_pairs - prev_pairs
    for a, b in formed:
        sym_a = mol.GetAtomWithIdx(a).GetSymbol() if a < mol.GetNumAtoms() else '?'
        sym_b = mol.GetAtomWithIdx(b).GetSymbol() if b < mol.GetNumAtoms() else '?'
        changes.append({
            "type": "formed",
            "atoms": [a, b],
            "symbols": [sym_a, sym_b],
            "old_order": 0,
            "new_order": post_order_map.get((a, b), 1),
        })

    # 断裂的键（反应物中有但产物中没有的）
    broken = prev_pairs - post_pairs
    for a, b in broken:
        sym_a = mol.GetAtomWithIdx(a).GetSymbol() if a < mol.GetNumAtoms() else '?'
        sym_b = mol.GetAtomWithIdx(b).GetSymbol() if b < mol.GetNumAtoms() else '?'
        changes.append({
            "type": "broken",
            "atoms": [a, b],
            "symbols": [sym_a, sym_b],
            "old_order": prev_order_map.get((a, b), 1),
            "new_order": 0,
        })

    # 键级变化（两侧都有但键级不同的）
    common = prev_pairs & post_pairs
    for a, b in common:
        old_o = prev_order_map.get((a, b), 1)
        new_o = post_order_map.get((a, b), 1)
        if abs(old_o - new_o) > 0.01:
            sym_a = mol.GetAtomWithIdx(a).GetSymbol() if a < mol.GetNumAtoms() else '?'
            sym_b = mol.GetAtomWithIdx(b).GetSymbol() if b < mol.GetNumAtoms() else '?'
            changes.append({
                "type": "order_changed",
                "atoms": [a, b],
                "symbols": [sym_a, sym_b],
                "old_order": old_o,
                "new_order": new_o,
            })

    return changes


def parse_reaction_to_arrows(reaction_smarts: str) -> List[Dict[str, Any]]:
    """从反应 SMARTS 推断电子推动箭头序列。

    通过比较反应物和产物的键集合变化，推断出可能的电子推动箭头。
    这是一种近似推断方法——论文中的实际箭头数据来自 M-CSA XML 文件。

    推断规则：
    - 键断裂 → 2-电子箭头：从断裂的键到保留电子的原子
    - 键形成 → 2-电子箭头：从亲核原子到亲电原子
    - 键级变化 → 2-电子箭头：在两个键合原子之间

    Parameters
    ----------
    reaction_smarts : str
        反应 SMARTS 字符串（格式：反应物 >> 产物）。

    Returns
    -------
    List[Dict[str, Any]]
        推断的箭头描述列表，每个箭头为：
        ``{"source_atoms": [int, ...], "target_atoms": [int, ...], "electrons": int}``

    示例
    -----
    >>> arrows = parse_reaction_to_arrows('[C:1][O:2]>>[C:1].[O:2]')
    """
    if not reaction_smarts or '>>' not in reaction_smarts:
        return []

    try:
        parts = reaction_smarts.split('>>', 1)
        reactant_smarts = parts[0].strip()
        product_smarts = parts[1].strip()

        if not reactant_smarts or not product_smarts:
            return []

        # 解析反应物和产物为带原子映射编号的分子
        reactant_mol = Chem.MolFromSmarts(reactant_smarts)
        product_mol = Chem.MolFromSmarts(product_smarts)

        if reactant_mol is None or product_mol is None:
            # 尝试作为普通 SMILES 解析
            reactant_mol = Chem.MolFromSmiles(reactant_smarts)
            product_mol = Chem.MolFromSmiles(product_smarts)
            if reactant_mol is None or product_mol is None:
                return []

        # 构建反应物的原子映射编号 → 索引查找表
        r_am_to_idx: Dict[int, int] = {}
        for atom in reactant_mol.GetAtoms():
            amn = atom.GetAtomMapNum()
            if amn != 0:
                r_am_to_idx[amn] = atom.GetIdx()

        # 构建产物的原子映射编号 → 索引查找表
        p_am_to_idx: Dict[int, int] = {}
        for atom in product_mol.GetAtoms():
            amn = atom.GetAtomMapNum()
            if amn != 0:
                p_am_to_idx[amn] = atom.GetIdx()

        # 检查是否有有效的原子映射编号
        if not r_am_to_idx or not p_am_to_idx:
            return []

        # 获取反应物和产物的键集合（使用原子映射编号）
        r_bonds: Set[Tuple[int, int, float]] = set()
        for bond in reactant_mol.GetBonds():
            am1 = bond.GetBeginAtom().GetAtomMapNum()
            am2 = bond.GetEndAtom().GetAtomMapNum()
            if am1 != 0 and am2 != 0:
                order = bond.GetBondTypeAsDouble()
                r_bonds.add((min(am1, am2), max(am1, am2), order))

        p_bonds: Set[Tuple[int, int, float]] = set()
        for bond in product_mol.GetBonds():
            am1 = bond.GetBeginAtom().GetAtomMapNum()
            am2 = bond.GetEndAtom().GetAtomMapNum()
            if am1 != 0 and am2 != 0:
                order = bond.GetBondTypeAsDouble()
                p_bonds.add((min(am1, am2), max(am1, am2), order))

        r_pairs = {(a, b) for a, b, _ in r_bonds}
        p_pairs = {(a, b) for a, b, _ in p_bonds}
        r_order_map = {(a, b): o for a, b, o in r_bonds}
        p_order_map = {(a, b): o for a, b, o in p_bonds}

        arrows: List[Dict[str, Any]] = []

        # 推断断裂的键 → 键断裂箭头
        broken_bonds = r_pairs - p_pairs
        for a, b in broken_bonds:
            # 键断裂：2-电子箭头，电子留在电负性更强的原子上
            # 使用原子索引（0-based）而非映射编号
            a_idx = r_am_to_idx.get(a, a - 1)
            b_idx = r_am_to_idx.get(b, b - 1)
            # 两个源原子（键上的两个原子），电子流向一个目标
            arrows.append({
                "source_atoms": [a_idx, b_idx],
                "target_atoms": [a_idx],  # 电子留在第一个原子上
                "electrons": 2,
                "description": f"键断裂 {a}-{b}",
            })

        # 推断新键形成 → 成键箭头
        formed_bonds = p_pairs - r_pairs
        for a, b in formed_bonds:
            a_idx = r_am_to_idx.get(a, a - 1)
            b_idx = r_am_to_idx.get(b, b - 1)
            # 新键形成：亲核原子向亲电原子提供电子
            arrows.append({
                "source_atoms": [a_idx],
                "target_atoms": [b_idx],
                "electrons": 2,
                "description": f"键形成 {a}-{b}",
            })

        # 推断键级变化
        common_bonds = r_pairs & p_pairs
        for a, b in common_bonds:
            old_order = r_order_map.get((a, b), 1)
            new_order = p_order_map.get((a, b), 1)
            if abs(new_order - old_order) > 0.01:
                a_idx = r_am_to_idx.get(a, a - 1)
                b_idx = r_am_to_idx.get(b, b - 1)
                arrows.append({
                    "source_atoms": [a_idx],
                    "target_atoms": [b_idx],
                    "electrons": 2,
                    "description": f"键级变化 {a}-{b}: {old_order}→{new_order}",
                })

        return arrows

    except Exception:
        return []


def apply_arrows_to_smiles(smiles: str, arrows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """将箭头序列应用到 SMILES 字符串上的便捷函数。

    将 SMILES 转换为 Mol 对象，应用所有箭头变换，返回产物 SMILES。

    Parameters
    ----------
    smiles : str
        输入分子的 SMILES 字符串。
    arrows : List[Dict[str, Any]]
        箭头描述列表，格式同 ``simulate_reaction_mechanism``。

    Returns
    -------
    Dict[str, Any]
        结果字典：
        - ``product_smiles`` (str): 产物 SMILES 字符串
        - ``success`` (bool): 是否成功
        - ``error`` (Optional[str]): 错误信息

    示例
    -----
    >>> result = apply_arrows_to_smiles('CCO', [
    ...     {"source_atoms": [0, 1], "target_atoms": [1], "electrons": 2}
    ... ])
    """
    if not smiles:
        return {
            "product_smiles": "",
            "success": False,
            "error": "空 SMILES 字符串",
        }

    if not arrows:
        return {
            "product_smiles": smiles,
            "success": True,
            "error": None,
        }

    mol = _smiles_to_mol(smiles)
    if mol is None:
        return {
            "product_smiles": "",
            "success": False,
            "error": f"无法解析 SMILES: {smiles}",
        }

    # 分配原子映射编号
    mol = _assign_atom_map_nums(mol)

    try:
        for arrow in arrows:
            source_atoms = arrow.get("source_atoms", [])
            target_atoms = arrow.get("target_atoms", [])
            electrons = arrow.get("electrons", 2)

            source_map_nums = frozenset(idx + 1 for idx in source_atoms)
            target_map_nums = frozenset(idx + 1 for idx in target_atoms)

            # 使用安全的箭头变换函数（自动处理自键）
            mol = _apply_single_arrow_safe(
                mol, source_map_nums, target_map_nums, electrons
            )

        product_smi = _mol_to_smiles(mol)
        return {
            "product_smiles": product_smi,
            "success": True,
            "error": None,
        }

    except Exception as e:
        return {
            "product_smiles": _mol_to_smiles(mol),
            "success": False,
            "error": str(e),
        }


def validate_arrow_sequence(
    mol: Chem.Mol,
    arrows: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """验证箭头序列的化学合理性。

    检查项目包括：
    - 电子守恒（每步的总电子数变化）
    - 不可能的键级（超过 3 或为负数）
    - 电荷平衡
    - 原子索引越界

    Parameters
    ----------
    mol : Chem.Mol
        输入的 RDKit Mol 对象。
    arrows : List[Dict[str, Any]]
        待验证的箭头序列。

    Returns
    -------
    Dict[str, Any]
        验证结果：
        - ``valid`` (bool): 序列是否化学上合理
        - ``warnings`` (List[str]): 警告信息列表
        - ``charge_imbalance`` (float): 总电荷不平衡量
    """
    warnings: List[str] = []
    total_charge_change = 0.0
    n_atoms = mol.GetNumAtoms()

    if not arrows:
        return {
            "valid": True,
            "warnings": [],
            "charge_imbalance": 0.0,
        }

    for step_idx, arrow in enumerate(arrows):
        source_atoms = arrow.get("source_atoms", [])
        target_atoms = arrow.get("target_atoms", [])
        electrons = arrow.get("electrons", 2)

        # 检查原子索引越界
        for idx in source_atoms:
            if idx < 0 or idx >= n_atoms:
                warnings.append(
                    f"第 {step_idx + 1} 步：源原子索引 {idx} 越界 (0~{n_atoms - 1})"
                )
        for idx in target_atoms:
            if idx < 0 or idx >= n_atoms:
                warnings.append(
                    f"第 {step_idx + 1} 步：目标原子索引 {idx} 越界 (0~{n_atoms - 1})"
                )

        # 检查电子数合理性
        if electrons not in (1, 2):
            warnings.append(
                f"第 {step_idx + 1} 步：电子数 {electrons} 异常（应为 1 或 2）"
            )

        # 检查源和目标原子重叠
        overlap = set(source_atoms) & set(target_atoms)
        if overlap:
            warnings.append(
                f"第 {step_idx + 1} 步：源和目标原子有重叠 {overlap}"
            )

        # 模拟键级变化
        for s_idx in source_atoms:
            for t_idx in target_atoms:
                if s_idx < 0 or s_idx >= n_atoms or t_idx < 0 or t_idx >= n_atoms:
                    continue
                bond = mol.GetBondBetweenAtoms(s_idx, t_idx)
                if bond is not None:
                    new_order = bond.GetBondTypeAsDouble() + electrons
                    if new_order > 3.5:
                        warnings.append(
                            f"第 {step_idx + 1} 步：键级 {new_order:.1f} 超过三键上限"
                        )
                    if new_order < 0:
                        warnings.append(
                            f"第 {step_idx + 1} 步：键级 {new_order:.1f} 为负值"
                        )

        # 电荷变化追踪（2-电子箭头：源原子 +2，目标原子 -2）
        if electrons == 2:
            total_charge_change += len(source_atoms) * 2
            total_charge_change -= len(target_atoms) * 2

    # 电荷平衡分析
    charge_imbalance = abs(total_charge_change)
    if charge_imbalance > 0 and charge_imbalance != len(arrows) * 2:
        warnings.append(
            f"电荷不平衡: 总电荷变化为 {total_charge_change:+.0f}，"
            f"可能需要额外的质子转移步骤"
        )

    valid = len(warnings) == 0 or all(
        "电荷不平衡" in w for w in warnings  # 电荷不平衡仅为警告，不影响合法性
    )

    return {
        "valid": valid,
        "warnings": warnings,
        "charge_imbalance": charge_imbalance,
    }


def generate_reaction_scheme_svg(
    reactant_smiles: str,
    arrows: List[Dict[str, Any]],
    product_smiles: str = "",
) -> str:
    """生成带箭头标注的反应机理 SVG 图。

    使用 RDKit 的 2D 坐标生成原子位置，在反应物分子上绘制箭头指示器，
    展示电子推动箭头的方向。如果未提供产物 SMILES，则通过箭头变换计算。

    Parameters
    ----------
    reactant_smiles : str
        反应物 SMILES 字符串。
    arrows : List[Dict[str, Any]]
        箭头描述列表。
    product_smiles : str
        产物 SMILES 字符串（可选，为空时自动计算）。

    Returns
    -------
    str
        SVG 字符串，包含反应物结构图和箭头标注。
    """
    # 解析反应物
    mol = _smiles_to_mol(reactant_smiles)
    if mol is None:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="50">'
        '</svg>'

    # 分配原子映射编号并生成 2D 坐标
    mol = _assign_atom_map_nums(mol)
    _ensure_2d_coords(mol)

    # 计算产物（如果未提供）
    if not product_smiles:
        result = apply_arrows_to_smiles(reactant_smiles, arrows)
        product_smiles = result.get("product_smiles", "")

    # 解析产物
    prod_mol = _smiles_to_mol(product_smiles) if product_smiles else None

    # 获取反应物 2D 坐标
    reactant_svg = _mol_to_svg_string(mol)

    # 获取产物 SVG
    product_svg = ""
    if prod_mol is not None:
        _ensure_2d_coords(prod_mol)
        product_svg = _mol_to_svg_string(prod_mol)

    # 构建箭头 SVG 标注
    arrow_overlays = _build_arrow_svg_overlays(mol, arrows)

    # 组合完整 SVG
    if product_svg:
        # 双面板：反应物 → 产物
        full_svg = _compose_multi_mol_svg(
            [reactant_svg, arrow_overlays, product_svg],
            labels=["反应物", "箭头", "产物"],
            gap=60,
        )
    else:
        # 单面板：仅反应物 + 箭头标注
        full_svg = _compose_single_with_arrows(reactant_svg, arrow_overlays)

    return full_svg


# ===========================================================================
# SVG 辅助函数
# ===========================================================================

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


def _get_atom_2d_pos(mol: Chem.Mol, idx: int) -> Tuple[float, float]:
    """获取分子中原子的 2D 坐标。

    Parameters
    ----------
    mol : Chem.Mol
        具有构象的分子。
    idx : int
        原子索引。

    Returns
    -------
    Tuple[float, float]
        (x, y) 坐标，无构象时返回 (0, 0)。
    """
    try:
        if mol.GetNumConformers() == 0:
            return (0.0, 0.0)
        conf = mol.GetConformer()
        pos = conf.GetAtomPosition(idx)
        return (pos.x, pos.y)
    except Exception:
        return (0.0, 0.0)


def _build_arrow_svg_overlays(mol: Chem.Mol, arrows: List[Dict[str, Any]]) -> str:
    """构建箭头 SVG 覆盖层。

    在反应物分子上绘制从源原子到目标原子的弧形箭头。

    Parameters
    ----------
    mol : Chem.Mol
        具有 2D 坐标的反应物分子。
    arrows : List[Dict[str, Any]]
        箭头描述列表。

    Returns
    -------
    str
        包含箭头的 SVG 片段字符串。
    """
    if not arrows or mol.GetNumConformers() == 0:
        return ''

    svg_parts: List[str] = []
    svg_parts.append('<g class="arrow-overlays">')

    # 颜色方案：不同箭头使用不同颜色
    colors = [
        '#e74c3c',  # 红色
        '#2980b9',  # 蓝色
        '#27ae60',  # 绿色
        '#f39c12',  # 橙色
        '#8e44ad',  # 紫色
    ]

    # 获取分子边界（用于缩放坐标）
    xs = []
    ys = []
    for atom in mol.GetAtoms():
        x, y = _get_atom_2d_pos(mol, atom.GetIdx())
        xs.append(x)
        ys.append(y)

    if not xs or not ys:
        return ''

    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    range_x = max(max_x - min_x, 1.0)
    range_y = max(max_y - min_y, 1.0)

    # RDKit SVG 坐标映射参数
    svg_width = 300
    svg_height = 250
    padding = 30

    def mol_to_svg_coord(mx: float, my: float) -> Tuple[float, float]:
        """将分子坐标转换为 SVG 坐标。"""
        # RDKit 的坐标系 y 轴向上，SVG y 轴向下
        sx = padding + ((mx - min_x) / range_x) * (svg_width - 2 * padding)
        sy = padding + ((max_y - my) / range_y) * (svg_height - 2 * padding)
        return (sx, sy)

    for i, arrow in enumerate(arrows):
        source_atoms = arrow.get("source_atoms", [])
        target_atoms = arrow.get("target_atoms", [])
        color = colors[i % len(colors)]

        if not source_atoms or not target_atoms:
            continue

        # 计算源原子和目标原子的质心
        src_xs, src_ys = [], []
        for idx in source_atoms:
            x, y = _get_atom_2d_pos(mol, idx)
            src_xs.append(x)
            src_ys.append(y)

        tgt_xs, tgt_ys = [], []
        for idx in target_atoms:
            x, y = _get_atom_2d_pos(mol, idx)
            tgt_xs.append(x)
            tgt_ys.append(y)

        src_cx = sum(src_xs) / len(src_xs)
        src_cy = sum(src_ys) / len(src_ys)
        tgt_cx = sum(tgt_xs) / len(tgt_xs)
        tgt_cy = sum(tgt_ys) / len(tgt_ys)

        sx, sy = mol_to_svg_coord(src_cx, src_cy)
        tx, ty = mol_to_svg_coord(tgt_cx, tgt_cy)

        # 绘制弧形箭头
        electrons = arrow.get("electrons", 2)
        arrow_type = "↿" if electrons == 1 else "⇀"
        line_width = 1.5 if electrons == 2 else 1.0

        # 计算控制点（弧形偏移）
        mid_x = (sx + tx) / 2
        mid_y = (sy + ty) / 2
        dx = tx - sx
        dy = ty - sy
        length = math.sqrt(dx * dx + dy * dy)
        if length < 1:
            continue

        # 垂直偏移量（弧形弯曲方向）
        offset = 20
        perp_x = -dy / length * offset
        perp_y = dx / length * offset
        ctrl_x = mid_x + perp_x
        ctrl_y = mid_y + perp_y

        # SVG path 弧形箭头
        path_d = f"M {sx:.1f},{sy:.1f} Q {ctrl_x:.1f},{ctrl_y:.1f} {tx:.1f},{ty:.1f}"

        # 箭头端点（三角形）
        arrow_size = 6
        angle = math.atan2(ty - ctrl_y, tx - ctrl_x)
        ax1 = tx - arrow_size * math.cos(angle - 0.4)
        ay1 = ty - arrow_size * math.sin(angle - 0.4)
        ax2 = tx - arrow_size * math.cos(angle + 0.4)
        ay2 = ty - arrow_size * math.sin(angle + 0.4)

        svg_parts.append(
            f'<path d="{path_d}" fill="none" stroke="{color}" '
            f'stroke-width="{line_width}" opacity="0.8"/>'
        )
        svg_parts.append(
            f'<polygon points="{tx:.1f},{ty:.1f} {ax1:.1f},{ay1:.1f} '
            f'{ax2:.1f},{ay2:.1f}" fill="{color}" opacity="0.8"/>'
        )

        # 电子数标注
        label_x = ctrl_x + perp_x * 0.3
        label_y = ctrl_y + perp_y * 0.3
        svg_parts.append(
            f'<text x="{label_x:.1f}" y="{label_y:.1f}" '
            f'font-size="10" fill="{color}" text-anchor="middle" '
            f'font-weight="bold">{electrons}e⁻</text>'
        )

    svg_parts.append('</g>')
    return '\n'.join(svg_parts)


def _compose_multi_mol_svg(
    svg_parts: List[str],
    labels: Optional[List[str]] = None,
    gap: int = 60,
) -> str:
    """将多个分子 SVG 片段水平组合为一个完整的 SVG 图。

    Parameters
    ----------
    svg_parts : List[str]
        SVG 片段列表。
    labels : Optional[List[str]]
        每个片段的标签文字。
    gap : int
        片段之间的间距（像素）。

    Returns
    -------
    str
        组合后的完整 SVG 字符串。
    """
    # 过滤空片段
    valid_parts = [p for p in svg_parts if p]
    if not valid_parts:
        return '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="50"></svg>'

    # 提取各部分的宽度
    widths = []
    for part in valid_parts:
        w_match = _extract_svg_width(part)
        widths.append(w_match if w_match > 0 else 300)

    total_width = sum(widths) + gap * (len(valid_parts) - 1)
    total_width = max(total_width, 200)
    height = 300

    result_parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_width}" height="{height}">',
        f'<rect width="100%" height="100%" fill="white"/>',
    ]

    x_offset = 0
    for i, (part, w) in enumerate(zip(valid_parts, widths)):
        # 提取内部 SVG 内容（去掉外层 <svg> 标签）
        inner = _strip_svg_tags(part)
        result_parts.append(f'<g transform="translate({x_offset}, 0)">')
        result_parts.append(inner)
        result_parts.append('</g>')

        # 添加标签
        if labels and i < len(labels) and labels[i]:
            label_x = x_offset + w / 2
            result_parts.append(
                f'<text x="{label_x}" y="{height - 5}" font-size="12" '
                f'fill="#666" text-anchor="middle">{labels[i]}</text>'
            )

        # 添加箭头符号（片段之间）
        if i < len(valid_parts) - 1:
            arrow_x = x_offset + w + gap / 2
            result_parts.append(
                f'<text x="{arrow_x}" y="{height / 2}" font-size="24" '
                f'fill="#333" text-anchor="middle">→</text>'
            )

        x_offset += w + gap

    result_parts.append('</svg>')
    return '\n'.join(result_parts)


def _compose_single_with_arrows(mol_svg: str, arrow_svg: str) -> str:
    """将分子 SVG 和箭头覆盖层组合。

    Parameters
    ----------
    mol_svg : str
        分子 SVG 字符串。
    arrow_svg : str
        箭头覆盖层 SVG 字符串。

    Returns
    -------
    str
        组合后的 SVG 字符串。
    """
    width = _extract_svg_width(mol_svg) or 300
    height = 300

    inner_mol = _strip_svg_tags(mol_svg)

    result = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        f'<rect width="100%" height="100%" fill="white"/>',
        inner_mol,
    ]
    if arrow_svg:
        result.append(arrow_svg)
    result.append('</svg>')

    return '\n'.join(result)


def _extract_svg_width(svg_text: str) -> int:
    """从 SVG 文本中提取宽度。

    Parameters
    ----------
    svg_text : str
        SVG 字符串。

    Returns
    -------
    int
        SVG 宽度，提取失败返回 0。
    """
    import re
    match = re.search(r'width="(\d+)"', svg_text)
    if match:
        return int(match.group(1))
    return 0


def _strip_svg_tags(svg_text: str) -> str:
    """去除 SVG 的外层 <svg> 和 </svg> 标签，保留内部内容。

    Parameters
    ----------
    svg_text : str
        SVG 字符串。

    Returns
    -------
    str
        去除外层标签后的内部内容。
    """
    import re
    # 去除开头的 <svg ...> 标签
    text = re.sub(r'<svg[^>]*>', '', svg_text, count=1)
    # 去除结尾的 </svg> 标签
    text = re.sub(r'</svg>', '', text, count=1)
    return text.strip()


# ===========================================================================
# 模块级快捷入口（供 mechanism_search.py 集成使用）
# ===========================================================================

def get_arrow_info_for_step(reaction_smarts: str) -> Dict[str, Any]:
    """为机制搜索的某个步骤获取箭头信息（供 mechanism_search.py 集成使用）。

    Parameters
    ----------
    reaction_smarts : str
        该步骤的反应 SMARTS 字符串。

    Returns
    -------
    Dict[str, Any]
        包含 arrows 和 arrow_svg 的字典：
        - ``arrows`` (List[Dict]): 推断的箭头序列
        - ``arrow_svg`` (str): 箭头 SVG 图
    """
    arrows = parse_reaction_to_arrows(reaction_smarts)

    # 生成简要 SVG（仅显示箭头信息，无需完整分子图）
    arrow_svg = ''
    if arrows:
        arrow_svg = _build_arrow_summary_svg(arrows)

    return {
        "arrows": arrows,
        "arrow_svg": arrow_svg,
    }


def _build_arrow_summary_svg(arrows: List[Dict[str, Any]]) -> str:
    """构建箭头序列的摘要 SVG。

    Parameters
    ----------
    arrows : List[Dict[str, Any]]
        箭头描述列表。

    Returns
    -------
    str
        摘要 SVG 字符串。
    """
    width = 200
    height = max(40, 30 + len(arrows) * 22)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="#f8f9fa" rx="4"/>',
    ]

    colors = ['#e74c3c', '#2980b9', '#27ae60', '#f39c12', '#8e44ad']

    for i, arrow in enumerate(arrows):
        y = 20 + i * 22
        source = arrow.get("source_atoms", [])
        target = arrow.get("target_atoms", [])
        electrons = arrow.get("electrons", 2)
        desc = arrow.get("description", "")
        color = colors[i % len(colors)]

        # 箭头类型图标
        if len(source) > 1:
            icon = "⇌"  # 键断裂/重排
        else:
            icon = "→"  # 电子转移

        text = f"{icon} {desc or f'{source}→{target} ({electrons}e⁻)'}"
        parts.append(
            f'<text x="10" y="{y}" font-size="11" fill="{color}" '
            f'font-family="monospace">{text}</text>'
        )

    parts.append('</svg>')
    return '\n'.join(parts)
