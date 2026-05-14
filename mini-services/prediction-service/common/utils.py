"""
EzMechanism 预测系统通用工具函数

本模块从论文 EzMechanism 的源码 (common/utils.py) 中提取并改编而来，
去除了对 Django 框架的依赖，使所有函数可独立使用。

主要功能包括：
- 分子合并与标准化（combine_mols_list, make_mol_equivalent）
- 集合合并与 Union-Find（consolidate_sets）
- SMARTS 字符串处理（remove_atom_mapping_from_smarts）
- 分子环境提取（get_mol_env）
- 电子推动箭头变换（get_mol_product_from_arrow, get_mol_product_from_props）
- 原子查询重置与标签生成（reset_atom_query, label_int_to_letters）

依赖：rdkit, re, typing
"""

from typing import List, Set, Tuple, Optional, FrozenSet
from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem.rdchem import RWMol, BondType, Atom
import re


# ---------------------------------------------------------------------------
# 1. combine_mols_list — 合并多个 RDKit Mol 对象为一个分子
# ---------------------------------------------------------------------------

def combine_mols_list(mols: List[Chem.Mol]) -> Chem.Mol:
    """将多个 RDKit Mol 对象合并为一个分子。

    通过 ``Chem.CombineMols`` 迭代合并 *mols* 中的所有分子。
    若 *mols* 为空或仅含 ``None``，则返回空分子。

    Parameters
    ----------
    mols : List[Chem.Mol]
        待合并的分子列表。

    Returns
    -------
    Chem.Mol
        包含输入列表中所有原子与键的单一分子。

    示例
    -----
    >>> from rdkit import Chem
    >>> m1 = Chem.MolFromSmiles('C')
    >>> m2 = Chem.MolFromSmiles('O')
    >>> combined = combine_mols_list([m1, m2])
    >>> combined.GetNumAtoms()
    2
    """
    valid_mols = [m for m in mols if m is not None]
    if not valid_mols:
        return Chem.Mol()  # 空分子
    combined = valid_mols[0]
    for mol in valid_mols[1:]:
        combined = Chem.CombineMols(combined, mol)
    return combined


# ---------------------------------------------------------------------------
# 2. make_mol_equivalent — 标准化分子用于去重
# ---------------------------------------------------------------------------

def make_mol_equivalent(mol: Chem.Mol) -> Chem.Mol:
    """标准化分子，使连接性等价的分子产生相同的 SMILES 字符串。

    依次执行：清除原子映射编号、移除显式氢、分子净化（sanitize）。
    主要用于分子去重，而非保留化学语义。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子。

    Returns
    -------
    Chem.Mol
        去除显式氢与原子映射编号后的标准化副本。
        若处理失败则返回空分子。
    """
    try:
        # 先清除原子映射编号
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        # 移除显式氢
        mol_no_h = Chem.RemoveHs(mol, sanitize=False)
        # 净化分子使其处于合法状态
        Chem.SanitizeMol(mol_no_h)
        return mol_no_h
    except Exception:
        return Chem.Mol()


# ---------------------------------------------------------------------------
# 3. consolidate_sets — 基于并查集（Union-Find）的集合合并
# ---------------------------------------------------------------------------

def consolidate_sets(list_of_sets: List[Set[int]]) -> List[Set[int]]:
    """将重叠或相连的集合合并为不相交的组。

    使用并查集（Disjoint Set Union）数据结构，高效地将
    在任意输入集合中共同出现的元素归为一组。

    Parameters
    ----------
    list_of_sets : List[Set[int]]
        需按连接性分组的集合列表。

    Returns
    -------
    List[Set[int]]
        合并后的不相交集合列表，每个集合代表一个连通组。

    示例
    -----
    >>> consolidate_sets([{1, 2}, {2, 3}, {5}])
    [{1, 2, 3}, {5}]
    """
    if not list_of_sets:
        return []

    # -- 并查集实现 --
    parent: dict = {}

    def find(x: int) -> int:
        """带路径压缩的查找操作。"""
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        """合并操作。"""
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry

    # 收集所有元素并在每个集合内部执行合并
    for s in list_of_sets:
        elements = list(s)
        if not elements:
            continue
        for e in elements:
            parent.setdefault(e, e)
        first = elements[0]
        for e in elements[1:]:
            union(first, e)

    # 构建合并后的分组
    groups: dict = {}
    for elem in parent:
        root = find(elem)
        groups.setdefault(root, set()).add(elem)

    return list(groups.values())


# ---------------------------------------------------------------------------
# 4. remove_atom_mapping_from_smarts — 去除 SMARTS 中的原子映射编号
# ---------------------------------------------------------------------------

def remove_atom_mapping_from_smarts(smarts: str) -> str:
    """去除 SMARTS 字符串中的原子映射编号（``:N`` 模式）。

    仅移除出现在右括号 ``]`` 或原子符号之后的映射标签，
    ``[#6]`` 等原子查询模式不会被影响。

    Parameters
    ----------
    smarts : str
        输入的 SMARTS 字符串。

    Returns
    -------
    str
        去除原子映射编号后的 SMARTS 字符串。

    示例
    -----
    >>> remove_atom_mapping_from_smarts('[CH3:1][CH2:2][OH:3]')
    '[CH3][CH2][OH]'
    """
    # 去除 ':数字' 模式（原子映射编号）
    return re.sub(r':\d+', '', smarts)


# ---------------------------------------------------------------------------
# 5. get_mol_env — 通过 BFS 提取指定原子映射编号周围的分子环境
# ---------------------------------------------------------------------------

def get_mol_env(
    mol: Chem.Mol,
    atom_map_nums: Set[int],
    radius: int = 1,
) -> Tuple[Optional[Chem.Mol], Set[int]]:
    """提取指定原子映射编号周围的分子环境子结构。

    从 *atom_map_nums* 中的原子出发，通过 BFS 向外扩展 *radius* 个键，
    收集所有到达的原子索引，返回对应的子结构分子（含查询原子）。

    Parameters
    ----------
    mol : Chem.Mol
        输入分子（原子必须已设置映射编号）。
    atom_map_nums : Set[int]
        定义环境中心的原子映射编号集合。
    radius : int
        向外扩展的键的层数（默认为 1）。

    Returns
    -------
    Tuple[Optional[Chem.Mol], Set[int]]
        ``(env_mol, included_indices)``
        - *env_mol*：提取的子结构（原子为查询原子），若未找到则返回 None。
        - *included_indices*：被包含的原始原子索引集合。
    """
    if mol is None or not atom_map_nums:
        return None, set()

    # 构建原子映射编号 → 原子索引的映射
    am_to_idx: dict = {}
    for atom in mol.GetAtoms():
        amn = atom.GetAtomMapNum()
        if amn != 0:
            am_to_idx[amn] = atom.GetIdx()

    # 确定起始原子索引
    start_indices = {am_to_idx[am] for am in atom_map_nums if am in am_to_idx}
    if not start_indices:
        return None, set()

    # BFS 扩展
    visited: Set[int] = set(start_indices)
    frontier = set(start_indices)
    for _ in range(radius):
        next_frontier: Set[int] = set()
        for idx in frontier:
            atom = mol.GetAtomWithIdx(idx)
            for neighbour in atom.GetNeighbors():
                n_idx = neighbour.GetIdx()
                if n_idx not in visited:
                    visited.add(n_idx)
                    next_frontier.add(n_idx)
        frontier = next_frontier

    # 提取子结构（保留原子映射编号）
    env_indices = sorted(visited)
    env_mol = Chem.PathToSubmol(mol, env_indices)

    # 将所有原子转换为查询原子，以支持灵活的 SMARTS 匹配。
    # 常规 Atom 对象没有 SetQuery() 方法，
    # 因此需要使用 RWMol.ReplaceAtom() 将每个原子替换为 QueryAtom。
    rw_env = RWMol(env_mol)
    for atom in list(rw_env.GetAtoms()):
        atomic_num = atom.GetAtomicNum()
        map_num = atom.GetAtomMapNum()
        if atomic_num > 0:
            q_atom = Chem.rdqueries.AtomNumEqualsQueryAtom(atomic_num)
            q_atom.SetAtomMapNum(map_num)
            rw_env.ReplaceAtom(atom.GetIdx(), q_atom)

    try:
        Chem.SanitizeMol(rw_env)
    except Exception:
        pass

    return rw_env.GetMol(), visited


# ---------------------------------------------------------------------------
# 6. get_mol_product_from_arrow — 应用单个电子推动箭头变换
# ---------------------------------------------------------------------------

def get_mol_product_from_arrow(
    mol: Chem.Mol,
    arrow_info: Tuple[FrozenSet[int], FrozenSet[int], int],
) -> Chem.Mol:
    """应用单个电子推动箭头（curly arrow）变换，生成产物分子。

    Parameters
    ----------
    mol : Chem.Mol
        反应物分子（原子需设置映射编号）。
    arrow_info : Tuple[FrozenSet[int], FrozenSet[int], int]
        ``(source_atom_map_nums, target_atom_map_nums, electron_count)``
        - *source_atom_map_nums*：箭头源原子映射编号集合
        - *target_atom_map_nums*：箭头目标原子映射编号集合
        - *electron_count*：电子数（1 为自由基，2 为电子对）

    Returns
    -------
    Chem.Mol
        应用箭头变换后的产物分子。

    Notes
    -----
    对于 2-电子箭头（从源原子 A 到目标原子 B）：

    * A-B 键**已存在** → 键级增加 1，并调整形式电荷。
    * A-B 键**不存在** → 形成新的单键 A-B，并调整形式电荷。

    电荷守恒规则：
      - 源原子电荷增加 electron_count（失去电子）。
      - 目标原子电荷减少 electron_count（获得电子）。

    对于 electron_count == 1（自由基），形式电荷不变，
    但键级仍会调整。
    """
    source_ams, target_ams, electron_count = arrow_info

    rwmol = Chem.RWMol(mol)

    # 构建原子映射编号 → 原子索引的查找表
    am_to_idx: dict = {}
    for atom in rwmol.GetAtoms():
        amn = atom.GetAtomMapNum()
        if amn != 0:
            am_to_idx[amn] = atom.GetIdx()

    for s_am in source_ams:
        for t_am in target_ams:
            s_idx = am_to_idx.get(s_am)
            t_idx = am_to_idx.get(t_am)
            if s_idx is None or t_idx is None:
                continue

            bond = rwmol.GetBondBetweenAtoms(s_idx, t_idx)
            if bond is not None:
                # 已有键 → 修改键级
                current_order = bond.GetBondTypeAsDouble()
                new_order = current_order + electron_count
                if new_order <= 0:
                    # 键完全断裂
                    rwmol.RemoveBond(s_idx, t_idx)
                else:
                    # 将键级映射到最接近的标准键类型
                    if new_order >= 3:
                        bond.SetBondType(BondType.TRIPLE)
                    elif new_order >= 2:
                        bond.SetBondType(BondType.DOUBLE)
                    elif new_order >= 1.5:
                        bond.SetBondType(BondType.AROMATIC)
                    else:
                        bond.SetBondType(BondType.SINGLE)
            else:
                # 无键 → 形成新的单键
                rwmol.AddBond(s_idx, t_idx, BondType.SINGLE)

            # 调整形式电荷（仅对成对电子生效）
            if electron_count == 2:
                s_atom = rwmol.GetAtomWithIdx(s_idx)
                t_atom = rwmol.GetAtomWithIdx(t_idx)
                s_atom.SetFormalCharge(s_atom.GetFormalCharge() + electron_count)
                t_atom.SetFormalCharge(t_atom.GetFormalCharge() - electron_count)

    # 尝试净化；中间体结构可能不完全合法，静默跳过
    try:
        Chem.SanitizeMol(rwmol)
    except Exception:
        pass

    return rwmol.GetMol()


# ---------------------------------------------------------------------------
# 7. get_mol_product_from_props — 通过原子属性应用反应变换
# ---------------------------------------------------------------------------

def get_mol_product_from_props(
    mol: Chem.Mol,
) -> Tuple[Chem.Mol, List[Tuple[FrozenSet[int], FrozenSet[int], int]]]:
    """通过原子属性编码的反应变换，生成产物分子。

    分子中的原子携带箭头相关的整数属性，命名约定如下：

    * ``arrow_<N>_out``   — 箭头 N 的源原子标记（1 表示"此原子为源"）
    * ``arrow_<N>_in``    — 箭头 N 的目标原子标记
    * ``arrow_<N>_count`` — 箭头 N 的电子数（1 或 2）

    此外还读取以下属性（仅读取，不修改变换逻辑）：
    ``rc``, ``env_term``, ``charge_plus``, ``charge_minus``

    Parameters
    ----------
    mol : Chem.Mol
        带有反应属性注释的分子。

    Returns
    -------
    Tuple[Chem.Mol, List[Tuple[FrozenSet[int], FrozenSet[int], int]]]
        ``(product_mol, arrow_info_list)``
        - *product_mol*：变换后的产物分子
        - *arrow_info_list*：解析后的箭头描述列表，
          每项为 ``(source_ams, target_ams, electron_count)``
    """
    if mol is None:
        return Chem.Mol(), []

    # --- 发现箭头属性 ---
    # 扫描所有原子，查找匹配 arrow_N_out / arrow_N_in 的属性键
    arrow_numbers: Set[int] = set()
    for atom in mol.GetAtoms():
        for prop_name in atom.GetPropNames():
            if prop_name.startswith("arrow_") and prop_name.endswith("_out"):
                # 从 "arrow_N_out" 中提取 N
                parts = prop_name.split("_")
                if len(parts) == 3:
                    try:
                        arrow_numbers.add(int(parts[1]))
                    except ValueError:
                        pass

    # --- 构建箭头信息列表 ---
    arrow_info_list: List[Tuple[FrozenSet[int], FrozenSet[int], int]] = []
    for n in sorted(arrow_numbers):
        source_ams: Set[int] = set()
        target_ams: Set[int] = set()
        electron_count = 2  # 默认为 2 电子

        for atom in mol.GetAtoms():
            out_key = f"arrow_{n}_out"
            in_key = f"arrow_{n}_in"
            count_key = f"arrow_{n}_count"

            if out_key in atom.GetPropNames():
                val = atom.GetIntProp(out_key)
                if val:
                    amn = atom.GetAtomMapNum()
                    if amn != 0:
                        source_ams.add(amn)

            if in_key in atom.GetPropNames():
                val = atom.GetIntProp(in_key)
                if val:
                    amn = atom.GetAtomMapNum()
                    if amn != 0:
                        target_ams.add(amn)

            if count_key in atom.GetPropNames():
                electron_count = atom.GetIntProp(count_key)

        if source_ams and target_ams:
            arrow_info_list.append(
                (frozenset(source_ams), frozenset(target_ams), electron_count)
            )

    # --- 应用属性中的电荷调整 ---
    rwmol = Chem.RWMol(mol)
    for atom in rwmol.GetAtoms():
        props = atom.GetPropNames()
        if "charge_plus" in props:
            atom.SetFormalCharge(
                atom.GetFormalCharge() + atom.GetIntProp("charge_plus")
            )
        if "charge_minus" in props:
            atom.SetFormalCharge(
                atom.GetFormalCharge() - atom.GetIntProp("charge_minus")
            )

    # --- 依次应用每个箭头变换 ---
    for arrow_info in arrow_info_list:
        source_ams, target_ams, electron_count = arrow_info

        # 为 rwmol 的当前状态构建映射编号 → 索引查找表
        am_to_idx: dict = {}
        for atom in rwmol.GetAtoms():
            amn = atom.GetAtomMapNum()
            if amn != 0:
                am_to_idx[amn] = atom.GetIdx()

        for s_am in source_ams:
            for t_am in target_ams:
                s_idx = am_to_idx.get(s_am)
                t_idx = am_to_idx.get(t_am)
                if s_idx is None or t_idx is None:
                    continue

                # 在前一次键修改之后重新计算索引
                s_idx = am_to_idx.get(s_am, s_idx)
                t_idx = am_to_idx.get(t_am, t_idx)

                bond = rwmol.GetBondBetweenAtoms(s_idx, t_idx)
                if bond is not None:
                    current_order = bond.GetBondTypeAsDouble()
                    new_order = current_order + electron_count
                    if new_order <= 0:
                        rwmol.RemoveBond(s_idx, t_idx)
                    else:
                        if new_order >= 3:
                            bond.SetBondType(BondType.TRIPLE)
                        elif new_order >= 2:
                            bond.SetBondType(BondType.DOUBLE)
                        elif new_order >= 1.5:
                            bond.SetBondType(BondType.AROMATIC)
                        else:
                            bond.SetBondType(BondType.SINGLE)
                else:
                    rwmol.AddBond(s_idx, t_idx, BondType.SINGLE)

    try:
        Chem.SanitizeMol(rwmol)
    except Exception:
        pass

    return rwmol.GetMol(), arrow_info_list


# ---------------------------------------------------------------------------
# 8. reset_atom_query — 重置原子查询以匹配特定元素
# ---------------------------------------------------------------------------

def reset_atom_query(atom: Atom, atomic_num: int) -> None:
    """重置 RDKit 原子的查询条件，使其仅匹配指定的原子序数。

    当原子原本是查询原子（例如 ``GetAtomicNum() == 0``）时，
    可使用此函数将其固定为特定元素。

    Parameters
    ----------
    atom : Atom
        要修改的 RDKit 原子对象（就地修改，无返回值）。
    atomic_num : int
        要匹配的原子序数（例如 6 表示碳）。
    """
    atom.SetQuery(Chem.rdqueries.AtomNumEqualsQueryAtom(atomic_num))


# ---------------------------------------------------------------------------
# 9. label_int_to_letters — 将整数转换为 Excel 风格的列字母标签
# ---------------------------------------------------------------------------

def label_int_to_letters(n: int) -> str:
    """将正整数转换为 Excel 风格的列字母（小写）。

    Parameters
    ----------
    n : int
        1-based 索引。

    Returns
    -------
    str
        对应的小写字母标签。

    示例
    -----
    >>> label_int_to_letters(1)
    'a'
    >>> label_int_to_letters(26)
    'z'
    >>> label_int_to_letters(27)
    'aa'
    >>> label_int_to_letters(28)
    'ab'
    """
    if n <= 0:
        return ""
    result = ""
    while n > 0:
        n -= 1
        result = chr(97 + (n % 26)) + result
        n //= 26
    return result
