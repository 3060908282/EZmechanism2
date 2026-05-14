# -*- coding: utf-8 -*-
"""
化学检测工具模块 —— 氨基酸编码映射、分子指纹相似度、子结构匹配与反应类型检测。

本模块从论文原始代码（Django 版）中剥离而来，去除了所有 Django 依赖，
改为接收纯 RDKit Mol 对象和原子映射字典，可直接在独立 Python 环境中运行。

功能概览
--------
1. 氨基酸三字母 / 单字母编码双向映射（AA_1TO3, AA_3TO1, AA_20_CODE3）
2. MD5 哈希（get_md5）
3. Tanimoto 指纹相似度（fingerprint_similarity）
4. 最大模板子结构匹配（get_maximum_template_substructure）
5. 四种反应中心类型检测：
   - 异裂键断裂（_is_heterolytic_bond_cleavage）
   - 质子转移（_is_proton_transfer）
   - 亲核攻击（_is_nucleophilic_attack）
   - 自由基反应（_is_radical_reaction）
6. 方法 → 残基角色标签注册表（METHOD_TO_RESIDUE_FUNCTIONS）
"""

import hashlib
from typing import Callable, Dict, List, Optional, Tuple

from rdkit import Chem


# ===========================================================================
# 氨基酸编码映射
# ===========================================================================
# 包含 20 种标准氨基酸 + 5 种常见特殊编码（B/Z/X/U/O），用于将蛋白质序列
# 中的单字母缩写与三字母缩写进行相互转换。

AA_1TO3: Dict[str, str] = {
    "A": "Ala",  # 丙氨酸 Alanine
    "R": "Arg",  # 精氨酸 Arginine
    "N": "Asn",  # 天冬酰胺 Asparagine
    "D": "Asp",  # 天冬氨酸 Aspartic acid
    "C": "Cys",  # 半胱氨酸 Cysteine
    "E": "Glu",  # 谷氨酸 Glutamic acid
    "Q": "Gln",  # 谷氨酰胺 Glutamine
    "G": "Gly",  # 甘氨酸 Glycine
    "H": "His",  # 组氨酸 Histidine
    "I": "Ile",  # 异亮氨酸 Isoleucine
    "L": "Leu",  # 亮氨酸 Leucine
    "K": "Lys",  # 赖氨酸 Lysine
    "M": "Met",  # 甲硫氨酸 Methionine
    "F": "Phe",  # 苯丙氨酸 Phenylalanine
    "P": "Pro",  # 脯氨酸 Proline
    "S": "Ser",  # 丝氨酸 Serine
    "T": "Thr",  # 苏氨酸 Threonine
    "W": "Trp",  # 色氨酸 Tryptophan
    "Y": "Tyr",  # 酪氨酸 Tyrosine
    "V": "Val",  # 缬氨酸 Valine
    # ---- 常见特殊 / 模糊编码 ----
    "B": "Asx",  # 天冬氨酸或天冬酰胺（Asp 或 Asn）
    "Z": "Glx",  # 谷氨酸或谷氨酰胺（Glu 或 Gln）
    "X": "Xaa",  # 未知氨基酸（任意）
    "U": "Sec",  # 硒半胱氨酸 Selenocysteine
    "O": "Pyl",  # 吡咯赖氨酸 Pyrrolysine
}

AA_3TO1: Dict[str, str] = {v: k for k, v in AA_1TO3.items()}

AA_20_CODE3: List[str] = [
    "Ala", "Arg", "Asn", "Asp", "Cys",
    "Glu", "Gln", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tyr", "Val",
]


# ===========================================================================
# 通用辅助函数
# ===========================================================================


def get_md5(string: str) -> str:
    """返回字符串的 MD5 十六进制摘要值。

    用于为序列数据生成唯一标识符。

    参数
    ----
    string : str
        待哈希的输入字符串。

    返回
    ----
    str
        32 位十六进制 MD5 摘要。
    """
    return hashlib.md5(string.encode()).hexdigest()


def fingerprint_similarity(fp1, fp2) -> float:
    """计算两个分子指纹之间的 Tanimoto 相似度。

    支持三种输入类型：

    - **计数 / 列表向量**（``list[int]``）——使用 min/max Tanimoto 公式，
      适用于计数型指纹（如 ``bond_changes_fingerprint``）。
    - **RDKit 二进制指纹字符串** ——通过 ``DataStructs.cDataStructs``
      反序列化后调用 RDKit 原生 Tanimoto 实现。
    - **其他可迭代对象** ——按照计数型向量处理。

    当并集为零（即两个指纹全为 0）时返回 0.0。

    参数
    ----
    fp1, fp2 :
        待比较的指纹。可为列表、元组、NumPy 数组或 RDKit 二进制指纹字符串。

    返回
    ----
    float
        Tanimoto 相似度，范围 [0.0, 1.0]。

    异常
    ----
    TypeError
        当输入类型无法处理时抛出。
    """
    # ---- RDKit 序列化二进制指纹字符串 ----
    if isinstance(fp1, str) and isinstance(fp2, str):
        try:
            from rdkit import DataStructs
            rdkit_fp1 = DataStructs.cDataStructs.CreateFromBinaryText(fp1)
            rdkit_fp2 = DataStructs.cDataStructs.CreateFromBinaryText(fp2)
            return DataStructs.TanimotoSimilarity(rdkit_fp1, rdkit_fp2)
        except Exception:
            # 反序列化失败时回退到通用处理
            pass

    # ---- 通用计数型 Tanimoto ----
    # 适用于列表、元组、NumPy 数组等可迭代对象。
    try:
        iterator = zip(fp1, fp2)
    except TypeError:
        raise TypeError(
            "fingerprint_similarity 需要可迭代的指纹向量 "
            "或 RDKit 二进制指纹字符串作为输入"
        )

    intersection = sum(min(a, b) for a, b in iterator)
    union = sum(max(a, b) for a, b in zip(fp1, fp2))
    return intersection / union if union > 0 else 0.0


# ===========================================================================
# RDKit 子结构匹配辅助函数
# ===========================================================================


def get_maximum_template_substructure(
    smarts: str,
    mol: Chem.Mol,
) -> Optional[Tuple[int, ...]]:
    """返回 SMARTS 模板在分子中首次（最大）匹配的原子索引元组。

    RDKit 的 ``GetSubstructMatches`` 会按照匹配原子数降序返回结果，
    因此第一个匹配即为最大匹配。

    参数
    ----
    smarts : str
        代表残基模板的 SMARTS 模式字符串。
    mol : Chem.Mol
        待搜索的 RDKit 分子对象。

    返回
    ----
    tuple[int, ...] | None
        首次匹配的原子索引元组；若未找到匹配则返回 None。
    """
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return None
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return None
    return tuple(matches[0])


# ===========================================================================
# 箭头反应检测函数（无 Django 依赖版本）
# ===========================================================================
# 以下函数用于检测弯箭头（curly-arrow）机制图中的特定反应类型。
# 原始代码依赖 Django Step 对象（通过 step.am_to_mol 获取分子结构），
# 本适配版本改为直接接收 RDKit Mol 对象和原子映射字典，完全独立于 Django。
#
# 调用约定
# --------
# reaction_info : tuple
#     三元组 (source_atom_maps: List[int], target_atom_maps: List[int],
#             electron_count: int)，描述弯箭头的起点、终点和电子数。
# mol : Chem.Mol | None
#     RDKit 分子对象。若为 None，则尝试从 reaction_info 中提取或返回 None。
# am_to_idx : dict | None
#     原子映射编号 → 原子索引的字典。若为 None，则通过遍历 mol 的原子映射
#     编号来动态解析。
#
# 返回值
# -------
# tuple[str, str] | None
#     成功时返回 (源标签, 目标标签) 的二元组；失败时返回 None。


def _resolve_atom_map(mol: Chem.Mol, am_value: int) -> Optional[int]:
    """在分子中查找指定原子映射编号对应的原子索引。

    参数
    ----
    mol : Chem.Mol
        RDKit 分子对象。
    am_value : int
        原子映射编号。

    返回
    ----
    int | None
        对应的原子索引；未找到时返回 None。
    """
    for atom in mol.GetAtoms():
        if atom.GetAtomMapNum() == am_value:
            return atom.GetIdx()
    return None


def _extract_mol_from_reaction_info(reaction_info) -> Optional[Chem.Mol]:
    """尝试从 reaction_info 中提取 RDKit Mol 对象（兼容旧调用方式）。

    如果 reaction_info 是长度 > 3 的元组且最后一个是 Mol 对象，则返回之。
    否则返回 None。

    参数
    ----
    reaction_info : tuple
        反应信息元组。

    返回
    ----
    Chem.Mol | None
    """
    if isinstance(reaction_info, (list, tuple)) and len(reaction_info) > 3:
        last_item = reaction_info[-1]
        if isinstance(last_item, Chem.Mol):
            return last_item
    return None


def _is_heterolytic_bond_cleavage(
    reaction_info,
    mol: Optional[Chem.Mol] = None,
    am_to_idx: Optional[dict] = None,
) -> Optional[Tuple[str, str]]:
    """检测弯箭头反应中的异裂键断裂。

    异裂断裂是指共价键断裂后，两个电子全部归入其中一个原子（离去基团）。
    典型特征：电子数 = 2，且箭头起止两端存在化学键。

    参数
    ----
    reaction_info : tuple
        (source_atom_maps, target_atom_maps, electron_count) 三元组。
    mol : Chem.Mol | None
        RDKit 分子对象。若为 None，尝试从 reaction_info 中提取。
    am_to_idx : dict | None
        原子映射编号 → 原子索引的字典。若为 None，自动遍历 mol 查找。

    返回
    -------
    tuple[str, str] | None
        成功时返回 (电正性基团标签, 电负性基团标签)；失败返回 None。
    """
    if not isinstance(reaction_info, (list, tuple)) or len(reaction_info) < 3:
        return None

    source_ams, target_ams, e_count = reaction_info[0], reaction_info[1], reaction_info[2]

    # 异裂断裂必然涉及 2 个电子（一个共价键）
    if e_count != 2:
        return None

    # 尝试获取分子对象
    if mol is None:
        mol = _extract_mol_from_reaction_info(reaction_info)
    if mol is None:
        return None

    # 遍历 source → target 原子映射对，查找是否存在即将断裂的化学键
    try:
        for s_am in source_ams:
            for t_am in target_ams:
                # 解析原子映射编号 → 原子索引
                if am_to_idx is not None:
                    s_idx = am_to_idx.get(s_am)
                    t_idx = am_to_idx.get(t_am)
                else:
                    s_idx = _resolve_atom_map(mol, s_am)
                    t_idx = _resolve_atom_map(mol, t_am)

                if s_idx is not None and t_idx is not None:
                    bond = mol.GetBondBetweenAtoms(s_idx, t_idx)
                    if bond is not None:
                        return (f"a{s_am}", f"a{t_am}")
    except Exception:
        # 尽力检测：吞掉异常并返回未检测到
        pass

    return None


def _is_proton_transfer(
    reaction_info,
    mol: Optional[Chem.Mol] = None,
    am_to_idx: Optional[dict] = None,
) -> Optional[Tuple[str, str]]:
    """检测质子转移事件（H⁺ 在两个杂原子之间移动）。

    质子转移是一种 2 电子过程，特征为箭头一端是氢原子（原子序数 = 1），
    另一端是重原子（原子序数 > 1）。

    参数
    ----
    reaction_info : tuple
        (source_atom_maps, target_atom_maps, electron_count) 三元组。
    mol : Chem.Mol | None
        RDKit 分子对象。若为 None，尝试从 reaction_info 中提取。
    am_to_idx : dict | None
        原子映射编号 → 原子索引的字典。若为 None，自动遍历 mol 查找。

    返回
    -------
    tuple[str, str] | None
        成功时返回 (质子给体标签, 质子受体标签)；失败返回 None。
    """
    if not isinstance(reaction_info, (list, tuple)) or len(reaction_info) < 3:
        return None

    source_ams, target_ams, e_count = reaction_info[0], reaction_info[1], reaction_info[2]

    # 质子转移为 2 电子过程
    if e_count != 2:
        return None

    if mol is None:
        mol = _extract_mol_from_reaction_info(reaction_info)
    if mol is None:
        return None

    try:
        for s_am in source_ams:
            for t_am in target_ams:
                if am_to_idx is not None:
                    s_idx = am_to_idx.get(s_am)
                    t_idx = am_to_idx.get(t_am)
                else:
                    s_idx = _resolve_atom_map(mol, s_am)
                    t_idx = _resolve_atom_map(mol, t_am)

                if s_idx is None or t_idx is None:
                    continue

                s_atom = mol.GetAtomWithIdx(s_idx)
                t_atom = mol.GetAtomWithIdx(t_idx)
                s_is_h = s_atom.GetAtomicNum() == 1
                t_is_h = t_atom.GetAtomicNum() == 1

                # 一端是氢原子，另一端是重原子 —— 符合质子转移特征
                if s_is_h != t_is_h:
                    return (f"a{s_am}", f"a{t_am}")
    except Exception:
        pass

    return None


def _is_nucleophilic_attack(
    reaction_info,
    mol: Optional[Chem.Mol] = None,
    am_to_idx: Optional[dict] = None,
) -> Optional[Tuple[str, str]]:
    """检测亲核攻击（孤对电子形成新键）。

    亲核攻击是 2 电子过程，特征为富电子原子（N/O/S 等高电负性原子）
    的孤对电子推入缺电子中心（C/P/Si 等低电负性原子）。

    参数
    ----
    reaction_info : tuple
        (source_atom_maps, target_atom_maps, electron_count) 三元组。
    mol : Chem.Mol | None
        RDKit 分子对象。若为 None，尝试从 reaction_info 中提取。
    am_to_idx : dict | None
        原子映射编号 → 原子索引的字典。若为 None，自动遍历 mol 查找。

    返回
    -------
    tuple[str, str] | None
        成功时返回 (亲核试剂标签, 亲电试剂标签)；失败返回 None。
    """
    if not isinstance(reaction_info, (list, tuple)) or len(reaction_info) < 3:
        return None

    source_ams, target_ams, e_count = reaction_info[0], reaction_info[1], reaction_info[2]

    # 亲核攻击为 2 电子过程
    if e_count != 2:
        return None

    if mol is None:
        mol = _extract_mol_from_reaction_info(reaction_info)
    if mol is None:
        return None

    try:
        for s_am in source_ams:
            for t_am in target_ams:
                if am_to_idx is not None:
                    s_idx = am_to_idx.get(s_am)
                    t_idx = am_to_idx.get(t_am)
                else:
                    s_idx = _resolve_atom_map(mol, s_am)
                    t_idx = _resolve_atom_map(mol, t_am)

                if s_idx is None or t_idx is None:
                    continue

                s_atom = mol.GetAtomWithIdx(s_idx)
                t_atom = mol.GetAtomWithIdx(t_idx)

                # 亲核原子：高电负性元素（N=7, O=8, S=16）
                # 亲电原子：低电负性元素（C=6, P=15, Si=14）
                s_electroneg = s_atom.GetAtomicNum() in (7, 8, 16)
                t_electropos = t_atom.GetAtomicNum() in (6, 15, 14)

                if s_electroneg and t_electropos:
                    return (f"a{s_am}", f"a{t_am}")
    except Exception:
        pass

    return None


def _is_radical_reaction(
    reaction_info,
    mol: Optional[Chem.Mol] = None,
    am_to_idx: Optional[dict] = None,
) -> Optional[Tuple[str, str]]:
    """检测单电子转移（自由基反应）。

    自由基反应涉及 1 个电子的转移，反应中心通常存在未成对电子。
    由于自由基反应不需要检查分子结构细节，只需确认电子数 = 1。

    参数
    ----
    reaction_info : tuple
        (source_atom_maps, target_atom_maps, electron_count) 三元组。
    mol : Chem.Mol | None
        RDKit 分子对象（本函数不使用，保留接口一致性）。
    am_to_idx : dict | None
        原子映射字典（本函数不使用，保留接口一致性）。

    返回
    -------
    tuple[str, str] | None
        成功时返回 (自由基引发剂标签, 自由基受体标签)；失败返回 None。
    """
    if not isinstance(reaction_info, (list, tuple)) or len(reaction_info) < 3:
        return None

    source_ams, target_ams, e_count = reaction_info[0], reaction_info[1], reaction_info[2]

    # 自由基反应恰好涉及 1 个电子
    if e_count != 1:
        return None

    if source_ams and target_ams:
        return (f"a{source_ams[0]}", f"a{target_ams[0]}")

    return None


# 向后兼容别名（原始 Django 版代码中的公开接口名）
is_heterolytic_bond_cleavage = _is_heterolytic_bond_cleavage


# ===========================================================================
# 检测方法 → 残基角色标签注册表
# ===========================================================================
# 将每种反应检测函数映射到一组 (角色A标签, 角色B标签)。
# 当检测函数成功识别出反应类型时，可用此注册表查表获取对应的残基角色名称。

METHOD_TO_RESIDUE_FUNCTIONS: Dict[Callable, Tuple[str, str]] = {
    _is_heterolytic_bond_cleavage: ("electrofuge", "nucleofuge"),
    _is_proton_transfer: ("proton donor", "proton acceptor"),
    _is_nucleophilic_attack: ("nucleophile", "electrophile"),
    _is_radical_reaction: ("radical initiator", "radical acceptor"),
}
