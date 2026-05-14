#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
氨基酸残基角色识别与注释模块

本模块用于分析酶催化反应步骤中涉及的氨基酸残基角色，
包括质子供体/受体、亲核试剂/亲电试剂、自由基引发剂/受体等角色的识别。

主要功能：
- 从反应 SMARTS 中解析原子映射信息，构建合成反应信息元组
- 调用 METHOD_TO_RESIDUE_FUNCTIONS 中的检测函数进行残基角色识别
- 基于化学原理将角色映射到可能的氨基酸残基
- 检测催化三联体（Ser-His-Asp）和催化二联体（Ser-His）
- 检测金属离子参与

依赖：rdkit, common.others, reaction_type_analyzer
"""

import json
import re
import sys
from typing import Any, Dict, List, Optional, Tuple, Set

from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')

from common.others import (
    METHOD_TO_RESIDUE_FUNCTIONS,
    _is_proton_transfer,
    _is_nucleophilic_attack,
    _is_heterolytic_bond_cleavage,
    _is_radical_reaction,
    AA_1TO3,
    AA_3TO1,
)
from reaction_type_analyzer import classify_reaction_type


# ===========================================================================
# 常量定义
# ===========================================================================

# 常见金属元素符号（用于检测金属离子参与）
METAL_ELEMENTS: Set[str] = {
    'Fe', 'Zn', 'Mg', 'Mn', 'Cu', 'Ca',
    'Co', 'Ni', 'Mo', 'W', 'Se', 'V', 'Cr',
}

# 金属元素 → 原子序数映射（用于 SMARTS 中的原子检测）
_METAL_ATOMIC_NUMS: Set[int] = {
    26, 30, 12, 25, 29, 20,  # Fe, Zn, Mg, Mn, Cu, Ca
    27, 28, 42, 74, 34, 23, 24,  # Co, Ni, Mo, W, Se, V, Cr
}

# 角色 → 最可能的氨基酸残基列表
# 基于酶学中的常见催化残基化学性质
ROLE_TO_RESIDUES: Dict[str, List[str]] = {
    "proton donor": ["Ser", "His", "Tyr", "Cys", "Lys", "Asp", "Glu"],
    "proton acceptor": ["His", "Asp", "Glu", "Ser", "Tyr", "Cys", "Lys"],
    "nucleophile": ["Ser", "Cys", "Thr", "Tyr", "Lys", "Asp", "Glu", "His"],
    "electrophile": ["Gly", "Ala", "Val", "Leu"],  # 亲电中心通常在底物上
    "electrofuge": ["Gly", "Ala", "Val"],  # 离去基团通常在底物上
    "nucleofuge": ["Ser", "Cys", "Tyr", "His", "Asp", "Glu"],
    "radical initiator": ["Tyr", "Trp", "Cys", "Met", "His"],
    "radical acceptor": ["Phe", "Tyr", "Trp", "His", "Cys", "Met"],
}

# 角色 → 反转角色映射
ROLE_TO_COUNTER: Dict[str, str] = {
    "proton donor": "proton acceptor",
    "proton acceptor": "proton donor",
    "nucleophile": "electrophile",
    "electrophile": "nucleophile",
    "electrofuge": "nucleofuge",
    "nucleofuge": "electrofuge",
    "radical initiator": "radical acceptor",
    "radical acceptor": "radical initiator",
}

# 角色 → 检测方法名称映射
ROLE_TO_DETECTION_METHOD: Dict[str, str] = {
    "proton donor": "proton_transfer",
    "proton acceptor": "proton_transfer",
    "nucleophile": "nucleophilic_attack",
    "electrophile": "nucleophilic_attack",
    "electrofuge": "heterolytic_bond_cleavage",
    "nucleofuge": "heterolytic_bond_cleavage",
    "radical initiator": "radical_reaction",
    "radical acceptor": "radical_reaction",
}

# 角色 → 置信度（基于检测方法的可靠性）
ROLE_CONFIDENCE: Dict[str, float] = {
    "proton donor": 0.75,
    "proton acceptor": 0.75,
    "nucleophile": 0.80,
    "electrophile": 0.70,
    "electrofuge": 0.65,
    "nucleofuge": 0.65,
    "radical initiator": 0.60,
    "radical acceptor": 0.60,
}


# ===========================================================================
# 辅助函数：从反应 SMARTS 解析合成反应信息
# ===========================================================================

def _parse_atom_mappings_from_smarts(smarts: str) -> Dict[int, int]:
    """从 SMARTS 字符串中提取原子映射编号和对应的原子序数。

    通过解析 ``[原子符号:映射编号]`` 模式，构建
    映射编号 → 原子序数 的字典。

    参数
    ----
    smarts : str
        SMARTS 模式字符串（反应物侧或产物侧）。

    返回
    ----
    dict[int, int]
        原子映射编号 → 原子序数的字典。
    """
    am_to_atomic_num: Dict[int, int] = {}

    # 匹配 [元素符号:数字] 或 [#原子序数:数字] 模式
    # 例如 [CH3:1], [OH:2], [#6:3], [O-:4]
    pattern = re.compile(r'\[([^\]:]+):(\d+)\]')
    element_map = {
        'C': 6, 'N': 7, 'O': 8, 'S': 16, 'P': 15, 'H': 1,
        'F': 9, 'Cl': 17, 'Br': 35, 'I': 53, 'Se': 34,
        'Si': 14, 'B': 5, 'Na': 11, 'K': 19, 'Mg': 12,
        'Ca': 20, 'Fe': 26, 'Zn': 30, 'Mn': 25, 'Cu': 29,
        'Co': 27, 'Ni': 28, 'Mo': 42, 'W': 74, 'V': 23, 'Cr': 24,
    }

    for match in pattern.finditer(smarts):
        atom_desc = match.group(1)
        map_num = int(match.group(2))

        # 处理 [#原子序数:数字] 格式
        if atom_desc.startswith('#'):
            try:
                atomic_num = int(atom_desc[1:].split(',')[0].split(']')[0])
                am_to_atomic_num[map_num] = atomic_num
            except (ValueError, IndexError):
                pass
            continue

        # 提取元素符号（可能包含电荷、氢计数等标注）
        # 例如 "CH3", "OH-", "NH3+", "C"
        element = re.match(r'^([A-Z][a-z]?)', atom_desc)
        if element:
            symbol = element.group(1)
            if symbol in element_map:
                am_to_atomic_num[map_num] = element_map[symbol]
            else:
                # 尝试从 RDKit 获取原子序数
                try:
                    mol = Chem.MolFromSmarts(f'[{symbol}]')
                    if mol and mol.GetNumAtoms() == 1:
                        am_to_atomic_num[map_num] = mol.GetAtomWithIdx(0).GetAtomicNum()
                except Exception:
                    pass

    return am_to_atomic_num


def _build_synthetic_reaction_info(
    reaction_smarts: str,
) -> Tuple[Optional[Tuple[List[int], List[int], int]], Optional[Chem.Mol], Optional[Dict[int, int]]]:
    """从反应 SMARTS 构建合成的 reaction_info 元组和分子对象。

    由于本模块没有完整的弯箭头机制图（来自 M-CSA XML），
    因此从反应 SMARTS 中推断可能的电子转移事件。

    构建策略：
    1. 解析反应物和产物两侧的原子映射
    2. 识别键的变化（断裂、形成、键级变化）
    3. 为每个变化生成一个合成的 (source_atom_maps, target_atom_maps, electron_count) 元组
    4. 同时尝试解析反应物侧为 RDKit Mol 对象

    参数
    ----
    reaction_smarts : str
        完整的反应 SMARTS 字符串（格式为 ``reactant>>product``）。

    返回
    ----
    tuple
        (reaction_info, mol, am_to_idx)
        - reaction_info: 合成的反应信息元组，或 None
        - mol: 反应物侧的 RDKit Mol 对象，或 None
        - am_to_idx: 原子映射编号 → 原子索引的字典，或 None
    """
    if not reaction_smarts or '>>' not in reaction_smarts:
        return None, None, None

    try:
        parts = reaction_smarts.split('>>', 1)
        reactant_smarts = parts[0].strip()
        product_smarts = parts[1].strip()

        if not reactant_smarts or not product_smarts:
            return None, None, None

        # 解析两侧原子映射
        r_am_map = _parse_atom_mappings_from_smarts(reactant_smarts)
        p_am_map = _parse_atom_mappings_from_smarts(product_smarts)

        if not r_am_map and not p_am_map:
            return None, None, None

        # 尝试解析反应物为 Mol 对象
        mol = None
        am_to_idx: Optional[Dict[int, int]] = None
        try:
            mol = Chem.MolFromSmarts(reactant_smarts)
            if mol is not None:
                am_to_idx = {}
                for atom in mol.GetAtoms():
                    amn = atom.GetAtomMapNum()
                    if amn != 0:
                        am_to_idx[amn] = atom.GetIdx()
        except Exception:
            pass

        # 识别键变化
        # 解析反应物和产物的键集合（使用原子映射编号）
        r_bonds: Set[Tuple[int, int]] = set()
        p_bonds: Set[Tuple[int, int]] = set()

        if mol is not None:
            for bond in mol.GetBonds():
                am1 = bond.GetBeginAtom().GetAtomMapNum()
                am2 = bond.GetEndAtom().GetAtomMapNum()
                if am1 != 0 and am2 != 0:
                    r_bonds.add((min(am1, am2), max(am1, am2)))

        try:
            prod_mol = Chem.MolFromSmarts(product_smarts)
            if prod_mol is not None:
                for bond in prod_mol.GetBonds():
                    am1 = bond.GetBeginAtom().GetAtomMapNum()
                    am2 = bond.GetEndAtom().GetAtomMapNum()
                    if am1 != 0 and am2 != 0:
                        p_bonds.add((min(am1, am2), max(am1, am2)))
        except Exception:
            pass

        # 找出变化的键对，构建合成 reaction_info
        # 键断裂：反应物中有但产物中没有
        broken_bonds = r_bonds - p_bonds
        # 键形成：产物中有但反应物中没有
        formed_bonds = p_bonds - r_bonds

        # 为第一个键变化构建反应信息
        # 默认使用 2 电子过程
        electron_count = 2
        source_ams: List[int] = []
        target_ams: List[int] = []

        if broken_bonds:
            # 取第一个断裂的键
            a, b = list(broken_bonds)[0]
            source_ams = [a, b]
            target_ams = [a]  # 电子留在其中一个原子上
        elif formed_bonds:
            # 取第一个形成的键
            a, b = list(formed_bonds)[0]
            source_ams = [a]
            target_ams = [b]

        if source_ams and target_ams:
            reaction_info = (source_ams, target_ams, electron_count)
            return reaction_info, mol, am_to_idx

        # 如果没有键变化，但两侧都有原子映射，尝试基于电荷变化推断
        # 检查是否有电荷标注
        r_charges = len(re.findall(r'[+-]', reactant_smarts))
        p_charges = len(re.findall(r'[+-]', product_smarts))

        if r_am_map and p_am_map and r_charges != p_charges:
            # 可能是质子转移或自由基反应
            # 使用反应物和产物中共享的原子映射编号
            common_ams = set(r_am_map.keys()) & set(p_am_map.keys())
            if common_ams:
                am_list = sorted(common_ams)
                if len(am_list) >= 2:
                    reaction_info = ([am_list[0]], [am_list[1]], 2)
                    return reaction_info, mol, am_to_idx

        return None, mol, am_to_idx

    except Exception:
        return None, None, None


def _detect_metal_involvement(smarts: str) -> bool:
    """检测反应 SMARTS 中是否涉及金属元素。

    参数
    ----
    smarts : str
        反应 SMARTS 字符串。

    返回
    ----
    bool
        是否检测到金属元素参与。
    """
    # 直接检查金属元素符号
    for metal in METAL_ELEMENTS:
        # 使用正则匹配元素符号，避免匹配到含有金属符号的更长字符串
        if re.search(r'\[' + metal + r'[^a-z]', smarts):
            return True
        if re.search(r'\[#\d+' + r'\]' , smarts):
            # 检查原子序数是否为金属
            nums = re.findall(r'#(\d+)', smarts)
            for num_str in nums:
                try:
                    if int(num_str) in _METAL_ATOMIC_NUMS:
                        return True
                except ValueError:
                    pass

    # 尝试通过 RDKit 解析
    try:
        if '>>' in smarts:
            sides = smarts.split('>>', 1)
            for side in sides:
                mol = Chem.MolFromSmarts(side.strip())
                if mol is not None:
                    for atom in mol.GetAtoms():
                        if atom.GetAtomicNum() in _METAL_ATOMIC_NUMS:
                            return True
    except Exception:
        pass

    return False


def _check_catalytic_triad(residues: List[Dict[str, Any]]) -> bool:
    """检测是否形成催化三联体（Ser-His-Asp/Glu）。

    催化三联体是丝氨酸蛋白酶等酶类中常见的催化残基组合，
    由丝氨酸（亲核试剂）、组氨酸（碱）和天冬氨酸/谷氨酸（酸）组成。

    参数
    ----
    residues : list[dict]
        已识别的残基角色列表。

    返回
    ----
    bool
        是否检测到催化三联体模式。
    """
    suggested = {r.get('suggested_residue', '') for r in residues if r.get('suggested_residue')}
    roles = {r.get('role', '') for r in residues}

    has_ser = 'Ser' in suggested
    has_his = 'His' in suggested
    has_acid = 'Asp' in suggested or 'Glu' in suggested

    # 需要同时存在亲核角色和碱/酸角色
    has_nucleophile = 'nucleophile' in roles
    has_proton_transfer = 'proton donor' in roles or 'proton acceptor' in roles

    return has_ser and has_his and has_acid and has_nucleophile and has_proton_transfer


def _check_catalytic_dyad(residues: List[Dict[str, Any]]) -> bool:
    """检测是否形成催化二联体（Ser-His）。

    催化二联体是催化三联体的简化版本，存在于某些丝氨酸蛋白酶中。

    参数
    ----
    residues : list[dict]
        已识别的残基角色列表。

    返回
    ----
    bool
        是否检测到催化二联体模式。
    """
    suggested = {r.get('suggested_residue', '') for r in residues if r.get('suggested_residue')}

    has_ser = 'Ser' in suggested
    has_his = 'His' in suggested

    return has_ser and has_his


def _build_residue_description(
    role: str,
    suggested_residue: str,
    position: Optional[int],
    detection_method: str,
) -> str:
    """构建残基角色描述文本。

    参数
    ----
    role : str
        残基角色名称。
    suggested_residue : str
        建议的氨基酸残基三字母缩写。
    position : int | None
        残基位置编号。
    detection_method : str
        检测方法名称。

    返回
    ----
    str
        中文描述文本。
    """
    # 三字母 → 中文氨基酸名
    aa_names = {
        'Ser': '丝氨酸', 'His': '组氨酸', 'Asp': '天冬氨酸', 'Glu': '谷氨酸',
        'Cys': '半胱氨酸', 'Tyr': '酪氨酸', 'Lys': '赖氨酸', 'Thr': '苏氨酸',
        'Trp': '色氨酸', 'Met': '甲硫氨酸', 'Phe': '苯丙氨酸', 'Ala': '丙氨酸',
        'Gly': '甘氨酸', 'Val': '缬氨酸', 'Leu': '亮氨酸', 'Asn': '天冬酰胺',
        'Gln': '谷氨酰胺', 'Arg': '精氨酸', 'Pro': '脯氨酸', 'Ile': '异亮氨酸',
    }

    # 角色名 → 中文描述
    role_names = {
        'proton donor': '质子供体',
        'proton acceptor': '质子受体',
        'nucleophile': '亲核试剂',
        'electrophile': '亲电试剂',
        'electrofuge': '电正离去基团',
        'nucleofuge': '电负离去基团',
        'radical initiator': '自由基引发剂',
        'radical acceptor': '自由基受体',
    }

    aa_cn = aa_names.get(suggested_residue, suggested_residue)
    role_cn = role_names.get(role, role)

    pos_str = f'-{position}' if position else ''
    method_map = {
        'proton_transfer': '质子转移检测',
        'nucleophilic_attack': '亲核攻击检测',
        'heterolytic_bond_cleavage': '异裂键断裂检测',
        'radical_reaction': '自由基反应检测',
        'reaction_type_classification': '反应类型分类',
    }
    method_cn = method_map.get(detection_method, detection_method)

    return (
        f"{suggested_residue}{pos_str}（{aa_cn}）可能作为{role_cn}参与催化反应。"
        f"该推断基于{method_cn}方法。"
    )


# ===========================================================================
# 核心公共接口
# ===========================================================================

def analyze_residue_roles(
    reaction_smarts: str,
    substrate_smiles: str = "",
    product_smiles: str = "",
) -> Dict[str, Any]:
    """分析反应步骤中涉及的催化残基角色。

    通过以下步骤进行分析：
    1. 从反应 SMARTS 中解析原子映射信息
    2. 构建合成的 reaction_info 元组
    3. 调用 METHOD_TO_RESIDUE_FUNCTIONS 中的检测函数
    4. 使用反应类型分类器进行辅助判断
    5. 将检测到的角色映射到最可能的氨基酸残基
    6. 检测催化三联体/二联体和金属离子参与

    参数
    ----
    reaction_smarts : str
        反应 SMARTS 字符串（格式为 ``reactant>>product``）。
    substrate_smiles : str
        底物分子 SMILES 字符串（可选，用于辅助分析）。
    product_smiles : str
        产物分子 SMILES 字符串（可选，用于辅助分析）。

    返回
    ----
    dict
        包含以下键的字典：
        - ``residues`` (list[dict]): 识别的残基角色列表
        - ``residue_summary`` (dict): 残基角色汇总信息
    """
    result_residues: List[Dict[str, Any]] = []
    seen_roles: Set[str] = set()

    # ---- 步骤1：使用反应类型分类器获取初步判断 ----
    reaction_type_result = classify_reaction_type(
        reaction_smarts, substrate_smiles, product_smiles
    )
    rtype = reaction_type_result.get('type', 'general_reaction')
    rtype_confidence = reaction_type_result.get('confidence', 0.0)

    # 基于反应类型推断残基角色
    rtype_role_map = {
        'proton_transfer': [
            ('proton donor', 'proton acceptor'),
        ],
        'nucleophilic_attack': [
            ('nucleophile', 'electrophile'),
        ],
        'heterolytic_bond_cleavage': [
            ('electrofuge', 'nucleofuge'),
        ],
        'radical_reaction': [
            ('radical initiator', 'radical acceptor'),
        ],
    }

    if rtype in rtype_role_map and rtype_confidence >= 0.3:
        for role_a, role_b in rtype_role_map[rtype]:
            for role in [role_a, role_b]:
                if role not in seen_roles:
                    suggested_residues = ROLE_TO_RESIDUES.get(role, ['His', 'Ser'])
                    counter_role = ROLE_TO_COUNTER.get(role, '')
                    confidence = round(ROLE_CONFIDENCE.get(role, 0.5) * rtype_confidence, 2)

                    result_residues.append({
                        "role": role,
                        "counter_role": counter_role,
                        "confidence": confidence,
                        "detection_method": "reaction_type_classification",
                        "description": _build_residue_description(
                            role, suggested_residues[0], None, "reaction_type_classification"
                        ),
                        "suggested_residue": suggested_residues[0],
                        "suggested_position": None,
                    })
                    seen_roles.add(role)

    # ---- 步骤2：使用 METHOD_TO_RESIDUE_FUNCTIONS 检测函数进行精细检测 ----
    reaction_info, mol, am_to_idx = _build_synthetic_reaction_info(reaction_smarts)

    if reaction_info is not None and mol is not None:
        for detect_func, (role_a, role_b) in METHOD_TO_RESIDUE_FUNCTIONS.items():
            try:
                detection_result = detect_func(reaction_info, mol=mol, am_to_idx=am_to_idx)
                if detection_result is not None:
                    # 检测成功，提取角色
                    for role in [role_a, role_b]:
                        if role not in seen_roles:
                            suggested_residues = ROLE_TO_RESIDUES.get(role, ['His'])
                            counter_role = ROLE_TO_COUNTER.get(role, '')
                            base_confidence = ROLE_CONFIDENCE.get(role, 0.6)
                            # 如果反应类型分类器也支持此角色，提高置信度
                            if rtype in rtype_role_map:
                                for ra, rb in rtype_role_map[rtype]:
                                    if role == ra or role == rb:
                                        base_confidence = min(1.0, base_confidence + 0.15)
                                        break

                            detection_method_name = ROLE_TO_DETECTION_METHOD.get(
                                role, detect_func.__name__
                            )

                            result_residues.append({
                                "role": role,
                                "counter_role": counter_role,
                                "confidence": round(base_confidence, 2),
                                "detection_method": detection_method_name,
                                "description": _build_residue_description(
                                    role, suggested_residues[0], None, detection_method_name
                                ),
                                "suggested_residue": suggested_residues[0],
                                "suggested_position": None,
                            })
                            seen_roles.add(role)
            except Exception:
                # 检测函数调用失败，跳过
                pass

    # ---- 步骤3：如果没有任何检测结果，生成默认分析 ----
    if not result_residues and reaction_smarts:
        # 根据反应 SMARTS 中的元素组成推断可能的角色
        has_nitrogen = 'N' in reaction_smarts or '[#7' in reaction_smarts
        has_oxygen = 'O' in reaction_smarts or '[#8' in reaction_smarts
        has_sulfur = 'S' in reaction_smarts or '[#16' in reaction_smarts

        if has_nitrogen or has_oxygen:
            result_residues.append({
                "role": "proton donor",
                "counter_role": "proton acceptor",
                "confidence": 0.3,
                "detection_method": "element_analysis",
                "description": (
                    "基于反应物中含氮/氧杂原子的分析，"
                    "反应可能涉及质子转移步骤。"
                ),
                "suggested_residue": "His" if has_nitrogen else "Ser",
                "suggested_position": None,
            })
            if has_nitrogen:
                result_residues.append({
                    "role": "nucleophile",
                    "counter_role": "electrophile",
                    "confidence": 0.25,
                    "detection_method": "element_analysis",
                    "description": (
                        "基于反应物中含氮原子的分析，"
                        "可能存在含氮亲核试剂参与。"
                    ),
                    "suggested_residue": "His",
                    "suggested_position": None,
                })
        if has_sulfur:
            result_residues.append({
                "role": "nucleophile",
                "counter_role": "electrophile",
                "confidence": 0.3,
                "detection_method": "element_analysis",
                "description": (
                    "基于反应物中含硫原子的分析，"
                    "半胱氨酸可能作为亲核试剂参与反应。"
                ),
                "suggested_residue": "Cys",
                "suggested_position": None,
            })

    # ---- 步骤4：构建残基汇总信息 ----
    roles_breakdown: Dict[str, int] = {}
    for r in result_residues:
        role = r.get('role', 'unknown')
        roles_breakdown[role] = roles_breakdown.get(role, 0) + 1

    metal_involved = _detect_metal_involvement(reaction_smarts) if reaction_smarts else False
    catalytic_triad = _check_catalytic_triad(result_residues)
    catalytic_dyad = _check_catalytic_dyad(result_residues) and not catalytic_triad

    residue_summary = {
        "total_residues": len(result_residues),
        "roles_breakdown": roles_breakdown,
        "catalytic_triad": catalytic_triad,
        "catalytic_dyad": catalytic_dyad,
        "metal_ion": metal_involved,
    }

    return {
        "residues": result_residues,
        "residue_summary": residue_summary,
    }


# ===========================================================================
# CLI 接口
# ===========================================================================

def main():
    """命令行入口：分析反应步骤中的残基角色。

    用法
    ----
    python3 residue_info.py <reaction_smarts> [substrate_smiles] [product_smiles]

    输出
    ----
    JSON 格式到 stdout。
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "error": "用法: python3 residue_info.py <reaction_smarts> [substrate_smiles] [product_smiles]"
        }, ensure_ascii=False))
        sys.exit(1)

    reaction_smarts = sys.argv[1]
    substrate_smiles = sys.argv[2] if len(sys.argv) > 2 else ""
    product_smiles = sys.argv[3] if len(sys.argv) > 3 else ""

    result = analyze_residue_roles(reaction_smarts, substrate_smiles, product_smiles)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
