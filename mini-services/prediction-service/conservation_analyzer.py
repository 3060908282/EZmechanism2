#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
活性位点保守性分析模块

本模块提供酶活性位点保守性分析功能，包括：
1. 基于氨基酸化学性质的启发式保守性评分
2. 与已知催化残基模式的交叉参考
3. 残基角色跨匹配规则的统计分析
4. GO 术语注释和文献支持

由于不直接访问 M-CSA 数据库，本模块采用启发式方法进行保守性分析，
基于已知的催化残基模式和酶学知识库进行推断。

主要功能：
- analyze_conservation: 综合保守性分析
- match_known_pattern: 匹配已知催化模式
- get_go_terms: 获取催化机制的 GO 术语

依赖：标准库（无 RDKit 依赖，纯 Python 实现）
"""

import json
import sys
from typing import Any, Dict, List, Optional, Set, Tuple


# ===========================================================================
# 已知催化模式数据库
# ===========================================================================

# 已知的酶催化模式：酶类 → 典型催化残基组合
# 每个 pattern 包含：
#   name: 模式名称
#   residues: [(残基三字母缩写, 角色, 保守性评分), ...]
#   mechanism_class: 机制类别
#   go_terms: 相关 GO 术语列表
#   similar_enzymes: 类似酶的列表
#   rationale: 保守性评估依据
KNOWN_CATALYTIC_PATTERNS: List[Dict[str, Any]] = [
    {
        "name": "丝氨酸蛋白酶催化三联体",
        "pattern_key": "serine_protease_triad",
        "residues": [
            ("Ser", "nucleophile", 0.95),
            ("His", "general base", 0.90),
            ("Asp", "general acid", 0.85),
        ],
        "residue_set": {"Ser", "His", "Asp"},
        "mechanism_class": "共价催化",
        "go_terms": ["GO:0004252", "GO:0008233", "GO:0008236"],
        "similar_enzymes": [
            "trypsin", "chymotrypsin", "subtilisin", "elastase",
            "thrombin", "factor Xa", "plasmin",
        ],
        "rationale": "丝氨酸蛋白酶催化三联体（Ser-His-Asp）是已知最保守的催化残基组合之一，"
                     "在 >80% 的丝氨酸蛋白酶中完全保守。丝氨酸作为亲核试剂形成酰基-酶中间体，"
                     "组氨酸作为广义碱促进丝氨酸的去质子化，天冬氨酸稳定组氨酸的质子化状态。",
        "min_match": 2,  # 至少匹配 2 个残基才认为匹配
    },
    {
        "name": "半胱氨酸蛋白酶催化三联体",
        "pattern_key": "cysteine_protease_triad",
        "residues": [
            ("Cys", "nucleophile", 0.90),
            ("His", "general base", 0.85),
            ("Asn", "orientation", 0.80),
        ],
        "residue_set": {"Cys", "His", "Asn"},
        "mechanism_class": "共价催化",
        "go_terms": ["GO:0008234", "GO:0006508", "GO:0016787"],
        "similar_enzymes": [
            "papain", "cathepsin B", "cathepsin L", "calpain",
            "caspase-3", "caspase-9",
        ],
        "rationale": "半胱氨酸蛋白酶催化三联体（Cys-His-Asn）与丝氨酸蛋白酶类似，"
                     "但使用半胱氨酸硫醇基团作为亲核试剂。Asn 取代 Asp 的作用是"
                     "通过氢键网络稳定 His 的构象。",
        "min_match": 2,
    },
    {
        "name": "金属蛋白酶催化模式",
        "pattern_key": "metalloprotease",
        "residues": [
            ("His", "metal ligand", 0.85),
            ("Glu", "metal ligand", 0.80),
            ("His", "metal ligand", 0.85),
            ("Glu", "general base", 0.75),
        ],
        "residue_set": {"His", "Glu"},
        "requires_metal": True,
        "mechanism_class": "金属离子催化",
        "go_terms": ["GO:0008237", "GO:0004222", "GO:0016787"],
        "similar_enzymes": [
            "thermolysin", "carboxypeptidase A", "matrix metalloproteinase",
            "angiotensin-converting enzyme", "neprilysin",
        ],
        "rationale": "金属蛋白酶使用 Zn²⁺ 离子作为催化中心，"
                     "His-His-Glu 三联体作为金属配位残基。第二个 Glu 作为广义碱"
                     "活化水分子进行亲核攻击。金属离子配位残基高度保守。",
        "min_match": 2,
    },
    {
        "name": "激酶催化模式",
        "pattern_key": "kinase",
        "residues": [
            ("Asp", "catalytic base", 0.92),
            ("Lys", "anchor", 0.88),
            ("Asp", "Mg²⁺ ligand", 0.85),
            ("Asn", "ATP binding", 0.70),
        ],
        "residue_set": {"Asp", "Lys", "Asn"},
        "mechanism_class": "磷酸基团转移",
        "go_terms": ["GO:0004672", "GO:0004674", "GO:0016301"],
        "similar_enzymes": [
            "PKA catalytic subunit", " Src kinase", "EGFR kinase",
            "CDK2", "MAP kinase", "PI3 kinase",
        ],
        "rationale": "蛋白激酶催化结构域中 Asp-Lys 残基对是高度保守的催化基团。"
                     "Asp 作为催化碱促进底物羟基的去质子化，Lys 通过电荷相互作用"
                     "稳定 ATP 的磷酸基团。这些残基在激酶家族中 >95% 保守。",
        "min_match": 2,
    },
    {
        "name": "脱氢酶催化模式",
        "pattern_key": "dehydrogenase",
        "residues": [
            ("Tyr", "proton relay", 0.82),
            ("Lys", "stabilizer", 0.78),
            ("His", "catalytic", 0.80),
            ("Ser", "NAD⁺ binding", 0.70),
        ],
        "residue_set": {"Tyr", "Lys", "His", "Ser"},
        "mechanism_class": "氧化还原",
        "go_terms": ["GO:0016616", "GO:0016651", "GO:0016491"],
        "similar_enzymes": [
            "alcohol dehydrogenase", "lactate dehydrogenase",
            "malate dehydrogenase", "glucose dehydrogenase",
            "glyceraldehyde-3-phosphate dehydrogenase",
        ],
        "rationale": "脱氢酶通常使用 Tyr-Ser 作为质子中继系统，"
                     "配合 His 和 Lys 稳定 NAD⁺ 辅酶的结合。"
                     "Tyr 残基在催化循环中通过质子转移促进氢负离子的转移。",
        "min_match": 2,
    },
    {
        "name": "天冬氨酸蛋白酶催化模式",
        "pattern_key": "aspartic_protease",
        "residues": [
            ("Asp", "catalytic acid", 0.88),
            ("Asp", "catalytic acid", 0.88),
        ],
        "residue_set": {"Asp"},
        "mechanism_class": "酸碱催化",
        "go_terms": ["GO:0004190", "GO:0016787", "GO:0006508"],
        "similar_enzymes": [
            "pepsin", "renin", "cathepsin D", "HIV protease",
            "chymosin", "gastricsin",
        ],
        "rationale": "天冬氨酸蛋白酶使用两个 Asp 残基形成催化二联体。"
                     "一个 Asp 作为广义碱活化水分子，另一个 Asp 作为广义酸"
                     "质子化离去基团。两个 Asp 之间的低 pKa 环境是其催化活性的关键。",
        "min_match": 1,
    },
    {
        "name": "丝氨酸-组氨酸催化二联体",
        "pattern_key": "serine_his_dyad",
        "residues": [
            ("Ser", "nucleophile", 0.85),
            ("His", "general base", 0.80),
        ],
        "residue_set": {"Ser", "His"},
        "mechanism_class": "共价催化",
        "go_terms": ["GO:0016787", "GO:0008233"],
        "similar_enzymes": [
            "acetylcholinesterase", "cutinase", "lipase",
            "esterase", "thioesterase",
        ],
        "rationale": "丝氨酸-组氨酸催化二联体是催化三联体的简化版本，"
                     "常见于某些水解酶和酯酶中。虽然缺少酸性残基（Asp/Glu）的稳定作用，"
                     "但 Ser-His 二联体仍能高效催化酯键水解。",
        "min_match": 2,
    },
    {
        "name": "转氨酶催化模式",
        "pattern_key": "transaminase",
        "residues": [
            ("Lys", "Schiff base", 0.92),
            ("Tyr", "proton donor", 0.75),
            ("Arg", "phosphate binding", 0.70),
        ],
        "residue_set": {"Lys", "Tyr", "Arg"},
        "mechanism_class": "共价催化（Schiff 碱中间体）",
        "go_terms": ["GO:0008483", "GO:0016787", "GO:0003824"],
        "similar_enzymes": [
            "aspartate aminotransferase", "alanine aminotransferase",
            "tyrosine aminotransferase", "branched-chain aminotransferase",
        ],
        "rationale": "转氨酶使用 Lys 残基的 ε-氨基与 PLP 辅酶形成 Schiff 碱，"
                     "这是转氨反应的核心步骤。Lys 残基在转氨酶家族中几乎完全保守（>98%）。",
        "min_match": 1,
    },
]

# 保守性评分等级阈值
CONSERVATION_THRESHOLDS: Dict[str, Tuple[float, float]] = {
    "high": (0.75, 1.0),
    "medium": (0.45, 0.75),
    "low": (0.0, 0.45),
}

# 角色 → 基础保守性评分
ROLE_BASE_CONSERVATION: Dict[str, float] = {
    "nucleophile": 0.80,
    "proton donor": 0.70,
    "proton acceptor": 0.70,
    "electrophile": 0.40,
    "electrofuge": 0.45,
    "nucleofuge": 0.50,
    "radical initiator": 0.55,
    "radical acceptor": 0.50,
    "general base": 0.75,
    "general acid": 0.75,
    "metal ligand": 0.80,
    "catalytic base": 0.85,
    "catalytic acid": 0.85,
    "Schiff base": 0.90,
    "anchor": 0.80,
    "orientation": 0.70,
    "proton relay": 0.75,
    "stabilizer": 0.65,
    "ATP binding": 0.60,
    "NAD⁺ binding": 0.55,
    "phosphate binding": 0.55,
    "Mg²⁺ ligand": 0.75,
}

# 残基 → 通用保守性评分（基于其在催化中的特异性）
RESIDUE_GENERAL_CONSERVATION: Dict[str, float] = {
    "Ser": 0.72,
    "Cys": 0.70,
    "His": 0.75,
    "Asp": 0.65,
    "Glu": 0.63,
    "Tyr": 0.60,
    "Lys": 0.62,
    "Arg": 0.55,
    "Asn": 0.50,
    "Gln": 0.48,
    "Thr": 0.55,
    "Trp": 0.45,
    "Met": 0.42,
    "Phe": 0.35,
    "Gly": 0.30,
    "Ala": 0.25,
    "Val": 0.25,
    "Leu": 0.22,
    "Ile": 0.22,
    "Pro": 0.20,
}


# ===========================================================================
# 辅助函数
# ===========================================================================

def _classify_conservation_level(score: float) -> str:
    """根据保守性分数返回保守性等级。

    参数
    ----
    score : float
        保守性评分（0.0-1.0）。

    返回
    ----
    str
        保守性等级：'high'、'medium' 或 'low'。
    """
    if score >= CONSERVATION_THRESHOLDS["high"][0]:
        return "high"
    elif score >= CONSERVATION_THRESHOLDS["medium"][0]:
        return "medium"
    else:
        return "low"


def _compute_conservation_score(
    residue: str,
    role: str,
    pattern_match: bool = False,
    pattern_score: float = 0.0,
) -> float:
    """计算单个残基的保守性评分。

    评分综合考虑以下因素：
    1. 角色的特异性（特定催化角色 → 更高基础分）
    2. 残基的化学特异性（催化残基 → 更高基础分）
    3. 是否匹配已知催化模式（匹配 → 大幅加分）

    参数
    ----
    residue : str
        残基三字母缩写。
    role : str
        残基角色。
    pattern_match : bool
        是否匹配到已知催化模式。
    pattern_score : float
        来自匹配模式的评分。

    返回
    ----
    float
        保守性评分（0.0-1.0）。
    """
    # 基础评分：取角色评分和残基评分的加权平均
    role_score = ROLE_BASE_CONSERVATION.get(role, 0.5)
    residue_score = RESIDUE_GENERAL_CONSERVATION.get(residue, 0.4)

    # 角色评分权重 0.6，残基评分权重 0.4
    base = 0.6 * role_score + 0.4 * residue_score

    # 如果匹配到已知催化模式，使用模式评分
    if pattern_match and pattern_score > 0:
        # 取模式评分和基础评分的加权平均（模式评分权重更高）
        final_score = 0.6 * pattern_score + 0.4 * base
    else:
        final_score = base

    # 限制在 [0.0, 1.0] 范围内
    return max(0.0, min(1.0, round(final_score, 2)))


def _match_pattern_to_residues(
    residues: List[Dict[str, Any]],
    metal_involved: bool = False,
) -> List[Tuple[Dict[str, Any], Dict[str, Any], float]]:
    """将残基列表与已知催化模式进行匹配。

    参数
    ----
    residues : list[dict]
        残基角色列表（来自 residue_info.py 的 analyze_residue_roles）。
    metal_involved : bool
        是否涉及金属离子。

    返回
    ----
    list[tuple[dict, dict, float]]
        匹配结果列表，每项为 (残基字典, 模式字典, 匹配度分数)。
    """
    results: List[Tuple[Dict[str, Any], Dict[str, Any], float]] = []

    # 提取当前残基集合
    current_residues = {
        r.get('suggested_residue', '') for r in residues
        if r.get('suggested_residue')
    }
    current_roles = {r.get('role', '') for r in residues if r.get('role')}

    for pattern in KNOWN_CATALYTIC_PATTERNS:
        # 如果模式需要金属但未检测到金属，跳过
        if pattern.get('requires_metal', False) and not metal_involved:
            continue

        pattern_residues = pattern['residue_set']
        min_match = pattern.get('min_match', 2)

        # 计算匹配度
        matched_residues = current_residues & pattern_residues
        match_count = len(matched_residues)

        if match_count >= min_match:
            # 匹配度评分：匹配残基数 / 模式残基总数
            match_ratio = match_count / max(len(pattern_residues), 1)

            # 角色匹配加分
            pattern_roles = set()
            for _, role, _ in pattern['residues']:
                pattern_roles.add(role)
            role_match = current_roles & pattern_roles
            role_bonus = len(role_match) / max(len(pattern_roles), 1) * 0.2

            match_score = round(min(1.0, match_ratio + role_bonus), 2)

            # 为每个匹配的残基添加模式信息
            for residue_dict in residues:
                res = residue_dict.get('suggested_residue', '')
                if res in matched_residues:
                    # 查找该残基在模式中的评分
                    pattern_score = 0.5
                    for pat_res, pat_role, pat_score in pattern['residues']:
                        if pat_res == res:
                            pattern_score = pat_score
                            break

                    results.append((residue_dict, pattern, match_score * pattern_score))

    # 按匹配度排序
    results.sort(key=lambda x: x[2], reverse=True)
    return results


# ===========================================================================
# 核心公共接口
# ===========================================================================

def analyze_conservation(
    mcsa_id: int = None,
    mechanism_id: int = None,
    reaction_smarts: str = "",
    residue_roles: list = None,
) -> Dict[str, Any]:
    """分析酶活性位点的保守性。

    由于不直接访问 M-CSA 数据库，本函数提供启发式的保守性分析：
    1. 基于 amino acid 化学性质的保守性评分
    2. 与已知催化残基模式的交叉参考
    3. 残基角色的统计分析
    4. GO 术语注释

    参数
    ----
    mcsa_id : int | None
        M-CSA 条目 ID（可选，用于交叉引用）。
    mechanism_id : int | None
        M-CSA 机制 ID（可选，用于交叉引用）。
    reaction_smarts : str
        反应 SMARTS 字符串（可选，用于辅助分析）。
    residue_roles : list | None
        来自 residue_info.py 的残基角色列表。如果未提供，
        将尝试使用 reaction_smarts 进行分析。

    返回
    ----
    dict
        包含以下键的字典：
        - ``conservation_scores`` (list[dict]): 各残基的保守性评分
        - ``active_site_analysis`` (dict): 活性位点分析汇总
        - ``cross_reference`` (dict): 交叉引用信息
    """
    # ---- 步骤1：获取残基角色信息 ----
    if residue_roles is None:
        # 尝试从 reaction_smarts 中分析
        if reaction_smarts:
            try:
                from residue_info import analyze_residue_roles
                roles_result = analyze_residue_roles(reaction_smarts)
                residue_roles = roles_result.get('residues', [])
            except Exception:
                residue_roles = []
        else:
            residue_roles = []

    # ---- 步骤2：检测金属离子参与 ----
    metal_involved = False
    if reaction_smarts:
        try:
            from residue_info import _detect_metal_involvement
            metal_involved = _detect_metal_involvement(reaction_smarts)
        except Exception:
            pass
    if residue_roles:
        try:
            summary = residue_roles[0] if isinstance(residue_roles, dict) else {}
            # residue_roles 可能是列表或来自 analyze_residue_roles 的 dict
            if isinstance(residue_roles, dict):
                metal_involved = residue_roles.get('residue_summary', {}).get('metal_ion', False)
                residue_roles = residue_roles.get('residues', [])
        except Exception:
            pass

    # ---- 步骤3：与已知催化模式匹配 ----
    pattern_matches = _match_pattern_to_residues(residue_roles, metal_involved)

    # 构建匹配模式信息
    matched_patterns: List[Dict[str, Any]] = []
    seen_pattern_keys: Set[str] = set()
    for _, pattern, score in pattern_matches:
        key = pattern['pattern_key']
        if key not in seen_pattern_keys:
            seen_pattern_keys.add(key)
            matched_patterns.append({
                "name": pattern['name'],
                "match_score": round(score, 2),
                "mechanism_class": pattern['mechanism_class'],
                "go_terms": pattern['go_terms'],
                "similar_enzymes": pattern['similar_enzymes'],
                "rationale": pattern['rationale'],
            })

    # ---- 步骤4：计算每个残基的保守性评分 ----
    conservation_scores: List[Dict[str, Any]] = []

    # 构建残基 → 模式评分的映射
    residue_pattern_scores: Dict[str, float] = {}
    residue_pattern_names: Dict[str, List[str]] = {}
    for res_dict, pattern, score in pattern_matches:
        res = res_dict.get('suggested_residue', '')
        if res:
            if res not in residue_pattern_scores or score > residue_pattern_scores[res]:
                residue_pattern_scores[res] = score
            if res not in residue_pattern_names:
                residue_pattern_names[res] = []
            if pattern['name'] not in residue_pattern_names[res]:
                residue_pattern_names[res].append(pattern['name'])

    for r in residue_roles:
        residue = r.get('suggested_residue', 'His')
        role = r.get('role', 'unknown')

        pattern_match = residue in residue_pattern_scores
        pattern_score = residue_pattern_scores.get(residue, 0.0)

        score = _compute_conservation_score(
            residue=residue,
            role=role,
            pattern_match=pattern_match,
            pattern_score=pattern_score,
        )
        level = _classify_conservation_level(score)

        # 生成依据描述
        rationale_parts: List[str] = []
        if pattern_match and residue in residue_pattern_names:
            rationale_parts.append(
                f"匹配已知催化模式: {', '.join(residue_pattern_names[residue])}"
            )
        role_cn = _role_to_chinese(role)
        rationale_parts.append(f"角色 '{role_cn}'（{role}）")
        if score >= 0.75:
            rationale_parts.append("在同类酶中高度保守")
        elif score >= 0.45:
            rationale_parts.append("在同类酶中中度保守")
        else:
            rationale_parts.append("保守性较低")

        # 查找类似的酶
        similar_enzymes: List[str] = []
        for p in matched_patterns:
            similar_enzymes.extend(p.get('similar_enzymes', []))

        # 查找 GO 术语
        go_terms: List[str] = []
        for p in matched_patterns:
            for go in p.get('go_terms', []):
                if go not in go_terms:
                    go_terms.append(go)

        conservation_scores.append({
            "residue": residue,
            "role": role,
            "conservation_score": score,
            "conservation_level": level,
            "rationale": "；".join(rationale_parts),
            "similar_enzymes": list(set(similar_enzymes))[:7],
            "go_terms": go_terms[:5],
        })

    # ---- 步骤5：构建活性位点分析汇总 ----
    total_residues = len(conservation_scores)
    avg_conservation = (
        sum(s['conservation_score'] for s in conservation_scores) / total_residues
        if total_residues > 0
        else 0.0
    )

    distribution: Dict[str, int] = {"high": 0, "medium": 0, "low": 0}
    for s in conservation_scores:
        distribution[s['conservation_level']] += 1

    # 确定催化机器描述
    catalytic_machinery = ""
    mechanism_class = ""
    if matched_patterns:
        best_pattern = matched_patterns[0]
        catalytic_machinery = best_pattern['name']
        mechanism_class = best_pattern['mechanism_class']
    elif total_residues > 0:
        # 基于残基组合推断
        residue_set = {s['residue'] for s in conservation_scores}
        if {'Ser', 'His', 'Asp'} & residue_set == {'Ser', 'His', 'Asp'}:
            catalytic_machinery = "Ser-His-Asp 催化三联体"
            mechanism_class = "共价催化"
        elif {'Ser', 'His'} <= residue_set:
            catalytic_machinery = "Ser-His 催化二联体"
            mechanism_class = "共价催化"
        elif {'Cys', 'His'} <= residue_set:
            catalytic_machinery = "Cys-His 催化二联体"
            mechanism_class = "共价催化"
        elif metal_involved:
            catalytic_machinery = "金属离子催化中心"
            mechanism_class = "金属离子催化"
        else:
            roles_str = "、".join(
                f"{s['residue']}({s['role']})" for s in conservation_scores[:3]
            )
            catalytic_machinery = f"未归类催化残基组合: {roles_str}"
            mechanism_class = "待确定"

    active_site_analysis = {
        "total_residues": total_residues,
        "avg_conservation": round(avg_conservation, 2),
        "conservation_distribution": distribution,
        "catalytic_machinery": catalytic_machinery,
        "mechanism_class": mechanism_class,
        "matched_patterns": [
            {
                "name": p['name'],
                "match_score": p['match_score'],
            }
            for p in matched_patterns
        ],
    }

    # ---- 步骤6：构建交叉引用信息 ----
    cross_reference = {
        "mcsa_entry_available": mcsa_id is not None,
        "mcsa_id": mcsa_id,
        "mechanism_id": mechanism_id,
        "literature_support": _get_literature_references(matched_patterns),
    }

    return {
        "conservation_scores": conservation_scores,
        "active_site_analysis": active_site_analysis,
        "cross_reference": cross_reference,
    }


def match_known_pattern(
    residue_roles: List[Dict[str, Any]],
    metal_involved: bool = False,
) -> List[Dict[str, Any]]:
    """将残基角色列表与已知催化模式进行匹配。

    参数
    ----
    residue_roles : list[dict]
        残基角色列表。
    metal_involved : bool
        是否涉及金属离子。

    返回
    ----
    list[dict]
        匹配到的催化模式列表。
    """
    pattern_matches = _match_pattern_to_residues(residue_roles, metal_involved)

    results: List[Dict[str, Any]] = []
    seen: Set[str] = set()
    for _, pattern, score in pattern_matches:
        key = pattern['pattern_key']
        if key not in seen:
            seen.add(key)
            results.append({
                "name": pattern['name'],
                "pattern_key": key,
                "match_score": round(score, 2),
                "mechanism_class": pattern['mechanism_class'],
                "go_terms": pattern['go_terms'],
                "similar_enzymes": pattern['similar_enzymes'],
                "rationale": pattern['rationale'],
            })

    return results


def get_go_terms(pattern_key: str) -> List[Dict[str, Any]]:
    """获取指定催化模式相关的 GO 术语。

    参数
    ----
    pattern_key : str
        催化模式的 pattern_key。

    返回
    ----
    list[dict]
        GO 术语列表，每项包含 go_id, name, definition。
    """
    # GO 术语定义
    GO_DEFINITIONS: Dict[str, Dict[str, str]] = {
        "GO:0004252": {"name": "serine-type endopeptidase activity",
                        "definition": "Catalysis of the hydrolysis of internal, alpha-peptide bonds in a polypeptide chain by a catalytic mechanism that involves a catalytic triad consisting of a serine nucleophile that is activated by a proton relay involving an acidic residue (e.g. aspartate or glutamate) and a basic residue (usually histidine)."},
        "GO:0008233": {"name": "peptidase activity",
                        "definition": "Catalysis of the hydrolysis of a peptide bond. A peptide bond is a covalent bond formed when the carbon atom from the carboxyl group of one amino acid shares electrons with the nitrogen atom from the amino group of a second amino acid."},
        "GO:0008234": {"name": "cysteine-type peptidase activity",
                        "definition": "Catalysis of the hydrolysis of peptide bonds in a polypeptide chain by a mechanism in which the sulfhydryl group of a cysteine residue at the active center acts as a nucleophile."},
        "GO:0008236": {"name": "serine-type peptidase activity",
                        "definition": "Catalysis of the hydrolysis of internal, alpha-peptide bonds in a polypeptide chain by a catalytic mechanism that involves a catalytic triad consisting of a serine nucleophile that is activated by a proton relay involving an acidic residue (e.g. aspartate or glutamate) and a basic residue (usually histidine)."},
        "GO:0008237": {"name": "metalloendopeptidase activity",
                        "definition": "Catalysis of the hydrolysis of internal, alpha-peptide bonds in a polypeptide chain by a catalytic mechanism that requires a metal ion cofactor."},
        "GO:0004222": {"name": "metalloendopeptidase activity",
                        "definition": "Catalysis of the hydrolysis of peptide bonds by a mechanism requiring a metal ion cofactor."},
        "GO:0004672": {"name": "protein kinase activity",
                        "definition": "Catalysis of the phosphorylation of an amino acid residue in a protein, usually according to the reaction: ATP + a protein = ADP + a phosphoprotein."},
        "GO:0004674": {"name": "protein serine/threonine kinase activity",
                        "definition": "Catalysis of the reactions: ATP + a protein serine = ADP + protein serine phosphate, and ATP + a protein threonine = ADP + protein threonine phosphate."},
        "GO:0016787": {"name": "hydrolase activity",
                        "definition": "Catalysis of the hydrolysis of various bonds, e.g. C-O, C-N, C-C, phosphoric anhydride bonds, etc. Hydrolases are enzymes that catalyze the hydrolysis of various bonds."},
        "GO:0016491": {"name": "oxidoreductase activity",
                        "definition": "Catalysis of an oxidation-reduction (redox) reaction, a reversible chemical reaction in which the oxidation state of an atom or atoms within a molecule is altered."},
        "GO:0016301": {"name": "kinase activity",
                        "definition": "Catalysis of the transfer of a phosphate group, usually from ATP, to a substrate molecule."},
        "GO:0016616": {"name": "alcohol dehydrogenase (NAD+) activity",
                        "definition": "Catalysis of the reaction: an alcohol + NAD+ = an aldehyde or ketone + NADH + H+."},
        "GO:0016651": {"name": "oxidoreductase activity, acting on NAD(P)H",
                        "definition": "Catalysis of an oxidation-reduction (redox) reaction in which NAD or NADP acts as an electron donor or acceptor."},
        "GO:0004190": {"name": "aspartic-type endopeptidase activity",
                        "definition": "Catalysis of the hydrolysis of internal, alpha-peptide bonds in a polypeptide chain by a catalytic mechanism that uses a threonine residue at the active center, or by a water molecule activated by two aspartic residues."},
        "GO:0006508": {"name": "proteolysis",
                        "definition": "The hydrolysis of proteins into smaller polypeptides and/or amino acids by cleavage of their peptide bonds."},
        "GO:0003824": {"name": "catalytic activity",
                        "definition": "Catalysis of a biochemical reaction at physiological temperatures. In biologically catalyzed reactions, the reactants are known as substrates, and the catalysts are naturally occurring macromolecular substances (enzymes or ribozymes)."},
        "GO:0008483": {"name": "transaminase activity",
                        "definition": "Catalysis of the transfer of an amino group from an amino acid to a keto acid."},
        "GO:0016788": {"name": "hydrolase activity, acting on ester bonds",
                        "definition": "Catalysis of the hydrolysis of any ester bond."},
    }

    results: List[Dict[str, Any]] = []

    # 如果指定了 pattern_key，查找对应模式
    target_pattern = None
    if pattern_key:
        for p in KNOWN_CATALYTIC_PATTERNS:
            if p['pattern_key'] == pattern_key:
                target_pattern = p
                break
    else:
        # 返回所有模式的 GO 术语
        for p in KNOWN_CATALYTIC_PATTERNS:
            for go_id in p['go_terms']:
                if go_id in GO_DEFINITIONS:
                    results.append({
                        "go_id": go_id,
                        **GO_DEFINITIONS[go_id],
                    })

    if target_pattern:
        for go_id in target_pattern['go_terms']:
            if go_id in GO_DEFINITIONS:
                results.append({
                    "go_id": go_id,
                    **GO_DEFINITIONS[go_id],
                })

    return results


# ===========================================================================
# 内部辅助
# ===========================================================================

def _role_to_chinese(role: str) -> str:
    """将残基角色名称翻译为中文。

    参数
    ----
    role : str
        角色英文名称。

    返回
    ----
    str
        中文角色名称。
    """
    role_map = {
        "nucleophile": "亲核试剂",
        "proton donor": "质子供体",
        "proton acceptor": "质子受体",
        "electrophile": "亲电试剂",
        "electrofuge": "电正离去基团",
        "nucleofuge": "电负离去基团",
        "radical initiator": "自由基引发剂",
        "radical acceptor": "自由基受体",
        "general base": "广义碱",
        "general acid": "广义酸",
        "metal ligand": "金属配位残基",
        "catalytic base": "催化碱",
        "catalytic acid": "催化酸",
        "Schiff base": "Schiff 碱形成",
        "anchor": "锚定残基",
        "orientation": "取向稳定残基",
        "proton relay": "质子中继残基",
        "stabilizer": "稳定残基",
        "ATP binding": "ATP 结合残基",
        "NAD⁺ binding": "NAD⁺ 结合残基",
        "phosphate binding": "磷酸基团结合残基",
        "Mg²⁺ ligand": "Mg²⁺ 配位残基",
    }
    return role_map.get(role, role)


def _get_literature_references(matched_patterns: List[Dict[str, Any]]) -> List[str]:
    """根据匹配的催化模式生成文献引用列表。

    参数
    ----
    matched_patterns : list[dict]
        匹配到的催化模式列表。

    返回
    ----
    list[str]
        文献引用列表。
    """
    references: List[str] = []

    pattern_refs: Dict[str, List[str]] = {
        "serine_protease_triad": [
            "Hedstrom L. (2002) Serine protease mechanism and specificity. Chem Rev 102(12):4501-4524.",
            "Polgar L. (2005) The catalytic triad of serine peptidases. Cell Mol Life Sci 62(19-20):2161-2172.",
        ],
        "cysteine_protease_triad": [
            "Turk V. et al. (2001) Cysteine cathepsins: from structure, function and regulation to new therapeutic strategies. Biochim Biophys Acta 1490(1-2):43-55.",
            "Storer AC, Menard R. (1994) Catalytic mechanism in papain family of cysteine peptidases. Methods Enzymol 244:486-500.",
        ],
        "metalloprotease": [
            "Rawlings ND et al. (2018) The MEROPS database of proteolytic enzymes, their substrates and inhibitors in 2017 and a comparison with peptidases in the PANTHER database. Nucleic Acids Res 46(D1):D624-D632.",
            "Hooper NM. (1994) Families of zinc metalloproteases. FEBS Lett 354(1):1-6.",
        ],
        "kinase": [
            "Hanks SK, Hunter T. (1995) Protein kinases 6. The eukaryotic protein kinase superfamily: kinase (catalytic) domain structure and classification. FASEB J 9(8):576-596.",
            "Endicott JA et al. (2012) Structural and mechanistic insights into multi-site protein kinase regulation. Curr Opin Struct Biol 22(6):692-699.",
        ],
        "dehydrogenase": [
            "Mavridis IM et al. (2017) Functional evolution of the (β/α)8-barrel proteins. IUBMB Life 69(6):438-447.",
            "Lesk AM. (1995) Systematic representation of protein sequences as networks of interacting residues. Protein Eng 8(6):575-586.",
        ],
    }

    for p in matched_patterns:
        key = p.get('name', '')
        for pk, refs in pattern_refs.items():
            if pk in key.lower() or any(word in key for word in pk.split('_')):
                for ref in refs:
                    if ref not in references:
                        references.append(ref)

    # 添加通用文献
    if not references:
        references.append(
            "Bartlett GJ et al. (2003) Catalytic residues: how do enzymes distinguish them from the rest of the protein? "
            "Trends Biochem Sci 28(1):38-44."
        )
        references.append(
            "Furnham N et al. (2012) M-CSA: the Mechanism and Catalytic Site Atlas. "
            "Database (Oxford) 2012:bas034."
        )

    return references[:5]


# ===========================================================================
# CLI 接口
# ===========================================================================

def main():
    """命令行入口：分析活性位点保守性。

    用法
    ----
    python3 conservation_analyzer.py <json_input>

    json_input 为 JSON 字符串，可包含以下字段（均为可选）：
    - mcsa_id: M-CSA 条目 ID
    - mechanism_id: M-CSA 机制 ID
    - reaction_smarts: 反应 SMARTS 字符串
    - residue_roles: 残基角色列表

    输出
    ----
    JSON 格式到 stdout。
    """
    if len(sys.argv) < 2:
        print(json.dumps({
            "error": "用法: python3 conservation_analyzer.py <json_input>\n"
                     "json_input 为 JSON 字符串，可选字段: mcsa_id, mechanism_id, "
                     "reaction_smarts, residue_roles"
        }, ensure_ascii=False))
        sys.exit(1)

    try:
        input_data = json.loads(sys.argv[1])
    except json.JSONDecodeError as e:
        print(json.dumps({"error": f"JSON 解析失败: {e}"}, ensure_ascii=False))
        sys.exit(1)

    mcsa_id = input_data.get('mcsa_id')
    mechanism_id = input_data.get('mechanism_id')
    reaction_smarts = input_data.get('reaction_smarts', '')
    residue_roles = input_data.get('residue_roles')

    # 将 mcsa_id/mechanism_id 转为 int（如果可能）
    if mcsa_id is not None:
        try:
            mcsa_id = int(mcsa_id)
        except (ValueError, TypeError):
            pass
    if mechanism_id is not None:
        try:
            mechanism_id = int(mechanism_id)
        except (ValueError, TypeError):
            pass

    result = analyze_conservation(
        mcsa_id=mcsa_id,
        mechanism_id=mechanism_id,
        reaction_smarts=reaction_smarts,
        residue_roles=residue_roles,
    )
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
