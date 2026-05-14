# -*- coding: utf-8 -*-
"""
反应类型分析模块 —— 统一的反应类型检测与评分接口。

本模块在机制搜索管线（mechanism search pipeline）中提供反应类型识别能力，
基于 SMARTS 层面的键级变化分析，将规则匹配到的反应步骤分类为以下类型：

- 异裂键断裂（heterolytic_bond_cleavage）
- 质子转移（proton_transfer）
- 亲核攻击（nucleophilic_attack）
- 自由基反应（radical_reaction）
- 一般反应（general_reaction）

模块设计原则
------------
1. **SMARTS 级分析**：由于机制搜索中并不总是拥有原子映射的分子，
   因此本模块在 SMARTS 模式层面进行比较，而不是依赖完整的 RDKit Mol 对象。
2. **独立运行**：仅依赖 RDKit，不依赖 Django 或其他框架。
3. **可组合**：classify_reaction_type / score_reaction_by_type /
   get_bond_change_fingerprint 三个公开接口可独立调用。

函数概览
--------
- get_bond_change_fingerprint : 对比反应物和产物 SMARTS，提取键级变化指纹
- classify_reaction_type      : 基于键变化指纹对反应进行分类
- score_reaction_by_type      : 根据反应类型返回 Dijkstra 分数调节值
"""

import re
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog('rdApp.*')


# ===========================================================================
# 常量
# ===========================================================================

# 亲核原子（高电负性元素）原子序数集合
NUCLEOPHILIC_ELEMENTS: Set[int] = {7, 8, 16}  # N, O, S

# 亲电原子（低电负性目标）原子序数集合
ELECTROPHILIC_ELEMENTS: Set[int] = {6, 14, 15}  # C, Si, P

# 原子序数 → 元素符号映射表（涵盖有机化学常见元素）
_ATOMIC_NUM_TO_SYMBOL: Dict[int, str] = {
    1: 'H', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F',
    14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 35: 'Br', 53: 'I',
}

# 反应类型 → 分数调节值（负值 = 奖励，正值 = 惩罚）
REACTION_TYPE_SCORES: Dict[str, float] = {
    'proton_transfer': -0.1,           # 常见反应，微小奖励
    'nucleophilic_attack': -0.2,       # 特异性反应，较好奖励
    'heterolytic_bond_cleavage': -0.15,  # 特异性反应，中等奖励
    'radical_reaction': 0.2,           # 稀有反应，不稳定性惩罚
    'general_reaction': 0.0,           # 默认，无调节
}

# 反应类型 → 中文描述
REACTION_TYPE_DESCRIPTIONS: Dict[str, str] = {
    'heterolytic_bond_cleavage': '异裂键断裂：共价键的 2-电子断裂',
    'proton_transfer': '质子转移：H⁺ 在两个杂原子之间移动',
    'nucleophilic_attack': '亲核攻击：N/O/S 孤对电子向 C/P/Si 的电子对给予',
    'radical_reaction': '自由基反应：单电子转移过程',
    'general_reaction': '一般反应：未归类为特定类型',
}


# ===========================================================================
# SMARTS 解析辅助
# ===========================================================================

def _parse_smarts_to_bond_set(smarts: str) -> Dict[str, Any]:
    """将 SMARTS 字符串解析为键集合和原子集合。

    从 SMARTS 模式中提取：
    - bonds: 集合 {(原子序数1, 原子序数2, 键级)}，表示识别到的化学键
    - atoms: 集合 {原子序数}，表示出现的元素种类
    - has_hydrogen: 是否涉及氢原子

    由于 SMARTS 模式可能使用原子映射（:N）、电荷（+/-）等标注，
    本函数在解析时会剥离这些标注，仅关注原子元素类型和键拓扑。

    参数
    ----
    smarts : str
        SMARTS 模式字符串（反应物侧或产物侧）。

    返回
    ----
    dict
        包含 'bonds'、'atoms'、'has_hydrogen' 三个键的字典。
    """
    result: Dict[str, Any] = {
        'bonds': set(),
        'atoms': set(),
        'has_hydrogen': False,
    }

    if not smarts:
        return result

    try:
        mol = Chem.MolFromSmarts(smarts)
        if mol is None:
            # 尝试简化后解析
            simplified = _simplify_smarts_for_parsing(smarts)
            if simplified:
                mol = Chem.MolFromSmarts(simplified)
            if mol is None:
                return result
    except Exception:
        return result

    # 提取原子元素信息
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        result['atoms'].add(atomic_num)
        if atomic_num == 1:
            result['has_hydrogen'] = True

    # 提取键信息
    for bond in mol.GetBonds():
        begin_num = bond.GetBeginAtom().GetAtomicNum()
        end_num = bond.GetEndAtom().GetAtomicNum()
        # 用排序后的元组表示无序键对，附加键级
        bond_order = bond.GetBondTypeAsDouble()
        result['bonds'].add((tuple(sorted((begin_num, end_num))), bond_order))

    return result


def _simplify_smarts_for_parsing(smarts: str) -> str:
    """简化 SMARTS 字符串以便解析。

    移除原子映射（:N）、电荷（+/-）、同位素数字等标注，
    但保留原子类型和键拓扑信息。

    参数
    ----
    smarts : str
        原始 SMARTS 字符串。

    返回
    ----
    str
        简化后的 SMARTS 字符串。
    """
    s = re.sub(r':\d+', '', smarts)           # 移除原子映射编号
    s = re.sub(r'\+0', '', s)                 # 移除显式 +0 电荷
    s = re.sub(r'-0', '', s)                  # 移除显式 -0 电荷
    s = re.sub(r'[+-]\d*', '', s)             # 移除电荷标注
    s = re.sub(r'\d+(?=[A-Z#\[\(])', '', s)  # 移除同位素数字
    s = re.sub(r'[Rr]\d*', '', s)             # 移除环成员标注
    s = re.sub(r'\.+', '.', s)                # 合并多个点
    return s.strip('.')


def _parse_smiles_to_bond_set(smiles: str) -> Dict[str, Any]:
    """将 SMILES 字符串解析为键集合和原子集合。

    与 _parse_smarts_to_bond_set 类似，但使用 MolFromSmiles 解析，
    适用于 SMILES 而非 SMARTS 输入。

    参数
    ----
    smiles : str
        SMILES 字符串。

    返回
    ----
    dict
        包含 'bonds'、'atoms'、'has_hydrogen' 三个键的字典。
    """
    result: Dict[str, Any] = {
        'bonds': set(),
        'atoms': set(),
        'has_hydrogen': False,
    }

    if not smiles:
        return result

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result
        mol = Chem.AddHs(mol)  # 添加隐式氢以便检测 H 的变化
    except Exception:
        return result

    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        result['atoms'].add(atomic_num)
        if atomic_num == 1:
            result['has_hydrogen'] = True

    for bond in mol.GetBonds():
        begin_num = bond.GetBeginAtom().GetAtomicNum()
        end_num = bond.GetEndAtom().GetAtomicNum()
        bond_order = bond.GetBondTypeAsDouble()
        result['bonds'].add((tuple(sorted((begin_num, end_num))), bond_order))

    return result


# ===========================================================================
# 核心：键变化指纹提取
# ===========================================================================

def get_bond_change_fingerprint(
    reactant_smarts: str,
    product_smarts: str,
) -> Dict[str, Any]:
    """对比反应物和产物 SMARTS，提取键级变化指纹。

    本函数是反应类型分类的基础。它通过对比反应物和产物的键集合，
    识别出哪些键被断裂（在反应物中存在但在产物中不存在）和
    哪些键被形成（在产物中存在但在反应物中不存在）。

    同时统计电荷变化和元素组成变化。

    参数
    ----
    reactant_smarts : str
        反应物侧的 SMARTS 字符串。
    product_smarts : str
        产物侧的 SMARTS 字符串。

    返回
    ----
    dict
        包含以下键的字典：
        - ``bonds_formed`` : list[tuple[str, str]] — 新形成的键列表，
          每项为 (元素1符号, 元素2符号)。
        - ``bonds_broken`` : list[tuple[str, str]] — 断裂的键列表，
          每项为 (元素1符号, 元素2符号)。
        - ``charges_changed`` : int — 电荷变化的绝对值估计。
        - ``reactant_atoms`` : set[int] — 反应物中的原子序数集合。
        - ``product_atoms`` : set[int] — 产物中的原子序数集合。
        - ``has_hydrogen_in_reactant`` : bool — 反应物是否包含氢。
        - ``has_hydrogen_in_product`` : bool — 产物是否包含氢。
    """
    r_info = _parse_smarts_to_bond_set(reactant_smarts)
    p_info = _parse_smarts_to_bond_set(product_smarts)

    r_bonds: Set[tuple] = r_info['bonds']
    p_bonds: Set[tuple] = p_info['bonds']

    # 键的对称差集：识别变化
    bonds_formed_set = p_bonds - r_bonds  # 产物新增
    bonds_broken_set = r_bonds - p_bonds  # 反应物中消失

    # 将原子序数转换为元素符号
    def _bond_to_elements(bond_tuple: tuple) -> Tuple[str, str]:
        elements = list(bond_tuple[0])
        # 处理同元素键（如 C-C），frozenset 会折叠，此处用元组避免
        if len(elements) == 1:
            sym = _ATOMIC_NUM_TO_SYMBOL.get(elements[0], '?')
            return (sym, sym)
        elem1 = _ATOMIC_NUM_TO_SYMBOL.get(elements[0], '?')
        elem2 = _ATOMIC_NUM_TO_SYMBOL.get(elements[1], '?')
        return (elem1, elem2)

    bonds_formed = sorted(
        [_bond_to_elements(b) for b in bonds_formed_set]
    )
    bonds_broken = sorted(
        [_bond_to_elements(b) for b in bonds_broken_set]
    )

    # 估算电荷变化：通过统计反应物和产物中带电原子标注数
    r_charges = len(re.findall(r'[+-]', reactant_smarts))
    p_charges = len(re.findall(r'[+-]', product_smarts))
    charges_changed = abs(r_charges - p_charges)

    return {
        'bonds_formed': bonds_formed,
        'bonds_broken': bonds_broken,
        'charges_changed': charges_changed,
        'reactant_atoms': r_info['atoms'],
        'product_atoms': p_info['atoms'],
        'has_hydrogen_in_reactant': r_info['has_hydrogen'],
        'has_hydrogen_in_product': p_info['has_hydrogen'],
    }


# ===========================================================================
# 核心：反应类型分类
# ===========================================================================

def classify_reaction_type(
    reaction_smarts: str,
    substrate_smiles: str,
    product_smiles: str,
) -> Dict[str, Any]:
    """对反应步骤进行类型分类。

    综合使用 SMARTS 级键变化指纹分析和 SMILES 级分子结构比较，
    判定反应类型。分类优先级如下：

    1. **异裂键断裂**：反应物和产物之间存在键断裂（非 H 相关）
    2. **质子转移**：键变化涉及氢原子在两个杂原子间移动
    3. **亲核攻击**：新键在高电负性原子（N/O/S）与低电负性原子（C/P/Si）之间形成
    4. **自由基反应**：存在不配对电子变化特征（电荷标记数为奇数等）
    5. **一般反应**：以上均不匹配时的默认分类

    参数
    ----
    reaction_smarts : str
        完整的反应 SMARTS 字符串（格式为 ``reactant>>product``）。
    substrate_smiles : str
        底物分子 SMILES 字符串。
    product_smiles : str
        产物分子 SMILES 字符串。

    返回
    ----
    dict
        包含以下键的字典：
        - ``type`` : str — 反应类型名称。
        - ``confidence`` : float — 分类置信度 (0.0–1.0)。
        - ``details`` : str — 中文描述。
        - ``bond_changes`` : list[dict] — 每个键变化的详细信息列表。
    """
    # 默认结果
    default_result = {
        'type': 'general_reaction',
        'confidence': 0.3,
        'details': REACTION_TYPE_DESCRIPTIONS['general_reaction'],
        'bond_changes': [],
    }

    if not reaction_smarts or '>>' not in reaction_smarts:
        return default_result

    try:
        parts = reaction_smarts.split('>>', 1)
        reactant_smarts = parts[0].strip()
        product_smarts = parts[1].strip()

        if not reactant_smarts or not product_smarts:
            return default_result

        # ---- 步骤1：SMARTS 级键变化指纹 ----
        smarts_fingerprint = get_bond_change_fingerprint(reactant_smarts, product_smarts)

        # ---- 步骤2：SMILES 级键变化指纹（更精确） ----
        smiles_fingerprint = _get_smiles_bond_changes(substrate_smiles, product_smiles)

        # 合并两层分析结果：优先使用 SMILES 级（更精确），回退到 SMARTS 级
        bonds_formed = smiles_fingerprint.get('bonds_formed', []) or smarts_fingerprint.get('bonds_formed', [])
        bonds_broken = smiles_fingerprint.get('bonds_broken', []) or smarts_fingerprint.get('bonds_broken', [])

        # 如果 SMILES 分析失败，回退到 SMARTS 分析
        if not bonds_formed and not bonds_broken:
            bonds_formed = smarts_fingerprint['bonds_formed']
            bonds_broken = smarts_fingerprint['bonds_broken']

        charges_changed = smarts_fingerprint['charges_changed']
        # 氢原子存在性：使用 SMILES 级结果（更可靠），若失败则回退到 SMARTS 级
        has_h_reactant = smiles_fingerprint.get('has_hydrogen_in_reactant', False) or smarts_fingerprint['has_hydrogen_in_reactant']
        has_h_product = smiles_fingerprint.get('has_hydrogen_in_product', False) or smarts_fingerprint['has_hydrogen_in_product']

        # 构建键变化详情列表
        bond_changes: List[Dict[str, str]] = []
        for elem1, elem2 in bonds_formed:
            bond_changes.append({
                'type': 'formed',
                'elements': f"{elem1}-{elem2}",
            })
        for elem1, elem2 in bonds_broken:
            bond_changes.append({
                'type': 'broken',
                'elements': f"{elem1}-{elem2}",
            })

        # ---- 步骤3：反应类型判定 ----
        reaction_type, confidence, details = _classify_from_fingerprint(
            bonds_formed=bonds_formed,
            bonds_broken=bonds_broken,
            charges_changed=charges_changed,
            has_h_reactant=has_h_reactant,
            has_h_product=has_h_product,
            reaction_smarts=reaction_smarts,
        )

        return {
            'type': reaction_type,
            'confidence': round(confidence, 2),
            'details': details,
            'bond_changes': bond_changes,
        }
    except Exception:
        return default_result


def _get_smiles_bond_changes(
    substrate_smiles: str,
    product_smiles: str,
) -> Dict[str, Any]:
    """对比底物和产物 SMILES，提取键级变化指纹。

    参数
    ----
    substrate_smiles : str
        底物 SMILES。
    product_smiles : str
        产物 SMILES。

    返回
    ----
    dict
        与 get_bond_change_fingerprint 返回格式一致的字典。
        解析失败时返回空键集合。
    """
    empty: Dict[str, Any] = {
        'bonds_formed': [],
        'bonds_broken': [],
        'charges_changed': 0,
        'has_hydrogen_in_reactant': False,
        'has_hydrogen_in_product': False,
    }

    if not substrate_smiles or not product_smiles:
        return empty

    try:
        sub_info = _parse_smiles_to_bond_set(substrate_smiles)
        prod_info = _parse_smiles_to_bond_set(product_smiles)

        sub_bonds: Set[tuple] = sub_info['bonds']
        prod_bonds: Set[tuple] = prod_info['bonds']

        bonds_formed_set = prod_bonds - sub_bonds
        bonds_broken_set = sub_bonds - prod_bonds

        def _bond_to_elements_local(bond_tuple: tuple) -> Tuple[str, str]:
            elements = list(bond_tuple[0])
            if len(elements) >= 2:
                return (
                    _ATOMIC_NUM_TO_SYMBOL.get(elements[0], '?'),
                    _ATOMIC_NUM_TO_SYMBOL.get(elements[1], '?'),
                )
            sym = _ATOMIC_NUM_TO_SYMBOL.get(elements[0], '?') if elements else '?'
            return (sym, sym)

        bonds_formed = sorted([_bond_to_elements_local(b) for b in bonds_formed_set])
        bonds_broken = sorted([_bond_to_elements_local(b) for b in bonds_broken_set])

        return {
            'bonds_formed': bonds_formed,
            'bonds_broken': bonds_broken,
            'charges_changed': 0,
            'has_hydrogen_in_reactant': sub_info['has_hydrogen'],
            'has_hydrogen_in_product': prod_info['has_hydrogen'],
        }
    except Exception:
        return empty


def _classify_from_fingerprint(
    bonds_formed: List[Tuple[str, str]],
    bonds_broken: List[Tuple[str, str]],
    charges_changed: int,
    has_h_reactant: bool,
    has_h_product: bool,
    reaction_smarts: str,
) -> Tuple[str, float, str]:
    """根据键变化指纹判定反应类型。

    判定优先级：
    1. 自由基反应（电荷变化为奇数）
    2. 质子转移（H 在两个杂原子间移动）
    3. 亲核攻击（N/O/S → C/P/Si 新键形成）
    4. 异裂键断裂（非 H 键断裂）
    5. 一般反应（默认）

    参数
    ----
    bonds_formed : list[tuple[str, str]]
        新形成的键列表。
    bonds_broken : list[tuple[str, str]]
        断裂的键列表。
    charges_changed : int
        电荷变化数量。
    has_h_reactant : bool
        反应物是否包含氢。
    has_h_product : bool
        产物是否包含氢。
    reaction_smarts : str
        原始反应 SMARTS 字符串。

    返回
    ----
    tuple[str, float, str]
        (反应类型, 置信度, 描述) 三元组。
    """
    # 将元素符号映射回原子序数
    _SYMBOL_TO_NUM: Dict[str, int] = {v: k for k, v in _ATOMIC_NUM_TO_SYMBOL.items()}

    # ---- 检测1：自由基反应 ----
    # 自由基特征：反应物侧和产物侧的电荷标注数均为奇数（不配对电子保留）
    # 注意：单侧电荷变化（如 [O-] → [OH]）是质子转移的特征，不是自由基。
    # 只有两侧电荷数均为奇数时才可能是自由基反应。
    # 在酶学反应中自由基反应极其罕见，因此采用非常保守的判定标准。
    sides = reaction_smarts.split('>>', 1)
    r_charges_count = len(re.findall(r'[+-]', sides[0]))
    p_charges_count = len(re.findall(r'[+-]', sides[1]))
    if (r_charges_count % 2 == 1 and p_charges_count % 2 == 1
            and r_charges_count > 0 and p_charges_count > 0):
        return (
            'radical_reaction',
            0.6,
            REACTION_TYPE_DESCRIPTIONS['radical_reaction'],
        )

    # ---- 检测2：质子转移 ----
    # 特征：键变化涉及 H 原子，且断裂/形成的键都包含 H
    if has_h_reactant and has_h_product:
        h_bonds_formed = [b for b in bonds_formed if 'H' in b]
        h_bonds_broken = [b for b in bonds_broken if 'H' in b]

        # 如果键变化主要是 H 相关的
        if h_bonds_formed and h_bonds_broken:
            # 检查 H 是否从一个杂原子转移到另一个杂原子
            non_h_in_broken = [e for b in h_bonds_broken for e in b if e != 'H']
            non_h_in_formed = [e for b in h_bonds_formed for e in b if e != 'H']
            heteroatoms = {'N', 'O', 'S', 'P'}

            if (any(e in heteroatoms for e in non_h_in_broken) and
                    any(e in heteroatoms for e in non_h_in_formed)):
                confidence = 0.8 if len(h_bonds_formed) <= 2 else 0.6
                return (
                    'proton_transfer',
                    confidence,
                    REACTION_TYPE_DESCRIPTIONS['proton_transfer'],
                )

    # ---- 检测3：亲核攻击 ----
    # 特征：新键在亲核原子（N/O/S）和亲电原子（C/Si/P）之间形成
    for elem1, elem2 in bonds_formed:
        num1 = _SYMBOL_TO_NUM.get(elem1, 0)
        num2 = _SYMBOL_TO_NUM.get(elem2, 0)
        is_nucleophilic = num1 in NUCLEOPHILIC_ELEMENTS
        is_electrophilic = num2 in ELECTROPHILIC_ELEMENTS
        if not is_nucleophilic and not is_electrophilic:
            is_nucleophilic = num2 in NUCLEOPHILIC_ELEMENTS
            is_electrophilic = num1 in ELECTROPHILIC_ELEMENTS
        if is_nucleophilic and is_electrophilic:
            return (
                'nucleophilic_attack',
                0.75,
                REACTION_TYPE_DESCRIPTIONS['nucleophilic_attack'],
            )

    # ---- 检测4：异裂键断裂 ----
    # 特征：存在非 H 键断裂
    non_h_broken = [b for b in bonds_broken if 'H' not in b]
    if non_h_broken:
        confidence = 0.7 if len(non_h_broken) == 1 else 0.5
        return (
            'heterolytic_bond_cleavage',
            confidence,
            REACTION_TYPE_DESCRIPTIONS['heterolytic_bond_cleavage'],
        )

    # ---- 检测5：如果存在 H 键变化但不符合质子转移条件 ----
    h_broken = [b for b in bonds_broken if 'H' in b]
    h_formed = [b for b in bonds_formed if 'H' in b]
    if h_broken or h_formed:
        # H 相关的键变化，可能是质子转移的低置信度变体
        return (
            'proton_transfer',
            0.4,
            REACTION_TYPE_DESCRIPTIONS['proton_transfer'] + '（低置信度）',
        )

    # ---- 默认：一般反应 ----
    if bonds_formed or bonds_broken:
        return (
            'general_reaction',
            0.3,
            REACTION_TYPE_DESCRIPTIONS['general_reaction'],
        )

    return (
        'general_reaction',
        0.1,
        '无法识别的键变化模式',
    )


# ===========================================================================
# 核心：反应类型评分
# ===========================================================================

def score_reaction_by_type(
    reaction_smarts: str,
    substrate_smiles: str,
    product_smiles: str,
) -> float:
    """根据反应类型返回 Dijkstra 分数调节值。

    在机制搜索管线中，每当规则匹配成功并生成产物时，
    调用此函数获取基于反应类型的分数调节值，叠加到 Dijkstra 分数上。

    调节值范围：[-0.3, +0.3]
    - 负值：奖励（更优路径优先探索）
    - 正值：惩罚（降低探索优先级）
    - 零值：无调节

    评分策略
    --------
    - ``proton_transfer``        : -0.1（常见反应，微小奖励）
    - ``nucleophilic_attack``    : -0.2（特异性反应，较好奖励）
    - ``heterolytic_bond_cleavage`` : -0.15（特异性反应，中等奖励）
    - ``radical_reaction``       : +0.2（稀有反应，不稳定性惩罚）
    - ``general_reaction``       : 0.0（默认，无调节）

    参数
    ----
    reaction_smarts : str
        完整的反应 SMARTS 字符串（格式为 ``reactant>>product``）。
    substrate_smiles : str
        底物分子 SMILES 字符串。
    product_smiles : str
        产物分子 SMILES 字符串。

    返回
    ----
    float
        分数调节值，范围 [-0.3, +0.3]。
    """
    try:
        result = classify_reaction_type(reaction_smarts, substrate_smiles, product_smiles)
        rtype = result.get('type', 'general_reaction')
        base_score = REACTION_TYPE_SCORES.get(rtype, 0.0)

        # 置信度加权：低置信度时减弱分数调节
        confidence = result.get('confidence', 0.5)
        weighted_score = base_score * confidence

        # 限制在 [-0.3, +0.3] 范围内
        return max(-0.3, min(0.3, round(weighted_score, 3)))
    except Exception:
        return 0.0
