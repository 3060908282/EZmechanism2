#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
外部数据库交叉引用链接生成模块

本模块为匹配到的反应规则和分子生成外部数据库的超链接，
方便研究人员在 M-CSA、PDBe、PDB、UniProt、KEGG 等数据库中
进一步查看详细的酶学机制、蛋白质结构和分子信息。

主要功能：
- get_rule_links: 为反应规则生成外部数据库链接
- get_molecule_links: 为分子 SMILES 生成化学数据库搜索链接
- get_ec_links: 为 EC 编号生成 BRENDA/Expasy 链接
- get_all_links: 综合生成规则和分子的所有外部链接

数据来源参考论文中的 M-CSA 数据库（https://www.ebi.ac.uk/mcsa/）。

依赖：标准库（urllib.parse）
"""

import json
import re
import sys
from typing import Any, Dict, List, Optional
from urllib.parse import quote, urlencode


# ===========================================================================
# 数据库 URL 常量
# ===========================================================================

# 基础 URL 模板
MCSA_BASE = "https://www.ebi.ac.uk/mcsa"
PDBE_BASE = "https://www.ebi.ac.uk/pdbe/entry/pdb"
CATH_BASE = "http://www.cathdb.info/domain"
PDB_BASE = "https://www.rcsb.org/structure"
PDBSUM_BASE = "https://www.ebi.ac.uk/pdbsum"
UNIPROT_BASE = "https://www.uniprot.org/uniprot"
KEGG_COMPOUND_BASE = "https://www.genome.jp/dbget-bin/www_bget"
INTERPRO_BASE = "https://www.ebi.ac.uk/interpro/entry/InterPro"
FUNTREE_BASE = "https://www.ebi.ac.uk/funtree"
MACIE_BASE = "https://www.ebi.ac.uk/thornton-srv/databases/MACiE"
BRENDA_BASE = "https://www.brenda-enzymes.org"
EXPASY_BASE = "https://enzyme.expasy.org"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov"
CHEMSPIDER_BASE = "https://www.chemspider.com"
CHEBI_BASE = "https://www.ebi.ac.uk/chebi"

# 数据库名称映射（中文）
DB_LABELS: Dict[str, str] = {
    "mcsa": "M-CSA",
    "pdbe": "PDBe",
    "cath": "CATH",
    "pdb": "PDB",
    "pdbsum": "PDBsum",
    "uniprot": "UniProt",
    "kegg": "KEGG",
    "integrity": "数据完整性",
    "pubchem": "PubChem",
    "chemspider": "ChemSpider",
    "chebi_search": "ChEBI 搜索",
    "brenda": "BRENDA",
    "expasy": "ExPASy",
    "interpro": "InterPro",
    "funtree": "FunTree",
    "macie": "MACiE",
}


# ===========================================================================
# 辅助函数
# ===========================================================================

def _extract_pdb_id(mcsa_id: Any) -> Optional[str]:
    """从 M-CSA 数据中提取 PDB ID。

    M-CSA ID 通常为整数，其对应的 PDB ID 需要通过数据库查询获取。
    此处提供启发式提取：如果 mcsa_id 包含 PDB ID 格式（4 位字母数字），
    则直接提取。

    参数
    ----
    mcsa_id : any
        M-CSA 标识符（整数或字符串）。

    返回
    ----
    str | None
        提取的 PDB ID，或 None。
    """
    if mcsa_id is None:
        return None

    s = str(mcsa_id)

    # 如果本身就是 4 位字母数字的 PDB ID 格式
    if re.match(r'^[0-9][a-zA-Z0-9]{3}$', s):
        return s.lower()

    # 尝试从字符串中查找 PDB ID 模式
    match = re.search(r'\b([0-9][a-zA-Z0-9]{3})\b', s)
    if match:
        return match.group(1).lower()

    return None


def _extract_uniprot_id(rule: Dict[str, Any]) -> Optional[str]:
    """从规则字典中尝试提取 UniProt ID。

    参数
    ----
    rule : dict
        规则字典。

    返回
    ----
    str | None
        UniProt ID，或 None。
    """
    # 常见 UniProt ID 格式：以字母开头，后跟 5-10 位字母数字
    for key in ['uniprot_id', 'uniprot', 'accession']:
        if key in rule and rule[key]:
            val = str(rule[key]).strip()
            if re.match(r'^[A-Z][A-Z0-9]{5,10}$', val):
                return val

    return None


def _extract_ec_number(rule: Dict[str, Any]) -> Optional[str]:
    """从规则字典中提取 EC 编号。

    参数
    ----
    rule : dict
        规则字典。

    返回
    ----
    str | None
        EC 编号字符串，或 None。
    """
    for key in ['ec_number', 'ec', 'ec_no']:
        if key in rule and rule[key]:
            val = str(rule[key]).strip()
            if val and val != 'None' and re.match(r'EC\s*', val, re.IGNORECASE):
                return val
            # 纯数字格式：1.1.1.1
            if re.match(r'^\d+\.\d+[\.\-]*\d*[\.\-]*\d*$', val):
                return f"EC {val}"

    return None


def _extract_chain(rule: Dict[str, Any]) -> Optional[str]:
    """从规则字典中提取蛋白质链 ID。

    参数
    ----
    rule : dict
        规则字典。

    返回
    ----
    str | None
        链 ID（单字母），或 None。
    """
    for key in ['chain', 'chain_id', 'pdb_chain']:
        if key in rule and rule[key]:
            val = str(rule[key]).strip()
            if val and len(val) <= 2 and val != 'None':
                return val.upper()

    return None


def _make_link(url: str, label: str, available: bool = True) -> Dict[str, Any]:
    """创建标准化的链接字典。

    参数
    ----
    url : str
        链接 URL。
    label : str
        显示标签。
    available : bool
        链接是否可用。

    返回
    ----
    dict
        包含 url, label, available 的字典。
    """
    return {
        "url": url,
        "label": label,
        "available": available,
    }


# ===========================================================================
# 核心公共接口
# ===========================================================================

def get_rule_links(rule: Dict[str, Any]) -> Dict[str, Any]:
    """为反应规则生成外部数据库交叉引用链接。

    从规则字典中提取 M-CSA ID、机制 ID、步骤 ID、EC 编号等信息，
    构建指向各外部数据库的超链接。

    参数
    ----
    rule : dict
        规则字典，应包含以下字段（可选）：
        - mcsa_id: M-CSA 条目 ID
        - mechanism_id: M-CSA 机制 ID
        - step_id: M-CSA 步骤 ID
        - rule_id: 规则 ID
        - ec_number: EC 编号
        - enzyme: 酶名称
        - reaction_smarts: 反应 SMARTS

    返回
    ----
    dict
        包含各数据库链接信息的字典。每个链接格式为：
        ``{url, label, available}``。
        额外包含 ``integrity`` 字段，记录数据的完整性检查结果。
    """
    mcsa_id = rule.get('mcsa_id')
    mechanism_id = rule.get('mechanism_id')
    step_id = rule.get('step_id')
    rule_id = rule.get('rule_id', '')

    pdb_id = _extract_pdb_id(mcsa_id)
    uniprot_id = _extract_uniprot_id(rule)
    ec_number = _extract_ec_number(rule)
    chain = _extract_chain(rule)

    # ---- M-CSA 链接 ----
    mcsa_available = mcsa_id is not None and str(mcsa_id).strip() not in ('None', '', '0')
    mcsa_url = ""
    if mcsa_available:
        mcsa_url = f"{MCSA_BASE}/{mcsa_id}"

    links: Dict[str, Any] = {
        "mcsa": _make_link(mcsa_url, "M-CSA 条目", mcsa_available),
    }

    # ---- M-CSA 机制链接 ----
    mechanism_available = mechanism_id is not None and str(mechanism_id).strip() not in ('None', '', '0')
    mechanism_url = ""
    if mechanism_available:
        mechanism_url = f"{MCSA_BASE}/mechanisms/{mechanism_id}"

    links["mcsa_mechanism"] = _make_link(
        mechanism_url,
        "M-CSA 机制详情",
        mechanism_available,
    )

    # ---- M-CSA 步骤链接 ----
    step_available = (
        step_id is not None and str(step_id).strip() not in ('None', '', '0')
        and mechanism_available
    )
    step_url = ""
    if step_available:
        step_url = f"{MCSA_BASE}/mechanisms/{mechanism_id}/steps/{step_id}"

    links["mcsa_step"] = _make_link(
        step_url,
        "M-CSA 步骤详情",
        step_available,
    )

    # ---- PDB 链接 ----
    pdb_available = pdb_id is not None
    pdb_url = ""
    if pdb_available:
        pdb_url = f"{PDB_BASE}/{pdb_id}"

    links["pdb"] = _make_link(pdb_url, "RCSB PDB", pdb_available)

    # ---- PDBe 链接 ----
    pdbe_url = ""
    if pdb_available:
        pdbe_url = f"{PDBE_BASE}/{pdb_id}"

    links["pdbe"] = _make_link(pdbe_url, "PDBe", pdb_available)

    # ---- PDBsum 链接 ----
    pdbsum_url = ""
    if pdb_available:
        pdbsum_url = f"{PDBSUM_BASE}/{pdb_id}"

    links["pdbsum"] = _make_link(pdbsum_url, "PDBsum", pdb_available)

    # ---- CATH 链接 ----
    cath_available = pdb_id is not None and chain is not None
    cath_url = ""
    if cath_available:
        cath_url = f"{CATH_BASE}/{pdb_id}/{chain}"

    links["cath"] = _make_link(cath_url, "CATH", cath_available)

    # ---- UniProt 链接 ----
    uniprot_available = uniprot_id is not None
    uniprot_url = ""
    if uniprot_available:
        uniprot_url = f"{UNIPROT_BASE}/{uniprot_id}"

    links["uniprot"] = _make_link(uniprot_url, "UniProt", uniprot_available)

    # ---- KEGG 链接 ----
    kegg_id = rule.get('kegg_id', rule.get('kegg_compound_id', ''))
    kegg_available = kegg_id and str(kegg_id).strip() not in ('None', '')
    kegg_url = ""
    if kegg_available:
        kegg_url = f"{KEGG_COMPOUND_BASE}?cpd:{kegg_id}"

    links["kegg"] = _make_link(kegg_url, "KEGG Compound", kegg_available)

    # ---- BRENDA 链接 ----
    brenda_available = ec_number is not None
    brenda_url = ""
    if brenda_available:
        # BRENDA URL 使用 EC 编号
        ec_clean = ec_number.replace('EC ', '').strip()
        brenda_url = f"{BRENDA_BASE}/enzyme.php?ecno={quote(ec_clean)}"

    links["brenda"] = _make_link(brenda_url, "BRENDA", brenda_available)

    # ---- ExPASy 链接 ----
    expasy_available = ec_number is not None
    expasy_url = ""
    if expasy_available:
        ec_clean = ec_number.replace('EC ', '').strip()
        expasy_url = f"{EXPASY_BASE}/EC/{ec_clean}"

    links["expasy"] = _make_link(expasy_url, "ExPASy", expasy_available)

    # ---- MACiE 链接（参考） ----
    macie_id = rule.get('macie_id', '')
    macie_available = macie_id and str(macie_id).strip() not in ('None', '')
    macie_url = ""
    if macie_available:
        macie_url = f"{MACIE_BASE}/"

    links["macie"] = _make_link(macie_url, "MACiE", macie_available)

    # ---- FunTree 链接（参考） ----
    links["funtree"] = _make_link(FUNTREE_BASE, "FunTree", False)

    # ---- InterPro 链接 ----
    interpro_id = rule.get('interpro_id', '')
    interpro_available = interpro_id and str(interpro_id).strip() not in ('None', '')
    interpro_url = ""
    if interpro_available:
        interpro_url = f"{INTERPRO_BASE}/{interpro_id}"

    links["interpro"] = _make_link(interpro_url, "InterPro", interpro_available)

    # ---- 数据完整性检查 ----
    integrity = {
        "has_mcsa_id": bool(mcsa_id and str(mcsa_id).strip() not in ('None', '', '0')),
        "has_mechanism_id": bool(
            mechanism_id and str(mechanism_id).strip() not in ('None', '', '0')
        ),
        "has_ec": ec_number is not None,
        "has_pdb": pdb_id is not None,
        "has_uniprot": uniprot_id is not None,
        "has_kegg": bool(kegg_id and str(kegg_id).strip() not in ('None', '')),
        "completeness_score": 0,  # 0-5，表示链接完整度
    }

    # 计算完整度分数
    score = 0
    if integrity["has_mcsa_id"]:
        score += 1
    if integrity["has_mechanism_id"]:
        score += 1
    if integrity["has_ec"]:
        score += 1
    if integrity["has_pdb"]:
        score += 1
    if integrity["has_uniprot"]:
        score += 1
    integrity["completeness_score"] = score

    links["integrity"] = integrity

    return links


def get_molecule_links(smiles: str) -> Dict[str, Any]:
    """为分子 SMILES 生成化学数据库搜索链接。

    生成指向 PubChem、ChemSpider、ChEBI 等数据库的搜索链接，
    供研究人员进一步查看分子的详细化学信息。

    参数
    ----
    smiles : str
        分子 SMILES 字符串。

    返回
    ----
    dict
        包含各数据库链接信息的字典。每个链接格式为：
        ``{url, label, available}``。
    """
    if not smiles:
        return {
            "error": "空 SMILES 字符串",
            "pubchem": _make_link("", "PubChem", False),
            "chemspider": _make_link("", "ChemSpider", False),
            "chebi_search": _make_link("", "ChEBI 搜索", False),
        }

    encoded_smiles = quote(smiles)

    return {
        "pubchem": _make_link(
            f"{PUBCHEM_BASE}/rest/pug/compound/smiles/{encoded_smiles}/PNG",
            "PubChem",
            True,
        ),
        "chemspider": _make_link(
            f"{CHEMSPIDER_BASE}/Search.aspx?q={encoded_smiles}",
            "ChemSpider",
            True,
        ),
        "chebi_search": _make_link(
            f"{CHEBI_BASE}/searchString.do?searchString={encoded_smiles}",
            "ChEBI 搜索",
            True,
        ),
    }


def get_ec_links(ec_number: str) -> Dict[str, Any]:
    """为 EC 编号生成酶学数据库链接。

    参数
    ----
    ec_number : str
        EC 编号字符串（例如 "EC 1.1.1.1"）。

    返回
    ----
    dict
        包含 BRENDA 和 ExPASy 链接的字典。
    """
    if not ec_number:
        return {
            "brenda": _make_link("", "BRENDA", False),
            "expasy": _make_link("", "ExPASy", False),
        }

    ec_clean = ec_number.replace('EC ', '').strip()

    return {
        "brenda": _make_link(
            f"{BRENDA_BASE}/enzyme.php?ecno={quote(ec_clean)}",
            "BRENDA",
            True,
        ),
        "expasy": _make_link(
            f"{EXPASY_BASE}/EC/{ec_clean}",
            "ExPASy",
            True,
        ),
    }


def get_all_links(
    rule: Dict[str, Any],
    substrate_smiles: str = "",
    product_smiles: str = "",
) -> Dict[str, Any]:
    """综合生成规则和分子的所有外部链接。

    参数
    ----
    rule : dict
        规则字典。
    substrate_smiles : str
        底物分子 SMILES 字符串。
    product_smiles : str
        产物分子 SMILES 字符串。

    返回
    ----
    dict
        包含 rule_links 和 molecule_links 的字典。
    """
    result: Dict[str, Any] = {
        "rule_links": get_rule_links(rule),
    }

    if substrate_smiles:
        result["substrate_links"] = get_molecule_links(substrate_smiles)

    if product_smiles:
        result["product_links"] = get_molecule_links(product_smiles)

    return result


# ===========================================================================
# CLI 接口
# ===========================================================================

def main():
    """命令行入口：生成外部数据库链接。

    用法
    ----
    # 为规则生成链接（json_input 为 JSON 编码的规则字典）
    python3 database_links.py rule '{"mcsa_id": 123, "mechanism_id": 1, "step_id": 1, "ec_number": "EC 1.1.1.1"}'

    # 为分子生成链接（json_input 为 SMILES 字符串）
    python3 database_links.py molecule "CCO"

    输出
    ----
    JSON 格式到 stdout。
    """
    if len(sys.argv) < 3:
        print(json.dumps({
            "error": "用法: python3 database_links.py <action> <json_input>\n"
                     "  action 'rule': json_input 为 JSON 编码的规则字典\n"
                     "  action 'molecule': json_input 为 SMILES 字符串"
        }, ensure_ascii=False))
        sys.exit(1)

    action = sys.argv[1].lower().strip()
    json_input = sys.argv[2]

    if action == 'rule':
        try:
            rule = json.loads(json_input)
            if not isinstance(rule, dict):
                print(json.dumps({"error": "JSON 输入必须是一个字典"}, ensure_ascii=False))
                sys.exit(1)
            result = get_rule_links(rule)
        except json.JSONDecodeError as e:
            print(json.dumps({"error": f"JSON 解析失败: {e}"}, ensure_ascii=False))
            sys.exit(1)

    elif action == 'molecule':
        result = get_molecule_links(json_input)

    elif action == 'ec':
        result = get_ec_links(json_input)

    elif action == 'all':
        try:
            rule = json.loads(json_input)
            if not isinstance(rule, dict):
                print(json.dumps({"error": "JSON 输入必须是一个字典"}, ensure_ascii=False))
                sys.exit(1)
            result = get_all_links(rule)
        except json.JSONDecodeError as e:
            print(json.dumps({"error": f"JSON 解析失败: {e}"}, ensure_ascii=False))
            sys.exit(1)

    else:
        print(json.dumps({
            "error": f"未知操作: {action}。支持的操作: rule, molecule, ec, all"
        }, ensure_ascii=False))
        sys.exit(1)

    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
