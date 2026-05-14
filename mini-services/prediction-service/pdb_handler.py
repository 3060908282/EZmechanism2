# -*- coding: utf-8 -*-
"""
PDB 结构文件解析与活性位点预测模块。

本模块提供从 RCSB 蛋白质数据库获取 PDB 文件、解析 PDB 格式文本、
预测催化活性位点残基等核心功能。模块可作为独立 CLI 工具运行，
也可作为 Python 库导入使用。

主要功能
--------
1. ``fetch_pdb_info`` — 从 RCSB 在线获取并解析 PDB 结构
2. ``parse_pdb_text``  — 直接解析 PDB 格式文本（适用于本地文件）
3. ``get_active_site_prediction`` — 基于启发式规则的催化残基预测

依赖
----
- urllib.request（标准库，用于 RCSB 网络请求）
- common.others.AA_3TO1（氨基酸三字母→单字母编码映射）
- math（标准库，距离计算）

用法示例
--------
命令行::

    python3 pdb_handler.py fetch 1TSG
    python3 pdb_handler.py parse < pdb_file.txt
    python3 pdb_handler.py active-site 1TSG

Python API::

    from pdb_handler import fetch_pdb_info, get_active_site_prediction
    info = fetch_pdb_info("1TSG")
    active = get_active_site_prediction(info)
"""

import json
import math
import sys
import urllib.request
import urllib.error
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

# 尝试从 common.others 导入氨基酸映射；若不可用则使用内联备份
try:
    from common.others import AA_3TO1 as _AA_3TO1_RAW
    # common.others 中的键为 Title Case（如 "Ala"），需要统一为大写
    AA_3TO1: Dict[str, str] = {k.upper(): v for k, v in _AA_3TO1_RAW.items()}
except ImportError:
    AA_3TO1: Dict[str, str] = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        "ASX": "B", "GLX": "Z", "XAA": "X", "SEC": "U", "PYL": "O",
    }

# ===========================================================================
# 常量定义
# ===========================================================================

#: RCSB PDB 文件下载基础 URL
RCSB_DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

#: 常见催化残基类型及其在催化三联体/二联体中的角色映射
CATALYTIC_RESIDUE_TYPES = {
    "SER": {"role": "nucleophile", "weight": 1.0},
    "CYS": {"role": "nucleophile", "weight": 0.9},
    "HIS": {"role": "general_base", "weight": 1.0},
    "ASP": {"role": "general_acid", "weight": 0.9},
    "GLU": {"role": "general_acid", "weight": 0.85},
    "LYS": {"role": "general_base", "weight": 0.7},
    "TYR": {"role": "proton_donor", "weight": 0.6},
    "THR": {"role": "nucleophile", "weight": 0.4},
    "ASN": {"role": "general_base", "weight": 0.5},
}

#: 催化三联体 / 二联体模式（残基类型组合 → 模式名称）
TRIAD_PATTERNS = {
    frozenset({"SER", "HIS", "ASP"}): "Ser-His-Asp",
    frozenset({"SER", "HIS", "GLU"}): "Ser-His-Glu",
    frozenset({"CYS", "HIS", "ASN"}): "Cys-His-Asn",
    frozenset({"CYS", "HIS", "ASP"}): "Cys-His-Asp",
    frozenset({"SER", "HIS"}): "Ser-His",
    frozenset({"CYS", "HIS"}): "Cys-His",
    frozenset({"ASP", "HIS"}): "Asp-His",
    frozenset({"GLU", "HIS"}): "Glu-His",
}

#: 水分子残基名称
WATER_RESNAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "TIP4"}

#: 禁止纳入配体列表的 HETATM 残基名称（金属离子等小分子）
EXCLUDED_HETATM = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "TIP4"}

#: 活性位点预测中，残基与配体之间的最大距离阈值（Å）
ACTIVE_SITE_DISTANCE_THRESHOLD = 5.0


# ===========================================================================
# 内部辅助函数
# ===========================================================================


def _safe_float(val: str, default: float = 0.0) -> float:
    """安全地将字符串转换为浮点数，转换失败时返回默认值。

    参数
    ----
    val : str
        待转换的字符串。
    default : float
        转换失败时的默认返回值。

    返回
    ----
    float
        转换后的浮点数或默认值。
    """
    try:
        return float(val.strip())
    except (ValueError, AttributeError):
        return default


def _safe_int(val: str, default: int = 0) -> int:
    """安全地将字符串转换为整数，转换失败时返回默认值。

    参数
    ----
    val : str
        待转换的字符串。
    default : int
        转换失败时的默认返回值。

    返回
    ----
    int
        转换后的整数或默认值。
    """
    try:
        return int(val.strip())
    except (ValueError, AttributeError):
        return default


def _extract_int_prefix(val: str, default: int = 0) -> int:
    """从字符串中提取前缀整数（处理插入码等附加字符）。

    PDB 残基编号可能附带插入码（如 ``"128A"``, ``"1B"``），
    本函数提取数值部分并忽略后缀字符。

    参数
    ----
    val : str
        待提取的字符串。
    default : int
        未找到有效整数时的默认返回值。

    返回
    ----
    int
        提取的整数值或默认值。
    """
    stripped = val.strip()
    num_str = ""
    for ch in stripped:
        if ch.isdigit() or (ch == "-" and not num_str):
            num_str += ch
        else:
            break
    return int(num_str) if num_str else default


def _calc_distance(
    x1: float, y1: float, z1: float,
    x2: float, y2: float, z2: float,
) -> float:
    """计算两个三维坐标之间的欧几里得距离。

    参数
    ----
    x1, y1, z1 : float
        第一个点的三维坐标。
    x2, y2, z2 : float
        第二个点的三维坐标。

    返回
    ----
    float
        两点之间的欧几里得距离（Å）。
    """
    return math.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)


def _build_molecular_formula(atom_counts: Dict[str, int]) -> str:
    """根据原子计数字典构建分子式字符串（Hill 排序）。

    Hill 排序规则：先 C，再 H（若存在），其余按字母序排列。
    若分子不含碳原子，则所有元素按字母序排列。

    参数
    ----
    atom_counts : dict[str, int]
        元素符号 → 原子数量的映射。

    返回
    ----
    str
        Hill 排序的分子式字符串（如 ``"C8H15NO6"``）。
    """
    if not atom_counts:
        return ""

    parts = []
    has_carbon = "C" in atom_counts

    if has_carbon:
        # 碳优先
        c_count = atom_counts.pop("C")
        parts.append(f"C{c_count}" if c_count > 1 else "C")
        # 氢其次
        if "H" in atom_counts:
            h_count = atom_counts.pop("H")
            parts.append(f"H{h_count}" if h_count > 1 else "H")

    # 剩余元素按字母序排列
    for elem in sorted(atom_counts.keys()):
        count = atom_counts[elem]
        parts.append(f"{elem}{count}" if count > 1 else elem)

    return "".join(parts)


def _parse_dbref_records(pdb_text: str) -> Dict[str, dict]:
    """解析 DBREF / DBREF1 / DBREF2 记录，构建链→UniProt 映射。

    PDB 文件中的 DBREF 记录将 PDB 链残基编号映射到 UniProt 序列位置。
    本函数提取所有 UniProt（UNP / SWS）相关的映射记录。

    对于同一链的多条 DBREF1/DBREF2 记录，自动合并为连续区间。

    解析策略
    --------
    由于不同年代的 PDB 文件 DBREF 记录的列位置存在差异，本函数
    采用「空格分割 + 字段位置推断」的方式提取关键字段，比固定列
    宽解析更具鲁棒性。

    标准字段顺序（空格分割后）：

    ``[DBREF, PDB_ID, CHAIN, pdb_begin, pdb_end, DB_NAME, ACCESSION, ..., db_begin, db_end]``

    参数
    ----
    pdb_text : str
        PDB 格式文本。

    返回
    ----
    dict[str, dict]
        链 ID → 映射信息字典，每个字典包含：

        - ``pdb_start`` : int — PDB 残基编号起始
        - ``pdb_end`` : int — PDB 残基编号终止
        - ``uniprot_start`` : int — UniProt 序列位置起始
        - ``uniprot_end`` : int — UniProt 序列位置终止
        - ``uniprot_accession`` : str — UniProt 登录号
    """
    dbref_map: Dict[str, dict] = {}

    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec not in ("DBREF", "DBREF1", "DBREF2"):
            continue

        # 使用空格分割提取字段
        parts = line.split()
        if len(parts) < 10:
            continue

        # 字段位置说明（空格分割后）:
        # [0] = 记录类型 (DBREF/DBREF1/DBREF2)
        # [1] = PDB ID
        # [2] = 链 ID
        # [3] = PDB 序列起始（可能含插入码，如 "1A"）
        # [4] = PDB 序列终止
        # [5] = 数据库名称（如 "UNP", "SWS", "PDB"）
        # [6] = 数据库登录号
        # [7:] = 描述、dbChainID、数据库序列范围等
        # 最后两个字段应为数据库序列起始和终止

        chain_id = parts[2]
        pdb_begin = _extract_int_prefix(parts[3])
        pdb_end = _extract_int_prefix(parts[4])

        if pdb_begin == 0 and pdb_end == 0:
            continue

        db_name = parts[5].upper()
        db_id = parts[6]

        # 从尾部字段中提取数据库序列范围（最后两个整数字段）
        # parts[7:] 可能包含: 描述名、dbChainID、db_seqBegin、db_seqEnd
        db_int_fields: List[int] = []
        for p in parts[7:]:
            val = _safe_int(p)
            if val > 0:
                db_int_fields.append(val)
            elif db_int_fields:
                # 已找到数字后又遇到非数字，停止收集
                break

        if len(db_int_fields) >= 2:
            db_begin = db_int_fields[-2]
            db_end = db_int_fields[-1]
        else:
            db_begin = 0
            db_end = 0

        # 判断是否为 UniProt 数据库
        is_uniprot = db_name in ("UNP", "SWS")
        if not is_uniprot and db_id:
            first_char = db_id[0].upper()
            # UniProtKB 旧格式: O/P/Q + 5 位数字
            # UniProtKB 新格式: A0A 开头的 10 位字符
            is_uniprot = (
                first_char in ("O", "P", "Q")
                or (len(db_id) >= 6 and db_id[:3].upper() == "A0A")
            )

        if not is_uniprot or pdb_begin == 0 or db_begin == 0:
            continue

        # 同一链可能有多条记录（DBREF1/DBREF2），合并区间
        if chain_id not in dbref_map:
            dbref_map[chain_id] = {
                "pdb_start": pdb_begin,
                "pdb_end": pdb_end,
                "uniprot_start": db_begin,
                "uniprot_end": db_end,
                "uniprot_accession": db_id,
            }
        else:
            existing = dbref_map[chain_id]
            existing["pdb_start"] = min(existing["pdb_start"], pdb_begin)
            existing["pdb_end"] = max(existing["pdb_end"], pdb_end)
            existing["uniprot_start"] = min(existing["uniprot_start"], db_begin)
            existing["uniprot_end"] = max(existing["uniprot_end"], db_end)
            # 保留非空的登录号
            if not existing.get("uniprot_accession") and db_id:
                existing["uniprot_accession"] = db_id

    return dbref_map


def _get_residue_center(residue_atoms: List[dict]) -> Optional[Tuple[float, float, float]]:
    """计算残基所有重原子的几何中心坐标。

    参数
    ----
    residue_atoms : list[dict]
        残基的原子信息列表，每个原子需包含 ``x``, ``y``, ``z`` 键。

    返回
    ----
    tuple[float, float, float] | None
        几何中心坐标；若残基无有效原子则返回 None。
    """
    if not residue_atoms:
        return None

    total = len(residue_atoms)
    cx = sum(a["x"] for a in residue_atoms) / total
    cy = sum(a["y"] for a in residue_atoms) / total
    cz = sum(a["z"] for a in residue_atoms) / total
    return (cx, cy, cz)


# ===========================================================================
# 核心解析函数
# ===========================================================================


def parse_pdb_text(pdb_text: str, pdb_id: str = "UNKNOWN") -> dict:
    """解析 PDB 格式文本，提取链、残基、配体、水分子和元数据。

    本函数解析标准 PDB 文本格式，支持 ATOM 和 HETATM 记录。
    返回结构化的字典，包含蛋白质链序列、配体信息、活性位点预测等。

    参数
    ----
    pdb_text : str
        PDB 格式文本（完整的 PDB 文件内容）。
    pdb_id : str
        PDB 标识符，默认为 ``"UNKNOWN"``。用于标识输出结果。

    返回
    ----
    dict
        结构化解析结果，包含以下字段：

        - ``pdb_id`` : str — PDB 标识符
        - ``title`` : str — 结构标题
        - ``resolution`` : float | None — 分辨率（Å），X 射线晶体结构
        - ``deposition_date`` : str — 沉积日期（YYYY-MM-DD）
        - ``num_chains`` : int — 蛋白质链数量
        - ``chains`` : list[dict] — 每条链的详细信息
        - ``ligands`` : list[dict] — 配体信息列表
        - ``water_count`` : int — 水分子数量
        - ``total_atoms`` : int — 总原子数
        - ``active_site_residues`` : list[dict] — 活性位点残基预测

    异常
    ----
    ValueError
        当 PDB 文本为空或格式严重不合法时抛出。
    """
    if not pdb_text or not pdb_text.strip():
        raise ValueError("PDB 文本为空，无法解析")

    # ---- 提取元数据 ----
    title = ""
    resolution = None
    deposition_date = ""

    for line in pdb_text.splitlines():
        record_name = line[:6].strip()

        # HEADER 记录：deposition_date 在第 51-59 列（DD-MMM-YY 格式）
        if record_name == "HEADER" and len(line) >= 59:
            dep = line[50:59].strip()
            if dep and len(dep) >= 9:
                day = dep[0:2]
                mon = dep[3:6]
                year = dep[7:9]
                # 月份缩写 → 数字
                MONTH_MAP = {
                    "JAN": "01", "FEB": "02", "MAR": "03",
                    "APR": "04", "MAY": "05", "JUN": "06",
                    "JUL": "07", "AUG": "08", "SEP": "09",
                    "OCT": "10", "NOV": "11", "DEC": "12",
                }
                m = MONTH_MAP.get(mon.upper(), "01")
                # 处理两位数年份（RCSB 使用 1900s 基线，PDB 档案范围通常是 1971-1999 用 19xx，
                # 2000+ 用 20xx —— 简单判断：>=70 为 19xx，<70 为 20xx）
                y = int(year) if year.isdigit() else 0
                full_year = (1900 + y) if y >= 70 else (2000 + y)
                deposition_date = f"{full_year}-{m}-{day.zfill(2)}"

        # TITLE 记录：续行也追加
        if record_name == "TITLE":
            t = line[10:].strip()
            if t:
                title += (" " if title else "") + t

        # REMARK 2 中包含分辨率信息（X 射线晶体结构）
        if record_name == "REMARK" and line[6:10].strip() == "2":
            if "RESOLUTION." in line.upper() or "ANGSTROM" in line.upper():
                parts = line.split()
                for p in parts:
                    try:
                        val = float(p)
                        if 0.5 < val < 20.0:
                            resolution = val
                            break
                    except ValueError:
                        continue

    title = title.strip() or "UNTITLED"

    # ---- 解析 DBREF 记录（链→UniProt 映射） ----
    dbref_map = _parse_dbref_records(pdb_text)

    # ---- 解析 ATOM / HETATM 记录 ----
    # 数据结构
    chain_atoms: Dict[str, List[dict]] = defaultdict(list)       # 链 → 原子列表
    chain_residues: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: defaultdict(list))  # 链 → (res_num, icode) → 原子列表
    # 配体按 (chain, res_name, res_num) 三级键分组，确保同一 res_name 的不同实例（如多个 HOH）
    # 各自独立，可以被前端分别选中映射
    ligand_instance_groups: Dict[str, Dict[str, Dict[int, List[dict]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )  # chain → res_name → res_num → [atoms]
    water_count = 0
    total_atoms = 0

    for line in pdb_text.splitlines():
        record_name = line[:6].strip()

        if record_name not in ("ATOM", "HETATM"):
            continue

        total_atoms += 1

        # PDB 固定列宽解析
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21].strip() if len(line) > 21 else " "
        res_num = line[22:26].strip()
        insertion = line[26].strip() if len(line) > 26 else " "
        x = _safe_float(line[30:38])
        y = _safe_float(line[38:46])
        z = _safe_float(line[46:54])
        element = line[76:78].strip() if len(line) > 77 else ""

        # 推断元素符号（若 PDB 文件未提供）
        if not element and atom_name:
            for ch in atom_name:
                if ch.isalpha() and ch.upper() != "H":
                    element = ch.upper()
                    break

        atom_info = {
            "atom_name": atom_name,
            "res_name": res_name,
            "chain": chain_id,
            "res_num": res_num,
            "insertion": insertion,
            "x": x,
            "y": y,
            "z": z,
            "element": element,
        }

        if record_name == "ATOM":
            # 蛋白质原子
            chain_atoms[chain_id].append(atom_info)
            residue_key = f"{res_num}{insertion}"
            chain_residues[chain_id][residue_key].append(atom_info)

        elif record_name == "HETATM":
            # 非聚合物（配体 / 水）
            is_water = res_name in WATER_RESNAMES
            if is_water:
                water_count += 1
            # 所有 HETATM（包括水）都按 (chain, res_name, res_num) 归入配体实例表
            # 水分子也需要独立条目，以便前端映射特定 HOH 位点到水分子的 substrate
            res_num_int = _safe_int(res_num)
            ligand_instance_groups[chain_id][res_name][res_num_int].append(atom_info)

    # ---- 构建链信息 ----
    chains = []
    for cid in sorted(chain_atoms.keys()):
        residues_list = []
        seq_parts = []
        seq_counter = 0  # SEQRES 序列位置计数器（无 DBREF 时的回退方案）

        # 获取该链的 DBREF 映射（若有）
        dbref = dbref_map.get(cid)
        db_offset = None
        if dbref:
            db_offset = dbref["uniprot_start"] - dbref["pdb_start"]

        for rkey in sorted(chain_residues[cid].keys(), key=lambda k: (int(k[:-1]) if k[:-1].lstrip("-").isdigit() else 0, k[-1] if k else " ")):
            atoms = chain_residues[cid][rkey]
            if not atoms:
                continue
            first = atoms[0]
            res_name = first["res_name"]
            rnum = first["res_num"]
            ins = first["insertion"]

            rnum_int = _safe_int(rnum) if rnum.lstrip("-").isdigit() else None

            # 跳过非标准残基（非 20 种标准氨基酸），但仍记录
            one_letter = AA_3TO1.get(res_name.upper(), "X")
            is_standard = one_letter != "X"
            if is_standard:
                seq_counter += 1

            seq_parts.append(one_letter)

            # 计算 UniProt 序列位置（seq_pos）
            if (
                dbref
                and rnum_int is not None
                and dbref["pdb_start"] <= rnum_int <= dbref["pdb_end"]
            ):
                # 有 DBREF 映射且残基在映射范围内：使用 UniProt 位置
                seq_pos = rnum_int + db_offset
                res_db_offset = db_offset
            else:
                # 无 DBREF 映射或超出映射范围或非标准残基：
                # seq_pos=0 表示没有 UniProt 映射，前端将显示 "—"
                seq_pos = 0
                res_db_offset = None

            residue_info = {
                "res_name": res_name,
                "res_num": _safe_int(rnum) if rnum.lstrip("-").isdigit() else rnum,
                "insertion": ins if ins else " ",
                "chain": cid,
                "seq_pos": seq_pos,
            }
            if res_db_offset is not None:
                residue_info["seq_db_offset"] = res_db_offset

            residues_list.append(residue_info)

        chain_info = {
            "chain_id": cid,
            "num_residues": len(residues_list),
            "sequence": "".join(seq_parts),
            "residues": residues_list,
        }
        if dbref:
            chain_info["uniprot_accession"] = dbref["uniprot_accession"]

        chains.append(chain_info)

    # ---- 构建配体信息 ----
    # 每个 (chain, res_name, res_num) 实例独立输出，确保同一 res_name 的多个实例
    # （如 HOH_A_402 / HOH_A_403）可分别被前端选中映射
    ligands = []
    for cid in sorted(ligand_instance_groups.keys()):
        for res_name in sorted(ligand_instance_groups[cid].keys()):
            is_water = res_name in WATER_RESNAMES
            for rn in sorted(ligand_instance_groups[cid][res_name].keys()):
                atoms = ligand_instance_groups[cid][res_name][rn]
                if not atoms:
                    continue

                # 统计非氢原子类型
                elem_counts: Dict[str, int] = defaultdict(int)
                for a in atoms:
                    elem = a.get("element", "")
                    if elem and elem.upper() != "H":
                        elem_counts[elem.upper()] += 1
                    elif not elem:
                        # 尝试从原子名称推断（跳过 H）
                        for ch in a["atom_name"]:
                            if ch.isalpha() and ch.upper() != "H":
                                elem_counts[ch.upper()] += 1
                                break

                formula = _build_molecular_formula(dict(elem_counts))

                # 取第一个原子坐标作为代表
                rep_atom = atoms[0]

                ligand_entry = {
                    "res_name": res_name,
                    "res_num": rn,
                    "chain": cid,
                    "num_atoms": len(atoms),
                    "formula": formula,
                    "x": rep_atom["x"],
                    "y": rep_atom["y"],
                    "z": rep_atom["z"],
                    "is_water": is_water,
                }
                ligands.append(ligand_entry)

    # ---- 活性位点预测（内嵌简化版本） ----
    # 活性位点预测只用非水配体，水分子不应影响位点判定
    non_water_ligands = [lg for lg in ligands if not lg.get("is_water")]
    active_site_residues = _predict_active_site_internal(
        chains, non_water_ligands, chain_residues
    )

    return {
        "pdb_id": pdb_id.upper(),
        "title": title,
        "resolution": resolution,
        "deposition_date": deposition_date,
        "num_chains": len(chains),
        "chains": chains,
        "ligands": ligands,
        "water_count": water_count,
        "total_atoms": total_atoms,
        "active_site_residues": active_site_residues,
    }


# ===========================================================================
# 配体到活性位点几何中心的距离计算
# （原版 models.py Prediction.center_of_geometry / distance_to_center_of_geometry 的等效实现）
# ===========================================================================


def _extract_pdb_internals(pdb_text: str) -> Tuple[
    Dict[str, Dict[str, Dict[int, List[dict]]]],
    Dict[str, Dict[str, List[dict]]],
]:
    """从 PDB 文本提取内部数据结构（ligand_instance_groups + chain_residues）。

    用于 ligand-distances 计算，避免完整 parse_pdb_text 的开销。
    仅解析 ATOM/HETATM 行，构建分组结构。

    返回 (ligand_instance_groups, chain_residues) 二元组。
    """
    chain_residues: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: defaultdict(list))
    ligand_instance_groups: Dict[str, Dict[str, Dict[int, List[dict]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )

    for line in pdb_text.splitlines():
        record_name = line[:6].strip()
        if record_name not in ("ATOM", "HETATM"):
            continue

        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21].strip() if len(line) > 21 else " "
        res_num = line[22:26].strip()
        insertion = line[26].strip() if len(line) > 26 else " "
        x = _safe_float(line[30:38])
        y = _safe_float(line[38:46])
        z = _safe_float(line[46:54])
        element = line[76:78].strip() if len(line) > 77 else ""

        if not element and atom_name:
            for ch in atom_name:
                if ch.isalpha() and ch.upper() != "H":
                    element = ch.upper()
                    break

        atom_info = {
            "atom_name": atom_name,
            "res_name": res_name,
            "chain": chain_id,
            "res_num": res_num,
            "insertion": insertion,
            "x": x, "y": y, "z": z,
            "element": element,
        }

        if record_name == "ATOM":
            residue_key = f"{res_num}{insertion}"
            chain_residues[chain_id][residue_key].append(atom_info)
        elif record_name == "HETATM":
            res_num_int = _safe_int(res_num)
            ligand_instance_groups[chain_id][res_name][res_num_int].append(atom_info)

    return ligand_instance_groups, chain_residues


def compute_center_of_geometry(
    atoms: List[Dict],
) -> Optional[Tuple[float, float, float]]:
    """计算重原子坐标的几何中心。

    原版 ``models.py Prediction.center_of_geometry`` 的等效实现。
    取所有非氢原子的平均位置作为活性位点的几何中心。

    参数
    ----
    atoms : list[dict]
        原子坐标列表，每个原子含有 ``x``, ``y``, ``z``, ``element`` 字段。

    返回
    ----
    tuple[float, float, float] | None
        几何中心坐标 (cx, cy, cz)，若无重原子则返回 None。
    """
    coords = []
    for atom in atoms:
        elem = atom.get('element', '').upper()
        if elem and elem != 'H':
            coords.append((atom['x'], atom['y'], atom['z']))
    if not coords:
        return None
    n = len(coords)
    return (
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n,
    )


def compute_ligand_distances(
    ligand_instance_groups: Dict[str, Dict[str, Dict[int, List[dict]]]],
    chain_residues: Dict[str, Dict[str, List[dict]]],
    selected_residues: List[Dict],
) -> Dict[str, float]:
    """计算每个配体实例到活性位点几何中心的最小距离。

    原版 ``distance_to_center_of_geometry`` 的等效实现，但应用于所有配体（含水）。
    对每个配体残基，取其所有重原子到几何中心的最小欧几里得距离。

    参数
    ----
    ligand_instance_groups : dict
        chain → res_name → res_num → [atoms]（``parse_pdb_text`` 的中间数据结构）。
    chain_residues : dict
        chain → residue_key → [atoms]（``parse_pdb_text`` 的中间数据结构）。
    selected_residues : list[dict]
        用户选中的催化残基列表，每项含 ``chain`` 和 ``res_num`` 字段。

    返回
    ----
    dict[str, float]
        ``{res_name_chain_res_num: distance_Å}`` 映射。
    """
    # 1. 收集选中残基的所有重原子
    active_atoms: List[Dict] = []
    for sr in selected_residues:
        sr_chain = (sr.get('chain', '') or '').strip()
        sr_res_num = sr.get('res_num', 0)
        if not sr_chain:
            continue
        for rkey, atoms in chain_residues.get(sr_chain, {}).items():
            for a in atoms:
                try:
                    a_res_num = _safe_int(a.get('res_num', '0'))
                except (ValueError, TypeError):
                    continue
                if a_res_num == sr_res_num:
                    active_atoms.append(a)

    # 2. 计算几何中心
    center = compute_center_of_geometry(active_atoms)
    if center is None:
        return {}
    cx, cy, cz = center

    # 3. 遍历所有配体实例，计算最小原子距离
    result: Dict[str, float] = {}
    for chain, res_names in ligand_instance_groups.items():
        for res_name, res_nums in res_names.items():
            for res_num, atoms in res_nums.items():
                min_dist = float('inf')
                for atom in atoms:
                    elem = atom.get('element', '').upper()
                    if elem and elem != 'H':
                        d = _calc_distance(atom['x'], atom['y'], atom['z'], cx, cy, cz)
                        if d < min_dist:
                            min_dist = d
                key = f"{res_name}_{chain}_{res_num}"
                result[key] = round(min_dist, 2) if min_dist != float('inf') else 999.0

    return result


# ===========================================================================
# 活性位点预测（内部实现）
# ===========================================================================


def _predict_active_site_internal(
    chains: List[dict],
    ligands: List[dict],
    chain_residues: Dict[str, Dict[str, List[dict]]],
) -> List[dict]:
    """基于配体邻近性和催化残基类型预测活性位点。

    本函数为内部实现，被 ``parse_pdb_text`` 调用以填充
    ``active_site_residues`` 字段。

    算法步骤
    --------
    1. 找出所有潜在催化残基（Ser, His, Asp, Cys, Glu, Lys, Tyr, Thr, Asn）
    2. 计算每个催化残基到最近配体的距离
    3. 距离阈值内的残基按距离排序，赋予置信度分数
    4. 检查催化三联体 / 二联体模式匹配

    参数
    ----
    chains : list[dict]
        链信息列表（来自 ``parse_pdb_text`` 的输出）。
    ligands : list[dict]
        配体信息列表。
    chain_residues : dict
        链 → 残基 → 原子列表的嵌套字典。

    返回
    ----
    list[dict]
        预测的活性位点残基列表。
    """
    if not ligands:
        return []

    # 收集所有催化残基候选
    candidates: List[dict] = []

    for chain_info in chains:
        cid = chain_info["chain_id"]
        for res in chain_info["residues"]:
            rname = res["res_name"].upper()
            if rname in CATALYTIC_RESIDUE_TYPES:
                rkey = f"{res['res_num']}{res['insertion'] if res['insertion'] != ' ' else ''}"
                atoms = chain_residues.get(cid, {}).get(rkey, [])
                center = _get_residue_center(atoms)
                if center is not None:
                    candidates.append({
                        "res_name": rname,
                        "res_num": res["res_num"],
                        "chain": cid,
                        "center": center,
                        "atoms": atoms,
                    })

    if not candidates:
        return []

    # 计算每个候选残基到最近配体的距离
    ligand_coords = [(lg["x"], lg["y"], lg["z"]) for lg in ligands]

    for cand in candidates:
        cx, cy, cz = cand["center"]
        min_dist = float("inf")
        nearest_ligand_idx = 0

        for idx, (lx, ly, lz) in enumerate(ligand_coords):
            d = _calc_distance(cx, cy, cz, lx, ly, lz)
            if d < min_dist:
                min_dist = d
                nearest_ligand_idx = idx

        cand["nearest_ligand_dist"] = min_dist
        cand["nearest_ligand_idx"] = nearest_ligand_idx

    # 筛选距离在阈值内的残基
    nearby = [
        c for c in candidates
        if c["nearest_ligand_dist"] <= ACTIVE_SITE_DISTANCE_THRESHOLD
    ]

    if not nearby:
        # 如果没有足够近的残基，取距离最近的几个
        candidates_sorted = sorted(candidates, key=lambda c: c["nearest_ligand_dist"])
        nearby = candidates_sorted[:5]

    # 计算置信度分数
    # 基础分 = 类型权重，距离加成 = 距离越近分数越高
    for cand in nearby:
        type_info = CATALYTIC_RESIDUE_TYPES.get(cand["res_name"], {})
        base_weight = type_info.get("weight", 0.5)
        dist = cand["nearest_ligand_dist"]

        # 距离因子：0Å → 1.0, 5Å → 0.0，线性递减
        dist_factor = max(0.0, 1.0 - dist / ACTIVE_SITE_DISTANCE_THRESHOLD)
        confidence = round(base_weight * (0.5 + 0.5 * dist_factor), 2)

        cand["confidence"] = confidence
        cand["role"] = type_info.get("role", "unknown")

    # 按置信度降序排列
    nearby.sort(key=lambda c: c["confidence"], reverse=True)

    # 构建输出
    results = []
    for cand in nearby:
        lg_idx = cand["nearest_ligand_idx"]
        lg = ligands[lg_idx] if lg_idx < len(ligands) else None

        results.append({
            "res_name": cand["res_name"],
            "res_num": cand["res_num"],
            "chain": cand["chain"],
            "role": cand["role"],
            "confidence": cand["confidence"],
            "nearest_ligand_dist": round(cand["nearest_ligand_dist"], 1),
            "nearest_ligand": (
                {"res_name": lg["res_name"], "res_num": lg["res_num"], "chain": lg["chain"]}
                if lg else None
            ),
        })

    return results


# ===========================================================================
# 公开 API 函数
# ===========================================================================


def fetch_pdb_info(pdb_id: str) -> dict:
    """从 RCSB 蛋白质数据库在线获取并解析 PDB 结构。

    使用 ``urllib.request`` 从 RCSB 下载指定 PDB ID 的结构文件，
    然后调用 ``parse_pdb_text`` 进行解析。

    参数
    ----
    pdb_id : str
        4 字符 PDB 标识符（例如 ``"1TSG"``）。

    返回
    ----
    dict
        结构化解析结果，格式与 ``parse_pdb_text`` 相同。

    异常
    ----
    ValueError
        当 PDB ID 格式不合法时抛出。
    ConnectionError
        当网络请求失败时抛出。
    RuntimeError
        当 RCSB 返回错误状态码时抛出。
    """
    pdb_id = pdb_id.strip().upper()
    if not pdb_id or len(pdb_id) != 4 or not pdb_id[0].isdigit():
        raise ValueError(
            f"无效的 PDB ID: '{pdb_id}'。"
            f"PDB ID 应为 4 个字符且首字符为数字（例如 '1TSG'）。"
        )

    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id.lower())

    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "EzMechanism-PDBHandler/1.0"
        })
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                raise RuntimeError(
                    f"RCSB 返回状态码 {response.status}，"
                    f"无法获取 PDB 文件 '{pdb_id}'"
                )
            pdb_text = response.read().decode("utf-8", errors="replace")

    except urllib.error.HTTPError as e:
        if e.code == 404:
            raise RuntimeError(f"PDB ID '{pdb_id}' 在 RCSB 中不存在")
        raise ConnectionError(
            f"网络错误（HTTP {e.code}）：无法下载 '{pdb_id}' — {e.reason}"
        )
    except urllib.error.URLError as e:
        raise ConnectionError(
            f"网络连接失败：无法访问 RCSB（{e.reason}）。"
            f"请检查网络连接状态。"
        )
    except Exception as e:
        if "timed out" in str(e).lower() or "timeout" in str(e).lower():
            raise ConnectionError("请求超时（30秒），请检查网络连接后重试")
        raise

    if not pdb_text.strip():
        raise RuntimeError(f"RCSB 返回了空的 PDB 文件 '{pdb_id}'")

    return parse_pdb_text(pdb_text, pdb_id=pdb_id)


def get_active_site_prediction(pdb_info: dict) -> dict:
    """基于已解析的 PDB 信息预测催化活性位点。

    使用多层启发式规则进行预测：

    1. **配体邻近性筛选**：找出距离配体 5Å 以内的催化残基候选
    2. **残基类型评分**：根据 Ser/His/Asp/Cys/Glu/Lys/Tyr 的催化
       活性赋予基础权重
    3. **三联体模式检测**：检查预测的残基组合是否匹配已知催化
       三联体或二联体模式（如 Ser-His-Asp 丝氨酸蛋白酶三联体）

    参数
    ----
    pdb_info : dict
        ``fetch_pdb_info`` 或 ``parse_pdb_text`` 的返回值。

    返回
    ----
    dict
        活性位点预测结果，包含以下字段：

        - ``predicted_residues`` : list[dict] — 预测的催化残基列表
        - ``predicted_triad`` : str | None — 匹配的催化三联体 / 二联体名称
        - ``nearest_ligand`` : dict | None — 最近配体信息及距离

    异常
    ----
    ValueError
        当输入字典缺少必要字段时抛出。
    """
    if not pdb_info or not isinstance(pdb_info, dict):
        raise ValueError("输入的 pdb_info 必须为非空字典")

    chains = pdb_info.get("chains", [])
    ligands = pdb_info.get("ligands", [])
    active_site_residues = pdb_info.get("active_site_residues", [])

    # 若已有 active_site_residues 结果，直接使用
    predicted_residues = []

    if active_site_residues:
        for res in active_site_residues:
            predicted_residues.append({
                "res_name": res["res_name"],
                "res_num": res["res_num"],
                "chain": res["chain"],
                "role": res.get("role", "unknown"),
                "confidence": res.get("confidence", 0.0),
            })
    else:
        # 回退：重新计算（需要原始原子坐标）
        # 仅基于链的残基信息进行简化预测
        for chain_info in chains:
            for res in chain_info["residues"]:
                rname = res["res_name"].upper()
                if rname in CATALYTIC_RESIDUE_TYPES:
                    type_info = CATALYTIC_RESIDUE_TYPES[rname]
                    predicted_residues.append({
                        "res_name": rname,
                        "res_num": res["res_num"],
                        "chain": res["chain"],
                        "role": type_info["role"],
                        "confidence": round(type_info["weight"] * 0.5, 2),
                    })

    # 按置信度降序排列
    predicted_residues.sort(key=lambda r: r["confidence"], reverse=True)

    # ---- 检测催化三联体 / 二联体 ----
    predicted_triad = None
    if predicted_residues:
        # 收集前 N 个残基的类型
        top_types = set()
        for r in predicted_residues[:10]:
            top_types.add(r["res_name"])

        # 从大到小匹配模式
        for pattern, name in sorted(
            TRIAD_PATTERNS.items(),
            key=lambda item: len(item[0]),
            reverse=True,
        ):
            if pattern.issubset(top_types):
                predicted_triad = name
                break

    # ---- 最近配体信息 ----
    nearest_ligand = None
    if ligands and active_site_residues:
        # 从 active_site_residues 找距离最近的配体
        best_dist = float("inf")
        for res in active_site_residues:
            dist = res.get("nearest_ligand_dist", float("inf"))
            if dist < best_dist:
                best_dist = dist
                lg = res.get("nearest_ligand")
                if lg:
                    nearest_ligand = {
                        "res_name": lg["res_name"],
                        "res_num": lg["res_num"],
                        "chain": lg["chain"],
                        "distance": round(best_dist, 1),
                    }
    elif ligands:
        # 无法计算精确距离时，仅列出第一个配体
        lg = ligands[0]
        nearest_ligand = {
            "res_name": lg["res_name"],
            "res_num": lg["res_num"],
            "chain": lg["chain"],
            "distance": None,
        }

    return {
        "predicted_residues": predicted_residues,
        "predicted_triad": predicted_triad,
        "nearest_ligand": nearest_ligand,
    }


# ===========================================================================
# CLI 命令行接口
# ===========================================================================


def _cmd_fetch(pdb_id: str) -> None:
    """CLI 子命令：从 RCSB 获取并打印 PDB 信息。

    参数
    ----
    pdb_id : str
        PDB 标识符。
    """
    print(f"正在从 RCSB 获取 {pdb_id.upper()} ...", file=sys.stderr)
    try:
        info = fetch_pdb_info(pdb_id)
        print(json.dumps(info, indent=2, ensure_ascii=False))
    except (ValueError, ConnectionError, RuntimeError) as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


def _cmd_fetch_raw(pdb_id: str) -> None:
    """CLI 子命令：从 RCSB 获取原始 PDB 文本（不解析）。

    用于前端 3D 分子查看器，直接输出原始 PDB 文本到 stdout。

    参数
    ----
    pdb_id : str
        PDB 标识符。
    """
    pdb_id = pdb_id.strip().upper()
    url = RCSB_DOWNLOAD_URL.format(pdb_id=pdb_id.lower())
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "EzMechanism-PDBHandler/1.0"
        })
        with urllib.request.urlopen(req, timeout=30) as response:
            if response.status != 200:
                raise RuntimeError(f"RCSB returned status {response.status}")
            pdb_text = response.read().decode("utf-8", errors="replace")
            # 直接输出原始 PDB 文本（不 JSON 编码，前端接收为字符串）
            print(pdb_text, end="")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


def _cmd_parse() -> None:
    """CLI 子命令：从标准输入读取 PDB 文本并解析。"""
    pdb_text = sys.stdin.read()
    try:
        info = parse_pdb_text(pdb_text)
        print(json.dumps(info, indent=2, ensure_ascii=False))
    except ValueError as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


def _cmd_active_site(pdb_id: Optional[str] = None, info_json: Optional[str] = None) -> None:
    """CLI 子命令：预测活性位点。

    参数
    ----
    pdb_id : str | None
        PDB 标识符，用于从 RCSB 获取结构。
    info_json : str | None
        已有的 PDB 信息 JSON 字符串。
    """
    if info_json:
        try:
            info = json.loads(info_json)
        except json.JSONDecodeError as e:
            print(f"错误：JSON 解析失败 — {e}", file=sys.stderr)
            sys.exit(1)
    elif pdb_id:
        print(f"正在从 RCSB 获取 {pdb_id.upper()} ...", file=sys.stderr)
        info = fetch_pdb_info(pdb_id)
    else:
        print("错误：请提供 PDB ID 或 --info-json 参数", file=sys.stderr)
        sys.exit(1)

    try:
        result = get_active_site_prediction(info)
        print(json.dumps(result, indent=2, ensure_ascii=False))
    except (ValueError, RuntimeError) as e:
        print(f"错误：{e}", file=sys.stderr)
        sys.exit(1)


def _cmd_ligand_distances(selected_residues_json: str) -> None:
    """CLI 子命令：计算配体到活性位点几何中心的距离。

    从标准输入读取 PDB 文本，从命令行参数读取选中的残基列表。

    用法::

        python3 pdb_handler.py ligand-distances '<json>' < pdb_file.txt

    参数
    ----
    selected_residues_json : str
        JSON 字符串，选中的残基列表，如 ``[{"chain":"A","res_num":70}]``。
    """
    # 解析选中的残基
    try:
        selected_residues = json.loads(selected_residues_json)
    except json.JSONDecodeError as e:
        print(f"错误：JSON 解析失败 — {e}", file=sys.stderr)
        sys.exit(1)
    if not isinstance(selected_residues, list) or len(selected_residues) == 0:
        print(json.dumps({"error": "selected_residues must be a non-empty list", "distances": {}}))
        sys.exit(1)

    # 从标准输入读取 PDB 文本
    pdb_text = sys.stdin.read()
    if not pdb_text.strip():
        print(json.dumps({"error": "No PDB text provided (read from stdin)", "distances": {}}))
        sys.exit(1)

    # 提取内部数据结构
    ligand_instance_groups, chain_residues = _extract_pdb_internals(pdb_text)

    # 计算距离
    distances = compute_ligand_distances(ligand_instance_groups, chain_residues, selected_residues)

    # 输出结果
    result = {
        "num_selected_residues": len(selected_residues),
        "num_ligands": len(distances),
        "distances": distances,
    }
    print(json.dumps(result, ensure_ascii=False))


def main() -> None:
    """命令行入口函数。

    用法
    ----
    ::

        python3 pdb_handler.py fetch <pdb_id>
        python3 pdb_handler.py parse < pdb_file.txt
        python3 pdb_handler.py active-site <pdb_id>
        python3 pdb_handler.py active-site --info-json '<json>'
    """
    if len(sys.argv) < 2:
        print(
            "用法:\n"
            "  python3 pdb_handler.py fetch <pdb_id>           从 RCSB 获取 PDB 结构\n"
            "  python3 pdb_handler.py fetch-raw <pdb_id>       获取原始 PDB 文本（用于 3D 查看）\n"
            "  python3 pdb_handler.py parse < pdb_file.txt     解析本地 PDB 文件\n"
            "  python3 pdb_handler.py active-site <pdb_id>     预测活性位点\n"
            "  python3 pdb_handler.py active-site --info-json '<json>'  基于已有数据预测\n"
            "  python3 pdb_handler.py ligand-distances '<json>' < pdb  计算配体到活性位点距离\n",
            file=sys.stderr,
        )
        sys.exit(1)

    command = sys.argv[1]

    if command == "fetch":
        if len(sys.argv) < 3:
            print("错误：请提供 PDB ID", file=sys.stderr)
            sys.exit(1)
        _cmd_fetch(sys.argv[2])

    elif command == "fetch-raw":
        if len(sys.argv) < 3:
            print("错误：请提供 PDB ID", file=sys.stderr)
            sys.exit(1)
        _cmd_fetch_raw(sys.argv[2])

    elif command == "parse":
        _cmd_parse()

    elif command == "active-site":
        if len(sys.argv) >= 4 and sys.argv[2] == "--info-json":
            _cmd_active_site(info_json=sys.argv[3])
        elif len(sys.argv) >= 3:
            _cmd_active_site(pdb_id=sys.argv[2])
        else:
            print("错误：请提供 PDB ID 或 --info-json 参数", file=sys.stderr)
            sys.exit(1)

    elif command == "ligand-distances":
        if len(sys.argv) < 3:
            print("错误：请提供 selected_residues JSON", file=sys.stderr)
            sys.exit(1)
        _cmd_ligand_distances(sys.argv[2])

    else:
        print(f"错误：未知命令 '{command}'", file=sys.stderr)
        print(
            "可用命令: fetch, fetch-raw, parse, active-site, ligand-distances",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
