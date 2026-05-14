#!/usr/bin/env python3
"""
batch_predict.py — 批量预测模块

本模块提供了批量 SMILES 预测功能，支持并行执行多个分子的反应规则匹配和产物预测。

主要功能：
- 并行批量预测（run_batch_prediction）
- 逐 SMILES 预测（_predict_single_smiles）
- 规则匹配和产物生成（核心逻辑复用 shared.py）
- 统计摘要生成（_build_summary）
- CLI 接口支持命令行参数和 stdin 输入

依赖：rdkit, shared.py（仅限 RDKit，无 Django 框架依赖）
"""

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Any, Dict, List, Optional

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

from shared import (
    is_valid_smiles,
    simplify_smarts,
    load_rules,
    RULES_FILE,
)


# ---------------------------------------------------------------------------
# 核心逻辑：单 SMILES 预测
# ---------------------------------------------------------------------------

def _match_rules_for_smiles(
    smiles: str,
    rules: List[Dict[str, Any]],
    max_products_per_rule: int = 5,
) -> List[Dict[str, Any]]:
    """对单个 SMILES 执行规则匹配，返回匹配到的规则和产物列表。

    使用两阶段匹配策略：
    1. 先尝试简化后的 SMARTS 匹配（更灵活）
    2. 再尝试原始 SMARTS 匹配（更精确）

    Parameters
    ----------
    smiles : str
        输入分子的 SMILES 字符串。
    rules : List[Dict[str, Any]]
        预加载的反应规则列表。
    max_products_per_rule : int
        每条规则最多生成的产物数量。

    Returns
    -------
    List[Dict[str, Any]]
        匹配结果列表，每项包含 rule_id, match_score, products 等字段。
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []

    matches: List[Dict[str, Any]] = []

    for rule in rules:
        matched = False
        score = 0.0
        simplified = simplify_smarts(rule['reactant_smarts'])

        # 第一阶段：简化 SMARTS 匹配
        for comp in simplified.split('.'):
            comp = comp.strip()
            if not comp or len(comp) < 3 or len(comp) > 300:
                continue
            try:
                pat = Chem.MolFromSmarts(comp)
                if pat and mol.HasSubstructMatch(pat):
                    matched, score = True, 0.8
                    break
            except Exception:
                pass

        # 第二阶段：原始 SMARTS 匹配
        if not matched:
            for comp in rule['reactant_smarts'].split('.'):
                comp = comp.strip()
                if not comp or len(comp) > 300:
                    continue
                try:
                    pat = Chem.MolFromSmarts(comp)
                    if pat and mol.HasSubstructMatch(pat):
                        matched, score = True, 1.0
                        break
                except Exception:
                    pass

        if matched:
            products: List[str] = []
            # 尝试用反应 SMARTS 生成产物
            try:
                rxn = AllChem.ReactionFromSmarts(rule['reaction_smarts'])
                if rxn:
                    results = rxn.RunReactants([mol], maxProducts=max_products_per_rule)
                    seen: set = set()
                    for ptuple in results:
                        for p in ptuple:
                            try:
                                Chem.SanitizeMol(p)
                                smi = Chem.MolToSmiles(p)
                                if smi not in seen:
                                    seen.add(smi)
                                    products.append(smi)
                            except Exception:
                                pass
            except Exception:
                pass

            matches.append({
                'rule_id': rule.get('rule_id', ''),
                'rule_name': rule.get('rule_name', ''),
                'mcsa_id': rule.get('mcsa_id'),
                'mechanism_id': rule.get('mechanism_id'),
                'step_id': rule.get('step_id'),
                'reaction_smarts': rule.get('reaction_smarts', ''),
                'match_score': score,
                'category': rule.get('category', ''),
                'source': rule.get('source', 'mcsa'),
                'products': products,
            })

    return matches


def _predict_single_smiles(
    smiles: str,
    rules: List[Dict[str, Any]],
    max_steps: int = 3,
) -> Dict[str, Any]:
    """对单个 SMILES 执行多步预测。

    Parameters
    ----------
    smiles : str
        输入分子的 SMILES 字符串。
    rules : List[Dict[str, Any]]
        预加载的反应规则列表。
    max_steps : int
        最大预测步数。

    Returns
    -------
    Dict[str, Any]
        预测结果，包含 smiles, status, matches_count, products,
        prediction_steps, elapsed 等字段。
    """
    t0 = time.time()

    # 验证 SMILES
    if not is_valid_smiles(smiles):
        return {
            "smiles": smiles,
            "status": "error",
            "error": "Invalid SMILES characters",
            "elapsed": round(time.time() - t0, 3),
        }

    # 解析分子
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "smiles": smiles,
            "status": "error",
            "error": "Invalid SMILES: RDKit cannot parse",
            "elapsed": round(time.time() - t0, 3),
        }

    # 多步预测
    all_products: List[str] = []
    all_matches: List[Dict[str, Any]] = []
    visited = {smiles}
    current = [smiles]
    steps_completed = 0

    for step_num in range(1, max_steps + 1):
        step_products = []
        for sub in current:
            matches = _match_rules_for_smiles(sub, rules)
            seen = set()
            for m in matches:
                all_matches.append(m)
                for prod_smi in m.get('products', []):
                    if prod_smi not in seen and prod_smi not in visited:
                        seen.add(prod_smi)
                        step_products.append(prod_smi)
        if not step_products:
            break
        for p in step_products:
            visited.add(p)
            if p not in all_products:
                all_products.append(p)
        current = step_products
        steps_completed = step_num

    elapsed = time.time() - t0
    return {
        "smiles": smiles,
        "status": "success",
        "matches_count": len(all_matches),
        "products": all_products[:20],  # 限制产物数量避免输出过大
        "prediction_steps": steps_completed,
        "elapsed": round(elapsed, 3),
    }


# ---------------------------------------------------------------------------
# 批量预测
# ---------------------------------------------------------------------------

def run_batch_prediction(
    smiles_list: List[str],
    max_steps: int = 3,
    max_workers: int = 4,
) -> Dict[str, Any]:
    """并行批量预测多个 SMILES 的反应产物。

    使用 ThreadPoolExecutor 并行执行多个 SMILES 的预测任务，
    收集结果并构建统计摘要。

    Parameters
    ----------
    smiles_list : List[str]
        SMILES 字符串列表。
    max_steps : int
        每个 SMILES 的最大预测步数（1-10）。
    max_workers : int
        并行工作线程数（1-8）。

    Returns
    -------
    Dict[str, Any]
        批量预测结果，包含：
        - total_inputs: 输入总数
        - successful: 成功数
        - failed: 失败数
        - max_workers: 使用的线程数
        - elapsed_seconds: 总耗时（秒）
        - results: 每个输入的预测结果列表
        - summary: 统计摘要（总匹配数、总产物数、平均耗时等）
    """
    # 限制参数范围
    max_workers = max(1, min(8, max_workers))
    max_steps = max(1, min(10, max_steps))

    t_total = time.time()

    # 加载规则（仅加载一次，所有线程共享）
    t_load = time.time()
    rules = load_rules()
    load_time = time.time() - t_load

    results: List[Dict[str, Any]] = []

    # 使用线程池并行执行
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_smiles = {}
        for smi in smiles_list:
            future = executor.submit(
                _predict_single_smiles, smi, rules, max_steps
            )
            future_to_smiles[future] = smi

        for future in as_completed(future_to_smiles):
            try:
                result = future.result(timeout=120)
                results.append(result)
            except Exception as e:
                smi = future_to_smiles[future]
                results.append({
                    "smiles": smi,
                    "status": "error",
                    "error": str(e),
                    "elapsed": round(time.time() - t_total, 3),
                })

    # 按输入顺序排序结果
    smiles_order = {smi: i for i, smi in enumerate(smiles_list)}
    results.sort(key=lambda r: smiles_order.get(r.get("smiles", ""), 999))

    total_elapsed = time.time() - t_total

    # 构建统计摘要
    summary = _build_summary(results)

    return {
        "total_inputs": len(smiles_list),
        "successful": summary["successful"],
        "failed": summary["failed"],
        "max_workers": max_workers,
        "elapsed_seconds": round(total_elapsed, 3),
        "results": results,
        "summary": summary,
    }


# ---------------------------------------------------------------------------
# 统计摘要
# ---------------------------------------------------------------------------

def _build_summary(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """从批量预测结果中构建统计摘要。

    Parameters
    ----------
    results : List[Dict[str, Any]]
        每个输入的预测结果列表。

    Returns
    -------
    Dict[str, Any]
        统计摘要，包含：
        - successful: 成功数
        - failed: 失败数
        - total_matches: 总匹配规则数
        - total_products: 总产物数
        - avg_time_per_smiles: 每个 SMILES 平均耗时
        - most_matches: 匹配最多的 SMILES
        - most_products: 产物最多的 SMILES
    """
    successful = [r for r in results if r.get("status") == "success"]
    failed = [r for r in results if r.get("status") != "success"]

    total_matches = sum(r.get("matches_count", 0) for r in successful)
    total_products = sum(len(r.get("products", [])) for r in successful)
    total_elapsed = sum(r.get("elapsed", 0) for r in results)
    avg_time = total_elapsed / len(results) if results else 0

    # 匹配最多的 SMILES
    most_matches = {"smiles": "", "count": 0}
    for r in successful:
        count = r.get("matches_count", 0)
        if count > most_matches["count"]:
            most_matches = {"smiles": r.get("smiles", ""), "count": count}

    # 产物最多的 SMILES
    most_products = {"smiles": "", "count": 0}
    for r in successful:
        count = len(r.get("products", []))
        if count > most_products["count"]:
            most_products = {"smiles": r.get("smiles", ""), "count": count}

    return {
        "successful": len(successful),
        "failed": len(failed),
        "total_matches": total_matches,
        "total_products": total_products,
        "avg_time_per_smiles": round(avg_time, 3),
        "most_matches": most_matches,
        "most_products": most_products,
    }


# ---------------------------------------------------------------------------
# 表格输出
# ---------------------------------------------------------------------------

def _format_table(results: List[Dict[str, Any]], summary: Dict[str, Any]) -> str:
    """将批量预测结果格式化为 ASCII 表格。

    Parameters
    ----------
    results : List[Dict[str, Any]]
        预测结果列表。
    summary : Dict[str, Any]
        统计摘要。

    Returns
    -------
    str
        格式化后的 ASCII 表格字符串。
    """
    lines = []
    lines.append("=" * 90)
    lines.append("BATCH PREDICTION RESULTS")
    lines.append("=" * 90)
    lines.append(f"Total: {summary['successful']} successful, {summary['failed']} failed | "
                 f"Matches: {summary['total_matches']} | Products: {summary['total_products']} | "
                 f"Avg: {summary['avg_time_per_smiles']}s/SMILES")
    lines.append("-" * 90)

    header = f"{'#':>3}  {'SMILES':<20}  {'Status':<10}  {'Matches':>8}  {'Products':>9}  {'Steps':>6}  {'Time(s)':>8}"
    lines.append(header)
    lines.append("-" * 90)

    for i, r in enumerate(results, 1):
        smiles = r.get("smiles", "")
        if len(smiles) > 18:
            smiles = smiles[:15] + "..."
        status = r.get("status", "unknown")
        matches = r.get("matches_count", 0)
        products = len(r.get("products", []))
        steps = r.get("prediction_steps", 0)
        elapsed = r.get("elapsed", 0)

        if status == "error":
            lines.append(f"{i:>3}  {smiles:<20}  {'ERROR':<10}  {'—':>8}  {'—':>9}  {'—':>6}  {elapsed:>8.3f}")
        else:
            lines.append(f"{i:>3}  {smiles:<20}  {'OK':<10}  {matches:>8}  {products:>9}  {steps:>6}  {elapsed:>8.3f}")

    lines.append("=" * 90)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    """命令行入口：支持直接传参和 stdin 输入。

    用法：
        python3 batch_predict.py <smiles1> <smiles2> ...
        echo "CCO\\nCC(=O)O" | python3 batch_predict.py --stdin
        python3 batch_predict.py --stdin --max-steps 2 --workers 4 --output table
    """
    parser = argparse.ArgumentParser(
        description="批量预测模块 — 并行执行多个 SMILES 的反应规则匹配和产物预测",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例：
  python3 batch_predict.py CCO CC(=O)O C1=CC=CC=C1
  echo -e "CCO\\nCC(=O)O" | python3 batch_predict.py --stdin
  python3 batch_predict.py --stdin --max-steps 2 --workers 4 --output table
  python3 batch_predict.py CCO INVALID_SMILES "C1=CC=CC=C1" --output json
        """,
    )
    parser.add_argument('smiles', nargs='*', default=[],
                        help='SMILES 字符串列表（或使用 --stdin 从标准输入读取）')
    parser.add_argument('--stdin', action='store_true',
                        help='从标准输入读取 SMILES（每行一个）')
    parser.add_argument('--max-steps', type=int, default=3,
                        help='每个 SMILES 的最大预测步数（1-10，默认 3）')
    parser.add_argument('--workers', type=int, default=4,
                        help='并行工作线程数（1-8，默认 4）')
    parser.add_argument('--output', choices=['json', 'table'], default='json',
                        help='输出格式：json（默认）或 table')

    args = parser.parse_args()

    # 收集 SMILES 列表
    smiles_list: List[str] = []

    if args.stdin:
        # 从标准输入逐行读取
        for line in sys.stdin:
            line = line.strip()
            if line and not line.startswith('#'):
                smiles_list.append(line)
    else:
        smiles_list = args.smiles

    if not smiles_list:
        print(json.dumps({
            'error': '未提供 SMILES 输入。使用 --stdin 或直接传入 SMILES 参数。'
        }))
        sys.exit(1)

    # 执行批量预测
    result = run_batch_prediction(
        smiles_list=smiles_list,
        max_steps=args.max_steps,
        max_workers=args.workers,
    )

    # 格式化输出
    if args.output == 'table':
        summary = result.get('summary', {})
        print(_format_table(result.get('results', []), summary))
    else:
        print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
