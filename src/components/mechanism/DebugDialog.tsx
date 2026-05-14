"use client";

import React, { useState, useEffect, useCallback } from "react";
import { Bug, Dna, Loader2 } from "lucide-react";
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogTrigger,
} from "@/components/ui/dialog";
import { Button } from "@/components/ui/button";
import type { SelectedResidue, MechanismSearchResult } from "@/lib/types";
import { RESIDUE_SMILES_MAP } from "@/lib/constants";

interface DebugDialogProps {
  selectedResidues: SelectedResidue[];
  reactionSmiles: string;
  searchResult: MechanismSearchResult | null;
}

export default function DebugDialog({
  selectedResidues,
  reactionSmiles,
  searchResult,
}: DebugDialogProps) {
  const [open, setOpen] = useState(false);
  const [explicitHData, setExplicitHData] = useState<{
    reactants: string;
    products: string;
  } | null>(null);
  const [loadingH, setLoadingH] = useState(false);

  // Split reaction SMILES by ">>" to get reactant/product parts
  const reactionParts = reactionSmiles.split(">>");
  const reactantSmiles = (reactionParts[0] || "").trim();
  const productSmiles = reactionParts.length > 1 ? (reactionParts[1] || "").trim() : "";

  // Build the complete config molecules list (what actually gets sent to backend)
  // Substrates from Step 4 + ALL residue SMILES from Step 2 (no dedup)
  const reactantConfigMols = React.useMemo(() => {
    const mols: { smi: string; label: string; color: string }[] = [];
    // Substrates from Step 4
    if (reactantSmiles.trim()) {
      reactantSmiles.trim().split('.').forEach((s, i) => {
        if (s.trim()) mols.push({ smi: s.trim(), label: `Step 4 substrate #${i + 1}`, color: "text-gray-800" });
      });
    }
    // ALL residue SMILES from Step 2 (no dedup)
    selectedResidues.forEach((sr, i) => {
      const smi = RESIDUE_SMILES_MAP[sr.res_name?.toUpperCase()];
      if (smi) {
        mols.push({ smi, label: `Step 2 ${sr.res_name}${sr.res_num}${sr.chain}`, color: "text-violet-700" });
      }
    });
    return mols;
  }, [reactantSmiles, selectedResidues]);

  const productConfigMols = React.useMemo(() => {
    const mols: { smi: string; label: string; color: string }[] = [];
    // Products from Step 4
    if (productSmiles.trim()) {
      productSmiles.trim().split('.').forEach((s, i) => {
        if (s.trim()) mols.push({ smi: s.trim(), label: `Step 4 product #${i + 1}`, color: "text-gray-800" });
      });
    }
    // ALL residue SMILES from Step 2 (no dedup)
    selectedResidues.forEach((sr, i) => {
      const smi = RESIDUE_SMILES_MAP[sr.res_name?.toUpperCase()];
      if (smi) {
        mols.push({ smi, label: `Step 2 ${sr.res_name}${sr.res_num}${sr.chain}`, color: "text-violet-700" });
      }
    });
    return mols;
  }, [productSmiles, selectedResidues]);

  // Fetch explicit-H SMILES when dialog opens
  const fetchExplicitH = useCallback(async () => {
    // Build full SMILES strings for both configs
    const rSmiles = reactantConfigMols.map(m => m.smi).join('.');
    const pSmiles = productConfigMols.map(m => m.smi).join('.');
    if (!rSmiles && !pSmiles) return;

    setLoadingH(true);
    try {
      const results = await Promise.allSettled([
        rSmiles ? fetch('/api', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ action: 'smiles-explicit-h', smiles: rSmiles }),
        }).then(r => r.json()) : Promise.resolve(null),
        pSmiles ? fetch('/api', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ action: 'smiles-explicit-h', smiles: pSmiles }),
        }).then(r => r.json()) : Promise.resolve(null),
      ]);
      setExplicitHData({
        reactants: results[0].status === 'fulfilled' && results[0].value?.explicit_h_smiles
          ? results[0].value.explicit_h_smiles : rSmiles,
        products: results[1].status === 'fulfilled' && results[1].value?.explicit_h_smiles
          ? results[1].value.explicit_h_smiles : pSmiles,
      });
    } catch {
      setExplicitHData({ reactants: rSmiles, products: pSmiles });
    } finally {
      setLoadingH(false);
    }
  }, [reactantConfigMols, productConfigMols]);

  useEffect(() => {
    if (open) {
      setExplicitHData(null);
      fetchExplicitH();
    }
  }, [open, fetchExplicitH]);

  // Count unique residue types vs total selected
  const uniqueResTypes = new Set(selectedResidues.map(r => r.res_name?.toUpperCase())).size;
  const unknownResidues = selectedResidues.filter(
    sr => !RESIDUE_SMILES_MAP[sr.res_name?.toUpperCase()]
  );

  return (
    <Dialog open={open} onOpenChange={setOpen}>
      <DialogTrigger asChild>
        <Button variant="outline" size="sm" className="h-7 text-xs border-amber-300 text-amber-700 hover:bg-amber-50 gap-1">
          <Bug className="w-3 h-3" />Debug
        </Button>
      </DialogTrigger>
      <DialogContent className="max-w-3xl max-h-[85vh] overflow-y-auto">
        <DialogHeader>
          <DialogTitle className="flex items-center gap-2 text-sm">
            <Bug className="w-4 h-4 text-amber-600" />完整匹配数据调试
          </DialogTitle>
        </DialogHeader>
        <div className="space-y-3 text-xs">
          {/* ===== Section 1: Step 2 选择的残基 ===== */}
          <div className="rounded-lg border border-violet-200 bg-violet-50 p-3 space-y-1.5">
            <p className="font-medium text-violet-800 text-xs flex items-center gap-1">
              <Dna className="w-3 h-3" />
              Step 2 选择的催化残基（{selectedResidues.length} 个，{uniqueResTypes} 种类型）
              <span className="text-violet-500 font-normal">→ 每个残基对应一个侧链 SMILES 参与匹配</span>
            </p>
            {selectedResidues.length === 0 ? (
              <p className="text-amber-700">⚠ 未选择任何残基！残基作为催化剂参与反应，没有残基大部分规则无法匹配。</p>
            ) : (
              <>
                {unknownResidues.length > 0 && (
                  <p className="text-red-600">⚠ {unknownResidues.length} 个残基类型未知（无对应侧链 SMILES）：
                    {unknownResidues.map((r, i) => <span key={i} className="ml-1 font-bold">{r.res_name}</span>)}
                  </p>
                )}
                {selectedResidues.map((sr, i) => {
                  const smi = RESIDUE_SMILES_MAP[sr.res_name?.toUpperCase()];
                  return (
                    <div key={i} className="ml-1 space-y-0.5">
                      <span className="text-gray-600">
                        #{i + 1} [{sr.res_name}] chain={sr.chain} num={sr.res_num} part={sr.part}
                      </span>
                      {smi ? (
                        <div className="text-gray-800 font-mono break-all">→ {smi}</div>
                      ) : (
                        <div className="text-red-500 font-mono">→ (未知残基类型，无法生成 SMILES)</div>
                      )}
                    </div>
                  );
                })}
              </>
            )}
          </div>

          {/* ===== Section 2: Step 4 画的反应物/产物 ===== */}
          <div className="rounded-lg border border-amber-200 bg-amber-50 p-3 space-y-1.5">
            <p className="font-medium text-gray-700 text-xs">Step 4 绘制的底物分子（化学编辑器输出的 SMILES）</p>
            <div>
              <span className="text-gray-500">反应物 SMILES：</span>
              <code className="break-all block mt-0.5 bg-white rounded p-1">{reactantSmiles || "（空）"}</code>
            </div>
            <div>
              <span className="text-gray-500">产物 SMILES：</span>
              <code className={`break-all block mt-0.5 rounded p-1 ${!productSmiles.trim() ? "bg-red-50 text-red-700 font-bold" : "bg-white"}`}>
                {productSmiles || "（空）"}
              </code>
              {!productSmiles.trim() && (
                <p className="text-amber-600 mt-1">⚠ 产物为空！将以 forward-only 模式搜索（无双向 meeting point）。</p>
              )}
            </div>
          </div>

          {/* ===== Section 3: 将发给后端的完整 Config（含显式 H） ===== */}
          <div className="rounded-lg border border-blue-200 bg-blue-50 p-3 space-y-2">
            <div className="flex items-center justify-between">
              <p className="font-medium text-blue-800 text-xs">📋 后端实际匹配使用的完整 Config</p>
              {loadingH && <Loader2 className="w-3 h-3 animate-spin text-blue-500" />}
            </div>
            <p className="text-gray-500">后端将 Step 4 底物 + Step 2 所有残基侧链 SMILES 组合为 Config，不自动添加水。</p>

            {/* Reactant Config */}
            <div className="space-y-1">
              <p className="text-gray-600 font-medium">Reactant Config（{reactantConfigMols.length} 个分子）：</p>
              {reactantConfigMols.map((m, i) => (
                <div key={`r${i}`} className={`ml-2 font-mono break-all ${m.color}`}>
                  [{i + 1}] {m.smi}
                  <span className="text-gray-400 ml-2 text-[10px]">← {m.label}</span>
                </div>
              ))}
              {reactantConfigMols.length === 0 && (
                <div className="ml-2 text-gray-400">（空）</div>
              )}
              {/* Explicit H version */}
              {explicitHData && (
                <div className="ml-2 mt-1 p-1.5 bg-white rounded border border-blue-100">
                  <span className="text-gray-400 text-[10px]">显式 H SMILES：</span>
                  <code className="break-all block mt-0.5 text-[11px] text-blue-800">{explicitHData.reactants}</code>
                </div>
              )}
            </div>

            {/* Product Config */}
            {productSmiles.trim() && (
              <div className="space-y-1">
                <p className="text-gray-600 font-medium">Product Config（{productConfigMols.length} 个分子）：</p>
                {productConfigMols.map((m, i) => (
                  <div key={`p${i}`} className={`ml-2 font-mono break-all ${m.color}`}>
                    [{i + 1}] {m.smi}
                    <span className="text-gray-400 ml-2 text-[10px]">← {m.label}</span>
                  </div>
                ))}
                {/* Explicit H version */}
                {explicitHData && (
                  <div className="ml-2 mt-1 p-1.5 bg-white rounded border border-blue-100">
                    <span className="text-gray-400 text-[10px]">显式 H SMILES：</span>
                    <code className="break-all block mt-0.5 text-[11px] text-blue-800">{explicitHData.products}</code>
                  </div>
                )}
              </div>
            )}
          </div>

          {/* ===== Section 4: 后端实际返回的 debug_input（搜索执行后） ===== */}
          {searchResult?.debug_input ? (
            <div className="rounded-lg border border-gray-200 bg-gray-50 p-3 space-y-1.5">
              <p className="font-medium text-gray-700 text-xs">后端实际处理结果</p>

              {/* Show actual parsed molecule lists from backend */}
              {searchResult.debug_input.reactant_mols_parsed && (
                <div className="mt-1 space-y-1">
                  <p className="text-gray-500">后端实际 reactant_mols（{searchResult.debug_input.reactant_mols_parsed.length} 个）：</p>
                  <div className="ml-2 bg-white rounded p-1.5 font-mono text-[11px] break-all border">
                    {searchResult.debug_input.reactant_mols_parsed.join(' . ')}
                  </div>
                </div>
              )}
              {searchResult.debug_input.product_mols_parsed && (
                <div className="mt-1 space-y-1">
                  <p className="text-gray-500">后端实际 product_mols（{searchResult.debug_input.product_mols_parsed.length} 个）：</p>
                  <div className="ml-2 bg-white rounded p-1.5 font-mono text-[11px] break-all border">
                    {searchResult.debug_input.product_mols_parsed.join(' . ')}
                  </div>
                </div>
              )}
              {searchResult.debug_input.residue_smiles && searchResult.debug_input.residue_smiles.length > 0 && (
                <div className="mt-1 space-y-1">
                  <p className="text-violet-600">后端添加的残基 SMILES（{searchResult.debug_input.residue_smiles.length} 个）：</p>
                  {searchResult.debug_input.residue_smiles.map((smi: string, i: number) => (
                    <div key={i} className="ml-2 text-violet-700 font-mono text-[11px] break-all">[{i + 1}] {smi}</div>
                  ))}
                </div>
              )}

              {searchResult.debug_input.react_key_eq_prod_key === true && (
                <div className="rounded border border-red-200 bg-red-50 p-2 mt-2">
                  <p className="font-bold text-red-700">⚠ 反应物 = 产物，搜索未执行！</p>
                </div>
              )}
              <div className="grid grid-cols-2 gap-1 mt-2">
                <div className="flex justify-between"><span className="text-gray-500">规则数</span><span className="font-mono">{searchResult.debug_input.total_rules}</span></div>
                <div className="flex justify-between"><span className="text-gray-500">反应物=产物</span><span className={`font-mono ${searchResult.debug_input.react_key_eq_prod_key === true ? "text-red-600 font-bold" : "text-emerald-600"}`}>{String(searchResult.debug_input.react_key_eq_prod_key)}</span></div>
                <div className="flex justify-between"><span className="text-gray-500">正向探索</span><span className="font-mono">{searchResult.debug_input.forward_configs_explored} configs</span></div>
                <div className="flex justify-between"><span className="text-gray-500">反向探索</span><span className="font-mono">{searchResult.debug_input.backward_configs_explored} configs</span></div>
              </div>
              {searchResult.stats && (
                <pre className="mt-2 bg-white rounded p-2 text-[11px] overflow-x-auto border">{JSON.stringify(searchResult.stats, null, 2)}</pre>
              )}
            </div>
          ) : (
            <div className="text-center py-3 text-gray-400 text-xs border rounded-lg border-dashed">
              搜索尚未执行
            </div>
          )}
        </div>
      </DialogContent>
    </Dialog>
  );
}
