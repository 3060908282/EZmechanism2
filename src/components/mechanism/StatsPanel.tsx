"use client";

import { Database, CheckCircle2, ArrowRight, Clock, Bug, AlertTriangle, Network, Atom } from "lucide-react";
import { Progress } from "@/components/ui/progress";
import { Table, TableBody, TableCell, TableRow } from "@/components/ui/table";
import { Collapsible, CollapsibleContent, CollapsibleTrigger } from "@/components/ui/collapsible";
import { MechanismSearchResult } from "@/lib/types";
import StatCard from "@/components/mechanism/StatCard";

export default function StatsPanel({
  searchResult,
}: {
  searchResult: MechanismSearchResult | null;
}) {
  const graphNodes = searchResult?.graph?.nodes?.length ?? 0;
  const graphEdges = searchResult?.graph?.edges?.length ?? 0;
  const pathsCount = searchResult?.paths?.length ?? 0;
  const searchTime = (searchResult?.search_time ?? searchResult?.total_time ?? 0) as number;
  const totalRules = searchResult?.debug_input?.total_rules ?? 0;
  const fwdConfigs = searchResult?.debug_input?.forward_configs_explored ?? 0;
  const bwdConfigs = searchResult?.debug_input?.backward_configs_explored ?? 0;

  const debug = searchResult?.debug_input;
  const isTrivialCase = debug?.react_key_eq_prod_key === true;

  return (
    <div className="space-y-4">
      {/* Debug Input — always visible when searchResult exists */}
      {debug && (
        <div className={`rounded-lg border p-3.5 space-y-2 ${isTrivialCase ? "border-red-200 bg-red-50" : "border-amber-200 bg-amber-50"}`}>
          <div className="flex items-center gap-2">
            <Bug className="w-3.5 h-3.5 text-amber-600" />
            <span className="text-xs font-semibold text-gray-700">搜索输入数据</span>
            {isTrivialCase && (
              <span className="flex items-center gap-1 text-[10px] font-medium text-red-600 bg-red-100 px-1.5 py-0.5 rounded">
                <AlertTriangle className="w-3 h-3" />
                反应物 = 产物，搜索未执行！
              </span>
            )}
          </div>

          <div className="space-y-1.5 text-xs">
            <div>
              <span className="text-gray-500">原始反应物 SMILES：</span>
              <code className="ml-1 text-gray-800 break-all">{debug.reactants_smiles_raw}</code>
            </div>
            <div>
              <span className="text-gray-500">原始产物 SMILES：</span>
              <code className="ml-1 text-gray-800 break-all">{debug.products_smiles_raw}</code>
            </div>

            <Collapsible>
              <CollapsibleTrigger className="text-[11px] text-amber-700 hover:underline cursor-pointer">
                展开解析后的分子列表 ▾
              </CollapsibleTrigger>
              <CollapsibleContent className="mt-1.5 space-y-1">
                <div className="text-gray-500">反应物分子 ({debug.reactant_mols_parsed.length})：</div>
                {debug.reactant_mols_parsed.map((smi, i) => (
                  <div key={i} className="ml-3 text-gray-700 break-all">
                    <span className="text-gray-400">[{i}]</span> {smi}
                  </div>
                ))}
                <div className="text-gray-500 mt-1">产物分子 ({debug.product_mols_parsed.length})：</div>
                {debug.product_mols_parsed.map((smi, i) => (
                  <div key={i} className="ml-3 text-gray-700 break-all">
                    <span className="text-gray-400">[{i}]</span> {smi}
                  </div>
                ))}
                {debug.residue_smiles.length > 0 && (
                  <>
                    <div className="text-gray-500 mt-1">残基 SMILES ({debug.residue_smiles.length})：</div>
                    {debug.residue_smiles.map((smi, i) => (
                      <div key={i} className="ml-3 text-gray-700 break-all">
                        <span className="text-gray-400">[{i}]</span> {smi}
                      </div>
                    ))}
                  </>
                )}
              </CollapsibleContent>
            </Collapsible>
          </div>

          <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs pt-1 border-t border-amber-200/50">
            <div className="flex justify-between">
              <span className="text-gray-500">规则数</span>
              <span className="font-mono">{debug.total_rules}</span>
            </div>
            <div className="flex justify-between">
              <span className="text-gray-500">反应物=产物</span>
              <span className={`font-mono ${isTrivialCase ? "text-red-600 font-bold" : "text-emerald-600"}`}>
                {String(debug.react_key_eq_prod_key)}
              </span>
            </div>
            <div className="flex justify-between">
              <span className="text-gray-500">正向探索</span>
              <span className="font-mono">{debug.forward_configs_explored} configs</span>
            </div>
            <div className="flex justify-between">
              <span className="text-gray-500">反向探索</span>
              <span className="font-mono">{debug.backward_configs_explored} configs</span>
            </div>
          </div>
        </div>
      )}

      <div className="grid grid-cols-2 sm:grid-cols-4 gap-2.5">
        <StatCard icon={<Network className="w-3.5 h-3.5 text-violet-600" />} label="Pathways" value={pathsCount} subtext="Found paths" />
        <StatCard icon={<Atom className="w-3.5 h-3.5 text-teal-600" />} label="Configs" value={graphNodes} subtext="Explored nodes" />
        <StatCard icon={<ArrowRight className="w-3.5 h-3.5 text-amber-600" />} label="Reactions" value={graphEdges} subtext="Reaction edges" />
        <StatCard icon={<Clock className="w-3.5 h-3.5 text-gray-500" />} label="Search Time" value={`${searchTime.toFixed(1)}s`} subtext="Total time" />
      </div>

      <div className="rounded-lg border border-gray-100 p-3.5 space-y-2.5">
        <div className="flex items-center justify-between text-xs">
          <span className="font-medium text-gray-700">Search Coverage</span>
          <span className="text-emerald-600 font-semibold text-sm">{fwdConfigs + bwdConfigs} configs</span>
        </div>
        <Progress value={Math.min((fwdConfigs + bwdConfigs) / Math.max(totalRules, 1) * 100, 100)} className="h-2 [&_[data-slot=progress-indicator]]:bg-emerald-500" />
      </div>

      <div className="rounded-lg border border-gray-100 p-3.5">
        <h4 className="text-xs font-medium text-gray-700 mb-2">Performance Details</h4>
        <Table>
          <TableBody>
            {[
              ["Total Rules", totalRules.toLocaleString()],
              ["Forward Configs", fwdConfigs.toLocaleString()],
              ["Backward Configs", bwdConfigs.toLocaleString()],
              ["Graph Nodes", graphNodes.toLocaleString()],
              ["Graph Edges", graphEdges.toLocaleString()],
              ["Pathways Found", pathsCount.toLocaleString()],
              ["Search Time", `${searchTime.toFixed(1)}s`],
            ].map(([label, value]) => (
              <TableRow key={label}>
                <TableCell className="text-xs text-gray-500 py-1.5">{label}</TableCell>
                <TableCell className="text-xs font-mono text-right py-1.5">{value}</TableCell>
              </TableRow>
            ))}
          </TableBody>
        </Table>
      </div>
    </div>
  );
}
