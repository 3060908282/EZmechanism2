"use client";

import { useState } from "react";
import { Atom, ChevronRight, CheckCircle2, Clock, AlertCircle, ExternalLink, FlaskConical, Network, ArrowRight, Bug, AlertTriangle, ChevronDown, ChevronUp, ChevronLeft, Maximize2, Ruler } from "lucide-react";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Separator } from "@/components/ui/separator";
import { Dialog, DialogContent, DialogHeader, DialogTitle } from "@/components/ui/dialog";
import { MechanismSearchResult } from "@/lib/types";
import SmilesDisplay from "@/components/mechanism/SmilesDisplay";
import MoleculeSVG from "@/components/mechanism/MoleculeSVG";

/** Debug input panel — shows what data was actually sent to the search engine */
function DebugInputPanel({ searchResult }: { searchResult: MechanismSearchResult }) {
  const debug = searchResult.debug_input;
  const [expanded, setExpanded] = useState(false);
  if (!debug) return null;

  const isTrivialCase = debug.react_key_eq_prod_key === true;

  return (
    <div className={`rounded-lg border p-3 space-y-2 ${isTrivialCase ? "border-red-200 bg-red-50" : "border-amber-200 bg-amber-50"}`}>
      <div className="flex items-center justify-between">
        <div className="flex items-center gap-2">
          <Bug className="w-3.5 h-3.5 text-amber-600" />
          <span className="text-xs font-semibold text-gray-700">搜索输入数据</span>
          {isTrivialCase && (
            <span className="flex items-center gap-1 text-[10px] font-bold text-red-600 bg-red-100 px-1.5 py-0.5 rounded">
              <AlertTriangle className="w-3 h-3" />
              反应物 = 产物，搜索未执行！
            </span>
          )}
        </div>
        <button onClick={() => setExpanded(!expanded)} className="text-[11px] text-amber-700 hover:underline flex items-center gap-0.5">
          {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
          {expanded ? "收起" : "展开"}
        </button>
      </div>

      {/* Always visible: key info */}
      <div className="text-xs space-y-1">
        <div className="flex gap-2">
          <span className="text-gray-500 shrink-0">反应物：</span>
          <code className="text-gray-800 break-all">{debug.reactants_smiles_raw}</code>
        </div>
        <div className="flex gap-2">
          <span className="text-gray-500 shrink-0">产物：</span>
          <code className={`break-all ${isTrivialCase ? "text-red-700 font-bold" : "text-gray-800"}`}>{debug.products_smiles_raw}</code>
        </div>
      </div>

      {/* Expandable details */}
      {expanded && (
        <div className="space-y-2 pt-1 border-t border-amber-200/50 text-xs">
          <div className="text-gray-500">解析后反应物分子 ({debug.reactant_mols_parsed.length})：</div>
          {debug.reactant_mols_parsed.map((smi, i) => (
            <div key={i} className="ml-3 text-gray-700 break-all"><span className="text-gray-400">[{i}]</span> {smi}</div>
          ))}
          <div className="text-gray-500 mt-1">解析后产物分子 ({debug.product_mols_parsed.length})：</div>
          {debug.product_mols_parsed.map((smi, i) => (
            <div key={i} className="ml-3 text-gray-700 break-all"><span className="text-gray-400">[{i}]</span> {smi}</div>
          ))}
          {debug.residue_smiles.length > 0 && (
            <>
              <div className="text-gray-500 mt-1">残基 SMILES ({debug.residue_smiles.length})：</div>
              {debug.residue_smiles.map((smi, i) => (
                <div key={i} className="ml-3 text-gray-700 break-all"><span className="text-gray-400">[{i}]</span> {smi}</div>
              ))}
            </>
          )}
          <div className="grid grid-cols-2 gap-x-4 gap-y-1 text-xs pt-1 border-t border-amber-200/50">
            <div className="flex justify-between"><span className="text-gray-500">规则数</span><span className="font-mono">{debug.total_rules}</span></div>
            <div className="flex justify-between"><span className="text-gray-500">反应物=产物</span><span className={`font-mono ${isTrivialCase ? "text-red-600 font-bold" : "text-emerald-600"}`}>{String(debug.react_key_eq_prod_key)}</span></div>
            <div className="flex justify-between"><span className="text-gray-500">正向探索</span><span className="font-mono">{debug.forward_configs_explored} configs</span></div>
            <div className="flex justify-between"><span className="text-gray-500">反向探索</span><span className="font-mono">{debug.backward_configs_explored} configs</span></div>
          </div>
        </div>
      )}
    </div>
  );
}

/** A single reaction step rendered as a visual reaction diagram */
function ReactionStepCard({
  stepNum,
  fromMols,
  toMols,
  ruleName,
  direction,
  enlarged = false,
  maxBondLength,
  score,
}: {
  stepNum: number;
  fromMols: string[];
  toMols: string[];
  ruleName?: string;
  direction?: string;
  /** Render in enlarged (zoomed) mode with bigger molecules */
  enlarged?: boolean;
  maxBondLength?: number;
  score?: number;
}) {
  const hasFrom = fromMols.length > 0;
  const hasTo = toMols.length > 0;
  const molSize = enlarged ? undefined : 135;
  const molWidth = enlarged ? 250 : undefined;
  const molHeight = enlarged ? 180 : undefined;

  return (
    <div className={`rounded-lg border bg-white ${enlarged ? "border-violet-200 p-5" : "border-gray-100 p-3"}`}>
      {/* Step header */}
      <div className="flex items-center gap-2 mb-2">
        <Badge variant="outline" className={`font-mono ${enlarged ? "text-sm px-2.5 py-0.5" : "text-[10px] px-1.5 py-0"} border-violet-200 text-violet-600`}>
          {stepNum}
        </Badge>
        {direction && direction !== "forward" && (
          <Badge variant="secondary" className={`bg-violet-50 text-violet-700 ${enlarged ? "text-xs" : "text-[10px]"}`}>
            {direction === "reverse" ? "BWD" : direction}
          </Badge>
        )}
        {ruleName && (
          <span className={`text-gray-400 truncate ${enlarged ? "text-sm" : "text-[10px]"}`}>{ruleName}</span>
        )}
        {/* Bond length badge — always visible */}
        {maxBondLength != null && maxBondLength > 0 && (
          <Badge
            variant="outline"
            className={`ml-auto shrink-0 flex items-center gap-1 ${enlarged ? "text-xs px-2 py-0.5" : "text-[10px] px-1.5 py-0"} ${
              maxBondLength > 6
                ? "border-red-300 text-red-600 bg-red-50"
                : maxBondLength > 4
                  ? "border-amber-300 text-amber-700 bg-amber-50"
                  : "border-emerald-300 text-emerald-700 bg-emerald-50"
            }`}
            title={`Maximum new bond length: ${maxBondLength.toFixed(2)} Å`}
          >
            <Ruler className={enlarged ? "w-3 h-3" : "w-2.5 h-2.5"} />
            {maxBondLength.toFixed(2)} Å
          </Badge>
        )}
        {/* Score badge */}
        {score != null && score > 0 && (
          <Badge
            variant="outline"
            className={`shrink-0 ${enlarged ? "text-xs px-2 py-0.5" : "text-[10px] px-1.5 py-0"} border-gray-300 text-gray-600 bg-gray-50`}
            title={`Step score: ${score.toFixed(3)}`}
          >
            score: {score.toFixed(1)}
          </Badge>
        )}
      </div>

      {/* Reaction diagram: From → To */}
      <div className="flex items-center gap-4" style={{ minHeight: enlarged ? 190 : 135 }}>
        {/* From molecules */}
        {hasFrom ? (
          <div className={`flex-1 flex flex-wrap ${enlarged ? "gap-3" : "gap-2"} justify-end`}>
            {fromMols.map((m, i) => (
              <MoleculeSVG key={`f-${i}`} smiles={m} size={molSize} width={molWidth} height={molHeight} className="border border-emerald-100 rounded-md" />
            ))}
          </div>
        ) : (
          <div className="flex-1 text-right text-xs text-gray-300 italic">—</div>
        )}

        {/* Arrow + rule name */}
        <div className="shrink-0 flex flex-col items-center gap-1 px-2">
          <ArrowRight className={`text-gray-400 ${enlarged ? "w-8 h-8" : "w-4 h-4"}`} />
          {ruleName && (
            <span className={`text-gray-400 text-center leading-tight line-clamp-2 ${enlarged ? "text-sm max-w-[160px]" : "text-[8px] max-w-[80px]"}`}>{ruleName}</span>
          )}
        </div>

        {/* To molecules */}
        {hasTo ? (
          <div className={`flex-1 flex flex-wrap ${enlarged ? "gap-3" : "gap-2"} justify-start`}>
            {toMols.map((m, i) => (
              <MoleculeSVG key={`t-${i}`} smiles={m} size={molSize} width={molWidth} height={molHeight} className="border border-violet-100 rounded-md" />
            ))}
          </div>
        ) : (
          <div className="flex-1 text-xs text-gray-300 italic">—</div>
        )}
      </div>
    </div>
  );
}

/** Zoomed path dialog — shows all steps of a single path at larger size */
function PathZoomDialog({
  open,
  onClose,
  pathIndex,
  steps,
  totalScore,
}: {
  open: boolean;
  onClose: () => void;
  pathIndex: number;
  steps: Array<{
    step_num?: number;
    from_config?: { molecules?: string[] };
    to_config?: { molecules?: string[] };
    rule_name?: string;
    direction?: string;
    max_bond_length?: number;
    score?: number;
  }>;
  totalScore?: number;
}) {
  const filterMols = (mols?: string[]) =>
    (mols?.filter((m) => !m.startsWith("__WATER") && !m.startsWith("__PROTON") && !m.startsWith("<rdkit")) || []);

  return (
    <Dialog open={open} onOpenChange={(o) => !o && onClose()}>
      <DialogContent className="!max-w-[96vw] w-full max-h-[90vh] overflow-y-auto p-6">
        <DialogHeader className="mb-3">
          <DialogTitle className="text-base text-violet-800">
            Path {pathIndex + 1} — {steps.length} step{steps.length > 1 ? "s" : ""}
            {totalScore != null && (
              <span className="text-sm text-gray-400 font-normal ml-2">score: {totalScore}</span>
            )}
          </DialogTitle>
        </DialogHeader>
        <div className="space-y-4">
          {steps.map((step, sIdx) => {
            const fromMols = filterMols(step.from_config?.molecules);
            const toMols = filterMols(step.to_config?.molecules);
            const hasMolData = fromMols.length > 0 || toMols.length > 0;
            if (!hasMolData) {
              return (
                <div key={sIdx} className="text-xs text-gray-400 italic p-3 bg-gray-50 rounded-lg">
                  Step {step.step_num ?? sIdx + 1}: No molecule data
                </div>
              );
            }
            return (
              <ReactionStepCard
                key={sIdx}
                stepNum={step.step_num ?? sIdx + 1}
                fromMols={fromMols}
                toMols={toMols}
                ruleName={step.rule_name}
                direction={step.direction}
                enlarged
                maxBondLength={step.max_bond_length}
                score={step.score}
              />
            );
          })}
        </div>
      </DialogContent>
    </Dialog>
  );
}

export default function MechanismPanel({
  searchResult,
  onExploreProduct,
  onSwitchToGraph,
}: {
  searchResult: MechanismSearchResult | null;
  onExploreProduct: (smiles: string) => void;
  onSwitchToGraph?: () => void;
}) {
  // Track collapsed state per path
  const [collapsedPaths, setCollapsedPaths] = useState<Set<number>>(new Set());
  // Zoom dialog
  const [zoomedPath, setZoomedPath] = useState<number | null>(null);

  const togglePath = (idx: number) => {
    setCollapsedPaths((prev) => {
      const next = new Set(prev);
      if (next.has(idx)) next.delete(idx);
      else next.add(idx);
      return next;
    });
  };

  // Show mechanism-search paths if predict has no results but search found paths
  const hasSearchPaths = searchResult && searchResult.paths && searchResult.paths.length > 0;
  const hasSearchError = searchResult && (searchResult as unknown as Record<string, unknown>).error;
  const hasSearchGraph = searchResult && searchResult.graph &&
    searchResult.graph.nodes && searchResult.graph.nodes.length > 0;

  if (hasSearchPaths) {
    const zoomedSteps = zoomedPath !== null && searchResult.paths[zoomedPath]
      ? searchResult.paths[zoomedPath].steps
      : null;

    return (
      <div className="space-y-3">
        {/* Debug panel — always at top */}
        <DebugInputPanel searchResult={searchResult} />

        <div className="flex items-center gap-2 px-3 py-2 rounded-lg bg-violet-50 border border-violet-100">
          <Badge variant="secondary" className="bg-violet-100 text-violet-800 text-xs shrink-0">
            Search
          </Badge>
          <span className="text-xs text-violet-700">
            Bidirectional mechanism search found {searchResult.paths.length} pathway(s)
          </span>
        </div>

        {searchResult.paths.map((path, pIdx) => {
          const isCollapsed = collapsedPaths.has(pIdx);
          const stepCount = path.total_steps ?? path.num_steps ?? path.steps?.length ?? 0;
          return (
            <div key={path.path_id ?? pIdx} className="rounded-lg border border-violet-100 bg-white overflow-hidden">
              {/* Path header with collapse + zoom buttons */}
              <div className="flex items-center justify-between px-3 py-1.5 bg-violet-50/50 border-b border-violet-50">
                <span className="text-xs font-medium text-violet-700">Path {(path.path_id ?? pIdx) + 1}</span>
                <div className="flex items-center gap-1">
                  <span className="text-xs text-violet-500 mr-1">{stepCount} step(s)</span>
                  <Button
                    variant="ghost"
                    size="sm"
                    className="h-6 w-6 p-0 text-gray-400 hover:text-violet-600"
                    onClick={() => setZoomedPath(pIdx)}
                    title="放大查看"
                  >
                    <Maximize2 className="w-3.5 h-3.5" />
                  </Button>
                  <Button
                    variant="ghost"
                    size="sm"
                    className="h-6 w-6 p-0 text-gray-400 hover:text-violet-600"
                    onClick={() => togglePath(pIdx)}
                    title={isCollapsed ? "展开" : "收起"}
                  >
                    {isCollapsed
                      ? <ChevronLeft className="w-3.5 h-3.5 rotate-90" />
                      : <ChevronUp className="w-3.5 h-3.5" />
                    }
                  </Button>
                </div>
              </div>
              {/* Collapsible step content */}
              {!isCollapsed && (
                <div className="p-2.5 space-y-2">
                  {path.steps.map((step, sIdx) => {
                    const fromMols = step.from_config?.molecules?.filter(
                      (m) => !m.startsWith("__WATER") && !m.startsWith("__PROTON") && !m.startsWith("<rdkit")
                    ) || [];
                    const toMols = step.to_config?.molecules?.filter(
                      (m) => !m.startsWith("__WATER") && !m.startsWith("__PROTON") && !m.startsWith("<rdkit")
                    ) || [];
                    const hasMolData = fromMols.length > 0 || toMols.length > 0;

                    if (!hasMolData) {
                      return (
                        <div key={`${pIdx}-${sIdx}`} className="flex items-start gap-2">
                          <div className="flex items-center gap-1 mt-0.5 shrink-0">
                            <Badge variant="outline" className="text-[10px] font-mono px-1.5 py-0 border-violet-200 text-violet-600">
                              {step.step_num ?? sIdx + 1}
                            </Badge>
                            <ChevronRight className="w-3 h-3 text-violet-400" />
                          </div>
                          <div className="text-xs text-gray-400 italic">No molecule data</div>
                        </div>
                      );
                    }

                    return (
                      <ReactionStepCard
                        key={`${pIdx}-${sIdx}`}
                        stepNum={step.step_num ?? sIdx + 1}
                        fromMols={fromMols}
                        toMols={toMols}
                        ruleName={step.rule_name}
                        direction={step.direction}
                        maxBondLength={step.max_bond_length}
                        score={step.score}
                      />
                    );
                  })}
                </div>
              )}
            </div>
          );
        })}

        {/* Zoom dialog */}
        {zoomedSteps && (
          <PathZoomDialog
            open={true}
            onClose={() => setZoomedPath(null)}
            pathIndex={zoomedPath!}
            steps={zoomedSteps}
            totalScore={searchResult.paths[zoomedPath!]?.total_score ?? 0}
          />
        )}

        <Separator />
        <div className="flex items-center justify-between">
          <div className="flex items-center gap-4 text-xs text-gray-500">
            <span className="flex items-center gap-1">
              <Network className="w-3 h-3 text-violet-500" />
              {searchResult.paths.length} pathway(s) found
            </span>
            <span className="flex items-center gap-1">
              <Atom className="w-3 h-3 text-teal-500" />
              {searchResult.graph?.nodes?.length || 0} configs explored
            </span>
            <span className="flex items-center gap-1">
              <Clock className="w-3 h-3 text-gray-400" />
              {((searchResult.search_time || searchResult.total_time) || 0).toFixed(1)}s
            </span>
          </div>
          {onSwitchToGraph && (
            <Button variant="ghost" size="sm" className="text-xs text-violet-600 hover:text-violet-800" onClick={onSwitchToGraph}>
              View Graph <ArrowRight className="w-3 h-3 ml-1" />
            </Button>
          )}
        </div>
      </div>
    );
  }

  // No paths found — show diagnostic info from search graph
  const graphNodes = searchResult?.graph?.nodes?.length ?? 0;
  const graphEdges = searchResult?.graph?.edges?.length ?? 0;

  return (
    <div className="space-y-3">
      {/* Debug panel — always visible when searchResult exists */}
      {searchResult && <DebugInputPanel searchResult={searchResult} />}

      <div className="text-center py-8">
        {hasSearchError ? (
          <Alert className="max-w-md mx-auto border-red-200 bg-red-50">
            <AlertCircle className="h-4 w-4 text-red-500" />
            <AlertTitle className="text-sm">Search Error</AlertTitle>
            <AlertDescription className="text-xs text-red-700">
              {(() => {
                const errStr = String((searchResult as unknown as Record<string, unknown>).error || "");
                if (errStr.includes("Unexpected token") || errStr.includes("is not valid JSON")) {
                  return "The server returned an invalid response (likely a timeout). Please try reducing Max Rules and Max Configs, then run the search again.";
                }
                return errStr || "An unknown error occurred during the mechanism search.";
              })()}
            </AlertDescription>
          </Alert>
        ) : hasSearchGraph ? (
          <Alert className="max-w-md mx-auto">
            <AlertCircle className="h-4 w-4 text-amber-500" />
            <AlertTitle className="text-sm">No Multi-Step Pathway Found</AlertTitle>
            <AlertDescription className="text-xs space-y-1">
              <p>
                The bidirectional search explored {graphNodes} configuration(s) and {graphEdges} reaction edge(s),
                but could not find a complete pathway connecting reactants to products.
              </p>
              <p className="text-gray-400">
                M-CSA rules describe multi-molecule enzymatic reactions requiring enzyme residues and
                cofactors. Check the &quot;Graph&quot; tab for the explored configuration space.
              </p>
            </AlertDescription>
          </Alert>
        ) : (
          <div>
            <AlertCircle className="w-7 h-7 text-gray-300 mx-auto mb-2" />
            <p className="text-gray-500 text-sm">No matching rules found for this substrate.</p>
          </div>
        )}
      </div>
    </div>
  );
}
