"use client";

import React, { useState } from "react";
import { Search, Loader2, Dna, Link2, FlaskConical, ExternalLink } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { apiCall } from "@/lib/api";
import type { MatchRulesResponse, RuleMatch, ResidueInfoResponse, RuleLinksResponse, DbLink } from "@/lib/types";

export default function RulesPanel({
  matchResult,
  filter,
  onFilterChange,
  filteredMatches,
}: {
  matchResult: MatchRulesResponse | null;
  filter: "all" | "mcsa" | "builtin";
  onFilterChange: (f: "all" | "mcsa" | "builtin") => void;
  filteredMatches: RuleMatch[];
}) {
  const [residueCache, setResidueCache] = useState<Record<string, ResidueInfoResponse | null>>({});
  const [linksCache, setLinksCache] = useState<Record<string, RuleLinksResponse | null>>({});
  const [loadingResidue, setLoadingResidue] = useState<string | null>(null);
  const [loadingLinks, setLoadingLinks] = useState<string | null>(null);
  const [expandedResidue, setExpandedResidue] = useState<string | null>(null);

  const handleFetchResidueInfo = async (ruleId: string | number, reactionSmarts: string) => {
    const key = String(ruleId);
    if (residueCache[key] !== undefined) { setExpandedResidue(expandedResidue === key ? null : key); return; }
    setLoadingResidue(key);
    setExpandedResidue(key);
    try {
      const r = await apiCall("residue-info", { reaction_smarts: reactionSmarts }, 30000);
      if (!r.ok) throw new Error("Failed");
      const data = await r.json();
      setResidueCache((prev) => ({ ...prev, [key]: data }));
    } catch {
      setResidueCache((prev) => ({ ...prev, [key]: null }));
    } finally { setLoadingResidue(null); }
  };

  const handleFetchLinks = async (ruleId: string | number, rule: RuleMatch) => {
    const key = String(ruleId);
    if (linksCache[key] !== undefined) return;
    setLoadingLinks(key);
    try {
      const ruleData = { mcsa_id: rule.mcsa_id, mechanism_id: rule.mechanism_id, step_id: rule.step_id, enzyme: rule.enzyme, reaction_smarts: rule.reaction_smarts };
      const r = await apiCall("rule-links", { rule: ruleData }, 10000);
      if (!r.ok) throw new Error("Failed");
      const data = await r.json();
      setLinksCache((prev) => ({ ...prev, [key]: data }));
    } catch {
      setLinksCache((prev) => ({ ...prev, [key]: null }));
    } finally { setLoadingLinks(null); }
  };

  const ROLE_COLORS: Record<string, string> = {
    nucleophile: "bg-amber-100 text-amber-800 border-amber-200",
    "proton donor": "bg-emerald-100 text-emerald-800 border-emerald-200",
    "proton acceptor": "bg-blue-100 text-blue-800 border-blue-200",
    electrophile: "bg-red-100 text-red-800 border-red-200",
    "metal ion": "bg-violet-100 text-violet-800 border-violet-200",
    "general base": "bg-cyan-100 text-cyan-800 border-cyan-200",
    "general acid": "bg-orange-100 text-orange-800 border-orange-200",
    "transition state stabilizer": "bg-pink-100 text-pink-800 border-pink-200",
  };

  if (!matchResult || matchResult.matches.length === 0) {
    return (
      <div className="text-center py-8">
        <Search className="w-7 h-7 text-gray-300 mx-auto mb-2" />
        <p className="text-gray-500 text-sm">No matching rules found.</p>
      </div>
    );
  }

  const mcsaCount = matchResult.matches.filter((m) => m.source === "mcsa").length;
  const builtinCount = matchResult.matches.length - mcsaCount;

  return (
    <div className="space-y-3">
      <div className="flex items-center justify-between flex-wrap gap-2">
        <div className="flex items-center gap-2">
          <p className="text-xs text-gray-600">
            <span className="font-semibold text-emerald-700">{matchResult.matches_count.toLocaleString()}</span> of {matchResult.total_rules.toLocaleString()} rules matched
          </p>
          <span className="text-xs text-gray-400">in {matchResult.elapsed_seconds}s</span>
        </div>
        <div className="flex gap-1">
          {(["all", "mcsa", "builtin"] as const).map((f) => (
            <Button key={f} variant={filter === f ? "default" : "outline"} size="sm"
              className={`h-6 text-xs px-2 ${filter === f ? "bg-emerald-600 hover:bg-emerald-700 text-white" : "border-gray-200 hover:border-emerald-300"}`}
              onClick={() => onFilterChange(f)}>
              {f === "all" ? `All (${matchResult.matches_count})` : f === "mcsa" ? `M-CSA (${mcsaCount})` : `Built-in (${builtinCount})`}
            </Button>
          ))}
        </div>
      </div>

      <div className="max-h-[500px] overflow-y-auto space-y-1.5 pr-1">
        {filteredMatches.slice(0, 50).map((match, idx) => {
          const ruleKey = `${match.rule_id}-${idx}`;
          const cacheKey = String(match.rule_id);
          const residueInfo = residueCache[cacheKey];
          const linksInfo = linksCache[cacheKey];

          return (
            <div key={ruleKey} className="rounded-lg border border-gray-100 p-2.5 hover:border-emerald-200 hover:bg-emerald-50/30 transition-colors">
              <div className="flex items-start justify-between gap-2 mb-1.5">
                <div className="flex items-center gap-1.5 min-w-0">
                  <Badge variant="secondary" className="bg-emerald-50 text-emerald-700 text-xs shrink-0">#{match.rule_id}</Badge>
                  <span className="text-xs font-medium text-gray-800 truncate">{match.rule_name || `Rule ${match.rule_id}`}</span>
                </div>
                <div className="flex items-center gap-1 shrink-0">
                  {match.source === "mcsa" && (
                    <TooltipProvider>
                      <Tooltip>
                        <TooltipTrigger asChild>
                          <button onClick={() => handleFetchResidueInfo(match.rule_id, match.reaction_smarts)}
                            className={`p-1 rounded hover:bg-amber-50 transition-colors ${expandedResidue === cacheKey ? "text-amber-600" : "text-gray-400"}`}
                            disabled={loadingResidue === cacheKey}>
                            {loadingResidue === cacheKey ? <Loader2 className="w-3 h-3 animate-spin" /> : <Dna className="w-3 h-3" />}
                          </button>
                        </TooltipTrigger>
                        <TooltipContent>Residue Info</TooltipContent>
                      </Tooltip>
                    </TooltipProvider>
                  )}
                  <TooltipProvider>
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <button onClick={() => handleFetchLinks(match.rule_id, match)}
                          className="p-1 rounded hover:bg-blue-50 transition-colors text-gray-400 hover:text-blue-600"
                          disabled={loadingLinks === cacheKey}>
                          {loadingLinks === cacheKey ? <Loader2 className="w-3 h-3 animate-spin" /> : <Link2 className="w-3 h-3" />}
                        </button>
                      </TooltipTrigger>
                      <TooltipContent>Database Links</TooltipContent>
                    </Tooltip>
                  </TooltipProvider>
                </div>
              </div>

              {match.enzyme && (
                <div className="flex flex-wrap gap-1 mb-1">
                  <Badge variant="secondary" className="bg-teal-50 text-teal-700 text-xs">
                    <FlaskConical className="w-2 h-2 mr-0.5" />{match.enzyme}
                  </Badge>
                  {match.category && <Badge variant="outline" className="text-xs text-gray-500">{match.category}</Badge>}
                </div>
              )}

              <div className="text-xs font-mono text-gray-500 mb-1.5 break-all">
                {match.reaction_smarts.length > 80 ? match.reaction_smarts.slice(0, 80) + "..." : match.reaction_smarts}
              </div>

              {match.products && match.products.length > 0 && (
                <div className="flex flex-wrap gap-1">
                  {match.products.slice(0, 3).map((p, pi) => (
                    <code key={pi} className="text-xs font-mono bg-gray-50 text-gray-600 px-1.5 py-0.5 rounded truncate max-w-[150px]">{p}</code>
                  ))}
                  {match.products.length > 3 && <span className="text-xs text-gray-400">+{match.products.length - 3}</span>}
                </div>
              )}

              {expandedResidue === cacheKey && residueInfo && (
                <div className="mt-2 p-2 rounded-md bg-amber-50/50 border border-amber-100 space-y-1.5">
                  <div className="flex items-center gap-1"><Dna className="w-3 h-3 text-amber-600" /><span className="text-xs font-medium text-amber-700">Residue Analysis</span></div>
                  {residueInfo.residues.map((res, ri) => (
                    <div key={ri} className="flex items-center gap-1.5">
                      <Badge variant="outline" className={`text-xs ${ROLE_COLORS[res.role] || "bg-gray-100 text-gray-600"}`}>{res.role}</Badge>
                      <span className="text-xs text-gray-600">{res.suggested_residue}</span>
                      <span className="text-xs text-gray-400 ml-auto">{(res.confidence * 100).toFixed(0)}%</span>
                    </div>
                  ))}
                  {residueInfo.residue_summary && (
                    <div className="flex gap-1 flex-wrap mt-1">
                      {residueInfo.residue_summary.catalytic_triad && <Badge variant="secondary" className="text-xs bg-amber-100 text-amber-800">Triad: {residueInfo.residue_summary.catalytic_triad.join("-")}</Badge>}
                      {residueInfo.residue_summary.metal_ion && <Badge variant="secondary" className="text-xs bg-violet-100 text-violet-800">Metal: {residueInfo.residue_summary.metal_ion}</Badge>}
                    </div>
                  )}
                </div>
              )}

              {linksInfo && (
                <div className="mt-2 p-2 rounded-md bg-blue-50/50 border border-blue-100">
                  <div className="flex items-center gap-1 mb-1"><Link2 className="w-3 h-3 text-blue-600" /><span className="text-xs font-medium text-blue-700">Cross References</span></div>
                  <div className="flex flex-wrap gap-1">
                    {Object.entries(linksInfo).filter(([k]) => !["integrity"].includes(k)).map(([dbKey, dbVal]) => {
                      if (!dbVal) return null;
                      const link = dbVal as DbLink;
                      if (!link.available) return null;
                      return (
                        <a key={dbKey} href={link.url} target="_blank" rel="noopener noreferrer"
                          className="text-xs text-blue-600 hover:text-blue-800 hover:underline flex items-center gap-0.5">
                          <ExternalLink className="w-2.5 h-2.5" />{link.label || dbKey}
                        </a>
                      );
                    })}
                  </div>
                </div>
              )}
            </div>
          );
        })}
        {filteredMatches.length > 50 && (
          <p className="text-center text-xs text-gray-400 py-2">Showing 50 of {filteredMatches.length} rules. Use filters to narrow results.</p>
        )}
      </div>
    </div>
  );
}
