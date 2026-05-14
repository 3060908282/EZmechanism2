"use client";

import React, { useRef } from "react";
import {
  Beaker, ArrowRight, Sparkles, AlertCircle, RotateCcw, Eye, Microscope,
  GitBranch, ExternalLink, Search, Loader2, FlaskRound, Layers, Dna, Zap, Network,
} from "lucide-react";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Slider } from "@/components/ui/slider";
import { Separator } from "@/components/ui/separator";
import { useToast } from "@/hooks/use-toast";
import { QUICK_SUBSTRATES } from "@/lib/constants";
import type { MolInfo, MechanismSearchResult, HealthInfo, RuleMatch } from "@/lib/types";
import MoleculeSVG from "@/components/mechanism/MoleculeSVG";
import SmilesDisplay from "@/components/mechanism/SmilesDisplay";
import MechanismPanel from "@/components/mechanism/MechanismPanel";
import MechanismGraphPanel from "@/components/mechanism/MechanismGraphPanel";
import MoleculesPanel from "@/components/mechanism/MoleculesPanel";
import RulesPanel from "@/components/mechanism/RulesPanel";
import StatsPanel from "@/components/mechanism/StatsPanel";
import EmptyState from "@/components/mechanism/EmptyState";
import LoadingSkeleton from "@/components/mechanism/LoadingSkeleton";

export default function QuickPredictMode({
  smiles, setSmiles, maxSteps, setMaxSteps, molInfo,
  isLoading, productSmiles, setProductSmiles, maxConfigs, setMaxConfigs,
  mechSearchResult, isMechSearching, selectedSmiles, setSelectedSmiles,
  rulesFilter, setRulesFilter, activeTab, setActiveTab, healthInfo,
  onPredict, onMechSearch, onReset, onExploreProduct, filteredMatches,
}: {
  smiles: string; setSmiles: (s: string) => void;
  maxSteps: number; setMaxSteps: (n: number) => void;
  molInfo: MolInfo | null;
  isLoading: boolean;
  productSmiles: string; setProductSmiles: (s: string) => void;
  maxConfigs: number; setMaxConfigs: (n: number) => void;
  mechSearchResult: MechanismSearchResult | null;
  isMechSearching: boolean;
  selectedSmiles: string | null; setSelectedSmiles: (s: string | null) => void;
  rulesFilter: "all" | "mcsa" | "builtin"; setRulesFilter: (f: "all" | "mcsa" | "builtin") => void;
  activeTab: string; setActiveTab: (t: string) => void;
  healthInfo: HealthInfo | null;
  onPredict: () => void; onMechSearch: () => void; onReset: () => void;
  onExploreProduct: (s: string) => void;
  filteredMatches: RuleMatch[];
}) {
  const { toast } = useToast();
  const debounceRef = useRef<ReturnType<typeof setTimeout>>(undefined);

  // Debounced molInfo is handled by parent
  const hasResults = mechSearchResult;

  return (
    <div className="grid grid-cols-1 lg:grid-cols-12 gap-5">
      {/* Left Panel */}
      <div className="lg:col-span-4 space-y-4">
        {/* Input Card */}
        <Card className="border-emerald-100 shadow-sm">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center gap-2 text-sm">
              <Beaker className="w-4 h-4 text-emerald-600" />
              Substrate Input
            </CardTitle>
            <CardDescription className="text-xs">Enter a SMILES string or select a substrate</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="space-y-1.5">
              <Label htmlFor="qp-smiles" className="text-xs font-medium">Substrate SMILES</Label>
              <div className="relative">
                <Input id="qp-smiles" placeholder="e.g., CCO for ethanol" value={smiles} onChange={(e) => setSmiles(e.target.value)}
                  onKeyDown={(e) => e.key === "Enter" && onPredict()}
                  className="pr-9 font-mono text-xs border-emerald-200 focus:border-emerald-400 focus:ring-emerald-400/20 h-9" />
                {smiles && <button onClick={() => setSmiles("")} className="absolute right-2.5 top-1/2 -translate-y-1/2 text-gray-400 hover:text-gray-600"><RotateCcw className="w-3 h-3" /></button>}
              </div>
              {molInfo && molInfo.valid === false && (
                <p className="text-xs text-red-500 flex items-center gap-1"><AlertCircle className="w-3 h-3" />{molInfo.error || "Invalid SMILES"}</p>
              )}
            </div>

            <div className="space-y-1.5">
              <Label className="text-xs text-gray-400 uppercase tracking-wider">Quick Fill</Label>
              <div className="grid grid-cols-2 gap-1.5">
                {QUICK_SUBSTRATES.map((sub) => (
                  <Button key={sub.name} variant="outline" size="sm"
                    className={`justify-start text-xs h-7 px-2 border-emerald-100 hover:bg-emerald-50 hover:border-emerald-300 hover:text-emerald-700 ${smiles === sub.smiles ? "bg-emerald-50 border-emerald-300 text-emerald-700" : ""}`}
                    onClick={() => setSmiles(sub.smiles)}>
                    <span className="mr-1 text-sm">{sub.icon}</span>
                    <div className="text-left min-w-0">
                      <div className="truncate font-medium">{sub.name}</div>
                      <div className="text-xs text-gray-400 truncate">{sub.desc}</div>
                    </div>
                  </Button>
                ))}
              </div>
            </div>

            {smiles && molInfo?.valid && (
              <div className="rounded-lg border border-emerald-100 bg-white p-3">
                <div className="flex items-center gap-1.5 mb-2">
                  <Eye className="w-3 h-3 text-emerald-600" />
                  <span className="text-xs font-medium text-emerald-700 uppercase tracking-wider">Structure Preview</span>
                </div>
                <div className="flex justify-center">
                  <MoleculeSVG smiles={smiles} size={180} zoomable />
                </div>
              </div>
            )}

            {smiles && molInfo?.valid && (
              <div className="rounded-lg border border-gray-100 bg-gray-50/50 p-3">
                <div className="flex items-center gap-1.5 mb-2">
                  <Microscope className="w-3 h-3 text-gray-500" />
                  <span className="text-xs font-medium text-gray-600 uppercase tracking-wider">Molecular Properties</span>
                </div>
                <div className="grid grid-cols-2 gap-x-4 gap-y-1.5 text-xs">
                  {[
                    ["Formula", molInfo.formula], ["MW", molInfo.molecular_weight?.toFixed(2)],
                    ["Atoms", molInfo.num_atoms], ["Heavy", molInfo.num_heavy_atoms],
                    ["HBD", molInfo.num_h_bond_donors], ["HBA", molInfo.num_h_bond_acceptors],
                    ["LogP", molInfo.logp?.toFixed(2)], ["TPSA", molInfo.tpsa?.toFixed(1)],
                    ["Rings", molInfo.ring_count], ["Rot.B.", molInfo.num_rotatable_bonds],
                  ].map(([label, value]) => (
                    <div key={label} className="flex justify-between">
                      <span className="text-gray-500">{label}</span>
                      <span className="font-mono font-medium">{value}</span>
                    </div>
                  ))}
                </div>
              </div>
            )}

            <Separator />

            <div className="space-y-2">
              <div className="flex items-center justify-between">
                <Label className="text-xs font-medium">Max Prediction Steps</Label>
                <Badge variant="secondary" className="bg-emerald-50 text-emerald-700 font-mono text-xs">{maxSteps}</Badge>
              </div>
              <Slider value={[maxSteps]} onValueChange={(v) => setMaxSteps(v[0])} min={1} max={10} step={1}
                className="[&_[data-slot=slider-range]]:bg-emerald-500 [&_[data-slot=slider-thumb]]:border-emerald-500" />
            </div>

            <div className="flex gap-2">
              <Button onClick={onPredict} disabled={isLoading || !smiles.trim()} className="flex-1 bg-emerald-600 hover:bg-emerald-700 text-white shadow-sm h-9">
                {isLoading ? <Loader2 className="w-3.5 h-3.5 mr-1.5 animate-spin" /> : <Sparkles className="w-3.5 h-3.5 mr-1.5" />}
                {isLoading ? "Predicting..." : "Predict"}
              </Button>
              <Button variant="outline" onClick={onReset} className="border-emerald-200 hover:bg-emerald-50 h-9 w-9 p-0"><RotateCcw className="w-3.5 h-3.5" /></Button>
            </div>
          </CardContent>
        </Card>

        {/* Bidirectional Search Card */}
        <Card className="border-violet-100 shadow-sm">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center gap-2 text-sm">
              <GitBranch className="w-4 h-4 text-violet-600" />
              Bidirectional Search
            </CardTitle>
            <CardDescription className="text-xs">Find connecting mechanism paths between reactant and product</CardDescription>
          </CardHeader>
          <CardContent className="space-y-3">
            <div className="space-y-1.5">
              <Label className="text-xs font-medium">Product SMILES</Label>
              <Input placeholder="e.g., CC(=O)O for acetic acid" value={productSmiles} onChange={(e) => setProductSmiles(e.target.value)}
                className="font-mono text-xs border-violet-200 focus:border-violet-400 focus:ring-violet-400/20 h-9" />
            </div>
            <div className="space-y-2">
              <div className="flex items-center justify-between">
                <Label className="text-xs font-medium">Max Configs</Label>
                <Badge variant="secondary" className="bg-violet-50 text-violet-700 font-mono text-xs">{maxConfigs}</Badge>
              </div>
              <Slider value={[maxConfigs]} onValueChange={(v) => setMaxConfigs(v[0])} min={50} max={500} step={50}
                className="[&_[data-slot=slider-range]]:bg-violet-500 [&_[data-slot=slider-thumb]]:border-violet-500" />
            </div>
            {productSmiles && (
              <div className="flex justify-center">
                <MoleculeSVG smiles={productSmiles} size={120} zoomable />
              </div>
            )}
            <Button onClick={onMechSearch} disabled={isMechSearching || !smiles.trim() || !productSmiles.trim()}
              className="w-full bg-violet-600 hover:bg-violet-700 text-white shadow-sm h-9">
              {isMechSearching ? <Loader2 className="w-3.5 h-3.5 mr-1.5 animate-spin" /> : <GitBranch className="w-3.5 h-3.5 mr-1.5" />}
              {isMechSearching ? "Searching..." : "Bi-Search"}
            </Button>
            {mechSearchResult && (
              <div className="rounded-lg bg-violet-50/50 border border-violet-100 p-3 space-y-2 text-xs">
                <div className="flex justify-between"><span className="text-gray-500">Paths Found</span><span className="font-mono font-medium">{mechSearchResult.paths?.length ?? 0}</span></div>
                <div className="flex justify-between"><span className="text-gray-500">Time</span><span className="font-mono font-medium">{(mechSearchResult.total_time || mechSearchResult.search_time || 0).toFixed(1)}s</span></div>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Product Explorer */}
        {selectedSmiles && (
          <Card className="border-amber-100 shadow-sm bg-amber-50/30">
            <CardHeader className="pb-2 pt-3 px-4">
              <div className="flex items-center justify-between">
                <CardTitle className="flex items-center gap-2 text-xs text-amber-800">
                  <ExternalLink className="w-3.5 h-3.5" />Exploring Product
                </CardTitle>
                <Button variant="ghost" size="sm" className="h-6 w-6 p-0 text-amber-600 hover:text-amber-800 hover:bg-amber-100" onClick={() => setSelectedSmiles(null)}>
                  <RotateCcw className="w-3 h-3" />
                </Button>
              </div>
            </CardHeader>
            <CardContent className="px-4 pb-3">
              <div className="flex gap-3">
                <MoleculeSVG smiles={selectedSmiles} size={80} zoomable />
                <div className="flex-1 min-w-0">
                  <SmilesDisplay smiles={selectedSmiles} maxLen={25} />
                  <Button variant="outline" size="sm" className="mt-2 h-7 text-xs border-amber-200 hover:bg-amber-100 hover:border-amber-300 text-amber-700"
                    onClick={() => { setSmiles(selectedSmiles); setSelectedSmiles(null); }}>
                    <ArrowRight className="w-3 h-3 mr-1" />Use as New Substrate
                  </Button>
                </div>
              </div>
            </CardContent>
          </Card>
        )}
      </div>

      {/* Right Panel: Results */}
      <div className="lg:col-span-8">
        {isLoading || isMechSearching ? (
          <LoadingSkeleton />
        ) : hasResults ? (
          <Tabs value={activeTab} onValueChange={setActiveTab}>
            <Card className="border-emerald-100 shadow-sm">
              <CardHeader className="pb-2">
                <div className="flex items-center justify-between flex-wrap gap-2">
                  <CardTitle className="flex items-center gap-2 text-sm">
                    <Search className="w-4 h-4 text-emerald-600" />
                    Results
                  </CardTitle>
                  <TabsList className="bg-emerald-50 p-0.5 h-8">
                    <TabsTrigger value="mechanism" className="data-[state=active]:bg-white data-[state=active]:text-emerald-700 data-[state=active]:shadow-sm text-xs h-7 px-2.5">
                      <FlaskRound className="w-3 h-3 mr-1" />Mechanism
                    </TabsTrigger>
                    <TabsTrigger value="graph" className="data-[state=active]:bg-white data-[state=active]:text-violet-700 data-[state=active]:shadow-sm text-xs h-7 px-2.5">
                      <Network className="w-3 h-3 mr-1" />Graph
                    </TabsTrigger>
                    <TabsTrigger value="molecules" className="data-[state=active]:bg-white data-[state=active]:text-emerald-700 data-[state=active]:shadow-sm text-xs h-7 px-2.5">
                      <Layers className="w-3 h-3 mr-1" />Molecules
                    </TabsTrigger>
                    <TabsTrigger value="rules" className="data-[state=active]:bg-white data-[state=active]:text-emerald-700 data-[state=active]:shadow-sm text-xs h-7 px-2.5">
                      <Dna className="w-3 h-3 mr-1" />Rules
                    </TabsTrigger>
                    <TabsTrigger value="stats" className="data-[state=active]:bg-white data-[state=active]:text-emerald-700 data-[state=active]:shadow-sm text-xs h-7 px-2.5">
                      <Zap className="w-3 h-3 mr-1" />Stats
                    </TabsTrigger>
                  </TabsList>
                </div>
              </CardHeader>
              <CardContent>
                <TabsContent value="mechanism"><MechanismPanel searchResult={mechSearchResult} onExploreProduct={onExploreProduct} /></TabsContent>
                <TabsContent value="graph"><MechanismGraphPanel result={mechSearchResult} /></TabsContent>
                <TabsContent value="molecules"><MoleculesPanel searchResult={mechSearchResult} onExploreProduct={onExploreProduct} /></TabsContent>
                <TabsContent value="rules"><RulesPanel matchResult={null} filter={rulesFilter} onFilterChange={setRulesFilter} filteredMatches={filteredMatches} /></TabsContent>
                <TabsContent value="stats"><StatsPanel searchResult={mechSearchResult} /></TabsContent>
              </CardContent>
            </Card>
          </Tabs>
        ) : (
          <EmptyState mode="quick" />
        )}
      </div>
    </div>
  );
}
