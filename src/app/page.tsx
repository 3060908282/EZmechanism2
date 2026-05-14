"use client";

import React, { useState, useEffect } from "react";
import {
  Hexagon, Database, Activity,
} from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { API_BASE } from "@/lib/constants";
import type {
  HealthInfo,
  MechanismSearchResult, PdbInfo, SelectedResidue,
  SubstrateCofactor,
} from "@/lib/types";
import PdbWorkflowMode from "@/components/mechanism/PdbWorkflowMode";
import DebugDialog from "@/components/mechanism/DebugDialog";



export default function Home() {
  // PDB Workflow persistent state
  const [wfStep, setWfStep] = useState(1);
  const [wfPdbId, setWfPdbId] = useState("");
  const [wfPdbText, setWfPdbText] = useState("");
  const [wfPdbRawText, setWfPdbRawText] = useState("");
  const [wfPdbLoading, setWfPdbLoading] = useState(false);
  const [wfPdbInfo, setWfPdbInfo] = useState<PdbInfo | null>(null);
  const [wfPdbEntries, setWfPdbEntries] = useState<PdbInfo[]>([]);
  const [wfSelectedChain, setWfSelectedChain] = useState<string | null>("");
  const [wfResidueFilter, setWfResidueFilter] = useState("");
  const [wfSelectedResidues, setWfSelectedResidues] = useState<SelectedResidue[]>([]);
  const [wfSubstrates, setWfSubstrates] = useState<SubstrateCofactor[]>([]);
  const [wfReactantSmiles, setWfReactantSmiles] = useState("");
  const [wfProductSmiles, setWfProductSmiles] = useState("");
  const [wfMaxConfigs, setWfMaxConfigs] = useState(300);
  const [wfMaxBondLength, setWfMaxBondLength] = useState(11); // Max Bond Length (0-11 Å)
  const [wfComment, setWfComment] = useState("");
  const [wfRunning, setWfRunning] = useState(false);
  const [wfResult, setWfResult] = useState<MechanismSearchResult | null>(null);
  const [wfActiveTab, setWfActiveTab] = useState("mechanism");
  const [wfRulesFilter, setWfRulesFilter] = useState<"all" | "mcsa" | "builtin">("all");
  const [wfMcsaId, setWfMcsaId] = useState(0);
  const [wfSelectedSmiles, setWfSelectedSmiles] = useState<string | null>(null);
  const [wfTestId, setWfTestId] = useState<string | null>(null);
  const [healthInfo, setHealthInfo] = useState<HealthInfo | null>(null);

  // Fetch health info on mount
  useEffect(() => {
    fetch(API_BASE)
      .then((r) => r.json())
      .then((d) => setHealthInfo({ status: d.status, uptime: "", rules: { total_rules: 51647, mcsa_rules_test: 0, mcsa_rules_natmet: 51637, builtin_rules: 10 } }))
      .catch(() => { /* ignore */ });
  }, []);

  return (
    <div className="min-h-screen flex flex-col">
      {/* Background */}
      <div className="fixed inset-0 -z-10 bg-gradient-to-br from-emerald-50/80 via-white to-teal-50/60" />
      <div
        className="fixed inset-0 -z-10 opacity-[0.03]"
        style={{
          backgroundImage: `url("data:image/svg+xml,%3Csvg width='60' height='60' viewBox='0 0 60 60' xmlns='http://www.w3.org/2000/svg'%3E%3Cg fill='none' fill-rule='evenodd'%3E%3Cg fill='%2310b981' fill-opacity='1'%3E%3Cpath d='M36 34v-4h-2v4h-4v2h4v4h2v-4h4v-2h-4zm0-30V0h-2v4h-4v2h4v4h2V6h4V4h-4zM6 34v-4H4v4H0v2h4v4h2v-4h4v-2H6zM6 4V0H4v4H0v2h4v4h2V6h4V4H6z'/%3E%3C/g%3E%3C/g%3E%3C/svg%3E")`,
        }}
      />

      {/* Header */}
      <header className="border-b border-emerald-100 bg-white/80 backdrop-blur-sm sticky top-0 z-50">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8">
          <div className="flex items-center justify-between h-14">
            <div className="flex items-center gap-3">
              <div className="w-9 h-9 rounded-xl bg-gradient-to-br from-emerald-500 to-teal-600 flex items-center justify-center shadow-md">
                <Hexagon className="w-4 h-4 text-white" />
              </div>
              <div>
                <h1 className="text-base font-bold tracking-tight text-gray-900">M-CSA Mechanism Predictor</h1>
                <p className="text-[10px] text-gray-400 -mt-0.5">Enzyme Reaction Mechanism Prediction</p>
              </div>
            </div>
            <div className="flex items-center gap-2">
              {healthInfo && (
                <Badge variant="secondary" className="bg-emerald-50 text-emerald-700 border-emerald-200 text-xs">
                  <Database className="w-3 h-3 mr-1" />{healthInfo.rules.total_rules.toLocaleString()} Rules
                </Badge>
              )}
              <Badge variant="outline" className="text-gray-500 text-xs hidden sm:flex">
                <Activity className="w-3 h-3 mr-1" />{healthInfo?.uptime || "Online"}
              </Badge>
              <DebugDialog
                selectedResidues={wfSelectedResidues}
                reactionSmiles={wfReactantSmiles && wfProductSmiles ? `${wfReactantSmiles}>>${wfProductSmiles}` : wfReactantSmiles ? `${wfReactantSmiles}>>` : ""}
                searchResult={wfResult}
              />
            </div>
          </div>
        </div>
      </header>

      {/* Main Content */}
      <main className="flex-1 max-w-7xl mx-auto w-full px-4 sm:px-6 lg:px-8 py-5">
        {/* Hero Banner */}
        <div className="mb-5 rounded-2xl overflow-hidden relative h-36 sm:h-44">
          <img src="/hero-banner.png" alt="Enzyme Mechanism Prediction" className="w-full h-full object-cover" />
          <div className="absolute inset-0 bg-gradient-to-r from-emerald-900/85 via-emerald-900/60 to-teal-900/30 flex items-center">
            <div className="px-6 sm:px-10 max-w-xl">
              <h2 className="text-xl sm:text-2xl font-bold text-white mb-1.5">EzMechanism Workflow</h2>
              <p className="text-emerald-100/90 text-xs sm:text-sm leading-relaxed">
                Predict enzyme reaction mechanisms using M-CSA curated data with {healthInfo?.rules.total_rules.toLocaleString() || "50,000+"} reaction rules.
                Structure-based enzyme reaction mechanism prediction.
              </p>
            </div>
          </div>
        </div>

        <PdbWorkflowMode
          healthInfo={healthInfo}
          currentStep={wfStep} setCurrentStep={setWfStep}
          pdbId={wfPdbId} setPdbId={setWfPdbId}
          pdbText={wfPdbText} setPdbText={setWfPdbText}
          pdbRawText={wfPdbRawText} setPdbRawText={setWfPdbRawText}
          pdbLoading={wfPdbLoading} setPdbLoading={setWfPdbLoading}
          pdbInfo={wfPdbInfo} setPdbInfo={setWfPdbInfo}
          pdbEntries={wfPdbEntries} setPdbEntries={setWfPdbEntries}
          selectedChain={wfSelectedChain} setSelectedChain={setWfSelectedChain}
          residueFilter={wfResidueFilter} setResidueFilter={setWfResidueFilter}
          selectedResidues={wfSelectedResidues} setSelectedResidues={setWfSelectedResidues}
          substrates={wfSubstrates} setSubstrates={setWfSubstrates}
          reactantSmiles={wfReactantSmiles} setReactantSmiles={setWfReactantSmiles}
          productSmilesWf={wfProductSmiles} setProductSmilesWf={setWfProductSmiles}
          wfMaxConfigs={wfMaxConfigs} setWfMaxConfigs={setWfMaxConfigs}
          wfMaxBondLength={wfMaxBondLength} setWfMaxBondLength={setWfMaxBondLength}
          wfComment={wfComment} setWfComment={setWfComment}
          wfRunning={wfRunning} setWfRunning={setWfRunning}
          wfResult={wfResult} setWfResult={setWfResult}
          wfActiveTab={wfActiveTab} setWfActiveTab={setWfActiveTab}
          wfRulesFilter={wfRulesFilter} setWfRulesFilter={setWfRulesFilter}
          wfMcsaId={wfMcsaId} setWfMcsaId={setWfMcsaId}
          wfSelectedSmiles={wfSelectedSmiles} setWfSelectedSmiles={setWfSelectedSmiles}
          testId={wfTestId} onTestIdChange={setWfTestId}
        />
      </main>

      {/* Footer */}
      <footer className="border-t border-emerald-100 bg-white/80 backdrop-blur-sm mt-auto sticky bottom-0">
        <div className="max-w-7xl mx-auto px-4 sm:px-6 lg:px-8 py-3">
          <div className="flex items-center justify-between text-[11px] text-gray-400">
            <p>M-CSA Mechanism Predictor — EzMechanism</p>
            <p>
              Powered by <span className="text-emerald-600 font-medium">RDKit</span> &{" "}
              <span className="text-emerald-600 font-medium">M-CSA Database</span>
            </p>
          </div>
        </div>
      </footer>
    </div>
  );
}
