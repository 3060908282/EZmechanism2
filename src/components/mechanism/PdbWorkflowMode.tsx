"use client";

import React, { useState, useEffect, useCallback, useRef, useMemo } from "react";
import {
  Dna, Sparkles, FileText, ChevronRight, ChevronLeft, Download, Beaker,
  ArrowRight, CheckCircle2, AlertCircle, Trash2, Pencil, Plus, Loader2, X,
  FlaskRound, FlaskConical, Layers, Network, Zap, Activity, GitBranch,
  Search, Pencil as PencilLucide, Crosshair, FolderOpen, Eye, Copy, ZoomIn,
} from "lucide-react";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Textarea } from "@/components/ui/textarea";
import { Alert, AlertTitle, AlertDescription } from "@/components/ui/alert";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Slider } from "@/components/ui/slider";
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import {
  AlertDialog, AlertDialogAction, AlertDialogCancel, AlertDialogContent,
  AlertDialogDescription, AlertDialogFooter, AlertDialogHeader, AlertDialogTitle,
} from "@/components/ui/alert-dialog";
import {
  Dialog, DialogContent, DialogHeader, DialogTitle,
} from "@/components/ui/dialog";
import { useToast } from "@/hooks/use-toast";
import ChemicalEditor, { type ChemicalEditorHandle } from "@/components/ChemicalEditor";
import { apiCall, apiCallJson, genId, normalizeSearchResult } from "@/lib/api";
import { RESIDUE_SMILES_MAP, VALID_RESIDUE_NAMES } from "@/lib/constants";
import type {
  HealthInfo, PdbInfo, PdbResidue, SelectedResidue,
  SubstrateCofactor, MechanismSearchResult, SearchProgress,
} from "@/lib/types";
import StepWizard from "@/components/mechanism/StepWizard";
import MechanismPanel from "@/components/mechanism/MechanismPanel";
import MechanismGraphPanel from "@/components/mechanism/MechanismGraphPanel";
import MoleculesPanel from "@/components/mechanism/MoleculesPanel";
import RulesPanel from "@/components/mechanism/RulesPanel";
import StatsPanel from "@/components/mechanism/StatsPanel";
import MoleculeSVG from "@/components/mechanism/MoleculeSVG";
import EditMoleculeDialog, { type OverlayResult } from "@/components/mechanism/EditMoleculeDialog";
import PdbStructureViewer, { type PdbStructureViewerHandle, type OverlayMolecule } from "@/components/mechanism/PdbStructureViewer";
import ZoomModal from "@/components/mechanism/ZoomModal";

// Amino acid data for R-Group style rendering in reaction preview
// Side chain SMILES: represents only the side chain starting from Cβ
// atom 0 = attachment point (Cβ), labeled with residue name (e.g. "Cys145A")
// sideFormula: short formula for the side chain functional group (shown in label)
interface AminoAcidData {
  sideChain: string;   // Side chain SMILES (atom 0 = Cβ, the attachment point)
  sideFormula: string;  // Short formula for label, e.g. "OH" for Ser
}
const AMINO_ACID_SIDECHAIN: Record<string, AminoAcidData> = {
  ALA: { sideChain: "CC",            sideFormula: "CH₃" },    // Alanine: -CH(CH₃)
  ARG: { sideChain: "NCCC(N)",       sideFormula: "guanidino" }, // Arginine: long
  ASN: { sideChain: "CC(=O)N",       sideFormula: "CONH₂" },   // Asparagine
  ASP: { sideChain: "CC(=O)O",       sideFormula: "COOH" },    // Aspartic acid
  CYS: { sideChain: "CS",            sideFormula: "SH" },       // Cysteine
  GLN: { sideChain: "CCC(=O)N",      sideFormula: "CONH₂" },   // Glutamine
  GLU: { sideChain: "CCC(=O)O",      sideFormula: "COOH" },    // Glutamic acid
  GLY: { sideChain: "C",             sideFormula: "H" },        // Glycine (only H as side chain)
  HIS: { sideChain: "Cc1[nH]cnc1",   sideFormula: "imidazole" },// Histidine
  ILE: { sideChain: "CC(C)C",        sideFormula: "CH(CH₃)₂" },// Isoleucine
  LEU: { sideChain: "CC(C)C",        sideFormula: "CH₂CH(CH₃)₂" }, // Leucine
  LYS: { sideChain: "CCCCN",         sideFormula: "NH₃⁺" },    // Lysine
  MET: { sideChain: "CCSC",          sideFormula: "SCH₃" },    // Methionine
  PHE: { sideChain: "Cc1ccccc1",     sideFormula: "Ph" },      // Phenylalanine
  PRO: { sideChain: "C1CCC1",        sideFormula: "pyrrolidine" }, // Proline
  SER: { sideChain: "CO",            sideFormula: "OH" },       // Serine
  THR: { sideChain: "C(C)O",         sideFormula: "OH" },       // Threonine
  TRP: { sideChain: "Cc1c[nH]c2ccccc12", sideFormula: "indole" }, // Tryptophan
  TYR: { sideChain: "COc1ccccc1",    sideFormula: "OH" },       // Tyrosine
  VAL: { sideChain: "C(C)C",         sideFormula: "CH(CH₃)₂" },// Valine
};

// Directory test item type
interface PdbTestItem {
  id: string;
  testId: string;
  pdbId: string | null;
  state: string;
  comment: string;
  currentStep: number;
  createdAt: string;
  updatedAt: string;
}

export default function PdbWorkflowMode({ healthInfo,
  currentStep, setCurrentStep,
  pdbId, setPdbId, pdbText, setPdbText, pdbRawText, setPdbRawText, pdbLoading, setPdbLoading,
  pdbInfo, setPdbInfo, pdbEntries, setPdbEntries,
  selectedChain, setSelectedChain, residueFilter, setResidueFilter,
  selectedResidues, setSelectedResidues,
  substrates, setSubstrates, reactantSmiles, setReactantSmiles,
  productSmilesWf, setProductSmilesWf,
  wfMaxConfigs, setWfMaxConfigs, wfMaxBondLength, setWfMaxBondLength,
  wfComment, setWfComment, wfRunning, setWfRunning,
  wfResult, setWfResult, wfActiveTab, setWfActiveTab,
  wfRulesFilter, setWfRulesFilter, wfSelectedSmiles, setWfSelectedSmiles,
  wfMcsaId, setWfMcsaId,
  testId: testIdProp, onTestIdChange: onTestIdChangeProp,
}: {
  healthInfo: HealthInfo | null;
  currentStep: number; setCurrentStep: (s: number) => void;
  pdbId: string; setPdbId: (v: string) => void;
  pdbText: string; setPdbText: (v: string) => void;
  pdbRawText: string; setPdbRawText: (v: string) => void;
  pdbLoading: boolean; setPdbLoading: (v: boolean) => void;
  pdbInfo: PdbInfo | null; setPdbInfo: (v: PdbInfo | null) => void;
  pdbEntries: PdbInfo[]; setPdbEntries: React.Dispatch<React.SetStateAction<PdbInfo[]>>;
  selectedChain: string | null; setSelectedChain: (v: string | null) => void;
  residueFilter: string; setResidueFilter: (v: string) => void;
  selectedResidues: SelectedResidue[]; setSelectedResidues: React.Dispatch<React.SetStateAction<SelectedResidue[]>>;
  substrates: SubstrateCofactor[]; setSubstrates: React.Dispatch<React.SetStateAction<SubstrateCofactor[]>>;
  reactantSmiles: string; setReactantSmiles: (v: string) => void;
  productSmilesWf: string; setProductSmilesWf: (v: string) => void;
  wfMaxConfigs: number; setWfMaxConfigs: (v: number) => void;
  wfMaxBondLength: number; setWfMaxBondLength: (v: number) => void;
  wfComment: string; setWfComment: (v: string) => void;
  wfRunning: boolean; setWfRunning: (v: boolean) => void;
  wfResult: MechanismSearchResult | null; setWfResult: (v: MechanismSearchResult | null) => void;
  wfActiveTab: string; setWfActiveTab: (v: string) => void;

  wfRulesFilter: "all" | "mcsa" | "builtin"; setWfRulesFilter: (v: "all" | "mcsa" | "builtin") => void;
  wfMcsaId: number; setWfMcsaId: (v: number) => void;
  wfSelectedSmiles: string | null; setWfSelectedSmiles: (v: string | null) => void;
  testId?: string | null;
  onTestIdChange?: (id: string | null) => void;
}) {
  const { toast } = useToast();
  const [zoomModalOpen, setZoomModalOpen] = useState(false);
  // Async search state
  const [wfJobId, setWfJobId] = useState<string | null>(null);
  const [wfSearchProgress, setWfSearchProgress] = useState<SearchProgress | null>(null);
  const [wfShowProgress, setWfShowProgress] = useState(false); // Controls whether to show the Progress/Results view
  const [zoomModalData, setZoomModalData] = useState<string>("");
  const [zoomModalLabel, setZoomModalLabel] = useState<string>("");
  const [zoomModalAtomLabel, setZoomModalAtomLabel] = useState<string>("");
  const [rxnPreviewZoom, setRxnPreviewZoom] = useState(false);
  const [automapping, setAutomapping] = useState(false);
  const [verifying, setVerifying] = useState(false);
  const [ligandDistances, setLigandDistances] = useState<Record<string, number>>({});
  const [distancesLoading, setDistancesLoading] = useState(false);
  const viewerStep2Ref = useRef<PdbStructureViewerHandle>(null);
  const reactionEditorRef = useRef<ChemicalEditorHandle>(null);

  // ---- Overlay State ----
  const [overlayResults, setOverlayResults] = useState<OverlayResult[]>([]);
  const OVERLAY_COLORS = ["#FF66FF", "#66FFFF", "#FFAA00", "#66FF66", "#FF6666", "#6688FF", "#FFFF66"];

  // Convert overlay results to viewer-compatible format
  const overlayMolecules = useMemo(() => {
    return overlayResults.map((r, i) => ({
      substrateId: r.substrateId,
      molblock: r.molblock,
      smiles: r.smiles,
      rmsd: r.rmsd,
      resName: r.resName,
      chain: r.chain,
      resNum: r.resNum,
      color: OVERLAY_COLORS[i % OVERLAY_COLORS.length],
    }));
  }, [overlayResults]);

  const handleOverlayResult = useCallback((result: OverlayResult) => {
    setOverlayResults(prev => {
      // Replace existing overlay for same substrate, or add new
      const filtered = prev.filter(r => r.substrateId !== result.substrateId);
      return [...filtered, result];
    });
  }, []);

  // ---- Test Directory State ----
  const [internalTestId, setInternalTestId] = useState<string | null>(null);
  const activeTestId = testIdProp ?? internalTestId;
  const setActiveTestId = onTestIdChangeProp ?? setInternalTestId;

  const [testsList, setTestsList] = useState<PdbTestItem[]>([]);
  const [testsLoading, setTestsLoading] = useState(false);
  const [deleteTarget, setDeleteTarget] = useState<string | null>(null);
  const [editingCommentId, setEditingCommentId] = useState<string | null>(null);
  const [editingCommentText, setEditingCommentText] = useState("");
  const [currentTestIdLabel, setCurrentTestIdLabel] = useState<string | null>(null);

  // Ref to always have the latest saveable state snapshot
  const stateRef = useRef({
    currentStep, pdbId, pdbRawText, pdbInfo, pdbEntries, selectedChain,
    selectedResidues, substrates, reactantSmiles,
    productSmilesWf, wfMaxConfigs, wfMaxBondLength, wfComment,
    wfResult,
  });
  useEffect(() => {
    stateRef.current = {
      currentStep, pdbId, pdbRawText, pdbInfo, pdbEntries, selectedChain,
      selectedResidues, substrates, reactantSmiles,
      productSmilesWf, wfMaxConfigs, wfMaxBondLength, wfComment,
      wfResult,
    };
  });

  const hasResults = wfResult;
  const showProgressOrResults = wfShowProgress || hasResults;

  // ---- Fetch Tests for Directory ----
  const fetchTests = useCallback(async () => {
    setTestsLoading(true);
    try {
      const r = await fetch("/api/pdb-tests");
      if (r.ok) {
        const data = await r.json();
        setTestsList(data);
      }
    } catch { /* silent */ }
    setTestsLoading(false);
  }, []);

  useEffect(() => {
    if (!activeTestId) {
      fetchTests();
    }
  }, [activeTestId, fetchTests]);

  // ---- Auto-refresh Directory when any test is in "run" state ----
  useEffect(() => {
    if (activeTestId) return; // Only refresh when viewing Directory
    const hasRunning = testsList.some(t => t.state === "run");
    if (!hasRunning) return;
    const interval = setInterval(() => {
      fetchTests();
    }, 10000); // Refresh every 10s
    return () => clearInterval(interval);
  }, [activeTestId, testsList, fetchTests]);

  // ---- Build save payload from state ref ----
  const buildSavePayload = useCallback((step: number, extraData?: Record<string, unknown>) => {
    const s = stateRef.current;
    return {
      currentStep: step,
      pdbId: s.pdbInfo?.pdb_id || s.pdbId || null,
      pdbRawText: s.pdbRawText || null,
      pdbInfoJson: s.pdbInfo ? JSON.stringify(s.pdbInfo) : null,
      pdbEntriesJson: JSON.stringify(s.pdbEntries),
      selectedChain: s.selectedChain,
      selectedResiduesJson: JSON.stringify(s.selectedResidues),
      substratesJson: JSON.stringify(s.substrates),
      reactantSmiles: s.reactantSmiles,
      productSmiles: s.productSmilesWf,
      maxConfigs: s.wfMaxConfigs,
      maxRules: 0,  // 0 = load all rules
      maxBondLength: s.wfMaxBondLength,
      comment: s.wfComment,
      state: "edit",
      ...extraData,
    };
  }, []);

  // ---- Auto-save ----
  const handleAutoSave = useCallback(async (step: number, extraData?: Record<string, unknown>) => {
    if (!activeTestId) return;
    try {
      await fetch(`/api/pdb-tests/${activeTestId}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify(buildSavePayload(step, extraData)),
      });
      toast({ title: "Progress saved", description: `Step ${step} saved successfully` });
    } catch { /* silent fail — fire-and-forget */ }
  }, [activeTestId, toast, buildSavePayload]);

  // ---- Before-unload: save via sendBeacon to prevent data loss on browser close ----
  useEffect(() => {
    if (!activeTestId) return;
    const handler = () => {
      const payload = buildSavePayload(currentStep);
      const blob = new Blob([JSON.stringify(payload)], { type: "application/json" });
      navigator.sendBeacon(`/api/pdb-tests/${activeTestId}`, blob);
    };
    window.addEventListener("beforeunload", handler);
    return () => window.removeEventListener("beforeunload", handler);
  }, [activeTestId, currentStep, buildSavePayload]);

  // ---- Go-to-step with auto-save ----
  const goToStep = useCallback((targetStep: number) => {
    if (activeTestId) {
      handleAutoSave(currentStep); // fire-and-forget
    }
    setCurrentStep(targetStep);
  }, [activeTestId, currentStep, handleAutoSave, setCurrentStep]);

  // ---- Go to Directory with auto-save ----
  const goToDirectory = useCallback(() => {
    if (activeTestId) {
      handleAutoSave(currentStep); // save current step before leaving
    }
    setActiveTestId(null);
    setWfShowProgress(false);
  }, [activeTestId, currentStep, handleAutoSave]);

  // ---- Create New Test ----
  const handleCreateNew = useCallback(async () => {
    try {
      const r = await fetch("/api/pdb-tests", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({}),
      });
      if (!r.ok) throw new Error("Failed to create test");
      const test = await r.json();
      setActiveTestId(test.id);
      setCurrentTestIdLabel(test.testId);
      setCurrentStep(1);
      // Clear previous workflow state
      setPdbInfo(null);
      setPdbRawText("");
      setPdbEntries([]);
      setSelectedChain(null);
      setSelectedResidues([]);
      setSubstrates([]);
      setReactantSmiles("");
      setProductSmilesWf("");
      setWfComment("");
      setWfResult(null);
      setOverlayResults([]);
      toast({ title: "Test Created", description: `${test.testId} — start from Step 1` });
    } catch {
      toast({ title: "Create Failed", description: "Could not create new test", variant: "destructive" });
    }
  }, [toast, setActiveTestId, setCurrentStep, setPdbInfo, setPdbRawText, setPdbEntries, setSelectedChain, setSelectedResidues, setSubstrates, setReactantSmiles, setProductSmilesWf, setWfComment, setWfResult]);

  // ---- Resume polling for a running search (used when viewing a test in "run" state) ----
  const resumePolling = useCallback(async (jobId: string, testId: string) => {
    const POLL_INTERVAL = 3000;
    const MAX_POLL_TIME = 30 * 60 * 1000; // 30 min max
    const pollStart = Date.now();
    let localSearchResult: Record<string, any> | null = null;

    while (Date.now() - pollStart < MAX_POLL_TIME) {
      await new Promise(r => setTimeout(r, POLL_INTERVAL));
      try {
        const status = await apiCallJson<SearchProgress>(
          "mechanism-search-status", { jobId }, 15000
        );
        setWfSearchProgress(status);

        if (status.state === "DONE") {
          if (status.result) {
            const normalized = normalizeSearchResult(status.result as unknown as Record<string, unknown>);
            localSearchResult = normalized as unknown as MechanismSearchResult;
            setWfResult(normalized as unknown as MechanismSearchResult);
            setWfShowProgress(true);
            if (Array.isArray(normalized.paths) && normalized.paths.length > 0) {
              setWfActiveTab("graph");
            }
          }
          break;
        } else if (status.state === "ERROR") {
          const errorMsg = status.error || "Search failed with an unknown error";
          localSearchResult = { error: errorMsg } as unknown as MechanismSearchResult;
          setWfResult(localSearchResult as unknown as MechanismSearchResult);
          setWfShowProgress(true);
          break;
        } else if (status.state === "CANCELLED") {
          localSearchResult = { error: "Search was cancelled" } as unknown as MechanismSearchResult;
          setWfResult(localSearchResult as unknown as MechanismSearchResult);
          setWfShowProgress(true);
          break;
        }
      } catch (pollErr) {
        console.warn("Resume-poll error (will retry):", pollErr);
      }
    }

    // Timeout or completion — save results to DB
    if (!localSearchResult && Date.now() - pollStart >= MAX_POLL_TIME) {
      localSearchResult = { error: "Search timed out (30 minutes). The search is still running on the server." } as unknown as MechanismSearchResult;
      setWfResult(localSearchResult as unknown as MechanismSearchResult);
      setWfShowProgress(true);
    }

    const searchRes = localSearchResult as Record<string, any> | null;
    const hasAnyResult = (searchRes?.paths as any[])?.length > 0;
    const finalState = hasAnyResult ? "done" : "fail";

    if (testId) {
      fetch(`/api/pdb-tests/${testId}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          ...buildSavePayload(5, { state: finalState }),
          searchResultJson: localSearchResult ? JSON.stringify(localSearchResult) : null,
        }),
      }).then(() => {
        fetchTests(); // Refresh directory list after DB save
      }).catch(() => { /* silent */ });
    }

    if (finalState === "done") {
      toast({ title: "Search Complete", description: "Results ready" });
    } else if (searchRes?.error) {
      toast({ title: "Search Error", description: String(searchRes.error).slice(0, 150), variant: "destructive" });
    } else {
      toast({ title: "No Results Found", description: "No matching mechanisms or products were found", variant: "destructive" });
    }

    setWfRunning(false);
    setWfSearchProgress(null);
  }, [buildSavePayload, toast, setWfActiveTab, fetchTests]);

  // ---- View (Restore) Test ----
  const handleViewTest = useCallback(async (dbId: string) => {
    try {
      const r = await fetch(`/api/pdb-tests/${dbId}`);
      if (!r.ok) throw new Error("Failed to load test");
      const test = await r.json();

      // Restore all state from DB record
      if (test.pdbId) setPdbId(test.pdbId); else setPdbId("");
      if (test.pdbRawText) setPdbRawText(test.pdbRawText); else setPdbRawText("");
      if (test.pdbInfoJson) {
        try { setPdbInfo(JSON.parse(test.pdbInfoJson)); } catch { setPdbInfo(null); }
      } else { setPdbInfo(null); }
      if (test.pdbEntriesJson) {
        try { setPdbEntries(JSON.parse(test.pdbEntriesJson)); } catch { setPdbEntries([]); }
      } else { setPdbEntries([]); }
      setSelectedChain(test.selectedChain || null);
      if (test.selectedResiduesJson) {
        try { setSelectedResidues(JSON.parse(test.selectedResiduesJson)); } catch { setSelectedResidues([]); }
      } else { setSelectedResidues([]); }
      if (test.substratesJson) {
        try { setSubstrates(JSON.parse(test.substratesJson)); } catch { setSubstrates([]); }
      } else { setSubstrates([]); }
      setReactantSmiles(test.reactantSmiles || "");
      setProductSmilesWf(test.productSmiles || "");
      setWfMaxConfigs(test.maxConfigs || 100);
      setWfMaxBondLength(test.maxBondLength ?? 11);
      setWfComment(test.comment || "");

      setActiveTestId(dbId);
      setCurrentTestIdLabel(test.testId);
      setCurrentStep(test.currentStep || 1);

      // If test was completed, restore results
      if (test.state === "done" || test.state === "fail") {
        if (test.searchResultJson) {
          try {
            const parsed = JSON.parse(test.searchResultJson);
            setWfResult(normalizeSearchResult(parsed) as unknown as MechanismSearchResult);
            setWfShowProgress(true);
          } catch { /* ignore */ }
        } else if (test.state === "fail") {
          // fail test with no searchResultJson — show a synthetic error
          setWfResult({ error: "No result data available for this failed test." } as unknown as MechanismSearchResult);
          setWfShowProgress(true);
        } else {
          setWfResult(null);
        }
      } else if (test.state === "run") {
        // Test is currently running — resume polling if we have a jobId
        setWfResult(null);
        setWfShowProgress(true);
        const savedJobId = test.jobId || null;
        if (savedJobId) {
          setWfJobId(savedJobId);
          setWfRunning(true);
          setWfSearchProgress({
            state: "STARTING",
            explored_nodes: 0,
            total_nodes: 0,
            current_iteration: 0,
            elapsed_seconds: 0,
          });
          // Start resume-polling in the background (fire-and-forget)
          resumePolling(savedJobId, dbId);
        } else {
          // No jobId — search was started but we can't resume polling
          // Show as if searching but with a note
          setWfRunning(true);
          setWfSearchProgress({
            state: "STARTING",
            explored_nodes: 0,
            total_nodes: 0,
            current_iteration: 0,
            elapsed_seconds: 0,
          });
        }
      } else {
        setWfResult(null);
      }

      toast({ title: "Test Loaded", description: `${test.testId} — Step ${test.currentStep}` });
    } catch {
      toast({ title: "Load Failed", description: "Could not load test data", variant: "destructive" });
    }
  }, [toast, setActiveTestId, setCurrentStep, setPdbId, setPdbRawText, setPdbInfo, setPdbEntries, setSelectedChain, setSelectedResidues, setSubstrates, setReactantSmiles, setProductSmilesWf, setWfMaxConfigs, setWfMaxBondLength, setWfComment, setWfResult, resumePolling]);

  // ---- Delete Test (with confirmation) ----
  const confirmDeleteTest = useCallback((dbId: string) => {
    setDeleteTarget(dbId);
  }, []);

  const handleDeleteTest = useCallback(async () => {
    if (!deleteTarget) return;
    const id = deleteTarget;
    setDeleteTarget(null);
    try {
      const r = await fetch(`/api/pdb-tests/${id}`, { method: "DELETE" });
      if (!r.ok) throw new Error("Failed to delete");
      setTestsList((prev) => prev.filter((t) => t.id !== id));
      toast({ title: "Test Deleted", description: "Test removed from directory" });
    } catch {
      toast({ title: "Delete Failed", description: "Could not delete test", variant: "destructive" });
    }
  }, [deleteTarget, toast]);

  // ---- Save Comment Edit ----
  const startEditComment = useCallback((testId: string, currentComment: string) => {
    setEditingCommentId(testId);
    setEditingCommentText(currentComment || "");
  }, []);

  const saveEditComment = useCallback(async () => {
    if (!editingCommentId) return;
    const id = editingCommentId;
    const newComment = editingCommentText;
    setEditingCommentId(null);
    setEditingCommentText("");
    try {
      const r = await fetch(`/api/pdb-tests/${id}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ comment: newComment }),
      });
      if (!r.ok) throw new Error("Failed to save");
      setTestsList((prev) => prev.map((t) => t.id === id ? { ...t, comment: newComment } : t));
      toast({ title: "Comment Updated" });
    } catch {
      toast({ title: "Save Failed", variant: "destructive" });
    }
  }, [editingCommentId, editingCommentText, toast]);

  // ---- Relative Date Helper ----
  const relativeDate = (dateStr: string) => {
    const now = new Date();
    const date = new Date(dateStr);
    const diffMs = now.getTime() - date.getTime();
    const diffMin = Math.floor(diffMs / 60000);
    if (diffMin < 1) return "just now";
    if (diffMin < 60) return `${diffMin}m ago`;
    const diffHrs = Math.floor(diffMin / 60);
    if (diffHrs < 24) return `${diffHrs}h ago`;
    const diffDays = Math.floor(diffHrs / 24);
    if (diffDays < 30) return `${diffDays}d ago`;
    return date.toLocaleDateString();
  };

  // Generation counter to discard stale async results when PDB ID changes rapidly
  const fetchGenRef = useRef(0);

  // ---- Step Handlers ----
  const handleFetchPdb = useCallback(async () => {
    const handleFetchError = async (r: Response) => {
      let msg = `Server error: ${r.status}`;
      try {
        const detail = await r.json();
        if (detail?.error) msg = detail.error;
        if (detail?.hint) msg += `. ${detail.hint}`;
      } catch { /* use default message */ }
      throw new Error(msg);
    };
    if (pdbId.trim()) {
      const gen = ++fetchGenRef.current; // Capture generation for this fetch
      setPdbLoading(true);
      try {
        const r = await apiCall("pdb-fetch", { pdb_id: pdbId.trim().toUpperCase() }, 30000);
        if (!r.ok) await handleFetchError(r);
        const data = await r.json();
        // Discard if user changed PDB ID while we were fetching
        if (gen !== fetchGenRef.current) return;
        setPdbInfo(data);
        if (data.chains?.length > 0) setSelectedChain(data.chains[0].chain_id);
        setPdbEntries((prev) => {
          const exists = prev.find((p) => p.pdb_id === data.pdb_id);
          return exists ? prev : [...prev, data];
        });
        // Also fetch raw PDB text for 3D viewer
        try {
          const rawR = await apiCall("pdb-raw-fetch", { pdb_id: pdbId.trim().toUpperCase() }, 30000);
          if (rawR.ok) {
            const rawData = await rawR.json();
            // Discard stale raw fetch result
            if (gen !== fetchGenRef.current) return;
            if (rawData.pdb_text) setPdbRawText(rawData.pdb_text);
          }
        } catch { /* ignore raw fetch failure */ }
        setOverlayResults([]); // Clear overlays when PDB changes
        toast({ title: "PDB Loaded", description: `${data.title?.slice(0, 60) || data.pdb_id} (${data.chains?.length || 0} chains, ${data.ligands?.length || 0} ligands)` });
      } catch (error: unknown) {
        if (gen !== fetchGenRef.current) return; // Discard stale error
        const msg = error instanceof Error ? error.message : "Unknown error";
        toast({ title: "Fetch Failed", description: msg, variant: "destructive" });
      } finally {
        if (gen === fetchGenRef.current) setPdbLoading(false);
      }
    } else if (pdbText.trim().length > 100) {
      setPdbLoading(true);
      try {
        const r = await apiCall("pdb-parse", { pdb_text: pdbText.trim() }, 30000);
        if (!r.ok) await handleFetchError(r);
        const data = await r.json();
        setPdbInfo(data);
        if (data.chains?.length > 0) setSelectedChain(data.chains[0].chain_id);
        // Pasted text is already raw PDB text for the 3D viewer
        setPdbRawText(pdbText.trim());
        toast({ title: "PDB Parsed", description: `${data.chains?.length} chains detected` });
      } catch (error: unknown) {
        const msg = error instanceof Error ? error.message : "Unknown error";
        toast({ title: "Parse Failed", description: msg, variant: "destructive" });
      } finally { setPdbLoading(false); }
    }
  }, [pdbId, pdbText, toast]);

  // Auto-load Step3 substrates into Step4 single reaction editor (only once per entry)
  const step4InitializedRef = useRef(false);
  // onChange handler: split reaction SMILES by ">>" to update reactant/product state
  const handleReactionEditorChange = useCallback((fullSmiles: string) => {
    const parts = fullSmiles.split(">>");
    setReactantSmiles((parts[0] || "").trim());
    setProductSmilesWf(parts.length > 1 ? (parts[1] || "").trim() : "");
  }, [setReactantSmiles, setProductSmilesWf]);
  useEffect(() => {
    if (currentStep !== 4) {
      step4InitializedRef.current = false;
      return;
    }
    if (step4InitializedRef.current) return;
    step4InitializedRef.current = true;
    const subSmiles = substrates.filter((s) => s.smiles.trim()).map((s) => s.smiles.trim()).join(".");
    if (subSmiles) {
      setReactantSmiles(subSmiles);
      // Initialize single reaction editor: "substrates>>" (empty product side)
      setTimeout(() => {
        reactionEditorRef.current?.setMolecule(`${subSmiles}>>`);
      }, 500);
    }
  }, [currentStep]);

  // ---- Fetch ligand distances when selected residues or PDB text changes ----
  useEffect(() => {
    if (!pdbRawText || selectedResidues.length === 0) {
      setLigandDistances({});
      return;
    }
    let cancelled = false;
    (async () => {
      setDistancesLoading(true);
      try {
        const data = await apiCallJson<{distances?: Record<string, number>; error?: string}>(
          "ligand-distances",
          {
            pdb_text: pdbRawText,
            selected_residues: selectedResidues.map(sr => ({
              chain: sr.chain,
              res_num: sr.res_num,
            })),
          },
          30000,
        );
        if (!cancelled && data.distances && !data.error) {
          setLigandDistances(data.distances);
        }
      } catch {
        // Silently fail — fallback to rough distance estimation
      } finally {
        if (!cancelled) setDistancesLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, [selectedResidues, pdbRawText]);

  // Copy from Reactants to Products: builds "reactants>>reactants" in the single editor
  const handleCopyFromReactants = useCallback(async () => {
    if (!reactantSmiles) return;
    const newReaction = `${reactantSmiles}>>${reactantSmiles}`;
    const ok = await reactionEditorRef.current?.setMolecule(newReaction);
    if (ok) {
      setProductSmilesWf(reactantSmiles);
      toast({ title: "Copied", description: "Reactants copied to product side. Any existing products are replaced." });
    } else {
      toast({ title: "Copy Failed", description: "Could not copy reactants to product side", variant: "destructive" });
    }
  }, [reactantSmiles, setProductSmilesWf, toast]);

  // Auto Map: single editor — directly call Indigo automap on the reaction
  const handleAutomap = useCallback(async () => {
    setAutomapping(true);
    try {
      const result = await reactionEditorRef.current?.automap();
      if (result) {
        // Re-read mapped SMILES to update state
        const mappedSmiles = await reactionEditorRef.current?.getSmiles();
        if (mappedSmiles && mappedSmiles.includes(">>")) {
          const mappedParts = mappedSmiles.split(">>");
          setReactantSmiles((mappedParts[0] || "").trim());
          setProductSmilesWf((mappedParts.length > 1 ? mappedParts[1] : "").trim());
        }
        toast({ title: "Atom Mapping Applied", description: "Indigo has assigned atom map numbers to the reaction" });
      } else {
        toast({ title: "Automap Failed", description: "Could not assign atom mapping. Make sure both reactants and products are drawn.", variant: "destructive" });
      }
    } catch {
      toast({ title: "Automap Error", description: "Failed to run automap", variant: "destructive" });
    } finally {
      setAutomapping(false);
    }
  }, [toast, setReactantSmiles, setProductSmilesWf]);

  // Verify & Correct Atom Mapping: calls MCS-based remap to check Indigo's mapping quality
  const handleVerifyMapping = useCallback(async () => {
    // Get current reaction SMILES from editor
    const currentSmiles = await reactionEditorRef.current?.getSmiles();
    if (!currentSmiles || !currentSmiles.includes(">>")) {
      toast({ title: "No Reaction", description: "Draw both reactants and products before verifying", variant: "destructive" });
      return;
    }
    const parts = currentSmiles.split(">>");
    const reactants = (parts[0] || "").trim();
    const products = (parts[1] || "").trim();
    if (!reactants || !products) {
      toast({ title: "Incomplete Reaction", description: "Both reactants and products must be defined", variant: "destructive" });
      return;
    }

    setVerifying(true);
    try {
      const data = await apiCallJson<{
        error?: string;
        original_bond_changes?: number;
        corrected_bond_changes?: number;
        corrected?: boolean;
        corrected_reactants?: string;
        corrected_products?: string;
      }>("verify-atom-mapping", { reactants, products });

      if (data.error) {
        toast({ title: "Verification Error", description: String(data.error), variant: "destructive" });
        return;
      }

      const origChanges = data.original_bond_changes ?? 0;
      const corrChanges = data.corrected_bond_changes ?? 0;
      const corrected = data.corrected === true;

      if (corrected && data.corrected_reactants && data.corrected_products) {
        // Apply corrected mapping to the editor
        const correctedRxn = `${data.corrected_reactants}>>${data.corrected_products}`;
        const ok = await reactionEditorRef.current?.setMolecule(correctedRxn);
        if (ok) {
          setReactantSmiles(data.corrected_reactants);
          setProductSmilesWf(data.corrected_products);
        }
        toast({
          title: "Mapping Corrected ✓",
          description: `Bond changes: ${origChanges} → ${corrChanges} (MCS re-mapping applied). Reaction center is now smaller and more chemically accurate.`,
          duration: 6000,
        });
      } else {
        toast({
          title: "Mapping Verified ✓",
          description: `Current mapping is optimal (${origChanges} bond changes). No correction needed.`,
          duration: 4000,
        });
      }
    } catch {
      toast({ title: "Verification Error", description: "Failed to verify atom mapping", variant: "destructive" });
    } finally {
      setVerifying(false);
    }
  }, [toast, setReactantSmiles, setProductSmilesWf]);

  const handleNewSchemeDraft = useCallback(() => {
    const subSmiles = substrates.filter((s) => s.smiles.trim()).map((s) => s.smiles.trim()).join(".");
    setReactantSmiles(subSmiles);
    setProductSmilesWf("");
    step4InitializedRef.current = false;
    setTimeout(() => {
      reactionEditorRef.current?.setMolecule(`${subSmiles}>>`);
    }, 300);
    toast({ title: "Scheme Created", description: `Reactants populated from ${substrates.length} molecule(s)` });
  }, [substrates, setProductSmilesWf, toast]);

  const handleWfRun = useCallback(async () => {
    if (!reactantSmiles.trim()) {
      toast({ title: "Input Required", description: "Define reactants in Step 4", variant: "destructive" });
      return;
    }
    setWfRunning(true);
    setWfShowProgress(true);
    setWfSearchProgress(null);
    setWfJobId(null);

    // Mark state as "run" when search starts
    if (activeTestId) {
      fetch(`/api/pdb-tests/${activeTestId}`, {
        method: "PUT",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ state: "run" }),
      }).catch(() => { /* silent */ });
    }
    try {
      // --- Read current reaction SMILES from single editor & detect forward-only mode ---
      let currentReactantSmiles = reactantSmiles.trim();
      let currentProductSmiles = productSmilesWf.trim();
      let hasProducts = currentProductSmiles.length > 0;
      // Double-check from editor: the editor might have products that state hasn't caught yet
      if (!hasProducts) {
        try {
          const fullSmiles = await reactionEditorRef.current?.getSmiles();
          if (fullSmiles && fullSmiles.includes(">>")) {
            const parts = fullSmiles.split(">>");
            currentReactantSmiles = (parts[0] || "").trim() || currentReactantSmiles;
            currentProductSmiles = (parts.length > 1 ? parts[1] : "").trim();
            hasProducts = currentProductSmiles.length > 0;
          }
        } catch { /* use state values */ }
      }

      // --- Auto-call Indigo automap (only when products exist) ---
      if (hasProducts) {
        try {
          await reactionEditorRef.current?.automap();
          const mappedSmiles = await reactionEditorRef.current?.getSmiles();
          if (mappedSmiles && mappedSmiles.includes(">>")) {
            const mappedParts = mappedSmiles.split(">>");
            currentReactantSmiles = (mappedParts[0] || "").trim() || currentReactantSmiles;
            currentProductSmiles = (mappedParts.length > 1 ? mappedParts[1] : "").trim() || currentProductSmiles;
            setReactantSmiles(currentReactantSmiles);
            setProductSmilesWf(currentProductSmiles);
            toast({ title: "Auto Map Applied", description: "Atom mapping automatically applied before search" });
          }
        } catch {
          // automap error — continue with unmapped SMILES
        }
      } else {
        toast({ title: "Forward-Only Mode", description: "No products drawn — running forward-only mechanism search. For bidirectional search, draw expected products in Step 4." });
      }

      // Build FULL active site configuration
      const reqBody: Record<string, unknown> = {
        reactants: currentReactantSmiles,
        products: currentProductSmiles || '',
        max_configs: wfMaxConfigs,
        max_bond_length: wfMaxBondLength,
        ...(wfMcsaId > 0 ? { mcsa_id: wfMcsaId } : {}),
      };

      // Export atom-mapped RXN from single editor (only when products exist)
      if (hasProducts) {
        try {
          const rxnText = await reactionEditorRef.current?.getRxn();
          if (rxnText && rxnText.includes('$RXN')) {
            reqBody.rxn_text = rxnText;
          }
        } catch { /* getRxn failed — fall through to SMILES */ }
      }

      // --- Residue validation ---
      const unknownResidues = selectedResidues.filter(
        sr => !VALID_RESIDUE_NAMES.has(sr.res_name?.toUpperCase())
      );
      if (unknownResidues.length > 0) {
        const names = [...new Set(unknownResidues.map(r => r.res_name))].join(', ');
        toast({
          title: "Unknown Residue Types",
          description: `The following residue types are not recognized and will be excluded: ${names}. Most M-CSA rules require standard catalytic residues.`,
          variant: "destructive",
        });
      }
      const validResidues = unknownResidues.length > 0
        ? selectedResidues.filter(sr => VALID_RESIDUE_NAMES.has(sr.res_name?.toUpperCase()))
        : selectedResidues;
      if (validResidues.length > 0) {
        reqBody.residues = validResidues.map((sr) => ({
          res_name: sr.res_name,
          res_num: sr.res_num,
          chain: sr.chain,
          part: sr.part,
        }));
      }
      if (pdbInfo?.pdb_id) {
        reqBody.pdb_id = pdbInfo.pdb_id;
      }
      if (pdbRawText && pdbRawText.length > 100) {
        reqBody.pdb_text = pdbRawText;
      }

      // Pass ligand_mappings
      // Include all types (substrate, cofactor, water) that have a mapping —
      // water molecules mapped to HOH residues need their own independent PDB
      // coordinates, and the SMILES-keyed cache can't distinguish multiple O's.
      const ligandMappings = substrates
        .filter(s => s.mapped_ligand && s.smiles.trim())
        .map(s => {
          const parts = s.mapped_ligand!.split('_');
          return {
            smiles: s.smiles,
            res_name: parts[0] || '',
            chain: parts[1] || '',
            res_num: parseInt(parts[2]) || 0,
            ligand_type: s.type,  // 'substrate' | 'cofactor' | 'water'
          };
        });
      if (ligandMappings.length > 0) {
        reqBody.ligand_mappings = ligandMappings;
      }

      // ─── ASYNC SEARCH: Start + Poll ───
      // Step 1: Start the async search
      const startResult = await apiCallJson<{ jobId: string; state: string }>(
        "mechanism-search-start", reqBody, 30000
      );
      const jobId = startResult.jobId;
      setWfJobId(jobId);
      setWfSearchProgress({
        state: "STARTING",
        explored_nodes: 0,
        total_nodes: 0,
        current_iteration: 0,
        elapsed_seconds: 0,
      });

      // Save jobId to DB so it can be resumed if user navigates away
      if (activeTestId) {
        fetch(`/api/pdb-tests/${activeTestId}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ jobId }),
        }).catch(() => { /* silent */ });
      }

      // Step 2: Poll for progress every 3 seconds
      let localSearchResult: Record<string, any> | null = null;
      const POLL_INTERVAL = 3000;
      const MAX_POLL_TIME = 30 * 60 * 1000; // 30 min max
      const pollStart = Date.now();

      while (Date.now() - pollStart < MAX_POLL_TIME) {
        // Wait before polling (except the first time)
        await new Promise(r => setTimeout(r, POLL_INTERVAL));

        try {
          const status = await apiCallJson<SearchProgress>(
            "mechanism-search-status", { jobId }, 15000
          );
          setWfSearchProgress(status);

          if (status.state === "DONE") {
            // Search completed — extract result
            if (status.result) {
              const normalized = normalizeSearchResult(status.result as unknown as Record<string, unknown>);
              localSearchResult = normalized as unknown as MechanismSearchResult;
              setWfResult(normalized as unknown as MechanismSearchResult);
              setWfShowProgress(true); // Ensure results view is shown even if user went back to Step 5
              if (Array.isArray(normalized.paths) && normalized.paths.length > 0) {
                setWfActiveTab("graph");
              }
            }
            break;
          } else if (status.state === "ERROR") {
            const errorMsg = status.error || "Search failed with an unknown error";
            localSearchResult = { error: errorMsg } as unknown as MechanismSearchResult;
            setWfResult(localSearchResult as unknown as MechanismSearchResult);
            setWfShowProgress(true); // Show error results view
            break;
          } else if (status.state === "CANCELLED") {
            localSearchResult = { error: "Search was cancelled" } as unknown as MechanismSearchResult;
            setWfResult(localSearchResult as unknown as MechanismSearchResult);
            setWfShowProgress(true); // Show cancelled results view
            break;
          }
          // Still RUNNING or STARTING — continue polling
        } catch (pollErr) {
          // Network hiccup during polling — log and retry
          console.warn("Poll error (will retry):", pollErr);
        }
      }

      // Timeout check
      if (!localSearchResult && Date.now() - pollStart >= MAX_POLL_TIME) {
        localSearchResult = { error: "Search timed out (30 minutes). The search is still running on the server — results may appear later." } as unknown as MechanismSearchResult;
        setWfResult(localSearchResult as unknown as MechanismSearchResult);
        setWfShowProgress(true); // Show timeout results view
      }

      // Determine result state: "done" if search found paths, "fail" otherwise
      const searchRes = localSearchResult as Record<string, any> | null;
      const hasAnyResult = (searchRes?.paths as any[])?.length > 0;
      const finalState = hasAnyResult ? "done" : "fail";

      // Save results to DB after completion
      if (activeTestId) {
        fetch(`/api/pdb-tests/${activeTestId}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({
            ...buildSavePayload(5, { state: finalState }),
            searchResultJson: localSearchResult ? JSON.stringify(localSearchResult) : null,
          }),
        }).then(() => {
          // Refresh directory list after DB save so state is up-to-date
          fetchTests();
        }).catch(() => { /* silent */ });
      }

      if (finalState === "done") {
        toast({ title: "Search Complete", description: "Results ready" });
      } else if (searchRes?.error) {
        toast({ title: "Search Error", description: String(searchRes.error).slice(0, 150), variant: "destructive" });
      } else {
        toast({ title: "No Results Found", description: "No matching mechanisms or products were found", variant: "destructive" });
      }
    } catch (error: unknown) {
      const msg = error instanceof Error ? error.message : "Unknown error";
      toast({ title: "Search Failed", description: msg, variant: "destructive" });
      setWfShowProgress(false); // Return to Step 5 on failure
      // Revert state to "edit" on failure
      if (activeTestId) {
        fetch(`/api/pdb-tests/${activeTestId}`, {
          method: "PUT",
          headers: { "Content-Type": "application/json" },
          body: JSON.stringify({ state: "edit" }),
        }).catch(() => { /* silent */ });
      }
    } finally {
      setWfRunning(false);
      setWfSearchProgress(null);
    }
  }, [reactantSmiles, productSmilesWf, wfMaxConfigs, wfMaxBondLength, toast, selectedResidues, pdbInfo, activeTestId, wfActiveTab, buildSavePayload, substrates, pdbRawText, fetchTests]);

  const currentChainResidues = pdbInfo?.chains?.find((c) => c.chain_id === selectedChain)?.residues || [];
  const filteredResidues = currentChainResidues.filter((r) => {
    if (!residueFilter) return true;
    const q = residueFilter.toLowerCase();
    return r.res_name.toLowerCase().includes(q) || String(r.res_num).includes(q);
  });

  const wfFilteredMatches: any[] = [];

  // ==== DIRECTORY VIEW (when testId is null) ====
  if (!activeTestId) {
    return (
      <Card className="border-violet-100 shadow-sm">
        <CardHeader className="pb-3">
          <div className="flex items-center justify-between">
            <CardTitle className="flex items-center gap-2 text-sm">
              <FolderOpen className="w-4 h-4 text-violet-600" />PDB Test Directory
            </CardTitle>
            <Button
              onClick={handleCreateNew}
              className="bg-violet-600 hover:bg-violet-700 text-white h-8 text-xs"
            >
              <Plus className="w-3 h-3 mr-1" />Create New PDB Test
            </Button>
          </div>
          <CardDescription className="text-xs">Manage your saved PDB workflow tests</CardDescription>
        </CardHeader>
        <CardContent>
          {testsLoading ? (
            <div className="flex items-center justify-center py-12">
              <Loader2 className="w-5 h-5 text-violet-500 animate-spin mr-2" />
              <span className="text-xs text-gray-500">Loading tests...</span>
            </div>
          ) : testsList.length === 0 ? (
            <div className="flex flex-col items-center justify-center py-16 space-y-4">
              <div className="w-16 h-16 rounded-2xl bg-violet-50 flex items-center justify-center">
                <FolderOpen className="w-8 h-8 text-violet-300" />
              </div>
              <div className="text-center space-y-1">
                <p className="text-sm font-medium text-gray-600">No tests yet</p>
                <p className="text-xs text-gray-400">Create a new PDB test to start your enzyme mechanism search.</p>
              </div>
              <Button
                onClick={handleCreateNew}
                className="bg-violet-600 hover:bg-violet-700 text-white h-9 text-xs"
              >
                <Plus className="w-3 h-3 mr-1" />Create New PDB Test
              </Button>
            </div>
          ) : (
            <div className="max-h-[480px] overflow-y-auto rounded-lg border border-gray-100">
              <Table>
                <TableHeader>
                  <TableRow className="bg-gray-50/80">
                    <TableHead className="text-xs h-8 font-semibold">Test ID</TableHead>
                    <TableHead className="text-xs h-8 font-semibold">PDB ID</TableHead>
                    <TableHead className="text-xs h-8 font-semibold">State</TableHead>
                    <TableHead className="text-xs h-8 font-semibold">Step</TableHead>
                    <TableHead className="text-xs h-8 font-semibold">Comment</TableHead>
                    <TableHead className="text-xs h-8 font-semibold">Created</TableHead>
                    <TableHead className="text-xs h-8 font-semibold w-24">Actions</TableHead>
                  </TableRow>
                </TableHeader>
                <TableBody>
                  {testsList.map((test) => (
                    <TableRow key={test.id} className="hover:bg-violet-50/50">
                      <TableCell className="py-2">
                        <span className="text-xs font-mono font-medium text-violet-700">{test.testId}</span>
                      </TableCell>
                      <TableCell className="py-2">
                        {test.pdbId ? (
                          <Badge variant="outline" className="text-xs border-violet-200 text-violet-600 font-mono">{test.pdbId}</Badge>
                        ) : (
                          <span className="text-xs text-gray-300 italic">—</span>
                        )}
                      </TableCell>
                      <TableCell className="py-2">
                        <Badge
                          variant="secondary"
                          className={`text-xs ${
                            test.state === "done"
                              ? "bg-emerald-50 text-emerald-700 border border-emerald-200"
                              : test.state === "fail"
                                ? "bg-red-50 text-red-700 border border-red-200"
                                : test.state === "run"
                                  ? "bg-amber-50 text-amber-700 border border-amber-200"
                                  : "bg-blue-50 text-blue-700 border border-blue-200"
                          }`}
                        >
                          {test.state === "done" ? "✓ done" : test.state === "fail" ? "✗ fail" : test.state === "run" ? "● run" : "✎ edit"}
                        </Badge>
                      </TableCell>
                      <TableCell className="py-2">
                        <span className="text-xs font-mono text-gray-600">{test.currentStep}/5</span>
                      </TableCell>
                      <TableCell className="py-2 max-w-[200px]">
                        {editingCommentId === test.id ? (
                          <div className="flex items-center gap-1">
                            <Input
                              value={editingCommentText}
                              onChange={(e) => setEditingCommentText(e.target.value)}
                              onKeyDown={(e) => { if (e.key === "Enter") saveEditComment(); if (e.key === "Escape") { setEditingCommentId(null); setEditingCommentText(""); } }}
                              className="h-6 text-xs min-w-[100px] px-1.5 py-0"
                              autoFocus
                              placeholder="Add comment..."
                            />
                            <Button size="sm" variant="ghost" className="h-5 w-5 p-0 text-emerald-600 hover:bg-emerald-50"
                              onClick={saveEditComment}>
                              <CheckCircle2 className="w-3 h-3" />
                            </Button>
                            <Button size="sm" variant="ghost" className="h-5 w-5 p-0 text-gray-400 hover:bg-gray-100"
                              onClick={() => { setEditingCommentId(null); setEditingCommentText(""); }}>
                              <X className="w-3 h-3" />
                            </Button>
                          </div>
                        ) : (
                          <span
                            className="text-xs text-gray-500 truncate block cursor-pointer hover:text-violet-600 hover:bg-violet-50/50 rounded px-1 py-0.5 -mx-1 transition-colors"
                            onClick={() => startEditComment(test.id, test.comment)}
                            title="Click to edit"
                          >
                            {test.comment || <span className="italic text-gray-300">Click to add...</span>}
                          </span>
                        )}
                      </TableCell>
                      <TableCell className="py-2">
                        <span className="text-xs text-gray-400">{relativeDate(test.createdAt)}</span>
                      </TableCell>
                      <TableCell className="py-2">
                        <div className="flex items-center gap-1">
                          <Button
                            variant="ghost"
                            size="sm"
                            className="h-7 text-xs text-violet-600 hover:text-violet-800 hover:bg-violet-50 px-2"
                            onClick={() => handleViewTest(test.id)}
                          >
                            <Eye className="w-3 h-3 mr-1" />View
                          </Button>
                          <Button
                            variant="ghost"
                            size="sm"
                            className="h-7 text-xs text-red-400 hover:text-red-600 hover:bg-red-50 px-1.5"
                            onClick={() => confirmDeleteTest(test.id)}
                          >
                            <Trash2 className="w-3 h-3" />
                          </Button>
                        </div>
                      </TableCell>
                    </TableRow>
                  ))}
                </TableBody>
              </Table>
            </div>
          )}
        </CardContent>

        {/* Delete Confirmation Dialog */}
        <AlertDialog open={!!deleteTarget} onOpenChange={(open) => { if (!open) setDeleteTarget(null); }}>
          <AlertDialogContent>
            <AlertDialogHeader>
              <AlertDialogTitle>Delete Test</AlertDialogTitle>
              <AlertDialogDescription>
                Are you sure you want to delete this PDB test? This action cannot be undone. All saved progress and results will be permanently removed.
              </AlertDialogDescription>
            </AlertDialogHeader>
            <AlertDialogFooter>
              <AlertDialogCancel>Cancel</AlertDialogCancel>
              <AlertDialogAction
                onClick={handleDeleteTest}
                className="bg-red-600 hover:bg-red-700 text-white"
              >
                Delete
              </AlertDialogAction>
            </AlertDialogFooter>
          </AlertDialogContent>
        </AlertDialog>
      </Card>
    );
  }

  // ==== PROGRESS / RESULTS VIEW (shown when search is running OR results exist) ====
  if (showProgressOrResults) {
    return (
      <div className="space-y-4">
        <div className="flex items-center gap-3">
          {wfRunning ? (
            <Button variant="outline" size="sm" className="h-7 text-xs border-amber-200 text-amber-600 hover:bg-amber-50"
              onClick={goToDirectory}>
              <ChevronLeft className="w-3 h-3 mr-1" />Back to Directory
            </Button>
          ) : (
            <>
              <Button variant="outline" size="sm" className="h-7 text-xs border-violet-200 text-violet-600 hover:bg-violet-50"
                onClick={goToDirectory}>
                <FolderOpen className="w-3 h-3 mr-1" />Directory
              </Button>
              <Button variant="outline" size="sm" className="h-7 text-xs border-gray-200 text-gray-600 hover:bg-gray-50"
                onClick={() => { setWfResult(null); setWfShowProgress(false); setCurrentStep(5); }}>
                <ChevronLeft className="w-3 h-3 mr-1" />Back to Step 5
              </Button>
            </>
          )}
          <div className="flex items-center gap-2">
            {currentTestIdLabel && (
              <Badge variant="outline" className="bg-violet-50 text-violet-600 text-xs border-violet-200 font-mono">
                {currentTestIdLabel}
              </Badge>
            )}
            {pdbInfo && <Badge variant="secondary" className="bg-violet-50 text-violet-700 text-xs"><FileText className="w-3 h-3 mr-1" />{pdbInfo.pdb_id}</Badge>}
            {selectedResidues.length > 0 && <Badge variant="secondary" className="bg-emerald-50 text-emerald-700 text-xs"><Dna className="w-3 h-3 mr-1" />{selectedResidues.length} residues</Badge>}
          </div>
        </div>

        {/* Search Progress — shown while search is running */}
        {wfRunning && (
          <Card className="border-amber-100 shadow-sm">
            <CardHeader className="pb-3">
              <CardTitle className="flex items-center gap-2 text-sm">
                <Loader2 className="w-4 h-4 text-amber-500 animate-spin" />
                Search in Progress
              </CardTitle>
              <CardDescription className="text-xs">
                {wfSearchProgress?.state === "STARTING"
                  ? "Initializing search engine..."
                  : "Exploring mechanism pathways..."}
              </CardDescription>
            </CardHeader>
            <CardContent className="space-y-4">
              {/* Progress Stats */}
              {wfSearchProgress && wfSearchProgress.current_iteration > 0 && (
                <>
                  <div className="grid grid-cols-2 sm:grid-cols-4 gap-3">
                    <div className="rounded-lg bg-amber-50 p-3 text-center">
                      <div className="text-lg font-bold font-mono text-amber-700">
                        {wfSearchProgress.current_iteration}
                      </div>
                      <div className="text-xs text-amber-600">Iteration</div>
                    </div>
                    <div className="rounded-lg bg-violet-50 p-3 text-center">
                      <div className="text-lg font-bold font-mono text-violet-700">
                        {wfSearchProgress.explored_nodes}
                      </div>
                      <div className="text-xs text-violet-600">Explored Nodes</div>
                    </div>
                    <div className="rounded-lg bg-emerald-50 p-3 text-center">
                      <div className="text-lg font-bold font-mono text-emerald-700">
                        {wfSearchProgress.total_nodes}
                      </div>
                      <div className="text-xs text-emerald-600">Total Nodes</div>
                    </div>
                    <div className="rounded-lg bg-blue-50 p-3 text-center">
                      <div className="text-lg font-bold font-mono text-blue-700">
                        {wfSearchProgress.elapsed_seconds >= 60
                          ? `${Math.floor(wfSearchProgress.elapsed_seconds / 60)}m ${Math.round(wfSearchProgress.elapsed_seconds % 60)}s`
                          : `${Math.round(wfSearchProgress.elapsed_seconds)}s`}
                      </div>
                      <div className="text-xs text-blue-600">Elapsed</div>
                    </div>
                  </div>

                  {/* Progress Bar */}
                  <div className="space-y-1.5">
                    <div className="flex items-center justify-between text-xs text-gray-500">
                      <span>Exploring configurations...</span>
                      <span>{Math.round((wfSearchProgress.explored_nodes / Math.max(wfSearchProgress.total_nodes, 1)) * 100)}%</span>
                    </div>
                    <div className="w-full bg-gray-100 rounded-full h-2.5">
                      <div
                        className="bg-gradient-to-r from-amber-400 via-violet-500 to-emerald-500 h-2.5 rounded-full transition-all duration-700"
                        style={{ width: `${Math.min(100, Math.round((wfSearchProgress.explored_nodes / Math.max(wfSearchProgress.total_nodes, 1)) * 100))}%` }}
                      />
                    </div>
                  </div>
                </>
              )}

              {/* Starting indicator */}
              {(!wfSearchProgress || wfSearchProgress.current_iteration === 0) && (
                <div className="flex flex-col items-center justify-center py-8 space-y-3">
                  <div className="relative">
                    <div className="w-16 h-16 rounded-full border-4 border-amber-100" />
                    <div className="absolute inset-0 w-16 h-16 rounded-full border-4 border-amber-400 border-t-transparent animate-spin" />
                  </div>
                  <p className="text-sm text-gray-500">Initializing search engine...</p>
                  <p className="text-xs text-gray-400">This may take a moment to load rules and set up the search graph</p>
                </div>
              )}

              {/* Info about background search */}
              <div className="rounded-lg bg-amber-50/50 border border-amber-100 p-3">
                <p className="text-xs text-amber-700">
                  💡 The search runs on the server. You can click "Back to Step 5" and return later — the search will continue in the background.
                  Click "View" from the Directory to check progress.
                </p>
              </div>
            </CardContent>
          </Card>
        )}

        {/* Results — shown when search is complete */}
        {hasResults && !wfRunning && (
          <Tabs value={wfActiveTab} onValueChange={setWfActiveTab}>
            <Card className="border-violet-100 shadow-sm">
              <CardHeader className="pb-2">
                <div className="flex items-center justify-between flex-wrap gap-2">
                  <CardTitle className="flex items-center gap-2 text-sm">
                    <Sparkles className="w-4 h-4 text-violet-600" />EzMechanism Results
                  </CardTitle>
                  <TabsList className="bg-violet-50 p-0.5 h-8">
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
                <TabsContent value="mechanism"><MechanismPanel searchResult={wfResult} onExploreProduct={setWfSelectedSmiles} onSwitchToGraph={() => setWfActiveTab("graph")} /></TabsContent>
                <TabsContent value="graph"><MechanismGraphPanel result={wfResult} /></TabsContent>
                <TabsContent value="molecules"><MoleculesPanel searchResult={wfResult} onExploreProduct={setWfSelectedSmiles} /></TabsContent>
                <TabsContent value="rules"><RulesPanel matchResult={null} filter={wfRulesFilter} onFilterChange={setWfRulesFilter} filteredMatches={wfFilteredMatches} /></TabsContent>
                <TabsContent value="stats"><StatsPanel searchResult={wfResult} /></TabsContent>
              </CardContent>
            </Card>
          </Tabs>
        )}
      </div>
    );
  }

  return (
    <div className="space-y-4">
      <StepWizard currentStep={currentStep} onStepClick={goToStep} />

      {/* Step 1: PDB Structure */}
      {currentStep === 1 && (
        <Card className="border-violet-100 shadow-sm">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center gap-2 text-sm"><FileText className="w-4 h-4 text-violet-600" />Step 1: Choose PDB Structure</CardTitle>
            <CardDescription className="text-xs">Fetch a PDB structure by ID or paste PDB format text</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
              <div className="space-y-2">
                <Label className="text-xs font-medium">PDB ID</Label>
                <div className="flex gap-2">
                  <Input placeholder="e.g., 1EJG" value={pdbId} onChange={(e) => setPdbId(e.target.value.toUpperCase().slice(0, 4))}
                    className="font-mono text-xs border-violet-200 focus:border-violet-400 h-9 uppercase" maxLength={4} />
                  <Button onClick={handleFetchPdb} disabled={pdbLoading || !pdbId.trim()} className="bg-violet-600 hover:bg-violet-700 text-white h-9 text-xs shrink-0">
                    {pdbLoading ? <Loader2 className="w-3 h-3 mr-1 animate-spin" /> : <Download className="w-3 h-3 mr-1" />}Fetch
                  </Button>
                </div>
              </div>
              <div className="space-y-2">
                <Label className="text-xs font-medium">Or Paste PDB Text</Label>
                <Textarea placeholder="Paste PDB format text here..." value={pdbText} onChange={(e) => setPdbText(e.target.value)}
                  className="font-mono text-xs border-violet-200 focus:border-violet-400 min-h-[80px] resize-y" />
                <Button onClick={handleFetchPdb} disabled={pdbLoading || pdbText.trim().length < 100}
                  variant="outline" className="border-violet-200 text-violet-600 hover:bg-violet-50 h-8 text-xs">
                  {pdbLoading ? <Loader2 className="w-3 h-3 mr-1 animate-spin" /> : <FileText className="w-3 h-3 mr-1" />}Parse PDB
                </Button>
              </div>
            </div>

            {pdbRawText && (
              <PdbStructureViewer
                pdbData={pdbRawText}
                activeChain={selectedChain}
                chains={pdbInfo?.chains?.map((c) => c.chain_id) || []}
                onChainChange={(chain) => setSelectedChain(chain)}
                height={600}
                className="mb-4"
                autoZoom="all"
                overlayMolecules={overlayMolecules}
              />
            )}

            {pdbInfo && (
              <div className="rounded-lg border border-emerald-100 bg-emerald-50/30 p-4 space-y-3">
                <div className="flex items-center gap-2">
                  <CheckCircle2 className="w-4 h-4 text-emerald-600" />
                  <span className="text-sm font-semibold text-emerald-800">{pdbInfo.pdb_id} — {pdbInfo.title?.slice(0, 80)}</span>
                </div>
                <div className="grid grid-cols-2 sm:grid-cols-4 gap-2 text-xs">
                  <div className="rounded-md bg-white p-2 text-center"><div className="text-sm font-bold font-mono text-violet-700">{pdbInfo.chains?.length || 0}</div><div className="text-xs text-gray-500">Chains</div></div>
                  <div className="rounded-md bg-white p-2 text-center"><div className="text-sm font-bold font-mono text-violet-700">{pdbInfo.ligands?.length || 0}</div><div className="text-xs text-gray-500">Ligands</div></div>
                  <div className="rounded-md bg-white p-2 text-center"><div className="text-sm font-bold font-mono text-violet-700">{pdbInfo.resolution || "N/A"}</div><div className="text-xs text-gray-500">Resolution (Å)</div></div>
                  <div className="rounded-md bg-white p-2 text-center"><div className="text-sm font-bold font-mono text-violet-700">{pdbInfo.total_atoms || "N/A"}</div><div className="text-xs text-gray-500">Atoms</div></div>
                </div>
                {pdbInfo.chains && pdbInfo.chains.length > 0 && (
                  <div className="flex flex-wrap gap-1">
                    {pdbInfo.chains.map((c) => (
                      <Badge key={c.chain_id} variant="outline" className="text-xs border-violet-200 text-violet-700">Chain {c.chain_id} ({c.num_residues})</Badge>
                    ))}
                  </div>
                )}
                {pdbInfo.ligands && pdbInfo.ligands.length > 0 && (
                  <div className="flex flex-wrap gap-1">
                    {pdbInfo.ligands.slice(0, 8).map((l, i) => (
                      <Badge key={i} variant="secondary" className="text-xs bg-amber-50 text-amber-700">{l.res_name} ({l.chain}:{l.res_num})</Badge>
                    ))}
                    {pdbInfo.ligands.length > 8 && <span className="text-xs text-gray-400">+{pdbInfo.ligands.length - 8} more</span>}
                  </div>
                )}
                {pdbInfo.chains?.some((c) => c.uniprot_accession) && (
                  <div className="flex flex-wrap gap-1">
                    {pdbInfo.chains.filter((c) => c.uniprot_accession).map((c, i) => (
                      <Badge key={i} variant="outline" className="text-xs border-emerald-200 text-emerald-700">Chain {c.chain_id}: {c.uniprot_accession}</Badge>
                    ))}
                  </div>
                )}

                <div className="flex items-center gap-2 flex-wrap">
                  <Button variant="outline" size="sm" className="h-7 text-xs border-gray-200 text-gray-500 hover:bg-gray-50 hover:text-gray-700"
                    onClick={goToDirectory}>
                    <FolderOpen className="w-3 h-3 mr-1" />Directory
                  </Button>
                  <Button onClick={() => goToStep(2)} className="bg-violet-600 hover:bg-violet-700 text-white h-8 text-xs">
                    Proceed to Step 2 <ChevronRight className="w-3 h-3 ml-1" />
                  </Button>
                </div>
                {pdbInfo && (
                  <div className="flex items-center gap-2 text-xs text-emerald-600 pt-2 border-t border-gray-100 mt-2">
                    <CheckCircle2 className="w-3 h-3" />
                    <span>PDB structure data (chains, ligands, coordinates) will be available for residue selection in Step 2</span>
                  </div>
                )}
              </div>
            )}

            {pdbEntries.length > 1 && (
              <div className="space-y-2">
                <Label className="text-xs text-gray-400 uppercase tracking-wider">Previous PDB Entries</Label>
                <div className="max-h-32 overflow-y-auto">
                  <Table>
                    <TableHeader><TableRow><TableHead className="text-xs h-6">PDB ID</TableHead><TableHead className="text-xs h-6">Title</TableHead><TableHead className="text-xs h-6 w-16">Chains</TableHead></TableRow></TableHeader>
                    <TableBody>
                      {pdbEntries.map((e) => (
                        <TableRow key={e.pdb_id} className="cursor-pointer hover:bg-violet-50" onClick={() => { setPdbInfo(e); if (e.chains?.length > 0) setSelectedChain(e.chains[0].chain_id); }}>
                          <TableCell className="text-xs font-mono font-medium py-1">{e.pdb_id}</TableCell>
                          <TableCell className="text-xs py-1 max-w-[200px] truncate">{e.title}</TableCell>
                          <TableCell className="text-xs py-1 text-center">{e.chains?.length || 0}</TableCell>
                        </TableRow>
                      ))}
                    </TableBody>
                  </Table>
                </div>
              </div>
            )}
          </CardContent>
        </Card>
      )}

      {/* Step 2: Catalytic Residues */}
      {currentStep === 2 && (
        <Card className="border-violet-100 shadow-sm">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center gap-2 text-sm"><Dna className="w-4 h-4 text-violet-600" />Step 2: Choose Catalytic Residues</CardTitle>
            <CardDescription className="text-xs">Select chains, browse residues, or auto-predict active site</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {!pdbInfo ? (
              <Alert><AlertCircle className="h-4 w-4 text-amber-500" /><AlertTitle className="text-sm">No PDB Loaded</AlertTitle><AlertDescription className="text-xs">Go back to Step 1 and load a PDB structure first.</AlertDescription></Alert>
            ) : (
              <>
                {/* 3D Structure Viewer */}
                {pdbRawText && (
                  <PdbStructureViewer
                    ref={viewerStep2Ref}
                    pdbData={pdbRawText}
                    selectedResidues={selectedResidues}
                    activeChain={selectedChain}
                    chains={pdbInfo?.chains?.map((c) => c.chain_id) || []}
                    onChainChange={(chain) => setSelectedChain(chain)}
                    showActiveCenter={true}
                    autoZoom="chain"
                    onResidueClick={(chain, resNum) => {
                      setSelectedResidues((prev) => {
                        const exists = prev.some((sr) => sr.res_num === resNum && sr.chain === chain);
                        if (exists) return prev.filter((sr) => !(sr.res_num === resNum && sr.chain === chain));
                        const chainData = pdbInfo?.chains?.find((c) => c.chain_id === chain);
                        const res = chainData?.residues?.find((r) => r.res_num === resNum);
                        return [...prev, { res_name: res?.res_name || "UNK", res_num: resNum, chain, part: "side_chain" as const }];
                      });
                    }}
                    height={600}
                    overlayMolecules={overlayMolecules}
                  />
                )}

                {/* Chain selection + filter — below 3D viewer */}
                <div className="flex items-center gap-2 flex-wrap">
                  <Label className="text-xs font-medium">Chain:</Label>
                  <div className="flex gap-1">
                    {pdbInfo.chains?.map((c) => (
                      <Button key={c.chain_id} variant={selectedChain === c.chain_id ? "default" : "outline"} size="sm"
                        className={`h-7 text-xs px-3 ${selectedChain === c.chain_id ? "bg-violet-600 hover:bg-violet-700 text-white" : "border-violet-200 text-violet-600 hover:bg-violet-50"}`}
                        onClick={() => setSelectedChain(c.chain_id)}>
                        {c.chain_id} ({c.num_residues})
                      </Button>
                    ))}
                  </div>
                  <Input placeholder="Filter by name or number..." value={residueFilter} onChange={(e) => setResidueFilter(e.target.value)}
                    className="text-xs border-violet-200 h-8 flex-1 min-w-[140px]" />
                </div>

                <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
                  {/* Residue Browser */}
                  <div className="space-y-2">
                    <Label className="text-xs text-gray-400 uppercase tracking-wider">Residues ({filteredResidues.length})</Label>
                    <div className="max-h-80 overflow-y-auto rounded-md border border-gray-100">
                      <Table>
                        <TableHeader><TableRow><TableHead className="text-xs h-6">Name</TableHead><TableHead className="text-xs h-6 w-[72px]" title="Residue number in PDB chain">PDB Pos</TableHead><TableHead className="text-xs h-6 w-[72px]" title="Residue position in UniProt protein sequence">UniProt Pos</TableHead><TableHead className="text-xs h-6 w-16"></TableHead></TableRow></TableHeader>
                        <TableBody>
                          {filteredResidues.map((r) => {
                            const isSelected = selectedResidues.some((sr) => sr.res_num === r.res_num && sr.chain === r.chain);
                            const hasUniProtMapping = r.seq_pos > 0 && r.seq_db_offset !== undefined;
                            return (
                              <TableRow key={`${r.res_name}-${r.res_num}`} className={`text-xs ${isSelected ? "bg-violet-50" : ""}`}>
                                <TableCell className="font-mono font-medium py-1">{r.res_name}</TableCell>
                                <TableCell className="font-mono py-1 text-gray-800">{r.res_num}{r.insertion && r.insertion !== " " ? r.insertion : ""}</TableCell>
                                <TableCell className="font-mono py-1">
                                  {hasUniProtMapping ? (
                                    <span className="text-emerald-700" title={`UniProt position: ${r.seq_pos} (offset ${(r.seq_db_offset ?? 0) >= 0 ? "+" : ""}${r.seq_db_offset ?? 0})`}>{r.seq_pos}</span>
                                  ) : (
                                    <span className="text-gray-300 italic" title="No UniProt mapping — DBREF record not found or residue outside mapping range">—</span>
                                  )}
                                </TableCell>
                                <TableCell className="py-1">
                                  <Button size="sm" variant={isSelected ? "secondary" : "outline"}
                                    className={`h-6 text-xs px-2 ${isSelected ? "bg-violet-100 text-violet-700" : "border-gray-200 hover:border-violet-300"}`}
                                    onClick={() => {
                                      if (isSelected) {
                                        setSelectedResidues((prev) => prev.filter((sr) => !(sr.res_num === r.res_num && sr.chain === r.chain)));
                                      } else {
                                        setSelectedResidues((prev) => [...prev, { res_name: r.res_name, res_num: r.res_num, chain: r.chain, part: "side_chain" }]);
                                      }
                                    }}>
                                    {isSelected ? "✓" : "+"}
                                  </Button>
                                </TableCell>
                              </TableRow>
                            );
                          })}
                        </TableBody>
                      </Table>
                    </div>
                  </div>

                  {/* Selected Residues */}
                  <div className="space-y-2">
                    <div className="flex items-center justify-between">
                      <Label className="text-xs text-gray-400 uppercase tracking-wider">Selected Residues ({selectedResidues.length})</Label>
                      {selectedResidues.length > 1 && (
                        <Button
                          variant="ghost"
                          size="sm"
                          className="h-5 text-xs text-red-400 hover:text-red-600 hover:bg-red-50 px-1.5 gap-1"
                          onClick={() => setSelectedResidues([])}
                        >
                          <X className="w-2.5 h-2.5" />Clear All
                        </Button>
                      )}
                    </div>
                    {selectedResidues.length === 0 ? (
                      <div className="rounded-md border border-dashed border-gray-200 p-4 text-center">
                        <p className="text-xs text-gray-400">Click &quot;+&quot; to add catalytic residues</p>
                      </div>
                    ) : (
                      <div className="max-h-64 overflow-y-auto space-y-1">
                        {selectedResidues.map((sr, i) => {
                          // Look up UniProt position from PDB residue data
                          const chainData = pdbInfo?.chains?.find((c) => c.chain_id === sr.chain);
                          const resData = chainData?.residues?.find((r) => r.res_num === sr.res_num);
                          const uniProtPos = resData && resData.seq_pos > 0 ? resData.seq_pos : null;
                          return (
                          <div key={i} className="flex items-center gap-2 rounded-md border border-gray-100 bg-white px-3 py-2">
                            <Badge variant="outline" className="text-xs border-violet-200 text-violet-700 font-mono shrink-0">{sr.res_name}</Badge>
                            <span className="text-xs font-mono" title="PDB position">{sr.res_num}</span>
                            {uniProtPos !== null && (
                              <span className="text-xs font-mono text-emerald-600" title="UniProt position">{uniProtPos}</span>
                            )}
                            <Badge variant="secondary" className="text-xs bg-gray-50 shrink-0">Chain {sr.chain}</Badge>
                            <button
                              onClick={() => setSelectedResidues((prev) => prev.map((r, ri) => ri === i ? { ...r, part: r.part === "side_chain" ? "main_chain" : "side_chain" } : r))}
                              className={`text-xs px-1.5 py-0.5 rounded border ${sr.part === "side_chain" ? "border-emerald-200 text-emerald-700 bg-emerald-50" : "border-gray-200 text-gray-500 bg-gray-50"}`}>
                              {sr.part === "side_chain" ? "SC" : "MC"}
                            </button>
                            <button onClick={() => setSelectedResidues((prev) => prev.filter((_, ri) => ri !== i))} className="ml-auto text-gray-400 hover:text-red-500 shrink-0">
                              <Trash2 className="w-3 h-3" />
                            </button>
                          </div>
                          );
                        })}
                      </div>
                    )}
                  </div>
                </div>

                <div className="flex gap-2 flex-wrap">
                  <Button variant="outline" size="sm" className="h-7 text-xs border-gray-200 text-gray-500 hover:bg-gray-50 hover:text-gray-700"
                    onClick={goToDirectory}>
                    <FolderOpen className="w-3 h-3 mr-1" />Directory
                  </Button>
                  <Button variant="outline" onClick={() => goToStep(1)} className="border-gray-200 text-gray-600 hover:bg-gray-50 h-8 text-xs">
                    <ChevronLeft className="w-3 h-3 mr-1" />Step 1
                  </Button>
                  {selectedResidues.length > 0 && (
                    <Button
                      variant="outline"
                      onClick={() => viewerStep2Ref.current?.focusOnActiveCenter()}
                      className="border-amber-200 text-amber-600 hover:bg-amber-50 h-8 text-xs"
                    >
                      <Crosshair className="w-3 h-3 mr-1" />Show Active Center
                    </Button>
                  )}
                  <Button onClick={() => goToStep(3)} className="bg-violet-600 hover:bg-violet-700 text-white h-8 text-xs">
                    Proceed to Step 3 <ChevronRight className="w-3 h-3 ml-1" />
                  </Button>
                  {selectedResidues.length > 0 && (
                    <span className="text-xs text-emerald-600 flex items-center gap-1">
                      <CheckCircle2 className="w-3 h-3" />{selectedResidues.length} residue(s) ready
                    </span>
                  )}
                </div>
              </>
            )}
          </CardContent>
        </Card>
      )}

      {/* Step 3: Substrates & Cofactors */}
      {currentStep === 3 && (
        <Card className="border-violet-100 shadow-sm">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center gap-2 text-sm">
              <Beaker className="w-4 h-4 text-violet-600" />
              Step 3: Define Substrates & Cofactors
            </CardTitle>
            <CardDescription className="text-xs">
              Substrates and co-factors are treated in the same manner. Click &quot;+ Substrate/Cofactor&quot; to add a new molecule, then click the edit button to define it.
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            <Button onClick={() => setSubstrates(prev => [...prev, { id: genId(), name: "", smiles: "", type: "substrate", isEditing: false }])}
              variant="outline" className="border-violet-200 text-violet-600 hover:bg-violet-50 h-8 text-xs">
              <Plus className="w-3 h-3 mr-1" />+ Substrate/Cofactor
            </Button>

            {substrates.length === 0 ? (
              <div className="rounded-md border border-dashed border-gray-200 p-6 text-center">
                <p className="text-xs text-gray-400">No molecules added yet. Click the button above to add one.</p>
              </div>
            ) : (
              <div className="rounded-lg border border-gray-100 overflow-hidden">
                <Table>
                  <TableHeader>
                    <TableRow className="bg-gray-50/50">
                      <TableHead className="text-xs h-7 w-20">Type</TableHead>
                      <TableHead className="text-xs h-7">Name</TableHead>
                      <TableHead className="text-xs h-7">SMILES / Structure</TableHead>
                      <TableHead className="text-xs h-7 w-32">Mapped Ligand</TableHead>
                      <TableHead className="text-xs h-7 w-20"></TableHead>
                    </TableRow>
                  </TableHeader>
                  <TableBody>
                    {substrates.map((sub) => (
                      <TableRow key={sub.id} className="hover:bg-gray-50/50">
                        <TableCell className="py-2">
                          <select
                            value={sub.type}
                            onChange={(e) => setSubstrates(prev => prev.map(s => s.id === sub.id ? { ...s, type: e.target.value as SubstrateCofactor["type"] } : s))}
                            className={`text-xs border rounded px-1.5 py-0.5 bg-white h-6 ${
                              sub.type === "substrate" ? "border-emerald-200 text-emerald-700" :
                              sub.type === "cofactor" ? "border-violet-200 text-violet-600" : "border-cyan-200 text-cyan-700"
                            }`}
                          >
                            <option value="substrate">Substrate</option>
                            <option value="cofactor">Cofactor</option>
                            <option value="water">Water</option>
                          </select>
                        </TableCell>
                        <TableCell className="py-2">
                          <span className="text-xs text-gray-700 font-medium">{sub.name || <span className="text-gray-300 italic">unnamed</span>}</span>
                        </TableCell>
                        <TableCell className="py-2">
                          <div className="flex items-center gap-2">
                            {sub.smiles ? (
                              <>
                                <MoleculeSVG smiles={sub.smiles} size={32} />
                                <code className="text-xs font-mono text-gray-500 truncate max-w-[150px]">{sub.smiles}</code>
                              </>
                            ) : (
                              <span className="text-xs text-gray-300 italic">Not defined</span>
                            )}
                          </div>
                        </TableCell>
                        <TableCell className="py-2">
                          {sub.mapped_ligand ? (
                            <Badge variant="secondary" className="text-xs bg-amber-50 text-amber-700">{sub.mapped_ligand}</Badge>
                          ) : (
                            <span className="text-xs text-gray-300">None</span>
                          )}
                        </TableCell>
                        <TableCell className="py-2">
                          <div className="flex items-center gap-1">
                            <button
                              onClick={() => setSubstrates(prev => prev.map(s => s.id === sub.id ? { ...s, isEditing: !s.isEditing } : s))}
                              className="p-1 rounded hover:bg-violet-50 text-gray-400 hover:text-violet-600 transition-colors"
                              title="Edit molecule"
                            >
                              <PencilLucide className="w-3 h-3" />
                            </button>
                            <button
                              onClick={() => setSubstrates(prev => prev.filter(s => s.id !== sub.id))}
                              className="p-1 rounded hover:bg-red-50 text-gray-400 hover:text-red-500 transition-colors"
                              title="Remove"
                            >
                              <Trash2 className="w-3 h-3" />
                            </button>
                          </div>
                        </TableCell>
                      </TableRow>
                    ))}
                  </TableBody>
                </Table>
              </div>
            )}

            {/* Edit Dialog for Substrate/Cofactor */}
            {substrates.some(s => s.isEditing) && (() => {
              const editingSub = substrates.find(s => s.isEditing);
              if (!editingSub) return null;
              return (
                <EditMoleculeDialog
                  substrate={editingSub}
                  pdbLigands={pdbInfo?.ligands || []}
                  pdbInfo={pdbInfo}
                  selectedChain={selectedChain}
                  selectedResidues={selectedResidues}
                  pdbRawText={pdbRawText || ""}
                  ligandDistances={ligandDistances}
                  onSaveMol={(smiles) => setSubstrates(prev => prev.map(s => s.id === editingSub.id ? { ...s, smiles, isEditing: false } : s))}
                  onSaveName={(name) => setSubstrates(prev => prev.map(s => s.id === editingSub.id ? { ...s, name } : s))}
                  onSaveLigand={(ligand) => setSubstrates(prev => prev.map(s => s.id === editingSub.id ? { ...s, mapped_ligand: ligand } : s))}
                  onOverlay={handleOverlayResult}
                  onClose={() => setSubstrates(prev => prev.map(s => s.id === editingSub.id ? { ...s, isEditing: false } : s))}
                  toast={toast}
                />
              );
            })()}

            {/* Overlay Results Summary */}
            {overlayResults.length > 0 && (
              <div className="rounded-lg border border-violet-100 bg-violet-50/30 p-3 space-y-2">
                <div className="flex items-center gap-2">
                  <Layers className="w-3.5 h-3.5 text-violet-600" />
                  <span className="text-xs font-semibold text-violet-700">3D Overlay Results</span>
                  <Badge variant="secondary" className="text-xs bg-violet-100 text-violet-600 ml-auto">
                    {overlayResults.length} molecule(s)
                  </Badge>
                  <button
                    onClick={() => setOverlayResults([])}
                    className="p-0.5 rounded hover:bg-red-100 text-gray-400 hover:text-red-500 ml-1"
                    title="Clear all overlays"
                  >
                    <X className="w-3 h-3" />
                  </button>
                </div>
                <div className="grid grid-cols-1 sm:grid-cols-2 gap-2">
                  {overlayResults.map((or, i) => {
                    const sub = substrates.find(s => s.id === or.substrateId);
                    const color = OVERLAY_COLORS[i % OVERLAY_COLORS.length];
                    return (
                      <div key={or.substrateId} className="flex items-center gap-2 p-2 rounded-md bg-white border border-gray-100">
                        <div className="w-3 h-3 rounded-full shrink-0" style={{ backgroundColor: color }} />
                        <div className="flex-1 min-w-0">
                          <div className="text-xs font-medium truncate">
                            {sub?.name || `Molecule ${i + 1}`}
                          </div>
                          <div className="text-xs text-gray-400 font-mono">
                            {or.resName} {or.chain}:{or.resNum} · RMSD: {or.rmsd.toFixed(2)}Å · {or.numMapped} atoms
                          </div>
                        </div>
                        <Badge variant="secondary" className={`text-xs shrink-0 ${or.rmsd < 1.0 ? "bg-emerald-100 text-emerald-700" : or.rmsd < 2.0 ? "bg-amber-100 text-amber-700" : "bg-red-100 text-red-700"}`}>
                          {or.rmsd < 1.0 ? "Good" : or.rmsd < 2.0 ? "OK" : "Poor"}
                        </Badge>
                      </div>
                    );
                  })}
                </div>
                <p className="text-xs text-gray-400 leading-relaxed">
                  Overlaid molecules are shown in the 3D viewer on Step 1 and Step 2. Use the "Overlay" toggle button to show/hide them.
                </p>
              </div>
            )}

            <div className="flex gap-2">
              <Button variant="outline" size="sm" className="h-7 text-xs border-gray-200 text-gray-500 hover:bg-gray-50 hover:text-gray-700"
                onClick={goToDirectory}>
                <FolderOpen className="w-3 h-3 mr-1" />Directory
              </Button>
              <Button variant="outline" onClick={() => goToStep(2)} className="border-gray-200 text-gray-600 hover:bg-gray-50 h-8 text-xs">
                <ChevronLeft className="w-3 h-3 mr-1" />Step 2
              </Button>
              <Button onClick={() => goToStep(4)} className="bg-violet-600 hover:bg-violet-700 text-white h-8 text-xs">
                Proceed to Step 4 <ChevronRight className="w-3 h-3 ml-1" />
              </Button>
              {substrates.filter(s => s.smiles).length > 0 && (
                <span className="text-xs text-emerald-600 flex items-center gap-1 ml-2">
                  <CheckCircle2 className="w-3 h-3" />{substrates.filter(s => s.smiles).length} molecule(s) ready
                </span>
              )}
            </div>
          </CardContent>
        </Card>
      )}

      {/* Step 4: Active Site Configuration & Overall Reaction */}
      {currentStep === 4 && (
        <Card className="border-violet-100 shadow-sm">
          <CardHeader className="pb-1">
            <div className="flex items-center justify-between flex-wrap gap-2">
              <div>
                <CardTitle className="flex items-center gap-2 text-sm"><ArrowRight className="w-4 h-4 text-violet-600" />Step 4: Active Site Configuration & Reaction</CardTitle>
              </div>
              <div className="flex items-center gap-2">
                <Button onClick={handleNewSchemeDraft} variant="outline" className="border-violet-200 text-violet-600 hover:bg-violet-50 h-8 text-xs">
                  <Sparkles className="w-3 h-3 mr-1" />Re-fill from Step 3
                </Button>
              </div>
            </div>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* === ACTIVE SITE CONTEXT BAR === */}
            <div className="rounded-xl border border-amber-100 bg-gradient-to-r from-amber-50/60 via-white to-violet-50/60 p-4 space-y-3">
              <div className="flex items-center gap-2 mb-1">
                <Dna className="w-4 h-4 text-amber-600" />
                <span className="text-xs font-bold text-gray-700">Catalytic Residues</span>
                <span className="text-xs text-gray-400 ml-1">← Selected in Step 2</span>
              </div>
              {selectedResidues.length > 0 ? (
                <>
                  <div className="flex flex-wrap gap-2">
                    {selectedResidues.map((sr, i) => {
                      const aa = AMINO_ACID_SIDECHAIN[sr.res_name.toUpperCase()];
                      const resLabel = `${sr.res_name.charAt(0).toUpperCase()}${sr.res_name.slice(1).toLowerCase()}${sr.res_num}${sr.chain}`;
                      // Label for the image: "Ser70A—OH" format
                      const imgLabel = aa ? `${resLabel}—${aa.sideFormula}` : resLabel;
                      return (
                        <div key={i} className="flex flex-col items-center gap-1.5 rounded-lg border border-amber-200 bg-white p-2 cursor-pointer hover:border-amber-300 hover:shadow-md transition-all"
                          onClick={() => {
                            if (!aa) return;
                            setZoomModalOpen(true);
                            setZoomModalData(aa.sideChain);
                            setZoomModalLabel(imgLabel);
                            setZoomModalAtomLabel(resLabel);
                          }}>
                          {aa && <MoleculeSVG smiles={aa.sideChain} label={resLabel} idx={0} size={100} />}
                          <span className="text-xs font-mono text-amber-900 font-semibold">{imgLabel}</span>
                        </div>
                      );
                    })}
                  </div>
                </>
              ) : (
                <p className="text-xs text-amber-500 italic">No residues selected — go back to Step 2</p>
              )}
            </div>

            {/* === SINGLE REACTION EDITOR === */}
            <div className="space-y-2">
              <div className="flex items-center justify-between flex-wrap gap-2">
                <div className="flex items-center gap-2">
                  <div className="w-2.5 h-2.5 rounded-full bg-teal-400" />
                  <span className="text-xs font-semibold text-teal-700">Reaction Editor</span>
                  <span className="text-xs text-gray-400">(draw reactants → products)</span>
                </div>
                <div className="flex items-center gap-1">
                  <Button onClick={handleAutomap} disabled={automapping} variant="outline"
                    size="sm" className="h-7 text-xs border-teal-200 text-teal-600 hover:bg-teal-50 hover:text-teal-800"
                    title="Automatically assign atom map numbers using Indigo (also auto-called when running Step 5)">
                    {automapping ? <Loader2 className="w-3 h-3 mr-1 animate-spin" /> : <Zap className="w-3 h-3 mr-1" />}
                    Auto Map
                  </Button>
                  <Button onClick={handleVerifyMapping} disabled={verifying} variant="outline"
                    size="sm" className="h-7 text-xs border-amber-200 text-amber-600 hover:bg-amber-50 hover:text-amber-800"
                    title="Verify atom mapping quality using MCS — corrects Indigo errors like swapping carboxyl =O with water">
                    {verifying ? <Loader2 className="w-3 h-3 mr-1 animate-spin" /> : <Search className="w-3 h-3 mr-1" />}
                    Verify
                  </Button>
                  <Button onClick={handleCopyFromReactants} variant="outline"
                    size="sm" className="h-7 text-xs border-violet-200 text-violet-600 hover:bg-violet-50 hover:text-violet-800"
                    title="Copy reactant structures to the product side (clears existing products)">
                    <Copy className="w-3 h-3 mr-1" />Copy to Products
                  </Button>
                  <Button onClick={handleNewSchemeDraft} variant="outline"
                    size="sm" className="h-7 text-xs border-gray-200 text-gray-600 hover:bg-gray-50">
                    <Plus className="w-3 h-3 mr-1" />New Scheme
                  </Button>
                </div>
              </div>
              <ChemicalEditor
                ref={reactionEditorRef}
                height={520}
                onChange={handleReactionEditorChange}
              />
              <div className="flex items-center gap-2">
                <Label className="text-xs text-gray-400 uppercase tracking-wider shrink-0">Reaction SMILES:</Label>
                <Input
                  placeholder="e.g. CCO.N>>CCN.O"
                  value={reactantSmiles && productSmilesWf ? `${reactantSmiles}>>${productSmilesWf}` : reactantSmiles ? `${reactantSmiles}>>` : ""}
                  onChange={(e) => {
                    const val = e.target.value;
                    const parts = val.split(">>");
                    setReactantSmiles((parts[0] || "").trim());
                    setProductSmilesWf(parts.length > 1 ? (parts[1] || "").trim() : "");
                  }}
                  onBlur={() => {
                    const full = reactantSmiles && productSmilesWf ? `${reactantSmiles}>>${productSmilesWf}` : reactantSmiles ? `${reactantSmiles}>>` : "";
                    if (full.trim()) reactionEditorRef.current?.setMolecule(full.trim());
                  }}
                  className="font-mono text-xs h-8 border-gray-200"
                />
                {reactantSmiles && (
                  <MoleculeSVG smiles={reactantSmiles.split(".")[0]} size={28} />
                )}
              </div>
            </div>

            {/* Reaction Summary */}
            {reactantSmiles && productSmilesWf && (
              <div className="rounded-lg border border-gray-100 bg-gray-50/50 p-3 space-y-2">
                <div className="flex items-center justify-between">
                  <span className="text-xs font-medium text-gray-600">Reaction Summary</span>
                  <div className="flex items-center gap-1">
                    {selectedResidues.length > 0 && (
                      <Badge variant="secondary" className="bg-amber-50 text-amber-700 text-xs">
                        <Dna className="w-2 h-2 mr-0.5" />{selectedResidues.length} residues
                      </Badge>
                    )}
                    {pdbInfo?.pdb_id && (
                      <Badge variant="secondary" className="bg-violet-50 text-violet-700 text-xs">
                        <FileText className="w-2 h-2 mr-0.5" />{pdbInfo.pdb_id}
                      </Badge>
                    )}
                  </div>
                </div>
                <div className="flex items-center gap-3">
                  <div className="flex items-center gap-1.5 flex-1 min-w-0">
                    <Badge variant="secondary" className="bg-emerald-50 text-emerald-700 text-xs shrink-0">
                      {reactantSmiles.split(".").filter(Boolean).length} reactant(s)
                    </Badge>
                    <code className="text-xs font-mono text-gray-500 truncate">{reactantSmiles.length > 30 ? reactantSmiles.slice(0, 30) + "..." : reactantSmiles}</code>
                  </div>
                  <ArrowRight className="w-4 h-4 text-gray-400 shrink-0" />
                  <div className="flex items-center gap-1.5 flex-1 min-w-0">
                    <Badge variant="secondary" className="bg-violet-50 text-violet-700 text-xs shrink-0">
                      {productSmilesWf.split(".").filter(Boolean).length} product(s)
                    </Badge>
                    <code className="text-xs font-mono text-gray-500 truncate">{productSmilesWf.length > 30 ? productSmilesWf.slice(0, 30) + "..." : productSmilesWf}</code>
                  </div>
                </div>
                <p className="text-xs text-emerald-600 leading-relaxed">
                  The search engine will use all {reactantSmiles.split(".").filter(Boolean).length} reactant molecule(s) + {selectedResidues.length} catalytic residue(s) to match M-CSA rules.
                  Rules involving similar residue types will receive higher priority scores.
                </p>
              </div>
            )}

            {/* === OVERALL REACTION PREVIEW === */}
            {reactantSmiles && (
              <div className="rounded-xl border border-violet-200 bg-gradient-to-br from-violet-50/40 via-white to-emerald-50/40 p-5 space-y-3">
                <div className="flex items-center gap-2">
                  <FlaskConical className="w-4 h-4 text-violet-600" />
                  <span className="text-xs font-bold text-gray-700">Overall Reaction Preview</span>
                  {pdbInfo?.pdb_id && (
                    <Badge variant="outline" className="text-xs border-violet-200 text-violet-600 font-mono ml-auto">{pdbInfo.pdb_id}</Badge>
                  )}
                  <Button variant="ghost" size="sm" className="ml-2 h-7 w-7 p-0 text-gray-400 hover:text-violet-600 hover:bg-violet-50 shrink-0" title="Zoom reaction preview"
                    onClick={() => setRxnPreviewZoom(true)}>
                    <ZoomIn className="w-3.5 h-3.5" />
                  </Button>
                </div>
                {/* 2D molecule diagrams: ALL Reactants → enzyme → ALL Products */}
                <div className="flex items-center gap-2 flex-wrap">
                  {/* Reactants */}
                  <div className="flex-1 min-w-[200px] space-y-1">
                    <span className="text-xs font-semibold text-emerald-700 block">Reactants ({reactantSmiles.split(".").filter(Boolean).length})</span>
                    <div className="rounded-lg border border-emerald-100 bg-white p-3 flex flex-wrap items-center justify-center gap-3 min-h-[180px]">
                      {reactantSmiles.split(".").filter(Boolean).map((smi, i) => (
                        <div key={i} className="flex flex-col items-center gap-1">
                          <MoleculeSVG smiles={smi.trim()} size={reactantSmiles.split(".").filter(Boolean).length === 1 ? 240 : 200} />
                          <span className="text-xs text-gray-500 font-medium">R{i + 1}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                  {/* enzyme arrow */}
                  <div className="flex flex-col items-center gap-1 shrink-0 py-4">
                    <div className="rounded-full border-2 border-dashed border-violet-300 bg-violet-50/50 p-1.5">
                      {selectedResidues.length > 0 ? (
                        <span className="text-xs font-bold text-violet-600">{selectedResidues.length} res</span>
                      ) : (
                        <ArrowRight className="w-5 h-5 text-violet-400" />
                      )}
                    </div>
                    <span className="text-xs text-violet-400">enzyme</span>
                  </div>
                  {/* Products */}
                  <div className="flex-1 min-w-[200px] space-y-1">
                    <span className="text-xs font-semibold text-violet-700 block">Products ({productSmilesWf ? productSmilesWf.split(".").filter(Boolean).length : 0})</span>
                    <div className="rounded-lg border border-violet-100 bg-white p-3 flex flex-wrap items-center justify-center gap-3 min-h-[180px]">
                      {productSmilesWf ? (
                        productSmilesWf.split(".").filter(Boolean).map((smi, i) => (
                          <div key={i} className="flex flex-col items-center gap-1">
                            <MoleculeSVG smiles={smi.trim()} size={productSmilesWf.split(".").filter(Boolean).length === 1 ? 240 : 200} />
                            <span className="text-xs text-gray-500 font-medium">P{i + 1}</span>
                          </div>
                        ))
                      ) : (
                        <span className="text-xs text-gray-300 italic">Draw product in editor →</span>
                      )}
                    </div>
                  </div>
                </div>
                {/* Catalytic Residues Annotation Bar — with R-Group style side chain thumbnails */}
                {selectedResidues.length > 0 && (
                  <div className="rounded-lg border border-amber-100 bg-gradient-to-r from-amber-50/30 to-white p-3 space-y-2">
                    <div className="flex items-center gap-1.5">
                      <Dna className="w-3.5 h-3.5 text-amber-600" />
                      <span className="text-xs font-semibold text-amber-800">Catalytic Residues ({selectedResidues.length})</span>
                    </div>
                    <div className="flex flex-wrap gap-2">
                      {selectedResidues.map((r, i) => {
                        const aa = AMINO_ACID_SIDECHAIN[r.res_name.toUpperCase()];
                        const resLabel = `${r.res_name.charAt(0).toUpperCase()}${r.res_name.slice(1).toLowerCase()}${r.res_num}${r.chain}`;
                        const imgLabel = aa ? `${resLabel}—${aa.sideFormula}` : resLabel;
                        return (
                          <div key={i} className="flex flex-col items-center gap-1 rounded-lg border border-amber-200 bg-white p-1.5 cursor-pointer hover:border-amber-300 hover:shadow-sm transition-all"
                            onClick={() => {
                              if (!aa) return;
                              setZoomModalOpen(true);
                              setZoomModalData(aa.sideChain);
                              setZoomModalLabel(imgLabel);
                              setZoomModalAtomLabel(resLabel);
                            }}>
                            {aa && <MoleculeSVG smiles={aa.sideChain} label={resLabel} idx={0} size={100} />}
                            <span className="text-xs font-mono text-amber-900 font-semibold">{imgLabel}</span>
                          </div>
                        );
                      })}
                    </div>
                  </div>
                )}
              </div>
            )}

            {/* Residue Zoom Modal */}
            <ZoomModal open={zoomModalOpen} onClose={() => setZoomModalOpen(false)} smiles={zoomModalData} title={zoomModalLabel} label={zoomModalAtomLabel} idx={0} />

            {/* Reaction Preview Zoom Modal — 2x enlarged Reactants + Products */}
            <Dialog open={rxnPreviewZoom} onOpenChange={(o) => !o && setRxnPreviewZoom(false)}>
              <DialogContent className="!max-w-[1224px] w-full p-4 space-y-4">
                <DialogHeader>
                  <DialogTitle className="text-sm font-semibold text-gray-700 flex items-center gap-2">
                    <FlaskConical className="w-4 h-4 text-violet-600" />
                    Overall Reaction Preview
                    {pdbInfo?.pdb_id && (
                      <Badge variant="outline" className="text-xs border-violet-200 text-violet-600 font-mono">{pdbInfo.pdb_id}</Badge>
                    )}
                  </DialogTitle>
                </DialogHeader>
                <div className="flex items-center gap-3 flex-wrap">
                  {/* Reactants */}
                  <div className="flex-1 min-w-[220px] space-y-1">
                    <span className="text-xs font-semibold text-emerald-700 block">Reactants ({reactantSmiles.split(".").filter(Boolean).length})</span>
                    <div className="rounded-lg border border-emerald-100 bg-white p-4 flex flex-wrap items-center justify-center gap-4 min-h-[280px]">
                      {reactantSmiles.split(".").filter(Boolean).map((smi, i) => (
                        <div key={i} className="flex flex-col items-center gap-1">
                          <MoleculeSVG smiles={smi.trim()} size={reactantSmiles.split(".").filter(Boolean).length === 1 ? 480 : 400} />
                          <span className="text-sm text-gray-500 font-medium">R{i + 1}</span>
                        </div>
                      ))}
                    </div>
                  </div>
                  {/* enzyme arrow */}
                  <div className="flex flex-col items-center gap-1 shrink-0 py-4">
                    <div className="rounded-full border-2 border-dashed border-violet-300 bg-violet-50/50 p-2">
                      {selectedResidues.length > 0 ? (
                        <span className="text-sm font-bold text-violet-600">{selectedResidues.length} res</span>
                      ) : (
                        <ArrowRight className="w-7 h-7 text-violet-400" />
                      )}
                    </div>
                    <span className="text-xs text-violet-400">enzyme</span>
                  </div>
                  {/* Products */}
                  <div className="flex-1 min-w-[220px] space-y-1">
                    <span className="text-xs font-semibold text-violet-700 block">Products ({productSmilesWf ? productSmilesWf.split(".").filter(Boolean).length : 0})</span>
                    <div className="rounded-lg border border-violet-100 bg-white p-4 flex flex-wrap items-center justify-center gap-4 min-h-[280px]">
                      {productSmilesWf ? (
                        productSmilesWf.split(".").filter(Boolean).map((smi, i) => (
                          <div key={i} className="flex flex-col items-center gap-1">
                            <MoleculeSVG smiles={smi.trim()} size={productSmilesWf.split(".").filter(Boolean).length === 1 ? 480 : 400} />
                            <span className="text-sm text-gray-500 font-medium">P{i + 1}</span>
                          </div>
                        ))
                      ) : (
                        <span className="text-sm text-gray-300 italic">No products defined</span>
                      )}
                    </div>
                  </div>
                </div>
              </DialogContent>
            </Dialog>

            {/* No product warning */}
            {!productSmilesWf.trim() && (
              <Alert className="border-amber-200 bg-amber-50/50">
                <AlertCircle className="h-3.5 w-3.5 text-amber-500" />
                <AlertTitle className="text-xs text-amber-800">Products Not Yet Defined</AlertTitle>
                <AlertDescription className="text-xs text-amber-700">
                  Without products, only forward prediction (not bidirectional search) will run in Step 5.
                  Define the expected products for the most accurate mechanism prediction.
                </AlertDescription>
              </Alert>
            )}

            <div className="flex gap-2">
              <Button variant="outline" size="sm" className="h-7 text-xs border-gray-200 text-gray-500 hover:bg-gray-50 hover:text-gray-700"
                onClick={goToDirectory}>
                <FolderOpen className="w-3 h-3 mr-1" />Directory
              </Button>
              <Button variant="outline" onClick={() => goToStep(3)} className="border-gray-200 text-gray-600 hover:bg-gray-50 h-8 text-xs">
                <ChevronLeft className="w-3 h-3 mr-1" />Step 3
              </Button>
              <Button onClick={() => goToStep(5)} disabled={!reactantSmiles.trim()} className="bg-violet-600 hover:bg-violet-700 text-white h-8 text-xs">
                Proceed to Step 5 <ChevronRight className="w-3 h-3 ml-1" />
              </Button>
            </div>
          </CardContent>
        </Card>
      )}

      {/* Step 5: Parameters & Run */}
      {currentStep === 5 && (
        <Card className="border-violet-100 shadow-sm">
          <CardHeader className="pb-3">
            <CardTitle className="flex items-center gap-2 text-sm"><Sparkles className="w-4 h-4 text-violet-600" />Step 5: Parameters & Run</CardTitle>
            <CardDescription className="text-xs">Configure search parameters and run mechanism prediction</CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Context Summary — Full Active Site Configuration */}
            <div className="rounded-lg border border-gray-100 bg-gray-50/50 p-3 space-y-3">
              <Label className="text-xs text-gray-400 uppercase tracking-wider">Active Site Configuration Summary</Label>
              <div className="grid grid-cols-2 sm:grid-cols-5 gap-2">
                <div className="rounded-md bg-white p-2 text-center">
                  <div className="text-sm font-bold font-mono text-violet-700">{pdbInfo?.pdb_id || "N/A"}</div>
                  <div className="text-xs text-gray-500">PDB</div>
                </div>
                <div className="rounded-md bg-white p-2 text-center">
                  <div className="text-sm font-bold font-mono text-amber-700">{selectedResidues.length}</div>
                  <div className="text-xs text-gray-500">Residues</div>
                </div>
                <div className="rounded-md bg-white p-2 text-center">
                  <div className="text-sm font-bold font-mono text-emerald-700">{substrates.length}</div>
                  <div className="text-xs text-gray-500">Molecules</div>
                </div>
                <div className="rounded-md bg-white p-2 text-center">
                  <div className="text-sm font-bold font-mono text-emerald-700">{reactantSmiles.split(".").filter(Boolean).length}</div>
                  <div className="text-xs text-gray-500">Reactants</div>
                </div>
                <div className="rounded-md bg-white p-2 text-center">
                  <div className="text-sm font-bold font-mono text-violet-700">{productSmilesWf.trim() ? productSmilesWf.split(".").filter(Boolean).length : 0}</div>
                  <div className="text-xs text-gray-500">Products</div>
                </div>
              </div>
              {/* Visual flow diagram */}
              <div className="flex items-center gap-1.5 text-xs pt-1 border-t border-gray-100">
                <Badge variant="outline" className="text-xs border-violet-200 text-violet-600"><FileText className="w-2 h-2 mr-0.5" />PDB</Badge>
                <ChevronRight className="w-3 h-3 text-gray-300" />
                <Badge variant="outline" className="text-xs border-amber-200 text-amber-600"><Dna className="w-2 h-2 mr-0.5" />Residues</Badge>
                <ChevronRight className="w-3 h-3 text-gray-300" />
                <Badge variant="outline" className="text-xs border-emerald-200 text-emerald-600"><Beaker className="w-2 h-2 mr-0.5" />Molecules</Badge>
                <ChevronRight className="w-3 h-3 text-gray-300" />
                <Badge variant="outline" className="text-xs border-violet-200 text-violet-600"><ArrowRight className="w-2 h-2 mr-0.5" />Reaction</Badge>
                <ChevronRight className="w-3 h-3 text-gray-300" />
                <Badge variant="secondary" className="text-xs bg-violet-100 text-violet-700"><Sparkles className="w-2 h-2 mr-0.5" />Search</Badge>
              </div>
              {selectedResidues.length > 0 && (
                <p className="text-xs text-emerald-600">
                  ✓ Residue-aware scoring enabled: {selectedResidues.map(r => r.res_name).join(", ")} residues will boost matching M-CSA rules
                </p>
              )}
            </div>

            <div className="grid grid-cols-1 sm:grid-cols-2 gap-4">
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label className="text-xs font-medium">Max Configurations</Label>
                  <Badge variant="secondary" className="bg-violet-50 text-violet-700 font-mono text-xs">{wfMaxConfigs}</Badge>
                </div>
                <Slider value={[wfMaxConfigs]} onValueChange={(v) => setWfMaxConfigs(v[0])} min={50} max={2000} step={50}
                  className="[&_[data-slot=slider-range]]:bg-violet-500 [&_[data-slot=slider-thumb]]:border-violet-500" />
              </div>
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label className="text-xs font-medium">Max Bond Length (Å)</Label>
                  <Badge variant="secondary" className="bg-violet-50 text-violet-700 font-mono text-xs">{wfMaxBondLength}</Badge>
                </div>
                <Slider value={[wfMaxBondLength]} onValueChange={(v) => setWfMaxBondLength(v[0])} min={0} max={11} step={1}
                  className="[&_[data-slot=slider-range]]:bg-violet-500 [&_[data-slot=slider-thumb]]:border-violet-500" />
              </div>
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label className="text-xs font-medium">M-CSA ID (Own Rules)</Label>
                  <Badge variant="secondary" className="bg-violet-50 text-violet-700 font-mono text-xs">{wfMcsaId || 'All'}</Badge>
                </div>
                <Slider value={[wfMcsaId]} onValueChange={(v) => setWfMcsaId(v[0])} min={0} max={1000} step={1}
                  className="[&_[data-slot=slider-range]]:bg-violet-500 [&_[data-slot=slider-thumb]]:border-violet-500" />
                <p className="text-[10px] text-gray-400">0 = use all rules. Set to M-CSA enzyme ID to filter rules for that enzyme only.</p>
              </div>
            </div>

            <div className="space-y-1.5">
              <Label className="text-xs font-medium">Comment (optional)</Label>
              <Textarea placeholder="Add notes about this prediction..." value={wfComment} onChange={(e) => setWfComment(e.target.value)}
                className="text-xs border-gray-200 min-h-[60px] resize-y" />
            </div>

            <div className="flex gap-2">
              <Button variant="outline" size="sm" className="h-7 text-xs border-gray-200 text-gray-500 hover:bg-gray-50 hover:text-gray-700"
                onClick={goToDirectory}>
                <FolderOpen className="w-3 h-3 mr-1" />Directory
              </Button>
              <Button variant="outline" onClick={() => goToStep(4)} className="border-gray-200 text-gray-600 hover:bg-gray-50 h-9 text-xs">
                <ChevronLeft className="w-3 h-3 mr-1" />Step 4
              </Button>
              <Button onClick={handleWfRun} disabled={wfRunning || !reactantSmiles.trim()}
                className="flex-1 bg-gradient-to-r from-violet-600 to-emerald-600 hover:from-violet-700 hover:to-emerald-700 text-white shadow-md h-9 text-sm">
                <Sparkles className="w-4 h-4 mr-2" />
                Run Mechanism Search
              </Button>
            </div>
            <p className="text-[11px] text-gray-400 text-center mt-1">
              Auto Map will be applied automatically before search to assign atom map numbers.
            </p>
          </CardContent>
        </Card>
      )}
    </div>
  );
}