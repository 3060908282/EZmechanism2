"use client";

import React, { useState, useEffect, useRef, useMemo, useCallback } from "react";
import {
  Pencil, Save, BookOpen, RotateCcw, CheckCircle2, Loader2, Filter,
  Layers, AlertTriangle,
} from "lucide-react";
import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Label as LabelUI } from "@/components/ui/label";
import { Badge } from "@/components/ui/badge";
import { Table, TableBody, TableCell, TableRow } from "@/components/ui/table";
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from "@/components/ui/dialog";
import { useToast } from "@/hooks/use-toast";
import ChemicalEditor, { type ChemicalEditorHandle } from "@/components/ChemicalEditor";
import { apiCall } from "@/lib/api";
import type { SubstrateCofactor, PdbLigand, PdbInfo, SelectedResidue } from "@/lib/types";

export interface OverlayResult {
  substrateId: string;
 smiles: string;
 molblock: string;
 rmsd: number;
 mcsNumAtoms: number;
 numMapped: number;
 pdbLigandSmiles: string | null;
 resName: string;
 chain: string;
  resNum: number;
 message: string;
 isWater?: boolean;
}

export default function EditMoleculeDialog({
  substrate,
  pdbLigands,
  pdbInfo,
  selectedChain,
  selectedResidues,
  pdbRawText,
  ligandDistances,
  onSaveMol,
  onSaveName,
  onSaveLigand,
  onOverlay,
  onClose,
  toast,
}: {
  substrate: SubstrateCofactor;
  pdbLigands: PdbLigand[];
  pdbInfo: PdbInfo | null;
  selectedChain: string | null;
  selectedResidues: SelectedResidue[];
  pdbRawText: string;
  ligandDistances?: Record<string, number>;
  onSaveMol: (smiles: string) => void;
  onSaveName: (name: string) => void;
  onSaveLigand: (ligand: string) => void;
  onOverlay?: (result: OverlayResult) => void;
  onClose: () => void;
  toast: ReturnType<typeof useToast>["toast"];
}) {
  const editorRef = useRef<ChemicalEditorHandle>(null);
  const [editorReady, setEditorReady] = useState(false);
  const [smilesText, setSmilesText] = useState(substrate.smiles || "");
  const [chebiId, setChebiId] = useState("");
  const [nameText, setNameText] = useState(substrate.name || "");
  const [selectedLigand, setSelectedLigand] = useState(substrate.mapped_ligand || "");
  const [chebiLoading, setChebiLoading] = useState(false);
  const [mcsInfo, setMcsInfo] = useState<string | null>(null);
  const [mcsSvgA, setMcsSvgA] = useState<string | null>(null);
  const [mcsSvgB, setMcsSvgB] = useState<string | null>(null);
  const [mcsStats, setMcsStats] = useState<{similarity: number; mcs_atoms: number; mcs_bonds: number; mcs_smarts: string} | null>(null);
  const [computingSimilarity, setComputingSimilarity] = useState(false);
  const [computingMcs, setComputingMcs] = useState(false);
  const [overlayLoading, setOverlayLoading] = useState(false);
  const [overlayResult, setOverlayResult] = useState<OverlayResult | null>(null);
  const [chainFilter, setChainFilter] = useState<"selected" | "all">("selected");

  // Use refs for caches to avoid infinite re-render loops.
  // These are read/written inside async useEffects — putting them in state
  // and dependency arrays causes infinite loops.
  const ligandSmilesRef = useRef<Record<string, string>>({});
  const similarityRef = useRef<Record<string, number>>({});
  const mcsSmartsRef = useRef<Record<string, string>>({});
  const distanceRef = useRef<Record<string, number>>({});
  // A simple version counter to trigger useMemo recompute after async work finishes
  const [cacheVersion, setCacheVersion] = useState(0);

  // ---- Filter ligands by chain ----
  const filteredLigands = useMemo(() => {
    if (chainFilter === "all" || !selectedChain) return pdbLigands;
    return pdbLigands.filter(l => l.chain === selectedChain);
  }, [pdbLigands, chainFilter, selectedChain]);

  // ---- Compute distances to active site center of geometry ----
  useEffect(() => {
    if (filteredLigands.length === 0) return;
    // Use real distances from backend API if available
    if (ligandDistances && Object.keys(ligandDistances).length > 0) {
      for (const lig of filteredLigands) {
        const key = `${lig.res_name}_${lig.chain}_${lig.res_num}`;
        if (ligandDistances[key] !== undefined) {
          distanceRef.current[key] = ligandDistances[key];
        }
      }
    } else {
      // Fallback: rough estimation based on chain + coordinate presence
      for (const lig of filteredLigands) {
        const key = `${lig.res_name}_${lig.chain}_${lig.res_num}`;
        if (distanceRef.current[key] !== undefined) continue;
        if (selectedChain && lig.chain === selectedChain) {
          distanceRef.current[key] = 1.0;
        } else {
          distanceRef.current[key] = 10.0;
        }
        if (lig.x !== undefined && lig.y !== undefined && lig.z !== undefined) {
          const hasValidCoords = lig.x !== 0 || lig.y !== 0 || lig.z !== 0;
          if (hasValidCoords) distanceRef.current[key] = Math.min(distanceRef.current[key], 0.5);
        }
      }
    }
  }, [filteredLigands, selectedChain, ligandDistances]);

  // ---- Compute similarities for all filtered ligands when smilesText changes ----
  useEffect(() => {
    if (!smilesText.trim() || filteredLigands.length === 0) return;
    let cancelled = false;
    (async () => {
      setComputingSimilarity(true);
      const smiles = smilesText.trim();

      for (const lig of filteredLigands) {
        const key = `${lig.res_name}_${lig.chain}_${lig.res_num}`;
        // Reuse cached similarity
        if (similarityRef.current[key] !== undefined) continue;

        // Skip similarity computation for water molecules — O vs O is meaningless
        // and hundreds of HOH entries would cause excessive API calls
        if (lig.is_water) {
          similarityRef.current[key] = 0;
          continue;
        }

        // Try to get ligand SMILES from ref cache first, then API
        let ligSmi = ligandSmilesRef.current[`${lig.res_name}_${lig.chain}`] || ligandSmilesRef.current[lig.res_name];
        if (!ligSmi) {
          try {
            const r = await apiCall("ligand-pdb-smiles", { res_name: lig.res_name }, 5000);
            if (r.ok) {
              const d = await r.json();
              if (d.smiles) {
                ligSmi = d.smiles;
                ligandSmilesRef.current[lig.res_name] = ligSmi;
              }
            }
          } catch { /* skip */ }
          if (!ligSmi && pdbRawText) {
            try {
              const r2 = await apiCall("ligand-extract-smiles", {
                res_name: lig.res_name,
                chain: lig.chain,
                res_num: lig.res_num,
                pdb_text: pdbRawText,
              }, 10000);
              if (r2.ok) {
                const d2 = await r2.json();
                if (d2.smiles) {
                  ligSmi = d2.smiles;
                  ligandSmilesRef.current[`${lig.res_name}_${lig.chain}`] = ligSmi;
                }
              }
            } catch { /* skip */ }
          }
        }
        if (ligSmi) {
          try {
            const r = await apiCall("ligand-similarity", { smiles_a: smiles, smiles_b: ligSmi }, 10000);
            if (r.ok) {
              const d = await r.json();
              similarityRef.current[key] = d.similarity || 0;
              // Cache mcs_smarts from ligand-similarity for potential reuse by ligand-compare
              if (d.mcs_smarts) {
                mcsSmartsRef.current[key] = d.mcs_smarts;
              }
            }
          } catch { similarityRef.current[key] = 0; }
        } else {
          similarityRef.current[key] = 0;
        }
        if (cancelled) break;
      }

      if (!cancelled) {
        setComputingSimilarity(false);
        setCacheVersion(v => v + 1); // trigger useMemo recompute
      }
    })();
    return () => { cancelled = true; };
  }, [smilesText, filteredLigands, pdbRawText]);

  // ---- Sort ligands: non-water first (by similarity → proximity), then water (by proximity) ----
  const sortedLigands = useMemo(() => {
    if (filteredLigands.length === 0) return [];
    const isSubstrateWater = substrate.type === "water";
    const mapped = [...filteredLigands].map(l => {
      const key = `${l.res_name}_${l.chain}_${l.res_num}`;
      return {
        ...l,
        similarity: similarityRef.current[key] ?? l.similarity ?? 0,
        distance: distanceRef.current[key] ?? 999,
      };
    });

    // Split into water and non-water groups
    const waterLigs = mapped.filter(l => l.is_water);
    const nonWaterLigs = mapped.filter(l => !l.is_water);

    // Sort non-water by similarity (desc) → distance (asc)
    nonWaterLigs.sort((a, b) => {
      if (b.similarity !== a.similarity) return b.similarity - a.similarity;
      return a.distance - b.distance;
    });
    // Sort water by distance (asc) — proximity to active site is the only useful metric
    waterLigs.sort((a, b) => a.distance - b.distance);

    // If substrate is water, show water molecules first; otherwise non-water first
    if (isSubstrateWater) {
      return [...waterLigs, ...nonWaterLigs];
    }
    return [...nonWaterLigs, ...waterLigs];
  }, [filteredLigands, cacheVersion, substrate.type]);

  // ---- When a ligand is selected, compute MCS comparison ----
  useEffect(() => {
    if (!selectedLigand || !smilesText.trim()) {
      setMcsSvgA(null); setMcsSvgB(null); setMcsStats(null); setMcsInfo(null);
      return;
    }
    let cancelled = false;
    (async () => {
      setComputingMcs(true);
      const parts = selectedLigand.split("_");
      const resName = parts[0];
      const chain = parts[1];

      // Water molecules: skip SMILES lookup & MCS — HOH has no meaningful SMILES
      // and O-vs-O MCS comparison is meaningless
      const selectedLigObj = pdbLigands.find(l =>
        `${l.res_name}_${l.chain}_${l.res_num}` === selectedLigand
      );
      if (selectedLigObj?.is_water) {
        setMcsInfo("💧 Water molecule — coordinate mapping only (no MCS comparison needed)");
        setMcsSvgA(null); setMcsSvgB(null); setMcsStats(null);
        setComputingMcs(false);
        return;
      }

      // Try to get ligand SMILES from ref cache or API
      let ligSmi = ligandSmilesRef.current[`${resName}_${chain}`] || ligandSmilesRef.current[resName];
      if (!ligSmi) {
        try {
          const r = await apiCall("ligand-pdb-smiles", { res_name: resName }, 5000);
          if (r.ok) {
            const d = await r.json();
            if (d.smiles) {
              ligSmi = d.smiles;
              ligandSmilesRef.current[resName] = ligSmi;
            }
          }
        } catch { /* skip */ }
        if (!ligSmi && pdbRawText) {
          try {
            const r2 = await apiCall("ligand-extract-smiles", {
              res_name: resName,
              chain: chain,
              res_num: parseInt(parts[2], 10),
              pdb_text: pdbRawText,
            }, 10000);
            if (r2.ok) {
              const d2 = await r2.json();
              if (d2.smiles) {
                ligSmi = d2.smiles;
                ligandSmilesRef.current[`${resName}_${chain}`] = ligSmi;
              }
            }
          } catch { /* skip */ }
        }
      }

      if (ligSmi) {
        try {
          const r = await apiCall("ligand-compare", {
            smiles_a: smilesText.trim(), smiles_b: ligSmi, width: 380, height: 280,
          }, 15000);
          if (r.ok) {
            const d = await r.json();
            if (!cancelled) {
              setMcsSvgA(d.svg_a || null);
              setMcsSvgB(d.svg_b || null);
              setMcsStats({
                similarity: d.similarity || 0,
                mcs_atoms: d.mcs_num_atoms || 0,
                mcs_bonds: d.mcs_num_bonds || 0,
                mcs_smarts: d.mcs_smarts || "",
              });
              setMcsInfo(`MCS: ${d.mcs_num_atoms} atoms, ${d.mcs_num_bonds} bonds | Similarity: ${((d.similarity || 0) * 100).toFixed(1)}%`);
            }
          }
        } catch { setMcsInfo("Failed to compute MCS"); }
      } else {
        setMcsInfo(`No SMILES available for ${resName} — comparison not possible. The ligand may not have known chemical structure.`);
      }
      if (!cancelled) setComputingMcs(false);
    })();
    return () => { cancelled = true; };
  }, [selectedLigand, smilesText, pdbRawText]);

  // ---- ChEBI Import ----
  const handleFromChebi = useCallback(async () => {
    if (!chebiId.trim()) return;
    setChebiLoading(true);
    try {
      const r = await apiCall("chebi-fetch", { chebi_id: chebiId.trim() }, 15000);
      if (!r.ok) throw new Error("Failed");
      const data = await r.json();
      if (data.smiles) {
        setSmilesText(data.smiles);
        editorRef.current?.setMolecule(data.smiles);
        toast({ title: "ChEBI Loaded", description: `${data.name || chebiId}: ${data.smiles.slice(0, 40)}...` });
      } else {
        toast({ title: "No SMILES Found", description: `ChEBI:${chebiId} has no associated SMILES`, variant: "destructive" });
      }
    } catch {
      toast({ title: "ChEBI Fetch Failed", description: "Could not retrieve molecule from ChEBI", variant: "destructive" });
    } finally { setChebiLoading(false); }
  }, [chebiId, toast]);

  const handleSaveMol = useCallback(async () => {
    if (!smilesText.trim()) return;
    onSaveMol(smilesText.trim());
    toast({ title: "Molecule Saved", description: `SMILES: ${smilesText.slice(0, 50)}${smilesText.length > 50 ? "..." : ""}` });
  }, [smilesText, onSaveMol, toast]);

  const handleSaveName = useCallback(() => {
    if (!nameText.trim()) return;
    onSaveName(nameText.trim());
    toast({ title: "Name Saved", description: nameText });
  }, [nameText, onSaveName, toast]);

  const handleSaveLigand = useCallback(() => {
    onSaveLigand(selectedLigand);
    if (selectedLigand) {
      toast({ title: "Ligand Mapped", description: `Mapped to PDB ligand: ${selectedLigand}` });
    }
  }, [selectedLigand, onSaveLigand, toast]);

  // ---- Overlay: align hand-drawn molecule to PDB ligand 3D position ----
  const handleOverlay = useCallback(async () => {
    if (!smilesText.trim() || !selectedLigand || !pdbRawText) return;
    setOverlayLoading(true);
    setOverlayResult(null);
    try {
      const parts = selectedLigand.split("_");
      const resName = parts[0];
      const chain = parts[1];
      const resNum = parseInt(parts[2], 10);

      const r = await apiCall("overlay-substrate", {
        smiles: smilesText.trim(),
        res_name: resName,
        chain,
        res_num: resNum,
        pdb_text: pdbRawText,
      }, 30000);

      if (!r.ok) throw new Error(`Server error: ${r.status}`);
      const d = await r.json();

      if (!d.success) {
        toast({ title: "Overlay Failed", description: d.message || "Could not align molecule to PDB ligand", variant: "destructive" });
        setOverlayLoading(false);
        return;
      }

      const result: OverlayResult = {
        substrateId: substrate.id,
        smiles: smilesText.trim(),
        molblock: d.molblock,
        rmsd: d.rmsd,
        mcsNumAtoms: d.mcs_num_atoms,
        numMapped: d.num_mapped,
        pdbLigandSmiles: d.pdb_ligand_smiles || null,
        resName, chain, resNum,
        message: d.message,
        isWater: d.is_water || false,
      };
      setOverlayResult(result);

      // Save ligand mapping if not already done
      onSaveLigand(selectedLigand);

      // Notify parent of overlay result
      onOverlay?.(result);

      if (d.is_water) {
        toast({
          title: "Water Coordinate Mapping Applied",
          description: d.message || "Water O atom placed at PDB coordinates.",
        });
      } else {
        const rmsdStr = d.rmsd < 1.0 ? `${(d.rmsd * 100).toFixed(0)}cm` : `${d.rmsd.toFixed(2)}Å`;
        toast({
          title: "Overlay Applied",
          description: `${d.num_mapped} atoms aligned (RMSD: ${rmsdStr}). Check the 3D viewer to verify the overlay.`,
        });
      }
    } catch (err: unknown) {
      const msg = err instanceof Error ? err.message : "Overlay failed";
      toast({ title: "Overlay Error", description: msg, variant: "destructive" });
    } finally {
      setOverlayLoading(false);
    }
  }, [smilesText, selectedLigand, pdbRawText, substrate.id, onSaveLigand, onOverlay, toast]);

  const selectedLigObj = pdbLigands.find(l =>
    `${l.res_name}_${l.chain}_${l.res_num}` === selectedLigand
  );
  const isSelectedLigandWater = selectedLigObj?.is_water || false;

  const selectedChainLigands = pdbLigands.filter(l => l.chain === selectedChain);
  const allChains = useMemo(() => {
    const chains = new Set(pdbLigands.map(l => l.chain));
    return Array.from(chains).sort();
  }, [pdbLigands]);

  return (
    <Dialog open={true} onOpenChange={(o) => !o && onClose()}>
      <DialogContent className="sm:max-w-[90rem] max-w-[calc(100%-1rem)] w-[98vw] max-h-[92vh] overflow-y-auto p-0">
        <DialogHeader className="px-6 pt-4 pb-2 border-b">
          <DialogTitle className="flex items-center gap-2 text-sm">
            <Pencil className="w-4 h-4 text-violet-600" />
            Editing {substrate.type === "cofactor" ? "Cofactor" : substrate.type === "water" ? "Water" : "Substrate"}
            {substrate.name && <span className="text-gray-400">— {substrate.name}</span>}
          </DialogTitle>
        </DialogHeader>

        <div className="grid grid-cols-1 lg:grid-cols-5 gap-0 lg:divide-x">
          {/* LEFT SECTION: Molecule Definition (3/5) */}
          <div className="lg:col-span-3 p-4 space-y-3">
            <div className="flex items-center gap-2 mb-2">
              <span className="text-xs font-semibold text-violet-700 uppercase tracking-wider">Molecule Definition</span>
            </div>

            {/* Action buttons row */}
            <div className="flex items-center gap-2 flex-wrap">
              <Button onClick={async () => { const smi = await editorRef.current?.getSmiles(); if (smi) setSmilesText(smi); }}
                variant="outline" size="sm" className="h-7 text-xs border-violet-200 text-violet-600 hover:bg-violet-50">
                <RotateCcw className="w-3 h-3 mr-1" />Reload Mol
              </Button>
              <Button onClick={handleSaveMol} disabled={!smilesText.trim()}
                size="sm" className="h-7 text-xs bg-violet-600 hover:bg-violet-700 text-white">
                <Save className="w-3 h-3 mr-1" />Save Mol
              </Button>
              <div className="flex-1" />
              <div className="flex items-center gap-1">
                <Input placeholder="ChEBI ID (e.g. 15377 for L-Alanine)" value={chebiId}
                  onChange={(e) => setChebiId(e.target.value)}
                  onKeyDown={(e) => e.key === "Enter" && handleFromChebi()}
                  className="text-xs h-7 w-72 font-mono" />
                <Button onClick={handleFromChebi} disabled={chebiLoading || !chebiId.trim()}
                  size="sm" variant="outline" className="h-7 text-xs border-teal-200 text-teal-600 hover:bg-teal-50">
                  {chebiLoading ? <Loader2 className="w-3 h-3 animate-spin" /> : <BookOpen className="w-3 h-3 mr-1" />}From ChEBI
                </Button>
              </div>
            </div>

            {/* Ketcher Editor */}
            <ChemicalEditor
              ref={editorRef}
              height={470}
              className="border-violet-100"
              onReady={() => {
                setEditorReady(true);
                if (substrate.smiles) {
                  setTimeout(() => editorRef.current?.setMolecule(substrate.smiles), 500);
                }
              }}
            />

            {/* SMILES text input */}
            <div className="flex items-center gap-2">
              <LabelUI className="text-xs text-gray-400 uppercase tracking-wider shrink-0">SMILES:</LabelUI>
              <Input
                placeholder="Enter SMILES or draw in the editor above"
                value={smilesText}
                onChange={(e) => {
                  setSmilesText(e.target.value);
                  if (e.target.value) editorRef.current?.setMolecule(e.target.value);
                }}
                className="font-mono text-xs h-8 flex-1 border-gray-200"
              />
            </div>

            {/* Name input */}
            <div className="flex items-center gap-2 pt-2 border-t">
              <LabelUI className="text-xs text-gray-400 uppercase tracking-wider shrink-0">Name:</LabelUI>
              <Input
                placeholder="Enter molecule name"
                value={nameText}
                onChange={(e) => setNameText(e.target.value)}
                className="text-xs h-8 flex-1 border-gray-200"
              />
              <Button onClick={handleSaveName} disabled={!nameText.trim()}
                size="sm" className="h-8 text-xs bg-gray-600 hover:bg-gray-700 text-white">
                Save Name
              </Button>
            </div>
          </div>

          {/* RIGHT SECTION: PDB Ligand Mapping (2/5) */}
          <div className="lg:col-span-2 p-4 space-y-3 bg-gray-50/30">
            <div className="flex items-center gap-2 mb-2">
              <span className="text-xs font-semibold text-amber-700 uppercase tracking-wider">Cognate Ligand Mapping</span>
            </div>

            {/* Chain filter info bar */}
            <div className="flex items-center gap-2 p-2 rounded-md bg-blue-50/60 border border-blue-100">
              <Filter className="w-3 h-3 text-blue-500" />
              <span className="text-xs text-blue-700">
                {selectedChain
                  ? <>Showing <strong>{chainFilter === "selected" ? `${filteredLigands.length} ligand(s) on Chain ${selectedChain}` : `${filteredLigands.length} ligand(s) across ${allChains.length} chain(s)`}</strong>
                    {chainFilter === "selected" && selectedChainLigands.length === 0 && (
                      <span className="text-amber-600 ml-1">(no ligands on Chain {selectedChain} — showing all)</span>
                    )}
                  </>
                  : <strong>{filteredLigands.length} ligand(s) — no chain selected</strong>
                }
              </span>
              {selectedChain && (
                <Button
                  variant="outline"
                  size="sm"
                  className="h-5 text-xs ml-auto border-blue-200 text-blue-600 hover:bg-blue-100"
                  onClick={() => setChainFilter(prev => prev === "selected" ? "all" : "selected")}
                >
                  {chainFilter === "selected" ? "Show All Chains" : "Filter to Chain " + selectedChain}
                </Button>
              )}
            </div>

            {pdbLigands.length === 0 ? (
              <div className="rounded-md border border-dashed border-gray-200 p-4 text-center space-y-2">
                <p className="text-xs text-gray-400">No non-water ligands found in the loaded PDB structure.</p>
                {pdbInfo?.pdb_id && (
                  <p className="text-xs text-gray-300">Loaded: {pdbInfo.pdb_id} — {pdbInfo.ligands?.length === 0 ? "This structure contains only protein atoms and water molecules." : "PDB info may not have loaded correctly."}</p>
                )}
                {!pdbInfo?.pdb_id && (
                  <p className="text-xs text-amber-500">No PDB structure loaded. Go to Step 1 to load one.</p>
                )}
                {pdbInfo?.pdb_id && pdbInfo.ligands?.length === 0 && (
                  <p className="text-xs text-teal-500">Try PDBs with ligands: <span className="font-mono">1AKE</span> (kinase+ligand), <span className="font-mono">2VQE</span> (large complex), or <span className="font-mono">3NIR</span></p>
                )}
              </div>
            ) : (
              <>
                <div className="space-y-2">
                  <div className="flex items-center justify-between">
                    <LabelUI className="text-xs text-gray-500">Choose PDB Ligand</LabelUI>
                    {computingSimilarity && <span className="text-xs text-gray-400 flex items-center gap-1"><Loader2 className="w-3 h-3 animate-spin" />Computing similarity...</span>}
                  </div>
                  <select
                    value={selectedLigand}
                    onChange={(e) => setSelectedLigand(e.target.value)}
                    className="w-full text-xs border rounded-md px-2 py-1.5 bg-white border-gray-200 h-8"
                  >
                    <option value="">-- Select PDB Ligand (sorted by similarity → proximity) --</option>
                    {sortedLigands.map((l, i) => {
                      // Insert a visual separator before water molecules
                      const prevLig = i > 0 ? sortedLigands[i - 1] : null;
                      const isWaterStart = l.is_water && (!prevLig || !prevLig.is_water);
                      const hasRealDist = l.distance !== undefined && l.distance < 900;
                      return (
                        <option key={i} value={`${l.res_name}_${l.chain}_${l.res_num}`}>
                          {isWaterStart ? "── Water Molecules ── " : ""}{l.res_name} (Chain {l.chain}:{l.res_num}, {l.num_atoms} atoms{l.formula ? `, ${l.formula}` : ""})
                          {l.similarity > 0 && !l.is_water ? ` — sim: ${(l.similarity * 100).toFixed(0)}%` : ""}
                          {hasRealDist ? ` — ${l.distance.toFixed(1)}Å` : ""}
                          {l.chain !== selectedChain ? " ⚠ different chain" : ""}
                        </option>
                      );
                    })}
                  </select>
                  <Button onClick={handleSaveLigand} disabled={!selectedLigand}
                    size="sm" className="w-full h-7 text-xs bg-amber-500 hover:bg-amber-600 text-white">
                    Map to Selected Ligand
                  </Button>
                </div>

                {/* MCS Comparison + Overlay Section */}
                {selectedLigand && (
                  <div className="rounded-md border border-emerald-100 bg-emerald-50/30 p-3 space-y-2">
                    <div className="flex items-center gap-1.5">
                      <CheckCircle2 className="w-3.5 h-3.5 text-emerald-600" />
                      <span className="text-xs font-medium text-emerald-800">
                        Mapped: <span className="font-mono">{selectedLigand.replace(/_/g, " ")}</span>
                      </span>
                    </div>
                    {computingMcs && (
                      <div className="flex items-center justify-center py-4">
                        <Loader2 className="w-4 h-4 animate-spin text-emerald-500 mr-2" />
                        <span className="text-xs text-emerald-600">Computing MCS comparison...</span>
                      </div>
                    )}
                    {mcsStats && !computingMcs && (
                      <div className="flex items-center gap-3 text-xs">
                        <Badge variant="secondary" className="bg-emerald-100 text-emerald-700 text-xs">
                          Similarity: {(mcsStats.similarity * 100).toFixed(1)}%
                        </Badge>
                        <span className="text-gray-500">MCS: {mcsStats.mcs_atoms} atoms, {mcsStats.mcs_bonds} bonds</span>
                      </div>
                    )}
                    {mcsSvgA && mcsSvgB && (
                      <div className="grid grid-cols-2 gap-2 pt-1">
                        <div className="text-center space-y-1">
                          <span className="text-xs font-semibold text-violet-600 uppercase tracking-wider">Drawn Molecule</span>
                          <div className="rounded-md border border-gray-100 bg-white p-1 [&>svg]:w-full [&>svg]:h-auto"
                            dangerouslySetInnerHTML={{ __html: mcsSvgA }} />
                        </div>
                        <div className="text-center space-y-1">
                          <span className="text-xs font-semibold text-amber-600 uppercase tracking-wider">PDB Ligand</span>
                          <div className="rounded-md border border-gray-100 bg-white p-1 [&>svg]:w-full [&>svg]:h-auto"
                            dangerouslySetInnerHTML={{ __html: mcsSvgB }} />
                        </div>
                      </div>
                    )}
                    {mcsInfo && !computingMcs && (
                      <div className="text-xs text-gray-500 font-mono bg-white rounded p-1.5 border border-emerald-100">
                        {mcsInfo}
                      </div>
                    )}
                    <div className="text-xs text-gray-400 leading-relaxed">
                      Green highlighted atoms/bonds indicate the Maximum Common Substructure (MCS).
                      Verify the mapping by comparing residue counts and names with the 3D viewer on Step 2.
                    </div>

                    {/* === OVERLAY BUTTON === */}
                    <div className="pt-2 border-t border-emerald-200/60">
                      <Button
                        onClick={handleOverlay}
                        disabled={overlayLoading || !smilesText.trim() || !pdbRawText}
                        size="sm"
                        className={`w-full h-8 text-xs text-white shadow-sm ${
                          isSelectedLigandWater
                            ? "bg-gradient-to-r from-blue-500 to-cyan-500 hover:from-blue-600 hover:to-cyan-600"
                            : "bg-gradient-to-r from-violet-600 to-emerald-600 hover:from-violet-700 hover:to-emerald-700"
                        }`}
                      >
                        {overlayLoading ? (
                          <><Loader2 className="w-3 h-3 mr-1.5 animate-spin" />{isSelectedLigandWater ? "Mapping Water Coordinates..." : "Computing 3D Overlay..."}</>
                        ) : (
                          <><Layers className="w-3 h-3 mr-1.5" />{isSelectedLigandWater ? "Map Water Coordinates to PDB Position" : "Overlay to PDB Ligand 3D Position"}</>
                        )}
                      </Button>
                      <p className="text-xs text-gray-400 mt-1.5 leading-relaxed">
                        {isSelectedLigandWater
                          ? "Maps the water molecule's O atom directly to the PDB water coordinates. No MCS/Kabsch needed — direct coordinate placement."
                          : "Aligns the drawn molecule onto the PDB ligand's 3D coordinates using MCS + Kabsch algorithm. The overlay will be shown in the 3D viewer for visual verification."
                        }
                      </p>
                    </div>

                    {/* Overlay Result */}
                    {overlayResult && (
                      overlayResult.isWater ? (
                        /* Water overlay result — blue border, hide RMSD/MCS */
                        <div className="rounded-md border border-blue-200 bg-blue-50/50 p-2.5 space-y-1.5">
                          <div className="flex items-center gap-1.5">
                            <CheckCircle2 className="w-3.5 h-3.5 text-blue-600" />
                            <span className="text-xs font-medium text-blue-800">💧 Water Coordinate Mapping Applied</span>
                          </div>
                          <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs">
                            <div className="flex justify-between">
                              <span className="text-gray-500">Target:</span>
                              <span className="font-mono text-blue-700">{overlayResult.resName} {overlayResult.chain}:{overlayResult.resNum}</span>
                            </div>
                            <div className="flex justify-between">
                              <span className="text-gray-500">Mapped:</span>
                              <span className="font-mono">O atom → PDB coords</span>
                            </div>
                          </div>
                          {overlayResult.message && (
                            <p className="text-xs text-blue-600 leading-relaxed">{overlayResult.message}</p>
                          )}
                        </div>
                      ) : (
                        /* Non-water overlay result */
                        <div className={`rounded-md border p-2.5 space-y-1.5 ${
                          overlayResult.rmsd > 2.0
                            ? "border-amber-200 bg-amber-50/50"
                            : "border-emerald-200 bg-emerald-50/50"
                        }`}>
                          <div className="flex items-center gap-1.5">
                            {overlayResult.rmsd > 2.0 ? (
                              <AlertTriangle className="w-3.5 h-3.5 text-amber-600" />
                            ) : (
                              <CheckCircle2 className="w-3.5 h-3.5 text-emerald-600" />
                            )}
                            <span className="text-xs font-medium">
                              {overlayResult.rmsd > 2.0 ? "Overlay Applied (Low Quality)" : "Overlay Applied"}
                            </span>
                          </div>
                          <div className="grid grid-cols-2 gap-x-4 gap-y-0.5 text-xs">
                            <div className="flex justify-between">
                              <span className="text-gray-500">RMSD:</span>
                              <span className={`font-mono font-medium ${overlayResult.rmsd < 1.0 ? "text-emerald-700" : overlayResult.rmsd < 2.0 ? "text-amber-600" : "text-red-600"}`}>
                                {overlayResult.rmsd.toFixed(3)} Å
                              </span>
                            </div>
                            <div className="flex justify-between">
                              <span className="text-gray-500">Mapped Atoms:</span>
                              <span className="font-mono font-medium">{overlayResult.numMapped}</span>
                            </div>
                            <div className="flex justify-between">
                              <span className="text-gray-500">MCS Size:</span>
                              <span className="font-mono">{overlayResult.mcsNumAtoms} atoms</span>
                            </div>
                            <div className="flex justify-between">
                              <span className="text-gray-500">Target:</span>
                              <span className="font-mono">{overlayResult.resName} {overlayResult.chain}:{overlayResult.resNum}</span>
                            </div>
                          </div>
                          {overlayResult.rmsd > 2.0 && (
                            <p className="text-xs text-amber-600 leading-relaxed">
                              High RMSD indicates significant structural differences. The overlay is approximate —
                              verify in the 3D viewer and consider whether this is the correct PDB/ligand mapping.
                            </p>
                          )}
                        </div>
                      )
                    )}
                  </div>
                )}

                {/* Ligand list overview — sorted by similarity */}
                <div className="space-y-1">
                  <LabelUI className="text-xs text-gray-400 uppercase tracking-wider">
                    Available Ligands ({sortedLigands.length})
                    {sortedLigands.length > 0 && ` — sorted by similarity → proximity`}
                  </LabelUI>
                  <div className="max-h-48 overflow-y-auto rounded-md border border-gray-100">
                    <Table>
                      <TableBody>
                        {sortedLigands.slice(0, 30).map((l, i) => {
                          const ligKey = `${l.res_name}_${l.chain}_${l.res_num}`;
                          const isSelected = selectedLigand === ligKey;
                          const isOnChain = l.chain === selectedChain;
                          return (
                            <TableRow key={i}
                              className={`text-xs cursor-pointer hover:bg-amber-50 ${isSelected ? "bg-amber-50" : ""} ${!isOnChain ? "opacity-50" : ""} ${l.is_water ? "bg-blue-50/30" : ""}`}
                              onClick={() => setSelectedLigand(ligKey)}>
                              <TableCell className="font-mono font-medium py-1">
                                {l.is_water ? "💧" : ""}{l.res_name}
                              </TableCell>
                              <TableCell className="py-1">
                                <Badge variant="outline" className={`text-[7px] px-1 py-0 ${isOnChain ? "border-blue-200 text-blue-600" : "border-gray-200 text-gray-400"}`}>
                                  {l.chain}
                                </Badge>
                                :{l.res_num}
                              </TableCell>
                              <TableCell className="py-1 text-right">{l.num_atoms}</TableCell>
                              <TableCell className="py-1 text-right font-mono">
                                {l.is_water ? (
                                  <Badge variant="outline" className="text-[7px] px-1 py-0 border-blue-200 text-blue-500">💧 Water</Badge>
                                ) : (
                                  <span className={l.similarity > 0.5 ? "text-emerald-600 font-medium" : l.similarity > 0.2 ? "text-amber-600" : "text-gray-400"}>
                                    {l.similarity > 0 ? `${(l.similarity * 100).toFixed(0)}%` : "—"}
                                  </span>
                                )}
                              </TableCell>
                            </TableRow>
                          );
                        })}
                      </TableBody>
                    </Table>
                  </div>
                  {sortedLigands.length > 30 && (
                    <p className="text-xs text-gray-400 text-center">+{sortedLigands.length - 30} more ligands</p>
                  )}
                </div>

                {/* Help text when no SMILES available for many ligands */}
                {smilesText.trim() && sortedLigands.filter(l => l.similarity === 0).length > sortedLigands.length / 2 && (
                  <div className="rounded-md border border-amber-100 bg-amber-50/30 p-2">
                    <p className="text-xs text-amber-700 leading-relaxed">
                      <strong>Note:</strong> Many ligands show 0% similarity because their chemical structure is unknown
                      (not in the built-in dictionary and could not be extracted from the PDB file).
                      The ligands are sorted by proximity to the selected chain as a fallback.
                      Select the ligand manually by its name and chain position.
                    </p>
                  </div>
                )}
              </>
            )}
          </div>
        </div>
      </DialogContent>
    </Dialog>
  );
}
