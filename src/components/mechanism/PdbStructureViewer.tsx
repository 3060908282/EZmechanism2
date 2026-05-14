"use client";

import React, { useEffect, useRef, useState, useCallback, forwardRef, useImperativeHandle } from "react";
import { Loader2, Maximize2, Minimize2, RotateCcw, Eye, Crosshair, Box, Rotate3d, Atom, Sun, Moon, ChevronDown, Layers } from "lucide-react";
import { Button } from "@/components/ui/button";
import { Badge } from "@/components/ui/badge";
import type { SelectedResidue } from "@/lib/types";

// ---- Chain color palette ----
const CHAIN_COLORS: Record<string, string> = {
  A: "#4A90D9", B: "#E85D75", C: "#50C878", D: "#F5A623",
  E: "#9B59B6", F: "#1ABC9C", G: "#E67E22", H: "#3498DB",
  I: "#E74C3C", J: "#2ECC71", K: "#F39C12", L: "#8E44AD",
  M: "#16A085", N: "#D35400", O: "#2980B9", P: "#C0392B",
  Q: "#27AE60", R: "#F1C40F", S: "#7F8C8D", T: "#1D3557",
  U: "#457B9D", V: "#A8DADC", W: "#E63946", X: "#2A9D8F",
  Y: "#264653", Z: "#E9C46A",
};

function getChainColor(chainId: string): string {
  return CHAIN_COLORS[chainId] || CHAIN_COLORS[Object.keys(CHAIN_COLORS)[Object.keys(CHAIN_COLORS).length - 1]] || "#888888";
}

const RESIDUE_COLORS = ["#FF4444", "#FF8800", "#FFDD00", "#44FF44", "#4488FF", "#AA44FF", "#FF44AA", "#44FFFF"];

type DisplayMode = "cartoon" | "stick" | "ribbon";

interface CachedAtom { x: number; y: number; z: number; resn: string; resi: number; chain: string; }

export interface OverlayMolecule {
  substrateId: string;
  molblock: string;
  smiles: string;
  rmsd: number;
  resName: string;
  chain: string;
  resNum: number;
  color: string;
}

export interface PdbStructureViewerProps {
  pdbData: string;
  selectedResidues?: SelectedResidue[];
  activeChain?: string | null;
  chains?: string[];
  onChainChange?: (chain: string | null) => void;
  showActiveCenter?: boolean;
  onResidueClick?: (chain: string, resNum: number) => void;
  className?: string;
  height?: number;
  autoZoom?: "all" | "chain";
  overlayMolecules?: OverlayMolecule[];
  showOverlay?: boolean;
}

export interface PdbStructureViewerHandle {
  focusOnActiveCenter: () => void;
  focusOnResidue: (residue: SelectedResidue) => void;
  resetView: () => void;
  zoomToOverlay: (mol: OverlayMolecule) => void;
}

function ToggleButton({
  icon: Icon, label, active, onClick, showOn = false,
}: {
  icon: React.ElementType; label: string; active: boolean; onClick: () => void; showOn?: boolean;
}) {
  return (
    <Button
      size="sm"
      variant={active ? "default" : "ghost"}
      className={`h-6 px-2 text-xs gap-1 transition-all duration-150 ${
        active
          ? "bg-emerald-600/80 text-white hover:bg-emerald-500/80 shadow-sm shadow-emerald-500/30"
          : "text-gray-400 hover:text-white hover:bg-white/10"
      }`}
      onClick={onClick}
      title={label}
    >
      <Icon className="w-3 h-3" />
      <span className="hidden sm:inline">{label}</span>
      {showOn && active && <span className="text-xs font-bold ml-0.5 opacity-90">on</span>}
    </Button>
  );
}

const PdbStructureViewer = forwardRef<PdbStructureViewerHandle, PdbStructureViewerProps>(
  function PdbStructureViewerInner(
    { pdbData, selectedResidues = [], activeChain, chains = [], onChainChange, showActiveCenter = false, onResidueClick, className = "", height = 600, autoZoom, overlayMolecules = [], showOverlay = true },
    ref
  ) {
    const containerRef = useRef<HTMLDivElement>(null);
    const viewerRef = useRef<any>(null);
    const isInitialized = useRef(false);
    const modelRef = useRef<any>(null);
    const spinIntervalRef = useRef<ReturnType<typeof setInterval> | null>(null);
    const renderQueuedRef = useRef(false);

    // CA atom cache per chain
    const caAtomsCacheRef = useRef<Map<string, CachedAtom[]>>(new Map());

    // ---- Independent RES label management ----
    const resLabelRefsRef = useRef<any[]>([]);

    // Auto-zoom tracking
    const autoZoomDoneRef = useRef(false);
    const lastAutoZoomPdbRef = useRef<string>("");

    // Render state in refs
    const displayModeRef = useRef<DisplayMode>("cartoon");
    const resModeRef = useRef(false);
    const darkBgRef = useRef(true);
    const activeChainRef = useRef<string | null | undefined>(undefined);
    const selectedResiduesRef = useRef<SelectedResidue[]>([]);
    const onResidueClickRef = useRef<typeof onResidueClick>(undefined);
    const overlayModelRefsRef = useRef<Map<string, any>>(new Map());
    const [overlayVisible, setOverlayVisible] = useState(true);

    // React state
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);
    const [isFullscreen, setIsFullscreen] = useState(false);
    const [viewerReady, setViewerReady] = useState(false);
    const [displayMode, setDisplayModeState] = useState<DisplayMode>("cartoon");
    const [spinOn, setSpinOn] = useState(false);
    const [resMode, setResModeState] = useState(false);
    const [darkBg, setDarkBgState] = useState(true);
    const [viewerChains, setViewerChains] = useState<string[]>([]);
    const [renderTick, setRenderTick] = useState(0);
    const [activeCenterOpen, setActiveCenterOpen] = useState(false);

    const allChains = chains.length > 0 ? chains : viewerChains;

    // Sync React state → refs
    useEffect(() => { displayModeRef.current = displayMode; }, [displayMode]);
    useEffect(() => { resModeRef.current = resMode; }, [resMode]);
    useEffect(() => { darkBgRef.current = darkBg; }, [darkBg]);
    useEffect(() => { activeChainRef.current = activeChain; }, [activeChain]);
    useEffect(() => { selectedResiduesRef.current = selectedResidues; }, [selectedResidues]);
    useEffect(() => { onResidueClickRef.current = onResidueClick; }, [onResidueClick]);

    // Close dropdown on outside click
    useEffect(() => {
      if (!activeCenterOpen) return;
      const handler = () => setActiveCenterOpen(false);
      document.addEventListener("click", handler);
      return () => document.removeEventListener("click", handler);
    }, [activeCenterOpen]);

    // ---- Queue a single re-render (debounced via microtask) ----
    const queueRender = useCallback(() => {
      if (renderQueuedRef.current) return;
      renderQueuedRef.current = true;
      queueMicrotask(() => {
        renderQueuedRef.current = false;
        setRenderTick(t => t + 1);
      });
    }, []);

    // ---- Remove all stored RES labels ----
    const clearResLabels = useCallback(() => {
      const viewer = viewerRef.current;
      if (!viewer) return;
      for (const label of resLabelRefsRef.current) {
        try { viewer.removeLabel(label); } catch { /* */ }
      }
      resLabelRefsRef.current = [];
    }, []);

    // ---- Add RES labels for selected residues (color-coded to match atom highlight colors) ----
    const addResLabels = useCallback(() => {
      const viewer = viewerRef.current;
      if (!viewer) return;
      for (const label of resLabelRefsRef.current) {
        try { viewer.removeLabel(label); } catch { /* */ }
      }
      resLabelRefsRef.current = [];

      const isDark = darkBgRef.current;
      const residues = selectedResiduesRef.current;

      for (let i = 0; i < residues.length; i++) {
        const r = residues[i];
        const cache = caAtomsCacheRef.current.get(r.chain);
        const caAtom = cache?.find(a => a.resi === r.res_num && a.resn === r.res_name);
        if (!caAtom) continue;

        const color = RESIDUE_COLORS[i % RESIDUE_COLORS.length];
        const label = viewer.addLabel(`${r.res_name} ${r.res_num} (${r.chain})`, {
          position: { x: caAtom.x, y: caAtom.y + 2.0, z: caAtom.z },
          fontSize: 11,
          fontColor: "#ffffff",
          fontWeight: "bold",
          backgroundColor: color,
          backgroundOpacity: 0.92,
          borderColor: "#ffffff",
          borderThickness: 1.5,
          inFront: true,
          showBackground: true,
          borderRadius: 4,
        });
        if (label) resLabelRefsRef.current.push(label);
      }

      viewer.render();
    }, []);

    // ---- State setters ----
    const setDisplayMode = useCallback((v: DisplayMode) => {
      displayModeRef.current = v;
      setDisplayModeState(v);
      queueRender();
    }, [queueRender]);

    const toggleResMode = useCallback(() => {
      const next = !resModeRef.current;
      resModeRef.current = next;
      setResModeState(next);
    }, []);

    const toggleDarkBg = useCallback(() => {
      const next = !darkBgRef.current;
      darkBgRef.current = next;
      setDarkBgState(next);
      queueRender();
    }, [queueRender]);

    // ---- Chain cycling ----
    const handleChainClick = useCallback(() => {
      if (allChains.length <= 1) return;
      const current = activeChain;
      if (!current || current === "") {
        onChainChange?.(allChains[0]);
      } else {
        const idx = allChains.indexOf(current);
        if (idx === allChains.length - 1) {
          onChainChange?.(null);
        } else {
          onChainChange?.(allChains[idx + 1]);
        }
      }
    }, [allChains, onChainChange, activeChain]);

    // ---- RES effect ----
    useEffect(() => {
      if (!viewerReady || !pdbData || loading) return;
      if (resMode) {
        addResLabels();
      } else {
        clearResLabels();
        viewerRef.current?.render();
      }
    }, [resMode, selectedResidues, viewerReady, pdbData, loading, addResLabels, clearResLabels]);

    // ---- Initialize 3Dmol viewer ----
    useEffect(() => {
      if (!containerRef.current || isInitialized.current) return;
      const initViewer = () => {
        try {
          const $3Dmol = (window as any).$3Dmol;
          if (!$3Dmol) {
            setError("3Dmol library not loaded — check your internet connection and refresh");
            return;
          }
          if (!containerRef.current) return;
          const viewer = $3Dmol.createViewer(containerRef.current, {
            backgroundColor: "#0f172a", backgroundAlpha: 1.0, antialias: true,
            disableFog: false, cartoonQuality: 5, lowerZoomLimit: 5, upperZoomLimit: 500,
          });
          viewerRef.current = viewer;
          isInitialized.current = true;
          setViewerReady(true);
        } catch (err) {
          console.error("Failed to init 3Dmol:", err);
          setError("3D viewer initialization failed");
        }
      };

      // Wait for the global $3Dmol to be available (script tag may load after component mount)
      if ((window as any).$3Dmol) {
        initViewer();
      } else {
        const checkInterval = setInterval(() => {
          if ((window as any).$3Dmol) {
            clearInterval(checkInterval);
            initViewer();
          }
        }, 200);
        // Timeout after 10 seconds
        const timeout = setTimeout(() => {
          clearInterval(checkInterval);
          if (!isInitialized.current) {
            setError("3Dmol library loading timed out — please refresh the page");
          }
        }, 10000);
        return () => {
          clearInterval(checkInterval);
          clearTimeout(timeout);
        };
      }

      return () => {
        clearResLabels();
        if (spinIntervalRef.current) clearInterval(spinIntervalRef.current);
        if (viewerRef.current) {
          try { viewerRef.current.removeAllModels(); viewerRef.current.removeAllSurfaces(); viewerRef.current.removeAllShapes(); viewerRef.current.removeAllLabels(); } catch { /* */ }
        }
      };
    }, []);

    // ---- Load PDB data ----
    useEffect(() => {
      if (!viewerReady || !viewerRef.current || !pdbData) return;
      const viewer = viewerRef.current;
      React.startTransition(() => { setLoading(true); setError(null); });
      resLabelRefsRef.current = [];
      autoZoomDoneRef.current = false;

      try {
        viewer.removeAllModels();
        viewer.removeAllSurfaces();
        viewer.removeAllShapes();
        viewer.removeAllLabels();
        if (spinIntervalRef.current) { clearInterval(spinIntervalRef.current); spinIntervalRef.current = null; }

        // Validate PDB data has ATOM/HETATM records
        const dataLines = pdbData.split('\n');
        const atomLineCount = dataLines.filter(l => l.startsWith('ATOM') || l.startsWith('HETATM')).length;
        if (atomLineCount === 0) {
          React.startTransition(() => {
            setError(`PDB data has no ATOM records (${pdbData.length} chars, ${dataLines.length} lines)`);
            setLoading(false);
          });
          return;
        }

        modelRef.current = viewer.addModel(pdbData, "pdb");
        const chainIds: string[] = viewer.getUniqueValues("chain");
        React.startTransition(() => setViewerChains(chainIds));

        const caCache = new Map<string, CachedAtom[]>();
        chainIds.forEach((chainId: string) => {
          const caAtoms = viewer.selectedAtoms({ chain: chainId, hetflag: false, atom: "CA" });
          const mapped = caAtoms.map((a: any) => ({
            x: a.x, y: a.y, z: a.z, resn: a.resn, resi: a.resi, chain: a.chain,
          }));
          caCache.set(chainId, mapped);
        });
        caAtomsCacheRef.current = caCache;

        console.log(`[PdbViewer] Loaded: ${atomLineCount} atoms, ${chainIds.length} chains: [${chainIds.join(',')}]`);
        React.startTransition(() => setLoading(false));
      } catch (err) {
        console.error("Failed to load PDB data:", err);
        React.startTransition(() => {
          setError(`Failed to parse PDB structure: ${err instanceof Error ? err.message : String(err)}`);
          setLoading(false);
        });
      }
    }, [pdbData, viewerReady, clearResLabels]);

    // ---- Main render effect ----
    useEffect(() => {
      if (!viewerReady || !viewerRef.current || !pdbData || loading) return;
      const viewer = viewerRef.current;
      resLabelRefsRef.current = [];

      try {
        viewer.removeAllShapes();
        viewer.removeAllLabels();

        const mode = displayModeRef.current;
        const chain = activeChainRef.current;
        const isDark = darkBgRef.current;
        const residues = selectedResiduesRef.current;
        const onResClick = onResidueClickRef.current;
        const hasActiveFilter = chain != null && chain !== "";

        const chainIds: string[] = viewer.getUniqueValues("chain");
        if (mode === "cartoon") {
          chainIds.forEach((chainId: string) => {
            viewer.setStyle({ chain: chainId, hetflag: false }, {
              cartoon: {
                color: getChainColor(chainId),
                opacity: !hasActiveFilter ? 1.0 : chainId === chain ? 1.0 : 0.35,
                arrows: true, tubes: true, thickness: 0.4,
              },
            });
          });
        } else if (mode === "ribbon") {
          chainIds.forEach((chainId: string) => {
            viewer.setStyle({ chain: chainId, hetflag: false }, {
              cartoon: {
                style: "rectangle",
                color: getChainColor(chainId),
                opacity: !hasActiveFilter ? 1.0 : chainId === chain ? 1.0 : 0.35,
                arrows: true, thickness: 0.2,
              },
            });
          });
        } else if (mode === "stick") {
          chainIds.forEach((chainId: string) => {
            const isSelected = !hasActiveFilter || chainId === chain;
            viewer.setStyle({ chain: chainId, hetflag: false }, {
              stick: { colorscheme: "Jmol", radius: isSelected ? 0.12 : 0.08, opacity: isSelected ? 1.0 : 0.2 },
              sphere: { colorscheme: "Jmol", scale: isSelected ? 0.25 : 0.12, opacity: isSelected ? 1.0 : 0.2 },
            });
          });
        }
        viewer.setStyle({ hetflag: true, resn: ["HOH", "WAT"] }, { cross: { linewidth: 1, hidden: true } });
        viewer.setStyle({ hetflag: true }, { stick: { radius: 0.15, colorscheme: "Jmol", opacity: 0.9 } });

        // Highlight selected residues
        if (residues.length > 0) {
          residues.forEach((residue, idx) => {
            const color = RESIDUE_COLORS[idx % RESIDUE_COLORS.length];
            viewer.addStyle(
              { chain: residue.chain, resi: residue.res_num, resn: residue.res_name },
              {
                stick: { radius: mode === "stick" ? 0.22 : 0.25, color, opacity: 1.0 },
                sphere: { scale: mode === "stick" ? 0.4 : 0.45, color, opacity: 1.0 },
              }
            );
          });
        }

        // Interactions
        viewer.setClickable({ hetflag: false }, true, (atom: any) => {
          if (atom?.chain && atom.resi !== undefined) onResClick?.(atom.chain, atom.resi);
        });
        viewer.setHoverable(
          { hetflag: false }, true,
          (atom: any) => {
            if (atom?.label) {
              const lbl = viewer.addLabel(`${atom.resn}${atom.resi} (${atom.chain})`, {
                position: atom, fontSize: 13, fontWeight: "bold",
                fontColor: isDark ? "#ffffff" : "#1a1a1a",
                backgroundColor: isDark ? "rgba(0,0,0,0.8)" : "rgba(255,255,255,0.9)",
                backgroundOpacity: 0.9,
                borderColor: isDark ? "#60a5fa" : "#2563eb", borderThickness: 1.2, inFront: true,
              });
              atom.hoverLabel = lbl;
            }
          },
          (atom: any) => {
            if (atom?.hoverLabel) { viewer.removeLabel(atom.hoverLabel); atom.hoverLabel = null; }
          }
        );

        viewer.render();

        if (resModeRef.current) {
          addResLabels();
        }
      } catch (err) {
        console.error("Render error:", err);
      }
    }, [renderTick, activeChain, viewerReady, pdbData, loading, addResLabels]);

    // ---- Render overlay molecules from alignment pipeline ----
    useEffect(() => {
      const viewer = viewerRef.current;
      if (!viewer || !viewerReady || !pdbData || loading) return;

      // Remove all previous overlay models
      overlayModelRefsRef.current.forEach((model) => {
        try { if (model) viewer.removeModel(model); } catch { /* */ }
      });
      overlayModelRefsRef.current.clear();

      if (!showOverlay || !overlayMolecules || overlayMolecules.length === 0) {
        viewer.render();
        return;
      }

      // Add each overlay molecule as a separate model
      for (const mol of overlayMolecules) {
        if (!mol.molblock) continue;
        try {
          const overlayModel = viewer.addModel(mol.molblock, "sdf");
          // Style: colored sticks with spheres
          viewer.setStyle({ model: overlayModel }, {
            stick: { radius: 0.18, color: mol.color, opacity: 0.95 },
            sphere: { scale: 0.35, color: mol.color, opacity: 0.95 },
          });
          overlayModelRefsRef.current.set(mol.substrateId, overlayModel);
        } catch (err) {
          console.error("Failed to add overlay model:", err);
        }
      }

      viewer.render();
    }, [overlayMolecules, showOverlay, viewerReady, pdbData, loading]);

    // ---- Auto-zoom after PDB load ----
    useEffect(() => {
      if (!viewerReady || !viewerRef.current || !pdbData || loading) return;
      if (autoZoomDoneRef.current && lastAutoZoomPdbRef.current === pdbData) return;

      const timer = setTimeout(() => {
        const viewer = viewerRef.current;
        if (!viewer) return;
        autoZoomDoneRef.current = true;
        lastAutoZoomPdbRef.current = pdbData;

        if (autoZoom === "chain" && activeChainRef.current) {
          viewer.zoomTo({ chain: activeChainRef.current, hetflag: false }, 800);
        } else if (autoZoom === "all" || !autoZoom) {
          viewer.zoomTo(undefined, 800);
        }
        viewer.render();
      }, 300);

      return () => clearTimeout(timer);
    }, [pdbData, viewerReady, loading, autoZoom]);

    // ---- Background color change ----
    useEffect(() => {
      if (!viewerRef.current || !pdbData) return;
      try {
        viewerRef.current.setBackgroundColor(darkBgRef.current ? "#0f172a" : "#ffffff");
        viewerRef.current.render();
      } catch { /* */ }
    }, [darkBg, pdbData]);

    // ---- Spin toggle ----
    useEffect(() => {
      if (!viewerReady || !viewerRef.current || !pdbData) return;
      if (spinOn) {
        spinIntervalRef.current = setInterval(() => {
          if (viewerRef.current) { viewerRef.current.rotate(0.5, "y"); viewerRef.current.render(); }
        }, 50);
      } else {
        if (spinIntervalRef.current) { clearInterval(spinIntervalRef.current); spinIntervalRef.current = null; }
      }
      return () => {
        if (spinIntervalRef.current) { clearInterval(spinIntervalRef.current); spinIntervalRef.current = null; }
      };
    }, [spinOn, viewerReady, pdbData]);

    // ---- Handle resize ----
    useEffect(() => {
      const handleResize = () => { if (viewerRef.current) try { viewerRef.current.resize(); } catch { /* */ } };
      window.addEventListener("resize", handleResize);
      const timer = setTimeout(handleResize, 100);
      return () => { window.removeEventListener("resize", handleResize); clearTimeout(timer); };
    }, [isFullscreen]);

    // ---- Focus on a single specific residue ----
    const focusOnResidue = useCallback((residue: SelectedResidue) => {
      const viewer = viewerRef.current;
      if (!viewer) return;
      viewer.zoomTo(
        { chain: residue.chain, resi: residue.res_num, resn: residue.res_name, hetflag: false },
        800
      );
      viewer.render();
      setActiveCenterOpen(false);
    }, []);

    const focusOnActiveCenter = useCallback(() => {
      // Just toggles the dropdown
      setActiveCenterOpen(v => !v);
    }, []);

    const resetView = useCallback(() => {
      if (!viewerRef.current) return;
      setSpinOn(false);
      if (spinIntervalRef.current) { clearInterval(spinIntervalRef.current); spinIntervalRef.current = null; }
      displayModeRef.current = "cartoon";
      setDisplayModeState("cartoon");
      resModeRef.current = false;
      setResModeState(false);
      darkBgRef.current = true;
      setDarkBgState(true);
      try {
        viewerRef.current.setBackgroundColor("#0f172a");
        viewerRef.current.zoomTo(undefined, 800);
        viewerRef.current.render();
      } catch { /* */ }
      queueRender();
    }, [queueRender]);

    const zoomToOverlay = useCallback((mol: OverlayMolecule) => {
      const viewer = viewerRef.current;
      const overlayModel = overlayModelRefsRef.current.get(mol.substrateId);
      if (viewer && overlayModel) {
        viewer.zoomTo({ model: overlayModel }, 800);
        viewer.render();
      }
    }, []);

    const toggleOverlay = useCallback(() => {
      const next = !overlayVisible;
      setOverlayVisible(next);
      const viewer = viewerRef.current;
      if (!viewer) return;
      // Show/hide overlay models by setting style
      overlayModelRefsRef.current.forEach((model) => {
        if (model) {
          viewer.setStyle({ model }, {
            stick: { radius: 0.15, opacity: next ? 1.0 : 0.0 },
            sphere: { scale: 0.3, opacity: next ? 1.0 : 0.0 },
          });
        }
      });
      viewer.render();
    }, [overlayVisible]);

    useImperativeHandle(ref, () => ({ focusOnActiveCenter, focusOnResidue, resetView, zoomToOverlay }), [focusOnActiveCenter, focusOnResidue, resetView, zoomToOverlay]);

    // ---- No PDB data ----
    if (!pdbData) {
      return (
        <div>
          <div className="flex items-center gap-2 mb-2">
            <Badge variant="outline" className="text-xs h-5 border-gray-300 text-gray-500 bg-gray-50">3D Structure</Badge>
          </div>
          <div className={`rounded-xl border border-gray-200 bg-gray-50/80 flex flex-col items-center justify-center ${className}`} style={{ minHeight: height }}>
            <Box className="w-10 h-10 text-gray-300 mb-2" />
            <p className="text-xs text-gray-400">Load a PDB structure to view 3D preview</p>
            <p className="text-xs text-gray-300 mt-1">Step 1: Fetch or paste a PDB file</p>
          </div>
        </div>
      );
    }

    const bgColor = darkBg ? "#0f172a" : "#ffffff";
    const toolbarBgStart = darkBg ? "from-[#0f172a]/95" : "from-white/95";
    const toolbarBgEnd = darkBg ? "to-[#0f172a]/60" : "to-white/60";
    const legendBgStart = darkBg ? "from-[#0f172a]/90" : "from-white/90";
    const legendBgEnd = darkBg ? "to-[#0f172a]/40" : "to-white/40";
    const bottomBarBgStart = darkBg ? "from-[#0f172a]/80" : "from-white/80";
    const pillBg = darkBg ? "bg-black/30" : "bg-black/10";
    const resBadgeBorder = darkBg ? "border-amber-400/40 text-amber-300 bg-amber-900/30" : "border-amber-400 text-amber-700 bg-amber-50";
    const toggleIdleColor = darkBg ? "text-gray-400 hover:text-white hover:bg-white/10" : "text-gray-500 hover:text-gray-800 hover:bg-gray-200/60";
    const chainLabel = !activeChain ? "All" : `Chain ${activeChain}`;
    const fullscreenHeight = isFullscreen ? `calc(${height}px + 120px)` : height;

    return (
      <div className={className}>
        {/* Viewer container — no duplicate buttons above */}
        <div className="relative rounded-xl border border-gray-200 overflow-hidden" style={{ backgroundColor: bgColor }}>
          {/* Top Toolbar */}
          <div className={`absolute top-0 left-0 right-0 z-10 flex items-center justify-between px-3 py-1.5 bg-gradient-to-b ${toolbarBgStart} ${toolbarBgEnd}`}>
            <div className="flex items-center gap-1.5">
              {loading && <Loader2 className="w-3 h-3 text-gray-400 animate-spin" />}
              {!loading && !error && (
                <Badge variant="secondary" className="bg-emerald-900/50 text-emerald-400 text-xs h-5 border-emerald-700/30">3D Structure</Badge>
              )}
              {!loading && !error && allChains.length > 1 && (
                <button
                  onClick={(e) => { e.stopPropagation(); handleChainClick(); }}
                  className={`text-xs px-1.5 py-0.5 rounded-full border font-medium transition-all duration-150 cursor-pointer hover:opacity-80 ${
                    !activeChain ? (darkBg ? "border-gray-400/50 text-gray-200 bg-gray-600/30" : "border-gray-400 text-gray-600 bg-gray-100") : ""
                  }`}
                  style={activeChain ? {
                    borderColor: getChainColor(activeChain),
                    color: darkBg ? "#fff" : getChainColor(activeChain),
                    backgroundColor: darkBg ? getChainColor(activeChain) + "30" : getChainColor(activeChain) + "15",
                  } : undefined}
                  title="Click to switch chain"
                >
                  {chainLabel}
                </button>
              )}
              {selectedResidues.length > 0 && (
                <Badge variant="outline" className={`text-xs h-5 ${resBadgeBorder}`}>{selectedResidues.length} residue(s)</Badge>
              )}
            </div>

            <div className="flex items-center gap-1">
              {/* Show Active Center — dropdown list of selected residues */}
              {showActiveCenter && selectedResidues.length > 0 && (
                <div className="relative" onClick={(e) => e.stopPropagation()}>
                  <Button
                    size="sm"
                    variant="ghost"
                    className={`h-6 px-2 text-xs gap-1 text-amber-300 hover:text-amber-200 hover:bg-amber-900/40 ${activeCenterOpen ? "bg-amber-900/40" : ""}`}
                    onClick={() => setActiveCenterOpen(v => !v)}
                  >
                    <Crosshair className="w-3 h-3" />
                    <span className="hidden sm:inline">Active Center</span>
                    <ChevronDown className={`w-2.5 h-2.5 transition-transform ${activeCenterOpen ? "rotate-180" : ""}`} />
                  </Button>
                  {activeCenterOpen && (
                    <div className="absolute right-0 top-full mt-1 z-50 min-w-[180px] max-h-60 overflow-y-auto rounded-lg border border-gray-200 bg-white shadow-lg">
                      <div className="px-2 py-1.5 text-xs font-medium text-gray-500 border-b border-gray-100 bg-gray-50 rounded-t-lg sticky top-0">
                        Select residue to focus
                      </div>
                      {selectedResidues.map((r, idx) => {
                        const color = RESIDUE_COLORS[idx % RESIDUE_COLORS.length];
                        return (
                          <button
                            key={`${r.chain}-${r.res_num}-${idx}`}
                            className="w-full flex items-center gap-2 px-3 py-1.5 text-xs hover:bg-gray-50 transition-colors text-left border-b border-gray-50 last:border-0 cursor-pointer"
                            onClick={() => focusOnResidue(r)}
                          >
                            <div className="w-2 h-2 rounded-full shrink-0" style={{ backgroundColor: color }} />
                            <span className="font-mono font-medium text-gray-800">{r.res_name}</span>
                            <span className="font-mono text-gray-500">{r.res_num}</span>
                            <Badge variant="secondary" className="text-xs bg-gray-100 ml-auto shrink-0">Chain {r.chain}</Badge>
                          </button>
                        );
                      })}
                    </div>
                  )}
                </div>
              )}

              <Button
                size="sm" variant="outline"
                className={`h-6 px-2 text-xs gap-1 rounded-md font-medium transition-all duration-150 ${
                  darkBg ? "border-amber-400/60 text-amber-300 bg-amber-900/30 hover:bg-amber-800/40 hover:text-amber-200"
                    : "border-gray-400/60 text-gray-600 bg-gray-100 hover:bg-gray-200 hover:text-gray-700"
                }`}
                onClick={toggleDarkBg}
                title={darkBg ? "Switch to white background" : "Switch to dark background"}
              >
                {darkBg ? <Sun className="w-3 h-3" /> : <Moon className="w-3 h-3" />}
              </Button>
              <Button
                size="sm" variant="ghost"
                className={`h-6 px-2 text-xs ${toggleIdleColor}`}
                onClick={() => {
                  setIsFullscreen(v => !v);
                  setTimeout(() => { if (viewerRef.current) { viewerRef.current.resize(); viewerRef.current.render(); } }, 50);
                }}
                title={isFullscreen ? "Exit fullscreen" : "Fullscreen"}
              >
                {isFullscreen ? <Minimize2 className="w-3 h-3" /> : <Maximize2 className="w-3 h-3" />}
              </Button>
            </div>
          </div>

          {/* Bottom Toolbar */}
          <div className={`absolute bottom-8 left-0 right-0 z-10 flex flex-wrap items-center justify-center gap-1.5 px-3 py-1.5 bg-gradient-to-t ${bottomBarBgStart} to-transparent`}>
            <div className={`flex items-center gap-0.5 ${pillBg} rounded-lg p-0.5`}>
              <ToggleButton icon={Atom} label="Cartoon" active={displayMode === "cartoon"} onClick={() => setDisplayMode("cartoon")} />
              <ToggleButton icon={Rotate3d} label="Ribbon" active={displayMode === "ribbon"} onClick={() => setDisplayMode("ribbon")} />
              <ToggleButton icon={Box} label="Stick" active={displayMode === "stick"} onClick={() => setDisplayMode("stick")} />
            </div>
            <div className="w-px h-4 bg-white/10" />
            <div className={`flex items-center gap-0.5 ${pillBg} rounded-lg p-0.5`}>
              <ToggleButton icon={RotateCcw} label="Spin" active={spinOn} onClick={() => setSpinOn(v => !v)} showOn />
              <ToggleButton icon={Eye} label="Res" active={resMode} onClick={toggleResMode} showOn />
              <ToggleButton icon={Maximize2} label="Reset" active={false} onClick={resetView} />
            </div>
            {/* Overlay toggle — only show when overlay molecules exist */}
            {overlayMolecules.length > 0 && (
              <>
                <div className="w-px h-4 bg-white/10" />
                <ToggleButton icon={Layers} label="Overlay" active={overlayVisible} onClick={toggleOverlay} showOn />
              </>
            )}
          </div>

          {/* 3Dmol container */}
          <div ref={containerRef} className="w-full" style={{ height: fullscreenHeight, minHeight: 300 }} />

          {error && (
            <div className="absolute inset-0 flex flex-col items-center justify-center" style={{ backgroundColor: bgColor }}>
              <p className="text-xs text-red-400">{error}</p>
            </div>
          )}

          {/* Chain legend */}
          {!loading && !error && (
            <div className={`absolute bottom-0 left-0 right-0 z-10 flex flex-wrap items-center gap-1.5 px-3 py-1 bg-gradient-to-t ${legendBgStart} ${legendBgEnd}`}>
              {allChains.map((chainId: string) => (
                <div key={chainId} className="flex items-center gap-1">
                  <div className="w-2 h-2 rounded-full" style={{ backgroundColor: getChainColor(chainId) }} />
                  <span className={`text-xs font-mono ${darkBg ? "text-gray-400" : "text-gray-600"}`}>Chain {chainId}</span>
                </div>
              ))}
              <div className="ml-auto">
                {selectedResidues.length > 0 && <span className={`text-xs ${darkBg ? "text-gray-500" : "text-gray-400"}`}>Click atom to identify residue</span>}
                {overlayMolecules.length > 0 && (
                  <span className={`text-xs ${darkBg ? "text-violet-400" : "text-violet-600"} ml-2`}>
                    {overlayVisible ? "Overlay" : "Overlay hidden"}: {overlayMolecules.length} molecule(s)
                    {overlayMolecules.some(m => m.rmsd > 2) && <span className="text-amber-500 ml-1">⚠ high RMSD</span>}
                  </span>
                )}
              </div>
            </div>
          )}
        </div>
      </div>
    );
  }
);

export default PdbStructureViewer;
