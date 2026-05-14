"use client";

import React, { useState, useEffect, useRef, useCallback, useMemo } from "react";
import cytoscape, { type Core } from "cytoscape";
import { Network, GitBranch, Zap, Clock, Layers, ArrowRight, ArrowRightLeft, RotateCcw, Maximize2, Minimize2, Ruler } from "lucide-react";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { Badge } from "@/components/ui/badge";
import { Button } from "@/components/ui/button";
import { AlertCircle } from "lucide-react";
import { MechanismSearchResult } from "@/lib/types";
import MoleculeSVG from "@/components/mechanism/MoleculeSVG";

// ── Color palette per EzMechanism paper specification ──
const COLORS = {
  reactant:      { bg: "#f97316", border: "#ea580c", text: "#ffffff" },
  product:       { bg: "#ef4444", border: "#dc2626", text: "#ffffff" },
  exploration:   { bg: "#eab308", border: "#ca8a04", text: "#ffffff" },
  intermediate:  { bg: "#9ca3af", border: "#6b7280", text: "#ffffff" },
  edgeNormal:    "#9ca3af",
  edgeDbOnly:    "#ef4444",
  edgeCross:     "#f97316",
  highlighted:   "#f59e0b",
} as const;

/** Create a deterministic key from a list of SMILES */
function molKey(mols: string[] | undefined) {
  if (!mols || mols.length === 0) return "";
  return [...mols].sort().join("|");
}

export default function MechanismGraphPanel({ result }: { result: MechanismSearchResult | null }) {
  const cyContainerRef = useRef<HTMLDivElement>(null);
  const cyRef = useRef<Core | null>(null);
  const [selectedPath, setSelectedPath] = useState<number>(0);
  const [selectedNode, setSelectedNode] = useState<{ id: string; molecules: string[]; source: string; label: string } | null>(null);
  const [isExpanded, setIsExpanded] = useState(false);

  // ── Precompute the layout data from paths ──
  // Maps molecule-key → { id, row, col, label, source, molecules, colorKey }
  const layoutData = useMemo(() => {
    if (!result?.paths?.length || !result?.graph) return null;

    const graphNodes = result.graph.nodes;
    const graphEdges = result.graph.edges;

    // 1. Gather all unique config keys from ALL paths, and assign rows by step position
    const configKeyInfo = new Map<string, {
      nodeIds: string[];
      molecules: string[];
      source: string;
      label: string;
      smilesLabel: string;
      depth: number;
    }>();

    // Helper: find a graph node matching a molecule key
    const findGraphNodes = (mols: string[] | undefined) => {
      const key = molKey(mols);
      if (!key) return [];
      const matches = graphNodes.filter((n) => {
        const nk = molKey(n.molecules);
        return nk === key;
      });
      return matches.map((n) => n.id);
    };

    // Determine the maximum number of steps across all paths
    const maxSteps = Math.max(...result.paths.map((p) => p.total_steps ?? p.steps?.length ?? 0));

    // For each step row (0 = reactant side, 1..maxSteps-1 = intermediates, maxSteps = product side)
    const rowConfigs = new Map<number, Set<string>>();
    // Always add row 0 and row maxSteps
    rowConfigs.set(0, new Set());
    rowConfigs.set(maxSteps, new Set());
    for (let r = 1; r < maxSteps; r++) {
      rowConfigs.set(r, new Set());
    }

    // Reactant row (row 0): collect from_config of step 1 of every path
    const reactantKeys = new Set<string>();
    const productKeys = new Set<string>();

    result.paths.forEach((path) => {
      const steps = path.steps || [];
      if (steps.length === 0) return;

      // Reactant is the from_config of step 1
      const rKey = molKey(steps[0].from_config?.molecules);
      if (rKey) reactantKeys.add(rKey);

      // Product is the to_config of the last step
      const lastStep = steps[steps.length - 1];
      const pKey = molKey(lastStep.to_config?.molecules);
      if (pKey) productKeys.add(pKey);

      // Intermediate rows: to_config of step N goes to row N
      steps.forEach((step, idx) => {
        // from_config of step goes to row idx
        const fKey = molKey(step.from_config?.molecules);
        if (fKey && idx > 0) rowConfigs.get(idx)?.add(fKey);
        // to_config of step goes to row idx + 1
        const tKey = molKey(step.to_config?.molecules);
        if (tKey && idx + 1 < maxSteps) rowConfigs.get(idx + 1)?.add(tKey);
      });
    });

    // Ensure reactant and product are in their rows
    reactantKeys.forEach((k) => rowConfigs.get(0)?.add(k));
    productKeys.forEach((k) => rowConfigs.get(maxSteps)?.add(k));

    // 2. Register all configs with their info
    const allKeys = new Set<string>();
    rowConfigs.forEach((keys) => keys.forEach((k) => allKeys.add(k)));

    allKeys.forEach((key) => {
      if (configKeyInfo.has(key)) return;
      const mols = key.split("|").filter(Boolean);
      const nodeIds = findGraphNodes(mols);
      const firstNode = graphNodes.find((n) => molKey(n.molecules) === key);
      configKeyInfo.set(key, {
        nodeIds,
        molecules: mols,
        source: firstNode?.source || "exploration",
        label: firstNode?.label || "",
        smilesLabel: firstNode?.smiles_label || firstNode?.label || "",
        depth: firstNode?.depth ?? -1,
      });
    });

    // 3. Assign row and column positions
    // Find which row each key belongs to
    const keyToRow = new Map<string, number>();
    rowConfigs.forEach((keys, row) => {
      keys.forEach((k) => keyToRow.set(k, row));
    });

    // Assign columns within each row
    const keyToCol = new Map<string, number>();
    const rowEntries = Array.from(rowConfigs.entries()).sort(([a], [b]) => a - b);
    rowEntries.forEach(([, keys]) => {
      const sorted = Array.from(keys).sort();
      sorted.forEach((k, col) => keyToCol.set(k, col));
    });

    // 4. Compute node positions (x, y)
    // Layout: top-to-bottom, rows equally spaced
    const nodePositions = new Map<string, { x: number; y: number }>();
    const totalRows = maxSteps + 1; // row 0 to row maxSteps
    const rowSpacing = 120;
    const colSpacing = 100;

    // Find max columns per row to center
    const rowColCount = new Map<number, number>();
    rowEntries.forEach(([, keys]) => {
      rowColCount.set(0, 0);
    });

    keyToRow.forEach((row, key) => {
      const col = keyToCol.get(key) ?? 0;
      const currentMax = rowColCount.get(row) ?? 0;
      rowColCount.set(row, Math.max(currentMax, col + 1));
    });

    allKeys.forEach((key) => {
      const row = keyToRow.get(key) ?? 0;
      const col = keyToCol.get(key) ?? 0;
      const colsInRow = rowColCount.get(row) ?? 1;
      // Center columns horizontally
      const x = (col - (colsInRow - 1) / 2) * colSpacing;
      const y = row * rowSpacing;
      nodePositions.set(key, { x, y });
    });

    // 5. Build a canonical node-id mapping: use the first graph node ID for each key
    const keyToNodeId = new Map<string, string>();
    allKeys.forEach((key) => {
      const info = configKeyInfo.get(key);
      if (info && info.nodeIds.length > 0) {
        keyToNodeId.set(key, info.nodeIds[0]);
      } else {
        keyToNodeId.set(key, `config_${key.replace(/[^a-zA-Z0-9]/g, "_")}`);
      }
    });

    // 6. For path highlighting: precompute which edges belong to which path
    // Map: molecule-key → graph-edge-id for each step
    const pathEdgeMap = result.paths.map((path) => {
      const edgeIds: string[] = [];
      const steps = path.steps || [];
      steps.forEach((step) => {
        const fKey = molKey(step.from_config?.molecules);
        const tKey = molKey(step.to_config?.molecules);
        if (!fKey || !tKey) return;
        const fId = keyToNodeId.get(fKey);
        const tId = keyToNodeId.get(tKey);
        if (!fId || !tId) return;
        // Find matching graph edge
        const matchEdge = graphEdges.find((e) => {
          return (e.source === fId && e.target === tId) ||
                 (e.source === tId && e.target === fId);
        });
        if (matchEdge) edgeIds.push(matchEdge.id);
      });
      return edgeIds;
    });

    const pathNodeMap = result.paths.map((path) => {
      const nodeIds = new Set<string>();
      const steps = path.steps || [];
      steps.forEach((step) => {
        const fKey = molKey(step.from_config?.molecules);
        const tKey = molKey(step.to_config?.molecules);
        if (fKey) { const id = keyToNodeId.get(fKey); if (id) nodeIds.add(id); }
        if (tKey) { const id = keyToNodeId.get(tKey); if (id) nodeIds.add(id); }
      });
      return nodeIds;
    });

    return {
      configKeyInfo,
      keyToRow,
      keyToCol,
      nodePositions,
      keyToNodeId,
      pathEdgeMap,
      pathNodeMap,
      maxSteps,
      rowColCount,
      reactantKeys,
      productKeys,
    };
  }, [result]);

  // ── Build cytoscape graph ──
  const buildGraph = useCallback(() => {
    if (!result || !result.graph || result.graph.nodes.length === 0 || !cyContainerRef.current || !layoutData) return;
    const { nodes: graphNodes, edges: graphEdges } = result.graph;
    const { nodePositions, keyToNodeId, configKeyInfo, reactantKeys, productKeys } = layoutData;

    // Collect all node IDs that appear on any path
    const allPathNodeIds = new Set<string>();
    layoutData.pathNodeMap.forEach((ids) => ids.forEach((id) => allPathNodeIds.add(id)));

    // Build cytoscape nodes — only include nodes that are in our layout (on at least one path)
    const cyNodes: Array<{ data: Record<string, unknown>; classes?: string; position: { x: number; y: number } }> = [];
    const usedNodeIds = new Set<string>();

    configKeyInfo.forEach((info, key) => {
      const nodeId = keyToNodeId.get(key);
      if (!nodeId) return;
      if (usedNodeIds.has(nodeId)) return;
      usedNodeIds.add(nodeId);

      const isReactant = reactantKeys.has(key);
      const isProduct = productKeys.has(key);
      const isOnPath = allPathNodeIds.has(nodeId);

      let colorKey: "reactant" | "product" | "exploration" | "intermediate";
      if (isReactant) colorKey = "reactant";
      else if (isProduct) colorKey = "product";
      else if (isOnPath) colorKey = "exploration";
      else colorKey = "intermediate";

      // Determine label: "R" for reactant, "P" for product, distance_letters otherwise
      let displayLabel: string;
      if (isReactant) displayLabel = "R";
      else if (isProduct) displayLabel = "P";
      else displayLabel = info.label || info.smilesLabel || "";

      const pos = nodePositions.get(key) || { x: 0, y: 0 };

      cyNodes.push({
        data: {
          id: nodeId,
          label: displayLabel,
          smilesLabel: info.smilesLabel || info.label || "",
          molecules: info.molecules || [],
          source: info.source || (isReactant ? "reactant" : isProduct ? "product" : "exploration"),
          depth: info.depth,
          colorKey,
          molKey: key,
        },
        classes: `${colorKey}-node`,
        position: pos,
      });
    });

    // Build cytoscape edges — only edges between nodes that are in our layout
    const cyEdges: Array<{ data: Record<string, unknown> }> = [];
    graphEdges.forEach((e) => {
      if (!usedNodeIds.has(e.source) || !usedNodeIds.has(e.target)) return;

      const ruleStatus = e.rule_status || "normal";
      cyEdges.push({
        data: {
          id: e.id,
          source: e.source,
          target: e.target,
          label: e.label || "",
          ruleId: e.rule_id,
          direction: e.direction || "forward",
          score: e.score || 0,
          maxBondLength: e.max_bond_length || 0,
          ruleStatus,
        },
      });
    });

    // Destroy previous instance
    if (cyRef.current) {
      try { cyRef.current.destroy(); } catch { /* ignore */ }
      cyRef.current = null;
    }

    const cyInstance = cytoscape({
      container: cyContainerRef.current,
      elements: [...cyNodes, ...cyEdges],
      style: [
        // ── Base node style ──
        {
          selector: "node",
          style: {
            "shape": "ellipse",
            "width": 40,
            "height": 40,
            "background-color": COLORS.intermediate.bg,
            "border-color": COLORS.intermediate.border,
            "border-width": 3,
            "color": COLORS.intermediate.text,
            "font-size": "11px",
            "font-weight": "bold",
            "text-valign": "center",
            "text-halign": "center",
            "text-outline-color": "rgba(0,0,0,0.4)",
            "text-outline-width": 2,
            "text-wrap": "wrap",
            "text-max-width": "60px",
          },
        },
        // ── Reactant node: Orange ──
        {
          selector: "node.reactant-node",
          style: {
            "background-color": COLORS.reactant.bg,
            "border-color": COLORS.reactant.border,
            "border-width": 4,
            "color": COLORS.reactant.text,
            "width": 48,
            "height": 48,
            "font-size": "16px",
          },
        },
        // ── Product node: Red ──
        {
          selector: "node.product-node",
          style: {
            "background-color": COLORS.product.bg,
            "border-color": COLORS.product.border,
            "border-width": 4,
            "color": COLORS.product.text,
            "width": 48,
            "height": 48,
            "font-size": "16px",
          },
        },
        // ── Exploration node: Yellow ──
        {
          selector: "node.exploration-node",
          style: {
            "background-color": COLORS.exploration.bg,
            "border-color": COLORS.exploration.border,
            "border-width": 3,
            "color": COLORS.exploration.text,
            "font-size": "10px",
          },
        },
        // ── Intermediate node: Gray, faded ──
        {
          selector: "node.intermediate-node",
          style: {
            "background-color": COLORS.intermediate.bg,
            "border-color": COLORS.intermediate.border,
            "border-width": 2,
            "color": COLORS.intermediate.text,
            "font-size": "9px",
            "opacity": 0.35,
          },
        },
        // ── Base edge ──
        {
          selector: "edge",
          style: {
            "width": 2,
            "line-color": COLORS.edgeNormal,
            "target-arrow-color": COLORS.edgeNormal,
            "target-arrow-shape": "triangle",
            "curve-style": "bezier",
            "arrow-scale": 0.8,
            "label": "",
            "font-size": "8px",
            "text-rotation": "autorotate",
            "text-background-color": "#ffffff",
            "text-background-opacity": 0.85,
            "text-background-padding": "2px",
            "color": "#6b7280",
            "opacity": 0.5,
          },
        },
        // ── Normal edge: Gray ──
        {
          selector: "edge.edge-normal",
          style: {
            "line-color": COLORS.edgeNormal,
            "target-arrow-color": COLORS.edgeNormal,
            "opacity": 0.35,
          },
        },
        // ── Database-only edge: Red ──
        {
          selector: "edge.edge-database_only",
          style: {
            "line-color": COLORS.edgeDbOnly,
            "target-arrow-color": COLORS.edgeDbOnly,
            "width": 3,
            "opacity": 0.85,
          },
        },
        // ── Cross-mechanism edge: Orange ──
        {
          selector: "edge.edge-cross_mechanism",
          style: {
            "line-color": COLORS.edgeCross,
            "target-arrow-color": COLORS.edgeCross,
            "width": 3,
            "opacity": 0.85,
          },
        },
        // ── Highlighted edge (selected path) — amber, thick ──
        {
          selector: "edge.highlighted",
          style: {
            "line-color": COLORS.highlighted,
            "target-arrow-color": COLORS.highlighted,
            "width": 4,
            "opacity": 1,
            "z-index": 10,
          },
        },
        // ── Highlighted node ──
        {
          selector: "node.highlighted",
          style: {
            "border-width": 5,
            "border-color": COLORS.highlighted,
            "z-index": 10,
            "opacity": 1,
          },
        },
        // ── Dimmed (non-path elements when a path is selected) ──
        {
          selector: ".dimmed",
          style: {
            "opacity": 0.12,
          },
        },
        // ── Interaction ──
        { selector: "node:active", style: { "overlay-opacity": 0 } },
        { selector: "edge:active", style: { "overlay-opacity": 0 } },
      ],
      // Use preset positions (already set on nodes above)
      layout: { name: "preset" },
      userZoomingEnabled: true,
      userPanningEnabled: true,
      boxSelectionEnabled: false,
      autoungrabify: true,
      wheelSensitivity: 0.3,
      minZoom: 0.3,
      maxZoom: 3,
    });

    // Center the view on all nodes
    cyInstance.fit(undefined, 50);

    // Node tap handler
    cyInstance.on("tap", "node", (evt) => {
      const node = evt.target;
      setSelectedNode({
        id: node.id(),
        molecules: (node.data("molecules") as string[]) || [],
        source: (node.data("source") as string) || "",
        label: (node.data("label") as string) || (node.data("smilesLabel") as string) || "",
      });
    });

    // Highlight path function
    const highlightPath = (pathIdx: number) => {
      cyInstance.elements().removeClass("highlighted").removeClass("dimmed");

      if (pathIdx < 0 || !layoutData.pathEdgeMap[pathIdx]) {
        return;
      }

      const edgeIds = new Set(layoutData.pathEdgeMap[pathIdx]);
      const nodeIds = layoutData.pathNodeMap[pathIdx];

      // Dim everything
      cyInstance.elements().addClass("dimmed");

      // Highlight path edges
      cyInstance.edges().forEach((edge) => {
        if (edgeIds.has(edge.id())) {
          edge.removeClass("dimmed").addClass("highlighted");
        }
      });

      // Highlight path nodes
      cyInstance.nodes().forEach((node) => {
        if (nodeIds.has(node.id())) {
          node.removeClass("dimmed").addClass("highlighted");
        }
      });
    };

    (cyInstance as unknown as Record<string, unknown>)._highlightPath = highlightPath;
    cyRef.current = cyInstance;

    // Initial highlight
    highlightPath(selectedPath);
  }, [result, layoutData, selectedPath]);

  useEffect(() => {
    buildGraph();
    return () => {
      if (cyRef.current) {
        try { cyRef.current.destroy(); } catch { /* ignore */ }
        cyRef.current = null;
      }
    };
  }, [buildGraph]);

  // When selectedPath changes, re-highlight without rebuilding
  useEffect(() => {
    if (cyRef.current && layoutData) {
      const cy = cyRef.current;
      cy.elements().removeClass("highlighted").removeClass("dimmed");

      const pathIdx = selectedPath;
      if (pathIdx < 0 || !layoutData.pathEdgeMap[pathIdx]) return;

      const edgeIds = new Set(layoutData.pathEdgeMap[pathIdx]);
      const nodeIds = layoutData.pathNodeMap[pathIdx];

      cy.elements().addClass("dimmed");
      cy.edges().forEach((edge) => {
        if (edgeIds.has(edge.id())) edge.removeClass("dimmed").addClass("highlighted");
      });
      cy.nodes().forEach((node) => {
        if (nodeIds.has(node.id())) node.removeClass("dimmed").addClass("highlighted");
      });
    }
  }, [selectedPath, layoutData]);

  // ── Render ──
  if (!result) {
    return (
      <div className="text-center py-10">
        <Network className="w-8 h-8 text-gray-300 mx-auto mb-2" />
        <p className="text-gray-500 text-sm">Run a bidirectional search to visualize the mechanism graph.</p>
      </div>
    );
  }

  if (!result.paths || result.paths.length === 0) {
    const hasGraph = result.graph && result.graph.nodes && result.graph.nodes.length > 0;
    const hasStats = result.stats && (
      ((result.stats.forward_configs ?? 0) > 1 || (result.stats.backward_configs ?? 0) > 1)
    );

    return (
      <div className="space-y-4">
        <Alert>
          <AlertCircle className="h-4 w-4 text-amber-500" />
          <AlertTitle className="text-sm">No Complete Paths Found</AlertTitle>
          <AlertDescription className="text-xs">
            {hasStats
              ? "The search explored configurations but could not find a connecting mechanism path between reactants and products."
              : "The bidirectional search could not find a connecting mechanism path."}
          </AlertDescription>
        </Alert>
        {hasStats && (
          <div className="grid grid-cols-2 sm:grid-cols-4 gap-2">
            {[
              { label: "Forward Configs", value: result.stats.forward_configs },
              { label: "Backward Configs", value: result.stats.backward_configs },
              { label: "Graph Nodes", value: result.stats.total_nodes || 0 },
              { label: "Graph Edges", value: result.stats.total_edges || 0 },
            ].map((s) => (
              <div key={s.label} className="rounded-lg bg-gray-50 border border-gray-100 p-2.5 text-center">
                <div className="text-xs text-gray-500">{s.label}</div>
                <div className="text-sm font-mono font-semibold text-gray-800">{s.value}</div>
              </div>
            ))}
          </div>
        )}
        {result.stats && (result.stats.search_time || result.search_time) && (
          <div className="text-xs text-gray-400 text-center">
            Search completed in {(result.stats.search_time || result.search_time)?.toFixed(1)}s using {result.stats.total_rules} rules
          </div>
        )}
      </div>
    );
  }

  const currentPath = result.paths[selectedPath];

  return (
    <div className="space-y-4">
      {/* Path selector */}
      <div className="flex items-center gap-2 flex-wrap">
        <span className="text-xs text-gray-500 uppercase tracking-wider font-medium">Path</span>
        {result.paths.map((p, i) => (
          <Button key={p.path_id ?? i} variant={selectedPath === i ? "default" : "outline"} size="sm"
            className={`h-7 text-xs px-2.5 ${selectedPath === i ? "bg-amber-600 hover:bg-amber-700 text-white" : "border-amber-200 text-amber-700 hover:bg-amber-50"}`}
            onClick={() => setSelectedPath(i)}>
            Path {i + 1}
            <Badge variant="secondary" className={`ml-1 text-xs ${selectedPath === i ? "bg-amber-500 text-white" : "bg-amber-50 text-amber-700"}`}>{p.total_steps} steps</Badge>
          </Button>
        ))}
      </div>

      {/* Configuration Graph */}
      <div className="rounded-xl border border-gray-200 bg-white overflow-hidden">
        <div className="flex items-center justify-between px-3 py-2 bg-gray-50 border-b border-gray-200">
          <div className="flex items-center gap-2">
            <Network className="w-3.5 h-3.5 text-gray-600" />
            <span className="text-xs font-medium text-gray-700">Configuration Graph</span>
          </div>
          <div className="flex items-center gap-3 text-xs text-gray-600">
            <span className="flex items-center gap-1">
              <span className="w-3 h-3 rounded-full bg-orange-500 inline-block border-2 border-orange-600" /> Reactant (R)
            </span>
            <span className="flex items-center gap-1">
              <span className="w-3 h-3 rounded-full bg-red-500 inline-block border-2 border-red-600" /> Product (P)
            </span>
            <span className="flex items-center gap-1">
              <span className="w-3 h-3 rounded-full bg-yellow-500 inline-block border-2 border-yellow-600" /> Exploration
            </span>
            <span className="flex items-center gap-1">
              <span className="w-3 h-3 rounded-full bg-gray-400 inline-block border-2 border-gray-500" /> Other
            </span>
          </div>
        </div>
        {/* Edge color legend */}
        <div className="flex items-center gap-4 px-3 py-1.5 bg-gray-50/50 border-b border-gray-100 text-xs text-gray-500">
          <span className="font-medium">Edges:</span>
          <span className="flex items-center gap-1">
            <span className="w-5 h-0.5 bg-gray-400 inline-block rounded" /> Normal
          </span>
          <span className="flex items-center gap-1">
            <span className="w-5 h-0.5 bg-red-500 inline-block rounded" /> Database-only rule
          </span>
          <span className="flex items-center gap-1">
            <span className="w-5 h-0.5 bg-orange-500 inline-block rounded" /> Cross-mechanism rule
          </span>
          <span className="flex items-center gap-1">
            <span className="w-5 h-0.5 bg-amber-500 inline-block rounded" /> Selected path
          </span>
        </div>
        <div className="relative">
          <div
            ref={cyContainerRef}
            style={{
              width: "100%",
              height: isExpanded ? "520px" : "340px",
              minHeight: "280px",
            }}
          />
          <Button
            variant="outline"
            size="sm"
            className="absolute top-2 right-2 h-6 w-6 p-0 bg-white/80 border-gray-200 shadow-sm"
            onClick={() => setIsExpanded(!isExpanded)}
          >
            {isExpanded ? <Minimize2 className="w-3 h-3" /> : <Maximize2 className="w-3 h-3" />}
          </Button>
        </div>
      </div>

      {/* Path steps detail */}
      {currentPath && (
        <div className="space-y-2">
          <div className="flex items-center gap-2">
            <GitBranch className="w-3.5 h-3.5 text-amber-600" />
            <span className="text-xs font-medium text-gray-700">Path {selectedPath + 1} — {currentPath.total_steps ?? 0} step{(currentPath.total_steps ?? 0) > 1 ? "s" : ""}</span>
            {currentPath.total_score != null && currentPath.total_score > 0 && (
              <Badge variant="outline" className="text-[10px] px-1.5 py-0 border-gray-200 text-gray-500 bg-gray-50">
                total score: {currentPath.total_score.toFixed(1)}
              </Badge>
            )}
          </div>
          {currentPath.steps.map((step, sIdx) => (
            <div key={`${selectedPath}-${sIdx}`} className="relative ml-3 pl-5 border-l-2 border-amber-200">
              <div className="absolute -left-[5px] top-2 w-2.5 h-2.5 rounded-full bg-amber-400" />
              <div className="rounded-lg border border-gray-100 bg-white p-3 mb-2">
                <div className="flex items-center gap-2 mb-1.5 flex-wrap">
                  <Badge variant="outline" className="text-xs border-amber-200 text-amber-600 font-mono">Step {step.step_num}</Badge>
                  <Badge variant="secondary" className={`text-xs ${(step.direction || "forward") === "forward" ? "bg-orange-50 text-orange-700" : "bg-red-50 text-red-700"}`}>{(step.direction || "forward") === "forward" ? "FWD" : "BWD"}</Badge>
                  {step.reaction_type && step.reaction_type !== "general_reaction" && (
                    <Badge variant="secondary" className="bg-yellow-50 text-yellow-700 text-xs border border-yellow-200">
                      <Zap className="w-2.5 h-2.5 mr-0.5" />{step.reaction_type.replace(/_/g, " ")}
                    </Badge>
                  )}
                  {step.arrows && step.arrows.length > 0 && (
                    <Badge variant="secondary" className="bg-cyan-50 text-cyan-700 text-xs border border-cyan-200">
                      <ArrowRightLeft className="w-2.5 h-2.5 mr-0.5" />{step.arrows.length} arrow{step.arrows.length > 1 ? "s" : ""}
                    </Badge>
                  )}
                  {step.score != null && step.score > 0 && (
                    <Badge
                      variant="outline"
                      className="text-[10px] px-1.5 py-0 border-gray-200 text-gray-500 bg-gray-50"
                      title={`Step score: ${step.score}`}
                    >
                      score: {step.score.toFixed(1)}
                    </Badge>
                  )}
                  {/* Bond length badge */}
                  {step.max_bond_length != null && step.max_bond_length > 0 ? (
                    <Badge
                      variant="outline"
                      className={`flex items-center gap-1 text-[10px] px-1.5 py-0 ${
                        step.max_bond_length > 6
                          ? "border-red-300 text-red-600 bg-red-50"
                          : step.max_bond_length > 4
                            ? "border-amber-300 text-amber-700 bg-amber-50"
                            : "border-emerald-300 text-emerald-700 bg-emerald-50"
                      }`}
                      title={`Maximum new bond length: ${step.max_bond_length.toFixed(2)} Å`}
                    >
                      <Ruler className="w-2.5 h-2.5" />
                      {step.max_bond_length.toFixed(2)} Å
                    </Badge>
                  ) : (
                    <span className="text-[10px] text-gray-300 italic" title="No bond length data">— no bond data —</span>
                  )}
                </div>
                <div className="space-y-1.5">
                  {(step.from_config?.molecules?.length || step.to_config?.molecules?.length) ? (
                    <div>
                      <div className="flex items-center gap-2 min-h-[64px]">
                        {step.from_config?.molecules && step.from_config.molecules.length > 0 ? (
                          <div className="flex-1 flex flex-wrap gap-1 justify-end">
                            {step.from_config.molecules.map((mol, mi) => (
                              <MoleculeSVG key={mi} smiles={mol} size={60} className="border border-orange-100 rounded-md" />
                            ))}
                          </div>
                        ) : (
                          <div className="flex-1 text-right text-xs text-gray-300 italic">—</div>
                        )}
                        <div className="shrink-0 flex flex-col items-center gap-0.5 px-1">
                          <ArrowRight className="w-4 h-4 text-gray-400" />
                          <span className="text-[8px] text-gray-400 text-center max-w-[80px] leading-tight line-clamp-2">
                            {step.rule_name || step.reaction_smarts?.slice(0, 50) || "Rule applied"}
                          </span>
                        </div>
                        {step.to_config?.molecules && step.to_config.molecules.length > 0 ? (
                          <div className="flex-1 flex flex-wrap gap-1 justify-start">
                            {step.to_config.molecules.map((mol, mi) => (
                              <MoleculeSVG key={mi} smiles={mol} size={60} className="border border-red-100 rounded-md" />
                            ))}
                          </div>
                        ) : (
                          <div className="flex-1 text-xs text-gray-300 italic">—</div>
                        )}
                      </div>
                    </div>
                  ) : (
                    <div className="text-xs text-gray-400 italic">No molecule data available</div>
                  )}
                  {step.arrows && step.arrows.length > 0 && (
                    <div className="mt-2 rounded-md bg-cyan-50/50 border border-cyan-100 p-2">
                      <div className="flex items-center gap-1 mb-1">
                        <ArrowRightLeft className="w-3 h-3 text-cyan-600" />
                        <span className="text-xs font-medium text-cyan-700">Curly Arrow Mechanism</span>
                      </div>
                      <div className="space-y-0.5">
                        {step.arrows.map((arrow, ai) => (
                          <div key={ai} className="text-xs text-gray-600 flex items-center gap-1">
                            <span className="text-cyan-500">→</span>
                            <span>atoms [{arrow.source_atoms.join(", ")}] → [{arrow.target_atoms.join(", ")}]</span>
                            <Badge variant="outline" className="text-xs px-1 py-0 border-cyan-200 text-cyan-600">{arrow.electrons}e⁻</Badge>
                          </div>
                        ))}
                      </div>
                    </div>
                  )}
                </div>
              </div>
            </div>
          ))}
        </div>
      )}

      {/* Node detail panel */}
      {selectedNode && (
        <div className="rounded-lg border border-gray-100 bg-gray-50/50 p-3">
          <div className="flex items-center justify-between mb-2">
            <span className="text-xs font-medium text-gray-600 uppercase tracking-wider">Node Detail</span>
            <Button variant="ghost" size="sm" className="h-5 w-5 p-0 text-gray-400 hover:text-gray-600" onClick={() => setSelectedNode(null)}><RotateCcw className="w-3 h-3" /></Button>
          </div>
          <div className="flex items-center gap-2 mb-1">
            <Badge variant="outline" className="text-xs">{selectedNode.source}</Badge>
            <span className="text-xs text-gray-400 font-mono">{selectedNode.label || selectedNode.id}</span>
          </div>
          <div className="space-y-2 mt-2">
            <span className="text-xs text-gray-500 uppercase tracking-wider">Molecules ({selectedNode.molecules.length})</span>
            <div className="flex flex-wrap gap-1">
              {selectedNode.molecules.slice(0, 6).map((mol, i) => (
                <MoleculeSVG key={i} smiles={mol} size={70} className="border border-gray-100 rounded-md" zoomable />
              ))}
            </div>
            {selectedNode.molecules.length > 6 && <p className="text-xs text-gray-400">+{selectedNode.molecules.length - 6} more molecules</p>}
          </div>
        </div>
      )}

      {/* Stats footer */}
      {result.stats && (
        <div className="flex items-center gap-4 text-xs text-gray-500 flex-wrap">
          <span className="flex items-center gap-1"><GitBranch className="w-3 h-3" /> {result.stats.paths_found} path(s)</span>
          <span className="flex items-center gap-1"><Zap className="w-3 h-3" /> {result.stats.total_rules} rules</span>
          <span className="flex items-center gap-1"><Clock className="w-3 h-3" /> {(result.total_time || result.search_time || 0).toFixed(1)}s</span>
          {currentPath.merge_similarity != null && (
            <span className="flex items-center gap-1"><Layers className="w-3 h-3 text-amber-400" /> merge sim: {(currentPath.merge_similarity * 100).toFixed(0)}%</span>
          )}
        </div>
      )}
    </div>
  );
}
