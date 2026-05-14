'use client';

import React, { useEffect, useRef, useState, useCallback } from 'react';

// ── Types ────────────────────────────────────────────────────────────────
interface HighlightedResidue {
  chain: string;
  res_num: number;
  res_name?: string;
}

interface PdbViewerWebGLProps {
  pdbText: string;
  highlightedResidues?: HighlightedResidue[];
  selectedChain?: string;
  height?: string;
  onResidueClick?: (residue: { chain: string; res_num: number; res_name: string }) => void;
  className?: string;
}

// ── Display modes ────────────────────────────────────────────────────────
type DisplayMode = 'cartoon' | 'ribbon' | 'tube' | 'sphere' | 'stick' | 'line' | 'cross' | 'ballstick' | 'surface';

interface ModeOption {
  id: DisplayMode;
  label: string;
  icon: string;
  tooltip: string;
}

const DISPLAY_MODES: ModeOption[] = [
  { id: 'cartoon', label: 'Cartoon', icon: '🎞️', tooltip: 'Cartoon / ribbon diagram (PyMOL-style)' },
  { id: 'ribbon', label: 'Ribbon', icon: '🎀', tooltip: 'Flat ribbon representation' },
  { id: 'tube', label: 'Tube', icon: '🧬', tooltip: 'Thick tube trace' },
  { id: 'stick', label: 'Stick', icon: '⚛️', tooltip: 'Stick model' },
  { id: 'ballstick', label: 'Ball&Stick', icon: '🔮', tooltip: 'Ball and stick model' },
  { id: 'sphere', label: 'Sphere', icon: '🌐', tooltip: 'Spacefill / CPK spheres' },
  { id: 'line', label: 'Line', icon: '📐', tooltip: 'Wireframe lines' },
  { id: 'cross', label: 'Cross', icon: '✚', tooltip: 'Cross representation' },
  { id: 'surface', label: 'Surface', icon: '🫧', tooltip: 'Molecular surface (SES)' },
];

// ── Color schemes ────────────────────────────────────────────────────────
type ColorScheme = 'chain' | 'spectrum' | 'ss' | 'residue';

interface ColorOption {
  id: ColorScheme;
  label: string;
  icon: string;
  tooltip: string;
}

const COLOR_SCHEMES: ColorOption[] = [
  { id: 'chain', label: 'Chain', icon: '🔗', tooltip: 'Color by chain' },
  { id: 'spectrum', label: 'Spectrum', icon: '🌈', tooltip: 'Rainbow spectrum by residue position' },
  { id: 'ss', label: 'SS', icon: '🧱', tooltip: 'Color by secondary structure' },
  { id: 'residue', label: 'Residue', icon: '🧪', tooltip: 'Color by residue type' },
];

// ── 3Dmol.js dynamic import helper ───────────────────────────────────────
type $3DmolType = any;
type GLViewerType = $3DmolType;

let $3DmolPromise: Promise<$3DmolType> | null = null;

function load3Dmol(): Promise<$3DmolType> {
  if ($3DmolPromise) return $3DmolPromise;
  $3DmolPromise = new Promise<$3DmolType>((resolve, reject) => {
    const check = () => {
      const g = globalThis as any;
      if (g.$3Dmol && typeof g.$3Dmol.createViewer === 'function') {
        resolve(g.$3Dmol);
      } else {
        setTimeout(check, 200);
      }
    };
    check();
    setTimeout(() => reject(new Error('3Dmol.js loading timed out')), 10000);
  });
  return $3DmolPromise;
}

// ══════════════════════════════════════════════════════════════════════════
//  COMPONENT
// ══════════════════════════════════════════════════════════════════════════
export default function PdbViewerWebGL({
  pdbText,
  highlightedResidues = [],
  selectedChain,
  height = '650px',
  onResidueClick,
  className = '',
}: PdbViewerWebGLProps) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<GLViewerType>(null);
  const $3DmolRef = useRef<$3DmolType>(null);
  const modelRef = useRef<number | null>(null);

  const [status, setStatus] = useState<'loading' | 'ready' | 'error'>('loading');
  const [errorMsg, setErrorMsg] = useState('');
  const [displayMode, setDisplayMode] = useState<DisplayMode>('cartoon');
  const [colorScheme, setColorScheme] = useState<ColorScheme>('chain');
  const [spinning, setSpinning] = useState(false);
  const [showLabels, setShowLabels] = useState(false);
  const [darkBg, setDarkBg] = useState(false);

  // Refs for accessing latest values in callbacks (synced via useEffect)
  const displayModeRef = useRef(displayMode);
  const colorSchemeRef = useRef(colorScheme);
  const spinningRef = useRef(spinning);
  const showLabelsRef = useRef(showLabels);
  const darkBgRef = useRef(darkBg);
  const highlightedResiduesRef = useRef(highlightedResidues);
  const selectedChainRef = useRef(selectedChain);
  const onResidueClickRef = useRef(onResidueClick);

  // Keep refs in sync with state/props
  useEffect(() => {
    displayModeRef.current = displayMode;
    colorSchemeRef.current = colorScheme;
    spinningRef.current = spinning;
    showLabelsRef.current = showLabels;
    darkBgRef.current = darkBg;
    highlightedResiduesRef.current = highlightedResidues;
    selectedChainRef.current = selectedChain;
    onResidueClickRef.current = onResidueClick;
  }, [displayMode, colorScheme, spinning, showLabels, darkBg, highlightedResidues, selectedChain, onResidueClick]);

  // ── Apply style to the model ───────────────────────────────────────────
  const applyStyle = useCallback(() => {
    const viewer = viewerRef.current;
    const $3Dmol = $3DmolRef.current;
    const modelId = modelRef.current;
    if (!viewer || !$3Dmol || modelId === null) return;

    const model = viewer.getModel(modelId);
    if (!model) return;

    const currentScheme = colorSchemeRef.current;
    const currentMode = displayModeRef.current;
    const currentSelectedChain = selectedChainRef.current;
    const currentHighlight = highlightedResiduesRef.current;

    // Build colorscheme config
    let colorscheme: string | undefined = undefined;
    if (currentScheme === 'spectrum') {
      colorscheme = 'spectrum';
    } else if (currentScheme === 'ss') {
      colorscheme = 'ssPyMol';
    } else if (currentScheme === 'residue') {
      colorscheme = 'amino';
    }
    // 'chain' uses the default chain coloring in 3dmol (no explicit colorscheme needed)

    // Build the style spec based on display mode
    const styleSpec: Record<string, $3DmolType> = {};

    if (currentMode === 'cartoon') {
      styleSpec.cartoon = colorscheme ? { colorscheme } : {};
    } else if (currentMode === 'ribbon') {
      styleSpec.cartoon = { ...({ colorscheme } as object), style: 'ribbon' };
    } else if (currentMode === 'tube') {
      styleSpec.cartoon = { ...({ colorscheme } as object), style: 'tube' };
    } else if (currentMode === 'sphere') {
      styleSpec.sphere = colorscheme ? { colorscheme, scale: 0.5 } : { scale: 0.5 };
    } else if (currentMode === 'stick') {
      styleSpec.stick = colorscheme ? { colorscheme, radius: 0.1 } : { radius: 0.1 };
    } else if (currentMode === 'line') {
      styleSpec.line = colorscheme ? { colorscheme } : {};
    } else if (currentMode === 'cross') {
      styleSpec.cross = colorscheme ? { colorscheme, linewidth: 2 } : { linewidth: 2 };
    } else if (currentMode === 'ballstick') {
      styleSpec.stick = colorscheme ? { colorscheme, radius: 0.08 } : { radius: 0.08 };
      styleSpec.sphere = colorscheme ? { colorscheme, scale: 0.25 } : { scale: 0.25 };
    } else if (currentMode === 'surface') {
      // Surface mode: show cartoon + surface
      styleSpec.cartoon = colorscheme ? { colorscheme, opacity: 0.3 } : { opacity: 0.3 };
    }

    // Handle chain selection - dim non-selected chains
    if (currentSelectedChain) {
      model.setStyle({ chain: currentSelectedChain }, styleSpec);
      // Dim other chains
      const dimSpec: Record<string, $3DmolType> = {};
      if (currentMode === 'cartoon' || currentMode === 'ribbon' || currentMode === 'tube') {
        dimSpec.cartoon = { opacity: 0.15, colorscheme: 'chain' };
      } else if (currentMode === 'sphere') {
        dimSpec.sphere = { opacity: 0.1 };
      } else if (currentMode === 'stick' || currentMode === 'ballstick') {
        dimSpec.stick = { opacity: 0.1 };
      } else if (currentMode === 'line') {
        dimSpec.line = { opacity: 0.08 };
      } else if (currentMode === 'cross') {
        dimSpec.cross = { opacity: 0.1 };
      } else if (currentMode === 'surface') {
        dimSpec.cartoon = { opacity: 0.05 };
      }
      model.setStyle({ chain: { $ne: currentSelectedChain } }, dimSpec);
    } else {
      model.setStyle({}, styleSpec);
    }

    // Handle highlighted residues
    if (currentHighlight.length > 0) {
      for (const hl of currentHighlight) {
        const hlSel = { chain: hl.chain, resi: hl.res_num };
        const hlStyle: Record<string, $3DmolType> = {};
        if (currentMode === 'cartoon' || currentMode === 'ribbon' || currentMode === 'tube') {
          hlStyle.cartoon = { color: '#EF4444', opacity: 1.0 };
        } else if (currentMode === 'sphere' || currentMode === 'ballstick') {
          hlStyle.sphere = { color: '#EF4444' };
        } else if (currentMode === 'stick') {
          hlStyle.stick = { color: '#EF4444' };
        } else if (currentMode === 'line') {
          hlStyle.line = { color: '#EF4444' };
        } else if (currentMode === 'cross') {
          hlStyle.cross = { color: '#EF4444' };
        } else {
          hlStyle.cartoon = { color: '#EF4444' };
          hlStyle.sphere = { color: '#EF4444', scale: 0.3 };
        }
        model.setStyle(hlSel, hlStyle);
      }
    }

    // Surface mode: add molecular surface
    viewer.removeAllSurfaces();
    if (currentMode === 'surface') {
      const surfSel = currentSelectedChain ? { chain: currentSelectedChain } : {};
      viewer.addSurface($3Dmol.SurfaceType.SES, {
        opacity: 0.7,
        colorscheme: currentScheme === 'spectrum' ? 'spectrum'
          : currentScheme === 'ss' ? 'ssPyMol'
          : currentScheme === 'residue' ? 'amino'
          : 'chain',
      }, surfSel);
    }

    // Labels
    viewer.removeAllLabels();
    if (showLabelsRef.current) {
      const labelSel = currentSelectedChain ? { chain: currentSelectedChain, hetflag: false } : { hetflag: false };
      viewer.addResLabels(labelSel, {
        font: '12px sans-serif',
        fontColor: darkBgRef.current ? '#ffffff' : '#000000',
        backgroundColor: darkBgRef.current ? 'rgba(0,0,0,0.5)' : 'rgba(255,255,255,0.85)',
        backgroundOpacity: 0.8,
        showBackground: true,
        inFront: true,
      });
    }

    viewer.render();
  }, []);

  // ── Initialize 3Dmol viewer ────────────────────────────────────────────
  useEffect(() => {
    if (typeof window === 'undefined') return;
    const container = containerRef.current;
    if (!container) return;

    let destroyed = false;

    const init = async () => {
      try {
        const $3Dmol = await load3Dmol();
        if (destroyed || !container) return;
        $3DmolRef.current = $3Dmol;

        // Clean up any existing canvas children
        while (container.firstChild) {
          container.removeChild(container.firstChild);
        }

        const viewer = $3Dmol.createViewer(container, {
          backgroundColor: '#ffffff',
          antialias: true,
          disableFog: false,
          disableOutline: false,
          cartoonQuality: 10,
          lowerZoomLimit: 5,
          upperZoomLimit: 800,
        });

        if (destroyed) {
          // cleanup viewer
          viewer.removeAllModels();
          viewer.removeAllSurfaces();
          viewer.removeAllLabels();
          return;
        }

        viewerRef.current = viewer;

        // Load PDB data
        if (!pdbText.trim()) {
          setStatus('error');
          setErrorMsg('No PDB data provided');
          return;
        }

        const model = viewer.addModel(pdbText, 'pdb');
        modelRef.current = model.id;

        // Enable click handling for residue picking
        if (onResidueClickRef.current) {
          model.setClickable({}, true, (atom: $3DmolType) => {
            if (atom && atom.chain && atom.resi !== undefined && atom.resn) {
              onResidueClickRef.current?.({
                chain: atom.chain,
                res_num: atom.resi,
                res_name: atom.resn,
              });
            }
          });
        }

        viewer.zoomTo();
        applyStyle();
        setStatus('ready');
      } catch (err) {
        if (!destroyed) {
          setStatus('error');
          setErrorMsg(`Failed to initialize WebGL viewer: ${err instanceof Error ? err.message : String(err)}`);
        }
      }
    };

    init();

    return () => {
      destroyed = true;
      if (viewerRef.current) {
        try {
          viewerRef.current.spin(false);
          viewerRef.current.removeAllModels();
          viewerRef.current.removeAllSurfaces();
          viewerRef.current.removeAllLabels();
        } catch {
          // ignore
        }
        viewerRef.current = null;
      }
      modelRef.current = null;
      // Clean up container
      if (container) {
        while (container.firstChild) {
          container.removeChild(container.firstChild);
        }
      }
    };
  }, []);

  // ── Re-apply style on mode/scheme/highlight/chain changes ──────────────
  useEffect(() => {
    if (status !== 'ready') return;
    applyStyle();
  }, [displayMode, colorScheme, highlightedResidues, selectedChain, showLabels, darkBg, status, applyStyle]);

  // ── Reload on pdbText change ───────────────────────────────────────────
  useEffect(() => {
    const viewer = viewerRef.current;
    const $3Dmol = $3DmolRef.current;
    if (!viewer || !$3Dmol || !pdbText.trim()) return;

    // Remove old model
    viewer.removeAllSurfaces();
    viewer.removeAllLabels();

    if (modelRef.current !== null) {
      viewer.removeModel(modelRef.current);
      modelRef.current = null;
    }

    const model = viewer.addModel(pdbText, 'pdb');
    modelRef.current = model.id;

    // Re-enable click handling
    if (onResidueClickRef.current) {
      model.setClickable({}, true, (atom: $3DmolType) => {
        if (atom && atom.chain && atom.resi !== undefined && atom.resn) {
          onResidueClickRef.current?.({
            chain: atom.chain,
            res_num: atom.resi,
            res_name: atom.resn,
          });
        }
      });
    }

    viewer.zoomTo();
    applyStyle();
  }, [pdbText]);

  // ── Resize observer ────────────────────────────────────────────────────
  useEffect(() => {
    const container = containerRef.current;
    if (!container || typeof window === 'undefined') return;

    let resizeTimeout: ReturnType<typeof setTimeout> | null = null;

    const observer = new ResizeObserver(() => {
      if (resizeTimeout) clearTimeout(resizeTimeout);
      resizeTimeout = setTimeout(() => {
        const viewer = viewerRef.current;
        if (viewer) {
          try {
            viewer.resize();
            viewer.render();
          } catch {
            // ignore resize errors
          }
        }
      }, 100);
    });

    observer.observe(container);
    return () => {
      observer.disconnect();
      if (resizeTimeout) clearTimeout(resizeTimeout);
    };
  }, []);

  // ── Handlers ───────────────────────────────────────────────────────────
  const handleModeChange = useCallback((mode: DisplayMode) => {
    setDisplayMode(mode);
  }, []);

  const handleColorSchemeChange = useCallback((scheme: ColorScheme) => {
    setColorScheme(scheme);
  }, []);

  const handleToggleSpin = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    if (spinningRef.current) {
      viewer.spin(false);
      setSpinning(false);
    } else {
      viewer.spin(true);
      setSpinning(true);
    }
  }, []);

  const handleResetView = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    viewer.spin(false);
    setSpinning(false);
    viewer.zoomTo();
    viewer.render();
  }, []);

  const handleScreenshot = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    try {
      const imgData = viewer.pngURI();
      if (imgData) {
        const a = document.createElement('a');
        a.download = `structure_${Date.now()}.png`;
        a.href = imgData;
        a.click();
      }
    } catch {
      // fallback: capture the container's canvas
      const container = containerRef.current;
      if (container) {
        const canvas = container.querySelector('canvas');
        if (canvas) {
          const a = document.createElement('a');
          a.download = `structure_${Date.now()}.png`;
          a.href = canvas.toDataURL('image/png');
          a.click();
        }
      }
    }
  }, []);

  const handleToggleLabels = useCallback(() => {
    setShowLabels((p) => !p);
  }, []);

  const handleToggleBg = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    const newDark = !darkBgRef.current;
    setDarkBg(newDark);
    viewer.setBackgroundColor(newDark ? '#1a1a2e' : '#ffffff');
    applyStyle();
  }, [applyStyle]);

  const isReady = status === 'ready';

  // ── Button style helper ───────────────────────────────────────────────
  const btnBase: React.CSSProperties = {
    display: 'inline-flex', alignItems: 'center', gap: 3,
    padding: '4px 8px', fontSize: 11, lineHeight: 1,
    borderRadius: 6, cursor: 'pointer', whiteSpace: 'nowrap',
    backdropFilter: 'blur(8px)', boxShadow: '0 1px 3px rgba(0,0,0,0.08)',
    transition: 'all 0.15s ease',
    border: '1px solid rgba(0,0,0,0.1)',
    backgroundColor: 'rgba(255,255,255,0.92)',
    color: '#6b7280', fontWeight: 400,
  };

  const btnActive = (active: boolean, color: string, activeColor: string): React.CSSProperties =>
    active
      ? { ...btnBase, border: `1.5px solid ${color}`, backgroundColor: `${color}10`, color: activeColor, fontWeight: 600 }
      : btnBase;

  // ── Render ─────────────────────────────────────────────────────────────
  return (
    <div
      className={`pdb-viewer-webgl ${className}`}
      style={{
        width: '100%', height, position: 'relative',
        borderRadius: '0.5rem', overflow: 'hidden',
        border: '1px solid #e5e7eb',
        background: darkBg ? '#1a1a2e' : '#ffffff',
      }}
    >
      {/* 3Dmol container */}
      <div ref={containerRef} style={{ width: '100%', height: '100%' }} />

      {/* Loading / Error overlay */}
      {!isReady && (
        <div style={{
          position: 'absolute', inset: 0,
          display: 'flex', alignItems: 'center', justifyContent: 'center',
          background: 'rgba(255,255,255,0.95)', zIndex: 10,
        }}>
          {status === 'loading' ? (
            <div style={{ textAlign: 'center', color: '#6b7280' }}>
              <div style={{ fontSize: 28, marginBottom: 8 }}>🧬</div>
              <div style={{ fontSize: 13, fontWeight: 500 }}>Loading WebGL Viewer...</div>
              <div style={{ fontSize: 11, color: '#9ca3af', marginTop: 4 }}>Initializing 3Dmol.js</div>
            </div>
          ) : (
            <div style={{ textAlign: 'center', color: '#EF4444', fontSize: 13, maxWidth: 320, padding: 16 }}>
              <div style={{ fontSize: 28, marginBottom: 8 }}>⚠️</div>
              <div>{errorMsg}</div>
            </div>
          )}
        </div>
      )}

      {/* Controls toolbar */}
      {isReady && (
        <>
          {/* Top toolbar */}
          <div style={{
            position: 'absolute', top: 8, left: 8, right: 8,
            display: 'flex', alignItems: 'center', gap: 5, flexWrap: 'wrap', zIndex: 5,
          }}>
            {/* Display mode buttons */}
            {DISPLAY_MODES.map((m) => (
              <button
                key={m.id}
                title={m.tooltip}
                onClick={(e) => { e.stopPropagation(); handleModeChange(m.id); }}
                style={btnActive(displayMode === m.id, '#10B981', '#059669')}
              >
                <span style={{ fontSize: 12 }}>{m.icon}</span>
                <span>{m.label}</span>
              </button>
            ))}

            {/* Separator */}
            <div style={{ width: 1, height: 20, backgroundColor: 'rgba(0,0,0,0.1)', margin: '0 2px' }} />

            {/* Color scheme buttons */}
            {COLOR_SCHEMES.map((c) => (
              <button
                key={c.id}
                title={c.tooltip}
                onClick={(e) => { e.stopPropagation(); handleColorSchemeChange(c.id); }}
                style={btnActive(colorScheme === c.id, '#8B5CF6', '#7C3AED')}
              >
                <span style={{ fontSize: 12 }}>{c.icon}</span>
                <span>{c.label}</span>
              </button>
            ))}

            {/* Separator */}
            <div style={{ width: 1, height: 20, backgroundColor: 'rgba(0,0,0,0.1)', margin: '0 2px' }} />

            {/* Toggle buttons */}
            <button
              title={showLabels ? 'Hide residue labels' : 'Show residue labels'}
              onClick={(e) => { e.stopPropagation(); handleToggleLabels(); }}
              style={btnActive(showLabels, '#F59E0B', '#D97706')}
            >
              <span style={{ fontSize: 12 }}>🔢</span>
              <span>{showLabels ? 'Labels On' : 'Labels'}</span>
            </button>

            <button
              title={spinning ? 'Stop spinning' : 'Auto-rotate'}
              onClick={(e) => { e.stopPropagation(); handleToggleSpin(); }}
              style={btnActive(spinning, '#10B981', '#059669')}
            >
              <span style={{ fontSize: 12 }}>🔄</span>
              <span>Spin{spinning ? ' On' : ''}</span>
            </button>

            <button
              title="Reset view"
              onClick={(e) => { e.stopPropagation(); handleResetView(); }}
              style={btnBase}
            >
              <span style={{ fontSize: 12 }}>🎯</span><span>Reset</span>
            </button>

            <button
              title="Save screenshot (PNG)"
              onClick={(e) => { e.stopPropagation(); handleScreenshot(); }}
              style={btnBase}
            >
              <span style={{ fontSize: 12 }}>📷</span><span>PNG</span>
            </button>

            <button
              title={darkBg ? 'Switch to light background' : 'Switch to dark background'}
              onClick={(e) => { e.stopPropagation(); handleToggleBg(); }}
              style={btnActive(darkBg, '#6366F1', '#4F46E5')}
            >
              <span style={{ fontSize: 12 }}>{darkBg ? '☀️' : '🌙'}</span>
              <span>{darkBg ? 'Light' : 'Dark'}</span>
            </button>
          </div>

          {/* Bottom info bar */}
          <div style={{
            position: 'absolute', bottom: 6, left: 8,
            display: 'flex', alignItems: 'center', gap: 6, zIndex: 5,
          }}>
            <div style={{
              ...btnBase, cursor: 'default', opacity: 0.75, fontSize: 10,
              padding: '3px 6px',
            }}>
              <span style={{ color: '#10B981', fontWeight: 600 }}>3Dmol.js</span>
              <span style={{ color: '#9ca3af' }}>·</span>
              <span>{displayMode}</span>
              <span style={{ color: '#9ca3af' }}>·</span>
              <span>{colorScheme}</span>
            </div>
          </div>
        </>
      )}
    </div>
  );
}
