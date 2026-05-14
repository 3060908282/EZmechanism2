"use client";

import React, { useState, useEffect, useRef } from "react";
import { Loader2, ZoomIn } from "lucide-react";
import { apiCall } from "@/lib/api";
import ZoomModal from "@/components/mechanism/ZoomModal";

// ── Module-level SVG cache: SMILES → SVG string ──
// Avoids duplicate backend API calls for the same molecule.
const _svgCache = new Map<string, string>();
const _pendingCache = new Map<string, Promise<string | null>>();

export function clearMoleculeSvgCache() {
  _svgCache.clear();
  _pendingCache.clear();
}

// Cache version — bump to force re-fetch (e.g. when backend resolution changes)
const _CACHE_VERSION = "v5";

/** Build a stable cache key from SMILES + optional label/idx */
function buildCacheKey(smiles: string, label?: string, idx?: number): string {
  return JSON.stringify({ v: _CACHE_VERSION, s: smiles, l: label, i: idx });
}

/** Fetch SVG for a SMILES with cache dedup */
async function fetchSvg(smiles: string, label?: string, idx?: number): Promise<string | null> {
  if (!smiles) return null;

  const cacheKey = buildCacheKey(smiles, label, idx);
  if (_svgCache.has(cacheKey)) return _svgCache.get(cacheKey)!;

  // Dedup concurrent identical requests
  if (_pendingCache.has(cacheKey)) return _pendingCache.get(cacheKey) ?? null;

  const body: Record<string, unknown> = { smiles };
  if (label) body.label = label;
  if (idx !== undefined) body.idx = idx;
  // Always request default 300x300 from backend (vector SVG, scales perfectly)
  // so we get consistent high-quality rendering at any display size.

  const promise = apiCall("mol-svg", body, 10000)
    .then((r) => r.json())
    .then((d) => {
      const svg: string | null = d.svg || null;
      if (svg) _svgCache.set(cacheKey, svg);
      _pendingCache.delete(cacheKey);
      return svg;
    })
    .catch(() => {
      _pendingCache.delete(cacheKey);
      return null;
    });
  _pendingCache.set(cacheKey, promise);
  return promise;
}

export default function MoleculeSVG({
  smiles,
  size = 200,
  className = "",
  showLabel = false,
  zoomable = false,
  label,
  idx,
  width,
  height,
}: {
  smiles: string;
  /** Default square size (used when width/height not specified) */
  size?: number;
  className?: string;
  showLabel?: boolean;
  zoomable?: boolean;
  /** Optional atom label annotation (R-group style, e.g. "Cys145A") */
  label?: string;
  /** Atom index to apply label (default: 0) */
  idx?: number;
  /** Container display width in pixels (overrides size) */
  width?: number;
  /** Container display height in pixels (overrides size) */
  height?: number;
}) {
  const [svg, setSvg] = useState<string | null>(null);
  const [error, setError] = useState(false);
  const [zoomOpen, setZoomOpen] = useState(false);
  const mountedRef = useRef(true);

  useEffect(() => {
    if (!smiles) return;
    mountedRef.current = true;

    // Check cache synchronously first
    const cacheKey = buildCacheKey(smiles, label, idx);
    if (_svgCache.has(cacheKey)) {
      React.startTransition(() => setSvg(_svgCache.get(cacheKey)!));
      return;
    }

    fetchSvg(smiles, label, idx).then((s) => {
      if (mountedRef.current) {
        if (s) setSvg(s);
        else setError(true);
      }
    });

    return () => { mountedRef.current = false; };
  }, [smiles, label, idx]);

  // Container dimensions: prefer explicit width/height, fallback to square size
  const containerStyle = {
    width: width || size,
    height: height || size,
  };

  if (error) {
    return (
      <div
        className={`flex items-center justify-center bg-gray-50 rounded-lg ${className}`}
        style={containerStyle}
      >
        <span className="text-xs text-gray-400 text-center px-2">{smiles}</span>
      </div>
    );
  }

  if (!svg) {
    return (
      <div
        className={`flex items-center justify-center bg-gray-50 rounded-lg ${className}`}
        style={containerStyle}
      >
        <Loader2 className="w-4 h-4 text-gray-300 animate-spin" />
      </div>
    );
  }

  return (
    <>
      <div
        className={`relative overflow-hidden rounded-lg bg-white ${zoomable ? "cursor-pointer group" : ""} ${className}`}
        style={containerStyle}
        onClick={zoomable ? () => setZoomOpen(true) : undefined}
      >
        <div
          dangerouslySetInnerHTML={{ __html: svg }}
          className="w-full h-full [&>svg]:w-full [&>svg]:h-full [&>svg]:max-w-full [&>svg]:max-h-full"
        />
        {zoomable && (
          <div className="absolute top-1 right-1 opacity-0 group-hover:opacity-100 transition-opacity bg-white/80 rounded-md p-0.5 shadow-sm">
            <ZoomIn className="w-3 h-3 text-gray-500" />
          </div>
        )}
        {showLabel && (
          <div className="absolute bottom-0 left-0 right-0 bg-white/90 backdrop-blur-sm border-t border-gray-100 px-1.5 py-0.5 rounded-b-lg">
            <code className="text-xs font-mono text-gray-500 block truncate">{smiles}</code>
          </div>
        )}
      </div>
      <ZoomModal open={zoomOpen} onClose={() => setZoomOpen(false)} svgHtml={svg} smiles={smiles} title={label || smiles} label={label} idx={idx} />
    </>
  );
}
