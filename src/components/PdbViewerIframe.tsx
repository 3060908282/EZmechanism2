'use client';

import React, { useState, useEffect, useRef } from 'react';
import { Loader2, AlertCircle, ExternalLink } from 'lucide-react';

interface PdbViewerIframeProps {
  pdbId: string;
  height?: string;
  className?: string;
}

/**
 * Embedded RCSB PDB 3D viewer via iframe.
 * Uses RCSB's own WebGL-enabled Mol* viewer running on their servers,
 * providing professional PyMOL-quality cartoon rendering that matches
 * exactly what users see on rcsb.org.
 *
 * Requirements:
 * - Valid 4-character PDB ID (e.g., "1EJG", "3NIR")
 * - Internet access (loads from rcsb.org)
 *
 * Limitations:
 * - No custom residue highlighting (use PdbViewer3D for that)
 * - No chain selection/dimming
 * - No residue click callbacks
 * - Requires valid PDB ID in the RCSB database
 */
export default function PdbViewerIframe({
  pdbId,
  height = '650px',
  className = '',
}: PdbViewerIframeProps) {
  const [loaded, setLoaded] = useState(false);
  const [error, setError] = useState(false);
  const iframeRef = useRef<HTMLIFrameElement>(null);

  // Validate PDB ID format
  const isValidId = /^[A-Za-z0-9]{4}$/.test(pdbId.trim());

  // RCSB embed URL - this loads the Mol* viewer with the specified structure
  const embedUrl = isValidId
    ? `https://www.rcsb.org/structure/${pdbId.trim().toUpperCase()}/embedded`
    : '';

  const [prevPdbId, setPrevPdbId] = useState(pdbId);
  if (pdbId !== prevPdbId) {
    setPrevPdbId(pdbId);
    setLoaded(false);
    setError(false);
  }

  const handleLoad = () => {
    // Give the viewer a moment to initialize WebGL
    setTimeout(() => setLoaded(true), 800);
  };

  const handleError = () => {
    setError(true);
  };

  if (!isValidId) {
    return (
      <div
        className={`flex items-center justify-center bg-gray-50 rounded-md border border-gray-200 ${className}`}
        style={{ height }}
      >
        <div className="text-center text-gray-400">
          <AlertCircle className="w-8 h-8 mx-auto mb-2 text-gray-300" />
          <p className="text-sm">Invalid PDB ID: {pdbId}</p>
          <p className="text-xs text-gray-300 mt-1">Requires a valid 4-character PDB identifier</p>
        </div>
      </div>
    );
  }

  return (
    <div className={`relative ${className}`} style={{ height }}>
      {/* Loading overlay */}
      {!loaded && !error && (
        <div
          className="absolute inset-0 flex items-center justify-center bg-white rounded-md border border-gray-200 z-10"
          style={{ height }}
        >
          <div className="text-center">
            <Loader2 className="w-8 h-8 animate-spin text-emerald-500 mx-auto mb-3" />
            <p className="text-sm font-medium text-gray-600">Loading RCSB 3D Viewer</p>
            <p className="text-xs text-gray-400 mt-1">Initializing Mol* WebGL engine...</p>
            <p className="text-[10px] text-gray-300 mt-2 font-mono">{pdbId.trim().toUpperCase()}</p>
          </div>
        </div>
      )}

      {/* Error state */}
      {error && (
        <div
          className="absolute inset-0 flex items-center justify-center bg-gray-50 rounded-md border border-gray-200 z-10"
          style={{ height }}
        >
          <div className="text-center text-gray-500 px-6">
            <AlertCircle className="w-8 h-8 mx-auto mb-2 text-red-300" />
            <p className="text-sm font-medium">Failed to load RCSB viewer</p>
            <p className="text-xs text-gray-400 mt-1">The RCSB server may be temporarily unavailable</p>
            <a
              href={`https://www.rcsb.org/structure/${pdbId.trim().toUpperCase()}`}
              target="_blank"
              rel="noopener noreferrer"
              className="inline-flex items-center gap-1.5 mt-3 text-xs text-emerald-600 hover:text-emerald-700 hover:underline"
            >
              Open on RCSB website <ExternalLink className="w-3 h-3" />
            </a>
          </div>
        </div>
      )}

      {/* RCSB embedded viewer iframe */}
      <iframe
        ref={iframeRef}
        src={embedUrl}
        onLoad={handleLoad}
        onError={handleError}
        title={`RCSB PDB 3D Viewer - ${pdbId}`}
        className="w-full h-full rounded-md border border-gray-200"
        style={{
          minHeight: height,
          border: 'none',
        }}
        allow="fullscreen; clipboard-write"
        sandbox="allow-scripts allow-same-origin allow-popups allow-forms allow-top-navigation-by-user-activation"
      />
    </div>
  );
}
