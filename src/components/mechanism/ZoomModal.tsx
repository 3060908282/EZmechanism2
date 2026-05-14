"use client";

import React, { useState, useEffect } from "react";
import { Loader2 } from "lucide-react";
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from "@/components/ui/dialog";
import { apiCall } from "@/lib/api";

export default function ZoomModal({ 
  open, 
  onClose, 
  svgHtml, 
  smiles, 
  title,
  label,
  idx,
}: { 
  open: boolean; 
  onClose: () => void; 
  svgHtml?: string; 
  smiles?: string; 
  title?: string;
  /** Optional atom label (R-group style) */
  label?: string;
  /** Atom index to apply label */
  idx?: number;
}) {
  const [fetchedSvg, setFetchedSvg] = useState<string | null>(null);
  const [loading, setLoading] = useState(false);

  // When smiles is provided, fetch the SVG from the API (with optional label)
  useEffect(() => {
    if (!open || !smiles) {
      React.startTransition(() => setFetchedSvg(null));
      return;
    }
    let cancelled = false;
    React.startTransition(() => setLoading(true));
    const body: Record<string, unknown> = { smiles, width: 500, height: 400 };
    if (label) body.label = label;
    if (idx !== undefined) body.idx = idx;
    apiCall("mol-svg", body, 10000)
      .then((r) => r.json())
      .then((d) => {
        if (!cancelled && d.svg) setFetchedSvg(d.svg);
      })
      .catch(() => {})
      .finally(() => {
        if (!cancelled) setLoading(false);
      });
    return () => { cancelled = true; };
  }, [open, smiles, label, idx]);

  const displaySvg = svgHtml || fetchedSvg;

  return (
    <Dialog open={open} onOpenChange={(o) => !o && onClose()}>
      <DialogContent className="max-w-3xl w-full p-2">
        <DialogHeader className="sr-only">
          <DialogTitle>{title || "Molecule Structure"}</DialogTitle>
        </DialogHeader>
        {loading ? (
          <div className="bg-white rounded-lg overflow-hidden flex items-center justify-center min-h-[400px]">
            <div className="text-center">
              <Loader2 className="w-6 h-6 text-violet-400 animate-spin mx-auto mb-2" />
              <p className="text-xs text-gray-400">Loading structure...</p>
            </div>
          </div>
        ) : displaySvg ? (
          <div className="bg-white rounded-lg overflow-hidden flex items-center justify-center min-h-[400px]">
            <div
              dangerouslySetInnerHTML={{ __html: displaySvg }}
              className="[&>svg]:w-full [&>svg]:h-auto [&>svg]:max-h-[70vh] [&>svg]:max-w-full"
            />
          </div>
        ) : (
          <div className="bg-gray-50 rounded-lg overflow-hidden flex items-center justify-center min-h-[400px]">
            <p className="text-sm text-gray-400">Unable to render structure</p>
          </div>
        )}
        {title && (
          <div className="text-center -mt-1">
            <span className="text-sm font-mono font-semibold text-violet-700">{title}</span>
          </div>
        )}
      </DialogContent>
    </Dialog>
  );
}
