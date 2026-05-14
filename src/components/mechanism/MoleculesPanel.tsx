"use client";

import React, { useState } from "react";
import { Layers, Loader2, Link2, ExternalLink } from "lucide-react";
import { Card, CardContent } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Popover, PopoverContent, PopoverTrigger } from "@/components/ui/popover";
import { apiCall } from "@/lib/api";
import { MoleculeLinksResponse, DbLink, MechanismSearchResult } from "@/lib/types";
import MoleculeSVG from "@/components/mechanism/MoleculeSVG";

export default function MoleculesPanel({
  searchResult,
  onExploreProduct,
}: {
  searchResult: MechanismSearchResult | null;
  onExploreProduct: (smiles: string) => void;
}) {
  const [molLinksCache, setMolLinksCache] = useState<Record<string, MoleculeLinksResponse>>({});
  const [loadingLinks, setLoadingLinks] = useState<string | null>(null);
  const [linksOpen, setLinksOpen] = useState<string | null>(null);

  const handleFetchMolLinks = async (smi: string) => {
    if (molLinksCache[smi]) { setLinksOpen(linksOpen === smi ? null : smi); return; }
    setLoadingLinks(smi);
    setLinksOpen(smi);
    try {
      const r = await apiCall("molecule-links", { smiles: smi }, 10000);
      if (!r.ok) throw new Error("Failed");
      const data = await r.json();
      setMolLinksCache((prev) => ({ ...prev, [smi]: data }));
    } catch { /* ignore */ } finally {
      setLoadingLinks(null);
    }
  };

  // Extract unique molecules from search result graph nodes
  const allProducts = new Set<string>();
  if (searchResult?.graph?.nodes) {
    for (const node of searchResult.graph.nodes) {
      if (node.molecules) {
        for (const mol of node.molecules) {
          if (mol && !mol.startsWith("__WATER") && !mol.startsWith("__PROTON") && !mol.startsWith("<rdkit")) {
            allProducts.add(mol);
          }
        }
      }
    }
  }
  const products = Array.from(allProducts);

  if (products.length === 0) {
    return (
      <div className="text-center py-8">
        <Layers className="w-7 h-7 text-gray-300 mx-auto mb-2" />
        <p className="text-gray-500 text-sm">No products predicted yet.</p>
      </div>
    );
  }

  return (
    <div className="space-y-4">
      <div className="flex items-center justify-between">
        <p className="text-xs text-gray-500">
          <span className="font-semibold text-emerald-700">{products.length}</span> unique product(s) predicted
        </p>
      </div>
      <div className="grid grid-cols-2 sm:grid-cols-3 gap-3">
        {products.map((smi, idx) => {
          const molLinks = molLinksCache[smi];
          return (
            <Card
              key={`${smi}-${idx}`}
              className="border-gray-100 hover:border-emerald-200 hover:shadow-sm transition-all cursor-pointer group relative"
              onClick={() => onExploreProduct(smi)}
            >
              <CardContent className="p-3">
                <div className="absolute top-1.5 right-1.5 z-10 opacity-0 group-hover:opacity-100 transition-opacity">
                  <Popover open={linksOpen === smi} onOpenChange={(open) => { if (!open) setLinksOpen(null); }}>
                    <PopoverTrigger asChild>
                      <button
                        onClick={(e) => { e.stopPropagation(); handleFetchMolLinks(smi); }}
                        className="p-1 rounded-md bg-white shadow-sm border border-gray-100 hover:bg-blue-50 hover:border-blue-200 transition-colors text-gray-400 hover:text-blue-600"
                      >
                        {loadingLinks === smi ? <Loader2 className="w-3 h-3 animate-spin" /> : <Link2 className="w-3 h-3" />}
                      </button>
                    </PopoverTrigger>
                    <PopoverContent className="w-48 p-2" side="right" align="start" onClick={(e) => e.stopPropagation()}>
                      <div className="space-y-0.5">
                        <div className="text-xs font-medium text-gray-500 uppercase tracking-wider mb-1">Search In</div>
                        {molLinks && Object.entries(molLinks).map(([dbKey, dbVal]) => {
                          if (!dbVal) return null;
                          const link = dbVal as DbLink;
                          if (!link.available || !link.url) return null;
                          return (
                            <a key={dbKey} href={link.url} target="_blank" rel="noopener noreferrer" className="flex items-center gap-1.5 text-xs text-blue-600 hover:text-blue-800 hover:underline py-0.5">
                              <ExternalLink className="w-2.5 h-2.5 shrink-0" />
                              {link.label || dbKey}
                            </a>
                          );
                        })}
                        {loadingLinks === smi && !molLinks && <span className="text-xs text-gray-400">Loading...</span>}
                      </div>
                    </PopoverContent>
                  </Popover>
                </div>
                <div className="flex justify-center mb-2">
                  <MoleculeSVG smiles={smi} size={130} showLabel />
                </div>
                <div className="flex justify-center">
                  <Button variant="ghost" size="sm" className="h-6 text-xs text-emerald-600 hover:text-emerald-700 opacity-0 group-hover:opacity-100 transition-opacity">
                    <ExternalLink className="w-3 h-3 mr-1" />
                    Explore
                  </Button>
                </div>
              </CardContent>
            </Card>
          );
        })}
      </div>
    </div>
  );
}
