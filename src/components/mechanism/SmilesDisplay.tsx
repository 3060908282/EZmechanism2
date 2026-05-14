"use client";

import React, { useState } from "react";
import { ChevronDown, ChevronUp } from "lucide-react";
import CopyButton from "@/components/mechanism/CopyButton";

export default function SmilesDisplay({
  smiles,
  maxLen = 40,
  className = "",
}: {
  smiles: string;
  maxLen?: number;
  className?: string;
}) {
  const [expanded, setExpanded] = useState(false);
  const truncated = smiles.length > maxLen && !expanded;
  const display = truncated ? smiles.slice(0, maxLen) + "..." : smiles;

  return (
    <div className={`flex items-center gap-1.5 ${className}`}>
      <code className="text-sm font-mono text-gray-800 break-all">{display}</code>
      {smiles.length > maxLen && (
        <button
          onClick={() => setExpanded(!expanded)}
          className="text-gray-400 hover:text-emerald-600 shrink-0"
        >
          {expanded ? <ChevronUp className="w-3 h-3" /> : <ChevronDown className="w-3 h-3" />}
        </button>
      )}
      <CopyButton text={smiles} />
    </div>
  );
}
