"use client";

import { FlaskConical } from "lucide-react";
import { Card, CardContent } from "@/components/ui/card";

export default function EmptyState({ mode }: { mode?: string }) {
  return (
    <Card className="border-emerald-100 border-dashed shadow-none bg-white/50">
      <CardContent className="flex flex-col items-center justify-center py-14 text-center">
        <div className="w-14 h-14 rounded-2xl bg-gradient-to-br from-emerald-100 to-teal-100 flex items-center justify-center mb-4">
          <FlaskConical className="w-7 h-7 text-emerald-400" />
        </div>
        <h3 className="text-base font-semibold text-gray-600 mb-1">
          No Prediction Yet
        </h3>
        <p className="text-xs text-gray-400 max-w-sm leading-relaxed">
          {mode === "pdb"
            ? "Complete the EzMechanism workflow steps to predict a mechanism from PDB structure data."
            : "Enter a substrate SMILES and click \"Predict\" to discover matching reaction rules and predicted pathways."}
        </p>
      </CardContent>
    </Card>
  );
}
