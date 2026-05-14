"use client";

import { Copy } from "lucide-react";
import { useToast } from "@/hooks/use-toast";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";

export default function CopyButton({ text }: { text: string }) {
  const { toast } = useToast();
  return (
    <TooltipProvider>
      <Tooltip>
        <TooltipTrigger asChild>
          <button
            onClick={() => {
              navigator.clipboard.writeText(text);
              toast({ title: "Copied", description: "SMILES copied to clipboard" });
            }}
            className="text-gray-400 hover:text-emerald-600 transition-colors shrink-0"
          >
            <Copy className="w-3 h-3" />
          </button>
        </TooltipTrigger>
        <TooltipContent>Copy SMILES</TooltipContent>
      </Tooltip>
    </TooltipProvider>
  );
}
