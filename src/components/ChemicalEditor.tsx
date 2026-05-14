"use client";

import React, { useEffect, useRef, useImperativeHandle, forwardRef, useState } from "react";

export interface ChemicalEditorHandle {
  getSmiles: () => Promise<string>;
  setMolecule: (smiles: string) => Promise<boolean>;
  getMolfile: () => Promise<string>;
  getRxn: () => Promise<string>;
  automap: () => Promise<string>;
  clear: () => Promise<void>;
  generateImage: (width?: number, height?: number) => Promise<string>;
  isReady: () => boolean;
  waitForReady: (timeout?: number) => Promise<boolean>;
}

interface ChemicalEditorProps {
  height?: number;
  className?: string;
  onReady?: () => void;
  onChange?: (smiles: string) => void;
}

let msgIdCounter = 0;
function nextMsgId() {
  return `msg-${Date.now()}-${++msgIdCounter}`;
}

function waitForResponse(type: string, id: string, timeout = 10000): Promise<{ data?: string; error?: string; success?: boolean }> {
  return new Promise((resolve) => {
    const handler = (event: MessageEvent) => {
      const d = event.data;
      if (d && d.type === type && d.id === id) {
        window.removeEventListener("message", handler);
        resolve(d);
      }
    };
    window.addEventListener("message", handler);
    const timer = setTimeout(() => {
      window.removeEventListener("message", handler);
      resolve({ data: "", error: "timeout" });
    }, timeout);
  });
}

// Ketcher served as static files from Next.js public/ketcher/
// Must use explicit /ketcher/index.html because Next.js dev mode
// intercepts /ketcher/ and /ketcher as routes (308 redirect → 404)
const EDITOR_URL = "/ketcher/index.html";

const ChemicalEditor = forwardRef<ChemicalEditorHandle, ChemicalEditorProps>(
  function ChemicalEditor({ height = 420, className = "", onReady, onChange }, ref) {
    const iframeRef = useRef<HTMLIFrameElement>(null);
    const readyRef = useRef(false);
    const readyPromiseRef = useRef<(() => void)[]>([]);
    const [loading, setLoading] = useState(true);
    const [error, setError] = useState<string | null>(null);
    const changeTimerRef = useRef<ReturnType<typeof setTimeout> | null>(null);
    const lastSmilesRef = useRef<string>("");
    const mountedRef = useRef(true);

    // Wait for the editor to become ready (resolves existing promises)
    const notifyReady = useRef(() => {
      readyRef.current = true;
      setLoading(false);
      setError(null);
      // Resolve all pending waitForReady promises
      readyPromiseRef.current.forEach((resolve) => resolve());
      readyPromiseRef.current = [];
      onReady?.();
    });

    useImperativeHandle(ref, () => ({
      getSmiles: async () => {
        if (!readyRef.current) return "";
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "getSmiles", id }, "*");
          const res = await waitForResponse("getSmilesResult", id, 5000);
          if (res.error) return "";
          const smiles = res.data || "";
          lastSmilesRef.current = smiles;
          return smiles;
        } catch {
          return "";
        }
      },
      setMolecule: async (smiles: string): Promise<boolean> => {
        // Wait for editor to be ready with a timeout
        if (!readyRef.current) {
          const ready = await new Promise<boolean>((resolve) => {
            const timer = setTimeout(() => resolve(false), 15000);
            readyPromiseRef.current.push(() => {
              clearTimeout(timer);
              resolve(true);
            });
          });
          if (!ready || !mountedRef.current) return false;
        }
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "setMolecule", id, data: smiles }, "*");
          const res = await waitForResponse("setMoleculeResult", id, 10000);
          if (!res.error) {
            // After setMolecule, immediately fetch Ketcher's canonical SMILES
            // to prevent spurious onChange from polling detecting a SMILES diff
            try {
              const getId = nextMsgId();
              iframeRef.current?.contentWindow?.postMessage({ type: "getSmiles", id: getId }, "*");
              const getRes = await waitForResponse("getSmilesResult", getId, 5000);
              if (!getRes.error && getRes.data) {
                lastSmilesRef.current = getRes.data;
              } else {
                lastSmilesRef.current = smiles;
              }
            } catch {
              lastSmilesRef.current = smiles;
            }
            return true;
          }
          return false;
        } catch {
          return false;
        }
      },
      getMolfile: async () => {
        if (!readyRef.current) return "";
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "getMolfile", id }, "*");
          const res = await waitForResponse("getMolfileResult", id, 5000);
          return res.error ? "" : (res.data || "");
        } catch {
          return "";
        }
      },
      getRxn: async () => {
        if (!readyRef.current) return "";
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "getRxn", id }, "*");
          const res = await waitForResponse("getRxnResult", id, 10000);
          return res.error ? "" : (res.data || "");
        } catch {
          return "";
        }
      },
      automap: async () => {
        if (!readyRef.current) return "";
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "automap", id }, "*");
          const res = await waitForResponse("automapResult", id, 30000);
          return res.error ? "" : (res.data || "");
        } catch {
          return "";
        }
      },
      clear: async () => {
        if (!readyRef.current) return;
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "clear", id }, "*");
          await waitForResponse("clearResult", id, 5000);
          lastSmilesRef.current = "";
        } catch {
          // ignore
        }
      },
      generateImage: async (w = 500, h = 400) => {
        if (!readyRef.current) return "";
        try {
          const id = nextMsgId();
          iframeRef.current?.contentWindow?.postMessage({ type: "generateImage", id, data: { width: w, height: h } }, "*");
          const res = await waitForResponse("generateImageResult", id, 15000);
          return res.error ? "" : (res.data || "");
        } catch {
          return "";
        }
      },
      isReady: () => readyRef.current,
      waitForReady: async (timeout = 15000) => {
        if (readyRef.current) return true;
        return new Promise<boolean>((resolve) => {
          const timer = setTimeout(() => resolve(false), timeout);
          readyPromiseRef.current.push(() => {
            clearTimeout(timer);
            resolve(true);
          });
        });
      },
    }));

    useEffect(() => {
      mountedRef.current = true;
      const handler = (event: MessageEvent) => {
        const d = event.data;
        if (d && d.type === "ketcherReady") {
          notifyReady.current();
        }
        if (d && d.type === "ketcherError") {
          setError(d.error || "Failed to load editor");
          setLoading(false);
        }
      };
      window.addEventListener("message", handler);
      return () => {
        mountedRef.current = false;
        window.removeEventListener("message", handler);
      };
    }, [onReady]);

    // Auto-detect changes by polling
    useEffect(() => {
      if (!onChange) return;
      let polling = false;
      const interval = setInterval(async () => {
        if (!readyRef.current || !iframeRef.current || polling) return;
        polling = true;
        try {
          const id = nextMsgId();
          iframeRef.current.contentWindow?.postMessage({ type: "getSmiles", id }, "*");
          const res = await waitForResponse("getSmilesResult", id, 5000);
          const newSmiles = (res.data && !res.error) ? (res.data || "") : lastSmilesRef.current;
          if (newSmiles !== lastSmilesRef.current) {
            if (changeTimerRef.current) clearTimeout(changeTimerRef.current);
            changeTimerRef.current = setTimeout(() => {
              lastSmilesRef.current = newSmiles;
              onChange(newSmiles);
            }, 300);
          }
        } catch {
          // ignore polling errors
        } finally {
          polling = false;
        }
      }, 2000);
      return () => {
        clearInterval(interval);
        if (changeTimerRef.current) clearTimeout(changeTimerRef.current);
      };
    }, [onChange]);

    return (
      <div className={`relative rounded-lg overflow-hidden border border-gray-200 bg-white ${className}`}>
        <iframe
          ref={iframeRef}
          src={EDITOR_URL}
          style={{ width: "100%", height, border: "none" }}
          title="Chemical Structure Editor"
        />
        {loading && (
          <div className="absolute inset-0 flex items-center justify-center bg-gray-50/90 z-10">
            <div className="text-center">
              <div className="w-8 h-8 border-3 border-gray-200 border-t-teal-500 rounded-full animate-spin mx-auto mb-2" />
              <p className="text-xs text-gray-500">Loading chemical editor...</p>
            </div>
          </div>
        )}
        {error && (
          <div className="absolute inset-0 flex items-center justify-center bg-gray-50/95 z-10">
            <div className="text-center p-4">
              <p className="text-sm text-red-600 font-medium mb-1">Editor Error</p>
              <p className="text-xs text-gray-500">{error}</p>
            </div>
          </div>
        )}
      </div>
    );
  }
);

export default ChemicalEditor;
