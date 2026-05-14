import { API_BASE } from "@/lib/constants";

export async function apiCall(action: string, body: Record<string, unknown>, timeoutMs = 120000): Promise<Response> {
  return fetch(API_BASE, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ action, ...body }),
    signal: AbortSignal.timeout(timeoutMs),
  });
}

/**
 * Safely read a Response body as JSON.
 * Handles: HTML error pages, non-200 status codes, empty responses, network errors.
 * NEVER throws SyntaxError — always throws a readable Error message.
 */
async function safeParseJson(r: Response, context = "API"): Promise<Record<string, unknown>> {
  const statusMsg = `${r.status} ${r.statusText}`;

  // Read raw text first — never call r.json() directly
  let text: string;
  try {
    text = await r.text();
  } catch {
    throw new Error(`Network error reading ${context} response (connection lost or timeout)`);
  }

  if (!text.trim()) {
    throw new Error(`Empty response from ${context} (${statusMsg})`);
  }

  // Detect HTML response (even if status is 200, proxies can inject HTML)
  const trimmed = text.trimStart();
  if (trimmed.startsWith("<!DOCTYPE") || trimmed.startsWith("<html") || trimmed.startsWith("<HTML")) {
    // Extract visible text from HTML for a readable error message
    const bodyMatch = text.match(/<body[^>]*>([\s\S]*?)<\/body>/i) || text.match(/<h[1-6][^>]*>(.*?)<\/h[1-6]>/i);
    const visibleText = bodyMatch
      ? bodyMatch[1].replace(/<[^>]*>/g, " ").replace(/\s+/g, " ").trim().slice(0, 150)
      : text.replace(/<[^>]*>/g, " ").replace(/\s+/g, " ").trim().slice(0, 150);
    throw new Error(
      `${context} returned an HTML error page (${statusMsg}). ` +
      `This usually means the server or proxy timed out. ` +
      (visibleText ? `Message: "${visibleText}"` : "The request may have taken too long — try reducing max_configs.")
    );
  }

  // Try to parse as JSON
  try {
    return JSON.parse(text);
  } catch {
    // Not valid JSON — show a snippet
    const snippet = text.slice(0, 120).replace(/\s+/g, " ");
    throw new Error(
      `${context} returned invalid JSON (${statusMsg}). ` +
      `Response preview: "${snippet}..."`
    );
  }
}

/**
 * Like apiCall() but also validates response and parses JSON safely.
 * Returns the parsed JSON body on success.
 * NEVER throws SyntaxError — always throws a readable Error message.
 */
export async function apiCallJson<T = Record<string, unknown>>(
  action: string,
  body: Record<string, unknown>,
  timeoutMs = 120000,
): Promise<T> {
  let r: Response;
  try {
    r = await apiCall(action, body, timeoutMs);
  } catch (err) {
    if (err instanceof DOMException && err.name === "TimeoutError") {
      throw new Error(
        `Request to ${action} timed out after ${Math.round(timeoutMs / 1000)}s. ` +
        "The server is busy — try again with fewer rules or configs."
      );
    }
    throw new Error(`Network error calling ${action}: ${err instanceof Error ? err.message : String(err)}`);
  }

  // Check HTTP status
  if (!r.ok) {
    const parsed = await safeParseJson(r, action).catch(() => null);
    if (parsed && typeof parsed === "object" && parsed.error) {
      let msg = String(parsed.error);
      if (parsed.hint) msg += `. ${parsed.hint}`;
      throw new Error(msg);
    }
    // Couldn't parse JSON error — check for HTML
    const contentType = r.headers.get("content-type") || "";
    if (contentType.includes("text/html")) {
      throw new Error(
        `Server error ${r.status} — ${action} request failed. ` +
        "The server may be overloaded or the request timed out."
      );
    }
    throw new Error(`Server error ${r.status} ${r.statusText} for action "${action}"`);
  }

  // Status is OK — but still safely parse (proxy might have injected HTML)
  const data = await safeParseJson(r, action);
  return data as T;
}

let _idCounter = 0;
export function genId() {
  return `mol-${Date.now()}-${++_idCounter}`;
}

/**
 * Normalize a mechanism search result.
 *
 * Fixes legacy/incomplete data by:
 * 1. Adding path_id / total_steps if missing
 * 2. Adding step_num / direction if missing
 * 3. Filtering out invalid molecule strings (<rdkit garbage)
 * 4. Providing empty-array defaults for from_config / to_config
 */
export function normalizeSearchResult(data: Record<string, unknown>): Record<string, unknown> {
  const result = { ...data };

  // Normalize paths
  const paths = (result.paths as Array<Record<string, unknown>> | undefined) || [];
  const graph = result.graph as Record<string, unknown> | undefined;
  const graphNodes = (graph?.nodes as Array<Record<string, unknown>> | undefined) || [];

  // Build a lookup: node id → molecules array
  const nodeMolMap = new Map<string, string[]>();
  for (const node of graphNodes) {
    const id = String(node.id ?? "");
    const mols = (node.molecules as string[]) || [];
    if (id && mols.length > 0) {
      nodeMolMap.set(id, mols);
    }
  }

  const normalizedPaths = paths.map((path, pIdx) => {
    const nPath = { ...path };
    if (nPath.path_id == null) nPath.path_id = pIdx;
    if (nPath.total_steps == null && nPath.num_steps != null) nPath.total_steps = nPath.num_steps;

    const steps = (nPath.steps as Array<Record<string, unknown>> | undefined) || [];
    nPath.steps = steps.map((step, sIdx) => {
      const nStep = { ...step };
      if (nStep.step_num == null) nStep.step_num = sIdx + 1;
      if (!nStep.direction) nStep.direction = "forward";

      // Sanitize from_config / to_config molecules
      for (const configKey of ["from_config", "to_config"] as const) {
        const config = nStep[configKey] as Record<string, unknown> | undefined;
        if (config && Array.isArray(config.molecules)) {
          // Filter out invalid strings (Mol object reprs, empty strings)
          config.molecules = config.molecules
            .filter((m: unknown) =>
              typeof m === "string" &&
              m.length > 0 &&
              !m.startsWith("<rdkit") &&
              !m.startsWith("__WATER") &&
              !m.startsWith("__PROTON")
            );
        } else {
          // Provide empty default
          nStep[configKey] = { molecules: [], source: "unknown" };
        }
      }

      return nStep;
    });

    return nPath;
  });

  result.paths = normalizedPaths;

  // Ensure graph edges have direction
  const graphEdges = (graph?.edges as Array<Record<string, unknown>> | undefined) || [];
  for (const edge of graphEdges) {
    if (!edge.direction) edge.direction = "forward";
    if (!edge.label && edge.rule_name) edge.label = edge.rule_name;
  }

  return result;
}
