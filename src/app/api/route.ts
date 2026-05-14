import { NextRequest, NextResponse } from "next/server";

/**
 * Proxy API route that executes Python prediction scripts on-demand.
 *
 * Architecture:
 * - Heavy requests (mechanism-search) use ASYNC mode: mechanism-search-start
 *   spawns a background Python process and returns a jobId immediately.
 *   Frontend polls mechanism-search-status for progress + results.
 * - The old sync mechanism-search action is kept as fallback.
 * - Lightweight requests (mol-info, mol-svg) always use child_process.
 * - All user inputs are passed as CLI arguments (not string interpolation)
 *   to prevent code injection attacks.
 */

import path from "path";
import { randomUUID } from "crypto";

// Use relative paths so it works on any OS (Linux, Windows, Mac)
const PREDICTION_DIR = path.join(process.cwd(), "mini-services", "prediction-service");
const FLASK_BASE = `http://localhost:3003`;
// Use conda Python with RDKit on Linux sandbox, system python on Windows
const PYTHON_BIN = process.platform === "win32"
  ? "python"
  : "/home/z/.venv/bin/python3";

// Job directory for async mechanism search
const JOB_BASE_DIR = path.join(process.cwd(), ".mech-search-jobs");
const JOB_TTL_MS = 30 * 60 * 1000; // 30 minutes

async function runPython(args: string[], timeoutMs: number = 120000, stdinData?: string): Promise<string> {
  const { execFile, spawn } = await import("child_process");
  const { promisify } = await import("util");

  // If stdinData is provided, use spawn + pipe (execFile 'input' option is broken in Node v24)
  if (stdinData) {
    return new Promise<string>((resolve, reject) => {
      const child = spawn(PYTHON_BIN, args, {
        cwd: PREDICTION_DIR,
        stdio: ["pipe", "pipe", "pipe"],
      });
      const timer = setTimeout(() => { child.kill("SIGKILL"); reject(new Error("Python script timed out")); }, timeoutMs);

      let stdout = "";
      let stderr = "";
      child.stdout.on("data", (d: Buffer) => { stdout += d.toString(); });
      child.stderr.on("data", (d: Buffer) => { stderr += d.toString(); });

      child.stdin.write(stdinData);
      child.stdin.end();

      child.on("close", (code) => {
        clearTimeout(timer);
        if (code === 0 || (code !== 0 && stdout.trim())) {
          resolve(stdout);
        } else {
          reject(new Error(stderr.trim() || `Python script exited with code ${code}`));
        }
      });
      child.on("error", (err) => { clearTimeout(timer); reject(err); });
    });
  }

  // No stdin needed — use execFile (fast, simple)
  const execFileAsync = promisify(execFile);
  try {
    const result = await execFileAsync(PYTHON_BIN, args, {
      cwd: PREDICTION_DIR,
      timeout: timeoutMs,
      maxBuffer: 50 * 1024 * 1024, // 50MB buffer for large rule sets
      encoding: "utf-8", // Always return strings, not Buffers
    });
    return result.stdout;
  } catch (error: unknown) {
    const err = error as { stdout?: string; stderr?: string; message?: string };
    if (err.stdout) return err.stdout;
    throw new Error(err.stderr?.trim() || err.message || "Python script execution failed");
  }
}

/** Try to call the Flask backend. Returns null if Flask is not available. */
async function tryFlask(
  action: string,
  body: Record<string, unknown>,
  timeoutMs: number = 15000,
): Promise<NextResponse<unknown> | null> {
  try {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), timeoutMs);
    const res = await fetch(`${FLASK_BASE}/api/${action}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(body),
      signal: controller.signal,
    });
    clearTimeout(timer);
    if (!res.ok) return null;
    const data = await res.json();
    return NextResponse.json(data);
  } catch {
    // Flask not running or request failed — fall back
    return null;
  }
}

// SMILES validation: only allow safe characters
const SMILES_REGEX = /^[A-Za-z0-9@+\-[\]()\\/{}%=#.:;!?*]+$/;

function validateSmiles(smiles: unknown): string | null {
  if (typeof smiles !== "string" || !smiles.trim()) return null;
  const s = smiles.trim();
  if (!SMILES_REGEX.test(s)) return null;
  if (s.length > 500) return null; // unreasonably long
  return s;
}

// SMARTS validation: slightly more permissive than SMILES
const SMARTS_REGEX = /^[A-Za-z0-9@+\-[\]()\\/{}%=#.:;!?*&$~]+$/;

function validateSmarts(smarts: unknown): string | null {
  if (typeof smarts !== "string" || !smarts.trim()) return null;
  const s = smarts.trim();
  if (!SMARTS_REGEX.test(s)) return null;
  if (s.length > 2000) return null;
  return s;
}

export async function POST(request: NextRequest) {
  let actionName = "predict";
  try {
    const body = await request.json();
    const action = body.action || "predict";
    actionName = action;

    if (action === "health") {
      return NextResponse.json({
        status: "ok",
        service: "M-CSA Prediction (direct Python)",
        note: "Rules loaded on-demand from rules_nat_met.xlsx",
      });
    }

    if (action === "mol-info") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }

      const flaskMi = await tryFlask("mol-info", { smiles }, 10000);
      if (flaskMi) return flaskMi;

      const result = await runPython([
        path.join(PREDICTION_DIR, "mol_info.py"), smiles,
      ], 10000);

      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "smiles-explicit-h") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }

      // Try Flask backend first
      const flaskResult = await tryFlask("smiles-explicit-h", { smiles }, 10000);
      if (flaskResult) return flaskResult;

      // Fallback: direct Python script
      const result = await runPython([
        path.join(PREDICTION_DIR, "mol_info.py"), smiles, "--explicit-h",
      ], 10000);

      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "mol-svg") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }

      const args: string[] = [path.join(PREDICTION_DIR, "mol_svg.py"), smiles];
      // Optional label to annotate a specific atom (R-group style)
      if (typeof body.label === "string" && body.label.trim()) {
        args.push("--label", body.label.trim());
      }
      // Atom index to apply label (default: 0)
      if (typeof body.idx === "number") {
        args.push("--idx", String(body.idx));
      }
      // Optional width/height
      if (typeof body.width === "number") {
        args.push("--width", String(body.width));
      }
      if (typeof body.height === "number") {
        args.push("--height", String(body.height));
      }

      const result = await runPython(args, 10000);

      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "mechanism-search-start") {
      // ─── Async search: validate, spawn Python in background, return jobId ───
      const reactants = validateSmiles(body.reactants);
      const products = validateSmiles(body.products);
      if (!reactants) {
        return NextResponse.json({ error: "Valid reactants SMILES are required" }, { status: 400 });
      }

      const rxnText = typeof body.rxn_text === "string" ? body.rxn_text.trim() : "";
      if (!rxnText && !products) {
        return NextResponse.json({ error: "Both valid reactants and products SMILES are required" }, { status: 400 });
      }

      const maxConfigs = Math.min(Math.max(Math.floor(Number(body.max_configs) || 300), 50), 2000);
      const maxBondLength = Math.min(Math.max(Number(body.max_bond_length) ?? 11, 0), 20);
      const reactantPdb = typeof body.reactant_pdb === "string" ? body.reactant_pdb.trim() : "";
      const productPdb = typeof body.product_pdb === "string" ? body.product_pdb.trim() : "";

      // Validate residues if provided
      let residuesJson = '';
      if (Array.isArray(body.residues)) {
        if (body.residues.length > 20) {
          return NextResponse.json({ error: "Maximum 20 residues allowed" }, { status: 400 });
        }
        for (const r of body.residues) {
          if (typeof r !== 'object' || !r.res_name || typeof r.res_num !== 'number') {
            return NextResponse.json({ error: "Each residue must have res_name (string) and res_num (number)" }, { status: 400 });
          }
        }
        residuesJson = JSON.stringify(body.residues);
      }

      const pdbId = typeof body.pdb_id === 'string' ? body.pdb_id.trim().toUpperCase() : '';
      const pdbText = typeof body.pdb_text === 'string' ? body.pdb_text.trim() : '';

      // Create job directory
      const fs = await import("fs");
      const os = await import("os");
      const jobId = randomUUID();
      const jobDir = path.join(JOB_BASE_DIR, jobId);
      fs.mkdirSync(jobDir, { recursive: true });

      const progressFile = path.join(jobDir, "progress.json");
      const resultFile = path.join(jobDir, "result.json");

      // Write initial progress
      fs.writeFileSync(progressFile, JSON.stringify({
        state: "STARTING",
        explored_nodes: 0,
        total_nodes: 0,
        current_iteration: 0,
        elapsed_seconds: 0,
        created_at: new Date().toISOString(),
        updated_at: new Date().toISOString(),
      }), "utf-8");

      // Build script arguments (same as sync mode, plus --progress-file and --result-file)
      const scriptArgs: string[] = [
        path.join(PREDICTION_DIR, "mechanism_search.py"),
        "--reactants", reactants,
        "--products", products || "",
        "--max-configs", String(maxConfigs),
        "--max-rules", "0",  // 0 = load all rules
        "--max-bond-length", String(maxBondLength),
        "--progress-file", progressFile,
        "--result-file", resultFile,
      ];

      if (reactantPdb) scriptArgs.push("--reactant-pdb", reactantPdb);
      if (productPdb) scriptArgs.push("--product-pdb", productPdb);
      if (residuesJson) scriptArgs.push("--residues", residuesJson);
      if (pdbId) scriptArgs.push("--pdb-id", pdbId);

      // M-CSA enzyme ID filter (Own Rules mode)
      const mcsaId = Math.floor(Number(body.mcsa_id) || 0);
      if (mcsaId > 0) scriptArgs.push("--mcsa-id", String(mcsaId));

      // Write PDB text to temp file
      let pdbTempFile: string | null = null;
      if (pdbText && pdbText.length > 100) {
        pdbTempFile = path.join(jobDir, "input.pdb");
        fs.writeFileSync(pdbTempFile, pdbText, "utf-8");
        scriptArgs.push("--pdb-text-file", pdbTempFile);
      }

      // Write RXN text to temp file
      let rxnTempFile: string | null = null;
      if (rxnText) {
        rxnTempFile = path.join(jobDir, "input.rxn");
        fs.writeFileSync(rxnTempFile, rxnText, "utf-8");
        scriptArgs.push("--rxn-text-file", rxnTempFile);
      }

      // Write ligand_mappings to temp file
      const ligandMappings = body.ligand_mappings;
      if (Array.isArray(ligandMappings) && ligandMappings.length > 0) {
        const validMappings = ligandMappings.filter(
          (m: any) => m && typeof m.smiles === 'string' && m.smiles.trim()
        );
        if (validMappings.length > 0) {
          const ligandFile = path.join(jobDir, "ligand_mappings.json");
          fs.writeFileSync(ligandFile, JSON.stringify(validMappings), "utf-8");
          scriptArgs.push("--ligand-mappings-file", ligandFile);
        }
      }

      // Spawn Python process in background (fire-and-forget)
      const { spawn } = await import("child_process");
      const child = spawn(PYTHON_BIN, scriptArgs, {
        cwd: PREDICTION_DIR,
        stdio: ["ignore", "pipe", "pipe"],
        detached: true,
      });

      // Log stderr for diagnostics but don't block
      let stderrLog = "";
      child.stderr?.on("data", (d: Buffer) => { stderrLog += d.toString(); });
      child.on("close", (code, signal) => {
        if (code !== 0 || signal) {  // Also catch signal-based kills (SIGKILL, OOM, etc.)
          // If process exited with error and progress file still says RUNNING/STARTING,
          // update it to ERROR
          try {
            const progressData = JSON.parse(fs.readFileSync(progressFile, "utf-8"));
            if (progressData.state === "RUNNING" || progressData.state === "STARTING") {
              fs.writeFileSync(progressFile, JSON.stringify({
                ...progressData,
                state: "ERROR",
                error: `Python process exited (code=${code}, signal=${signal}): ${stderrLog.slice(-500)}`,
                updated_at: new Date().toISOString(),
              }), "utf-8");
            }
          } catch { /* progress file may not exist or be corrupt */ }
        }
        // Clean up temp input files (PDB, RXN) — no longer needed after process exits
        for (const f of [pdbTempFile, rxnTempFile]) {
          if (f) { try { fs.unlinkSync(f); } catch { /* ignore */ } }
        }
      });
      child.unref(); // Don't wait for the child process

      return NextResponse.json({ jobId, state: "STARTING" });
    }

    if (action === "mechanism-search-status") {
      // ─── Read progress + result for a jobId ───
      const fs = await import("fs");
      const jobId = typeof body.jobId === "string" ? body.jobId.trim() : "";
      if (!jobId) {
        return NextResponse.json({ error: "jobId is required" }, { status: 400 });
      }

      // Validate jobId format (UUID) to prevent path traversal
      if (!/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i.test(jobId)) {
        return NextResponse.json({ error: "Invalid jobId format" }, { status: 400 });
      }

      const jobDir = path.join(JOB_BASE_DIR, jobId);
      const progressFile = path.join(jobDir, "progress.json");

      // TTL cleanup: remove old job directories (skip RUNNING ones)
      try {
        const entries = fs.readdirSync(JOB_BASE_DIR, { withFileTypes: true });
        const now = Date.now();
        for (const entry of entries) {
          if (!entry.isDirectory()) continue;
          try {
            const ep = path.join(JOB_BASE_DIR, entry.name, "progress.json");
            const stat = fs.statSync(ep);
            const pData = JSON.parse(fs.readFileSync(ep, "utf-8"));
            // Only clean up: (1) older than TTL AND (2) not RUNNING
            if (now - stat.mtimeMs > JOB_TTL_MS && pData.state !== "RUNNING" && pData.state !== "STARTING") {
              fs.rmSync(path.join(JOB_BASE_DIR, entry.name), { recursive: true, force: true });
            }
          } catch { /* skip invalid directories */ }
        }
      } catch { /* JOB_BASE_DIR may not exist yet */ }

      // Read progress file
      if (!fs.existsSync(progressFile)) {
        return NextResponse.json({ error: "Job not found", state: "UNKNOWN" }, { status: 404 });
      }

      try {
        const progressData = JSON.parse(fs.readFileSync(progressFile, "utf-8"));

        // Stale job detection: if RUNNING but not updated for 30+ minutes, mark as ERROR
        if (progressData.state === "RUNNING" || progressData.state === "STARTING") {
          const lastUpdate = new Date(progressData.updated_at).getTime();
          const staleThreshold = 30 * 60 * 1000; // 30 minutes
          if (Date.now() - lastUpdate > staleThreshold) {
            progressData.state = "ERROR";
            progressData.error = "Search process appears to have died (no progress update for 30+ minutes)";
            progressData.updated_at = new Date().toISOString();
            // Update the file too
            try {
              fs.writeFileSync(progressFile, JSON.stringify(progressData, null, 2), "utf-8");
            } catch { /* ignore write failure */ }
          }
        }

        // If search is DONE, also read the result file and include it
        if (progressData.state === "DONE" && progressData.result_file) {
          const resultPath = progressData.result_file;
          if (fs.existsSync(resultPath)) {
            try {
              const resultData = JSON.parse(fs.readFileSync(resultPath, "utf-8"));
              return NextResponse.json({ ...progressData, result: resultData });
            } catch {
              // Result file corrupt — return progress without result
              return NextResponse.json({ ...progressData, error: "Result file is corrupt" });
            }
          }
        }

        return NextResponse.json(progressData);
      } catch {
        return NextResponse.json({ error: "Failed to read progress", state: "UNKNOWN" }, { status: 500 });
      }
    }

    if (action === "mechanism-search-cancel") {
      // ─── Cancel a running search job ───
      const fs = await import("fs");
      const jobId = typeof body.jobId === "string" ? body.jobId.trim() : "";
      if (!jobId || !/^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i.test(jobId)) {
        return NextResponse.json({ error: "Valid jobId is required" }, { status: 400 });
      }

      const jobDir = path.join(JOB_BASE_DIR, jobId);
      const progressFile = path.join(jobDir, "progress.json");

      if (!fs.existsSync(progressFile)) {
        return NextResponse.json({ error: "Job not found" }, { status: 404 });
      }

      try {
        const progressData = JSON.parse(fs.readFileSync(progressFile, "utf-8"));
        if (progressData.state !== "RUNNING" && progressData.state !== "STARTING") {
          return NextResponse.json({ message: "Job already finished", state: progressData.state });
        }

        // Write CANCELLED state — the Python process will still run but
        // the frontend will stop polling and discard the result
        fs.writeFileSync(progressFile, JSON.stringify({
          ...progressData,
          state: "CANCELLED",
          updated_at: new Date().toISOString(),
        }), "utf-8");

        return NextResponse.json({ message: "Job cancelled", jobId, state: "CANCELLED" });
      } catch {
        return NextResponse.json({ error: "Failed to cancel job" }, { status: 500 });
      }
    }

    if (action === "mechanism-search") {
      const reactants = validateSmiles(body.reactants);
      const products = validateSmiles(body.products);
      // Allow empty products (forward-only mode)
      if (!reactants) {
        return NextResponse.json({ error: "Valid reactants SMILES are required" }, { status: 400 });
      }

      // Check for RXN file input (from Ketcher editor with atom mapping)
      const rxnText = typeof body.rxn_text === "string" ? body.rxn_text.trim() : "";
      if (!rxnText && !products) {
        return NextResponse.json({ error: "Both valid reactants and products SMILES are required" }, { status: 400 });
      }

      const maxConfigs = Math.min(Math.max(Math.floor(Number(body.max_configs) || 300), 50), 2000);
      // Max Bond Length: user-controlled slider (0-11 Å), used for scoring only (no hard cutoff)
      const maxBondLength = Math.min(Math.max(Number(body.max_bond_length) ?? 11, 0), 20);
      const reactantPdb = typeof body.reactant_pdb === "string" ? body.reactant_pdb.trim() : "";
      const productPdb = typeof body.product_pdb === "string" ? body.product_pdb.trim() : "";

      // Validate residues if provided
      let residuesJson = '';
      if (Array.isArray(body.residues)) {
        if (body.residues.length > 20) {
          return NextResponse.json({ error: "Maximum 20 residues allowed" }, { status: 400 });
        }
        for (const r of body.residues) {
          if (typeof r !== 'object' || !r.res_name || typeof r.res_num !== 'number') {
            return NextResponse.json({ error: "Each residue must have res_name (string) and res_num (number)" }, { status: 400 });
          }
        }
        residuesJson = JSON.stringify(body.residues);
      }

      const pdbId = typeof body.pdb_id === 'string' ? body.pdb_id.trim().toUpperCase() : '';
      const pdbText = typeof body.pdb_text === 'string' ? body.pdb_text.trim() : '';

      // NOTE: Do NOT use Flask for mechanism-search — Flask preloads all 51647 rules
      // (~6GB RAM). Mechanism search spawns its own subprocess that also loads rules,
      // causing OOM kill on memory-constrained environments. Always use child_process.
      // Direct child_process
      const scriptArgs: string[] = [
        path.join(PREDICTION_DIR, "mechanism_search.py"),
        "--reactants", reactants,
        "--products", products || "",
        "--max-configs", String(maxConfigs),
        "--max-rules", "0",  // 0 = load all rules
        "--max-bond-length", String(maxBondLength),
      ];

      if (reactantPdb) {
        scriptArgs.push("--reactant-pdb", reactantPdb);
      }
      if (productPdb) {
        scriptArgs.push("--product-pdb", productPdb);
      }
      if (residuesJson) {
        scriptArgs.push("--residues", residuesJson);
      }
      if (pdbId) {
        scriptArgs.push("--pdb-id", pdbId);
      }
      // Write PDB text to temp file to avoid E2BIG (argument list too long)
      let pdbTempFile: string | null = null;
      if (pdbText && pdbText.length > 100) {
        const fs = await import("fs");
        const os = await import("os");
        const path = await import("path");
        pdbTempFile = path.join(os.tmpdir(), `pdb_${Date.now()}.pdb`);
        fs.writeFileSync(pdbTempFile, pdbText, "utf-8");
        scriptArgs.push("--pdb-text-file", pdbTempFile);
      }
      // Write RXN text to temp file (atom-mapped reaction from Ketcher)
      let rxnTempFile: string | null = null;
      if (rxnText) {
        const fs = await import("fs");
        const os = await import("os");
        const path = await import("path");
        rxnTempFile = path.join(os.tmpdir(), `rxn_${Date.now()}.rxn`);
        fs.writeFileSync(rxnTempFile, rxnText, "utf-8");
        scriptArgs.push("--rxn-text-file", rxnTempFile);
      }

      // Write ligand_mappings to temp JSON file (precise PDB ligand coordinates)
      // Declared before the if-block so the finally cleanup can safely reference it
      // even if writeFileSync throws (e.g. permission error).
      let ligandMappingsTempFile: string | null = null;
      const ligandMappings = body.ligand_mappings;
      if (Array.isArray(ligandMappings) && ligandMappings.length > 0) {
        // Validate each mapping has required fields
        const validMappings = ligandMappings.filter(
          (m: any) => m && typeof m.smiles === 'string' && m.smiles.trim()
        );
        if (validMappings.length > 0) {
          const fs = await import("fs");
          const os = await import("os");
          const path = await import("path");
          ligandMappingsTempFile = path.join(os.tmpdir(), `ligand_mappings_${Date.now()}.json`);
          fs.writeFileSync(ligandMappingsTempFile, JSON.stringify(validMappings), "utf-8");
          scriptArgs.push("--ligand-mappings-file", ligandMappingsTempFile);
        }
      }

      try {
        const raw = await runPython(scriptArgs, 600000);
        // mechanism_search.py logs to stdout; extract the last JSON line
        const lines = raw.trim().split('\n');
        const jsonLine = lines.find(l => l.trim().startsWith('{')) || raw;
        const data = JSON.parse(jsonLine);
        return NextResponse.json(data);
      } finally {
        // Clean up temp PDB file
        if (pdbTempFile) {
          try {
            const fs = await import("fs");
            fs.unlinkSync(pdbTempFile);
          } catch { /* ignore */ }
        }
        // Clean up temp RXN file
        if (rxnTempFile) {
          try {
            const fs = await import("fs");
            fs.unlinkSync(rxnTempFile);
          } catch { /* ignore */ }
        }
        // Clean up temp ligand mappings file
        if (ligandMappingsTempFile) {
          try {
            const fs = await import("fs");
            fs.unlinkSync(ligandMappingsTempFile);
          } catch { /* ignore */ }
        }
      }

    }

    // ---- Feature 1: Residue Info Analysis ----
    if (action === "residue-info") {
      const reactionSmarts = validateSmarts(body.reaction_smarts);
      if (!reactionSmarts) {
        return NextResponse.json({ error: "Valid reaction SMARTS is required" }, { status: 400 });
      }
      const substrateSmiles = validateSmiles(body.substrate_smiles) || "";
      const productSmiles = validateSmiles(body.product_smiles) || "";

      const args = [path.join(PREDICTION_DIR, "residue_info.py"), reactionSmarts];
      if (substrateSmiles) args.push(substrateSmiles);
      if (productSmiles) args.push(productSmiles);

      const result = await runPython(args, 30000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Feature 2: Database Cross-Reference Links ----
    if (action === "rule-links") {
      if (!body.rule || typeof body.rule !== "object") {
        return NextResponse.json({ error: "Valid rule object is required" }, { status: 400 });
      }
      const ruleJson = JSON.stringify(body.rule);
      if (ruleJson.length > 5000) {
        return NextResponse.json({ error: "Rule data too large" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "database_links.py"), "rule", ruleJson,
      ], 10000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "molecule-links") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "database_links.py"), "molecule", smiles,
      ], 10000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Feature 3: Conservation Analysis ----
    if (action === "conservation") {
      const payload: Record<string, unknown> = {};
      if (body.mcsa_id != null) payload.mcsa_id = Number(body.mcsa_id);
      if (body.mechanism_id != null) payload.mechanism_id = Number(body.mechanism_id);
      if (body.reaction_smarts) payload.reaction_smarts = String(body.reaction_smarts);
      if (Array.isArray(body.residue_roles)) payload.residue_roles = body.residue_roles;
      if (!Object.keys(payload).length) {
        return NextResponse.json({ error: "At least one parameter is required (mcsa_id, mechanism_id, reaction_smarts, residue_roles)" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "conservation_analyzer.py"), JSON.stringify(payload),
      ], 30000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Feature 4: Batch Prediction ----
    if (action === "batch-predict") {
      if (!Array.isArray(body.smiles_list) || body.smiles_list.length === 0) {
        return NextResponse.json({ error: "Non-empty smiles_list array is required" }, { status: 400 });
      }
      if (body.smiles_list.length > 50) {
        return NextResponse.json({ error: "Maximum 50 SMILES per batch" }, { status: 400 });
      }
      // Validate each SMILES
      const validatedList: string[] = [];
      for (const s of body.smiles_list) {
        const v = validateSmiles(s);
        if (!v) {
          return NextResponse.json({ error: `Invalid SMILES in list: ${String(s)}` }, { status: 400 });
        }
        validatedList.push(v);
      }
      const maxSteps = Math.min(Math.max(Math.floor(Number(body.max_steps) || 3), 1), 10);
      const workers = Math.min(Math.max(Math.floor(Number(body.workers) || 3), 1), 8);

      const result = await runPython([
        path.join(PREDICTION_DIR, "batch_predict.py"),
        ...validatedList,
        "--max-steps", String(maxSteps),
        "--workers", String(workers),
      ], 300000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- PDB Handler ----
    if (action === "pdb-fetch") {
      const pdbId = typeof body.pdb_id === "string" ? body.pdb_id.trim().toUpperCase() : "";
      if (!pdbId || !/^[A-Z0-9]{4}$/.test(pdbId)) {
        return NextResponse.json({ error: "Valid 4-character PDB ID is required" }, { status: 400 });
      }

      // Try Flask backend first (avoids Python path issues on local deploy)
      const flaskResult = await tryFlask("pdb-fetch", { pdb_id: pdbId }, 30000);
      if (flaskResult) return flaskResult;

      // Fallback: direct Python script
      const result = await runPython([
        path.join(PREDICTION_DIR, "pdb_handler.py"), "fetch", pdbId,
      ], 30000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "pdb-parse") {
      const pdbText = typeof body.pdb_text === "string" ? body.pdb_text.trim() : "";
      if (!pdbText || pdbText.length < 100) {
        return NextResponse.json({ error: "Valid PDB format text is required (min 100 chars)" }, { status: 400 });
      }
      if (pdbText.length > 500000) {
        return NextResponse.json({ error: "PDB text too large (max 500KB)" }, { status: 400 });
      }

      // Try Flask backend first
      const flaskResult = await tryFlask("pdb-parse", { pdb_text: pdbText }, 30000);
      if (flaskResult) return flaskResult;

      // Fallback: direct Python script
      const result = await runPython([
        path.join(PREDICTION_DIR, "pdb_handler.py"), "parse",
      ], 30000, pdbText);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "pdb-raw-fetch") {
      const pdbId = typeof body.pdb_id === "string" ? body.pdb_id.trim().toUpperCase() : "";
      if (!pdbId || !/^[A-Z0-9]{4}$/.test(pdbId)) {
        return NextResponse.json({ error: "Valid 4-character PDB ID is required" }, { status: 400 });
      }
      // Try Flask backend first
      try {
        const flaskResult = await tryFlask("pdb-raw-fetch", { pdb_id: pdbId }, 30000);
        if (flaskResult) return flaskResult;
      } catch { /* fall through */ }
      try {
        const result = await runPython([
          path.join(PREDICTION_DIR, "pdb_handler.py"), "fetch-raw", pdbId,
        ], 30000);
        // Ensure result is a string (runPython may return Buffer)
        const pdbText = typeof result === "string" ? result : String(result);
        return NextResponse.json({ pdb_text: pdbText });
      } catch (error: unknown) {
        const msg = error instanceof Error ? error.message : "Failed to fetch PDB file";
        return NextResponse.json({ error: msg }, { status: 500 });
      }
    }

    if (action === "pdb-active-site") {
      const pdbId = typeof body.pdb_id === "string" ? body.pdb_id.trim().toUpperCase() : "";
      if (!pdbId || !/^[A-Z0-9]{4}$/.test(pdbId)) {
        return NextResponse.json({ error: "Valid 4-character PDB ID is required" }, { status: 400 });
      }

      // Try Flask backend first
      const flaskResult = await tryFlask("pdb-active-site", { pdb_id: pdbId }, 45000);
      if (flaskResult) return flaskResult;

      // Fallback: direct Python script
      const result = await runPython([
        path.join(PREDICTION_DIR, "pdb_handler.py"), "active-site", pdbId,
      ], 45000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Feature 5: Arrow Editor ----
    if (action === "arrow-atoms") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "arrow_editor.py"), "atoms", smiles,
      ], 15000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "arrow-apply") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }
      if (!Array.isArray(body.arrows) || body.arrows.length === 0) {
        return NextResponse.json({ error: "Non-empty arrows array is required" }, { status: 400 });
      }
      // Validate arrows structure
      for (const arrow of body.arrows) {
        if (!Array.isArray(arrow.source_atoms) || !Array.isArray(arrow.target_atoms) || typeof arrow.electrons !== "number") {
          return NextResponse.json({ error: "Each arrow must have source_atoms, target_atoms, and electrons" }, { status: 400 });
        }
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "arrow_editor.py"), "apply", smiles, JSON.stringify(body.arrows),
      ], 15000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "arrow-infer") {
      const reactant = validateSmiles(body.reactant_smiles);
      const product = validateSmiles(body.product_smiles);
      if (!reactant || !product) {
        return NextResponse.json({ error: "Both reactant_smiles and product_smiles are required" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "arrow_editor.py"), "infer", reactant, product,
      ], 30000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- ChEBI Fetch ----
    if (action === "chebi-fetch") {
      const chebiId = typeof body.chebi_id === "string" ? body.chebi_id.trim() : "";
      if (!chebiId) {
        return NextResponse.json({ error: "Valid ChEBI ID is required" }, { status: 400 });
      }
      // Accept formats: CHEBI:15377 or 15377
      const chebiNum = chebiId.replace(/^CHEBI:/i, "");
      if (!/^\d+$/.test(chebiNum)) {
        return NextResponse.json({ error: "ChEBI ID must be numeric (e.g. CHEBI:15377 or 15377)" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "chebi_fetch.py"), chebiNum,
      ], 15000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Ligand Comparison (MCS Highlight) ----
    if (action === "ligand-compare") {
      const smilesA = validateSmiles(body.smiles_a);
      const smilesB = validateSmiles(body.smiles_b);
      if (!smilesA || !smilesB) {
        return NextResponse.json({ error: "Both smiles_a and smiles_b are required" }, { status: 400 });
      }
      const w = Math.min(Math.max(Number(body.width) || 350, 100), 800);
      const h = Math.min(Math.max(Number(body.height) || 300, 100), 600);

      const result = await runPython([
        path.join(PREDICTION_DIR, "ligand_compare.py"), "compare", smilesA, smilesB,
        "--width", String(w), "--height", String(h),
      ], 15000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "ligand-pdb-smiles") {
      const resName = typeof body.res_name === "string" ? body.res_name.trim().toUpperCase() : "";
      if (!resName || resName.length > 10) {
        return NextResponse.json({ error: "Valid res_name (max 10 chars) is required" }, { status: 400 });
      }
      if (!/^[A-Z0-9]+$/.test(resName)) {
        return NextResponse.json({ error: "res_name must be alphanumeric" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "ligand_compare.py"), "pdb_smiles", resName,
      ], 5000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    if (action === "ligand-similarity") {
      const smilesA = validateSmiles(body.smiles_a);
      const smilesB = validateSmiles(body.smiles_b);
      if (!smilesA || !smilesB) {
        return NextResponse.json({ error: "Both smiles_a and smiles_b are required" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "ligand_compare.py"), "similarity", smilesA, smilesB,
      ], 10000);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Ligand SMILES Extraction from PDB Structure ----
    if (action === "ligand-extract-smiles") {
      const resName = typeof body.res_name === "string" ? body.res_name.trim().toUpperCase() : "";
      if (!resName || resName.length > 10) {
        return NextResponse.json({ error: "Valid res_name (max 10 chars) is required" }, { status: 400 });
      }
      if (!/^[A-Z0-9]+$/.test(resName)) {
        return NextResponse.json({ error: "res_name must be alphanumeric" }, { status: 400 });
      }
      const pdbText = typeof body.pdb_text === "string" ? body.pdb_text.trim() : "";
      const chain = typeof body.chain === "string" ? body.chain.trim() : "";
      const resNum = typeof body.res_num === "number" ? body.res_num : undefined;

      const args = [path.join(PREDICTION_DIR, "ligand_compare.py"), "extract_smiles", resName];
      if (chain) args.push(chain);
      if (resNum !== undefined && chain) args.push(String(resNum));

      const result = await runPython(args, 15000, pdbText || undefined);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    // ---- Atom Mapping Verification & Correction ----
    if (action === "verify-atom-mapping") {
      const reactants = typeof body.reactants === "string" ? body.reactants.trim() : "";
      const products = typeof body.products === "string" ? body.products.trim() : "";
      if (!reactants || !products) {
        return NextResponse.json({ error: "Both reactants and products SMILES are required" }, { status: 400 });
      }

      const result = await runPython([
        path.join(PREDICTION_DIR, "mechanism_search.py"),
        "--reactants", reactants,
        "--products", products,
        "--verify-mapping",
      ], 15000);
      try {
        const data = JSON.parse(result.trim());
        return NextResponse.json(data);
      } catch {
        return NextResponse.json({ error: "Failed to parse verification result" }, { status: 500 });
      }
    }

    // ---- Ligand Distances to Active Site Center of Geometry ----
    if (action === "ligand-distances") {
      const pdbText = typeof body.pdb_text === "string" ? body.pdb_text.trim() : "";
      const selectedResidues = Array.isArray(body.selected_residues) ? body.selected_residues : [];
      if (!pdbText || pdbText.length < 100) {
        return NextResponse.json({ error: "Valid PDB text is required (min 100 chars)" }, { status: 400 });
      }
      if (selectedResidues.length === 0) {
        return NextResponse.json({ error: "At least one selected residue is required" }, { status: 400 });
      }
      // Validate each residue entry
      for (const sr of selectedResidues) {
        if (typeof sr.chain !== "string" || typeof sr.res_num !== "number") {
          return NextResponse.json({ error: "Each selected residue must have 'chain' (string) and 'res_num' (number)" }, { status: 400 });
        }
      }

      const residuesJson = JSON.stringify(selectedResidues);
      const result = await runPython([
        path.join(PREDICTION_DIR, "pdb_handler.py"), "ligand-distances", residuesJson,
      ], 30000, pdbText);
      try {
        const data = JSON.parse(result.trim());
        return NextResponse.json(data);
      } catch {
        return NextResponse.json({ error: "Failed to parse ligand distances result" }, { status: 500 });
      }
    }

    // ---- Substrate Overlay Pipeline ----
    if (action === "overlay-substrate") {
      const smiles = validateSmiles(body.smiles);
      if (!smiles) {
        return NextResponse.json({ error: "Valid SMILES string is required" }, { status: 400 });
      }
      const resName = typeof body.res_name === "string" ? body.res_name.trim().toUpperCase() : "";
      if (!resName || resName.length > 10) {
        return NextResponse.json({ error: "Valid res_name (max 10 chars) is required" }, { status: 400 });
      }
      if (!/^[A-Z0-9]+$/.test(resName)) {
        return NextResponse.json({ error: "res_name must be alphanumeric" }, { status: 400 });
      }
      const pdbText = typeof body.pdb_text === "string" ? body.pdb_text.trim() : "";
      if (!pdbText || pdbText.length < 100) {
        return NextResponse.json({ error: "Valid PDB text is required (min 100 chars)" }, { status: 400 });
      }
      if (pdbText.length > 500000) {
        return NextResponse.json({ error: "PDB text too large (max 500KB)" }, { status: 400 });
      }
      const chain = typeof body.chain === "string" ? body.chain.trim() : "";
      const resNum = typeof body.res_num === "number" ? body.res_num : undefined;
      const maxConformers = Math.min(Math.max(Number(body.max_conformers) || 50, 5), 200);
      const rmsdThreshold = Math.min(Math.max(Number(body.rmsd_threshold) || 2.0, 0.5), 10.0);

      const args = [
        path.join(PREDICTION_DIR, "overlay_substrate.py"), smiles, "-",
        "--res-name", resName,
        "--max-conformers", String(maxConformers),
        "--rmsd-threshold", String(rmsdThreshold),
      ];
      if (chain) args.push("--chain", chain);
      if (resNum !== undefined) args.push("--res-num", String(resNum));

      const result = await runPython(args, 30000, pdbText);
      const data = JSON.parse(result.trim());
      return NextResponse.json(data);
    }

    return NextResponse.json({ error: `Unknown action: ${String(action)}` }, { status: 400 });
  } catch (error: unknown) {
    const msg = error instanceof Error ? error.message : "Internal error";
    console.error(`[API] Action "${actionName}" failed:`, msg);
    return NextResponse.json({
      error: msg,
      hint: "Make sure the prediction service is running (python mini-services/prediction-service/index.py) and you have network access to RCSB (files.rcsb.org)",
    }, { status: 500 });
  }
}

export async function GET() {
  return NextResponse.json({
    status: "ok",
    service: "M-CSA Prediction API",
    architecture: {
      primary: "child_process (Python scripts on-demand)",
      secondary: "Flask backend (port 3003) for mechanism-search with rule caching",
      shared_module: "shared.py (simplify_smarts, get_builtin_rules, load_rules)",
    },
    endpoints: {
      "POST /api": {
        actions: [
          "health", "mol-info", "mol-svg",
          "mechanism-search-start", "mechanism-search-status", "mechanism-search-cancel",
          "mechanism-search",
          "residue-info", "rule-links", "molecule-links", "conservation",
          "batch-predict", "arrow-atoms", "arrow-apply", "arrow-infer",
          "pdb-fetch", "pdb-parse", "pdb-active-site",
          "chebi-fetch",
        ],
      },
      "residue-info": {
        description: "Analyze catalytic residue roles from reaction SMARTS",
        example: { action: "residue-info", reaction_smarts: "[CH2:1][OH:2]>>[CH:1]=[O:3]" },
      },
      "rule-links": {
        description: "Generate external database cross-reference links for a rule",
        example: { action: "rule-links", rule: { mcsa_id: 123, mechanism_id: 5 } },
      },
      "molecule-links": {
        description: "Generate molecule search links for external databases",
        example: { action: "molecule-links", smiles: "CCO" },
      },
      "conservation": {
        description: "Analyze active site conservation for a mechanism",
        example: { action: "conservation", reaction_smarts: "[CH2:1][OH:2]>>[CH:1]=[O:3]", residue_roles: [{ role: "nucleophile" }] },
      },
      "batch-predict": {
        description: "Batch prediction for multiple SMILES (up to 50)",
        example: { action: "batch-predict", smiles_list: ["CCO", "CC(=O)O"], workers: 3 },
      },
      "arrow-atoms": {
        description: "Get molecule atom indices for arrow editing",
        example: { action: "arrow-atoms", smiles: "CCO" },
      },
      "arrow-apply": {
        description: "Apply arrow edits to a molecule",
        example: { action: "arrow-apply", smiles: "CCO", arrows: [{ source_atoms: [2], target_atoms: [1], electrons: 2 }] },
      },
      "arrow-infer": {
        description: "Infer arrow sequence from reactant→product transformation",
        example: { action: "arrow-infer", reactant_smiles: "CCO", product_smiles: "CC=O" },
      },
    },
    rules_source: "rules_nat_met.xlsx (51,637 M-CSA rules + 10 built-in)",
  });
}
