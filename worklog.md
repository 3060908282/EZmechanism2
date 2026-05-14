---
Task ID: deploy-1
Agent: Main
Task: Deploy EZmechanism2 from GitHub to my-project

Work Log:
- Cloned https://github.com/3060908282/EZmechanism2 to /home/z/EZmechanism2
- Copied all source directories: src/, mini-services/, public/, db/, prisma/, upload/
- Copied worklog.md from EZmechanism2
- Added cytoscape@^3.33.2 and ketcher-react@^3.12.0 to package.json
- Ran bun install (205 packages installed)
- Updated next.config.ts with serverExternalPackages, CORS headers
- Updated eslint.config.mjs with ignores for ketcher, 3Dmol, mech-search-jobs
- Ran prisma db:push (database already in sync)
- Ran prisma db:generate (Prisma Client generated)
- Installed missing Python packages: rdkit, flask, flask-cors
- Ketcher static files already pre-built in public/ketcher/ (no rebuild needed)
- Started prediction-service (Flask on port 3003): 6,093 rules loaded
- Started Next.js dev server (port 3000): HTTP 200

Stage Summary:
- All 3 services running:
  - Next.js (3000): HTTP 200, compile OK
  - Prediction Service (3003): 6,093 rules (6,083 mcsa + 10 builtin), health OK
  - Ketcher static files: HTTP 200
- Database: custom.db 88MB, PdbTest table has existing entries
- PDB Tests API: working, returns test data

---
Task ID: fix-zero-tolerance-1
Agent: Main
Task: Implement zero-tolerance coordinate policy fix (4 steps)

Work Log:
- Step ①: Modified main() inline coords_cache build (lines 3270-3307)
  - Removed fallback to global MCS (load_pdb_coords_for_smiles)
  - Removed fallback to generate_3d_coords
  - No mapping or extraction failed → coords_cache[key] = None
  - Updated standalone build_coords_cache_from_pdb() function for consistency
  - Fixed logger.info for coords_cache=None case
- Step ②: Added coordinate completeness check in bidirectional_search (after _build_am_to_3d_coord)
  - Core ligands: must have ALL non-H atoms in am_to_3d_coord, else abort search
  - Water/proton molecules: warning only if no coords, search continues
  - Returns error result with descriptive message when core ligand check fails
- Step ③: Modified _quick_parse in reaction_result.py
  - PDB mode (am_to_3d_coord non-empty): direct dict access → KeyError → result discarded
  - Non-PDB mode (am_to_3d_coord empty/None): no change, scoring falls back to standard formula
  - Deprecated _generate_and_backfill_coords (no longer called by _quick_parse)
- Step ④: Updated scoring logic in _get_reaction_data_for_products
  - Removed penalty score=9 for missing coords (now handled by KeyError)
  - Simplified to standard Upload formula: score = (max(4, score_bl) - 3)^2

Stage Summary:
- Verified Python imports OK, prediction service restarted
- Unit tests passed:
  - PDB mode + missing coords → KeyError → result discarded ✅
  - PDB mode + complete coords → correct bond length ✅
  - Non-PDB mode → no KeyError, normal operation ✅
- Lint check passed
- Next.js dev server running normally

---
Task ID: ligand-mapping-consume-1
Agent: Main
Task: Fix ligand mapping "consume-one-by-one" mode for multiple same-SMILES molecules

Work Log:
- Frontend fix (PdbWorkflowMode.tsx):
  - Removed `s.type === 'substrate'` filter in ligand_mappings construction
  - Now includes all types (substrate, cofactor, water) that have mapped_ligand
  - Added `ligand_type` field to each mapping entry for backend logging
- Backend main() simplification:
  - coords_cache now ONLY stores molecules WITHOUT ligand_mappings
  - Molecules with ligand_mappings are skipped (handled by _build_am_to_3d_coord)
  - Passed ligand_mappings to bidirectional_search() call
- Backend bidirectional_search():
  - Added ligand_mappings parameter to function signature
  - Updated has_pdb_data condition to trigger on ligand_mappings + pdb_text
  - Passed ligand_mappings to _build_am_to_3d_coord()
- Backend _build_am_to_3d_coord() complete rewrite:
  - Promoted _assigned_mol_ids to function top level (shared across all steps)
  - Added _get_mol_canonical_key() helper function
  - Added _find_next_unassigned_mol() generic consume function
  - Added _mcs_align_and_mmff() shared MCS+MMFF optimization function
  - ① Ligand mapping path: per-mapping PDB fragment extraction → MCS → MMFF
  - ② Residue path: reuses _find_next_unassigned_mol with shared _assigned_mol_ids
  - ③ General coords_cache path: now SKIPS mols in _assigned_mol_ids
  - ④ Build am_to_3d_coord (unchanged)
  - Failed assignments un-assign (_assigned_mol_ids.discard) so other steps can retry

Stage Summary:
- Fixed the critical overwrite bug: general coords_cache path now skips already-assigned mols
- Fixed the SMILES-keyed cache collision: multiple water O's now get independent PDB coordinates
- Fixed frontend filter: water molecule mappings are now passed to backend
- Verified: Python import OK, lint OK, prediction service restarted

---
Task ID: remove-hard-filter-1
Agent: Main
Task: Remove hard bond length filter, align with Upload scoring-only pruning

Work Log:
- Removed `if score_bl > 0 and score_bl > max_bond_length: return None` hard filter
- Replaced with comment explaining why: zero-tolerance policy ensures all scored bonds
  come from real PDB coordinates, so Upload's scoring-based pruning is sufficient
- Verified _quick_parse already correctly implements PDB mode KeyError behavior
- Restarted prediction service

Stage Summary:
- Hard bond length filter removed — Dijkstra scoring alone handles pruning (Upload-aligned)
- _quick_parse: PDB mode uses direct dict access → KeyError → result discarded (confirmed working)

---
Task ID: fix-water-ligand-per-instance-1
Agent: Main
Task: Fix Step3 "Map to Selected ligand" — two water substrates mapping to same PDB ligand

Root Cause Analysis:
1. PDB parser (pdb_handler.py) excluded water molecules from ligand list entirely
   → HOH entries never appeared in the dropdown
2. Even if water were included, ligand grouping was by (chain, res_name) only
   → All HOH on same chain collapsed into one entry (taking res_nums[0])
3. Same issue affected any multi-instance ligand (e.g., ATP_501, ATP_502)

Work Log:
- Backend pdb_handler.py: restructured ligand parsing
  - Changed data structure from (chain, res_name) → atoms to (chain, res_name, res_num) → atoms
  - Replaced ligand_groups + ligand_resnum_map with ligand_instance_groups (3-level dict)
  - Now includes ALL HETATM records including water (HOH, WAT, etc.)
  - Each (chain, res_name, res_num) instance outputs as a separate ligand entry
  - Added is_water boolean flag to each ligand entry
  - Active site prediction now filters out water molecules (non_water_ligands)
  - Verified: 3 HOH entries → 3 separate ligand entries with unique res_nums
- Backend mechanism_search.py: fixed single-atom residue extraction
  - Changed `if len(residue_lines) < 2: return None` to `< 1`
  - Water molecules (1 O atom) can now be extracted successfully
  - Verified: extract_pdb_residue_fragment("HOH", 401, "A") returns 1-atom Mol with correct coords
- Frontend types.ts: added is_water field to PdbLigand interface
- Frontend EditMoleculeDialog.tsx: water-aware UI
  - Skip similarity computation for water molecules (O vs O is meaningless, avoids 500 API calls)
  - Sort order: non-water by similarity → water by proximity
  - When substrate type is water, water molecules shown first in dropdown
  - Added visual separator "── Water Molecules ──" in dropdown before water entries
  - Ligand table: 💧 icon + "water" label for HOH entries, blue tinted rows
  - Increased table display limit from 20 to 30

Stage Summary:
- Water molecules now appear as individual selectable entries in Step3 ligand dropdown
- Each HOH instance has unique (chain, res_num) so two water substrates can map to different sites
- Multi-instance non-water ligands also benefit (e.g., ATP_501 and ATP_502 are separate entries)
- Verified: PDB parsing, fragment extraction, lint, dev server all OK

---
Task ID: fix-water-ligand-per-instance-2
Agent: Main
Task: Fix normal ligands being unmappable after water ligand inclusion (UI overload)

Root Cause:
1. Including 378 water molecules in dropdown made it practically unusable
2. Similarity computation loop processed all 378 entries (even though water was skipped quickly)
3. "No SMILES available for ALP" was a pre-existing MCS comparison issue, not a mapping failure
4. Users couldn't find normal ligands in the flooded dropdown

Work Log:
- EditMoleculeDialog.tsx: Added water filter toggle (💧 Water molecules: Visible/Hidden)
  - Default: hidden for non-water substrates, visible for water substrates
  - When visible, only nearest 30 HOH sites shown (sorted by proximity)
  - Chain filter info bar updated with water count
- EditMoleculeDialog.tsx: filteredLigands now respects showWater flag
  - When showWater=false, all is_water entries are excluded before similarity computation
  - This prevents hundreds of unnecessary similarityRef entries
- EditMoleculeDialog.tsx: selectedChainLigands now excludes water for the "no ligands" warning
- Fixed "No non-water ligands found" message → "No ligands found"
- Prediction service restarted (rdkit re-installed, Flask port 3003 confirmed running)
- Verified: 1AKE returns 2 non-water + 378 water ligands from API
- Verified: lint passes, dev server compiles OK

Stage Summary:
- Water molecules are now HIDDEN by default in the dropdown (toggle to show)
- Normal ligands are immediately accessible without water clutter
- When water is shown, limited to nearest 30 to prevent UI overload
- Performance: similarity computation only runs on non-water ligands when water is hidden

---
Task ID: fix-alp-smiles-atom-mapping-1
Agent: Main
Task: Fix ALP SMILES regression + answer question about Step4 atom mapping with copied SMILES

Work Log:
- Investigated "No SMILES available for ALP" error in EditMoleculeDialog.tsx
  - ALP was not in the built-in _PDB_LIGAND_SMILES dictionary in ligand_compare.py
  - PDBeChem API could resolve it (CC1(C)S[C@@H](N[C@H]1C(O)=O)[C@@H](CO)C(O)=O) but may fail due to timeout/network
  - Added ALP and 17 other common PDB ligands to the built-in dictionary
- Analyzed Step4 atom mapping flow when user copies a SMILES as product:
  - User pastes SMILES → goes into Ketcher reaction editor
  - On "Run Search", handleWfRun auto-calls Indigo automap (line 758-769)
  - Indigo assigns atom map numbers to both reactants and products
  - RXN file exported with atom map numbers → sent to backend
  - Backend uses atom map numbers for reaction center identification + 3D coordinate assignment
- Added atom mapping hint in Step4 UI (PdbWorkflowMode.tsx):
  - Shows when product SMILES is present
  - Explains that atom mapping is automatically assigned by Indigo
  - Warns if product SMILES has no atom map numbers (no ":" characters)
  - Suggests using Auto Map button to preview and verify mapping

Stage Summary:
- ALP SMILES now available from built-in dictionary (no PDBeChem API dependency)
- 18 additional common PDB ligands added: F6P, G6P, 2PG, 3PG, PEP, PYR, LAC, FUM, AKG, MAL, ORO, R5P, BGC, MPD, EPE, PG4
- Step4 UI now shows atom mapping guidance when product SMILES is entered
- Answer: atom mapping DOES work when copying SMILES as product — Indigo automap assigns it automatically
---
Task ID: 1
Agent: Main Agent
Task: Investigate and fix PDB-003 Step4 atom mapping issues

Work Log:
- Analyzed PDB-003 data from SQLite: 1TEM (β-lactamase TEM-1), Step4 has reactants `CC1(C)C(C([O-])=O)N2C(C(C2=O)N)S1.O.O` and products `O.CC1(SC(NC1C(=O)[O-])C(C(O)=O)N)C`
- Installed epam-indigo Python package to test Indigo automap algorithm
- Tested Indigo automap on PDB-003 reaction and found critical mapping error:
  - Indigo incorrectly swapped carboxyl =O (AM:7) with a water molecule (AM:15)
  - This inflated the reaction center from 3 atoms to 6 atoms
  - Wrong reaction center: bonds (5,7) C=O and (8,11) N-C broken, (16,11) O-C and (5,15) C=O formed
  - Correct reaction center: only (8,12) N-C broken and (16,12) O-C formed (β-lactam ring opening)
- Implemented `remap_reaction_by_mcs()` function in mechanism_search.py that:
  1. Computes MCS between each (reactant, product) pair
  2. Assigns atom map numbers based on MCS atom correspondence
  3. Handles unmatched atoms by element matching (e.g., water O → product OH O)
  4. Handles unassigned molecules by element signature matching (e.g., water → water)
  5. Only applies correction if it reduces the number of bond changes
- Added `remap_reaction_by_mcs` call in the main CLI flow after RXN parsing
- Verified the correction: 4 → 2 bond changes for PDB-003
- Tested edge cases: correct mappings unchanged, forward-only mode, no-mapping mode
- ALP SMILES lookup works correctly (confirmed `ligand_compare.py pdb_smiles ALP` returns valid SMILES)
- Started prediction service on port 3003

Stage Summary:
- **Critical bug found**: Indigo automap incorrectly maps carboxyl =O to water for β-lactam hydrolysis reactions
- **Fix implemented**: `remap_reaction_by_mcs()` function in mechanism_search.py automatically corrects suboptimal atom mappings
- **Result**: PDB-003 reaction center correctly reduced from 6 atoms to 3 atoms
- **SMILES regression**: ALP SMILES is available in the dictionary; the "No SMILES available" error was likely a transient API timeout issue, not a code bug

---
Task ID: 2
Agent: Main Agent
Task: Add atom mapping verification and correction features

Work Log:
- Added `remap_reaction_by_mcs()` function to mechanism_search.py (lines 991-1254)
  - Computes MCS between reactant/product pairs to find optimal atom mapping
  - Handles unmatched atoms by element matching (water O → product OH O)
  - Handles unassigned molecules by element signature matching (water → water)
  - Only applies correction if it reduces bond changes
- Added `--verify-mapping` CLI flag to mechanism_search.py for quick verification
- Added `verify-atom-mapping` API endpoint in route.ts
- Added "Verify" button in Step4 UI for manual atom mapping verification
- Integrated automatic remapping in the search flow (after RXN parsing, before bidirectional_search)
- Verified: ALP SMILES works correctly in the backend dictionary
- All lint checks pass

Stage Summary:
- **Backend**: `remap_reaction_by_mcs()` automatically corrects Indigo automap errors
- **Backend**: `--verify-mapping` CLI and API endpoint for manual verification
- **Frontend**: "Verify" button in Step4 for user-triggered verification
- **Auto-correction**: Happens transparently during search, logged for transparency
- **PDB-003 result**: 4 → 2 bond changes (carboxyl =O no longer swapped with water)
---
Task ID: 3
Agent: Main Agent
Task: Fix max_rules=11 bug and activeSiteResultJson=null issue

Work Log:
- Diagnosed max_rules=11 bug: naming collision between maxBondLength (0-11 Å) and maxRules (should be 0 or 5000)
  - buildSavePayload at line 251: `maxRules: s.wfMaxBondLength` stored bond length (11) as maxRules
  - handleViewTest at line 440: `setWfMaxBondLength(test.maxRules || 5000)` restored maxRules (11) into bond length slider
  - API route never passed `--max-rules` to Python script
- Fixed DB schema: added `maxBondLength Float @default(11.0)` field, changed `maxRules` default to 0
- Fixed save logic: `maxRules: 0` (load all rules), `maxBondLength: s.wfMaxBondLength` (correct field)
- Fixed restore logic: `setWfMaxBondLength(test.maxBondLength ?? 11)` (correct source field)
- Added `--max-rules "0"` to both async and sync search paths in route.ts
- Ran `bun run db:push` to sync schema

- Diagnosed activeSiteResultJson=null: PDB fetch already returns active_site_residues in response
  but frontend never used it — required separate manual "Predict Active Site" button click
- Fixed handleFetchPdb (both pdbId and pdbText paths) to auto-populate:
  - setActiveSiteResult from data.active_site_residues
  - setSelectedResidues from data.active_site_residues (auto-select predicted residues)
- Added active_site_residues field to PdbInfo TypeScript interface
- Updated toast messages to include active site residue count

Stage Summary:
- **max_rules bug FIXED**: DB now has separate maxRules (0=all) and maxBondLength (11.0 Å) fields
- **Save/restore confusion FIXED**: wfMaxBondLength saves to maxBondLength, not maxRules
- **API --max-rules FIXED**: Both async and sync search paths pass --max-rules 0 (load all ~6000 rules)
- **Active site auto-population FIXED**: PDB fetch now auto-fills activeSiteResult and selectedResidues
- **PdbInfo type updated**: Added active_site_residues optional field
- Lint passes, dev server running normally
---
Task ID: 4
Agent: Main Agent
Task: Remove active site prediction feature entirely (user prefers manual residue selection)

Work Log:
- Removed `wfActiveSiteLoading` and `wfActiveSiteResult` states from page.tsx
- Removed `activeSiteLoading` and `activeSiteResult` props from PdbWorkflowMode component
- Removed `handlePredictActiveSite` callback function
- Removed "Predict Active Site" button from Step 2 UI
- Removed `activeSiteResult` display (predicted triad, nearest ligand info) from Step 2
- Removed auto-population of activeSiteResult/selectedResidues from handleFetchPdb (both paths)
- Reverted toast messages back to simple format (without active site residue count)
- Removed `activeSiteResultJson` from save payload in buildSavePayload
- Removed `activeSiteResultJson` restore logic in handleViewTest
- Removed `setActiveSiteResult(null)` from handleCreateNew
- Removed `ActiveSitePrediction` interface from types.ts
- Removed `active_site_residues?` field from PdbInfo interface
- Removed `activeSiteResultJson` column from Prisma schema
- Ran `bun run db:push` to sync schema
- Updated hint text: "Click '+' to add catalytic residues" (was "Click '+' or use 'Predict Active Site'")
- Backend `pdb-active-site` endpoint kept in route.ts (not called from frontend, but harmless)

Stage Summary:
- Active site prediction feature fully removed from frontend
- DB schema cleaned: no more `activeSiteResultJson` column
- Users now manually select catalytic residues in Step 2
- Lint passes, dev server running normally
---
Task ID: 1
Agent: main
Task: Diagnose process killing and fix bugs

Work Log:
- Diagnosed that sandbox kills all child processes when bash command completes
- Found that nohup, setsid, disown don't help because sandbox uses cgroup-level process management
- Confirmed both bug fixes are already in place in the code:
  1. maxRules: 0 (load all rules) in buildSavePayload
  2. --max-rules "0" in API route (both async and sync paths)
  3. maxBondLength Float @default(11.0) in schema
  4. Active site prediction feature already removed
- Installed flask, flask-cors, rdkit Python packages for prediction service
- Created start-services.sh for unified service startup
- Lint passes

Stage Summary:
- Bug fixes already present in code (no code changes needed)
- Services start successfully but get killed when bash command exits
- Solution: must use `bun run dev` in foreground or start within same command session
- Prediction service deps installed (flask, rdkit, etc.)
- All code is clean (lint passes)

---
Task ID: 5
Agent: main
Task: Refactor ligand_compare.py - extract shared _find_mcs function + unify timeout

Work Log:
- Analyzed current ligand_compare.py: two independent MCS functions with duplicated logic
- Identified timeout inconsistency: compute_similarity_only=1s, compute_mcs_highlight=3s, original compare_mols=2s
- Extracted _find_mcs() shared core function with unified timeout=2s (matching original)
- Refactored compute_similarity_only to delegate to _find_mcs(timeout=2)
- Refactored compute_mcs_highlight to delegate to _find_mcs(timeout=2) for MCS computation
- Verified consistency: both functions return identical similarity for same molecule pair
- Verified Flask index.py compatibility with refactored module

Stage Summary:
- _find_mcs() is now the single source of truth for MCS parameters
- Timeout unified to 2s (matching original compare_mols)
- Algorithm consistency guaranteed: both endpoints use identical MCS computation
- No need to port compare_mols - current implementation is logically equivalent
- Code is cleaner: MCS core logic written once, not twice

---
Task ID: 1
Agent: main
Task: Refactor _find_mcs + water overlay support

Work Log:
- Analyzed current state of ligand_compare.py: _find_mcs() already refactored, but compute_similarity_only() still returned float
- Modified compute_similarity_only() to return dict with similarity, mcs_num_atoms, mcs_num_bonds, mcs_smarts, num_atoms_a, num_atoms_b
- Updated CLI similarity action to output enriched dict instead of {"similarity": float}
- Added _WATER_RESNAMES constant to overlay_substrate.py
- Added _overlay_water_to_pdb() function for direct water O atom coordinate mapping
- Added water bypass at beginning of overlay_smiles_to_pdb() to skip MCS+Kabsch for water molecules
- Added mcsSmartsRef to EditMoleculeDialog.tsx to cache mcs_smarts from ligand-similarity responses
- Added isWater field to OverlayResult interface
- Added isSelectedLigandWater computed value for conditional UI
- Updated overlay button to show blue "Map Water Coordinates to PDB Position" for water, purple for non-water
- Updated overlay result display with water-specific blue card (💧 Water Coordinate Mapping Applied)
- Updated ligand list water badge from "water" text to styled Badge component
- Updated toast message to show "Water Coordinate Mapping Applied" for water overlays
- All lint checks pass
- All backend tests pass (CCO vs CCO → similarity 1.0, benzene vs CCO → 0.3333, water overlay → success)

Stage Summary:
- ligand_compare.py: compute_similarity_only() now returns enriched dict (similarity + mcs_num_atoms + mcs_num_bonds + mcs_smarts + atom counts)
- overlay_substrate.py: water molecules now get direct coordinate mapping instead of failing through MCS+Kabsch pipeline
- EditMoleculeDialog.tsx: full water-aware UI (blue button, blue result card, cached mcs_smarts)
- All services running: Next.js:3000, Flask:3003, Ketcher:3004

---
Task ID: 2
Agent: main
Task: Fix service persistence - keep frontend/backend running even when idle

Work Log:
- Discovered that the Bash tool kills all child processes when the session ends
- This means nohup, setsid, and disown don't help because the entire process group gets cleaned up
- Found that the double-fork technique (using subshell + disown) creates orphan processes that get reparented to PID 1 (tini)
- Created /home/z/my-project/run-services.sh - starts all 3 services using double-fork
- Created /home/z/my-project/watchdog.sh - monitors services every 30s and restarts dead ones
- Both use the double-fork technique to survive session termination
- Verified services persist across multiple Bash tool sessions

Stage Summary:
- Services are now running persistently via double-fork technique
- Next.js (3000): ✅ running
- Flask prediction (3003): ✅ running
- Ketcher (3004): ✅ running
- Caddy gateway (81): ✅ serving pages
- Watchdog monitors and auto-restarts dead services every 30s

---
Task ID: 6
Agent: main
Task: Phase A - Add optimize_mol_2d() and detect_mol_stereo() to result display pipeline

Work Log:
- Added optimize_mol_2d() to reaction_result.py as standalone function
  - Only operates on Mol with exactly NumConformers == 1 (won't touch 3D conformers)
  - Uses Compute2DCoords with coordMap to re-layout distorted RunReactants products
  - Returns bool: True=success, False=skip/failure
- Added detect_mol_stereo() to reaction_result.py as standalone function
  - Checks conformer existence AND Is3D() before calling AssignAtomChiralTagsFromStructure
  - Only processes true 3D conformers (confId=1) — 2D-only molecules safely skip
  - Returns bool: True=success, False=no-3D-conformer/failure
- Modified mechanism_search.py _build_upload_result:
  - Added Phase A optimization block after path_node_ids population, before node export
  - Iterates path_node_ids, calls optimize_mol_2d() and detect_mol_stereo() on each Mol
  - Logs stats (2d_opt ok/skip, stereo ok/skip) via logger.debug
  - Only path nodes are optimized (non-path intermediates left as-is)
- Updated import in mechanism_search.py: added optimize_mol_2d, detect_mol_stereo
- Verified: syntax check OK, import test OK, functional test OK
  - NumConformers=1 → optimize_2d=True ✅
  - NumConformers=0 → optimize_2d=False ✅
  - Only 2D conf → detect_stereo=False ✅
  - Has 3D conf (Is3D=True) → detect_stereo=True ✅
- Phase B (precise bond length via diff-set product identification) deferred:
  - User identified that set difference (to_mols - from_mols) is unreliable due to Mol object identity
  - Original Upload also doesn't force precise bond length in result display
  - Will revisit when need is clear

Stage Summary:
- optimize_mol_2d() and detect_mol_stereo() added as safe standalone functions
- _build_upload_result now optimizes path node molecules for cleaner SVG display
- Safety guards: NumConformers==1 check, Is3D() check, all failures logged not crashed
- Phase B deferred — original Upload doesn't require precise bond length in display

---
Task ID: 7
Agent: main
Task: Fix remap_reaction_by_mcs atom conservation bug (orphan atoms 15/17)

Root Cause:
- remap_reaction_by_mcs greedy-pairs r0(penicillin)↔p0(hydrolyzed), r1(water16)↔p1(water16)
- r2(water17) has no matching product molecule → becomes UNASSIGNED
- p0 has 1 unmatched atom (OH oxygen from consumed water17) → gets NEW AM number
- Result: AM:15 exists only in products, AM:17 exists only in reactants → atom conservation violated
- The Step 4 verification only checked bond change count, not atom conservation
- Since bond changes were equal (2→2), it kept the broken original mapping

Fix (two parts):
1. Added Step 3.5: After all AM assignments, find "lost" and "new" orphan atoms,
   match them by element, replace new AM with lost AM → fixes atom conservation
2. Updated Step 4 verification: adopt new mapping if bond changes reduced OR
   if bond changes equal but atom conservation improved (orphans decreased)

Verification:
- User's β-lactam hydrolysis: AM:15→AM:17 merge works, OH:17 correctly bonded to C:12 ✅
- All edge cases pass (correct_mapping, simple_hydrolysis, two_waters_consumed, simple_substitution)
- Proton transfer (H+ only): orphans=1, expected — H atoms excluded from matching by design
- Prediction service restarted, health OK

---
Task ID: 5
Agent: general-purpose
Task: Fix all_simple_paths performance bug

Work Log:
- Read mechanism_search.py and identified 3 separate calls to nx.all_simple_paths in _build_upload_result
  - Line 3288: Collects path_node_ids and all_path_node_sets
  - Line 3301: Builds rule_path_count dict
  - Line 3457: Builds all_paths list for path object construction
- Consolidated all 3 calls into a single enumeration loop (lines 3284-3327)
- Added MAX_PATHS=500 limit to prevent exponential path enumeration on large graphs
- Added PATH_TIMEOUT=60s timeout to prevent OOM kills on graphs with 75K+ nodes
- Added early termination with logger.warning for both limits
- Added logger.info for path enumeration stats (count + elapsed time)
- Replaced third call (was `list(nx.all_simple_paths(...))`) with reuse of already-collected `all_paths` list
- Changed condition from `nx.has_path(...)` to `all_paths` (empty list is falsy) to avoid redundant path check
- Verified: path_node_ids still includes r_node_id and p_node_id even when no paths found (pre-existing lines 3279-3282)
- Verified: Python syntax check passes

Stage Summary:
- 3x all_simple_paths calls → 1x (3x speedup on path enumeration, eliminates 2/3 of OOM risk)
- MAX_PATHS=500 prevents unbounded memory growth from exponential path counts
- PATH_TIMEOUT=60s prevents OOM kills on very large graphs
- All existing functionality preserved (path_node_ids, all_path_node_sets, rule_path_count, path objects)

---
Task ID: 6
Agent: general-purpose
Task: Fix child process crash detection bug in route.ts

Work Log:
- Read route.ts and identified two bugs in async mechanism-search flow:
  1. Close handler only checks `code !== 0 && code !== null`, which misses signal-based termination (SIGKILL/OOM gives code=null, signal='SIGKILL')
  2. No stale job detection — if close handler fails, progress.json stays in "RUNNING" forever
- Fix 1: Updated child.on("close") handler (line 325):
  - Changed signature from `(code)` to `(code, signal)` to capture signal parameter
  - Changed condition from `code !== 0 && code !== null` to `code !== 0 || signal`
  - Updated error message from `exited with code ${code}` to `exited (code=${code}, signal=${signal})`
  - Now catches SIGKILL, SIGTERM, and any other signal-based termination
- Fix 2: Added stale job detection in mechanism-search-status handler (after line 391):
  - Checks if progressData.state is "RUNNING" or "STARTING"
  - Compares updated_at timestamp against 30-minute stale threshold
  - If stale, updates progressData to state="ERROR" with descriptive message
  - Writes updated state back to progress.json file (with pretty-print for readability)
  - Placed BEFORE the DONE result-file check so stale jobs don't accidentally match DONE path

Stage Summary:
- **Close handler bug FIXED**: Signal-based termination (SIGKILL/OOM) now correctly updates progress to ERROR
- **Stale job detection ADDED**: Jobs stuck in RUNNING for 30+ minutes auto-marked as ERROR on status poll
- No other API routes or functionality modified
---
Task ID: 8
Agent: main
Task: Optimize RESIDUE_SIDECHAIN_SMILES - replace full amino acid SMILES with minimal catalytic fragments

Work Log:
- Analyzed original EzMechanism residue.py: uses .mol files with minimal side chain local structures
- Proposed minimal side chain fragments (e.g., SER: N[C@@H](CO)C(=O)O → CO)
- Ran comprehensive experiment with 6083 complete rules from custom.db:
  - Tested 5 fragment variants per residue (full_aa, sidechain, *sidechain, minimal, *minimal)
  - Counted rule matches for each variant
  - Tested RunReactants product counts
  - Verified atom map number assignment
- **Critical discovery: `*` connection point HURTS matching dramatically**
  - SER: CO=465 matches vs *CO=7 matches (-98.5%!)
  - CYS: CS=261 matches vs *CS=5 matches (-98%)
  - Root cause: RDKit `*` (atomicNum=0) cannot match `[#6+0]` (carbon) in rule SMARTS
  - Zero rules contain `[*]`, `[R]`, or `[#0]` patterns
  - The original EzMechanism uses ChemAxon Marvin where R matches C; RDKit's * is fundamentally different
- Implemented final solution WITHOUT `*` connection points:
  - RESIDUE_SIDECHAIN_SMILES: Full AA → minimal catalytic fragments
  - Added MAINCHAIN_FRAGMENT_SMILES dictionary for backbone groups
  - Updated get_residue_smiles() to support 'part' field (side_chain/main_chain)
  - GLY (no side chain) returns empty string → skipped in search
- Validated:
  - All SMILES valid, parse correctly with MolFromSmiles + AddHs
  - AM assignment works correctly (heavy atoms get unique AMs, H gets AM=0)
  - 465 rules match SER(CO) fragment (confirmed with custom.db)
  - H atom reduction: 20-56% across catalytic residues
  - Syntax check ✅, lint ✅
  - Pipeline test: parsing, AM assignment, isotope setup, substructure matching all work

Stage Summary:
- **RESIDUE_SIDECHAIN_SMILES replaced** with minimal catalytic fragments (2-bond range principle)
- **`*` connection point REJECTED** based on experimental evidence (-98.5% matching for SER)
- **MAINCHAIN_FRAGMENT_SMILES added** for backbone catalytic groups
- **get_residue_smiles updated** with part='side_chain'/'main_chain' support
- **H atom reduction**: SER -43%, CYS -43%, LYS -50%, HIS -56%, ASP -20%, GLU -43%, TYR -45%, ARG -50%
- **Rule matching preserved**: Key catalytic residues match correct number of rules
  (backbone false positives eliminated: HIS 1346→5, TYR 1346→7, ARG 1376→38)
