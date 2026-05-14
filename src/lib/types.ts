// ============================================================================
// TypeScript Interfaces for M-CSA Mechanism Predictor
// ============================================================================

export interface RuleMatch {
  rule_id: string | number;
  rule_name?: string;
  mcsa_id?: number | null;
  step_id?: number | null;
  mechanism_id?: number | null;
  reaction_smarts: string;
  reactant_smarts: string;
  product_smarts?: string;
  products: string[];
  match_score: number;
  category?: string;
  enzyme?: string;
  source?: string;
}

export interface PredictionStep {
  step: number;
  substrates: string[];
  products: string[];
  rules_applied: {
    rule_id: string | number;
    rule_name?: string;
    mcsa_id?: number | null;
    reaction_smarts: string;
    substrate: string;
    product: string;
    match_score: number;
    category?: string;
    enzyme?: string;
  }[];
  rules_count: number;
}

export interface MatchRulesResponse {
  smiles: string;
  matches: RuleMatch[];
  total_rules: number;
  matches_count: number;
  elapsed_seconds: number;
}

export interface PredictResponse {
  smiles: string;
  steps: PredictionStep[];
  max_steps: number;
  total_steps: number;
  total_rules_checked: number;
  total_rules_available: number;
  elapsed_seconds: number;
  unique_products: number;
}

export interface MolInfo {
  valid: boolean;
  smiles?: string;
  canonical_smiles?: string;
  formula?: string;
  num_atoms?: number;
  num_heavy_atoms?: number;
  molecular_weight?: number;
  num_rotatable_bonds?: number;
  num_h_bond_donors?: number;
  num_h_bond_acceptors?: number;
  tpsa?: number;
  logp?: number;
  ring_count?: number;
  aromatic_rings?: number;
  error?: string;
}

export interface HealthInfo {
  status: string;
  uptime: string;
  rules: {
    total_rules: number;
    mcsa_rules_test: number;
    mcsa_rules_natmet: number;
    builtin_rules: number;
  };
}

export interface BondLengthInfo {
  estimated_bond_length: number;
  bond_length_feasible: boolean;
  bond_length_favorable: boolean;
  closest_non_bonded_distance?: number;
  total_non_bonded_pairs?: number;
}

export interface GeometryStats {
  coords_loaded: number;
  max_bond_length: number;
  favorable_bond_length: number;
  bond_length_bonus: number;
  bond_length_penalty: number;
}

export interface ArrowInfo {
  source_atoms: number[];
  target_atoms: number[];
  electrons: number;
  description?: string;
}

export interface SearchDebugInput {
  reactants_smiles_raw: string;
  products_smiles_raw: string;
  reactant_mols_parsed: string[];
  product_mols_parsed: string[];
  residue_smiles: string[];
  total_rules: number;
  forward_start_key: boolean;
  backward_start_key: boolean;
  react_key_eq_prod_key: boolean | string;
  forward_configs_explored: number;
  backward_configs_explored: number;
}

export interface SearchProgress {
  state: "STARTING" | "RUNNING" | "DONE" | "ERROR" | "CANCELLED" | "UNKNOWN";
  explored_nodes: number;
  total_nodes: number;
  current_iteration: number;
  elapsed_seconds: number;
  created_at?: string;
  updated_at?: string;
  error?: string;
  result_file?: string;
  result?: MechanismSearchResult;
}

export interface MechanismSearchResult {
  reactants: string;
  products: string;
  debug_input?: SearchDebugInput;
  paths: {
    path_id?: number;
    path_type?: string;   // "score_optimal" | "step_minimal"
    path_label?: string;  // "Score-Optimal Path" | "Shortest-Step Path"
    total_steps?: number;
    num_steps?: number;
    total_score?: number;
    merge_similarity?: number;
    steps: {
      step_num?: number;
      from_config?: { molecules: string[]; source: string };
      to_config?: { molecules: string[]; source: string };
      rule_id?: string | number;
      rule_name?: string;
      reaction_smarts?: string;
      direction?: string;
      score?: number;
      max_bond_length?: number;
      reaction_type?: string;
      reaction_type_confidence?: number;
      arrows?: ArrowInfo[];
      arrow_svg?: string;
    }[];
  }[];
  graph: {
    nodes: { id: string; label: string; smiles_label?: string; molecules: string[]; source: string; depth: number }[];
    edges: { id: string; source: string; target: string; rule_id?: string; direction?: string; label?: string; score?: number; reaction_smarts?: string; rule_status?: string; max_bond_length?: number }[];
  };
  stats: {
    search_mode?: string;
    forward_explored?: number;
    backward_explored?: number;
    forward_configs?: number;
    backward_configs?: number;
    meeting_points?: number;
    total_rules?: number;
    paths_found?: number;
    total_nodes?: number;
    total_edges?: number;
    total_expanded?: number;
    max_depth_reached?: number;
    search_time?: number;
  };
  geometry_stats?: GeometryStats | null;
  load_time?: number;
  search_time?: number;
  total_time?: number;
  flask_elapsed?: number;
}

// ---- Residue Info Types ----
export interface ResidueRole {
  role: string;
  counter_role?: string;
  confidence: number;
  detection_method: string;
  description: string;
  suggested_residue: string;
  suggested_position?: string;
}

export interface ResidueSummary {
  total_residues: number;
  roles_breakdown: Record<string, number>;
  catalytic_triad?: string[];
  catalytic_dyad?: string[];
  metal_ion?: string;
}

export interface ResidueInfoResponse {
  residues: ResidueRole[];
  residue_summary: ResidueSummary;
}

// ---- Database Links Types ----
export interface DbLink {
  url: string;
  label: string;
  available: boolean;
}

export interface RuleLinksResponse {
  mcsa?: DbLink;
  pdbe?: DbLink;
  pdb?: DbLink;
  uniprot?: DbLink;
  kegg?: DbLink;
  brenda?: DbLink;
  expasy?: DbLink;
  [key: string]: unknown;
  integrity?: { score: number; max: number };
}

export interface MoleculeLinksResponse {
  pubchem?: DbLink;
  chemspider?: DbLink;
  chebi_search?: DbLink;
  [key: string]: unknown;
}

// ---- Arrow Editor Types ----
export interface ArrowAtom {
  index: number;
  symbol: string;
  atomic_num: number;
  is_aromatic: boolean;
  degree: number;
  formal_charge: number;
}

export interface ArrowBond {
  source: number;
  target: number;
  order: number;
  type: string;
}

export interface ArrowAtomsResponse {
  smiles: string;
  canonical_smiles: string;
  num_atoms: number;
  atoms: ArrowAtom[];
  bonds: ArrowBond[];
  svg?: string;
}

export interface UserArrow {
  source_atoms: number[];
  target_atoms: number[];
  electrons: number;
}

// ---- PDB Types ----
export interface PdbChain {
  chain_id: string;
  num_residues: number;
  sequence: string;
  residues: PdbResidue[];
  uniprot_accession?: string;
}

export interface PdbResidue {
  res_name: string;
  res_num: number;
  insertion: string;
  chain: string;
  seq_pos: number;
  seq_db_offset?: number;
}

export interface PdbLigand {
  res_name: string;
  res_num: number;
  chain: string;
  num_atoms: number;
  formula?: string;
  similarity?: number;
  x?: number;
  y?: number;
  z?: number;
  is_water?: boolean;
}

export interface PdbInfo {
  pdb_id: string;
  title: string;
  resolution?: number;
  deposition_date?: string;
  num_chains: number;
  chains: PdbChain[];
  ligands: PdbLigand[];
  water_count?: number;
  total_atoms?: number;
}

export interface SelectedResidue {
  res_name: string;
  res_num: number;
  chain: string;
  part: "side_chain" | "main_chain";
}

export interface SubstrateCofactor {
  id: string;
  name: string;
  smiles: string;
  type: "substrate" | "cofactor" | "water";
  mapped_ligand?: string;
  isEditing: boolean;
}
