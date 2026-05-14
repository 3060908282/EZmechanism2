export const API_BASE = "/api";

export const STEPS = [
  { id: 1, label: "PDB Structure", icon: "📄" },
  { id: 2, label: "Catalytic Residues", icon: "🧬" },
  { id: 3, label: "Substrates", icon: "⚗️" },
  { id: 4, label: "Reaction", icon: "🔄" },
  { id: 5, label: "Run Search", icon: "🚀" },
] as const;

export const QUICK_SUBSTRATES: Array<{ name: string; smiles: string; icon: string; desc: string }> = [
  { name: "Ethanol", smiles: "CCO", icon: "🍷", desc: "Simple alcohol" },
  { name: "Acetate", smiles: "CC(=O)O", icon: "🧪", desc: "Simple acid" },
  { name: "Glucose", smiles: "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", icon: "🍬", desc: "Sugar monomer" },
  { name: "Benzene", smiles: "c1ccccc1", icon: "⭕", desc: "Aromatic ring" },
];

/**
 * Residue sidechain SMILES — mirrors RESIDUE_SIDECHAIN_SMILES in mechanism_search.py
 * Used by DebugDialog and residue validation
 */
export const RESIDUE_SMILES_MAP: Record<string, string> = {
  SER: 'N[C@@H](CO)C(=O)O',
  CYS: 'N[C@@H](CS)C(=O)O',
  LYS: 'N[C@@H](CCCCN)C(=O)O',
  ARG: 'N[C@@H](CCCNC(N)=N)C(=O)O',
  HIS: 'N[C@@H](CC1=CNC=N1)C(=O)O',
  ASP: 'N[C@@H](CC(=O)[O-])C(=O)[O-]',
  GLU: 'N[C@@H](CCC(=O)[O-])C(=O)[O-]',
  ASN: 'N[C@@H](CC(=O)N)C(=O)O',
  GLN: 'N[C@@H](CCC(=O)N)C(=O)O',
  TYR: 'N[C@@H](CC1=CC=C(O)C=C1)C(=O)O',
  THR: 'N[C@@H](C(C)O)C(=O)O',
  MET: 'N[C@@H](CCSC)C(=O)O',
  TRP: 'N[C@@H](CC1=CNC2=C1C=CC=C2)C(=O)O',
  PHE: 'N[C@@H](CC1=CC=CC=C1)C(=O)O',
  LEU: 'N[C@@H](CC(C)C)C(=O)O',
  ILE: 'N[C@@H](C(C)CC)C(=O)O',
  VAL: 'N[C@@H](C(C)C)C(=O)O',
  PRO: 'N[C@@H]1CCCC(=O)O1',
  ALA: 'N[C@@H](C)C(=O)O',
  GLY: 'NCC(=O)O',
};

/** Valid standard amino acid residue names — used for residue input validation */
export const VALID_RESIDUE_NAMES: Set<string> = new Set(Object.keys(RESIDUE_SMILES_MAP));
