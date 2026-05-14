/**
 * PDB File Parser — parses ATOM/HETATM/DBREF/DBREF1/DBREF2 records
 * to extract chains, residues, ligands, and UniProt sequence mappings.
 *
 * PDB format reference: https://www.wwpdb.org/documentation/file-format
 */

// Standard 20 amino acids: 3-letter → 1-letter
const AA_3TO1: Record<string, string> = {
  ALA: "A", ARG: "R", ASN: "N", ASP: "D", CYS: "C",
  GLN: "Q", GLU: "E", GLY: "G", HIS: "H", ILE: "I",
  LEU: "L", LYS: "K", MET: "M", PHE: "F", PRO: "P",
  SER: "S", THR: "T", TRP: "W", TYR: "Y", VAL: "V",
  // Common modified / alternative residue names
  MSE: "M", SEC: "U", PYL: "O", ASX: "B", GLX: "Z",
  CSO: "C", HIP: "H", HIE: "H", HID: "H", HSD: "H",
  HSE: "H", HSP: "H", LYZ: "K", CYX: "C", CYM: "C",
};

const STANDARD_AAS = new Set(Object.keys(AA_3TO1));

/** A parsed residue key for deduplication */
interface ResidueKey {
  res_name: string;
  res_num: number;
  insertion: string;
  chain: string;
}

/** DBREF mapping info for a single chain */
interface DbrefMapping {
  chain: string;
  db_name: string;
  db_accession: string;
  pdb_start: number;
  pdb_end: number;
  ins_start: string;
  ins_end: string;
  db_start: number;
  db_end: number;
}

/** Intermediate residue accumulator */
interface ResidueAccum {
  res_name: string;
  res_num: number;
  insertion: string;
  chain: string;
  atom_count: number;
}

export interface ParsedPdbInfo {
  pdb_id: string;
  title: string;
  resolution?: number;
  deposition_date?: string;
  num_chains: number;
  chains: ParsedPdbChain[];
  ligands: ParsedPdbLigand[];
  water_count: number;
  total_atoms: number;
}

export interface ParsedPdbChain {
  chain_id: string;
  num_residues: number;
  sequence: string;
  residues: ParsedPdbResidue[];
  uniprot_accession?: string;
}

export interface ParsedPdbResidue {
  res_name: string;
  res_num: number;
  insertion: string;
  chain: string;
  seq_pos: number;
  seq_db_offset?: number;
}

export interface ParsedPdbLigand {
  res_name: string;
  res_num: number;
  chain: string;
  num_atoms: number;
  formula?: string;
}

/**
 * Parse a complete PDB file text and return structured info.
 */
export function parsePdbText(pdbText: string): ParsedPdbInfo {
  const lines = pdbText.split("\n");

  // --- Extract metadata ---
  let pdbId = "CUSTOM";
  let title = "Parsed PDB Structure";
  let resolution: number | undefined;
  let depositionDate: string | undefined;
  const titleParts: string[] = [];

  for (const line of lines) {
    const rec = line.slice(0, 6).trim();
    if (rec === "HEADER") {
      // HEADER    CLASSIFICATION DATE   IDCODE
      const idcode = line.slice(62, 66).trim();
      if (idcode && idcode.length === 4) pdbId = idcode.toUpperCase();
      const dateStr = line.slice(50, 59).trim();
      if (dateStr) depositionDate = formatDate(dateStr);
    } else if (rec === "TITLE") {
      const t = line.slice(10).trim();
      if (t) titleParts.push(t);
    } else if (rec === "REMARK" && line.slice(7, 10).trim() === "2" && line.includes("RESOLUTION.")) {
      // REMARK   2 RESOLUTION. ... <value> ANGSTROMS
      const resMatch = line.match(/RESOLUTION\.\s+(\d+\.?\d*)/);
      if (resMatch) resolution = parseFloat(resMatch[1]);
    }
  }
  if (titleParts.length > 0) {
    title = titleParts.join(" ").trim();
  }

  // --- Parse DBREF records ---
  const dbrefMap = new Map<string, DbrefMapping>();
  for (const line of lines) {
    const rec = line.slice(0, 6).trim();
    if (rec === "DBREF") {
      // Legacy DBREF (old format)
      // COLUMNS: 8-11 ID, 13 Chain, 15-18 initRes, 19 insInit, 21-24 endRes, 25 insEnd,
      //          27-32 dbName, 48-67 accession, 73-78 dbStart, 80-85 dbEnd
      const chain = line.slice(12, 13).trim();
      if (!chain) continue;
      const pdbStart = parseInt(line.slice(14, 18).trim(), 10);
      const insStart = (line.slice(18, 19) || " ").trim();
      const pdbEnd = parseInt(line.slice(22, 26).trim(), 10);
      const insEnd = (line.slice(26, 27) || " ").trim();
      const dbName = line.slice(26, 32).trim(); // note: overlaps with insEnd, but old format is tricky
      const accession = line.slice(47, 67).trim();
      const dbStart = parseInt(line.slice(72, 78).trim(), 10);
      const dbEnd = parseInt(line.slice(79, 85).trim(), 10);

      if (!isNaN(pdbStart) && !isNaN(pdbEnd) && !isNaN(dbStart) && !isNaN(dbEnd)) {
        // Fix dbName from proper column
        const realDbName = line.slice(26, 32).trim();
        dbrefMap.set(chain, {
          chain,
          db_name: realDbName,
          db_accession: accession,
          pdb_start: pdbStart,
          pdb_end: pdbEnd,
          ins_start: insStart,
          ins_end: insEnd,
          db_start: dbStart,
          db_end: dbEnd,
        });
      }
    } else if (rec === "DBREF1") {
      // New format DBREF1
      // COLUMNS: 8-11 ID, 13 Chain, 15-18 initRes, 19 insInit, 21-24 endRes, 25 insEnd,
      //          27-32 dbName, 48-67 accession
      const chain = line.slice(12, 13).trim();
      if (!chain) continue;
      const pdbStart = parseInt(line.slice(14, 18).trim(), 10);
      const insStart = (line.slice(18, 19) || " ").trim();
      const pdbEnd = parseInt(line.slice(22, 26).trim(), 10);
      const insEnd = (line.slice(26, 27) || " ").trim();
      const dbName = line.slice(26, 32).trim();
      const accession = line.slice(47, 67).trim();

      const existing = dbrefMap.get(chain) || {
        chain, db_name: dbName, db_accession: accession,
        pdb_start: pdbStart, pdb_end: pdbEnd,
        ins_start: insStart, ins_end: insEnd,
        db_start: 0, db_end: 0,
      };
      existing.pdb_start = pdbStart;
      existing.pdb_end = pdbEnd;
      existing.ins_start = insStart;
      existing.ins_end = insEnd;
      existing.db_name = dbName;
      existing.db_accession = accession;
      dbrefMap.set(chain, existing);
    } else if (rec === "DBREF2") {
      // New format DBREF2
      // COLUMNS: 8-11 ID, 13 Chain, 19-40 accession, 46-55 dbStart, 58-67 dbEnd
      const chain = line.slice(12, 13).trim();
      if (!chain) continue;
      const accession = line.slice(18, 40).trim();
      const dbStart = parseInt(line.slice(45, 55).trim(), 10);
      const dbEnd = parseInt(line.slice(57, 67).trim(), 10);

      const existing = dbrefMap.get(chain);
      if (existing) {
        if (accession) existing.db_accession = accession;
        if (!isNaN(dbStart)) existing.db_start = dbStart;
        if (!isNaN(dbEnd)) existing.db_end = dbEnd;
        dbrefMap.set(chain, existing);
      }
    }
  }

  // --- Parse ATOM & HETATM records ---
  const residueMap = new Map<string, ResidueAccum>();
  const ligandMap = new Map<string, { res_name: string; res_num: number; chain: string; atom_count: number }>();
  let totalAtoms = 0;
  let waterCount = 0;

  // Track water residues to exclude from ligands
  const waterResKeys = new Set<string>();

  for (const line of lines) {
    const rec = line.slice(0, 6).trim();
    if (rec !== "ATOM" && rec !== "HETATM") continue;

    totalAtoms++;

    // Parse columns
    const atomName = line.slice(12, 16).trim();
    const resName = line.slice(17, 20).trim();
    const chain = line.slice(21, 22).trim();
    const resNumStr = line.slice(22, 26).trim();
    const insertion = (line.slice(26, 27) || " ").trim();
    const element = line.slice(76, 78).trim();

    // Skip hydrogen atoms for counting (but still track residue existence)
    const isHydrogen = element === "H" || element === "D" || element === "T" ||
      (element === "" && atomName.startsWith("H"));

    const resNum = parseInt(resNumStr, 10);
    if (isNaN(resNum)) continue;

    // Skip HOH / WAT residues (water)
    if (resName === "HOH" || resName === "WAT" || resName === "H2O") {
      const wkey = `${chain}:${resName}:${resNum}:${insertion}`;
      if (!waterResKeys.has(wkey)) {
        waterResKeys.add(wkey);
        waterCount++;
      }
      continue;
    }

    const isStandardAA = STANDARD_AAS.has(resName);
    const isAtom = rec === "ATOM";

    if (isStandardAA && isAtom) {
      // Standard amino acid residue in ATOM record
      const key = `${chain}:${resName}:${resNum}:${insertion}`;
      const existing = residueMap.get(key);
      if (existing && !isHydrogen) {
        existing.atom_count++;
      } else if (!existing) {
        residueMap.set(key, {
          res_name: resName,
          res_num: resNum,
          insertion,
          chain,
          atom_count: isHydrogen ? 0 : 1,
        });
      }
    } else if (rec === "HETATM" && !isStandardAA) {
      // Heteroatom ligand (non-standard, non-water)
      const key = `${chain}:${resName}:${resNum}:${insertion}`;
      const existing = ligandMap.get(key);
      if (existing && !isHydrogen) {
        existing.atom_count++;
      } else if (!existing) {
        ligandMap.set(key, {
          res_name: resName,
          res_num: resNum,
          chain,
          atom_count: isHydrogen ? 0 : 1,
        });
      }
    }
  }

  // --- Build chains from residues ---
  // Group residues by chain, ordered by res_num
  const chainResidues = new Map<string, ResidueAccum[]>();
  for (const residue of residueMap.values()) {
    const list = chainResidues.get(residue.chain) || [];
    list.push(residue);
    chainResidues.set(residue.chain, list);
  }

  // Sort each chain's residues by res_num, then by insertion code
  for (const list of chainResidues.values()) {
    list.sort((a, b) => {
      if (a.res_num !== b.res_num) return a.res_num - b.res_num;
      return a.insertion.localeCompare(b.insertion);
    });
  }

  // --- Build ParsedPdbChain objects ---
  const chains: ParsedPdbChain[] = [];
  for (const [chainId, residues] of chainResidues) {
    const dbref = dbrefMap.get(chainId);

    const sequenceChars: string[] = [];
    const parsedResidues: ParsedPdbResidue[] = [];

    for (const r of residues) {
      // Compute UniProt sequence position from DBREF mapping
      let seqPos = 0;
      let seqDbOffset: number | undefined;

      if (dbref) {
        // offset = db_start - pdb_start
        const offset = dbref.db_start - dbref.pdb_start;
        seqPos = r.res_num + offset;
        // Check if this residue is within the mapped range
        if (r.res_num >= dbref.pdb_start && r.res_num <= dbref.pdb_end) {
          seqDbOffset = offset;
        } else {
          // Residue is outside DBREF mapping range
          seqPos = 0;
          seqDbOffset = undefined;
        }
      } else {
        // No DBREF mapping — seq_pos defaults to 0 (shown as "—" in UI)
        seqPos = 0;
        seqDbOffset = undefined;
      }

      const oneLetter = AA_3TO1[r.res_name] || "X";
      sequenceChars.push(oneLetter);

      parsedResidues.push({
        res_name: r.res_name,
        res_num: r.res_num,
        insertion: r.insertion,
        chain: r.chain,
        seq_pos: seqPos,
        seq_db_offset: seqDbOffset,
      });
    }

    chains.push({
      chain_id: chainId,
      num_residues: residues.length,
      sequence: sequenceChars.join(""),
      residues: parsedResidues,
      uniprot_accession: dbref?.db_name === "UNP" ? dbref.db_accession : undefined,
    });
  }

  // Sort chains alphabetically
  chains.sort((a, b) => a.chain_id.localeCompare(b.chain_id));

  // --- Build ligands ---
  const ligands: ParsedPdbLigand[] = [];
  for (const lig of ligandMap.values()) {
    ligands.push({
      res_name: lig.res_name,
      res_num: lig.res_num,
      chain: lig.chain,
      num_atoms: lig.atom_count,
    });
  }
  // Sort ligands by chain, then res_num
  ligands.sort((a, b) => {
    if (a.chain !== b.chain) return a.chain.localeCompare(b.chain);
    return a.res_num - b.res_num;
  });

  return {
    pdb_id: pdbId,
    title,
    resolution,
    deposition_date: depositionDate,
    num_chains: chains.length,
    chains,
    ligands,
    water_count: waterCount,
    total_atoms: totalAtoms,
  };
}

/**
 * Convert parsed PDB info to the format expected by the frontend (PdbInfo).
 */
export function toFrontendPdbInfo(parsed: ParsedPdbInfo) {
  return {
    pdb_id: parsed.pdb_id,
    title: parsed.title,
    resolution: parsed.resolution,
    deposition_date: parsed.deposition_date,
    num_chains: parsed.num_chains,
    chains: parsed.chains.map((c) => ({
      chain_id: c.chain_id,
      num_residues: c.num_residues,
      sequence: c.sequence,
      residues: c.residues.map((r) => ({
        res_name: r.res_name,
        res_num: r.res_num,
        insertion: r.insertion,
        chain: r.chain,
        seq_pos: r.seq_pos,
        seq_db_offset: r.seq_db_offset,
      })),
      uniprot_accession: c.uniprot_accession,
    })),
    ligands: parsed.ligands.map((l) => ({
      res_name: l.res_name,
      res_num: l.res_num,
      chain: l.chain,
      num_atoms: l.num_atoms,
    })),
    water_count: parsed.water_count,
    total_atoms: parsed.total_atoms,
  };
}

/**
 * Try to parse a date string from PDB HEADER into ISO format.
 */
function formatDate(dateStr: string): string | undefined {
  // Format: DD-MMM-YY (e.g., "15-JUN-23")
  const months: Record<string, string> = {
    JAN: "01", FEB: "02", MAR: "03", APR: "04", MAY: "05", JUN: "06",
    JUL: "07", AUG: "08", SEP: "09", OCT: "10", NOV: "11", DEC: "12",
  };
  const match = dateStr.match(/^(\d{1,2})-(\w{3})-(\d{2,4})$/);
  if (!match) return dateStr;
  const day = match[1].padStart(2, "0");
  const month = months[match[2].toUpperCase()];
  let year = match[3];
  if (year.length === 2) year = year >= "25" ? "19" + year : "20" + year;
  if (!month) return dateStr;
  return `${year}-${month}-${day}`;
}
