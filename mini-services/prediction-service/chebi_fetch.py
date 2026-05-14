#!/usr/bin/env python3
"""
ChEBI molecule fetcher — 通过 OLS4 API + PubChem 获取分子信息。

用法:
    python3 chebi_fetch.py <chebi_id>
    # chebi_id 可以是 "15377" 或 "CHEBI:15377" 格式

输出:
    JSON: { "chebi_id": "...", "name": "...", "smiles": "...", "definition": "...", "formula": "..." }
"""

import json
import sys
import urllib.request
import urllib.parse
import urllib.error


def _url_get(url: str, timeout: int = 10) -> str | None:
    """Simple GET request with timeout."""
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "EzMechanism/1.0",
            "Accept": "application/json",
        })
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8")
    except Exception:
        return None


def fetch_chebi(chebi_num: str) -> dict:
    """Fetch molecule info for a ChEBI ID using OLS4 + PubChem."""
    result = {
        "chebi_id": f"CHEBI:{chebi_num}",
        "name": None,
        "smiles": None,
        "definition": None,
        "formula": None,
        "inchi": None,
        "inchikey": None,
    }

    iri = f"http://purl.obolibrary.org/obo/CHEBI_{chebi_num}"
    iri_encoded = urllib.parse.quote(iri, safe="")

    # Step 1: Query OLS4 general endpoint for all ontology matches
    try:
        url = f"https://www.ebi.ac.uk/ols4/api/terms?iri={iri_encoded}"
        text = _url_get(url, timeout=10)
        if text:
            data = json.loads(text)
            terms = data.get("_embedded", {}).get("terms", [])

            # Find the chebi ontology entry first for name/definition
            for t in terms:
                if t.get("ontology_name") == "chebi":
                    result["name"] = t.get("label")
                    desc = t.get("description", [])
                    if desc:
                        result["definition"] = desc[0]
                    break

            # Also check cido (Chemical Information and Data Ontology) for annotations
            # and any other ontology that has smiles
            if not result["name"] and terms:
                result["name"] = terms[0].get("label")
                desc = terms[0].get("description", [])
                if desc:
                    result["definition"] = desc[0]

            # Look for smiles/inchi/formula across all ontology entries
            for t in terms:
                anno = t.get("annotation", {})
                if not result.get("smiles"):
                    smi = anno.get("smiles", [])
                    if smi and smi[0] and smi[0] != "N/A":
                        result["smiles"] = smi[0]
                if not result.get("inchi"):
                    inch = anno.get("inchi", [])
                    if inch and inch[0] and inch[0] != "N/A":
                        result["inchi"] = inch[0]
                if not result.get("inchikey"):
                    key = anno.get("inchikey", [])
                    if key and key[0] and key[0] != "N/A":
                        result["inchikey"] = key[0]
                if not result.get("formula"):
                    form = anno.get("formula", [])
                    if form and form[0] and form[0] != "N/A":
                        result["formula"] = form[0]
    except Exception as e:
        print(f"Warning: OLS4 lookup failed: {e}", file=sys.stderr)

    # Step 2: Try PubChem for SMILES using the name
    if not result.get("smiles") and result.get("name"):
        try:
            name_encoded = urllib.parse.quote(result["name"])
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name_encoded}/property/IsomericSMILES,CanonicalSMILES,MolecularFormula/JSON"
            text = _url_get(url, timeout=10)
            if text:
                data = json.loads(text)
                props = data.get("PropertyTable", {}).get("Properties", [])
                if props and len(props) > 0:
                    result["smiles"] = props[0].get("IsomericSMILES") or props[0].get("SMILES")
                    if not result.get("formula"):
                        result["formula"] = props[0].get("MolecularFormula")
        except Exception as e:
            print(f"Warning: PubChem lookup failed: {e}", file=sys.stderr)

    # Step 3: Try InChI → SMILES conversion via RDKit
    if not result.get("smiles") and result.get("inchi"):
        try:
            from rdkit import Chem
            mol = Chem.MolFromInchi(result["inchi"])
            if mol:
                result["smiles"] = Chem.MolToSmiles(mol)
        except Exception:
            pass

    # Step 4: Clean up SMILES (remove [H] tags for standard form)
    if result.get("smiles"):
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(result["smiles"])
            if mol:
                result["smiles"] = Chem.MolToSmiles(mol)
        except Exception:
            pass

    return result


def main():
    if len(sys.argv) < 2:
        print(json.dumps({"error": "Usage: chebi_fetch.py <chebi_id>"}))
        sys.exit(1)

    chebi_id = sys.argv[1].strip().replace("CHEBI:", "").replace("chebi:", "")

    if not chebi_id.isdigit():
        print(json.dumps({"error": "Invalid ChEBI ID: must be numeric"}))
        sys.exit(1)

    result = fetch_chebi(chebi_id)
    print(json.dumps(result))


if __name__ == "__main__":
    main()
