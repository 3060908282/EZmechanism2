#!/usr/bin/env python3
"""Generate SVG representation for a SMILES string.

Usage:
    python3 mol_svg.py <smiles> [--label <atom_alias>] [--width <w>] [--height <h>] [--idx <atom_idx>]

Options:
    --label <text>    Set display alias/label on the first atom (or atom at --idx).
                      Example: --label "Cys145A" will render the atom as "Cys145A"
                      instead of its element symbol. Used for residue-style rendering
                      (e.g., HS—Cys145A in reaction preview).
    --idx <n>         Atom index to apply the label to (default: 0 = first atom).
    --width <w>       SVG width in pixels (default: 300).
    --height <h>      SVG height in pixels (default: 300).

Output: JSON to stdout
"""

import json
import sys

from rdkit import RDLogger, Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D, MolDrawOptions
RDLogger.DisableLog('rdApp.*')

from shared import validate_smiles


def parse_args(argv):
    """Parse command-line arguments."""
    args = {
        'smiles': None,
        'label': None,
        'idx': 0,
        'width': 600,
        'height': 600,
    }
    i = 1
    while i < len(argv):
        if argv[i] == '--label' and i + 1 < len(argv):
            args['label'] = argv[i + 1]
            i += 2
        elif argv[i] == '--idx' and i + 1 < len(argv):
            try:
                args['idx'] = int(argv[i + 1])
            except ValueError:
                print(json.dumps({'error': 'Invalid --idx value'}))
                sys.exit(1)
            i += 2
        elif argv[i] == '--width' and i + 1 < len(argv):
            try:
                args['width'] = int(argv[i + 1])
            except ValueError:
                print(json.dumps({'error': 'Invalid --width value'}))
                sys.exit(1)
            i += 2
        elif argv[i] == '--height' and i + 1 < len(argv):
            try:
                args['height'] = int(argv[i + 1])
            except ValueError:
                print(json.dumps({'error': 'Invalid --height value'}))
                sys.exit(1)
            i += 2
        elif not args['smiles']:
            args['smiles'] = argv[i]
            i += 1
        else:
            i += 1
    return args


def main() -> None:
    args = parse_args(sys.argv)

    if not args['smiles']:
        print(json.dumps({'error': 'Usage: python3 mol_svg.py <smiles> [--label <text>] [--idx <n>]'}))
        sys.exit(1)

    smiles = args['smiles']
    validate_smiles(smiles)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(json.dumps({'error': 'Invalid SMILES'}))
        sys.exit(1)

    AllChem.Compute2DCoords(mol)

    w = args['width']
    h = args['height']

    # Set up draw options — minimal padding so molecule fills the canvas
    opts = MolDrawOptions()
    opts.padding = 0.02

    # Apply atom alias label using atomLabels map (rendered as custom text)
    if args['label']:
        idx = args['idx']
        if idx < mol.GetNumAtoms():
            opts.atomLabels[idx] = args['label']
            # Dynamic font scaling based on label length:
            #   >12 chars → much larger (e.g. "Glu1234AB" style)
            #   >8 chars  → larger (e.g. "Cys145A" style)
            label_len = len(args['label'])
            if label_len > 12:
                opts.atomLabelFontSize = 22
            elif label_len > 8:
                opts.atomLabelFontSize = 17

    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    drawer.SetDrawOptions(opts)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Force SVG to fill container without letterboxing
    svg = svg.replace('preserveAspectRatio="xMidYMid meet"', '')
    svg = svg.replace('preserveAspectRatio="xMidYMid"', '')
    if 'preserveAspectRatio' not in svg:
        svg = svg.replace('<svg ', '<svg preserveAspectRatio="xMidYMid slice" ', 1)

    result = {
        'smiles': smiles,
        'svg': svg,
        'width': w,
        'height': h,
    }
    if args['label']:
        result['label'] = args['label']
    print(json.dumps(result))


if __name__ == '__main__':
    main()
