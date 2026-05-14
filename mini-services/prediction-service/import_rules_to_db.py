#!/usr/bin/env python3
"""
One-time import script: xlsx reaction rules → SQLite database.
Usage: python import_rules_to_db.py
"""
import os
import sys
import sqlite3
import openpyxl

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'db', 'custom.db')
RULES_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'rules_nat_met.xlsx')
RULES_FILE_TEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 'rules_test.xlsx')


def import_xlsx(filepath: str, source_tag: str, conn: sqlite3.Connection) -> int:
    """Import rules from xlsx into SQLite."""
    if not os.path.exists(filepath):
        print(f"  File not found: {filepath}")
        return 0

    wb = openpyxl.load_workbook(filepath, read_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    wb.close()

    imported = 0
    skipped = 0
    for row in rows:
        if not row or len(row) < 9:
            continue
        rule_id = str(row[3]).strip() if row[3] else ""
        reaction_smarts = str(row[6]).strip() if row[6] else ""
        if not rule_id or not reaction_smarts:
            skipped += 1
            continue

        # Split reaction_smarts into reactant/product sides
        if ">>" in reaction_smarts:
            parts = reaction_smarts.split(">>", 1)
            reactant_smarts = parts[0].strip()
            product_smarts = parts[1].strip()
        else:
            reactant_smarts = reaction_smarts
            product_smarts = ""

        try:
            conn.execute(
                """INSERT OR IGNORE INTO ReactionRule 
                   (ruleId, mcsaId, mechanismId, stepId, isReversed,
                    ruleCompleteInStep, reactionSmarts, reactantSmarts, productSmarts,
                    radicalInStep, ruleArrows, sourceTag)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    rule_id,
                    str(row[0]).strip() if row[0] else "",
                    str(row[1]).strip() if row[1] else "",
                    str(row[2]).strip() if row[2] else "",
                    bool(row[4]),
                    str(row[5]).strip().lower() in ('true', '1') if row[5] else False,
                    reaction_smarts,
                    reactant_smarts,
                    product_smarts,
                    str(row[7]).strip().lower() in ('true', '1') if row[7] else False,
                    str(row[8]).strip() if row[8] else "",
                    source_tag,
                ),
            )
            imported += 1
        except Exception as e:
            skipped += 1
            if skipped <= 3:
                print(f"  Skip rule {rule_id}: {e}")

    return imported


def main():
    print("=" * 50)
    print("  Import Reaction Rules to SQLite")
    print("=" * 50)
    print(f"  DB: {DB_PATH}")
    print()

    if not os.path.exists(DB_PATH):
        print(f"ERROR: Database not found at {DB_PATH}")
        print("  Run 'bun run db:push' first from the project root.")
        sys.exit(1)

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")

    # Clear existing rules
    conn.execute("DELETE FROM ReactionRule")
    print("  Cleared existing rules.")

    # Import primary rules
    print(f"\n  Importing from: {os.path.basename(RULES_FILE)}")
    n1 = import_xlsx(RULES_FILE, "mcsa", conn)
    print(f"  Imported: {n1} rules")

    # Import test rules
    print(f"\n  Importing from: {os.path.basename(RULES_FILE_TEST)}")
    n2 = import_xlsx(RULES_FILE_TEST, "mcsa-test", conn)
    print(f"  Imported: {n2} rules")

    conn.commit()

    # Stats
    total = conn.execute("SELECT COUNT(*) FROM ReactionRule").fetchone()[0]
    mcsa = conn.execute("SELECT COUNT(*) FROM ReactionRule WHERE sourceTag='mcsa'").fetchone()[0]
    test = conn.execute("SELECT COUNT(*) FROM ReactionRule WHERE sourceTag='mcsa-test'").fetchone()[0]

    # DB file size
    db_size = os.path.getsize(DB_PATH) / 1024

    print()
    print(f"  Total: {total} rules (mcsa: {mcsa}, test: {test})")
    print(f"  DB size: {db_size:.1f} KB")
    print()
    print("  Done!")
    print("=" * 50)

    conn.close()


if __name__ == "__main__":
    main()
