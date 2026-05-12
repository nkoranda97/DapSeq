"""
Migrate pipeline_db.csv → pipeline_db.db

Safe to re-run: output_dirs already present in the SQLite database are
skipped, so no rows are duplicated.

Usage (from the pipeline root):
    pixi run python workflow/scripts/migrate_db.py
    pixi run python workflow/scripts/migrate_db.py --csv /other/path/pipeline_db.csv
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

# update_db lives in the same directory; import its helpers
sys.path.insert(0, str(Path(__file__).parent))
import update_db as db_mod


def parse_args():
    root = Path(__file__).parent.parent.parent  # DapSeq/
    parser = argparse.ArgumentParser(description="Migrate pipeline_db.csv to pipeline_db.db")
    parser.add_argument(
        "--csv",
        default=str(root / "pipeline_db.csv"),
        help="Path to the source CSV (default: pipeline_db.csv next to this repo)",
    )
    parser.add_argument(
        "--db",
        default=str(root / "pipeline_db.db"),
        help="Path to the destination SQLite database (default: pipeline_db.db next to this repo)",
    )
    return parser.parse_args()


def read_csv(csv_path):
    """Return rows grouped by output_dir: {output_dir: [tuple, ...]}"""
    path = Path(csv_path)
    if not path.exists():
        sys.exit(f"Error: CSV not found: {csv_path}")

    groups = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = set(db_mod.COLS) - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"Error: CSV is missing expected columns: {missing}")
        for row in reader:
            t = tuple(row.get(c, "") for c in db_mod.COLS)
            groups[row["output_dir"]].append(t)
    return groups


def existing_output_dirs(db_path):
    """Return the set of output_dirs already in the database."""
    import sqlite3
    path = Path(db_path)
    if not path.exists():
        return set()
    con = sqlite3.connect(str(path), timeout=60)
    try:
        rows = con.execute("SELECT DISTINCT output_dir FROM pipeline_runs").fetchall()
        return {r[0] for r in rows}
    except sqlite3.OperationalError:
        return set()
    finally:
        con.close()


def main():
    args = parse_args()

    groups = read_csv(args.csv)
    if not groups:
        print("CSV is empty — nothing to migrate.")
        return

    already_in_db = existing_output_dirs(args.db)
    to_migrate = {od: rows for od, rows in groups.items() if od not in already_in_db}
    skipped = set(groups) - set(to_migrate)

    if skipped:
        print(f"Skipping {len(skipped)} output_dir(s) already in the database:")
        for od in sorted(skipped):
            print(f"  {od}")

    if not to_migrate:
        print("Nothing new to migrate.")
        return

    print(f"Migrating {len(to_migrate)} output_dir(s) → {args.db}")
    total_rows = 0
    for output_dir, rows in sorted(to_migrate.items()):
        db_mod.write_rows(args.db, output_dir, rows)
        print(f"  {output_dir}: {len(rows)} row(s)")
        total_rows += len(rows)

    print(f"Done — {total_rows} row(s) written.")


if __name__ == "__main__":
    main()
