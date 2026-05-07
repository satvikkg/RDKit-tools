"""
load_sdf_to_postgres.py
-----------------------
Loads an SDF file into a PostgreSQL database with:
  - 2D mol column (RDKit MOL type) for substructure/similarity search
  - 3D molblock column (TEXT) preserving original 3D coordinates
  - Morgan fingerprint column (BFP) for similarity search
  - All SD properties as a JSONB column

Usage:
    python load_sdf_to_postgres.py \
        --sdf library.sdf \
        --dsn "postgresql://user:password@localhost:5432/screening" \
        --table compounds \
        --batch-size 1000
"""

import argparse
import json
import logging
import sys
import time

import psycopg2
import psycopg2.extras
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------

DDL = """
CREATE EXTENSION IF NOT EXISTS rdkit;

CREATE TABLE IF NOT EXISTS {table} (
    id          SERIAL PRIMARY KEY,
    mol_name    TEXT,
    smiles      TEXT,
    mol         MOL,            -- 2D RDKit mol (substructure / similarity)
    molblock    TEXT,           -- full molblock preserving 3D coordinates
    mfp2        BFP,            -- Morgan FP radius 2 (Tanimoto similarity)
    properties  JSONB           -- all SD tag key/value pairs
);

CREATE INDEX IF NOT EXISTS {table}_mol_gist  ON {table} USING gist(mol);
CREATE INDEX IF NOT EXISTS {table}_mfp2_gist ON {table} USING gist(mfp2);
CREATE INDEX IF NOT EXISTS {table}_props_gin ON {table} USING gin(properties);
"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def mol_to_2d(mol):
    """Return a copy of mol with 2D coords (removes 3D), for the MOL column."""
    mol2d = Chem.RWMol(mol)
    mol2d.RemoveAllConformers()
    AllChem.Compute2DCoords(mol2d)
    return mol2d


def compute_mfp2(mol):
    """Morgan bit vector FP radius 2, 2048 bits, as hex string for Postgres BFP."""
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    return fp.ToBitString()


def get_properties(mol):
    """Extract all SD properties into a plain dict."""
    props = {}
    for key in mol.GetPropNames():
        try:
            props[key] = mol.GetProp(key)
        except Exception:
            pass
    return props


def iter_mols(sdf_path, remove_hs=False):
    """Yield (mol, molblock) tuples from an SDF, skipping unreadable entries."""
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=remove_hs, sanitize=True)
    for i, mol in enumerate(suppl):
        if mol is None:
            log.warning("Skipping molecule at index %d (parse/sanitize failed)", i)
            continue
        molblock = Chem.MolToMolBlock(mol)   # preserves 3D if present
        yield mol, molblock


def build_row(mol, molblock):
    """Return a dict ready for psycopg2 execute."""
    mol2d = mol_to_2d(mol)
    smiles = Chem.MolToSmiles(mol)
    mfp2_bits = compute_mfp2(mol)

    return {
        "mol_name":   mol.GetProp("_Name") if mol.HasProp("_Name") else None,
        "smiles":     smiles,
        "mol":        smiles,           # Postgres will cast TEXT -> MOL
        "molblock":   molblock,
        "mfp2":       mfp2_bits,        # Postgres will cast TEXT -> BFP
        "properties": json.dumps(get_properties(mol)),
    }


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

INSERT_SQL = """
INSERT INTO {table}
    (mol_name, smiles, mol, molblock, mfp2, properties)
VALUES
    (%(mol_name)s, %(smiles)s,
     mol_from_smiles(%(mol)s::cstring),
     %(molblock)s,
     %(mfp2)s::bfp,
     %(properties)s::jsonb)
ON CONFLICT DO NOTHING;
"""


def load(sdf_path, dsn, table, batch_size):
    conn = psycopg2.connect(dsn)
    conn.autocommit = False
    cur = conn.cursor()

    log.info("Creating schema (if not exists) …")
    cur.execute(DDL.format(table=table))
    conn.commit()

    insert_sql = INSERT_SQL.format(table=table)

    batch = []
    total_ok = 0
    total_skip = 0
    t0 = time.time()

    log.info("Reading %s …", sdf_path)

    for mol, molblock in iter_mols(sdf_path):
        try:
            row = build_row(mol, molblock)
        except Exception as exc:
            log.warning("Skipping molecule — fingerprint/SMILES error: %s", exc)
            total_skip += 1
            continue

        batch.append(row)

        if len(batch) >= batch_size:
            psycopg2.extras.execute_batch(cur, insert_sql, batch, page_size=batch_size)
            conn.commit()
            total_ok += len(batch)
            batch = []
            elapsed = time.time() - t0
            log.info("  Loaded %d compounds  (%.1f/s)", total_ok, total_ok / elapsed)

    # flush remainder
    if batch:
        psycopg2.extras.execute_batch(cur, insert_sql, batch, page_size=batch_size)
        conn.commit()
        total_ok += len(batch)

    cur.close()
    conn.close()

    elapsed = time.time() - t0
    log.info("Done. %d loaded, %d skipped in %.1f s (%.0f cpd/s)",
             total_ok, total_skip, elapsed, total_ok / max(elapsed, 1))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Load SDF into RDKit/Postgres screening DB")
    p.add_argument("--sdf",        required=True,  help="Path to input SDF file")
    p.add_argument("--dsn",        required=True,  help="Postgres DSN, e.g. postgresql://user:pw@host/db")
    p.add_argument("--table",      default="compounds", help="Target table name (default: compounds)")
    p.add_argument("--batch-size", type=int, default=1000, help="Rows per commit batch (default: 1000)")
    p.add_argument("--keep-hs",    action="store_true",   help="Keep explicit hydrogens")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    try:
        load(
            sdf_path=args.sdf,
            dsn=args.dsn,
            table=args.table,
            batch_size=args.batch_size,
        )
    except KeyboardInterrupt:
        log.info("Interrupted.")
        sys.exit(1)
    except Exception as exc:
        log.error("Fatal: %s", exc)
        sys.exit(1)
