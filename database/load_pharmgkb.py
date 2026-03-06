"""
WEEK 2 — PharmGKB Data Loader
================================
Reads relationships.tsv and loads gene-drug relationships
into the pharmgkb_annotations table in genomics_db.

We filter to:
  - Only Gene → Chemical relationships
  - Only 'associated' (not 'not associated' or 'ambiguous')
  - Only rows with clinical or variant annotation evidence

Run from your project root:
  python database/load_pharmgkb.py
"""

import csv
import psycopg2
from pathlib import Path

# ── Database connection ──────────────────────────────────────
DB_CONFIG = {
    "host"     : "localhost",
    "database" : "genomics_db",
    "user"     : "genomics_user",
    "password" : "genomics123",
    "port"     : "5432"
}

# ── File path ────────────────────────────────────────────────
PHARMGKB_FILE = Path("data/raw/relationships.tsv")

# ── Only keep rows with these evidence types ─────────────────
KEEP_EVIDENCE = {"ClinicalAnnotation", "VariantAnnotation"}

# ── Only keep confirmed associations ─────────────────────────
KEEP_ASSOCIATION = {"associated"}

# ── CYP genes — the core pharmacogenomics genes ──────────────
# These get flagged specially in your pipeline
CYP_GENES = {
    "CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
    "CYP2D6", "CYP2E1", "CYP3A4", "CYP3A5", "CYP3A7",
    "DPYD", "TPMT", "UGT1A1", "SLCO1B1", "ABCG2",
    "G6PD", "IFNL3", "HLA-A", "HLA-B"
}


def connect():
    """Create and return a database connection."""
    return psycopg2.connect(**DB_CONFIG)


def clean_value(value):
    """Convert empty or missing values to None."""
    if value in ("", "-", "NA", "N/A", "na"):
        return None
    return value


def determine_phenotype_category(pk, pd, evidence):
    """
    Determine the phenotype category from PK/PD flags and evidence.
    PK = Pharmacokinetics (how body processes drug — metabolism)
    PD = Pharmacodynamics (how drug affects body — efficacy/toxicity)
    """
    if pk and pk.strip():
        return "Metabolism/PK"
    if pd and pd.strip():
        return "Efficacy/PD"
    if "ClinicalAnnotation" in evidence:
        return "Clinical"
    return "Other"


def build_effect_summary(gene, drug, association, pk, pd):
    """
    Build a human-readable effect summary string.
    This is what shows up on your dashboard.
    """
    category = ""
    if pk and pk.strip():
        category = "metabolism"
    elif pd and pd.strip():
        category = "efficacy/response"

    if gene in CYP_GENES:
        return f"{gene} variant affects {category or 'response'} of {drug}"
    return f"{gene} is {association} with {drug} {category}".strip()


def load_pharmgkb():
    """Load PharmGKB gene-drug relationships into the database."""
    print(f"Opening {PHARMGKB_FILE}...")

    if not PHARMGKB_FILE.exists():
        print(f"ERROR: File not found at {PHARMGKB_FILE}")
        print("Make sure you're running this from ~/genomics-dashboard")
        return

    conn   = connect()
    cursor = conn.cursor()

    # Clear existing data before reloading
    cursor.execute("TRUNCATE TABLE pharmgkb_annotations RESTART IDENTITY;")
    print("Cleared existing PharmGKB data.")

    inserted = 0
    skipped  = 0

    with open(PHARMGKB_FILE, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            # ── Filter 1: Must be Gene → Chemical relationship ──
            if row.get("Entity1_type") != "Gene":
                skipped += 1
                continue
            if row.get("Entity2_type") != "Chemical":
                skipped += 1
                continue

            # ── Filter 2: Must be a confirmed association ──
            association = row.get("Association", "").strip().lower()
            if association not in KEEP_ASSOCIATION:
                skipped += 1
                continue

            # ── Filter 3: Must have clinical or variant evidence ──
            evidence = row.get("Evidence", "")
            evidence_types = set(e.strip() for e in evidence.split(","))
            if not evidence_types.intersection(KEEP_EVIDENCE):
                skipped += 1
                continue

            # ── Extract fields ──
            gene       = clean_value(row.get("Entity1_name"))
            drug       = clean_value(row.get("Entity2_name"))
            pk         = clean_value(row.get("PK"))
            pd         = clean_value(row.get("PD"))
            pmids      = clean_value(row.get("PMIDs"))

            if not gene or not drug:
                skipped += 1
                continue

            # ── Build derived fields ──
            phenotype_category = determine_phenotype_category(pk, pd, evidence)
            effect_summary     = build_effect_summary(gene, drug, association, pk, pd)

            # ── Determine evidence level ──
            # ClinicalAnnotation = stronger evidence (like PharmGKB 1A/1B)
            # VariantAnnotation  = weaker evidence (like PharmGKB 3/4)
            if "ClinicalAnnotation" in evidence_types:
                evidence_level = "1A" if gene in CYP_GENES else "1B"
            else:
                evidence_level = "3"

            # ── Insert row ──
            cursor.execute("""
                INSERT INTO pharmgkb_annotations (
                    gene_name,
                    drug_name,
                    effect_summary,
                    evidence_level,
                    phenotype_category,
                    pmid
                ) VALUES (%s, %s, %s, %s, %s, %s)
            """, (
                gene,
                drug,
                effect_summary,
                evidence_level,
                phenotype_category,
                pmids,
            ))
            inserted += 1

    conn.commit()
    cursor.close()
    conn.close()

    print(f"\n✅ Done.")
    print(f"   Inserted : {inserted:,} gene-drug relationships")
    print(f"   Skipped  : {skipped:,} rows (not Gene-Chemical, ambiguous, or no evidence)")


def verify():
    """Quick check — print counts and sample rows."""
    conn   = connect()
    cursor = conn.cursor()

    # Total count
    cursor.execute("SELECT COUNT(*) FROM pharmgkb_annotations;")
    count = cursor.fetchone()[0]
    print(f"\nTotal rows in pharmgkb_annotations: {count:,}")

    # CYP genes specifically (most important for your project)
    print("\nCYP gene drug relationships (your pharmacogenomics core):")
    cursor.execute("""
        SELECT gene_name, drug_name, phenotype_category, evidence_level
        FROM pharmgkb_annotations
        WHERE gene_name LIKE 'CYP%'
        ORDER BY gene_name, evidence_level
        LIMIT 15;
    """)
    for row in cursor.fetchall():
        print(f"  {row[0]:<12} {row[1]:<35} {row[2]:<20} level {row[3]}")

    # Breakdown by category
    print("\nBreakdown by phenotype category:")
    cursor.execute("""
        SELECT phenotype_category, COUNT(*) as count
        FROM pharmgkb_annotations
        GROUP BY phenotype_category
        ORDER BY count DESC;
    """)
    for row in cursor.fetchall():
        print(f"  {row[0]:<25} : {row[1]:,}")

    # Top 10 most studied genes
    print("\nTop 10 most studied genes:")
    cursor.execute("""
        SELECT gene_name, COUNT(*) as drug_count
        FROM pharmgkb_annotations
        GROUP BY gene_name
        ORDER BY drug_count DESC
        LIMIT 10;
    """)
    for row in cursor.fetchall():
        print(f"  {row[0]:<15} : {row[1]} drugs")

    cursor.close()
    conn.close()


if __name__ == "__main__":
    print("=== LOADING PHARMGKB DATA ===\n")
    load_pharmgkb()
    verify()