import gzip
import csv
import psycopg2
from pathlib import Path

DB_CONFIG = {
    "host"     : "localhost",
    "database" : "genomics_db",
    "user"     : "genomics_user",
    "password" : "genomics123",
    "port"     : "5432"
}

CLINVAR_FILE = Path("data/raw/variant_summary.txt.gz")

KEEP_SIGNIFICANCE = {
    "Pathogenic",
    "Likely pathogenic",
    "Uncertain significance",
    "Pathogenic/Likely pathogenic",
}

KEEP_ASSEMBLY = "GRCh38"
BATCH_SIZE = 1000


def connect():
    """Create and return a database connection."""
    return psycopg2.connect(**DB_CONFIG)


def clean_value(value):
    if value in ("-", "", "na", "NA", "N/A"):
        return None
    return value


def load_clinvar(limit=None):
    print(f"Opening {CLINVAR_FILE}...")

    if not CLINVAR_FILE.exists():
        print(f"ERROR: File not found at {CLINVAR_FILE}")
        print("Make sure you're running this from ~/genomics-dashboard")
        return

    conn   = connect()
    cursor = conn.cursor()

    cursor.execute("TRUNCATE TABLE clinvar_annotations RESTART IDENTITY;")
    print("Cleared existing ClinVar data.")

    inserted  = 0
    skipped   = 0
    batch     = []

    with gzip.open(CLINVAR_FILE, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for i, row in enumerate(reader):

            if i % 500000 == 0 and i > 0:
                print(f"  Scanned {i:,} rows — inserted {inserted:,} so far...")

            if limit and inserted >= limit:
                break

            if row.get("Assembly") != KEEP_ASSEMBLY:
                skipped += 1
                continue

            significance = row.get("ClinicalSignificance", "")
            if not any(s in significance for s in KEEP_SIGNIFICANCE):
                skipped += 1
                continue

            gene = clean_value(row.get("GeneSymbol"))
            if gene is None:
                skipped += 1
                continue

            chromosome = clean_value(row.get("Chromosome"))
            if chromosome:
                chromosome = f"chr{chromosome}"  

            position   = clean_value(row.get("PositionVCF"))
            if position:
                try:
                    position = int(position)
                except ValueError:
                    skipped += 1
                    continue
                if position < 1:
                    skipped +=1
                    continue

            ref_allele = clean_value(row.get("ReferenceAlleleVCF"))
            alt_allele = clean_value(row.get("AlternateAlleleVCF"))
            condition  = clean_value(row.get("PhenotypeList"))
            if condition:
                condition = condition.split("|")[0].strip()
                if not condition:
                    condition = None
            review     = clean_value(row.get("ReviewStatus"))
            accession  = clean_value(row.get("RCVaccession"))

            # ── Add to batch ──
            batch.append((
                chromosome,
                position,
                ref_allele,
                alt_allele,
                gene,
                significance,
                condition,
                review,
                accession,
            ))

            if len(batch) >= BATCH_SIZE:
                insert_batch(cursor, batch)
                inserted += len(batch)
                batch = []

        if batch:
            insert_batch(cursor, batch)
            inserted += len(batch)

    conn.commit()
    cursor.close()
    conn.close()

    print(f"\n✅ Done.")
    print(f"   Inserted : {inserted:,} variants")
    print(f"   Skipped  : {skipped:,} rows (wrong assembly, benign, or no gene)")


def insert_batch(cursor, batch):
    """Insert a batch of rows into clinvar_annotations."""
    cursor.executemany("""
        INSERT INTO clinvar_annotations (
            chromosome,
            position,
            ref_allele,
            alt_allele,
            gene_name,
            clinical_significance,
            condition_name,
            review_status,
            clinvar_accession
        ) VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
    """, batch)


def verify():
    conn   = connect()
    cursor = conn.cursor()

    cursor.execute("SELECT COUNT(*) FROM clinvar_annotations;")
    count = cursor.fetchone()[0]
    print(f"\nTotal rows in clinvar_annotations: {count:,}")

    print("\nSample pathogenic variants:")
    cursor.execute("""
        SELECT gene_name, chromosome, position, clinical_significance, condition_name
        FROM clinvar_annotations
        WHERE clinical_significance = 'Pathogenic'
        LIMIT 10;
    """)
    rows = cursor.fetchall()
    for row in rows:
        print(f"  {row[0]:<15} {row[1]:<8} {row[2]:<12} {row[3]:<20} {str(row[4])[:50]}")

    print("\nVariants by significance:")
    cursor.execute("""
        SELECT clinical_significance, COUNT(*) as count
        FROM clinvar_annotations
        GROUP BY clinical_significance
        ORDER BY count DESC;
    """)
    for row in cursor.fetchall():
        print(f"  {row[0]:<40} : {row[1]:,}")

    cursor.close()
    conn.close()


if __name__ == "__main__":
    print("=== LOADING CLINVAR DATA (test: 5,000 rows) ===\n")
    load_clinvar(limit=None)
    verify()

    print("\n" + "=" * 50)
    print("Test load successful!")
    print("To load the full dataset, edit the last line to:")
    print("  load_clinvar(limit=None)")
    print("=" * 50)