"""
modules/ingestion/vcf_parser.py
================================
Parses a VCF file and inserts variants into the database.
This is the entry point for the entire pipeline.

Workflow:
    1. Create a user in the database (or reuse existing)
    2. Create a vcf_upload record
    3. Parse the VCF file line by line
    4. For each variant:
       a. Decode chromosome, position, REF, ALT, genotype
       b. Validate the data
       c. Check if gene is a CYP gene → flag as 'pharmacogenomics'
       d. Insert into variants table
    5. Update the upload status to 'complete'

Usage:
    python modules/ingestion/vcf_parser.py
"""

import sys
import gzip
from pathlib import Path

# Add project root to path so we can import database.connect
sys.path.append(str(Path(__file__).parent.parent.parent))

from database.connect import get_connection, execute_insert, execute_query


# ── CYP genes — flagged as pharmacogenomics ──────────────────
CYP_GENES = {
    "CYP1A2", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19",
    "CYP2D6", "CYP2E1", "CYP3A4", "CYP3A5", "CYP3A7",
    "DPYD", "TPMT", "UGT1A1", "SLCO1B1", "ABCG2",
    "G6PD", "IFNL3", "HLA-A", "HLA-B"
}

# ── Disease genes — flagged as clinical ──────────────────────
DISEASE_GENES = {
    "BRCA1", "BRCA2",           # breast/ovarian cancer
    "TP53",                      # Li-Fraumeni syndrome
    "PTEN",                      # Cowden syndrome
    "APC",                       # colorectal cancer
    "MLH1", "MSH2", "MSH6",     # Lynch syndrome
    "CFTR",                      # cystic fibrosis
    "HEXA",                      # Tay-Sachs
    "HBB",                       # sickle cell / thalassemia
    "FBN1",                      # Marfan syndrome
    "LDLR",                      # familial hypercholesterolemia
    "RB1",                       # retinoblastoma
    "VHL",                       # Von Hippel-Lindau
    "NF1", "NF2",               # neurofibromatosis
}

# ── Batch size for database inserts ──────────────────────────
BATCH_SIZE = 500


# ─────────────────────────────────────────────────────────────
# SECTION 1: USER AND UPLOAD MANAGEMENT
# ─────────────────────────────────────────────────────────────

def get_or_create_user(username, email):
    """
    Get existing user by email, or create a new one.
    Returns the user_id.
    """
    # Check if user already exists
    existing = execute_query(
        "SELECT user_id FROM users WHERE email = %s",
        params=(email,)
    )
    if existing:
        print(f"  Found existing user: {username} (id={existing[0]['user_id']})")
        return existing[0]['user_id']

    # Create new user
    user_id = execute_insert(
        "INSERT INTO users (username, email) VALUES (%s, %s) RETURNING user_id",
        params=(username, email)
    )
    print(f"  Created new user: {username} (id={user_id})")
    return user_id


def create_upload_record(user_id, filename):
    """
    Create a vcf_uploads record and return the upload_id.
    Status starts as 'processing'.
    """
    upload_id = execute_insert("""
        INSERT INTO vcf_uploads (user_id, filename, status)
        VALUES (%s, %s, 'processing')
        RETURNING upload_id
    """, params=(user_id, filename))
    return upload_id


def update_upload_status(upload_id, status, total_variants=0):
    """Update the upload record when processing is done."""
    execute_query("""
        UPDATE vcf_uploads
        SET status = %s, total_variants = %s
        WHERE upload_id = %s
    """, params=(status, total_variants, upload_id), fetch=False)


# ─────────────────────────────────────────────────────────────
# SECTION 2: VCF PARSING HELPERS
# ─────────────────────────────────────────────────────────────

def open_vcf(filepath):
    """Open a VCF file whether plain text or gzip compressed."""
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8")
    return open(filepath, "r", encoding="utf-8")


def decode_genotype(gt_raw):
    """
    Convert a raw genotype string to a zygosity label.

    GT values:
      0/0 or 0|0 → homozygous reference (no mutation)
      0/1 or 0|1 → heterozygous (one mutated copy)
      1/1 or 1|1 → homozygous alternate (both copies mutated)
    """
    gt_map = {
        "0/0": "homozygous_ref",
        "0|0": "homozygous_ref",
        "0/1": "heterozygous",
        "1/0": "heterozygous",
        "0|1": "heterozygous",
        "1|0": "heterozygous",
        "1/1": "homozygous_alt",
        "1|1": "homozygous_alt",
    }
    return gt_map.get(gt_raw, "unknown")


def determine_flag(gene_name):
    """
    Determine the pipeline flag for a variant based on its gene.

    Returns:
        'pharmacogenomics' → CYP gene, send to PGx module
        'clinical'         → known disease gene, send to pathogenicity module
        'both'             → gene is in both lists
        None               → no special flag
    """
    if gene_name is None:
        return None

    gene = gene_name.upper().strip()
    in_cyp     = gene in CYP_GENES
    in_disease = gene in DISEASE_GENES

    if in_cyp and in_disease:
        return "both"
    if in_cyp:
        return "pharmacogenomics"
    if in_disease:
        return "clinical"
    return None


def validate_variant(chrom, pos, ref, alt):
    """
    Basic validation before inserting a variant.
    Returns (is_valid, reason) tuple.
    """
    if not chrom:
        return False, "missing chromosome"

    if pos is None or pos < 1:
        return False, f"invalid position: {pos}"

    if not ref or ref == ".":
        return False, "missing REF allele"

    if not alt or alt == ".":
        return False, "missing ALT allele"

    # Check for valid DNA characters only
    valid_chars = set("ACGTNacgtn*<>[].,")
    if not all(c in valid_chars for c in ref):
        return False, f"invalid REF characters: {ref}"

    return True, "ok"


def parse_info_field(info_str):
    """
    Parse the INFO field into a dictionary.
    Example: 'DP=30;AF=0.5;Gene=BRCA1' → {'DP': '30', 'AF': '0.5', 'Gene': 'BRCA1'}
    """
    info = {}
    if not info_str or info_str == ".":
        return info
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key.strip()] = value.strip()
        else:
            info[item.strip()] = True   # flag fields with no value
    return info


def parse_variant_line(line):
    """
    Parse one VCF data line into a clean dictionary.
    Returns None if the line is malformed.
    """
    fields = line.strip().split("\t")

    # VCF requires at least 8 columns
    if len(fields) < 8:
        return None

    try:
        chrom  = fields[0].strip()
        pos    = int(fields[1])
        rsid   = fields[2] if fields[2] != "." else None
        ref    = fields[3].strip().upper()
        alt    = fields[4].strip().upper()
        qual   = fields[5] if fields[5] != "." else None
        filter_= fields[6]
        info   = parse_info_field(fields[7])

        # Normalize chromosome format
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        # Parse quality score
        qual_float = None
        if qual:
            try:
                qual_float = float(qual)
            except ValueError:
                qual_float = None

        # Parse depth and allele frequency from INFO
        depth = None
        if "DP" in info:
            try:
                depth = int(info["DP"])
            except (ValueError, TypeError):
                depth = None

        af = None
        if "AF" in info:
            try:
                af = float(info["AF"])
            except (ValueError, TypeError):
                af = None

        # Parse genotype from FORMAT + SAMPLE columns
        gt_raw   = None
        zygosity = None
        if len(fields) >= 10:
            fmt_keys    = fields[8].split(":")
            sample_vals = fields[9].split(":")
            genotype    = dict(zip(fmt_keys, sample_vals))
            gt_raw      = genotype.get("GT", None)
            if gt_raw:
                zygosity = decode_genotype(gt_raw)

            # Also try to get depth from FORMAT if not in INFO
            if depth is None and "DP" in genotype:
                try:
                    depth = int(genotype["DP"])
                except (ValueError, TypeError):
                    depth = None

        # Try to get gene name from INFO field
        gene_name = info.get("Gene") or info.get("GENE") or info.get("ANN", None)
        if gene_name and "," in str(gene_name):
            gene_name = gene_name.split(",")[0]

        return {
            "chrom"    : chrom,
            "pos"      : pos,
            "rsid"     : rsid,
            "ref"      : ref,
            "alt"      : alt,
            "qual"     : qual_float,
            "filter"   : filter_,
            "depth"    : depth,
            "af"       : af,
            "gene_name": gene_name,
            "zygosity" : zygosity,
        }

    except (ValueError, IndexError) as e:
        return None


# ─────────────────────────────────────────────────────────────
# SECTION 3: DATABASE INSERT
# ─────────────────────────────────────────────────────────────

def insert_variants_batch(cursor, upload_id, batch):
    """Insert a batch of parsed variants into the database."""
    cursor.executemany("""
        INSERT INTO variants (
            upload_id, chromosome, position, ref_allele, alt_allele,
            variant_id_rs, gene_name, zygosity, quality_score,
            depth, allele_freq, flag
        ) VALUES (
            %s, %s, %s, %s, %s,
            %s, %s, %s, %s,
            %s, %s, %s
        )
    """, batch)


# ─────────────────────────────────────────────────────────────
# SECTION 4: MAIN INGESTION FUNCTION
# ─────────────────────────────────────────────────────────────

def ingest_vcf(filepath, username, email):
    """
    Full pipeline: parse a VCF file and insert all variants.

    Args:
        filepath : path to the .vcf or .vcf.gz file
        username : name of the user uploading
        email    : email of the user uploading

    Returns:
        upload_id : the ID of the upload record created
    """
    filepath = Path(filepath)
    if not filepath.exists():
        print(f"❌ File not found: {filepath}")
        return None

    print(f"\n{'='*55}")
    print(f"  INGESTING: {filepath.name}")
    print(f"{'='*55}")

    # ── Step 1: Create user and upload record ──
    print("\n[1/4] Setting up user and upload record...")
    user_id   = get_or_create_user(username, email)
    upload_id = create_upload_record(user_id, filepath.name)
    print(f"  Upload record created (id={upload_id})")

    # ── Step 2: Parse the VCF ──
    print("\n[2/4] Parsing VCF file...")
    conn   = get_connection()
    cursor = conn.cursor()

    inserted  = 0
    skipped   = 0
    flagged   = {"pharmacogenomics": 0, "clinical": 0, "both": 0}
    batch     = []

    try:
        with open_vcf(filepath) as f:
            for line in f:
                # Skip header lines
                if line.startswith("#"):
                    continue

                # Parse the variant
                variant = parse_variant_line(line)
                if variant is None:
                    skipped += 1
                    continue

                # Validate
                is_valid, reason = validate_variant(
                    variant["chrom"], variant["pos"],
                    variant["ref"],   variant["alt"]
                )
                if not is_valid:
                    skipped += 1
                    continue

                # Determine flag
                flag = determine_flag(variant["gene_name"])
                if flag and flag in flagged:
                    flagged[flag] += 1

                # Add to batch
                batch.append((
                    upload_id,
                    variant["chrom"],
                    variant["pos"],
                    variant["ref"],
                    variant["alt"],
                    variant["rsid"],
                    variant["gene_name"],
                    variant["zygosity"],
                    variant["qual"],
                    variant["depth"],
                    variant["af"],
                    flag,
                ))

                # Insert when batch is full
                if len(batch) >= BATCH_SIZE:
                    insert_variants_batch(cursor, upload_id, batch)
                    inserted += len(batch)
                    batch = []

        # Insert remaining
        if batch:
            insert_variants_batch(cursor, upload_id, batch)
            inserted += len(batch)

        conn.commit()

        # ── Step 3: Update upload status ──
        print("\n[3/4] Updating upload status...")
        update_upload_status(upload_id, "complete", inserted)

        # ── Step 4: Print summary ──
        print("\n[4/4] Ingestion complete!\n")
        print(f"  File          : {filepath.name}")
        print(f"  User          : {username}")
        print(f"  Upload ID     : {upload_id}")
        print(f"  Inserted      : {inserted:,} variants")
        print(f"  Skipped       : {skipped:,} invalid lines")
        print(f"\n  Flagged variants:")
        print(f"    Pharmacogenomics : {flagged['pharmacogenomics']:,}")
        print(f"    Clinical         : {flagged['clinical']:,}")
        print(f"    Both             : {flagged['both']:,}")
        print(f"    Unflagged        : {inserted - sum(flagged.values()):,}")

        return upload_id

    except Exception as e:
        conn.rollback()
        update_upload_status(upload_id, "failed")
        print(f"\n❌ Ingestion failed: {e}")
        raise

    finally:
        cursor.close()
        conn.close()


# ─────────────────────────────────────────────────────────────
# SECTION 5: TEST RUN
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Use the sample VCF we created in Week 1
    # Make sure sample.vcf is in your data/raw folder
    sample_vcf = Path("data/raw/CEU.low_coverage.2010_09.genotypes.vcf.gz")

    # Create sample.vcf if it doesn't exist yet
    if not sample_vcf.exists():
        print("Creating sample VCF for testing...")
        sample_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t925952\trs1234567\tG\tA\t50\tPASS\tDP=30;AF=0.5;Gene=CFTR\tGT:DP\t0/1:30
chr2\t179415139\trs3456789\tA\tG\t75\tPASS\tDP=60;AF=0.75;Gene=CYP2D6\tGT:DP\t1/1:60
chr7\t117548628\trs4567890\tT\tC\t88\tPASS\tDP=22;AF=0.5;Gene=CYP2C19\tGT:DP\t0/1:22
chr13\t32339916\trs5678901\tG\tT\t92\tPASS\tDP=55;AF=0.5;Gene=BRCA2\tGT:DP\t0/1:55
chr17\t41244429\trs6789012\tA\tC\t95\tPASS\tDP=70;AF=1.0;Gene=BRCA1\tGT:DP\t1/1:70
chr17\t41246709\trs7890123\tC\tG\t61\tPASS\tDP=33;AF=0.5;Gene=BRCA1\tGT:DP\t0/1:33
chr19\t15879621\trs8901234\tT\tA\t44\tPASS\tDP=18;AF=0.5;Gene=CYP3A4\tGT:DP\t0/1:18
chr22\t19724571\trs9012345\tG\tC\t83\tPASS\tDP=40;AF=0.25;Gene=NF2\tGT:DP\t0/0:40
chr22\t19724772\trs1902345\tT\tG\t77\tPASS\tDP=28;AF=0.5;Gene=TP53\tGT:DP\t0/1:28
"""
        sample_vcf.write_text(sample_content)
        print(f"  Created {sample_vcf}")

    # Run the ingestion
    upload_id = ingest_vcf(
        filepath = sample_vcf,
        username = "CEU_sample",
        email    = "ceu@1000genomes.org"
    )

    if upload_id:
        # Verify what was inserted
        print("\n── Verifying inserted variants ──")
        variants = execute_query("""
            SELECT chromosome, position, gene_name, zygosity, flag
            FROM variants
            WHERE upload_id = %s
            ORDER BY chromosome, position
        """, params=(upload_id,))

        print(f"\n  {'CHROM':<8} {'POS':<12} {'GENE':<10} {'ZYGOSITY':<20} {'FLAG'}")
        print(f"  {'-'*65}")
        for v in variants:
            print(
                f"  {v['chromosome']:<8} "
                f"{v['position']:<12} "
                f"{str(v['gene_name']):<10} "
                f"{str(v['zygosity']):<20} "
                f"{str(v['flag'])}"
            )