"""
modules/ingestion/vcf_parser.py
================================
Updated VCF parser with:
  - Multi-sample support (one user per sample)
  - Auto-detect genome assembly (GRCh38 vs older)
  - Progress updates every 50,000 variants
  - Duplicate upload handling (replace old upload)
  - Skip homozygous_ref variants (no mutation = no risk)

Usage:
    python modules/ingestion/vcf_parser.py
"""

import sys
import gzip
from pathlib import Path
from datetime import datetime

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
    "BRCA1", "BRCA2", "TP53", "PTEN", "APC",
    "MLH1", "MSH2", "MSH6", "CFTR", "HEXA",
    "HBB", "FBN1", "LDLR", "RB1", "VHL",
    "NF1", "NF2",
}

# ── Known GRCh38 assembly strings found in VCF headers ───────
GRCH38_IDENTIFIERS = {
    "grch38","GRCh38", "hg38", "gca_000001405.15",
    "genome reference consortium human build 38"
}

BATCH_SIZE     = 500
PROGRESS_EVERY = 50000


# ─────────────────────────────────────────────────────────────
# SECTION 1: ASSEMBLY DETECTION
# ─────────────────────────────────────────────────────────────

def detect_assembly(header_lines):
    """
    Scan VCF header lines to detect the genome assembly version.

    Returns:
        'GRCh38'  → current standard, ClinVar compatible
        'GRCh37'  → previous standard (hg19)
        'older'   → build 36 or earlier
        'unknown' → no assembly info found in header
    """
    for line in header_lines:
        line_lower = line.lower() 
        if any(ident in line_lower for ident in GRCH38_IDENTIFIERS):
            return "GRCh38"

        if any(x in line_lower for x in ["grch37", "hg19", "b37", "build 37"]):
            return "GRCh37"

        if any(x in line_lower for x in ["b36", "hg18", "build 36", "ncbi36", "b35", "hg17"]):
            return "older"

    return "unknown"


def warn_assembly(assembly, filename):
    """Print a clear warning if the assembly is not GRCh38."""
    if assembly == "GRCh38":
        print(f"  ✅ Assembly: GRCh38 — fully compatible with ClinVar")
    elif assembly == "GRCh37":
        print(f"  ⚠️  Assembly: GRCh37 (hg19) — coordinates differ from ClinVar (GRCh38)")
        print(f"     Annotation matching will be unreliable.")
        print(f"     Consider lifting over to GRCh38 using UCSC liftOver.")
    elif assembly == "older":
        print(f"  ❌ Assembly: pre-GRCh37 — incompatible with ClinVar")
        print(f"     Variants will be stored but annotation matching will fail.")
        print(f"     Strongly recommend lifting over to GRCh38.")
    else:
        print(f"  ⚠️  Assembly: unknown — could not detect from header")
        print(f"     Annotation matching may be unreliable.")


# ─────────────────────────────────────────────────────────────
# SECTION 2: HEADER PARSING
# ─────────────────────────────────────────────────────────────

def open_vcf(filepath):
    """Open a VCF file whether plain text or gzip compressed."""
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8")
    return open(filepath, "r", encoding="utf-8")


def parse_vcf_header(filepath):
    """
    Read the VCF header and extract:
      - All metadata lines (##)
      - Sample names from the #CHROM line
      - Detected assembly version

    Returns:
        header_lines  : list of ## metadata lines
        sample_names  : list of sample column names
        assembly      : detected assembly string
    """
    header_lines = []
    sample_names = []

    with open_vcf(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("##"):
                header_lines.append(line)
            elif line.startswith("#CHROM"):
                fields = line.split("\t")
                sample_names = fields[9:] if len(fields) > 9 else ["SAMPLE"]
                break

    assembly = detect_assembly(header_lines)
    return header_lines, sample_names, assembly


# ─────────────────────────────────────────────────────────────
# SECTION 3: USER AND UPLOAD MANAGEMENT
# ─────────────────────────────────────────────────────────────

def get_or_create_user(username, email):
    """Get existing user by email or create a new one."""
    existing = execute_query(
        "SELECT user_id FROM users WHERE email = %s",
        params=(email,)
    )
    if existing:
        return existing[0]['user_id']
    return execute_insert(
        "INSERT INTO users (username, email) VALUES (%s, %s) RETURNING user_id",
        params=(username, email)
    )


def create_upload_record(user_id, filename, assembly):
    """Create a vcf_uploads record with assembly info in notes."""
    return execute_insert("""
        INSERT INTO vcf_uploads (user_id, filename, status, notes)
        VALUES (%s, %s, 'processing', %s)
        RETURNING upload_id
    """, params=(user_id, filename, f"assembly={assembly}"))


def update_upload_status(upload_id, status, total_variants=0):
    """Update the upload record when processing is done."""
    execute_query("""
        UPDATE vcf_uploads
        SET status = %s, total_variants = %s
        WHERE upload_id = %s
    """, params=(status, total_variants, upload_id), fetch=False)


def delete_existing_upload(user_id, filename):
    """
    If this user already has an upload with this filename,
    delete it and all its variants before re-ingesting.
    Implements the 'replace' duplicate strategy.
    """
    existing = execute_query("""
        SELECT upload_id FROM vcf_uploads
        WHERE user_id = %s AND filename = %s
    """, params=(user_id, filename))

    if not existing:
        return 0

    deleted = 0
    for row in existing:
        old_id = row['upload_id']

        # Delete variant_annotations first (foreign key)
        execute_query("""
            DELETE FROM variant_annotations
            WHERE variant_id IN (
                SELECT variant_id FROM variants WHERE upload_id = %s
            )
        """, params=(old_id,), fetch=False)

        # Delete variants
        execute_query(
            "DELETE FROM variants WHERE upload_id = %s",
            params=(old_id,), fetch=False
        )

        # Delete upload record
        execute_query(
            "DELETE FROM vcf_uploads WHERE upload_id = %s",
            params=(old_id,), fetch=False
        )

        print(f"  🗑️  Replaced previous upload (id={old_id})")
        deleted += 1

    return deleted


def generate_sample_email(filename, sample_name):
    """
    Generate a unique email for a sample from a multi-sample VCF.
    Format: samplename.filename@genomics.local
    """
    clean_filename = Path(filename).stem.replace(".", "_").replace("-", "_")[:30]
    clean_sample   = sample_name.replace(".", "_").replace("-", "_")
    return f"{clean_sample}.{clean_filename}@genomics.local"


# ─────────────────────────────────────────────────────────────
# SECTION 4: VCF LINE PARSING
# ─────────────────────────────────────────────────────────────

def decode_genotype(gt_raw):
    """Convert raw genotype string to zygosity label."""
    gt_map = {
        "0/0": "homozygous_ref",
        "0/1": "heterozygous",
        "1/0": "heterozygous",
        "1/1": "homozygous_alt",
    }
    # Normalize phased (|) to unphased (/)
    normalized = gt_raw.replace("|", "/") if gt_raw else ""
    return gt_map.get(normalized, "unknown")


def determine_flag(gene_name):
    """Determine the pipeline flag for a variant based on its gene."""
    if gene_name is None:
        return None
    gene       = gene_name.upper().strip()
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
    """Basic validation before inserting a variant."""
    if not chrom:
        return False, "missing chromosome"
    if pos is None or pos < 1:
        return False, f"invalid position: {pos}"
    if not ref or ref == ".":
        return False, "missing REF allele"
    if not alt or alt == ".":
        return False, "missing ALT allele"
    valid_chars = set("ACGTNacgtn*<>[].,")
    if not all(c in valid_chars for c in ref):
        return False, f"invalid REF: {ref}"
    return True, "ok"


def parse_info_field(info_str):
    """Parse the INFO field into a dictionary."""
    info = {}
    if not info_str or info_str == ".":
        return info
    for item in info_str.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key.strip()] = value.strip()
        else:
            info[item.strip()] = True
    return info


def parse_variant_line(line, sample_index=0):
    """
    Parse one VCF data line for a specific sample.

    Args:
        line         : raw VCF line string
        sample_index : which sample column to read (0 = first sample)

    Returns:
        dict of parsed fields, or None if malformed
    """
    fields = line.strip().split("\t")
    if len(fields) < 8:
        return None

    try:
        chrom = fields[0].strip()
        pos   = int(fields[1])
        rsid  = fields[2] if fields[2] != "." else None
        ref   = fields[3].strip().upper()
        alt   = fields[4].strip().upper()
        qual  = fields[5] if fields[5] != "." else None
        info  = parse_info_field(fields[7])

        # Normalize chromosome format
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        # Quality score
        qual_float = None
        if qual:
            try:
                qual_float = float(qual)
            except ValueError:
                pass

        # Depth from INFO
        depth = None
        if "DP" in info:
            try:
                depth = int(info["DP"])
            except (ValueError, TypeError):
                pass

        # Allele frequency from INFO
        af = None
        if "AF" in info:
            try:
                af = float(info["AF"])
            except (ValueError, TypeError):
                pass

        # Genotype for requested sample
        gt_raw     = None
        zygosity   = None
        sample_col = 9 + sample_index

        if len(fields) > sample_col and len(fields) > 8:
            fmt_keys    = fields[8].split(":")
            sample_vals = fields[sample_col].split(":")
            genotype    = dict(zip(fmt_keys, sample_vals))
            gt_raw      = genotype.get("GT")

            if gt_raw:
                zygosity = decode_genotype(gt_raw)

            # Depth from FORMAT if not in INFO
            if depth is None and "DP" in genotype:
                try:
                    depth = int(genotype["DP"])
                except (ValueError, TypeError):
                    pass

        # Gene name from INFO
        gene_name = info.get("Gene") or info.get("GENE") or info.get("ANN")
        if gene_name and "," in str(gene_name):
            gene_name = str(gene_name).split(",")[0]

        return {
            "chrom"    : chrom,
            "pos"      : pos,
            "rsid"     : rsid,
            "ref"      : ref,
            "alt"      : alt,
            "qual"     : qual_float,
            "depth"    : depth,
            "af"       : af,
            "gene_name": gene_name,
            "zygosity" : zygosity,
        }

    except (ValueError, IndexError):
        return None


# ─────────────────────────────────────────────────────────────
# SECTION 5: DATABASE INSERT
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
# SECTION 6: SINGLE SAMPLE INGESTION
# ─────────────────────────────────────────────────────────────

def ingest_single_sample(filepath, username, email,
                          sample_index=0, sample_name="SAMPLE",
                          assembly="unknown"):
    """Ingest one sample from a VCF file."""
    filepath  = Path(filepath)
    user_id   = get_or_create_user(username, email)
    filename  = f"{filepath.name}__sample_{sample_name}"

    # Handle duplicate — delete previous upload if exists
    delete_existing_upload(user_id, filename)

    upload_id = create_upload_record(user_id, filename, assembly)
    conn      = get_connection()
    cursor    = conn.cursor()

    inserted = 0
    skipped  = 0
    flagged  = {"pharmacogenomics": 0, "clinical": 0, "both": 0}
    batch    = []

    try:
        with open_vcf(filepath) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                variant = parse_variant_line(line, sample_index=sample_index)
                if variant is None:
                    skipped += 1
                    continue

                # Skip homozygous reference — no mutation, no risk
                if variant["zygosity"] == "homozygous_ref":
                    skipped += 1
                    continue

                # Validate
                is_valid, _ = validate_variant(
                    variant["chrom"], variant["pos"],
                    variant["ref"],   variant["alt"]
                )
                if not is_valid:
                    skipped += 1
                    continue

                # Flag
                flag = determine_flag(variant["gene_name"])
                if flag and flag in flagged:
                    flagged[flag] += 1

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

                if len(batch) >= BATCH_SIZE:
                    insert_variants_batch(cursor, upload_id, batch)
                    conn.commit()
                    inserted += len(batch)
                    batch = []
                    if inserted % PROGRESS_EVERY == 0:
                        print(f"    ↳ {inserted:,} variants inserted...")

        if batch:
            insert_variants_batch(cursor, upload_id, batch)
            conn.commit()
            inserted += len(batch)

        update_upload_status(upload_id, "complete", inserted)
        return upload_id, inserted, skipped, flagged

    except Exception as e:
        conn.rollback()
        update_upload_status(upload_id, "failed")
        print(f"  ❌ Failed on sample {sample_name}: {e}")
        raise

    finally:
        cursor.close()
        conn.close()


# ─────────────────────────────────────────────────────────────
# SECTION 7: MAIN INGESTION FUNCTION
# ─────────────────────────────────────────────────────────────

def ingest_vcf(filepath, username=None, email=None):
    """
    Main entry point. Handles both single and multi-sample VCFs.
    """
    filepath   = Path(filepath)
    start_time = datetime.now()

    if not filepath.exists():
        print(f"❌ File not found: {filepath}")
        return None

    print(f"\n{'='*60}")
    print(f"  INGESTING : {filepath.name}")
    print(f"  Started   : {start_time.strftime('%H:%M:%S')}")
    print(f"{'='*60}")

    # Parse header
    print("\n[1/3] Reading VCF header...")
    _, sample_names, assembly = parse_vcf_header(filepath)
    warn_assembly(assembly, filepath.name)
    print(f"  Samples found : {len(sample_names)}")

    if len(sample_names) > 1:
        preview = ', '.join(sample_names[:5])
        suffix  = f"... +{len(sample_names)-5} more" if len(sample_names) > 5 else ""
        print(f"  Sample names  : {preview}{suffix}")

    # Ingest each sample
    print(f"\n[2/3] Ingesting samples...")
    results = []

    for i, sample_name in enumerate(sample_names):
        print(f"\n  Sample {i+1}/{len(sample_names)}: {sample_name}")

        if len(sample_names) == 1 and username and email:
            s_username = username
            s_email    = email
        else:
            s_username = f"{Path(filepath.stem).stem}_{sample_name}"
            s_email    = generate_sample_email(filepath.name, sample_name)

        upload_id, inserted, skipped, flagged = ingest_single_sample(
            filepath     = filepath,
            username     = s_username,
            email        = s_email,
            sample_index = i,
            sample_name  = sample_name,
            assembly     = assembly,
        )

        results.append({
            "sample"   : sample_name,
            "upload_id": upload_id,
            "inserted" : inserted,
            "skipped"  : skipped,
            "flagged"  : flagged,
        })

        print(f"    ✅ Inserted: {inserted:,} | Skipped: {skipped:,} | "
              f"PGx: {flagged['pharmacogenomics']} | "
              f"Clinical: {flagged['clinical']}")

    # Summary
    elapsed        = (datetime.now() - start_time).seconds
    total_inserted = sum(r['inserted'] for r in results)
    total_skipped  = sum(r['skipped']  for r in results)

    print(f"\n[3/3] Complete!")
    print(f"{'='*60}")
    print(f"  File          : {filepath.name}")
    print(f"  Assembly      : {assembly}")
    print(f"  Samples       : {len(sample_names)}")
    print(f"  Total inserted: {total_inserted:,}")
    print(f"  Total skipped : {total_skipped:,}")
    print(f"  Time elapsed  : {elapsed}s")
    print(f"{'='*60}")

    return results


# ─────────────────────────────────────────────────────────────
# SECTION 8: TEST RUN
# ─────────────────────────────────────────────────────────────

if __name__ == "__main__":

    sample_vcf = Path("data/raw/multi_sample_test.vcf")

    if not sample_vcf.exists():
        print("Creating sample VCF for testing...")
        sample_content = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPATIENT_001
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
        print(f"  Created {sample_vcf}\n")

    results = ingest_vcf(
        filepath = sample_vcf,
        username = "younes",
        email    = "younesadi18@gmail.com"
    )

if results:
    print("\n── Verifying inserted variants ──")
    for result in results:
        print(f"\n  Sample: {result['sample']} (upload_id={result['upload_id']})")
        variants = execute_query("""
            SELECT chromosome, position, gene_name, zygosity, flag
            FROM variants
            WHERE upload_id = %s
            ORDER BY chromosome, position
        """, params=(result['upload_id'],))

        print(f"  {'CHROM':<8} {'POS':<12} {'GENE':<10} {'ZYGOSITY':<20} {'FLAG'}")
        print(f"  {'-'*65}")
        for v in variants:
            print(
                f"  {v['chromosome']:<8} "
                f"{v['position']:<12} "
                f"{str(v['gene_name']):<10} "
                f"{str(v['zygosity']):<20} "
                f"{str(v['flag'])}"
            )
        print(f"  Total: {len(variants)} variants")