import sys
import gzip
from pathlib import Path
from datetime import datetime
sys.path.append(str(Path(__file__).parent.parent.parent))
import re

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

BATCH_SIZE     = 50000
PROGRESS_EVERY = 50000
BATCH_SIZE = 5000



# ASSEMBLY DETECTION: FUNCTION TO DETECT IN WHICH ASSEMBLY THE VCF FILE FORMATED AND NOT ACCEPT OLD FORMATS

def detect_assembly(header_lines):
    for line in header_lines:
        if 'GRCh38' in line or 'hg38' in line:
            return "GRCh38"
        
    return None

# A simple open vcf file

def open_vcf(filepath):
    filepath = Path(filepath)
    if not (str(filepath).endswith('.vcf') or str(filepath).endswith('.vcf.gz')):
        print("Rejected, only .vcf and .vcf.gz files accepted")
        return None
    elif str(filepath).endswith('.vcf.gz'):
        return gzip.open(filepath,"rt", encoding="utf-8")
    else:
        return open(filepath, "r", encoding="utf-8")


def load_gene_index() -> dict:
    """Load all gene + exon coordinates into memory for fast lookup."""
    

    print("  Loading gene index into memory...")
    index = {}

    # Load exons first
    exons = execute_query(
        """
            SELECT chromosome, start_pos, end_pos, gene_name 
            FROM exon_coordinates
            ORDER BY chromosome, start_pos;
        """
    )

    for row in exons:
        chrom = row['chromosome']
        if chrom not in index:
            index[chrom] = []
        index[chrom].append((row['start_pos'], row['end_pos'], row['gene_name'], 'exonic'))
        # a tuple because fixed, protected data,less memory,Faster than lists
    
    genes = execute_query(
        """
            SELECT chromosome, start_pos, end_pos, gene_name
            from gene_coordinates
            ORDER BY chromosome, start_pos;
        """
    )

    for row in genes:
        chrom = row['chromosome']
        if chrom not in index:
            index[chrom] = []
        index[chrom].append((row['start_pos'], row['end_pos'], row['gene_name'], 'unknown'))

    total = sum(len(v) for v in index.values())
    print(f" Gene index ready {total:,} entries")
    return index

    # Load genes as fallback

def lookup_gene_memory(index: dict, chrom: str, pos: int):
    for start, end, gene_name, region in index.get(chrom, []):
        if start <= pos <= end:
            return gene_name, region
    return None, None

def parse_vcf_header(filepath):
    header_lines = []
    columns = []
    samples = []
    with open_vcf(filepath) as f:
        for line in f:
            if line.startswith("##"):
                header_lines.append(line)
    
            elif line.startswith("#C"):
                cols = line.strip().split("\t")
                cols[0] = cols[0].lstrip('#')
                columns = cols[:9]
                samples = cols[9:]

    return header_lines, columns, samples


#user management

def get_or_create_user(username, email) -> int:
    existing = execute_query(
        """
            SELECT user_id from users
                WHERE email = %s;
        """, (email,)
    )

    if existing:
        return existing[0]['user_id']
    else:
        user_id = execute_insert(
                """
                    INSERT INTO users (username, email)
                    VALUES (%s, %s)
                    RETURNING user_id; 
                """, (username, email,)

        )
        return user_id


def create_upload_record(user_id, filename, notes=None):
    return execute_insert(
        """
            INSERT INTO vcf_uploads (user_id, filename, notes)
            VALUES (%s, %s, %s)
            RETURNING upload_id;
        """, (user_id, filename, notes)
    )


def update_upload_status(upload_id, status, total_variants=0):
    execute_insert(
        """
            UPDATE vcf_uploads
            SET status = %s, total_variants = %s
            WHERE upload_id = %s;
        """,(status,total_variants,upload_id)
    )


def delete_existing_upload(user_id, filename):
    execute_query(
            """
                DELETE FROM vcf_uploads
                WHERE user_id = %s
                AND filename = %s
                RETURNING user_id, filename;
            """,(user_id,filename)
        )
        
def decode_genotype(gt_raw):
    gt_map = {
        "0/0": "homozygous_ref",
        "0/1": "heterozygous",
        "1/0": "heterozygous",
        "1/1": "homozygous_alt",
    }
    normalized = gt_raw.replace("|", "/") if gt_raw else ""
    return gt_map.get(normalized,"unknown")


def determine_flag(gene_name) -> str:
    if not gene_name:
        return None
    
    
    if gene_name in DISEASE_GENES:
        if gene_name in CYP_GENES:
            return "both"
        else:
            return "clinical"
    if gene_name in CYP_GENES:
        return "pharmacogenomics"
    
    return None

def validate_variant(chrom, pos, ref, alt):
    if not chrom.startswith("chr"):
        return False
    if not isinstance(pos, int) or pos <= 0:
        return False
    if not re.match(r'^[ACTG]+$',ref) or not re.match(r'^[ACTG]+$',alt):
        return False
    return True
    


def parse_info_field(info_str) -> dict:
    line = info_str.strip().split(';')

    info = {}
    for item in line:
        if '=' in item:
            key, value = item.split("=", 1)
            info[key] = value
    return info
  

def parse_variant_line(line)->dict:
    fields = line.strip().split("\t")
    
    if len(fields) < 8:
        return None
    
    try:
        chrom = fields[0].strip()
        pos = int(fields[1].strip())
        rs_id = fields[2].strip() if fields[2] != "." else None
        ref = fields[3].strip().upper()
        alt = fields[4].strip().upper()
        qual = float(fields[5].strip()) if fields[5] != "." else None
        info = parse_info_field(fields[7].strip())

        filter_status = fields[6].strip()
        if filter_status != "PASS" and filter_status != ".":
            return None
        
        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        if not validate_variant(chrom, pos, ref, alt):
            return None
        
        depth = int(info.get('DP', 0)) if info.get('DP') else None
        allele_freq = float(info.get('AF', 0)) if info.get('AF') else None
        
        return {
            'chrom'      : chrom,
            'pos'        : pos,
            'rs_id'      : rs_id,
            'ref'        : ref,
            'alt'        : alt,
            'qual'        : qual,
            'depth'      : depth,
            'allele_freq': allele_freq,
            'fields'     : fields,
        }
    except Exception as e:
        return None
        
def extract_sample_data(fields, sample_index, gene_index):
    sample = fields[9 + sample_index].strip().split(":")

    format_str = fields[8].strip().split(":")
    gt_index = format_str.index('GT') if 'GT'in format_str else 0

    genotype = sample[gt_index].strip()
    zygosity = decode_genotype(genotype)

    if zygosity == "homozygous_ref":
        return None
    
    gene_name, region = lookup_gene_memory(gene_index,fields[0].strip(),int(fields[1].strip()))
    
    flag = determine_flag(gene_name)

    return {
            'zygosity'  : zygosity,
            'gene_name' : gene_name,
            'region'    : region,
            'flag'      : flag,
    }

def insert_variants_batch(cursor, batch):
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


def ingest_vcf(filepath, username=None, email=None):
    if not filepath.exists():
        print(f"FILE NOT FOUND: {filepath}")
        return None
    
    conn = get_connection()
    cursor = conn.cursor()

    print("\n" + "-" * 60)
    print(f"INGESTING: {filepath.name}")
    print("-" * 60)

    print("\n[1/3] Reading VCF header...")
    header_lines, _, samples = parse_vcf_header(filepath)
    
    if not detect_assembly(header_lines):
        print("Rejected : only GRCh38 files accepted")
        return None
    
    if not samples:
        print("No samples found in VCF file")
        return None
    
    print(f"Samples found: {len(samples)}")
    
    gene_index = load_gene_index()

    user_id = get_or_create_user(username, email)

    print("ingesting samples...")

    upload_ids = {}
    batches = {}
    inserted = {}

    for sample_name in samples:
        upload_ids[sample_name] = create_upload_record(user_id, filepath.name, f"assembly=GRCh38")
        batches[sample_name] = []
        inserted[sample_name] = 0

    # one pass through the file
    with open_vcf(filepath) as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            base = parse_variant_line(line)
            if base is None:
                continue
            
            for i, sample_name in enumerate(samples):
                sample_data = extract_sample_data(base['fields'], i, gene_index)
                if sample_data is None:
                    continue
                
                batches[sample_name].append((
                    upload_ids[sample_name],
                    base['chrom'], base['pos'],
                    base['ref'],   base['alt'],
                    base['rs_id'],
                    sample_data['gene_name'],
                    sample_data['zygosity'],
                    base['qual'],
                    base['depth'],
                    base['allele_freq'],
                    sample_data['flag'],
                ))

                if len(batches[sample_name]) >= BATCH_SIZE:
                    insert_variants_batch(cursor, batches[sample_name])
                    conn.commit()
                    inserted[sample_name] += len(batches[sample_name])
                    batches[sample_name] = []
                    
        for sample_name in samples:
            if batches[sample_name]:
                insert_variants_batch(cursor, batches[sample_name])
                conn.commit()
                inserted[sample_name] += len(batches[sample_name])
                print(f"  Committed {len(batches[sample_name])} variants")
            
            
            update_upload_status(upload_ids[sample_name], 'complete', inserted[sample_name])
            print(f"  {sample_name}: {inserted[sample_name]:,} variants inserted")
    cursor.close()
    conn.close()

    return {
        'samples'    : samples,
        'upload_ids' : upload_ids,
        'inserted'   : inserted,
    }

if __name__ == "__main__":
    from pathlib import Path
    result = ingest_vcf(
        filepath = Path("data/raw/sample.vcf"),
        username = "test_user",
        email    = "test@genomics.com"
    )
    print(result)