import sys
import gzip
import re
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent))

from database.connect import get_connection

GTF_PATH   = Path("data/raw/Homo_sapiens.GRCh38.109.gtf.gz")
BATCH_SIZE = 1000



def parse_attributes(attr_string):

    attrs = {}
    for match in re.finditer(r'(\w+)\s+"([^"]+)"', attr_string):
        attrs[match.group(1)] = match.group(2)
    return attrs


def normalize_chrom(chrom):
    chrom = str(chrom).strip()
    # GTF uses '1', '2', 'X', 'MT' we need 'chr1', 'chr2', 'chrX', 'chrM'
    if chrom == "MT":
        return "chrM"
    if not chrom.startswith("chr"):
        return f"chr{chrom}"
    return chrom


def is_standard_chrom(chrom):
    standard = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
    return chrom in standard


def create_schema(conn):
    cursor = conn.cursor()
    cursor.execute(CREATE_TABLES_SQL)
    conn.commit()
    cursor.close()
    print("Schema created (gene_coordinates + exon_coordinates)")


def check_already_loaded(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT COUNT(*) FROM gene_coordinates")
    count = cursor.fetchone()[0]
    cursor.close()
    return count > 0


def load_gtf(gtf_path=GTF_PATH):
    gtf_path = Path(gtf_path)
    if not gtf_path.exists():
        print(f"❌ GTF file not found: {gtf_path}")
        print("   Download with:")
        print("   wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/"
              "Homo_sapiens.GRCh38.109.gtf.gz -P data/raw/")
        return

    conn   = get_connection()
    create_schema(conn)

    if check_already_loaded(conn):
        print("  ⚠️  gene_coordinates table already has data.")
        answer = input("  Re-load? This will delete existing data. (y/n): ").strip().lower()
        if answer != 'y':
            print("  Skipping load.")
            conn.close()
            return
        cursor = conn.cursor()
        cursor.execute("TRUNCATE TABLE exon_coordinates RESTART IDENTITY CASCADE")
        cursor.execute("TRUNCATE TABLE gene_coordinates RESTART IDENTITY CASCADE")
        conn.commit()
        cursor.close()
        print("  🗑️  Cleared existing data")

    cursor = conn.cursor()

    # ── Pass 1: Load genes ──
    print(f"\n[1/3] Loading genes from {gtf_path.name}...")

    gene_batch   = []
    genes_loaded = 0
    genes_seen   = set() 

    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature != "gene":
                continue

            chrom  = normalize_chrom(fields[0])
            if not is_standard_chrom(chrom):
                continue

            attrs   = parse_attributes(fields[8])
            gene_name = attrs.get("gene_name")
            gene_id   = attrs.get("gene_id")
            biotype   = attrs.get("gene_biotype", "")

            if not gene_name:
                continue

            key = (gene_name, chrom)
            if key in genes_seen:
                continue
            genes_seen.add(key)

            try:
                start  = int(fields[3])
                end    = int(fields[4])
                strand = fields[6] if fields[6] in ('+', '-') else None
            except ValueError:
                continue

            gene_batch.append((
                gene_name, gene_id, chrom,
                start, end, strand, biotype
            ))

            if len(gene_batch) >= BATCH_SIZE:
                cursor.executemany("""
                    INSERT INTO gene_coordinates
                        (gene_name, gene_id, chromosome, start_pos, end_pos, strand, biotype)
                    VALUES (%s, %s, %s, %s, %s, %s, %s)
                """, gene_batch)
                conn.commit()
                genes_loaded += len(gene_batch)
                gene_batch = []

                if genes_loaded % 10000 == 0:
                    print(f"    ↳ {genes_loaded:,} genes loaded...")

    if gene_batch:
        cursor.executemany("""
            INSERT INTO gene_coordinates
                (gene_name, gene_id, chromosome, start_pos, end_pos, strand, biotype)
            VALUES (%s, %s, %s, %s, %s, %s, %s)
        """, gene_batch)
        conn.commit()
        genes_loaded += len(gene_batch)

    print(f"  ✅ {genes_loaded:,} genes loaded")

    print("\n[2/3] Building gene lookup map...")
    cursor.execute("""
        SELECT gene_coord_id, gene_name, chromosome
        FROM gene_coordinates
    """)
    gene_id_map = {
        (row[1], row[2]): row[0]
        for row in cursor.fetchall()
    }
    print(f"  ✅ {len(gene_id_map):,} entries in lookup map")

    print(f"\n[3/3] Loading exons...")

    exon_batch   = []
    exons_loaded = 0
    exons_skipped = 0

    with gzip.open(gtf_path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature != "exon":
                continue

            chrom = normalize_chrom(fields[0])
            if not is_standard_chrom(chrom):
                continue

            attrs     = parse_attributes(fields[8])
            gene_name = attrs.get("gene_name")
            transcript_id = attrs.get("transcript_id")
            exon_number   = attrs.get("exon_number")

            if not gene_name:
                exons_skipped += 1
                continue

            gene_coord_id = gene_id_map.get((gene_name, chrom))
            if not gene_coord_id:
                exons_skipped += 1
                continue

            try:
                start  = int(fields[3])
                end    = int(fields[4])
                strand = fields[6] if fields[6] in ('+', '-') else None
                exon_n = int(exon_number) if exon_number else None
            except (ValueError, TypeError):
                exons_skipped += 1
                continue

            exon_batch.append((
                gene_coord_id, gene_name, chrom,
                start, end, strand, exon_n, transcript_id
            ))

            if len(exon_batch) >= BATCH_SIZE:
                cursor.executemany("""
                    INSERT INTO exon_coordinates
                        (gene_coord_id, gene_name, chromosome,
                         start_pos, end_pos, strand, exon_number, transcript_id)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                """, exon_batch)
                conn.commit()
                exons_loaded += len(exon_batch)
                exon_batch = []

                if exons_loaded % 100000 == 0:
                    print(f"    ↳ {exons_loaded:,} exons loaded...")

    if exon_batch:
        cursor.executemany("""
            INSERT INTO exon_coordinates
                (gene_coord_id, gene_name, chromosome,
                 start_pos, end_pos, strand, exon_number, transcript_id)
            VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
        """, exon_batch)
        conn.commit()
        exons_loaded += len(exon_batch)

    print(f"  ✅ {exons_loaded:,} exons loaded")
    print(f"  ○  {exons_skipped:,} exons skipped (no gene match)")

    cursor.close()
    conn.close()

    print(f"\n{'='*50}")
    print(f"  GTF LOAD COMPLETE")
    print(f"  Genes  : {genes_loaded:,}")
    print(f"  Exons  : {exons_loaded:,}")
    print(f"{'='*50}")



def lookup_gene(chromosome, position):
    from database.connect import execute_query

    exon_hit = execute_query("""
        SELECT gene_name, strand
        FROM exon_coordinates
        WHERE chromosome = %s
        AND   start_pos  <= %s
        AND   end_pos    >= %s
        LIMIT 1
    """, params=(chromosome, position, position))

    if exon_hit:
        return {
            "gene_name" : exon_hit[0]['gene_name'],
            "region"    : "exonic",
            "strand"    : exon_hit[0]['strand'],
        }


    gene_hit = execute_query("""
        SELECT gene_name, strand
        FROM gene_coordinates
        WHERE chromosome = %s
        AND   start_pos  <= %s
        AND   end_pos    >= %s
        LIMIT 1
    """, params=(chromosome, position, position))

    if gene_hit:
        return {
            "gene_name" : gene_hit[0]['gene_name'],
            "region"    : "intronic",
            "strand"    : gene_hit[0]['strand'],
        }

    return None


if __name__ == "__main__":
    load_gtf()


    print("\n── Testing coordinate lookup ──\n")

    test_cases = [
    ("chr17", 43094692,  "BRCA1"),
    ("chr13", 32339916,  "BRCA2"),
    ("chr17", 7674220,   "TP53"),
    ("chr22", 42126499,  "CYP2D6"),
    ("chr10", 94781859,  "CYP2C19"),
    ("chr7",  117480025, "CFTR"),
]

    for chrom, pos, label in test_cases:
        result = lookup_gene(chrom, pos)
        if result:
            print(f"  {label}")
            print(f"    {chrom}:{pos} → {result['gene_name']} "
                  f"({result['region']}, strand {result['strand']})")
        else:
            print(f"  {label}")
            print(f"    {chrom}:{pos} → No gene found")
        print()