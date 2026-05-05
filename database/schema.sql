-- GenomeDx database schema
-- 9 tables, create in order: users, vcf_uploads, variants, variant_annotations, risk_summary
-- clinvar_annotations, pharmgkb_annotations, gene_coordinates, exon_coordinates can be loaded in any order


-- 1. Who uploads files to the system
CREATE TABLE IF NOT EXISTS users (
    user_id    SERIAL PRIMARY KEY,
    username   VARCHAR(100) NOT NULL UNIQUE,
    email      VARCHAR(255) NOT NULL UNIQUE,
    created_at TIMESTAMP DEFAULT NOW()
);


-- 2. One row per sample in an uploaded VCF file
-- status: pending | processing | complete | failed
CREATE TABLE IF NOT EXISTS vcf_uploads (
    upload_id       SERIAL PRIMARY KEY,
    user_id         INT NOT NULL REFERENCES users(user_id) ON DELETE CASCADE,
    filename        VARCHAR(255) NOT NULL,
    uploaded_at     TIMESTAMP DEFAULT NOW(),
    status          VARCHAR(50) DEFAULT 'pending',
    total_variants  INT DEFAULT 0,
    notes           TEXT
);


-- 3. Every genetic variant found in a VCF file
-- flag: clinical | pharmacogenomics | both | NULL
-- zygosity: heterozygous | homozygous_alt
CREATE TABLE IF NOT EXISTS variants (
    variant_id    SERIAL PRIMARY KEY,
    upload_id     INT NOT NULL REFERENCES vcf_uploads(upload_id) ON DELETE CASCADE,
    chromosome    VARCHAR(10)  NOT NULL,
    position      BIGINT       NOT NULL,
    ref_allele    VARCHAR(500) NOT NULL,
    alt_allele    VARCHAR(500) NOT NULL,
    variant_id_rs VARCHAR(50),
    gene_name     VARCHAR(100),
    zygosity      VARCHAR(20),
    quality_score DOUBLE PRECISION,
    depth         INT,
    allele_freq   DOUBLE PRECISION,
    flag          VARCHAR(50),
    risk_score    DOUBLE PRECISION
);

CREATE INDEX IF NOT EXISTS idx_variants_gene     ON variants(gene_name);
CREATE INDEX IF NOT EXISTS idx_variants_upload   ON variants(upload_id);
CREATE INDEX IF NOT EXISTS idx_variants_position ON variants(chromosome, position);


-- 4. 2.6M known disease variants loaded from NIH ClinVar
CREATE TABLE IF NOT EXISTS clinvar_annotations (
    clinvar_id            SERIAL PRIMARY KEY,
    chromosome            VARCHAR(10),
    position              BIGINT,
    ref_allele            TEXT,
    alt_allele            TEXT,
    gene_name             TEXT,
    clinical_significance TEXT,
    condition_name        TEXT,
    review_status         TEXT,
    clinvar_accession     TEXT
);

CREATE INDEX IF NOT EXISTS idx_clinvar_position     ON clinvar_annotations(chromosome, position);
CREATE INDEX IF NOT EXISTS idx_clinvar_gene         ON clinvar_annotations(gene_name);
CREATE INDEX IF NOT EXISTS idx_clinvar_significance ON clinvar_annotations(clinical_significance);


-- 5. 5,120 gene-drug relationships loaded from Stanford PharmGKB
-- evidence_level: 1A (strongest) to 4 (weakest)
CREATE TABLE IF NOT EXISTS pharmgkb_annotations (
    pharmgkb_id        SERIAL PRIMARY KEY,
    gene_name          TEXT,
    drug_name          TEXT,
    effect_summary     TEXT,
    evidence_level     TEXT,
    phenotype_category TEXT,
    pmid               TEXT
);

CREATE INDEX IF NOT EXISTS idx_pharmgkb_gene ON pharmgkb_annotations(gene_name);
CREATE INDEX IF NOT EXISTS idx_pharmgkb_drug ON pharmgkb_annotations(drug_name);


-- 6. Results of matching a patient's variants against ClinVar and PharmGKB
-- source_id is a soft reference to either clinvar_id or pharmgkb_id
CREATE TABLE IF NOT EXISTS variant_annotations (
    annotation_id   SERIAL PRIMARY KEY,
    variant_id      INT NOT NULL REFERENCES variants(variant_id) ON DELETE CASCADE,
    source          VARCHAR(50) NOT NULL,
    source_id       INT,
    annotation_type VARCHAR(50),
    risk_score      DOUBLE PRECISION,
    notes           TEXT,
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_variant_annotations_variant ON variant_annotations(variant_id);
CREATE INDEX IF NOT EXISTS idx_variant_annotations_source  ON variant_annotations(source, annotation_type);


-- 7. One risk summary per upload shown on the dashboard
-- overall_score = pathogenicity x 0.6 + pharmacogenomics x 0.4
CREATE TABLE IF NOT EXISTS risk_summary (
    summary_id             SERIAL PRIMARY KEY,
    upload_id              INT NOT NULL REFERENCES vcf_uploads(upload_id) ON DELETE CASCADE,
    pathogenicity_score    DOUBLE PRECISION,
    pharmacogenomics_score DOUBLE PRECISION,
    overall_score          DOUBLE PRECISION,
    generated_at           TIMESTAMP DEFAULT NOW()
);


-- 8. 40,502 human genes from Ensembl GRCh38, used to resolve gene names from coordinates
CREATE TABLE IF NOT EXISTS gene_coordinates (
    gene_coord_id SERIAL PRIMARY KEY,
    gene_name     TEXT NOT NULL,
    gene_id       TEXT,
    chromosome    TEXT NOT NULL,
    start_pos     BIGINT NOT NULL,
    end_pos       BIGINT NOT NULL,
    strand        CHAR(1),
    biotype       TEXT
);

CREATE INDEX IF NOT EXISTS idx_gene_coord_lookup ON gene_coordinates(chromosome, start_pos, end_pos);
CREATE INDEX IF NOT EXISTS idx_gene_coord_name   ON gene_coordinates(gene_name);


-- 9. 1,552,149 exons from Ensembl GRCh38, used to distinguish exonic from intronic variants
CREATE TABLE IF NOT EXISTS exon_coordinates (
    exon_coord_id SERIAL PRIMARY KEY,
    gene_coord_id INT REFERENCES gene_coordinates(gene_coord_id),
    gene_name     TEXT NOT NULL,
    chromosome    TEXT NOT NULL,
    start_pos     BIGINT NOT NULL,
    end_pos       BIGINT NOT NULL,
    strand        CHAR(1),
    exon_number   INT,
    transcript_id TEXT
);

CREATE INDEX IF NOT EXISTS idx_exon_coord_lookup ON exon_coordinates(chromosome, start_pos, end_pos);
CREATE INDEX IF NOT EXISTS idx_exon_gene_name    ON exon_coordinates(gene_name);


SELECT table_name
FROM information_schema.tables
WHERE table_schema = 'public'
ORDER BY table_name;