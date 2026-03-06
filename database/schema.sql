-- =============================================================
-- GENOMICS DASHBOARD — DATABASE SCHEMA
-- =============================================================
-- Run this file with:
-- psql -h localhost -U genomics_user -d genomics_db -f schema.sql
-- =============================================================


-- =============================================================
-- SECTION 1: USERS
-- Who is using the system
-- =============================================================

CREATE TABLE IF NOT EXISTS users (
    user_id       SERIAL PRIMARY KEY,
    username      VARCHAR(100) NOT NULL UNIQUE,
    email         VARCHAR(255) NOT NULL UNIQUE,
    created_at    TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SECTION 2: VCF UPLOADS
-- One user can upload multiple VCF files over time
-- =============================================================

CREATE TABLE IF NOT EXISTS vcf_uploads (
    upload_id     SERIAL PRIMARY KEY,
    user_id       INT NOT NULL REFERENCES users(user_id) ON DELETE CASCADE,
    filename      VARCHAR(255) NOT NULL,
    uploaded_at   TIMESTAMP DEFAULT NOW(),
    status        VARCHAR(50) DEFAULT 'pending',
    -- status values: pending | processing | complete | failed
    total_variants INT DEFAULT 0,
    notes         TEXT
);


-- =============================================================
-- SECTION 3: VARIANTS
-- Every genetic variant found in a VCF file
-- This is the central table everything else connects to
-- =============================================================

CREATE TABLE IF NOT EXISTS variants (
    variant_id    SERIAL PRIMARY KEY,
    upload_id     INT NOT NULL REFERENCES vcf_uploads(upload_id) ON DELETE CASCADE,
    chromosome    VARCHAR(10)  NOT NULL,   -- e.g. 'chr17'
    position      BIGINT       NOT NULL,   -- genomic coordinate
    ref_allele    VARCHAR(500) NOT NULL,   -- reference base(s)
    alt_allele    VARCHAR(500) NOT NULL,   -- alternate base(s)
    variant_id_rs VARCHAR(50),             -- rsID from dbSNP e.g. 'rs1234567'
    gene_name     VARCHAR(100),            -- e.g. 'BRCA1', 'CYP2D6'
    zygosity      VARCHAR(20),             -- 'heterozygous' | 'homozygous_alt' | 'homozygous_ref'
    quality_score FLOAT,                   -- QUAL field from VCF
    depth         INT,                     -- read depth (DP)
    allele_freq   FLOAT,                   -- allele frequency (AF)
    flag          VARCHAR(50),             -- 'pharmacogenomics' | 'clinical' | 'both' | NULL
    risk_score    FLOAT                    -- normalized 0.0 to 1.0, filled in by modules
);

-- Index for fast lookups by gene (used constantly in queries)
CREATE INDEX IF NOT EXISTS idx_variants_gene
    ON variants(gene_name);

-- Index for fast lookups by upload (used to get all variants for a user)
CREATE INDEX IF NOT EXISTS idx_variants_upload
    ON variants(upload_id);

-- Index for chromosome + position (used in annotation matching)
CREATE INDEX IF NOT EXISTS idx_variants_position
    ON variants(chromosome, position);


-- =============================================================
-- SECTION 4: CLINVAR ANNOTATIONS
-- Loaded from ClinVar's free variant_summary.txt download
-- Maps known variants to diseases and clinical significance
-- =============================================================

CREATE TABLE IF NOT EXISTS clinvar_annotations (
    clinvar_id            SERIAL PRIMARY KEY,
    chromosome            VARCHAR(10),
    position              BIGINT,
    ref_allele            TEXT,
    alt_allele            TEXT,
    gene_name             TEXT,
    clinical_significance TEXT,  -- 'Pathogenic' | 'Likely pathogenic' | 'VUS' | 'Benign'
    condition_name        TEXT,          -- disease name e.g. 'Hereditary breast cancer'
    review_status         TEXT,  -- how well-evidenced: 'practice guideline' > 'reviewed by expert panel' etc.
    clinvar_accession     TEXT    -- ClinVar's own ID e.g. 'RCV000012345'
);

-- Index for fast annotation matching
CREATE INDEX IF NOT EXISTS idx_clinvar_position
    ON clinvar_annotations(chromosome, position);

CREATE INDEX IF NOT EXISTS idx_clinvar_gene
    ON clinvar_annotations(gene_name);

CREATE INDEX IF NOT EXISTS idx_clinvar_significance
    ON clinvar_annotations(clinical_significance);


-- =============================================================
-- SECTION 5: PHARMGKB ANNOTATIONS
-- Loaded from PharmGKB's free relationships download
-- Maps genes to drug metabolism effects
-- =============================================================

CREATE TABLE IF NOT EXISTS pharmgkb_annotations (
    pharmgkb_id      SERIAL PRIMARY KEY,
    gene_name        VARCHAR(100),   -- e.g. 'CYP2D6'
    drug_name        VARCHAR(255),   -- e.g. 'Codeine'
    effect_summary   TEXT,           -- e.g. 'Poor metabolizer — increased toxicity risk'
    evidence_level   VARCHAR(10),    -- PharmGKB scale: '1A' | '1B' | '2A' | '2B' | '3' | '4'
    phenotype_category VARCHAR(100), -- 'Metabolism/PK' | 'Efficacy' | 'Toxicity' | 'Dosage'
    pmid             VARCHAR(50)     -- PubMed ID of supporting study
);

CREATE INDEX IF NOT EXISTS idx_pharmgkb_gene
    ON pharmgkb_annotations(gene_name);

CREATE INDEX IF NOT EXISTS idx_pharmgkb_drug
    ON pharmgkb_annotations(drug_name);


-- =============================================================
-- SECTION 6: VARIANT ANNOTATIONS
-- The JOIN table — links YOUR variants to knowledge base findings
-- This is where ClinVar + PharmGKB results are attached to a user's variants
-- =============================================================

CREATE TABLE IF NOT EXISTS variant_annotations (
    annotation_id   SERIAL PRIMARY KEY,
    variant_id      INT NOT NULL REFERENCES variants(variant_id) ON DELETE CASCADE,
    source          VARCHAR(50)  NOT NULL,  -- 'clinvar' | 'pharmgkb'
    source_id       INT,                    -- ID from clinvar_annotations or pharmgkb_annotations
    annotation_type VARCHAR(50),            -- 'pathogenic' | 'likely_pathogenic' | 'VUS' | 'pharmacogenomics'
    risk_score      FLOAT,                  -- normalized 0.0 to 1.0
    notes           TEXT,
    created_at      TIMESTAMP DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_variant_annotations_variant
    ON variant_annotations(variant_id);

CREATE INDEX IF NOT EXISTS idx_variant_annotations_source
    ON variant_annotations(source, annotation_type);


-- =============================================================
-- SECTION 7: CRISPR SIMULATIONS
-- Stores results from the CRISPR off-target simulation module
-- Only created for pathogenic variants the user selects
-- =============================================================

CREATE TABLE IF NOT EXISTS crispr_simulations (
    simulation_id        SERIAL PRIMARY KEY,
    variant_id           INT NOT NULL REFERENCES variants(variant_id) ON DELETE CASCADE,
    guide_rna_seq        TEXT,    -- the designed gRNA sequence
    pam_sequence         VARCHAR(10),             -- PAM site e.g. 'NGG' for SpCas9
    on_target_efficiency FLOAT,                   -- 0.0 to 1.0
    off_target_score     FLOAT,                   -- 0.0 = safe, 1.0 = very risky
    off_target_sites     INT,                     -- number of predicted off-target locations
    safety_verdict       VARCHAR(50),             -- 'safe' | 'caution' | 'unsafe'
    simulated_at         TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SECTION 8: MICROBIOME SAMPLES
-- Stores gut microbiome diversity data (second input file)
-- Used to cross-reference with pharmacogenomics results
-- =============================================================

CREATE TABLE IF NOT EXISTS microbiome_samples (
    sample_id          SERIAL PRIMARY KEY,
    user_id            INT NOT NULL REFERENCES users(user_id) ON DELETE CASCADE,
    diversity_score    FLOAT,    -- Shannon diversity index (higher = more diverse = healthier)
    bacteroides_pct    FLOAT,    -- % Bacteroides (affects drug metabolism)
    firmicutes_pct     FLOAT,    -- % Firmicutes
    prevotella_pct     FLOAT,    -- % Prevotella
    other_pct          FLOAT,    -- % everything else
    metabolic_impact   VARCHAR(50),  -- 'low' | 'moderate' | 'high' impact on drug metabolism
    raw_summary        JSONB,    -- store full profile flexibly for future use
    uploaded_at        TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SECTION 9: RISK SCORES (DASHBOARD SUMMARY)
-- Stores the final normalized scores shown on the dashboard
-- One row per upload — the unified view across all modules
-- =============================================================

CREATE TABLE IF NOT EXISTS risk_summary (
    summary_id              SERIAL PRIMARY KEY,
    upload_id               INT NOT NULL REFERENCES vcf_uploads(upload_id) ON DELETE CASCADE,
    pathogenicity_score     FLOAT,   -- 0.0 to 1.0 — from ClinVar module
    pharmacogenomics_score  FLOAT,   -- 0.0 to 1.0 — from PharmGKB module
    crispr_safety_score     FLOAT,   -- 0.0 to 1.0 — from CRISPR module
    microbiome_score        FLOAT,   -- 0.0 to 1.0 — from microbiome module
    overall_score           FLOAT,   -- weighted average of all four
    generated_at            TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- QUICK SANITY CHECK — run after loading schema
-- Should print all 9 table names
-- =============================================================

SELECT table_name
FROM information_schema.tables
WHERE table_schema = 'public'
ORDER BY table_name;
