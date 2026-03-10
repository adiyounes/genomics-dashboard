-- =============================================================
-- GenomeDx — PostgreSQL Schema
-- Database : genomics_db
-- User     : genomics_user
-- Updated  : March 2026
--
-- Tables (11):
--   Core pipeline : users, vcf_uploads, variants
--   Annotations   : variant_annotations, risk_summary
--   Knowledge base: clinvar_annotations, pharmgkb_annotations
--   Reference     : gene_coordinates, exon_coordinates
--   Modules       : crispr_simulations, microbiome_samples
--
-- Load order matters — respect foreign key dependencies:
--   1. users
--   2. vcf_uploads
--   3. variants
--   4. variant_annotations, risk_summary, crispr_simulations
--   5. clinvar_annotations, pharmgkb_annotations (independent)
--   6. gene_coordinates
--   7. exon_coordinates
--   8. microbiome_samples
-- =============================================================


-- =============================================================
-- SECTION 1: USERS
-- Who is using the system
-- =============================================================

CREATE TABLE IF NOT EXISTS users (
    user_id    SERIAL PRIMARY KEY,
    username   VARCHAR(100) NOT NULL UNIQUE,
    email      VARCHAR(255) NOT NULL UNIQUE,
    created_at TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SECTION 2: VCF UPLOADS
-- One user can upload multiple VCF files over time
-- status: pending | processing | complete | failed
-- notes : stores assembly info e.g. "assembly=GRCh38"
-- =============================================================

CREATE TABLE IF NOT EXISTS vcf_uploads (
    upload_id       SERIAL PRIMARY KEY,
    user_id         INT NOT NULL REFERENCES users(user_id) ON DELETE CASCADE,
    filename        VARCHAR(255) NOT NULL,
    uploaded_at     TIMESTAMP DEFAULT NOW(),
    status          VARCHAR(50) DEFAULT 'pending',
    total_variants  INT DEFAULT 0,
    notes           TEXT
);


-- =============================================================
-- SECTION 3: VARIANTS
-- Every genetic variant found in a VCF file.
-- Central table — everything else connects to this.
--
-- flag values:
--   'clinical'         — variant is in a known disease gene
--   'pharmacogenomics' — variant is in a drug metabolism gene
--   'both'             — variant is in both
--   NULL               — variant is in an uncharacterised gene
--
-- zygosity values:
--   'heterozygous'     — one mutated copy (inherited from one parent)
--   'homozygous_alt'   — two mutated copies (inherited from both parents)
--   'homozygous_ref'   — no mutation (skipped during ingestion)
-- =============================================================

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

CREATE INDEX IF NOT EXISTS idx_variants_gene
    ON variants(gene_name);

CREATE INDEX IF NOT EXISTS idx_variants_upload
    ON variants(upload_id);

CREATE INDEX IF NOT EXISTS idx_variants_position
    ON variants(chromosome, position);


-- =============================================================
-- SECTION 4: CLINVAR ANNOTATIONS
-- Loaded from ClinVar variant_summary.txt (~2.6M rows)
-- Maps known variants to diseases and clinical significance
--
-- clinical_significance values:
--   'Pathogenic' | 'Likely pathogenic' | 'Uncertain significance'
--   'Likely benign' | 'Benign' | 'Conflicting...'
--
-- review_status (confidence):
--   'practice guideline' > 'reviewed by expert panel'
--   > 'criteria provided, multiple submitters'
--   > 'criteria provided, single submitter'
--   > 'no assertion criteria provided'
-- =============================================================

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

CREATE INDEX IF NOT EXISTS idx_clinvar_position
    ON clinvar_annotations(chromosome, position);

CREATE INDEX IF NOT EXISTS idx_clinvar_gene
    ON clinvar_annotations(gene_name);

CREATE INDEX IF NOT EXISTS idx_clinvar_significance
    ON clinvar_annotations(clinical_significance);


-- =============================================================
-- SECTION 5: PHARMGKB ANNOTATIONS
-- Loaded from PharmGKB relationships.tsv (~5,120 rows)
-- Maps genes to drug metabolism effects
--
-- evidence_level: '1A' (strongest) → '4' (weakest)
--   1A = clinical annotation with expert review
--   1B = clinical annotation
--   2A/2B = annotation with functional evidence
--   3/4   = case reports / preliminary evidence
-- =============================================================

CREATE TABLE IF NOT EXISTS pharmgkb_annotations (
    pharmgkb_id        SERIAL PRIMARY KEY,
    gene_name          TEXT,
    drug_name          TEXT,
    effect_summary     TEXT,
    evidence_level     TEXT,
    phenotype_category TEXT,
    pmid               TEXT
);

CREATE INDEX IF NOT EXISTS idx_pharmgkb_gene
    ON pharmgkb_annotations(gene_name);

CREATE INDEX IF NOT EXISTS idx_pharmgkb_drug
    ON pharmgkb_annotations(drug_name);


-- =============================================================
-- SECTION 6: VARIANT ANNOTATIONS
-- Links a patient's variants to knowledge base findings.
-- One variant can have many annotations (multiple ClinVar
-- submitters, multiple drugs per gene).
--
-- source values: 'clinvar' | 'pharmgkb'
-- source_id    : FK to clinvar_id or pharmgkb_id (soft reference)
--
-- annotation_type values:
--   'pathogenic' | 'likely_pathogenic' | 'VUS'
--   'likely_benign' | 'benign' | 'pharmacogenomics'
-- =============================================================

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

CREATE INDEX IF NOT EXISTS idx_variant_annotations_variant
    ON variant_annotations(variant_id);

CREATE INDEX IF NOT EXISTS idx_variant_annotations_source
    ON variant_annotations(source, annotation_type);


-- =============================================================
-- SECTION 7: RISK SUMMARY
-- One row per upload — the final unified risk scores shown
-- on the dashboard. Generated by the annotation engine.
--
-- Score weights (overall_score formula):
--   pathogenicity_score    × 0.50
--   pharmacogenomics_score × 0.30
--   crispr_safety_score    × 0.10
--   microbiome_score       × 0.10
-- =============================================================

CREATE TABLE IF NOT EXISTS risk_summary (
    summary_id             SERIAL PRIMARY KEY,
    upload_id              INT NOT NULL REFERENCES vcf_uploads(upload_id) ON DELETE CASCADE,
    pathogenicity_score    DOUBLE PRECISION,
    pharmacogenomics_score DOUBLE PRECISION,
    crispr_safety_score    DOUBLE PRECISION,
    microbiome_score       DOUBLE PRECISION,
    overall_score          DOUBLE PRECISION,
    generated_at           TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SECTION 8: GENE COORDINATES
-- Loaded from Ensembl GRCh38 GTF (~40,502 rows)
-- Maps every human gene to its chromosomal address.
-- Used by the VCF parser to resolve gene names from coordinates
-- when no Gene= tag is present in the VCF INFO field.
--
-- biotype examples:
--   'protein_coding' | 'lncRNA' | 'pseudogene' | 'miRNA'
-- strand: '+' (forward) | '-' (reverse)
-- =============================================================

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

CREATE INDEX IF NOT EXISTS idx_gene_coord_lookup
    ON gene_coordinates(chromosome, start_pos, end_pos);

CREATE INDEX IF NOT EXISTS idx_gene_coord_name
    ON gene_coordinates(gene_name);


-- =============================================================
-- SECTION 9: EXON COORDINATES
-- Loaded from Ensembl GRCh38 GTF (~1,552,149 rows)
-- One row per exon, linked to its parent gene.
-- Used to determine whether a variant falls in a coding
-- region (exonic) vs non-coding region (intronic).
-- Exonic variants carry higher clinical significance.
-- =============================================================

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

CREATE INDEX IF NOT EXISTS idx_exon_coord_lookup
    ON exon_coordinates(chromosome, start_pos, end_pos);

CREATE INDEX IF NOT EXISTS idx_exon_gene_name
    ON exon_coordinates(gene_name);


-- =============================================================
-- SECTION 10: CRISPR SIMULATIONS
-- Stores results from the CRISPR off-target prediction module.
-- Linked to a specific variant (the target site).
-- Built in Month 2.
--
-- safety_verdict: 'safe' | 'caution' | 'unsafe'
-- on_target_efficiency: 0.0 (won't cut) → 1.0 (cuts perfectly)
-- off_target_score    : 0.0 (safe) → 1.0 (very risky)
-- =============================================================

CREATE TABLE IF NOT EXISTS crispr_simulations (
    simulation_id        SERIAL PRIMARY KEY,
    variant_id           INT NOT NULL REFERENCES variants(variant_id) ON DELETE CASCADE,
    guide_rna_seq        TEXT,
    pam_sequence         VARCHAR(10),
    on_target_efficiency DOUBLE PRECISION,
    off_target_score     DOUBLE PRECISION,
    off_target_sites     INT,
    safety_verdict       VARCHAR(50),
    simulated_at         TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SECTION 11: MICROBIOME SAMPLES
-- Stores gut microbiome composition data.
-- Cross-referenced with pharmacogenomics results to refine
-- drug metabolism predictions. Built in Month 2.
--
-- diversity_score : Shannon diversity index
--   higher = more diverse = generally healthier gut
-- metabolic_impact: how much the microbiome affects drug metabolism
--   'low' | 'moderate' | 'high'
-- =============================================================

CREATE TABLE IF NOT EXISTS microbiome_samples (
    sample_id        SERIAL PRIMARY KEY,
    user_id          INT NOT NULL REFERENCES users(user_id) ON DELETE CASCADE,
    diversity_score  DOUBLE PRECISION,
    bacteroides_pct  DOUBLE PRECISION,
    firmicutes_pct   DOUBLE PRECISION,
    prevotella_pct   DOUBLE PRECISION,
    other_pct        DOUBLE PRECISION,
    metabolic_impact VARCHAR(50),
    raw_summary      JSONB,
    uploaded_at      TIMESTAMP DEFAULT NOW()
);


-- =============================================================
-- SANITY CHECK — run after loading schema
-- Should print all 11 table names
-- =============================================================

SELECT table_name
FROM information_schema.tables
WHERE table_schema = 'public'
ORDER BY table_name;