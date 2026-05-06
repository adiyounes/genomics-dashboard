import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent.parent.parent))

from database.connect import execute_query, execute_insert, get_connection


# Base risk scores by clinical significance
SIGNIFICANCE_SCORES = {
    "pathogenic"                        : 1.0,
    "pathogenic/likely pathogenic"      : 0.9,
    "pathogenic, low penetrance"        : 0.7,
    "likely pathogenic"                 : 0.8,
    "uncertain significance"            : 0.5,
    "likely benign"                     : 0.1,
    "benign"                            : 0.0,
}

# Zygosity multipliers
ZYGOSITY_MULTIPLIERS = {
    "homozygous_alt" : 1.0,   # both copies affected, full risk
    "heterozygous"   : 0.7,   # one working copy, reduced risk
    "unknown"        : 0.5,   # uncertain
}

# PharmGKB evidence level scores 
EVIDENCE_SCORES = {
    "1A": 0.9,
    "1B": 0.75,
    "2A": 0.6,
    "2B": 0.5,
    "3" : 0.4,
    "4" : 0.3,
}



def calculate_clinvar_score(clinical_significance, zygosity):
   
    if not clinical_significance:
        return 0.0
    
    multiplier = ZYGOSITY_MULTIPLIERS.get((zygosity or "unknown").lower(), 0.5)
    base_score = SIGNIFICANCE_SCORES.get(clinical_significance.lower().strip(), 0.1)
    return round(min(max(base_score * multiplier, 0.0), 1.0), 3)

def calculate_pharmgkb_score(evidence_level, zygosity):

    if not evidence_level:
        return 0
    
    base_score = EVIDENCE_SCORES.get(evidence_level.strip(), 0.1)
    multiplier = ZYGOSITY_MULTIPLIERS.get((zygosity or "unknown").lower(), 0.5)
    return round(min(max(base_score * multiplier, 0.0), 1.0), 3)

def match_clinvar(variant):
    annotations = []

    exact_match = execute_query(
        """
            SELECT
                clinvar_id,
                clinical_significance,
                condition_name,
                review_status,
                gene_name
            FROM clinvar_annotations
            WHERE chromosome = %s
            AND   position = %s
            AND   ref_allele = %s
            AND   alt_allele = %s
        """,(variant["chromosome"],
            variant["position"],
            variant["ref_allele"],
            variant["alt_allele"],
        )
    )

    for match in exact_match:
        score =calculate_clinvar_score(match["clinical_significance"],
                variant["zygosity"])
        annotations.append({
            "source"          : "clinvar",
            "source_id"       : match['clinvar_id'],
            "annotation_type" : classify_significance(match['clinical_significance']),
            "risk_score"      : score,
            "match_type"      : "exact",
            "notes"           : (
                f"{match['clinical_significance']} | "
                f"{match['condition_name'] or 'Unknown condition'} | "
                f"Review: {match['review_status'] or 'N/A'}"
            ),
        })

    if not annotations and variant.get('chromosome') and variant.get('position'):
        position_matches = execute_query("""
            SELECT
                clinvar_id,
                clinical_significance,
                condition_name,
                review_status,
                gene_name
            FROM clinvar_annotations
            WHERE chromosome = %s
            AND   position   = %s
            LIMIT 5
        """, params=(variant['chromosome'], variant['position']))

        for match in position_matches:
            score = calculate_clinvar_score(
                match['clinical_significance'],
                variant['zygosity']
            ) * 0.8   # reduce score for position only match
            annotations.append({
                "source"          : "clinvar",
                "source_id"       : match['clinvar_id'],
                "annotation_type" : classify_significance(match['clinical_significance']),
                "risk_score"      : score,
                "match_type"      : "position",
                "notes"           : (
                    f"{match['clinical_significance']} | "
                    f"{match['condition_name'] or 'Unknown condition'} | "
                    f"Position match only"
                ),
            })

    if not annotations and variant.get('gene_name'):
        gene_matches = execute_query("""
            SELECT
                clinvar_id,
                clinical_significance,
                condition_name,
                review_status
            FROM clinvar_annotations
            WHERE gene_name = %s
            AND   clinical_significance IN (
                'Pathogenic',
                'Likely pathogenic',
                'Pathogenic/Likely pathogenic'
            )
            ORDER BY
                CASE clinical_significance
                    WHEN 'Pathogenic' THEN 1
                    WHEN 'Pathogenic/Likely pathogenic' THEN 2
                    WHEN 'Likely pathogenic' THEN 3
                END
            LIMIT 3
        """, params=(variant['gene_name'],))

        for match in gene_matches:
            score = calculate_clinvar_score(
                match['clinical_significance'],
                variant['zygosity']
            ) * 0.5 # reduce score for gene name only match
            annotations.append({
                "source"          : "clinvar",
                "source_id"       : match['clinvar_id'],
                "annotation_type" : classify_significance(match['clinical_significance']),
                "risk_score"      : score,
                "match_type"      : "gene",
                "notes"           : (
                    f"{match['clinical_significance']} | "
                    f"{match['condition_name'] or 'Unknown condition'} | "
                    f"Gene-level match only"
                ),
            })

    return annotations

def classify_significance(significance):
    if not significance:
        return "unknown"
    sig = significance.lower()
    if "pathogenic" in sig and "likely" not in sig:
        return "pathogenic"
    if "likely pathogenic" in sig:
        return "likely_pathogenic"
    if "uncertain" in sig or "vus" in sig:
        return "vus"
    if "benign" in sig:
        return "benign"
    return "other"

def match_pharmgkb(variant):

    if not variant.get('gene_name'):
        return []

    matches = execute_query("""
        SELECT
            pharmgkb_id,
            drug_name,
            effect_summary,
            evidence_level,
            phenotype_category
        FROM pharmgkb_annotations
        WHERE gene_name = %s
        ORDER BY
            CASE evidence_level
                WHEN '1A' THEN 1
                WHEN '1B' THEN 2
                WHEN '2A' THEN 3
                WHEN '2B' THEN 4
                WHEN '3'  THEN 5
                WHEN '4'  THEN 6
                ELSE 7
            END
        LIMIT 10
    """, params=(variant['gene_name'],))

    annotations = []
    for match in matches:
        score = calculate_pharmgkb_score(
            match['evidence_level'],
            variant['zygosity']
        )
        annotations.append({
            "source"          : "pharmgkb",
            "source_id"       : match['pharmgkb_id'],
            "annotation_type" : "pharmacogenomics",
            "risk_score"      : score,
            "match_type"      : "gene",
            "notes"           : (
                f"{match['effect_summary']} | "
                f"Drug: {match['drug_name']} | "
                f"Evidence: {match['evidence_level']} | "
                f"Category: {match['phenotype_category']}"
            ),
        })

    return annotations

def save_annotations(variant_id, annotations):
    if not annotations:
        return 0.0

    max_score = 0.0

    for ann in annotations:
        execute_insert("""
            INSERT INTO variant_annotations (
                variant_id,
                source,
                source_id,
                annotation_type,
                risk_score,
                notes
            ) VALUES (%s, %s, %s, %s, %s, %s)
            RETURNING annotation_id
        """, params=(
            variant_id,
            ann['source'],
            ann['source_id'],
            ann['annotation_type'],
            ann['risk_score'],
            ann['notes'],
        ))

        if ann['risk_score'] > max_score:
            max_score = ann['risk_score']

    return max_score

def update_variant_risk_score(variant_id, risk_score):
    execute_insert("""
        UPDATE variants
        SET risk_score = %s
        WHERE variant_id = %s
    """, params=(risk_score, variant_id))

def save_risk_summary(upload_id):
    # Get all annotated variants for this upload
    variants = execute_query("""
        SELECT
            v.variant_id,
            v.flag,
            v.risk_score,
            v.zygosity
        FROM variants v
        WHERE v.upload_id = %s
        AND   v.risk_score IS NOT NULL
    """, params=(upload_id,))

    if not variants:
        return None

    # Separate scores by module
    pathogenicity_scores    = []
    pharmacogenomics_scores = []

    for v in variants:
        if v['flag'] in ('clinical', 'both') and v['risk_score']:
            pathogenicity_scores.append(v['risk_score'])
        if v['flag'] in ('pharmacogenomics', 'both') and v['risk_score']:
            pharmacogenomics_scores.append(v['risk_score'])

    # Calculate module scores — use max score per module
    # (worst case is most clinically relevant)
    path_score  = max(pathogenicity_scores)    if pathogenicity_scores    else 0.0
    pgx_score   = max(pharmacogenomics_scores) if pharmacogenomics_scores else 0.0

    # Overall = weighted average of available scores
    # Pathogenicity weighted higher (more clinically urgent)
    weights = [0.6, 0.4]
    scores  = [path_score, pgx_score]
    overall = round(sum(w * s for w, s in zip(weights, scores)), 3)

    # Delete existing summary for this upload if re-running
    execute_insert(
        "DELETE FROM risk_summary WHERE upload_id = %s",
        params=(upload_id,)
    )

    # Save new summary
    summary_id = execute_insert("""
        INSERT INTO risk_summary (
            upload_id,
            pathogenicity_score,
            pharmacogenomics_score,
            overall_score
        ) VALUES (%s, %s, %s, %s)
        RETURNING summary_id
    """, params=(
        upload_id,
        round(path_score, 3),
        round(pgx_score, 3),
        overall,
    ))

    return {
        "pathogenicity"    : round(path_score, 3),
        "pharmacogenomics" : round(pgx_score, 3),
        "overall"          : overall,
    }

def annotate_upload(upload_id):
    print(f"\n{'='*60}")
    print(f"  ANNOTATING upload_id = {upload_id}")
    print(f"{'='*60}")

    # Load variants
    variants = execute_query("""
        SELECT
            variant_id, chromosome, position,
            ref_allele, alt_allele, gene_name,
            zygosity, flag
        FROM variants
        WHERE upload_id = %s
        ORDER BY chromosome, position
    """, params=(upload_id,))

    if not variants:
        print(f"  ❌ No variants found for upload_id={upload_id}")
        return None

    print(f"\n  Found {len(variants)} variants to annotate\n")

    # Annotate each variant
    stats = {
        "total"          : len(variants),
        "clinvar_hits"   : 0,
        "pharmgkb_hits"  : 0,
        "no_match"       : 0,
        "annotations"    : 0,
    }

    for i, variant in enumerate(variants, 1):
        variant_id = variant['variant_id']
        gene       = variant['gene_name'] or "unknown"
        flag       = variant['flag']      or "none"

        print(f"  [{i}/{len(variants)}] {variant['chromosome']}:"
              f"{variant['position']} {gene} ({flag})")

        all_annotations = []

        # ClinVar matching
        if flag in ('clinical', 'both', None):
            clinvar_hits = match_clinvar(variant)
            if clinvar_hits:
                all_annotations.extend(clinvar_hits)
                stats['clinvar_hits'] += 1
                best = max(clinvar_hits, key=lambda x: x['risk_score'])
                print(f"    ✅ ClinVar: {best['annotation_type']} "
                      f"(score={best['risk_score']}, "
                      f"match={best['match_type']})")

        # PharmGKB matching 
        if flag in ('pharmacogenomics', 'both'):
            pgx_hits = match_pharmgkb(variant)
            if pgx_hits:
                all_annotations.extend(pgx_hits)
                stats['pharmgkb_hits'] += 1
                print(f"    💊 PharmGKB: {len(pgx_hits)} drug interactions found")
                for hit in pgx_hits[:3]:   # show top 3
                    print(f"       → {hit['notes'][:80]}")

        # Also run ClinVar on pharmacogenomics variants
        # CYP genes can also appear in ClinVar
        if flag == 'pharmacogenomics':
            clinvar_hits = match_clinvar(variant)
            if clinvar_hits:
                all_annotations.extend(clinvar_hits)
                stats['clinvar_hits'] += 1

        # Save annotations
        if all_annotations:
            max_score = save_annotations(variant_id, all_annotations)
            update_variant_risk_score(variant_id, max_score)
            stats['annotations'] += len(all_annotations)
        else:
            stats['no_match'] += 1
            print(f"    ○ No match found")

    #Calculate risk summary
    print(f"\n  Calculating risk summary...")
    summary = save_risk_summary(upload_id)

    # Print results
    print(f"  ANNOTATION COMPLETE")
    print(f"\n  Variants processed  : {stats['total']}")
    print(f"  ClinVar matches     : {stats['clinvar_hits']}")
    print(f"  PharmGKB matches    : {stats['pharmgkb_hits']}")
    print(f"  No match            : {stats['no_match']}")
    print(f"  Total annotations   : {stats['annotations']}")

    if summary:
        print(f"\n  Risk Summary ")
        print(f"  Pathogenicity score    : {summary['pathogenicity']}")
        print(f"  Pharmacogenomics score : {summary['pharmacogenomics']}")
        print(f"  Overall risk score     : {summary['overall']}")

    return {**stats, "summary": summary}

if __name__ == "__main__":

    # Get only uploads that haven't been annotated yet
    unannotated = execute_query("""
        SELECT
            vu.upload_id,
            vu.filename,
            vu.total_variants
        FROM vcf_uploads vu
        WHERE vu.status = 'complete'
        AND NOT EXISTS (
            SELECT 1
            FROM variants v
            JOIN variant_annotations va ON v.variant_id = va.variant_id
            WHERE v.upload_id = vu.upload_id
        )
        ORDER BY vu.upload_id
    """)

    if not unannotated:
        print("✅ All uploads are already annotated.")
        sys.exit(0)

    print(f"\nFound {len(unannotated)} unannotated upload(s):\n")
    for row in unannotated:
        print(f"  upload_id={row['upload_id']} | "
              f"{row['filename']} | "
              f"{row['total_variants']} variants")

    for row in unannotated:
        annotate_upload(row['upload_id'])