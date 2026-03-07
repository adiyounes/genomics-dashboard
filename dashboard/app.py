"""
dashboard/app.py
================
Genomics Dashboard — Streamlit Interface

Run with:
    cd ~/genomics-dashboard
    streamlit run dashboard/app.py
"""

import sys
import time
from pathlib import Path

import streamlit as st
import pandas as pd

sys.path.append(str(Path(__file__).parent.parent))

from modules.ingestion.vcf_parser  import ingest_vcf
from modules.annotation.annotator  import annotate_upload
from database.connect               import execute_query

# ── Page config ──────────────────────────────────────────────
st.set_page_config(
    page_title = "GenomeDx",
    page_icon  = "🧬",
    layout     = "wide",
    initial_sidebar_state = "expanded",
)

# ── Custom CSS ────────────────────────────────────────────────
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@300;400;500;600&display=swap');

html, body, [class*="css"] {
    font-family: 'IBM Plex Sans', sans-serif;
    background-color: #0a0e1a;
    color: #c8d6e5;
}
.main { background-color: #0a0e1a; }
.block-container { padding: 2rem 2.5rem; max-width: 1400px; }

.genome-header {
    display: flex;
    align-items: baseline;
    gap: 12px;
    margin-bottom: 2rem;
    padding-bottom: 1.5rem;
    border-bottom: 1px solid #1e2d42;
}
.genome-title {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 1.6rem;
    font-weight: 600;
    color: #e8f4fd;
}
.genome-subtitle {
    font-size: 0.85rem;
    color: #4a6080;
    letter-spacing: 0.08em;
    text-transform: uppercase;
}
.section-header {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.7rem;
    text-transform: uppercase;
    letter-spacing: 0.15em;
    color: #4a6080;
    margin-bottom: 1rem;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid #1e2d42;
}
.metric-card {
    background: #0f1729;
    border: 1px solid #1e2d42;
    border-radius: 6px;
    padding: 1.2rem 1.4rem;
}
.metric-label {
    font-size: 0.7rem;
    text-transform: uppercase;
    letter-spacing: 0.12em;
    color: #4a6080;
    margin-bottom: 0.5rem;
    font-family: 'IBM Plex Mono', monospace;
}
.metric-value {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 1.8rem;
    font-weight: 600;
    line-height: 1;
}
.metric-sub {
    font-size: 0.75rem;
    color: #4a6080;
    margin-top: 0.3rem;
}
.risk-bar-row {
    display: flex;
    align-items: center;
    margin-bottom: 0.8rem;
    gap: 1rem;
}
.risk-bar-name {
    font-size: 0.8rem;
    color: #7a90a8;
    min-width: 160px;
    font-family: 'IBM Plex Mono', monospace;
}
.risk-bar-track {
    flex: 1;
    background: #1e2d42;
    border-radius: 3px;
    height: 6px;
    overflow: hidden;
}
.risk-bar-fill {
    height: 100%;
    border-radius: 3px;
}
.risk-bar-score {
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.85rem;
    color: #e8f4fd;
    min-width: 45px;
    text-align: right;
}
section[data-testid="stSidebar"] {
    background: #0a0e1a;
    border-right: 1px solid #1e2d42;
}
.stButton > button {
    background: #1a3a5c;
    color: #e8f4fd;
    border: 1px solid #2563eb;
    border-radius: 4px;
    font-family: 'IBM Plex Mono', monospace;
    font-size: 0.8rem;
    padding: 0.5rem 1.2rem;
}
.stButton > button:hover { background: #2563eb; }
</style>
""", unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────────────────────

def get_risk_color(score):
    if score is None: return "#4a6080"
    if score >= 0.7:  return "#f87171"
    if score >= 0.4:  return "#fbbf24"
    if score >= 0.1:  return "#38bdf8"
    return "#4ade80"


def score_to_label(score):
    if score is None: return "—"
    if score >= 0.7:  return "HIGH"
    if score >= 0.4:  return "MODERATE"
    if score >= 0.1:  return "LOW"
    return "MINIMAL"


def get_all_uploads():
    return execute_query("""
        SELECT
            vu.upload_id,
            vu.filename,
            vu.total_variants,
            vu.status,
            vu.notes,
            u.username,
            u.email
        FROM vcf_uploads vu
        JOIN users u ON vu.user_id = u.user_id
        WHERE vu.status = 'complete'
        ORDER BY vu.upload_id DESC
    """)


def get_risk_summary(upload_id):
    result = execute_query("""
        SELECT * FROM risk_summary
        WHERE upload_id = %s
        ORDER BY summary_id DESC
        LIMIT 1
    """, params=(upload_id,))
    return result[0] if result else None


def get_variants_for_upload(upload_id):
    return execute_query("""
        SELECT
            v.variant_id,
            v.chromosome,
            v.position,
            v.ref_allele,
            v.alt_allele,
            v.gene_name,
            v.zygosity,
            v.flag,
            v.risk_score,
            v.quality_score,
            v.depth,
            (
                SELECT va.annotation_type
                FROM variant_annotations va
                WHERE va.variant_id = v.variant_id
                ORDER BY va.risk_score DESC
                LIMIT 1
            ) as top_annotation_type,
            (
                SELECT COUNT(*)
                FROM variant_annotations va
                WHERE va.variant_id = v.variant_id
            ) as annotation_count
        FROM variants v
        WHERE v.upload_id = %s
        ORDER BY v.risk_score DESC NULLS LAST
    """, params=(upload_id,))


def get_pharmgkb_for_upload(upload_id):
    return execute_query("""
        SELECT DISTINCT
            v.gene_name,
            v.zygosity,
            va.notes,
            va.risk_score
        FROM variants v
        JOIN variant_annotations va ON v.variant_id = va.variant_id
        WHERE v.upload_id = %s
        AND   va.source   = 'pharmgkb'
        ORDER BY va.risk_score DESC
    """, params=(upload_id,))


# ─────────────────────────────────────────────────────────────
# SIDEBAR
# ─────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("""
    <div style="margin-bottom:2rem;">
        <div style="font-family:'IBM Plex Mono',monospace;font-size:1.1rem;
                    font-weight:600;color:#e8f4fd;">🧬 GenomeDx</div>
        <div style="font-size:0.7rem;color:#4a6080;text-transform:uppercase;
                    letter-spacing:0.1em;margin-top:4px;">Genomic Analysis Platform</div>
    </div>
    """, unsafe_allow_html=True)

    page = st.radio(
        "Navigation",
        ["Upload & Analyse", "Browse Results", "Database Stats"],
        label_visibility="collapsed"
    )

    st.markdown("---")

    stats = execute_query("""
        SELECT
            (SELECT COUNT(*) FROM users)               as users,
            (SELECT COUNT(*) FROM vcf_uploads
             WHERE status='complete')                  as uploads,
            (SELECT COUNT(*) FROM variants)            as variants,
            (SELECT COUNT(*) FROM variant_annotations) as annotations
    """)
    if stats:
        s = stats[0]
        st.markdown(f"""
        <div style="font-family:'IBM Plex Mono',monospace;font-size:0.65rem;
                    text-transform:uppercase;letter-spacing:0.1em;color:#4a6080;
                    margin-bottom:0.8rem;">Database</div>
        <div style="font-size:0.8rem;color:#7a90a8;line-height:2.2;">
            {s['users']} users<br>
            {s['uploads']} uploads<br>
            {s['variants']:,} variants<br>
            {s['annotations']:,} annotations
        </div>
        """, unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────
# PAGE 1: UPLOAD & ANALYSE
# ─────────────────────────────────────────────────────────────

if page == "Upload & Analyse":

    st.markdown("""
    <div class="genome-header">
        <div class="genome-title">Upload & Analyse</div>
        <div class="genome-subtitle">VCF → Annotation → Risk Scores</div>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns([1, 1], gap="large")

    with col1:
        st.markdown('<div class="section-header">Patient Information</div>',
                    unsafe_allow_html=True)
        username = st.text_input("Username", placeholder="e.g. patient_001")
        email    = st.text_input("Email",    placeholder="e.g. patient@clinic.org")

        st.markdown('<div class="section-header" style="margin-top:1.5rem;">VCF File</div>',
                    unsafe_allow_html=True)
        uploaded_file = st.file_uploader(
            "Upload VCF", type=["vcf", "gz"],
            label_visibility="collapsed"
        )
        if uploaded_file:
            st.success(f"✓ {uploaded_file.name} ({uploaded_file.size/1024:.1f} KB)")

    with col2:
        st.markdown('<div class="section-header">Pipeline Steps</div>',
                    unsafe_allow_html=True)
        st.markdown("""
        <div style="font-size:0.82rem;color:#4a6080;line-height:2;
                    font-family:'IBM Plex Mono',monospace;">
            01 → Parse VCF file<br>
            02 → Detect assembly<br>
            03 → Extract variants<br>
            04 → Flag CYP / disease genes<br>
            05 → Match against ClinVar<br>
            06 → Match against PharmGKB<br>
            07 → Calculate risk scores<br>
            08 → Generate risk summary
        </div>
        """, unsafe_allow_html=True)

    st.markdown("<br>", unsafe_allow_html=True)

    if st.button("▶  Run Analysis"):
        if not uploaded_file:
            st.error("Please upload a VCF file first.")
        elif not username or not email:
            st.error("Please enter username and email.")
        else:
            tmp_path = Path(f"/tmp/{uploaded_file.name}")
            tmp_path.write_bytes(uploaded_file.getvalue())

            progress = st.progress(0, text="Starting pipeline...")
            status   = st.empty()

            try:
                status.info("🔬 Step 1/2 — Parsing VCF file...")
                progress.progress(20, text="Parsing VCF...")

                results = ingest_vcf(
                    filepath = tmp_path,
                    username = username,
                    email    = email,
                )

                if not results:
                    st.error("Failed to parse VCF file.")
                    st.stop()

                progress.progress(60, text="Annotating variants...")
                status.info("🧪 Step 2/2 — Running annotation engine...")

                for result in results:
                    annotate_upload(result['upload_id'])

                progress.progress(100, text="Complete!")
                status.empty()

                total_inserted = sum(r['inserted'] for r in results)
                st.success(
                    f"✅ Analysis complete — "
                    f"{total_inserted:,} variants ingested across "
                    f"{len(results)} sample(s)"
                )
                time.sleep(1)
                st.rerun()

            except Exception as e:
                st.error(f"Pipeline error: {e}")
            finally:
                if tmp_path.exists():
                    tmp_path.unlink()


# ─────────────────────────────────────────────────────────────
# PAGE 2: BROWSE RESULTS
# ─────────────────────────────────────────────────────────────

elif page == "Browse Results":

    st.markdown("""
    <div class="genome-header">
        <div class="genome-title">Browse Results</div>
        <div class="genome-subtitle">Variants · Annotations · Risk Scores</div>
    </div>
    """, unsafe_allow_html=True)

    uploads = get_all_uploads()

    if not uploads:
        st.info("No completed uploads found. Upload a VCF file first.")
        st.stop()

    # Upload selector
    upload_labels = {
        f"{u['username']} — {u['filename']} ({u['total_variants']} variants)": u['upload_id']
        for u in uploads
    }
    selected_label = st.selectbox(
        "Select upload",
        list(upload_labels.keys()),
        label_visibility="collapsed"
    )
    upload_id = upload_labels[selected_label]

    # Assembly badge
    selected = next(u for u in uploads if u['upload_id'] == upload_id)
    assembly = "unknown"
    if selected['notes']:
        for part in selected['notes'].split(";"):
            if "assembly=" in part:
                assembly = part.split("=")[1].strip()

    color = "#4ade80" if assembly == "GRCh38" else "#fbbf24" if assembly == "GRCh37" else "#f87171"
    st.markdown(
        f'<span style="font-family:IBM Plex Mono,monospace;font-size:0.8rem;'
        f'color:{color};">⬡ Assembly: {assembly}</span>',
        unsafe_allow_html=True
    )

    st.markdown("<br>", unsafe_allow_html=True)

    # ── Risk Summary Cards ──
    summary = get_risk_summary(upload_id)

    if summary:
        c1, c2, c3, c4 = st.columns(4)
        for col, label, score in [
            (c1, "Overall Risk",     summary['overall_score']),
            (c2, "Pathogenicity",    summary['pathogenicity_score']),
            (c3, "Pharmacogenomics", summary['pharmacogenomics_score']),
            (c4, "CRISPR/Microbiome", None),
        ]:
            rc = get_risk_color(score)
            with col:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-label">{label}</div>
                    <div class="metric-value" style="color:{rc};">
                        {f"{score:.2f}" if score is not None else "—"}
                    </div>
                    <div class="metric-sub">{score_to_label(score)}</div>
                </div>
                """, unsafe_allow_html=True)

        # Module score bars
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown('<div class="section-header">Module Scores</div>',
                    unsafe_allow_html=True)
        for label, score, color in [
            ("Pathogenicity",    summary['pathogenicity_score']    or 0, "#f87171"),
            ("Pharmacogenomics", summary['pharmacogenomics_score'] or 0, "#38bdf8"),
            ("CRISPR Safety",    0.0,                                    "#4a6080"),
            ("Microbiome",       0.0,                                    "#4a6080"),
        ]:
            pct = int((score or 0) * 100)
            st.markdown(f"""
            <div class="risk-bar-row">
                <div class="risk-bar-name">{label}</div>
                <div class="risk-bar-track">
                    <div class="risk-bar-fill"
                         style="width:{pct}%;background:{color};"></div>
                </div>
                <div class="risk-bar-score">{score:.3f}</div>
            </div>
            """, unsafe_allow_html=True)
    else:
        st.warning("No risk summary found. Run the annotation engine for this upload.")

    # ── Variants Table ──
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-header">Variants</div>', unsafe_allow_html=True)

    variants = get_variants_for_upload(upload_id)

    if variants:
        f1, f2, f3 = st.columns(3)
        with f1:
            flag_filter = st.multiselect(
                "Filter by flag",
                ["clinical", "pharmacogenomics", "both"],
                default=["clinical", "pharmacogenomics", "both"],
                placeholder="All flags..."
            )
        with f2:
            min_score = st.slider("Min risk score", 0.0, 1.0, 0.0, 0.05)
        with f3:
            show_unannotated = st.checkbox("Show unannotated variants", value=True)

        rows = []
        for v in variants:
            flag      = v['flag']      or "none"
            score     = v['risk_score']
            ann_type  = v['top_annotation_type'] or "—"
            ann_count = v['annotation_count']    or 0

            if flag_filter and flag not in flag_filter:
                continue
            if score is not None and score < min_score:
                continue
            if not show_unannotated and ann_count == 0:
                continue

            rows.append({
                "Chromosome"     : v['chromosome'],
                "Position"       : v['position'],
                "Allele"         : f"{v['ref_allele']} → {v['alt_allele']}",
                "Gene"           : v['gene_name'] or "—",
                "Zygosity"       : v['zygosity']  or "—",
                "Classification" : ann_type,
                "Risk Score"     : round(score, 3) if score is not None else None,
                "Annotations"    : ann_count,
            })

        if not rows:
            st.info("No variants match the current filters.")
        else:
            df = pd.DataFrame(rows)
            st.dataframe(df, use_container_width=True, hide_index=True)
            st.caption(f"{len(rows)} of {len(variants)} variants shown")

    else:
        st.info("No variants found for this upload.")

    # ── Drug Interactions ──
    pgx_data = get_pharmgkb_for_upload(upload_id)
    if pgx_data:
        st.markdown("<br>", unsafe_allow_html=True)
        st.markdown('<div class="section-header">Drug Interactions</div>',
                    unsafe_allow_html=True)

        drug_rows = []
        for row in pgx_data:
            parts = row['notes'].split(" | ")
            drug_rows.append({
                "Gene"     : row['gene_name'],
                "Drug"     : parts[1].replace("Drug: ", "")     if len(parts) > 1 else "—",
                "Effect"   : parts[0]                           if len(parts) > 0 else "—",
                "Evidence" : parts[2].replace("Evidence: ", "") if len(parts) > 2 else "—",
                "Category" : parts[3].replace("Category: ", "") if len(parts) > 3 else "—",
                "Score"    : round(row['risk_score'], 3),
            })

        df_pgx = pd.DataFrame(drug_rows)
        st.dataframe(df_pgx, use_container_width=True, hide_index=True)
        st.caption(f"{len(drug_rows)} drug interactions found")


# ─────────────────────────────────────────────────────────────
# PAGE 3: DATABASE STATS
# ─────────────────────────────────────────────────────────────

elif page == "Database Stats":

    st.markdown("""
    <div class="genome-header">
        <div class="genome-title">Database Stats</div>
        <div class="genome-subtitle">Knowledge Base · Pipeline Health</div>
    </div>
    """, unsafe_allow_html=True)

    col1, col2 = st.columns(2, gap="large")

    with col1:
        st.markdown('<div class="section-header">ClinVar Distribution</div>',
                    unsafe_allow_html=True)
        clinvar_dist = execute_query("""
            SELECT clinical_significance, COUNT(*) as count
            FROM clinvar_annotations
            GROUP BY clinical_significance
            ORDER BY count DESC
            LIMIT 8
        """)
        if clinvar_dist:
            df_cv = pd.DataFrame(clinvar_dist)
            df_cv.columns = ['Significance', 'Count']
            st.dataframe(df_cv, use_container_width=True, hide_index=True)

    with col2:
        st.markdown('<div class="section-header">PharmGKB Evidence Levels</div>',
                    unsafe_allow_html=True)
        pgx_dist = execute_query("""
            SELECT evidence_level, COUNT(*) as count
            FROM pharmgkb_annotations
            GROUP BY evidence_level
            ORDER BY evidence_level
        """)
        if pgx_dist:
            df_pgx2 = pd.DataFrame(pgx_dist)
            df_pgx2.columns = ['Evidence Level', 'Count']
            st.dataframe(df_pgx2, use_container_width=True, hide_index=True)

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<div class="section-header">Recent Uploads</div>',
                unsafe_allow_html=True)
    recent = execute_query("""
        SELECT
            vu.upload_id,
            u.username,
            vu.filename,
            vu.total_variants,
            vu.status
        FROM vcf_uploads vu
        JOIN users u ON vu.user_id = u.user_id
        ORDER BY vu.upload_id DESC
        LIMIT 10
    """)
    if recent:
        df_r = pd.DataFrame(recent)
        df_r.columns = ['ID', 'User', 'File', 'Variants', 'Status']
        st.dataframe(df_r, use_container_width=True, hide_index=True)

    st.markdown("""
    <div style="margin-top:3rem;padding:1rem;border:1px solid #1e2d42;
                border-radius:4px;font-size:0.75rem;color:#4a6080;line-height:1.6;">
        ⚠ This platform is intended for research purposes only.
        Results should not be used as a substitute for professional medical advice,
        diagnosis, or treatment.
    </div>
    """, unsafe_allow_html=True)