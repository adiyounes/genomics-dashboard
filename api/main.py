from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
import sys
from pathlib import Path
import httpx
from pydantic import BaseModel

sys.path.append(str(Path(__file__).parent.parent))

from database.connect import execute_query, execute_insert
from modules.ingestion.vcf_parser import ingest_vcf
from modules.annotation.annotator import annotate_upload

app = FastAPI(title="GenomeDX API")

AGENT_URL = "http://34.243.52.170:8000"

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000", 
                   "http://34.243.52.170:8000", 
                   "https://genomics-dashboard-ggin.vercel.app",
                   "genomics-dashboard-ggin-5h64x79xg-adiyounes-projects.vercel.app"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/health")
def helth():
    return {"status": "ok"}

@app.get("/uploads")
def get_uploads():
    uploads = execute_query(
        """
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
    """
    )
    return uploads

@app.get("/variants/{upload_id}")
def get_variants(upload_id: int):
    variants = execute_query(
        """
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
            v.depth
        FROM variants v
        WHERE v.upload_id = %s
    """
    , (upload_id,))
    return variants

@app.get("/risksummary/{upload_id}")
def get_risk_summary(upload_id: int):
    result = execute_query(
        """
        SELECT * FROM risk_summary
        WHERE upload_id = %s
        ORDER BY summary_id DESC
        LIMIT 1
    """,(upload_id,)
    )
    return result

@app.get("/drugs/{upload_id}")
def get_drugs(upload_id: int):
    result = execute_query(
        """
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
    """,(upload_id,))
    return result

@app.get("/stats")
def get_stats():
    result = execute_query(
        """
        SELECT
            (SELECT COUNT(*) FROM clinvar_annotations)  as clinvar_count,
            (SELECT COUNT(*) FROM pharmgkb_annotations) as pharmgkb_count,
            (SELECT COUNT(*) FROM variants)             as variant_count,
            (SELECT COUNT(*) FROM vcf_uploads
             WHERE status = 'complete')                 as upload_count
    """)
    return result[0] if result else {}

@app.post("/upload")
async def upload_vcf(
    file: UploadFile = File(...),
    username: str = Form(...),
    email: str = Form(...)
):
    temp_path = Path(f"/tmp/{file.filename}")

    try:
        contents = await file.read()
        temp_path.write_bytes(contents)

        results = ingest_vcf(
            filepath= temp_path,
            username= username,
            email= email
        )

        if not results:
            return {"error": "Failed tp parse vcf file"}
        
        for sample_name, upload_id in results['upload_ids'].items():
            annotate_upload(upload_id)

        return {
            "success"    : True,
            "samples"    : results['samples'],
            "upload_ids" : results['upload_ids'],
            "inserted"   : results['inserted']
        }
    except Exception as e:
        return {"error": str(e)}
    
    finally:
        if temp_path.exists():
            temp_path.unlink

@app.get("/stats/genes")
def get_top_genes():
    return execute_query("""
        SELECT gene_name, COUNT(*) as count
        FROM variants
        WHERE gene_name IS NOT NULL
        AND flag IS NOT NULL
        GROUP BY gene_name
        ORDER BY count DESC
        LIMIT 5
    """)

@app.get("/stats/risk-distribution")
def get_risk_distribution():
    return execute_query("""
        SELECT
            COUNT(CASE WHEN risk_score >= 0.7 THEN 1 END) as high,
            COUNT(CASE WHEN risk_score >= 0.4 AND risk_score < 0.7 THEN 1 END) as moderate,
            COUNT(CASE WHEN risk_score >= 0.1 AND risk_score < 0.4 THEN 1 END) as low,
            COUNT(CASE WHEN risk_score < 0.1 THEN 1 END) as minimal,
            COUNT(CASE WHEN risk_score IS NULL THEN 1 END) as unannotated
        FROM variants
    """)

@app.get("/stats/reference")
def get_reference_stats():
    return execute_query("""
        SELECT
            (SELECT COUNT(*) FROM gene_coordinates) as genes,
            (SELECT COUNT(*) FROM exon_coordinates) as exons
    """)

class ResearchRequest(BaseModel):
    gene: str
    condition: str

@app.post("/research")
async def get_research(request: ResearchRequest):
    try:
        async with httpx.AsyncClient(timeout=120.0) as client:
            response = await client.post(
                f"{AGENT_URL}/analyse",
                json={"gene": request.gene, "disease": request.condition}
            )
            if response.status_code == 200:
                return response.json()
            return {"error": "Agent unavailable"}
    except Exception as e:
        return {"error": str(e)}
    
@app.get("/conditions/{variant_id}")
def get_conditions(variant_id: int):
    return execute_query("""
        SELECT 
            condition_name,
            annotation_type,
            risk_score
        FROM (
            SELECT
                SPLIT_PART(notes, ' | ', 2) as condition_name,
                annotation_type,
                risk_score,
                ROW_NUMBER() OVER (
                    PARTITION BY SPLIT_PART(notes, ' | ', 2) 
                    ORDER BY risk_score DESC
                ) as rn
            FROM variant_annotations
            WHERE variant_id = %s
            AND source = 'clinvar'
            AND notes IS NOT NULL
            AND SPLIT_PART(notes, ' | ', 2) != 'not provided'
            AND SPLIT_PART(notes, ' | ', 2) != ''
        ) sub
        WHERE rn = 1
        ORDER BY risk_score DESC
        LIMIT 5
    """, (variant_id,))

@app.get("/debug-db")
def debug_db():
    import os
    return {
        "host": os.getenv("DB_HOST"),
        "db": os.getenv("DB_NAME"),
        "user": os.getenv("DB_USER"),
        "port": os.getenv("DB_PORT")
    }