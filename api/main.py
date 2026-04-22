from fastapi import FastAPI, UploadFile, File, Form
from fastapi.middleware.cors import CORSMiddleware
import sys
from pathlib import Path
import httpx

sys.path.append(str(Path(__file__).parent.parent))

from database.connect import execute_query, execute_insert
from modules.ingestion.vcf_parser import ingest_vcf
from modules.annotation.annotator import annotate_upload

app = FastAPI(title="GenomeDX API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"],
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

AGENT_API = "http://54.195.19.197:8000/analyse"

@app.post("/research")
async def get_research(gene: str, condition: str):
    try:
        async with httpx.AsyncClient(timeout=60.0) as client:
            response = await client.post(
                AGENT_API,
                json={"gene": gene, "disease":  condition}
            )
            if response.status_code == 200:
                return response.json()
            return {"error": "Agent unailable"}
    except Exception as e:
        return {"error": str(e)}