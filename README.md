# GenomeDx


[![CI](https://github.com/adiyounes/genomics-dashboard/actions/workflows/ci.yml/badge.svg)](https://github.com/adiyounes/genomics-dashboard/actions/workflows/ci.yml)


A genomics analysis platform that helps identify disease risks and drug interactions associated with a patient's genetic data. Researchers and clinicians upload a VCF file and the system returns a risk assessment across pathogenicity and pharmacogenomics domains, backed by 2.6M ClinVar variants and 5,120 PharmGKB gene-drug relationships.
 
---
 
## What it does
 
A researcher uploads a VCF file containing a patient's genetic mutations. The pipeline parses the file, resolves gene names from GRCh38 coordinates, and searches for associations in ClinVar and PharmGKB. It returns all flagged variants with risk scores, drug interactions, and a summary. The researcher can then explore the data or run the integrated literature agent to get AI-generated research summaries from PubMed and bioRxiv for any specific gene-disease pair.
 
---
 
## Features
 
- Parses real GRCh38 VCF files, single or multi-sample
- Resolves gene names from chromosome coordinates using 40,502 genes and 1.5M exons from Ensembl
- Matches variants against ClinVar with three confidence levels: exact, position, and gene
- Identifies drug interactions from PharmGKB ordered by evidence strength
- Calculates weighted risk scores per variant and per upload
- Integrates an autonomous literature mining agent that searches PubMed and bioRxiv and returns structured research summaries with citations
- React frontend with filters, risk score cards, and drug interaction tables
- FastAPI backend with 10+ endpoints
- 11 passing unit tests with CI via GitHub Actions
---
 
## Stack
 
| Layer | Technology |
|-------|-----------|
| Backend | Python, FastAPI |
| Database | PostgreSQL |
| Frontend | React |
| AI agent | Anthropic Claude API |
| Knowledge bases | NIH ClinVar, Stanford PharmGKB, Ensembl GTF |
| HTTP | httpx |
| Testing | pytest |
 
---
 
## Data sources
 
The pipeline uses three public datasets that need to be downloaded before running:
 
**ClinVar**  2.6M known disease variants from NIH
```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz -P data/raw/
```
 
**Ensembl GTF**  gene and exon coordinates for GRCh38
```bash
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz -P data/raw/
```
 
**PharmGKB**  gene-drug relationships from Stanford
```bash
wget https://api.pharmgkb.org/v1/download/file/data/relationships.zip -P data/raw/
```
 
---
 
## How to run
 
**Requirements:** Python 3.10, PostgreSQL 14, Node.js 20
 
```bash
# 1. Clone the repo
git clone git@github.com:adiyounes/genomics-dashboard.git
cd genomics-dashboard
 
# 2. Create virtual environment
python3 -m venv venv
source venv/bin/activate
 
# 3. Install dependencies
pip install -r requirements.txt
 
# 4. Set up database
createdb genomics_db
psql genomics_db < database/schema.sql
 
# 5. Configure credentials
cp database/credentials.example.py database/credentials.py
# Edit credentials.py with your database details
 
# 6. Load knowledge bases (takes 10-20 minutes)
python3 database/load_clinvar.py
python3 database/load_pharmgkb.py
python3 database/load_gtf.py
 
# 7. Start the API
uvicorn api.main:app --reload --port 8000
 
# 8. Start the frontend
cd frontend
npm install
npm start
```
 
Open `http://localhost:3000` in your browser.
 
---
 
## Testing
 
```bash
pytest tests/ -v
```
 
11 tests covering VCF parsing, annotation scoring, and database connection.
 
---
 
## Limitations
 
- Only accepts GRCh38 VCF files, older builds like GRCh37 are rejected
- Flags variants in 35 clinically established genes rare disease genes are not covered
- PharmGKB matching is gene-level only, not variant-level
- Full genome ingestion is slow for large files recommended to extract regions of interest first
---
 
## What's next
 
- **CRISPR Safety Prediction Agent** predicts off-target cut sites for gene therapy planning
- **Microbiome Drug Metabolism Agent** cross-references gut microbiome composition with pharmacogenomics to refine drug metabolism predictions
---
 
## License
 
MIT
