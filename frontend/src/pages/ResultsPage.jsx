import { useState, useEffect } from 'react'
import { getUploads, getVariants, getRiskSummary, getDrugs, getResearch, getConditions} from '../api/client'

function getRiskColor(score) {
    if (!score) return 'var(--text-dim)'
    if (score >= 0.7) return 'var(--risk-high)'
    if (score >= 0.4) return 'var(--risk-moderate)'
    if (score >= 0.1) return 'var(--risk-low)'
    return 'var(--risk-minimal)'
}

function getRiskLabel(score) {
    if (!score) return '—'
    if (score >= 0.7) return 'HIGH'
    if (score >= 0.4) return 'MODERATE'
    if (score >= 0.1) return 'LOW'
    return 'MINIMAL'
}



export default function ResultsPage() {
    const [uploads, setUploads] = useState([])
    const [selectedId, setSelectedId] = useState(null)
    const [variants, setVariant] = useState([])
    const [summary, setSummary] = useState(null)
    const [geneSearch, setGeneSearch] = useState('')
    const [flagFilter, setFlagFilter] = useState('all')
    const [minScore, setMinScore] = useState(0)
    const [flaggedOnly, setFlaggedOnly] = useState(false)
    const [drugs, setDrugs] = useState([])
    const [drugSearch, setDrugSearch] = useState('')
    const [research, setResearch] = useState(null)
    const [researchLoading, setResearchLoading] = useState(false)
    const [researchGene, setResearchGene] = useState(null)
    const [showModal, setShowModal] = useState(false)
    const [modalVariant, setModalVariant] = useState(null)
    const [conditions, setConditions] = useState([])

    useEffect(()=>{
        getUploads().then(data => setUploads(data))
    }, [])

    useEffect(() => {
    if (!selectedId) return
    getVariants(selectedId).then(data => setVariant(data))
    getRiskSummary(selectedId).then(data => setSummary(data[0]))
    getDrugs(selectedId).then(data => setDrugs(data))
    }, [selectedId])

    async function handleResearch(gene, condition) {
        setResearchGene(gene)
        setResearch(null)
        setResearchLoading(true)
        try {
            const data = await getResearch(gene, condition || 'disease')
            setResearch(data)
        } catch(e) {
            setResearch({error: "Agent unavailable"})
        } finally {
            setResearchLoading(false)
        }
    }

    async function openResearchModal(variant) {
        setModalVariant(variant)
        setShowModal(true)
        setConditions([])
        const data = await getConditions(variant.variant_id)
        setConditions(data)
    }

    const filtered = variants.filter(v => {
        if (geneSearch && !v.gene_name?.toLowerCase().includes(geneSearch.toLowerCase())) return false
        if (flagFilter !== 'all' && v.flag !== flagFilter) return false
        if (minScore > 0 && (!v.risk_score || v.risk_score < minScore)) return false
        if (flaggedOnly && !v.flag) return false
        return true
    })

    const filteredDrugs = drugs.filter(d => {
        if (!drugSearch) return true
        const parts = d.notes ? d.notes.split(' | ') : []
        const drug = parts[1]?.replace('Drug: ', '') || ''
        return drug.toLowerCase().includes(drugSearch.toLowerCase()) ||
           d.gene_name?.toLowerCase().includes(drugSearch.toLowerCase())
    })

    

    return (
     <div>
            <h1 className="page-title">Browse Results</h1>
            <p className="page-subtitle">Variants · Annotations · Risk Scores</p>

            <select onChange={e => setSelectedId(e.target.value)}>
                <option value="">Select an upload...</option>
                {uploads.map(u => (
                    <option key={u.upload_id} value={u.upload_id}>
                        {u.username} — {u.filename} ({u.total_variants} variants)
                    </option>
                ))}
            </select>

            {summary && (
                <>
                    <div className="cards-row">
                        <div className="score-card">
                            <div className="score-card-label">Overall Risk</div>
                            <div className="score-card-value" style={{color: getRiskColor(summary.overall_score)}}>
                                {summary.overall_score}
                            </div>
                            <div className="score-card-label-sub">{getRiskLabel(summary.overall_score)}</div>
                        </div>
                        <div className="score-card">
                            <div className="score-card-label">Pathogenicity</div>
                            <div className="score-card-value" style={{color: getRiskColor(summary.pathogenicity_score)}}>
                                {summary.pathogenicity_score}
                            </div>
                            <div className="score-card-label-sub">{getRiskLabel(summary.pathogenicity_score)}</div>
                        </div>
                        <div className="score-card">
                            <div className="score-card-label">Pharmacogenomics</div>
                            <div className="score-card-value" style={{color: getRiskColor(summary.pharmacogenomics_score)}}>
                                {summary.pharmacogenomics_score}
                            </div>
                            <div className="score-card-label-sub">{getRiskLabel(summary.pharmacogenomics_score)}</div>
                        </div>
                    </div>

                    
                    <div className="filters-row">
                        <input type="text" placeholder="Search gene..." onChange={e => setGeneSearch(e.target.value)}/>
                        <select onChange={e => setFlagFilter(e.target.value)}>
                            <option value="all">All flags</option>
                            <option value="clinical">Clinical</option>
                            <option value="pharmacogenomics">Pharmacogenomics</option>
                        </select>
                        <input type="number" placeholder="Min score" min="0" max="1" step="0.1" onChange={e => setMinScore(parseFloat(e.target.value) || 0)}/>
                        <label style={{display:'flex', alignItems:'center', gap:'0.5rem', color:'var(--text-muted)', fontSize:'0.8rem'}}>
                            <input type="checkbox" onChange={e => setFlaggedOnly(e.target.checked)}/>
                            Flagged only
                        </label>
                    </div>
                    <div style={{fontSize:'0.75rem', color:'var(--text-dim)', marginBottom:'0.5rem'}}>
                        {filtered.length} of {variants.length} variants
                    </div>
                    
                    {variants.length > 0 && (
                                            <div className="summary-bar">
                                                <span>Total: {variants.length}</span>
                                                <span style={{color: 'var(--risk-low)'}}>Clinical: {variants.filter(v => v.flag === 'clinical' || v.flag === 'both').length}</span>
                                                <span style={{color: 'var(--risk-minimal)'}}>Pharmacogenomics: {variants.filter(v => v.flag === 'pharmacogenomics' || v.flag === 'both').length}</span>
                                                <span style={{color: 'var(--text-dim)'}}>Unflagged: {variants.filter(v => !v.flag).length}</span>
                                                <span style={{color: 'var(--risk-high)'}}>High Risk: {variants.filter(v => v.risk_score >= 0.7).length}</span>
                                            </div>
                                            )}

                    <div className="table-wrapper">
                        <table>
                            <thead>
                                <tr>
                                    <th>Gene</th>
                                    <th>Chromosome</th>
                                    <th>Position</th>
                                    <th>Zygosity</th>
                                    <th>Flag</th>
                                    <th>Risk Score</th>
                                    <th>Agent Research</th>
                                </tr>
                            </thead>
                            <tbody>
                                {filtered.map(v => (
                                    <tr key={v.variant_id}>
                                        <td>{v.gene_name || '—'}</td>
                                        <td>{v.chromosome}</td>
                                        <td>{v.position}</td>
                                        <td>{v.zygosity}</td>
                                        <td>{v.flag || '—'}</td>
                                        <td style={{color: getRiskColor(v.risk_score), fontWeight: 600}}>
                                            {v.risk_score || '—'}
                                        </td>
                                        <td>
                                            {v.gene_name && v.flag === 'clinical' && (
                                                <button className='btn-research'
                                                    onClick={() => openResearchModal(v)}
                                                >
                                                    Research
                                                </button>
                                            )}
                                        </td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                    {researchLoading && (
                        <div className="research-panel">
                            <div className="research-panel-title">🔬 Searching PubMed + bioRxiv for {researchGene}...</div>
                        </div>
                    )}
                    
                    {research && !researchLoading && (
                        <div className="research-panel">
                            <div className="research-panel-title">
                                📄 Research Summary — {researchGene}
                                <span style={{fontSize: '0.7rem', color: 'var(--text-dim)', marginLeft: '1rem'}}>
                                    {research.papers_read} papers read
                                </span>
                            </div>
                            <div className="research-summary">
                                {research.summary}
                            </div>
                            {research.top_papers && (
                                <div style={{marginTop: '1rem'}}>
                                    <div className="score-card-label" style={{marginBottom: '0.5rem'}}>Top Papers</div>
                                    {research.top_papers.map(p => (
                                        <div key={p.pmid} className="research-paper">
                                            <a 
                                                href={`https://pubmed.ncbi.nlm.nih.gov/${p.pmid}`}
                                                target="_blank"
                                                rel="noreferrer"
                                            >
                                                {p.title}
                                            </a>
                                            <span style={{color: 'var(--text-dim)', fontSize: '0.7rem', marginLeft: '0.5rem'}}>
                                                {p.year}
                                            </span>
                                        </div>
                                    ))}
                                </div>
                            )}
                        </div>
                    )}
                </>
            )}

        {drugs.length > 0 && (
                <>
                    <h2 className="section-title">Drug Interactions</h2>
                    <div className="filters-row" style={{marginBottom: '0.75rem'}}>
                        <input
                            type="text"
                            placeholder="Search drug or gene..."
                            onChange={e => setDrugSearch(e.target.value)}
                        />
                        <div style={{fontSize:'0.75rem', color:'var(--text-dim)'}}>
                            {filteredDrugs.length} of {drugs.length} interactions
                        </div>
                    </div>
                    <div className="table-wrapper">
                        <table>
                            <thead>
                                <tr>
                                    <th>Gene</th>
                                    <th>Drug</th>
                                    <th>Effect</th>
                                    <th>Evidence</th>
                                    <th>Score</th>
                                </tr>
                            </thead>
                            <tbody>
                                {filteredDrugs.map((d, i) => {
                                    const parts = d.notes ? d.notes.split(' | ') : []
                                    return (
                                        <tr key={i}>
                                            <td>{d.gene_name}</td>
                                            <td>{parts[1]?.replace('Drug: ', '') || '—'}</td>
                                            <td>{parts[0] || '—'}</td>
                                            <td>{parts[2]?.replace('Evidence: ', '') || '—'}</td>
                                            <td style={{color: getRiskColor(d.risk_score), fontWeight: 600}}>{d.risk_score}</td>
                                        </tr>
                                    )
                                })}
                            </tbody>
                        </table>
                    </div>
                </>
            )}
        {showModal && modalVariant && (
            <div className="modal-overlay" onClick={() => setShowModal(false)}>
                <div className="modal" onClick={e => e.stopPropagation()}>
                    <div className="modal-header">
                        <span>🔬 Research {modalVariant.gene_name}</span>
                        <button className="modal-close" onClick={() => setShowModal(false)}>✕</button>
                    </div>
                    <div className="modal-body">
                        <p style={{color: 'var(--text-dim)', fontSize: '0.8rem', marginBottom: '1rem'}}>
                            Select a condition to research with the literature agent:
                        </p>
                        {conditions.length === 0 && (
                            <div className="loading">Loading conditions...</div>
                        )}
                        {conditions.map((c, i) => (
                            <button
                                key={i}
                                className="condition-btn"
                                onClick={() => {
                                    setShowModal(false)
                                    handleResearch(modalVariant.gene_name, c.condition_name)
                                }}
                            >
                                <span style={{flex: 1}}>{c.condition_name}</span>
                                <span style={{color: getRiskColor(c.risk_score), fontSize: '0.7rem', whiteSpace: 'nowrap'}}>
                                    {c.annotation_type}
                                </span>
                            </button>
                        ))}
                    </div>
                </div>
            </div>
        )}
        </div>
    )
}