import { useState, useEffect } from 'react'
import { getUploads, getVariants, getRiskSummary } from '../api/client'

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

    useEffect(()=>{
        getUploads().then(data => setUploads(data))
    }, [])

    useEffect(()=> {
        if(!selectedId) return
        getVariants(selectedId).then(data => setVariant(data))
        getRiskSummary(selectedId).then(data=>{
            setSummary(data[0])})
    },[selectedId])

    const filtered = variants.filter(v => {
        if (geneSearch && !v.gene_name?.toLowerCase().includes(geneSearch.toLowerCase())) return false
        if (flagFilter !== 'all' && v.flag !== flagFilter) return false
        if (minScore > 0 && (!v.risk_score || v.risk_score < minScore)) return false
        if (flaggedOnly && !v.flag) return false
        return true
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
                        <input type="text" placeholder="Search gene..." onChange={e => setGeneSearch(e.target.value)} style={{flex:1}}/>
                        <select onChange={e => setFlagFilter(e.target.value)} style={{width:'auto'}}>
                            <option value="all">All flags</option>
                            <option value="clinical">Clinical</option>
                            <option value="pharmacogenomics">Pharmacogenomics</option>
                            <option value="both">Both</option>
                        </select>
                        <input type="number" placeholder="Min score" min="0" max="1" step="0.1" onChange={e => setMinScore(parseFloat(e.target.value) || 0)} style={{width:'100px'}}/>
                        <label style={{display:'flex', alignItems:'center', gap:'0.5rem', color:'var(--text-muted)', fontSize:'0.8rem'}}>
                            <input type="checkbox" onChange={e => setFlaggedOnly(e.target.checked)}/>
                            Flagged only
                        </label>
                    </div>
                    <div style={{fontSize:'0.75rem', color:'var(--text-dim)', marginBottom:'0.5rem'}}>
                        {filtered.length} of {variants.length} variants
                    </div>

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
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </>
            )}
        </div>
    )
}