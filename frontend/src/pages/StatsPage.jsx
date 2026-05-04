import { useState, useEffect } from 'react'
import { getStats, getTopGenes, getRiskDistribution, getReferenceStats } from '../api/client'

export default function StatsPage() {
    const [stats, setStats] = useState(null)
    const [genes, setGenes] = useState([])
    const [riskDist, setRiskDist] = useState(null)
    const [reference, setReference] = useState(null)

    useEffect(() => {
        getStats().then(data => setStats(data))
        getTopGenes().then(data => setGenes(data))
        getRiskDistribution().then(data => setRiskDist(data[0]))
        getReferenceStats().then(data => setReference(data[0]))
    }, [])

    if (!stats) return <div className="loading">Loading...</div>

    return (
        <div>
            <h1 className="page-title">Database Stats</h1>
            <p className="page-subtitle">Knowledge Base · Pipeline Health</p>

            <div className="cards-row" style={{marginTop: '2rem'}}>
                <div className="score-card">
                    <div className="score-card-label">ClinVar Variants</div>
                    <div className="score-card-value" style={{color: 'var(--risk-low)'}}>
                        {stats.clinvar_count?.toLocaleString()}
                    </div>
                    <div className="score-card-label-sub">NIH disease variants</div>
                </div>
                <div className="score-card">
                    <div className="score-card-label">PharmGKB Annotations</div>
                    <div className="score-card-value" style={{color: 'var(--risk-minimal)'}}>
                        {stats.pharmgkb_count?.toLocaleString()}
                    </div>
                    <div className="score-card-label-sub">Gene-drug relationships</div>
                </div>
                <div className="score-card">
                    <div className="score-card-label">Variants Analysed</div>
                    <div className="score-card-value" style={{color: 'var(--risk-moderate)'}}>
                        {stats.variant_count?.toLocaleString()}
                    </div>
                    <div className="score-card-label-sub">Across all uploads</div>
                </div>
                <div className="score-card">
                    <div className="score-card-label">Total Uploads</div>
                    <div className="score-card-value" style={{color: 'var(--text-primary)'}}>
                        {stats.upload_count}
                    </div>
                    <div className="score-card-label-sub">Completed analyses</div>
                </div>
            </div>

            <div style={{display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1.5rem', marginTop: '2rem'}}>

                {reference && (
                    <div className="score-card">
                        <div className="score-card-label" style={{marginBottom: '1rem'}}>Genomic Reference</div>
                        <div style={{display: 'flex', flexDirection: 'column', gap: '0.75rem'}}>
                            <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
                                <span style={{color: 'var(--text-muted)', fontSize: '0.8rem', fontFamily: 'var(--font-mono)'}}>Gene coordinates</span>
                                <span style={{color: 'var(--text-primary)', fontFamily: 'var(--font-mono)', fontSize: '0.9rem'}}>{reference.genes?.toLocaleString()}</span>
                            </div>
                            <div style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
                                <span style={{color: 'var(--text-muted)', fontSize: '0.8rem', fontFamily: 'var(--font-mono)'}}>Exon coordinates</span>
                                <span style={{color: 'var(--text-primary)', fontFamily: 'var(--font-mono)', fontSize: '0.9rem'}}>{reference.exons?.toLocaleString()}</span>
                            </div>
                        </div>
                    </div>
                )}

                {riskDist && (
                    <div className="score-card">
                        <div className="score-card-label" style={{marginBottom: '1rem'}}>Risk Distribution</div>
                        <div style={{display: 'flex', flexDirection: 'column', gap: '0.75rem'}}>
                            {[
                                {label: 'High', value: riskDist.high, color: 'var(--risk-high)'},
                                {label: 'Moderate', value: riskDist.moderate, color: 'var(--risk-moderate)'},
                                {label: 'Low', value: riskDist.low, color: 'var(--risk-low)'},
                                {label: 'Minimal', value: riskDist.minimal, color: 'var(--risk-minimal)'},
                                {label: 'Unannotated', value: riskDist.unannotated, color: 'var(--text-dim)'},
                            ].map(item => (
                                <div key={item.label} style={{display: 'flex', justifyContent: 'space-between', alignItems: 'center'}}>
                                    <span style={{color: item.color, fontSize: '0.8rem', fontFamily: 'var(--font-mono)'}}>{item.label}</span>
                                    <span style={{color: 'var(--text-primary)', fontFamily: 'var(--font-mono)', fontSize: '0.9rem'}}>{item.value?.toLocaleString()}</span>
                                </div>
                            ))}
                        </div>
                    </div>
                )}
            </div>

            {genes.length > 0 && (
                <div style={{marginTop: '1.5rem'}}>
                    <h2 className="section-title">Top Flagged Genes</h2>
                    <div className="table-wrapper">
                        <table>
                            <thead>
                                <tr>
                                    <th>Gene</th>
                                    <th>Variant Count</th>
                                </tr>
                            </thead>
                            <tbody>
                                {genes.map(g => (
                                    <tr key={g.gene_name}>
                                        <td style={{fontFamily: 'var(--font-mono)'}}>{g.gene_name}</td>
                                        <td>{g.count?.toLocaleString()}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    </div>
                </div>
            )}

            <div className="disclaimer" style={{marginTop: '3rem'}}>
                ⚠ This platform is intended for research purposes only.
                Results should not be used as a substitute for professional
                medical advice, diagnosis, or treatment.
            </div>
        </div>
    )
}