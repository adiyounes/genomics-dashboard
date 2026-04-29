import { useState, useEffect } from 'react'
import { getUploads, getVariants, getRiskSummary } from '../api/client'

export default function ResultsPage() {
    const [uploads, setUploads] = useState([])
    const [selectedId, setSelectedId] = useState(null)
    const [variants, setVariant] = useState([])
    const [summary, setSummary] = useState(null)

    useEffect(()=>{
        getUploads().then(data => setUploads(data))
    }, [])

    useEffect(()=> {
        if(!selectedId) return
        getVariants(selectedId).then(data => setVariant(data))
        getRiskSummary(selectedId).then(data=>setSummary([0]))
    },[selectedId])

    return (
    <div>
        <h1>Browse Result</h1>
        <select onChange={e => setSelectedId(e.target.value)}>
            {uploads.map(u => (
                <option key={u.upload_id} value={u.upload_id}>
                    {u.username} — {u.filename} ({u.total_variants} variants)
                </option>
            ))}
        </select>
        {summary && (
            <div>
                <h2>Risk Summary</h2>
                <p>Overall: {summary.overall_score}</p>
                <p>Pathogenicity: {summary.pathogenicity_score}</p>
                <p>Pharmacogenomics: {summary.pharmacogenomics_score}</p>
            </div>
        )}
        {variants.length > 0 && (
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
                        {variants.map(v => (
                            <tr key={v.variant_id}>
                                <td>{v.gene_name}</td>
                                <td>{v.chromosome}</td>
                                <td>{v.position}</td>
                                <td>{v.zygosity}</td>
                                <td>{v.flag}</td>
                                <td>{v.risk_score}</td>
                            </tr>
                        ))}
                    </tbody>
                </table>
            )}
    </div>
    )
}