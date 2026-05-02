import { useState, useEffect } from 'react'
import { getStats } from '../api/client'

export default function StatsPage() {
    const [stats, setStats] = useState(null)

    useEffect(() => {
        getStats().then(data => setStats(data))
    }, [])

    if (!stats) return <div>Loading...</div>

    return (
        <div>
            <h1>Database Stats</h1>
            <p>ClinVar variants: {stats.clinvar_count?.toLocaleString()}</p>
            <p>PharmGKB annotations: {stats.pharmgkb_count?.toLocaleString()}</p>
            <p>Total variants: {stats.variant_count?.toLocaleString()}</p>
            <p>Total uploads: {stats.upload_count}</p>
        </div>
    )
}