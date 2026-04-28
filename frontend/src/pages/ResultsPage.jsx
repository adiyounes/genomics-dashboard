import { useState, useEffect } from 'react'
import { getUploads } from '../api/client'

export default function ResultsPage() {
    const [uploads, setUploads] = useState([])
    const [selectedId, setSelectedId] = useState(null)

    useEffect(()=>{
        getUploads().then(data => setUploads(data))
    }, [])
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
    </div>
    )
}