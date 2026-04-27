import { useState, useEffect } from 'react'
import { getUploads } from '../api/client'

export default function ResultsPage() {
    const [uploads, setUploads] = useState([])
    const [selectedId, setSelectedId] = useState(null)

    useEffect(()=>{

    })
    return (
    <div>
        <h1>Browse Result</h1>
    </div>
    )
}