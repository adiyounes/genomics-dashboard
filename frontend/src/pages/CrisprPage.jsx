import { useState, useEffect } from "react"
import { analyzeCrispr, getSimulations } from "../api/client"


export default function CrisprPage() {
    const [geneName, setGeneName] = useState('')
    const [variant, setVariant] = useState('')
    const [simulation, setSimulation] = useState(null)
    const [simulations, setSimulations] = useState([])
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)

    useEffect(() => {
        async function fetchSimulations() {
            try {
                const data = await getSimulations()
                setSimulations(data)
            } catch (err) {
                setError(err)
            }
        }
        fetchSimulations()
    }, [])

    const handleSubmit = async (e) => {
        e.preventDefault()
        setLoading(true)
        setError(null)
        try {
            const result = await analyzeCrispr(geneName, variant)
            setSimulation(result)
        } catch (err) {
            setError(err)
        } finally {
            setLoading(false)
        }
    }

    return (
        <div>
            <h1>CRISPR Analysis</h1>
            <form onSubmit={handleSubmit}>
                <input
                    type="text"
                    placeholder="Gene Name"
                    value={geneName}
                    onChange={(e) => setGeneName(e.target.value)}
                />
                <input
                    type="text"
                    placeholder="Variant"
                    value={variant}
                    onChange={(e) => setVariant(e.target.value)}
                />
                <button type="submit" disabled={loading}>
                    {loading ? 'Analyzing...' : 'Analyze'}
                </button>
            </form>
            {error && <p style={{color: 'red'}}>Error: {error.message}</p>}
            {simulation && (
                <div>
                    <h2>Result</h2>
                    <p>Verdict: <strong>{simulation.verdict}</strong></p>
                    <p>Risk Score: {simulation.risk_score}</p>
                    <p>gRNA: {simulation.grna_sequence}</p>
                    <p>Reasoning: {simulation.reasoning}</p>
                </div>
            )}
            <h2>Previous Simulations</h2>
            <ul>
                {simulations.map((sim) => (
                    <li key={sim.simulation_id}>{sim.gene} - {sim.variant} - {new Date(sim.analyzed_at).toLocaleString()}</li>
                ))}
            </ul>
        </div>
    )
}