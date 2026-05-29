import { useState, useEffect } from "react"
import { analyzeCrispr, getSimulations, getSimulation } from "../api/client"

export default function CrisprPage() {
    const [geneName, setGeneName] = useState('')
    const [variant, setVariant] = useState('')
    const [simulation, setSimulation] = useState(null)
    const [simulations, setSimulations] = useState([])
    const [loading, setLoading] = useState(false)
    const [error, setError] = useState(null)
    const [selectedSimulation, setSelectedSimulation] = useState(null)

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

    const handleSimulationClick = async (simulationId) => {
        try {
            const data = await getSimulation(simulationId)
            setSelectedSimulation(data)
        } catch (err) {
            setError(err)
            }
    }

    return (
    <div className="crispr-page">
        <h1>CRISPR Safety Analysis</h1>

        <form className="crispr-form" onSubmit={handleSubmit}>
            <input
                type="text"
                placeholder="Gene (e.g. BRCA1)"
                value={geneName}
                onChange={(e) => setGeneName(e.target.value)}
            />
            <input
                type="text"
                placeholder="Variant (e.g. c.5266dupC)"
                value={variant}
                onChange={(e) => setVariant(e.target.value)}
            />
            <button type="submit" disabled={loading}>
                {loading ? 'Analyzing...' : 'Analyze'}
            </button>
        </form>

        {error && <div className="crispr-error">Error: {error.message}</div>}

        {simulation && (
            <div className="crispr-result">
                <h2>Result</h2>
                <div className="crispr-result-row">
                    <span className="crispr-result-label">Verdict</span>
                    <span className={`verdict-badge verdict-${simulation.verdict}`}>
                        {simulation.verdict}
                    </span>
                </div>
                <div className="crispr-result-row">
                    <span className="crispr-result-label">Risk Score</span>
                    <span className="crispr-result-value">{simulation.risk_score}</span>
                </div>
                <div className="crispr-result-row">
                    <span className="crispr-result-label">gRNA</span>
                    <span className="crispr-result-value">{simulation.grna_sequence}</span>
                </div>
                <p className="crispr-reasoning">{simulation.reasoning}</p>
            </div>
        )}

        <h2>Previous Simulations</h2>
        <div className="simulations-list">
            {simulations.map((sim) => (
                <div className="simulation-row" 
                    key={sim.simulation_id}  
                    onClick={() => handleSimulationClick(sim.simulation_id)}>
                        <span className="sim-gene">{sim.gene}</span>
                        <span className="sim-variant">{sim.variant}</span>
                        <span className={`verdict-badge verdict-${sim.verdict}`}>{sim.verdict}</span>
                        <span className="sim-date">{new Date(sim.analyzed_at).toLocaleString()}</span>
                </div>
            ))}
        </div>
        {selectedSimulation && (
        <div className="crispr-result" style={{marginTop: 'var(--space-5)'}}>
            <h2>Simulation #{selectedSimulation.simulation_id} — {selectedSimulation.gene}</h2>
            <div className="crispr-result-row">
                <span className="crispr-result-label">Verdict</span>
                <span className={`verdict-badge verdict-${selectedSimulation.verdict}`}>
                    {selectedSimulation.verdict}
                </span>
            </div>
            <div className="crispr-result-row">
                <span className="crispr-result-label">gRNA</span>
                <span className="crispr-result-value">{selectedSimulation.grna_sequence}</span>
            </div>
            <p className="crispr-reasoning">{selectedSimulation.grna_reasoning}</p>
            <p className="crispr-reasoning">{selectedSimulation.reasoning}</p>
        </div>)
        }
    </div>
)
}