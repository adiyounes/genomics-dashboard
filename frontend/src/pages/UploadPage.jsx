import { useState } from 'react'
import { uploadVcf } from '../api/client'

export default function UploadPage() {
    const [username, setUsername] = useState('')
    const [email, setEmail] = useState('')
    const [file, setFile] = useState(null)
    const [loading, setLoading] = useState(false)
    const [result, setResult] = useState(null)
    const [error, setError] = useState(null)

    async function handleSubmit() {
        if (!file || !username || !email) {
            setError('Please fill all fields')
            return
        }
        setLoading(true)
        setError(null)

        try {
            const data = await uploadVcf(file, username, email)
            setResult(data)
        } catch (err) {
            setError('Upload failed')
        } finally {
            setLoading(false)
        }
    }

    

    return (
  <div className="upload-page">
    <div className="page-title">Upload & Analyse</div>
    <div className="page-subtitle">VCF → Annotation → Risk Scores</div>

    <div className="form-group">
      <label className="form-label">Username</label>
      <input type="text" placeholder="e.g. patient_001" onChange={e => setUsername(e.target.value)} />
    </div>

    <div className="form-group">
      <label className="form-label">Email</label>
      <input type="email" placeholder="e.g. patient@clinic.org" onChange={e => setEmail(e.target.value)} />
    </div>

    <div className="form-group">
      <label className="form-label">VCF File</label>
      <input type="file" onChange={e => setFile(e.target.files[0])} />
    </div>

    <button className="btn" onClick={handleSubmit} disabled={loading}>
      {loading ? 'Analysing...' : 'Analyse'}
    </button>

    {error && <div className="alert alert-error">{error}</div>}
    {result && <div className="alert alert-success">✅ {Object.values(result.inserted || {}).reduce((a,b) => a+b, 0)} variants inserted</div>}
  </div>
)
}