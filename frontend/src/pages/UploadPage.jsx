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
        <div>
            <h1>Upload VCF File</h1>
            <input type="text" placeholder="Username" onChange={e => setUsername(e.target.value)}/>
            <input type="text" placeholder='Email' onChange={e => setEmail(e.target.value)}/>
            <input type="file" onChange={e => setFile(e.target.files[0])}/>
            <button onClick={handleSubmit} disabled={loading}>
                {loading ? 'Analyse...' : 'Analyse'}</button>
            {error && <p style={{color: 'red'}}>{error}</p>}
            {result && <p style={{color: 'green'}}>✅ {result.inserted ? Object.values(result.inserted).reduce((a,b) => a+b, 0) : 0} variants inserted</p>}
        </div>
    )
}