import {useState} from 'react'

export default function UploadPage() {
    const [username, setUsername] = useState('')
    const [email, setEmail] = useState('')
    const [file, setFile] = useState(null)

    function handleSubmit() {
        console.log(username, email, file)
    }

    return (
        <div>
            <h1>Upload VCF File</h1>
            <input type="text" onChange={e => setUsername(e.target.value)}/>
            <input type="text" onChange={e => setEmail(e.target.value)}/>
            <input type="file" onChange={e => setFile(e.target.files[0])}/>
            <button onClick={handleSubmit()}>Analyse</button>
        </div>
    )
}