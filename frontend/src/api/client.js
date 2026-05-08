const BASE_URL = process.env.REACT_APP_API_URL || 'http://localhost:8000'

export async function getUploads() {
    const response = await fetch(`${BASE_URL}/uploads`)
    return response.json()
}

export async function getVariants(uploadId) {
    const response = await fetch(`${BASE_URL}/variants/${uploadId}`)
    return response.json()
}

export async function getRiskSummary(uploadId) {
    const response = await fetch(`${BASE_URL}/risksummary/${uploadId}`)
    return response.json()
}

export async function getDrugs(uploadId) {
    const response = await fetch(`${BASE_URL}/drugs/${uploadId}`)
    return response.json()
}

export async function getStats() {
    const response = await fetch(`${BASE_URL}/stats`)
    return response.json()
}

export async function uploadVcf(file, username, email) {
    const formData = new FormData()
    formData.append('file', file)
    formData.append('username', username)
    formData.append('email', email)

    const response = await fetch(`${BASE_URL}/upload`, {
        method: 'POST',
        body: formData
    })
    return response.json()
}

export async function getTopGenes() {
    const response = await fetch(`${BASE_URL}/stats/genes`)
    return response.json()
}

export async function getRiskDistribution() {
    const response = await fetch(`${BASE_URL}/stats/risk-distribution`)
    return response.json()
}

export async function getReferenceStats() {
    const response = await fetch(`${BASE_URL}/stats/reference`)
    return response.json()
}

export async function getResearch(gene, condition) {
    const response = await fetch(`${BASE_URL}/research`, {
        method: 'POST',
        headers: {'Content-Type': 'application/json'},
        body: JSON.stringify({gene, condition})
    })
    return response.json()
}

export async function getConditions(variantId) {
    const response = await fetch(`${BASE_URL}/conditions/${variantId}`)
    return response.json()
}