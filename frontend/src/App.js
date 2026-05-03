import {BrowserRouter, Routes, Route, NavLink} from 'react-router-dom'
import UploadPage from './pages/UploadPage'
import ResultsPage from './pages/ResultsPage'
import StatsPage from './pages/StatsPage'
import './App.css'

function App() {
  return (
    <BrowserRouter>
      <div className="app-layout">
        <nav>
          <NavLink to="/" className="nav-brand" e>
            🧬 GenomeDx <span>Genomic Analysis Platform</span>
          </NavLink>
          <NavLink to="/" end className={({isActive}) => isActive ? 'nav-link active' : 'nav-link'}>Upload</NavLink>
          <NavLink to="/results" className={({isActive}) => isActive ? 'nav-link active' : 'nav-link'}>Results</NavLink>
          <NavLink to="/stats" className={({isActive}) => isActive ? 'nav-link active' : 'nav-link'}>Stats</NavLink>
        </nav>
        <main>
          <Routes>
            <Route path="/" element={<UploadPage />} />
            <Route path="/results" element={<ResultsPage />} />
            <Route path="/stats" element={<StatsPage />} />
          </Routes>
        </main>
      </div>
    </BrowserRouter>
  )
}

export default App;
