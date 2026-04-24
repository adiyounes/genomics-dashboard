import {BrowserRouter, Routes, Route, NavLink} from 'react-router-dom'
import UploadPage from './pages/UploadPage'
import ResultsPage from './pages/ResultsPage'
import StatsPage from './pages/StatsPage'
import './App.css'

function App() {
  return (
    <BrowserRouter>
      <nav>
        <NavLink to="/">Upload</NavLink>
        <NavLink to="/results">Results</NavLink>
        <NavLink to="/stats">Stats</NavLink>
      </nav>
      <main>
        <Routes>
          <Route path="/" element={<UploadPage />} />
          <Route path="/results" element={<ResultsPage />} />
          <Route path="/stats" element={<StatsPage />} />
        </Routes>
      </main>
    </BrowserRouter>
  );
}

export default App;
