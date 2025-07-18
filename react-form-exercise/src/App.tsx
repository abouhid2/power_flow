import React, { useState } from 'react';
import './App.css';
import UserForm from './components/UserForm';
import ApiDataPage from './components/ApiDataPage';

function App() {
  const [activePage, setActivePage] = useState<'form' | 'api'>('form');

  return (
    <div className="App">
      <header className="App-header">
        <h1>React Exercise</h1>
        <p>Practice React with forms and API data fetching</p>
        <nav className="navigation">
          <button 
            className={activePage === 'form' ? 'active' : ''} 
            onClick={() => setActivePage('form')}
          >
            Form Exercise
          </button>
          <button 
            className={activePage === 'api' ? 'active' : ''} 
            onClick={() => setActivePage('api')}
          >
            API Data
          </button>
        </nav>
      </header>
      <main>
        {activePage === 'form' ? <UserForm /> : <ApiDataPage />}
      </main>
      <footer>
        <p>Tasks to complete:</p>
        <ul>
          {activePage === 'form' ? (
            <>
              <li>Implement form validation for all fields</li>
              <li>Display error messages for invalid inputs</li>
              <li>Add loading state during form submission</li>
              <li>Show success message after successful submission</li>
              <li>Add proper styling for form states</li>
            </>
          ) : (
            <>
              <li>Fetch data from an API</li>
              <li>Display loading state while fetching</li>
              <li>Handle and display errors</li>
              <li>Render the fetched data</li>
              <li>Add interactivity to the displayed data</li>
            </>
          )}
        </ul>
      </footer>
    </div>
  );
}

export default App;
