import React from 'react';
import './App.css';
import UserForm from './components/UserForm';

function App() {
  return (
    <div className="App">
      <header className="App-header">
        <h1>React Form Exercise</h1>
        <p>Complete the form implementation to practice React form handling</p>
      </header>
      <main>
        <UserForm />
      </main>
      <footer>
        <p>Tasks to complete:</p>
        <ul>
          <li>Implement form validation for all fields</li>
          <li>Display error messages for invalid inputs</li>
          <li>Add loading state during form submission</li>
          <li>Show success message after successful submission</li>
          <li>Add proper styling for form states</li>
        </ul>
      </footer>
    </div>
  );
}

export default App;
