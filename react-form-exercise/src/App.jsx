import './App.css';
import React, {useState, useRef} from 'react';

function App() {
    const [timer, setTimer] = useState(0);
    const intervalRef = useRef(null);

    const handleStart = () => {
        if (intervalRef.current !== null) return;
        
        intervalRef.current = setInterval(() => {
            setTimer(prev => prev + 1);
        }, 1000);
    }
    
    const handleStop = () => {
        if (intervalRef.current === null) return;
        
        clearInterval(intervalRef.current);
        intervalRef.current = null;
    }
    
    const handleReset = () => {
        handleStop();
        setTimer(0);
    }
    
  return (
    <div className="App">
     <div className="timer">{timer}</div>
     <div className="controls">
       <button onClick={handleStart}>Start</button>
       <button onClick={handleStop}>Stop</button>
       <button onClick={handleReset}>Reset</button>
     </div>
    </div>
  );
}

export default App;