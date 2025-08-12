import React, { useState, useEffect, useRef, useReducer, useCallback, useMemo, createContext, useContext } from 'react';
import './UserForm.css';

interface UserData {
  name: string;
  email: string;
  age: string;
  occupation: string;
}

// Example 7: useContext - Create a form context
const FormContext = createContext<{
  isDarkMode: boolean;
  toggleDarkMode: () => void;
}>({ isDarkMode: false, toggleDarkMode: () => {} });

// Example 8: Custom hook for form validation
const useFormValidation = (userData: UserData) => {
  const [errors, setErrors] = useState<Record<string, string>>({});
  
  // Validate a specific field
  const validateField = useCallback((name: string, value: string) => {
    let newErrors = { ...errors };
    
    if (name === 'name') {
      if (!value.trim()) {
        newErrors[name] = 'Name is required';
      } else if (value.length < 8) {
        newErrors[name] = 'Name must be at least 8 characters long';
      } else {
        delete newErrors[name];
      }
    }

    if (name === 'email') {
      const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
      if (!value.trim()) {
        newErrors[name] = 'Email is required';
      } else if (!emailRegex.test(value)) {
        newErrors[name] = 'Please enter a valid email address';
      } else {
        delete newErrors[name];
      }
    }

    if (name === 'age') {
      const ageValue = parseInt(value);
      if (!value.trim()) {
        newErrors[name] = 'Age is required';
      } else if (isNaN(ageValue) || ageValue < 18 || ageValue > 100) {
        newErrors[name] = 'Age must be between 18 and 100';
      } else {
        delete newErrors[name];
      }
    }

    if (name === 'occupation') {
      if (!value.trim()) {
        newErrors[name] = 'Occupation is required';
      } else if (value.length < 3) {
        newErrors[name] = 'Occupation must be at least 3 characters long';
      } else {
        delete newErrors[name];
      }
    }
    
    setErrors(newErrors);
    return newErrors;
  }, [errors]);
  
  // Validate the entire form
  const validateForm = useCallback(() => {
    let newErrors = { ...errors };
    
    // Validate all fields
    const fields = ['name', 'email', 'age', 'occupation'] as const;
    fields.forEach(field => {
      newErrors = { ...newErrors, ...validateField(field, userData[field]) };
    });
    
    setErrors(newErrors);
    return Object.keys(newErrors).length === 0;
  }, [userData, errors, validateField]);
  
  return { errors, validateField, validateForm };
};

// Example 4: useReducer - Form state reducer
type FormAction = 
  | { type: 'UPDATE_FIELD'; field: string; value: string }
  | { type: 'SET_SUBMITTING'; value: boolean }
  | { type: 'SET_SUBMITTED'; value: boolean }
  | { type: 'RESET_FORM' };

const formReducer = (state: { 
  userData: UserData; 
  isSubmitting: boolean; 
  isSubmitted: boolean 
}, action: FormAction) => {
  switch (action.type) {
    case 'UPDATE_FIELD':
      return {
        ...state,
        userData: {
          ...state.userData,
          [action.field]: action.value
        }
      };
    case 'SET_SUBMITTING':
      return {
        ...state,
        isSubmitting: action.value
      };
    case 'SET_SUBMITTED':
      return {
        ...state,
        isSubmitted: action.value
      };
    case 'RESET_FORM':
      return {
        ...state,
        userData: {
          name: '',
          email: '',
          age: '',
          occupation: ''
        },
        isSubmitted: false
      };
    default:
      return state;
  }
};

const UserForm: React.FC = () => {
  // Example 7: useContext - Using the context
  const [isDarkMode, setIsDarkMode] = useState(false);
  const toggleDarkMode = () => setIsDarkMode(!isDarkMode);
  const formContextValue = { isDarkMode, toggleDarkMode };
  
  // Example 4: useReducer - Replace multiple useState calls
  const [formState, dispatch] = useReducer(formReducer, {
    userData: {
      name: '',
      email: '',
      age: '',
      occupation: ''
    },
    isSubmitting: false,
    isSubmitted: false
  });
  
  const { userData, isSubmitting, isSubmitted } = formState;
  
  // Example 8: Custom hook - Use our form validation hook
  const { errors, validateField, validateForm } = useFormValidation(userData);
  
  // Example 3: useRef - Create refs for form elements
  const nameInputRef = useRef<HTMLInputElement>(null);
  const emailInputRef = useRef<HTMLInputElement>(null);
  const ageInputRef = useRef<HTMLInputElement>(null);
  const occupationInputRef = useRef<HTMLInputElement>(null);
  
  // Map field names to their refs for error focusing
  const fieldRefs = {
    name: nameInputRef,
    email: emailInputRef,
    age: ageInputRef,
    occupation: occupationInputRef
  };
  
  // Focus the first field with an error
  const focusFirstError = useCallback((currentErrors: Record<string, string>) => {
    const firstErrorField = Object.keys(currentErrors)[0];
    if (firstErrorField && fieldRefs[firstErrorField as keyof typeof fieldRefs]?.current) {
      fieldRefs[firstErrorField as keyof typeof fieldRefs].current?.focus();
    }
  }, [fieldRefs]);
  
  // Example 5: useCallback - Memoize event handlers
  const handleChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    
    // Update form state
    dispatch({ type: 'UPDATE_FIELD', field: name, value });
    
    // Validate the field
    validateField(name, value);
  }, [validateField]);
  
  // Handle form submission with useCallback
  const handleSubmit = useCallback((e: React.FormEvent) => {
    e.preventDefault();
    
    // Validate the form
    const isValid = validateForm();
    
    if (isValid) {
      dispatch({ type: 'SET_SUBMITTING', value: true });
      
      // Simulate API call with setTimeout
      setTimeout(() => {
        console.log('Form submitted with:', userData);
        dispatch({ type: 'SET_SUBMITTING', value: false });
        dispatch({ type: 'SET_SUBMITTED', value: true });
        
        // Reset form after successful submission
        dispatch({ type: 'RESET_FORM' });
      }, 1500);
    } else {
      dispatch({ type: 'SET_SUBMITTED', value: false });
      // Focus the first field with an error
      focusFirstError(errors);
    }
  }, [userData, validateForm, errors, focusFirstError]);
  
  // Example 6: useMemo - Calculate derived state
  const isFormValid = useMemo(() => {
    return Object.keys(errors).length === 0 &&
      userData.name.trim() !== '' &&
      userData.email.trim() !== '' &&
      userData.age.trim() !== '' &&
      userData.occupation.trim() !== '';
  }, [userData, errors]);
  
  // Format user data for submission
  const formattedUserData = useMemo(() => {
    return {
      ...userData,
      age: userData.age ? parseInt(userData.age) : 0,
      fullName: userData.name.trim(),
      emailAddress: userData.email.toLowerCase().trim()
    };
  }, [userData]);
  
  return (
    <FormContext.Provider value={formContextValue}>
      <div className={`form-container ${isDarkMode ? 'dark-mode' : ''}`}>
        <div className="theme-toggle">
          <button onClick={toggleDarkMode}>
            {isDarkMode ? 'Switch to Light Mode' : 'Switch to Dark Mode'}
          </button>
        </div>
        
        <h2>User Registration Form</h2>
        
        {isSubmitted && (
          <div className="success-message">
            Form submitted successfully!
          </div>
        )}
        
        <form onSubmit={handleSubmit}>
          <div className="form-group">
            <label htmlFor="name">Name:</label>
            <input
              ref={nameInputRef}
              type="text"
              id="name"
              name="name"
              value={userData.name}
              onChange={handleChange}
            />
            {errors["name"] && <p className="error-message">{errors["name"]}</p>}
          </div>

          <div className="form-group">
            <label htmlFor="email">Email:</label>
            <input
              ref={emailInputRef}
              type="email"
              id="email"
              name="email"
              value={userData.email}
              onChange={handleChange}
            />
            {errors["email"] && (
              <p className="error-message">{errors["email"]}</p>
            )}
          </div>

          <div className="form-group">
            <label htmlFor="age">Age:</label>
            <input
              ref={ageInputRef}
              type="number"
              id="age"
              name="age"
              value={userData.age}
              onChange={handleChange}
            />
            {errors["age"] && <p className="error-message">{errors["age"]}</p>}
          </div>

          <div className="form-group">
            <label htmlFor="occupation">Occupation:</label>
            <input
              ref={occupationInputRef}
              type="text"
              id="occupation"
              name="occupation"
              value={userData.occupation}
              onChange={handleChange}
            />
            {errors["occupation"] && <p className="error-message">{errors["occupation"]}</p>}
          </div>

          {/* Example 6: useMemo - Display submission status based on derived state */}
          <div className="form-status">
            <p>Form Status: {isFormValid ? 'Valid ✅' : 'Invalid ❌'}</p>
            <p>Fields completed: {
              Object.values(userData).filter(val => val.trim() !== '').length
            } of 4</p>
          </div>

          <button type="submit" disabled={isSubmitting || !isFormValid}>
            {isSubmitting ? 'Submitting...' : 'Submit'}
          </button>
        </form>
        
        {/* Example 6: useMemo - Display formatted data */}
        {userData.name && (
          <div className="preview-data">
            <h3>Form Preview:</h3>
            <pre>{JSON.stringify(formattedUserData, null, 2)}</pre>
          </div>
        )}
      </div>
    </FormContext.Provider>
  );
};

// Example 7: useContext - Create a component that uses the context
const ThemeToggle: React.FC = () => {
  const { isDarkMode, toggleDarkMode } = useContext(FormContext);
  
  return (
    <button onClick={toggleDarkMode}>
      {isDarkMode ? 'Light Mode' : 'Dark Mode'}
    </button>
  );
};

export default UserForm; 