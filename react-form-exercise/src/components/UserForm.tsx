import React, { useState, useEffect } from 'react';
import './UserForm.css';

interface UserData {
  name: string;
  email: string;
  age: string;
  occupation: string;
}

// This form is intentionally incomplete. 
// TODO: Complete the form implementation by:
// 1. Adding proper form validation
// 2. Implementing the handleSubmit function
// 3. Adding error messages for invalid inputs
// 4. Adding proper styling

const UserForm: React.FC = () => {
  // Initial state with empty values
  const [userData, setUserData] = useState<UserData>({
    name: '',
    email: '',
    age: '',
    occupation: ''
  });
  const [errors, setErrors] = useState<Record<string, string>>({});
  const [isSubmitting, setIsSubmitting] = useState(false);
  const [isSubmitted, setIsSubmitted] = useState(false);

  // Validate all fields on component mount
  useEffect(() => {
    const initialErrors: Record<string, string> = {};
    
    // Validate name
    if (!userData.name.trim()) {
      initialErrors.name = 'Name is required';
    } else if (userData.name.length < 8) {
      initialErrors.name = 'Name must be at least 8 characters long';
    }
    
    // Validate email
    const emailRegex = /^[^\s@]+@[^\s@]+\.[^\s@]+$/;
    if (!userData.email.trim()) {
      initialErrors.email = 'Email is required';
    } else if (!emailRegex.test(userData.email)) {
      initialErrors.email = 'Please enter a valid email address';
    }
    
    // Validate age
    const ageValue = parseInt(userData.age);
    if (!userData.age.trim()) {
      initialErrors.age = 'Age is required';
    } else if (isNaN(ageValue) || ageValue < 18 || ageValue > 100) {
      initialErrors.age = 'Age must be between 18 and 100';
    }
    
    // Validate occupation
    if (!userData.occupation.trim()) {
      initialErrors.occupation = 'Occupation is required';
    } else if (userData.occupation.length < 3) {
      initialErrors.occupation = 'Occupation must be at least 3 characters long';
    }
    
    setErrors(initialErrors);
  }, []);

  const validateField = (name: string, value: string) => {
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
  };
  
  // This function updates the state when inputs change
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;

    validateField(name, value);
    setUserData({
      ...userData,
      [name]: value
    });
  };
  
  // Validate the entire form
  const validateForm = (): boolean => {
    // Validate all fields
    validateField('name', userData.name);
    validateField('email', userData.email);
    validateField('age', userData.age);
    validateField('occupation', userData.occupation);
    
    // Check if there are any errors
    return Object.keys(errors).length === 0;
  };
  
  // Handle form submission
  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    
    // Validate the form
    const isValid = validateForm();
    
    if (isValid) {
      setIsSubmitting(true);
      
      // Simulate API call with setTimeout
      setTimeout(() => {
        console.log('Form submitted with:', userData);
        setIsSubmitting(false);
        setIsSubmitted(true);
        
        // Reset form after successful submission
        setUserData({
          name: '',
          email: '',
          age: '',
          occupation: ''
        });
      }, 1500);
    } else {
      setIsSubmitted(false)
    }
  };
  
  return (
    <div className="form-container">
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
            type="text"
            id="occupation"
            name="occupation"
            value={userData.occupation}
            onChange={handleChange}
          />
          {errors["occupation"] && <p className="error-message">{errors["occupation"]}</p>}
        </div>

        <button type="submit" disabled={isSubmitting}>
          {isSubmitting ? 'Submitting...' : 'Submit'}
        </button>
      </form>
    </div>
  );
};

export default UserForm; 