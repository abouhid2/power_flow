import React, { useState } from 'react';
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


  const validField = (name: string, value: string) => {

    return true
  }
  
  // TODO: Add form validation state
  const [errors, setErrors] = useState<Record<string, string>>({});
  
  // TODO: Add submission status state
  // const [isSubmitting, setIsSubmitting] = useState(false);
  // const [isSubmitted, setIsSubmitted] = useState(false);
  
  // This function updates the state when inputs change
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;

    if (validField(name, value)) {
      setUserData({
        ...userData,
        [name]: value
      });
      
    }
    // TODO: Clear errors when user types
  };
  
  // TODO: Implement form validation
  // const validateForm = (): boolean => {
  //   let isValid = true;
  //   const newErrors: Record<string, string> = {};
  //
  //   // Add validation logic here
  //
  //   // setErrors(newErrors);
  //   return isValid;
  // };
  
  // TODO: Implement form submission
  const handleSubmit = (e: React.FormEvent) => {
    e.preventDefault();
    
    // TODO: Add validation before submission
    
    // TODO: Simulate API call with setTimeout
    console.log('Form submitted with:', userData);
    
    // TODO: Reset form after successful submission
  };
  
  return (
    <div className="form-container">
      <h2>User Registration Form</h2>
      
      {/* TODO: Add success message when form is submitted successfully */}
      
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
          {/* TODO: Add error message for name */}
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
          {/* TODO: Add error message for email */}
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
          {/* TODO: Add error message for age */}
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
          {/* TODO: Add error message for occupation */}
        </div>
        
        <button type="submit">Submit</button>
      </form>
    </div>
  );
};

export default UserForm; 