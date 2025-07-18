import React from 'react';
import { render, screen, fireEvent } from '@testing-library/react';
import UserForm from './UserForm';

describe('UserForm Component', () => {
  test('renders the form with all fields', () => {
    render(<UserForm />);
    
    // Check if form title is rendered
    expect(screen.getByText('User Registration Form')).toBeInTheDocument();
    
    // Check if all form fields are rendered
    expect(screen.getByLabelText(/name:/i)).toBeInTheDocument();
    expect(screen.getByLabelText(/email:/i)).toBeInTheDocument();
    expect(screen.getByLabelText(/age:/i)).toBeInTheDocument();
    expect(screen.getByLabelText(/occupation:/i)).toBeInTheDocument();
    
    // Check if submit button is rendered
    expect(screen.getByRole('button', { name: /submit/i })).toBeInTheDocument();
  });
  
  test('allows entering values in all fields', () => {
    render(<UserForm />);
    
    // Get all input fields
    const nameInput = screen.getByLabelText(/name:/i);
    const emailInput = screen.getByLabelText(/email:/i);
    const ageInput = screen.getByLabelText(/age:/i);
    const occupationInput = screen.getByLabelText(/occupation:/i);
    
    // Enter values in all fields
    fireEvent.change(nameInput, { target: { value: 'John Doe' } });
    fireEvent.change(emailInput, { target: { value: 'john@example.com' } });
    fireEvent.change(ageInput, { target: { value: '30' } });
    fireEvent.change(occupationInput, { target: { value: 'Developer' } });
    
    // Check if values are updated
    expect(nameInput).toHaveValue('John Doe');
    expect(emailInput).toHaveValue('john@example.com');
    expect(ageInput).toHaveValue(30);
    expect(occupationInput).toHaveValue('Developer');
  });
  
  test('submits the form with user data', () => {
    // Mock console.log to check if form data is logged
    const consoleSpy = jest.spyOn(console, 'log');
    
    render(<UserForm />);
    
    // Get all input fields
    const nameInput = screen.getByLabelText(/name:/i);
    const emailInput = screen.getByLabelText(/email:/i);
    const ageInput = screen.getByLabelText(/age:/i);
    const occupationInput = screen.getByLabelText(/occupation:/i);
    
    // Enter values in all fields
    fireEvent.change(nameInput, { target: { value: 'John Doe' } });
    fireEvent.change(emailInput, { target: { value: 'john@example.com' } });
    fireEvent.change(ageInput, { target: { value: '30' } });
    fireEvent.change(occupationInput, { target: { value: 'Developer' } });
    
    // Submit the form
    fireEvent.click(screen.getByRole('button', { name: /submit/i }));
    
    // Check if form data is logged
    expect(consoleSpy).toHaveBeenCalledWith('Form submitted with:', {
      name: 'John Doe',
      email: 'john@example.com',
      age: '30',
      occupation: 'Developer'
    });
    
    // Clean up
    consoleSpy.mockRestore();
  });
  
  // TODO: Add more tests for form validation and error messages
  // once the form validation is implemented
}); 