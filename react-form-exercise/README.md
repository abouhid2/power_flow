# React Form Exercise

This project contains an incomplete React form component that you need to complete as an exercise to practice React form handling.

## Getting Started

1. Clone this repository
2. Run `npm install` to install dependencies
3. Run `npm start` to start the development server
4. Open [http://localhost:3000](http://localhost:3000) to view it in the browser

## Exercise Instructions

Your task is to complete the implementation of the `UserForm` component located in `src/components/UserForm.tsx`. The form is intentionally left incomplete for you to practice form handling in React.

### Tasks to Complete:

1. **Add Form Validation**

   - Implement the `validateForm` function to validate all form fields
   - Name should be required and at least 2 characters long
   - Email should be required and in a valid email format
   - Age should be required and between 18 and 100
   - Occupation should be required

2. **Display Error Messages**

   - Add error messages for invalid inputs
   - Show error messages below each input field
   - Style error messages appropriately

3. **Implement Form Submission**

   - Complete the `handleSubmit` function
   - Add validation before submission
   - Show loading state during submission (simulate API call with setTimeout)
   - Show success message after successful submission
   - Reset form after successful submission

4. **Add Proper Styling**
   - Style the form for different states (normal, error, success)
   - Style the submit button for different states (normal, disabled during submission)
   - Make the form responsive

## Tips

- Use the React useState hook to manage form state
- Use conditional rendering to show/hide error and success messages
- Use CSS classes to style different form states
- Test your form with different inputs to ensure validation works correctly

## Bonus Challenges

- Add a "Reset" button to clear the form
- Implement a confirmation dialog before submission
- Add client-side form validation using a library like Yup or Zod
- Add animations for form transitions

Good luck!
