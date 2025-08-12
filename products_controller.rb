class ProductsController < ApplicationController
  # TODO: Add authentication and authorization
  # Hint: Use before_action with appropriate methods
  
  # TODO: Add appropriate actions (index, show, new, create, edit, update, destroy)
  # Remember to implement proper error handling and redirects
  
  def index
    # TODO: Implement index action
    # Retrieve all products with proper eager loading to avoid N+1 queries
    # Implement caching for better performance
  end

  def show
    # TODO: Implement show action
    # Find the product by ID or slug
  end

  def new
    # TODO: Implement new action
    # Initialize a new product
  end

  def create
    # TODO: Implement create action
    # Create a new product with strong parameters
    # Handle success and failure cases
    # Enqueue background job for email notification
  end

  def edit
    # TODO: Implement edit action
    # Find the product by ID or slug
  end

  def update
    # TODO: Implement update action
    # Update the product with strong parameters
    # Handle success and failure cases
  end

  def destroy
    # TODO: Implement destroy action
    # Find and destroy the product
    # Ensure only admin users can perform this action
  end

  private

  # TODO: Implement strong parameters
  # Hint: Use params.require(:product).permit(...)
  
  # TODO: Implement authentication check
  # Hint: Redirect to login page if user is not authenticated
  
  # TODO: Implement authorization check
  # Hint: Check if current user has admin role for certain actions
  
  # TODO: Implement find_product method
  # Hint: Use find_by for more flexibility
end 