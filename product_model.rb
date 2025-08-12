class Product < ApplicationRecord
  # TODO: Add associations
  # A product belongs to a category
  # A product has many reviews
  
  # TODO: Add validations
  # Title must be present and unique
  # Price must be present and numeric
  # Description must be present with minimum length of 10 characters

  # TODO: Add callback to generate slug before saving
  # Hint: Use parameterize on the title to create the slug
  
  # TODO: Add a custom scope to return active products
  # Hint: Products with active field set to true

  # TODO: Implement a method to calculate the average rating from reviews
  # Hint: Use reviews association

  # TODO: Implement a method to check if the product is on sale
  # Hint: Compare sale_price with regular_price

  # TODO: Implement a method to format the price with currency symbol
  # Hint: Use number_to_currency helper or similar approach
end 