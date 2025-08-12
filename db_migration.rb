class CreateProductTables < ActiveRecord::Migration[6.1]
  def change
    # TODO: Create categories table
    # Hint: Should have name, description, and slug fields
    
    # TODO: Create products table
    # Hint: Should have title, description, price, sale_price, active, slug, and category_id fields
    
    # TODO: Create reviews table
    # Hint: Should have product_id, user_id, rating, comment, and approved fields
    
    # TODO: Add appropriate indexes
    # Hint: Consider indexing foreign keys, slug fields, and frequently queried columns
    
    # TODO: Add appropriate foreign key constraints
    # Hint: Use references with foreign_key: true
  end
end 