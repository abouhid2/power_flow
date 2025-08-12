# Ruby on Rails Interview Test

This repository contains an incomplete Ruby on Rails application that you need to complete as part of a technical interview assessment. The test covers common Rails topics that are frequently asked in job interviews.

## Tasks

1. Complete the RSpec tests in `rails_interview_test.rb`
2. Implement the Product model in `product_model.rb`
3. Implement the ProductsController in `products_controller.rb`
4. Implement the ProductNotificationJob in `product_notification_job.rb`

## Requirements

The application should have the following features:

1. **Product Management**:

   - CRUD operations for products
   - Products belong to categories and have many reviews
   - Validation for product attributes
   - Slug generation for SEO-friendly URLs

2. **Authentication & Authorization**:

   - User authentication system
   - Role-based access control (admin vs regular users)
   - Protect certain actions based on user roles

3. **Performance Optimization**:

   - Implement eager loading to prevent N+1 queries
   - Implement caching for frequently accessed pages
   - Optimize database queries

4. **Background Processing**:
   - Use background jobs for sending email notifications
   - Implement proper error handling and retries

## Evaluation Criteria

Your solution will be evaluated based on:

1. **Code Quality**: Clean, readable, and well-organized code
2. **Test Coverage**: Comprehensive tests for all features
3. **Rails Best Practices**: Following Rails conventions and best practices
4. **Performance**: Efficient database queries and optimizations
5. **Security**: Proper authentication, authorization, and parameter sanitization

## Getting Started

1. Review the files provided in this repository
2. Implement the missing functionality in each file
3. Run the tests to ensure your implementation passes all test cases
4. Document any assumptions or design decisions you made

## Bonus Points

- Implement API endpoints with proper JSON serialization
- Add pagination for the products index
- Implement search functionality
- Add sorting and filtering options
- Implement a simple admin dashboard

Good luck!
