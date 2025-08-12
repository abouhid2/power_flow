class ProductNotificationJob < ApplicationJob
  queue_as :default

  # TODO: Implement perform method
  # This job should send an email notification when a new product is created
  # Hint: Use ActionMailer to send emails
  def perform(product_id)
    # TODO: Find the product by ID
    
    # TODO: Send email notification to admin users
    # Hint: Use a mailer class
    
    # TODO: Implement retry mechanism for failed jobs
    # Hint: Use rescue with retry_job
    
    # TODO: Log job completion
    # Hint: Use Rails.logger
  end
end 