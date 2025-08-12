import React from 'react';
import { useCart } from '../context/CartContext';
import './CartPage.css';

// Sample product data - normally this would come from an API
const availableProducts = [
  { id: 1, name: 'Product 1', price: 19.99 },
  { id: 2, name: 'Product 2', price: 29.99 },
  { id: 3, name: 'Product 3', price: 39.99 },
];

const CartPage: React.FC = () => {
  // TODO: Use the useCart hook to access cart context
  // const { items, addItem, removeItem, updateQuantity } = useCart();

  // TODO: Implement function to add product to cart
  const handleAddToCart = (product: any) => {
    // TODO: Implement this function using the context
    console.log('Adding product to cart:', product);
  };

  // TODO: Implement function to remove item from cart
  const handleRemoveItem = (itemId: number) => {
    // TODO: Implement this function using the context
    console.log('Removing item:', itemId);
  };

  // TODO: Implement function to update item quantity
  const handleUpdateQuantity = (itemId: number, newQuantity: number) => {
    // TODO: Implement this function using the context
    console.log('Updating quantity for item:', itemId, 'to', newQuantity);
  };

  // TODO: Calculate cart total
  const calculateTotal = () => {
    // TODO: Implement this function using the context
    return 0;
  };

  return (
    <div className="cart-page">
      <h2>Shopping Cart</h2>
      
      <div className="cart-container">
        <div className="products-section">
          <h3>Available Products</h3>
          <div className="product-list">
            {availableProducts.map(product => (
              <div key={product.id} className="product-item">
                <h4>{product.name}</h4>
                <p>${product.price.toFixed(2)}</p>
                <button onClick={() => handleAddToCart(product)}>
                  Add to Cart
                </button>
              </div>
            ))}
          </div>
        </div>
        
        <div className="cart-section">
          <h3>Your Cart</h3>
          {/* TODO: Display cart items here using the context */}
          <div className="cart-items">
            {/* Replace this with actual cart items from context */}
            <p>Your cart is empty</p>
            
            {/* Example of how a cart item should look:
            <div className="cart-item">
              <div className="item-info">
                <h4>Product Name</h4>
                <p>$19.99</p>
              </div>
              <div className="item-controls">
                <button>-</button>
                <span>1</span>
                <button>+</button>
                <button className="remove-btn">Remove</button>
              </div>
            </div>
            */}
          </div>
          
          <div className="cart-summary">
            <h4>Total: ${calculateTotal().toFixed(2)}</h4>
            <button className="checkout-btn" disabled={true}>
              Checkout
            </button>
          </div>
        </div>
      </div>
    </div>
  );
};

export default CartPage; 