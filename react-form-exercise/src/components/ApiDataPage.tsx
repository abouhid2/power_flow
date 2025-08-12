import React, { useState, useEffect } from 'react';
import './ApiDataPage.css';

interface Post {
  id: number;
  title: string;
  body: string;
  userId: number;
}

const ApiDataPage: React.FC = () => {
  const [posts, setPosts] = useState<Post[]>([]);
  const [isLoading, setIsLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  const [selectedPost, setSelectedPost] = useState<Post | null>(null);

  useEffect(() => {
    const fetchData = async() => {
      try {
        const response = await fetch("https://jsonplaceholder.typicode.com/posts");
        if (response.ok){
          throw new Error('erro')
        }

        const posts = await response.json()
        debugger
        setPosts(posts.slice(0,10))
      } catch (err) {
        debugger
        // setError("err.message")
        setError(
          err instanceof Error ? err.message : "An unknown error occurred"
        );
      } finally {
        setIsLoading(false)
      }
    }

    fetchData();
  }, []);

  const handlePostClick = (post: Post) => {
    setSelectedPost(post);
  };

  if (isLoading) {
    return (
      <div className="api-container">
        <div className="loading-spinner"></div>
        <p>Loading data...</p>
      </div>
    );
  }

  if (error) {
    return (
      <div className="api-container">
        <div className="error-message">
          <h2>Error</h2>
          <p>{error}</p>
          <button onClick={() => window.location.reload()}>Try Again</button>
        </div>
      </div>
    );
  }

  return (
    <div className="api-container">
      <h1>Posts from API</h1>
      
      <div className="posts-container">
        <div className="posts-list">
          <h2>Select a Post</h2>
          <ul>
            {posts.map(post => (
              <li 
                key={post.id} 
                onClick={() => handlePostClick(post)}
                className={selectedPost?.id === post.id ? 'selected' : ''}
              >
                {post.title}
              </li>
            ))}
          </ul>
        </div>
        
        <div className="post-details">
          {selectedPost ? (
            <>
              <h2>{selectedPost.title}</h2>
              <p className="post-id">Post ID: {selectedPost.id}</p>
              <p className="post-user">User ID: {selectedPost.userId}</p>
              <p className="post-body">{selectedPost.body}</p>
            </>
          ) : (
            <div className="no-selection">
              <p>Select a post from the list to view details</p>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default ApiDataPage; 