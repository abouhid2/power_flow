export interface Todo {
  id: number;
  title: string;
  description: string;
  completed: boolean;
  createdAt: string;
}

// Mock data for development
export const mockTodos: Todo[] = [
  {
    id: 1,
    title: 'Learn Next.js',
    description: 'Study the fundamentals of Next.js framework',
    completed: false,
    createdAt: new Date().toISOString()
  },
  {
    id: 2,
    title: 'Build a Todo App',
    description: 'Create a full-stack todo application',
    completed: true,
    createdAt: new Date(Date.now() - 86400000).toISOString()
  },
  {
    id: 3,
    title: 'Deploy to Vercel',
    description: 'Deploy the application to Vercel platform',
    completed: false,
    createdAt: new Date(Date.now() - 172800000).toISOString()
  }
];

// API utility functions
export const api = {
  async getTodos(): Promise<Todo[]> {
    try {
      const response = await fetch('/api/todos');
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      return await response.json();
    } catch (error) {
      console.error('Failed to fetch todos:', error);
      throw error;
    }
  },

  async createTodo(todo: Omit<Todo, 'id' | 'createdAt'>): Promise<Todo> {
    try {
      const response = await fetch('/api/todos', {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(todo),
      });
      
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      
      return await response.json();
    } catch (error) {
      console.error('Failed to create todo:', error);
      throw error;
    }
  },

  async updateTodo(id: number, updates: Partial<Todo>): Promise<Todo> {
    try {
      const response = await fetch(`/api/todos/${id}`, {
        method: 'PUT',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(updates),
      });
      
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
      
      return await response.json();
    } catch (error) {
      console.error('Failed to update todo:', error);
      throw error;
    }
  },

  async deleteTodo(id: number): Promise<void> {
    try {
      const response = await fetch(`/api/todos/${id}`, {
        method: 'DELETE',
      });
      
      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }
    } catch (error) {
      console.error('Failed to delete todo:', error);
      throw error;
    }
  }
}; 