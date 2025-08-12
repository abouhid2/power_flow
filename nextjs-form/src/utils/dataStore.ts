import { Todo } from './api';

// Shared in-memory storage for todos (in a real app, this would be a database)
const todoData = {
  todos: [
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
  ] as Todo[]
};

// Helper functions to manage todos
export const todoStore = {
  getAll: () => {
    console.log('Getting all todos:', todoData.todos);
    return todoData.todos;
  },
  
  getById: (id: number) => {
    const todo = todoData.todos.find(todo => todo.id === id);
    console.log(`Getting todo by id ${id}:`, todo);
    return todo;
  },
  
  add: (todo: Omit<Todo, 'id' | 'createdAt'>) => {
    const newTodo: Todo = {
      id: Date.now(),
      title: todo.title,
      description: todo.description,
      completed: todo.completed,
      createdAt: new Date().toISOString()
    };
    console.log('Adding new todo:', newTodo);
    todoData.todos.push(newTodo);
    console.log('Todos after adding:', todoData.todos);
    return newTodo;
  },
  
  update: (id: number, updates: Partial<Todo>) => {
    console.log(`Updating todo ${id} with:`, updates);
    const index = todoData.todos.findIndex(todo => todo.id === id);
    if (index === -1) {
      console.log(`Todo ${id} not found`);
      return null;
    }
    
    todoData.todos[index] = { ...todoData.todos[index], ...updates };
    console.log('Updated todo:', todoData.todos[index]);
    return todoData.todos[index];
  },
  
  delete: (id: number) => {
    console.log(`Deleting todo ${id}`);
    const index = todoData.todos.findIndex(todo => todo.id === id);
    if (index === -1) {
      console.log(`Todo ${id} not found for deletion`);
      return false;
    }
    
    todoData.todos.splice(index, 1);
    console.log('Todos after deletion:', todoData.todos);
    return true;
  }
};

// Export todos for backward compatibility if needed
export const todos = todoData.todos; 