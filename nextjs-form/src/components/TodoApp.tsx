'use client';

import { useEffect, useState } from 'react';
import TodoList from './TodoList';
import TodoForm from './TodoForm';
import { Todo, api } from '../utils/api';

interface TodoAppProps {
  initialTodos: Todo[];
}

const TodoApp = ({ initialTodos }: TodoAppProps) => {
  const [todos, setTodos] = useState<Todo[]>(initialTodos);
  const [loading, setLoading] = useState(false);

  // Refresh todos from server
  const refreshTodos = async () => {
    setLoading(true);
    try {
      const data = await api.getTodos();
      setTodos(data);
    } catch (error) {
      console.error('Failed to refresh todos:', error);
    } finally {
      setLoading(false);
    }
  };

  const addTodo = async (newTodo: Todo) => {
    try {
      const createdTodo = await api.createTodo({
        title: newTodo.title,
        description: newTodo.description,
        completed: newTodo.completed
      });
      setTodos(prev => [...prev, createdTodo]);
    } catch (error) {
      console.error('Failed to add todo:', error);
      // Refresh todos to ensure consistency
      await refreshTodos();
    }
  };

  const updateTodo = async (id: number, updatedTodo: Todo) => {
    try {
      const result = await api.updateTodo(id, updatedTodo);
      setTodos(prev => prev.map(todo => todo.id === id ? result : todo));
    } catch (error) {
      console.error('Failed to update todo:', error);
      // Refresh todos to ensure consistency
      await refreshTodos();
    }
  };

  const deleteTodo = async (id: number) => {
    try {
      await api.deleteTodo(id);
      setTodos(prev => prev.filter(todo => todo.id !== id));
    } catch (error) {
      console.error('Failed to delete todo:', error);
      // Refresh todos to ensure consistency
      await refreshTodos();
    }
  };

  if (loading) {
    return (
      <div className="flex items-center justify-center py-8">
        <div className="text-lg">Refreshing todos...</div>
      </div>
    );
  }

  return (
    <>
      <TodoForm addTodo={addTodo} />
      <TodoList todos={todos} updateTodo={updateTodo} deleteTodo={deleteTodo} />
    </>
  );
};

export default TodoApp; 