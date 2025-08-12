"use client";

import { useEffect, useState } from "react";
import TodoList from "../components/TodoList";
import TodoForm from "../components/TodoForm";
import { Todo, api } from "../utils/api";

const Home = () => {
  const [todos, setTodos] = useState<Todo[]>([]);
  const [loading, setLoading] = useState(false);

  // Fetch todos on load
  useEffect(() => {
    const fetchTodos = async () => {
      setLoading(true);
      try {
        const data = await api.getTodos();
        setTodos(data);
      } catch (error) {
        console.error("Failed to fetch todos:", error);
      } finally {
        setLoading(false);
        // setLoading(false);
      }
    };
    fetchTodos();
  }, []);

  const addTodo = async (newTodo: Todo) => {
    try {
      const createdTodo = await api.createTodo({
        title: newTodo.title,
        description: newTodo.description,
        completed: newTodo.completed,
      });
      setTodos((prev) => [...prev, createdTodo]);
    } catch (error) {
      console.error("Failed to add todo:", error);
      // Optionally show user feedback here
    }
  };

  const updateTodo = async (id: number, updatedTodo: Todo) => {
    try {
      const result = await api.updateTodo(id, updatedTodo);
      setTodos((prev) => prev.map((todo) => (todo.id === id ? result : todo)));
    } catch (error) {
      console.error("Failed to update todo:", error);
      // Optionally show user feedback here
    }
  };

  const deleteTodo = async (id: number) => {
    try {
      await api.deleteTodo(id);
      setTodos((prev) => prev.filter((todo) => todo.id !== id));
    } catch (error) {
      console.error("Failed to delete todo:", error);
      // Optionally show user feedback here
    }
  };

  if (loading) {
    return (
      <div className="container mx-auto p-4">
        <div className="text-center">Loading todos...</div>
      </div>
    );
  }

  return (
    <div className="container mx-auto p-4">
      <h1 className="text-2xl font-bold">Todo Application</h1>
      <TodoForm addTodo={addTodo} />
      <TodoList todos={todos} updateTodo={updateTodo} deleteTodo={deleteTodo} />
    </div>
  );
};

export default Home;
