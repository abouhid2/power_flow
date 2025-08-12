'use client';

import { useState } from 'react';
import { Todo } from '../utils/api';

interface TodoListProps {
  todos: Todo[];
  updateTodo: (id: number, updatedTodo: Todo) => Promise<void>;
  deleteTodo: (id: number) => Promise<void>;
}

const TodoList = ({ todos, updateTodo, deleteTodo }: TodoListProps) => {
  const [editingId, setEditingId] = useState<number | null>(null);
  const [editTitle, setEditTitle] = useState('');
  const [editDescription, setEditDescription] = useState('');

  const handleEdit = (todo: Todo) => {
    setEditingId(todo.id);
    setEditTitle(todo.title);
    setEditDescription(todo.description);
  };

  const handleSave = async (todo: Todo) => {
    const updatedTodo = {
      ...todo,
      title: editTitle.trim(),
      description: editDescription.trim()
    };
    
    await updateTodo(todo.id, updatedTodo);
    setEditingId(null);
    setEditTitle('');
    setEditDescription('');
  };

  const handleCancel = () => {
    setEditingId(null);
    setEditTitle('');
    setEditDescription('');
  };

  const handleToggleComplete = async (todo: Todo) => {
    const updatedTodo = {
      ...todo,
      completed: !todo.completed
    };
    await updateTodo(todo.id, updatedTodo);
  };

  if (!todos || todos.length === 0) {
    return (
      <div className="bg-white rounded-lg shadow-md p-6">
        <p className="text-gray-500 text-center">No todos yet. Add one above!</p>
      </div>
    );
  }

  return (
    <div className="bg-white rounded-lg shadow-md p-6">
      <h2 className="text-xl font-semibold mb-4">Your Todos</h2>
      <div className="space-y-4">
        {todos.map((todo) => (
          <div
            key={todo.id}
            className={`border rounded-lg p-4 ${
              todo.completed ? 'bg-gray-50 border-gray-200' : 'bg-white border-gray-300'
            }`}
          >
            {editingId === todo.id ? (
              <div className="space-y-3">
                <input
                  type="text"
                  value={editTitle}
                  onChange={(e) => setEditTitle(e.target.value)}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                  placeholder="Todo title..."
                />
                <textarea
                  value={editDescription}
                  onChange={(e) => setEditDescription(e.target.value)}
                  className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-2 focus:ring-blue-500"
                  placeholder="Todo description..."
                  rows={2}
                />
                <div className="flex gap-2">
                  <button
                    onClick={() => handleSave(todo)}
                    className="px-4 py-2 bg-green-500 text-white rounded-md hover:bg-green-600 transition-colors"
                  >
                    Save
                  </button>
                  <button
                    onClick={handleCancel}
                    className="px-4 py-2 bg-gray-500 text-white rounded-md hover:bg-gray-600 transition-colors"
                  >
                    Cancel
                  </button>
                </div>
              </div>
            ) : (
              <div className="flex items-start justify-between">
                <div className="flex-1">
                  <div className="flex items-center gap-3">
                    <input
                      type="checkbox"
                      checked={todo.completed}
                      onChange={() => handleToggleComplete(todo)}
                      className="w-4 h-4 text-blue-600 rounded focus:ring-blue-500"
                    />
                    <h3
                      className={`text-lg font-medium ${
                        todo.completed ? 'line-through text-gray-500' : 'text-gray-900'
                      }`}
                    >
                      {todo.title}
                    </h3>
                  </div>
                  {todo.description && (
                    <p
                      className={`mt-2 text-gray-600 ${
                        todo.completed ? 'line-through' : ''
                      }`}
                    >
                      {todo.description}
                    </p>
                  )}
                  <p className="text-xs text-gray-400 mt-2">
                    Created: {new Date(todo.createdAt).toLocaleDateString()}
                  </p>
                </div>
                <div className="flex gap-2 ml-4">
                  <button
                    onClick={() => handleEdit(todo)}
                    className="px-3 py-1 bg-blue-500 text-white text-sm rounded-md hover:bg-blue-600 transition-colors"
                  >
                    Edit
                  </button>
                  <button
                    onClick={() => deleteTodo(todo.id)}
                    className="px-3 py-1 bg-red-500 text-white text-sm rounded-md hover:bg-red-600 transition-colors"
                  >
                    Delete
                  </button>
                </div>
              </div>
            )}
          </div>
        ))}
      </div>
    </div>
  );
};

export default TodoList; 