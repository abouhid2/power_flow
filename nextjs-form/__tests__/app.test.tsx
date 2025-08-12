import React from 'react'
import { render, screen, fireEvent, waitFor } from '@testing-library/react'
import userEvent from '@testing-library/user-event'
import '@testing-library/jest-dom'
import Home from '../src/app/page'
import TodoForm from '../src/components/TodoForm'
import TodoList from '../src/components/TodoList'
import { api } from '../src/utils/api'

// Mock the API module
jest.mock('../src/utils/api', () => ({
  api: {
    getTodos: jest.fn(),
    createTodo: jest.fn(),
    updateTodo: jest.fn(),
    deleteTodo: jest.fn(),
  },
}))

const mockApi = api as jest.Mocked<typeof api>

describe('Todo Application', () => {
  const mockTodos = [
    {
      id: 1,
      title: 'Test Todo 1',
      description: 'Test Description 1',
      completed: false,
      createdAt: '2024-01-01T00:00:00.000Z',
    },
    {
      id: 2,
      title: 'Test Todo 2',
      description: 'Test Description 2',
      completed: true,
      createdAt: '2024-01-02T00:00:00.000Z',
    },
  ]

  beforeEach(() => {
    jest.clearAllMocks()
    // Mock the initial fetch to return our test todos
    mockApi.getTodos.mockResolvedValue(mockTodos)
  })

  describe('Home Page Component', () => {
    it('renders with initial todos', async () => {
      render(<Home />)
      
      // Wait for the component to load and fetch todos
      await waitFor(() => {
        expect(screen.getByText('Todo Application')).toBeInTheDocument()
      })
      
      expect(screen.getByText('Add New Todo')).toBeInTheDocument()
      expect(screen.getByText('Your Todos')).toBeInTheDocument()
      expect(screen.getByText('Test Todo 1')).toBeInTheDocument()
      expect(screen.getByText('Test Todo 2')).toBeInTheDocument()
    })

    it('adds a new todo successfully', async () => {
      const user = userEvent.setup()
      const newTodo = {
        id: 3,
        title: 'New Todo',
        description: 'New Description',
        completed: false,
        createdAt: '2024-01-03T00:00:00.000Z',
      }
      
      mockApi.createTodo.mockResolvedValue(newTodo)
      
      render(<Home />)
      
      // Wait for initial load
      await waitFor(() => {
        expect(screen.getByText('Todo Application')).toBeInTheDocument()
      })
      
      const titleInput = screen.getByPlaceholderText('Enter todo title...')
      const descriptionInput = screen.getByPlaceholderText('Enter todo description...')
      const submitButton = screen.getByText('Add Todo')
      
      await user.type(titleInput, 'New Todo')
      await user.type(descriptionInput, 'New Description')
      await user.click(submitButton)
      
      await waitFor(() => {
        expect(mockApi.createTodo).toHaveBeenCalledWith({
          title: 'New Todo',
          description: 'New Description',
          completed: false,
        })
      })
    })

    it('updates a todo successfully', async () => {
      const user = userEvent.setup()
      const updatedTodo = { ...mockTodos[0], title: 'Updated Todo' }
      
      mockApi.updateTodo.mockResolvedValue(updatedTodo)
      
      render(<Home />)
      
      // Wait for initial load
      await waitFor(() => {
        expect(screen.getByText('Todo Application')).toBeInTheDocument()
      })
      
      const editButtons = screen.getAllByText('Edit')
      await user.click(editButtons[0])
      
      const titleInput = screen.getByDisplayValue('Test Todo 1')
      await user.clear(titleInput)
      await user.type(titleInput, 'Updated Todo')
      
      const saveButton = screen.getByText('Save')
      await user.click(saveButton)
      
      await waitFor(() => {
        expect(mockApi.updateTodo).toHaveBeenCalledWith(1, {
          title: 'Updated Todo',
          description: 'Test Description 1',
          completed: false,
          createdAt: '2024-01-01T00:00:00.000Z',
          id: 1,
        })
      })
    })

    it('deletes a todo successfully', async () => {
      const user = userEvent.setup()
      
      mockApi.deleteTodo.mockResolvedValue()
      
      render(<Home />)
      
      // Wait for initial load
      await waitFor(() => {
        expect(screen.getByText('Todo Application')).toBeInTheDocument()
      })
      
      const deleteButtons = screen.getAllByText('Delete')
      await user.click(deleteButtons[0])
      
      await waitFor(() => {
        expect(mockApi.deleteTodo).toHaveBeenCalledWith(1)
      })
    })

    it('toggles todo completion status', async () => {
      const user = userEvent.setup()
      const updatedTodo = { ...mockTodos[0], completed: true }
      
      mockApi.updateTodo.mockResolvedValue(updatedTodo)
      
      render(<Home />)
      
      // Wait for initial load
      await waitFor(() => {
        expect(screen.getByText('Todo Application')).toBeInTheDocument()
      })
      
      const checkboxes = screen.getAllByRole('checkbox')
      await user.click(checkboxes[0])
      
      await waitFor(() => {
        expect(mockApi.updateTodo).toHaveBeenCalledWith(1, {
          completed: true,
          createdAt: '2024-01-01T00:00:00.000Z',
          description: 'Test Description 1',
          id: 1,
          title: 'Test Todo 1',
        })
      })
    })
  })

  describe('TodoForm Component', () => {
    const mockAddTodo = jest.fn()

    beforeEach(() => {
      mockAddTodo.mockClear()
    })

    it('renders form fields correctly', () => {
      render(<TodoForm addTodo={mockAddTodo} />)
      
      expect(screen.getByLabelText('Title *')).toBeInTheDocument()
      expect(screen.getByLabelText('Description')).toBeInTheDocument()
      expect(screen.getByText('Add Todo')).toBeInTheDocument()
    })

    it('submits form with valid data', async () => {
      const user = userEvent.setup()
      
      render(<TodoForm addTodo={mockAddTodo} />)
      
      const titleInput = screen.getByLabelText('Title *')
      const descriptionInput = screen.getByLabelText('Description')
      const submitButton = screen.getByText('Add Todo')
      
      await user.type(titleInput, 'Test Todo')
      await user.type(descriptionInput, 'Test Description')
      await user.click(submitButton)
      
      expect(mockAddTodo).toHaveBeenCalledWith({
        id: expect.any(Number),
        title: 'Test Todo',
        description: 'Test Description',
        completed: false,
        createdAt: expect.any(String),
      })
    })

    it('does not submit form with empty title', async () => {
      const user = userEvent.setup()
      
      render(<TodoForm addTodo={mockAddTodo} />)
      
      const descriptionInput = screen.getByLabelText('Description')
      const submitButton = screen.getByText('Add Todo')
      
      await user.type(descriptionInput, 'Test Description')
      await user.click(submitButton)
      
      expect(mockAddTodo).not.toHaveBeenCalled()
    })

    it('clears form after successful submission', async () => {
      const user = userEvent.setup()
      
      render(<TodoForm addTodo={mockAddTodo} />)
      
      const titleInput = screen.getByLabelText('Title *')
      const descriptionInput = screen.getByLabelText('Description')
      const submitButton = screen.getByText('Add Todo')
      
      await user.type(titleInput, 'Test Todo')
      await user.type(descriptionInput, 'Test Description')
      await user.click(submitButton)
      
      expect(titleInput).toHaveValue('')
      expect(descriptionInput).toHaveValue('')
    })
  })

  describe('TodoList Component', () => {
    const mockUpdateTodo = jest.fn()
    const mockDeleteTodo = jest.fn()

    beforeEach(() => {
      mockUpdateTodo.mockClear()
      mockDeleteTodo.mockClear()
    })

    it('renders todos correctly', () => {
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      expect(screen.getByText('Test Todo 1')).toBeInTheDocument()
      expect(screen.getByText('Test Todo 2')).toBeInTheDocument()
      expect(screen.getByText('Test Description 1')).toBeInTheDocument()
      expect(screen.getByText('Test Description 2')).toBeInTheDocument()
    })

    it('shows empty state when no todos', () => {
      render(
        <TodoList
          todos={[]}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      expect(screen.getByText('No todos yet. Add one above!')).toBeInTheDocument()
    })

    it('handles edit mode correctly', async () => {
      const user = userEvent.setup()
      
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      const editButtons = screen.getAllByText('Edit')
      await user.click(editButtons[0])
      
      expect(screen.getByDisplayValue('Test Todo 1')).toBeInTheDocument()
      expect(screen.getByDisplayValue('Test Description 1')).toBeInTheDocument()
      expect(screen.getByText('Save')).toBeInTheDocument()
      expect(screen.getByText('Cancel')).toBeInTheDocument()
    })

    it('saves edited todo', async () => {
      const user = userEvent.setup()
      
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      const editButtons = screen.getAllByText('Edit')
      await user.click(editButtons[0])
      
      const titleInput = screen.getByDisplayValue('Test Todo 1')
      await user.clear(titleInput)
      await user.type(titleInput, 'Updated Todo')
      
      const saveButton = screen.getByText('Save')
      await user.click(saveButton)
      
      expect(mockUpdateTodo).toHaveBeenCalledWith(1, {
        title: 'Updated Todo',
        description: 'Test Description 1',
        completed: false,
        createdAt: '2024-01-01T00:00:00.000Z',
        id: 1,
      })
    })

    it('cancels edit mode', async () => {
      const user = userEvent.setup()
      
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      const editButtons = screen.getAllByText('Edit')
      await user.click(editButtons[0])
      
      const cancelButton = screen.getByText('Cancel')
      await user.click(cancelButton)
      
      expect(screen.queryByText('Save')).not.toBeInTheDocument()
      expect(screen.queryByText('Cancel')).not.toBeInTheDocument()
    })

    it('deletes todo when delete button is clicked', async () => {
      const user = userEvent.setup()
      
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      const deleteButtons = screen.getAllByText('Delete')
      await user.click(deleteButtons[0])
      
      expect(mockDeleteTodo).toHaveBeenCalledWith(1)
    })

    it('toggles todo completion', async () => {
      const user = userEvent.setup()
      
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      const checkboxes = screen.getAllByRole('checkbox')
      await user.click(checkboxes[0])
      
      expect(mockUpdateTodo).toHaveBeenCalledWith(1, {
        completed: true,
        createdAt: '2024-01-01T00:00:00.000Z',
        description: 'Test Description 1',
        id: 1,
        title: 'Test Todo 1',
      })
    })

    it('displays completed todos with strikethrough', () => {
      render(
        <TodoList
          todos={mockTodos}
          updateTodo={mockUpdateTodo}
          deleteTodo={mockDeleteTodo}
        />
      )
      
      const completedTodo = screen.getByText('Test Todo 2')
      expect(completedTodo).toHaveClass('line-through')
    })
  })

  describe('API Integration', () => {
    it('handles API errors gracefully', async () => {
      const consoleSpy = jest.spyOn(console, 'error').mockImplementation()
      mockApi.createTodo.mockRejectedValue(new Error('API Error'))
      
      const user = userEvent.setup()
      
      render(<Home />)
      
      // Wait for initial load
      await waitFor(() => {
        expect(screen.getByText('Todo Application')).toBeInTheDocument()
      })
      
      const titleInput = screen.getByPlaceholderText('Enter todo title...')
      const submitButton = screen.getByText('Add Todo')
      
      await user.type(titleInput, 'Test Todo')
      await user.click(submitButton)
      
      await waitFor(() => {
        expect(consoleSpy).toHaveBeenCalledWith('Failed to add todo:', expect.any(Error))
      })
      
      consoleSpy.mockRestore()
    })
  })
}) 