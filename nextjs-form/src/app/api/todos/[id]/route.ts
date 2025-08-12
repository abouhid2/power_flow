import { NextRequest, NextResponse } from 'next/server';
import { todoStore } from '../../../../utils/dataStore';

// PUT /api/todos/[id] - Update a specific todo
export async function PUT(
  request: NextRequest,
  { params }: { params: { id: string } }
) {
  try {
    const id = parseInt(params.id);
    const body = await request.json();
    const { title, description, completed } = body;

    // Check if todo exists
    const existingTodo = todoStore.getById(id);
    if (!existingTodo) {
      return NextResponse.json(
        { error: 'Todo not found' },
        { status: 404 }
      );
    }

    // Validate input
    if (title !== undefined && (typeof title !== 'string' || !title.trim())) {
      return NextResponse.json(
        { error: 'Title must be a non-empty string' },
        { status: 400 }
      );
    }

    // Update the todo
    const updatedTodo = todoStore.update(id, {
      title: title?.trim(),
      description: description?.trim(),
      completed
    });

    if (!updatedTodo) {
      return NextResponse.json(
        { error: 'Failed to update todo' },
        { status: 500 }
      );
    }

    // Simulate API delay
    await new Promise(resolve => setTimeout(resolve, 300));

    return NextResponse.json(updatedTodo, { status: 200 });
  } catch (error) {
    return NextResponse.json(
      { error: 'Failed to update todo' },
      { status: 500 }
    );
  }
}

// DELETE /api/todos/[id] - Delete a specific todo
export async function DELETE(
  request: NextRequest,
  { params }: { params: { id: string } }
) {
  try {
    const id = parseInt(params.id);
    
    // Check if todo exists
    const existingTodo = todoStore.getById(id);
    if (!existingTodo) {
      return NextResponse.json(
        { error: 'Todo not found' },
        { status: 404 }
      );
    }

    // Delete the todo
    const deleted = todoStore.delete(id);
    
    if (!deleted) {
      return NextResponse.json(
        { error: 'Failed to delete todo' },
        { status: 500 }
      );
    }

    // Simulate API delay
    await new Promise(resolve => setTimeout(resolve, 300));

    return NextResponse.json(
      { message: 'Todo deleted successfully' },
      { status: 200 }
    );
  } catch (error) {
    return NextResponse.json(
      { error: 'Failed to delete todo' },
      { status: 500 }
    );
  }
} 