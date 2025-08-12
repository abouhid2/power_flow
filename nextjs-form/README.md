# Fetchly Labs Technical Challenge - Todo Application

A modern, full-stack Todo application built with Next.js, TypeScript, and Tailwind CSS. This application demonstrates both client-side and server-side rendering, comprehensive testing, and responsive design principles.

## 🚀 Features

### Core Functionality
- ✅ **Create todos** with title and description
- ✅ **Update todos** with inline editing
- ✅ **Delete todos** with confirmation
- ✅ **Mark todos as complete/incomplete** with checkboxes
- ✅ **Real-time UI updates** with optimistic rendering

### Technical Features
- 🔄 **Server-Side Rendering (SSR)** for initial data loading
- 🎨 **Responsive design** using Tailwind CSS
- 🧪 **Comprehensive testing** with Jest and React Testing Library
- 📱 **Mobile-first design** with beautiful UI/UX
- 🔒 **Type safety** with TypeScript
- ⚡ **Fast performance** with Next.js optimizations

## 🛠️ Tech Stack

- **Framework**: Next.js 15 with App Router
- **Language**: TypeScript
- **Styling**: Tailwind CSS
- **Testing**: Jest + React Testing Library
- **API**: Next.js API Routes
- **State Management**: React Hooks

## 📁 Project Structure

```
nextjs-form/
├── src/
│   ├── app/
│   │   ├── api/
│   │   │   └── todos/
│   │   │       ├── route.ts              # GET, POST /api/todos
│   │   │       └── [id]/
│   │   │           └── route.ts          # PUT, DELETE /api/todos/[id]
│   │   ├── globals.css                   # Global styles
│   │   ├── layout.tsx                    # Root layout
│   │   └── page.tsx                      # Main page (SSR)
│   ├── components/
│   │   ├── TodoApp.tsx                   # Main app component
│   │   ├── TodoForm.tsx                  # Add todo form
│   │   └── TodoList.tsx                  # Todo list with CRUD
│   └── utils/
│       └── api.ts                        # API utilities and types
├── __tests__/
│   └── app.test.tsx                      # Comprehensive tests
├── jest.config.js                        # Jest configuration
├── jest.setup.js                         # Test setup
└── package.json
```

## 🚀 Getting Started

### Prerequisites
- Node.js 18+ 
- npm or yarn

### Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd nextjs-form
   ```

2. **Install dependencies**
   ```bash
   npm install
   ```

3. **Run the development server**
   ```bash
   npm run dev
   ```

4. **Open your browser**
   Navigate to [http://localhost:3000](http://localhost:3000)

### Available Scripts

- `npm run dev` - Start development server
- `npm run build` - Build for production
- `npm run start` - Start production server
- `npm run lint` - Run ESLint
- `npm test` - Run tests
- `npm run test:watch` - Run tests in watch mode
- `npm run test:coverage` - Run tests with coverage

## 🧪 Testing

The application includes comprehensive tests covering:

- **Component rendering** - All components render correctly
- **User interactions** - Form submissions, editing, deleting
- **API integration** - CRUD operations work correctly
- **Error handling** - Graceful error handling
- **Edge cases** - Empty states, validation

Run tests:
```bash
npm test
```

Run tests with coverage:
```bash
npm run test:coverage
```

## 📡 API Endpoints

The application includes a complete REST API:

### GET /api/todos
Retrieve all todos
```json
[
  {
    "id": 1,
    "title": "Learn Next.js",
    "description": "Study the fundamentals",
    "completed": false,
    "createdAt": "2024-01-01T00:00:00.000Z"
  }
]
```

### POST /api/todos
Create a new todo
```json
{
  "title": "New Todo",
  "description": "Optional description"
}
```

### PUT /api/todos/[id]
Update a specific todo
```json
{
  "title": "Updated Todo",
  "description": "Updated description",
  "completed": true
}
```

### DELETE /api/todos/[id]
Delete a specific todo

## 🎨 UI/UX Features

### Design Principles
- **Clean and modern** interface with gradient background
- **Responsive design** that works on all devices
- **Accessible** with proper ARIA labels and keyboard navigation
- **Smooth animations** and transitions
- **Intuitive interactions** with clear visual feedback

### Components
- **TodoForm**: Clean form for adding new todos
- **TodoList**: Card-based layout with inline editing
- **TodoApp**: Main container with state management

## 🔧 Architecture

### Server-Side Rendering
- Initial page load fetches todos server-side
- Improves SEO and initial load performance
- Hydrates to client-side interactivity

### State Management
- React hooks for local state
- Optimistic updates for better UX
- Error handling with fallback to server refresh

### API Design
- RESTful endpoints following best practices
- Proper error handling and validation
- Simulated delays for realistic testing

## 🚀 Deployment

### Vercel (Recommended)
1. Push code to GitHub
2. Connect repository to Vercel
3. Deploy automatically

### Other Platforms
```bash
npm run build
npm run start
```

## 📝 Future Enhancements

- [ ] **Authentication** - User accounts and sessions
- [ ] **Database integration** - PostgreSQL or MongoDB
- [ ] **Real-time updates** - WebSocket integration
- [ ] **Categories/tags** - Organize todos
- [ ] **Due dates** - Set deadlines and reminders
- [ ] **Search/filter** - Find specific todos
- [ ] **Dark mode** - Theme switching
- [ ] **Offline support** - PWA capabilities

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request

## 📄 License

This project is part of the Fetchly Labs Technical Challenge.

---

**Built with ❤️ using Next.js, TypeScript, and Tailwind CSS**
