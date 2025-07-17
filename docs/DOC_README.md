# CATChem Documentation

This directory contains the source files for CATChem's documentation website, built with [MkDocs](https://www.mkdocs.org/) and the [Material theme](https://squidfunk.github.io/mkdocs-material/).

## Quick Start

### Prerequisites

```bash
# Install Python 3.8+
python --version

# Install documentation dependencies
pip install -r docs/requirements.txt
```

### Development Server

```bash
# Start live development server
mkdocs serve

# Open browser to http://localhost:8000
```

### Building Documentation

```bash
# Build static site
mkdocs build

# Build with API documentation
mkdocs build --config-file mkdocs.yml
```

## Structure

```
docs/
├── mkdocs.yml              # Main configuration
├── requirements.txt        # Python dependencies
├── index.md               # Homepage
├── assets/                # Static assets
│   ├── stylesheets/       # NOAA theme CSS
│   ├── javascripts/       # Custom JS
│   └── images/           # Images and logos
├── quick-start/          # Quick start guide
├── user-guide/           # User documentation
│   └── processes/        # Process-specific docs
├── developer-guide/      # Developer documentation
│   ├── processes/        # Process development
│   ├── core/            # Core system docs
│   └── integration/     # Integration guides
├── api/                 # Auto-generated API docs
└── guides/              # Technical guides
```

## Features

### Material Theme with NOAA Colors
- Custom NOAA color palette (blue, teal, navy)
- Responsive design for all devices
- Dark/light mode toggle
- Search functionality

### MkDoxy Integration
- Automatic API documentation from Fortran source code
- Doxygen-style comment parsing
- Cross-references between code and docs

### Enhanced Markdown
- Mermaid diagrams for architecture visualization
- Code syntax highlighting with Fortran support
- Admonitions for notes, warnings, tips
- Tabbed content sections
- Mathematical equations with MathJax

## Contributing

Update documentation alongside code changes and test locally with `mkdocs serve`.

For detailed information, see the full documentation setup guide in the built documentation.
