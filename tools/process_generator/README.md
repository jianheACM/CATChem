# CATChem Process Generator

A powerful Python tool for generating standardized atmospheric chemistry process implementations for the CATChem framework.

## Features

- 🎯 **Template-driven generation** - Uses Jinja2 templates for consistent code structure
- ⚙️ **YAML configuration** - Simple configuration files define process parameters
- 🔍 **Validation** - Built-in validation for process configurations
- 🚀 **Modern Fortran** - Generates clean, modern Fortran 2008+ code
- 🔧 **CATChem integration** - Seamlessly integrates with CATChem infrastructure

## Quick Start

### 1. Set up development environment

```bash
# Quick setup - Full development environment (recommended for contributors)
./setup_dev_env.sh my-env

# Minimal setup - Just run the tool (for users)
./setup_dev_env.sh my-env 3.10 minimal

# Manual setup options:

# Option A: Full development setup with pyproject.toml
conda create -n catchem-process-gen python=3.10
conda activate catchem-process-gen
pip install -e ".[dev,test]"

# Option B: Minimal setup with requirements.txt
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt

# Option C: Manual venv with pyproject.toml
python3 -m venv venv
source venv/bin/activate
pip install -e .
```

### 2. Generate a process

```bash
# Basic usage
python process_generator.py configs/example_process.yaml

# With custom output directory
python process_generator.py configs/example_process.yaml --output-dir /path/to/output

# Validate configuration only
python process_generator.py configs/example_process.yaml --validate-only
```

### 3. Integration with CATChem

The generated processes automatically integrate with the CATChem framework:

```fortran
! The generated process includes:
! - Proper ProcessInterface implementation
! - Error handling integration
! - State variable management
! - Diagnostic capability
! - Registry-compatible factory pattern
```

## Configuration Format

Process configurations use YAML format:

```yaml
process:
  name: "MyChemicalProcess"
  description: "Custom atmospheric chemistry process"
  version: "1.0.0"

variables:
  - name: "temperature"
    type: "real"
    units: "K"
    description: "Temperature"

  - name: "concentration_o3"
    type: "real"
    units: "molecules/cm3"
    description: "Ozone concentration"

chemistry:
  reactions:
    - name: "R1"
      equation: "O + O2 + M -> O3 + M"
      rate_constant: "6.0e-34 * (300/T)**2.4"
```

## Dependencies

The process generator supports two installation approaches:

### Minimal Installation (requirements.txt)
For users who just want to run the tool:
```bash
pip install -r requirements.txt
```

**Core Dependencies:**
- **PyYAML** (≥6.0) - YAML configuration parsing
- **Jinja2** (≥3.1.0) - Template engine
- **click** (≥8.0.0) - Command-line interface
- **pathlib-mate** (≥1.0.0) - Enhanced path operations
- **colorama** (≥0.4.4) - Cross-platform colored terminal output
- **rich** (≥13.0.0) - Rich terminal formatting

### Full Development Installation (pyproject.toml)
For contributors and developers:
```bash
pip install -e ".[dev,test,docs,validation]"
```

**Additional Development Dependencies:**
- **pytest** - Testing framework
- **black** - Code formatting
- **isort** - Import sorting
- **flake8** - Linting
- **mypy** - Type checking
- **pre-commit** - Git hooks

**Optional Extensions:**
- **sphinx** - Documentation generation
- **jsonschema** - Configuration validation
- **yamllint** - YAML file validation
- **fortls** - Fortran language server support

## Project Structure

```
process_generator/
├── process_generator.py     # Main generator script
├── templates/              # Jinja2 templates
│   ├── process_interface.F90.j2
│   ├── process_creator.F90.j2
│   └── ...
├── configs/               # Example configurations
│   ├── example_process.yaml
│   └── ...
├── pyproject.toml        # Package configuration
├── setup_dev_env.sh      # Development setup script
└── README.md            # This file
```

## Development Workflow

### Testing
```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=process_generator

# Run specific test categories
pytest -m unit
pytest -m integration
pytest -m template
```

### Code Quality
```bash
# Format code
black .
isort .

# Check linting
flake8

# Type checking
mypy process_generator.py

# Pre-commit hooks
pre-commit install
pre-commit run --all-files
```

### Documentation
```bash
# Build documentation
sphinx-build docs docs/_build

# Live documentation server
sphinx-autobuild docs docs/_build
```

## Templates

The generator includes several Jinja2 templates:

- `process_interface.F90.j2` - Main process implementation
- `process_creator.F90.j2` - Factory pattern creator
- `process_config.yaml.j2` - Configuration template
- `process_test.F90.j2` - Unit test template

### Custom Templates

You can create custom templates by:

1. Adding `.j2` files to the `templates/` directory
2. Using the CATChem template variables
3. Following the established naming conventions

## Command Line Interface

```bash
# Basic usage
python process_generator.py CONFIG_FILE

# Options
  --output-dir DIR        Output directory for generated files
  --template-dir DIR      Custom template directory
  --validate-only         Only validate configuration
  --verbose              Enable verbose output
  --dry-run              Show what would be generated
  --force                Overwrite existing files
  --help                 Show help message
```

## Integration with CATChem Build System

Generated processes automatically integrate with the CATChem CMake build system:

```cmake
# Generated processes include proper CMake targets
# No manual build system modification required
```

## Troubleshooting

### Common Issues

1. **Import errors**: Ensure all dependencies are installed
   ```bash
   pip install -e ".[dev]"
   ```

2. **Template not found**: Check template directory path
   ```bash
   python process_generator.py --template-dir ./templates config.yaml
   ```

3. **YAML parsing errors**: Validate YAML syntax
   ```bash
   yamllint configs/your_config.yaml
   ```

### Getting Help

- Check the CATChem documentation: `../docs/`
- Review example configurations: `configs/`
- Run with `--verbose` for detailed output
- Use `--dry-run` to preview generation

## Contributing

1. Set up development environment: `./setup_dev_env.sh`
2. Make changes with tests: `pytest`
3. Ensure code quality: `pre-commit run --all-files`
4. Update documentation as needed

## License

Apache License 2.0 - See the main CATChem repository for details.
