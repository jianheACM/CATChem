#!/usr/bin/env python3
"""
CATChem Process Generator Tool

A comprehensive tool for generating standardized process implementations in CATChem.
Uses YAML configuration files and Jinja2 templates to create complete process
modules, schemes, documentation, tests, and CMake integration.

This tool follows the Process Infrastructure Guide and creates processes
that are compatible with the modern CATChem architecture.

Usage:
    python process_generator.py generate --config my_process.yaml
    python process_generator.py validate --config my_process.yaml
    python process_generator.py template --type process --output template.yaml

Author: CATChem Development Team
License: Apache 2.0
"""

import argparse
import logging
import os
import sys
import yaml
from pathlib import Path
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, field
from jinja2 import Environment, FileSystemLoader, select_autoescape
import json
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('ProcessGenerator')


@dataclass
class SchemeConfig:
    """Configuration for a process scheme."""
    name: str
    class_name: str
    description: str
    author: str = ""
    reference: str = ""
    parameters: Dict[str, Any] = field(default_factory=dict)
    required_met_fields: List[str] = field(default_factory=list)
    diagnostics: List[Dict[str, Any]] = field(default_factory=list)
    algorithm_type: str = "explicit"  # explicit, implicit, mixed


@dataclass
class ProcessConfig:
    """Configuration for a complete process."""
    # Basic metadata
    name: str
    class_name: str
    description: str
    author: str
    version: str = "1.0.0"
    license: str = "Apache 2.0"

    # Process characteristics
    process_type: str = "emission"  # emission, chemistry, transport, etc.
    is_multiphase: bool = False
    has_size_bins: bool = False
    supports_vectorization: bool = True

    # Species and chemistry
    species: List[str] = field(default_factory=list)
    size_bins: Optional[Dict[str, Any]] = None
    phases: List[str] = field(default_factory=lambda: ["gas"])

    # Schemes
    schemes: List[SchemeConfig] = field(default_factory=list)
    default_scheme: str = ""

    # Input/output
    required_met_fields: List[str] = field(default_factory=list)
    optional_met_fields: List[str] = field(default_factory=list)
    required_chem_fields: List[str] = field(default_factory=list)

    # Diagnostics
    diagnostics: List[Dict[str, Any]] = field(default_factory=list)

    # Integration settings
    timestep_dependency: str = "independent"  # independent, dependent, adaptive
    parallelization: str = "column"  # none, column, domain
    memory_requirements: str = "low"  # low, medium, high

    # File generation options
    generate_tests: bool = True
    generate_docs: bool = True
    generate_examples: bool = True

    # Output configuration
    output_dir: str = ""
    src_base_dir: str = "src/process"


class ProcessValidationError(Exception):
    """Exception raised for process configuration validation errors."""
    pass


class ProcessGenerator:
    """Main process generator class."""

    def __init__(self, template_dir: Optional[str] = None):
        """Initialize the process generator.

        Args:
            template_dir: Directory containing Jinja2 templates. If None,
                         uses default templates in same directory as this script.
        """
        if template_dir is None:
            template_dir = str(Path(__file__).parent / "templates")

        self.template_dir = Path(template_dir)
        self.env = Environment(
            loader=FileSystemLoader(str(self.template_dir)),
            autoescape=select_autoescape(['html', 'xml']),
            trim_blocks=True,
            lstrip_blocks=True
        )

        # Add custom filters
        self.env.filters['upper_snake'] = self._upper_snake_case
        self.env.filters['lower_snake'] = self._lower_snake_case
        self.env.filters['pascal_case'] = self._pascal_case
        self.env.filters['camel_case'] = self._camel_case
        self.env.filters['fortran_string'] = self._fortran_string
        self.env.filters['fortran_boolean'] = self._fortran_boolean

    @staticmethod
    def _upper_snake_case(s: str) -> str:
        """Convert string to UPPER_SNAKE_CASE."""
        return s.upper().replace(' ', '_').replace('-', '_')

    @staticmethod
    def _lower_snake_case(s: str) -> str:
        """Convert string to lower_snake_case."""
        return s.lower().replace(' ', '_').replace('-', '_')

    @staticmethod
    def _pascal_case(s: str) -> str:
        """Convert string to PascalCase."""
        return ''.join(word.capitalize() for word in s.replace('_', ' ').replace('-', ' ').split())

    @staticmethod
    def _camel_case(s: str) -> str:
        """Convert string to camelCase."""
        pascal = ProcessGenerator._pascal_case(s)
        return pascal[0].lower() + pascal[1:] if pascal else ""

    @staticmethod
    def _fortran_string(s: str, length: int = 64) -> str:
        """Format string for Fortran character declaration."""
        return f"'{s}'"

    @staticmethod
    def _fortran_boolean(b: bool) -> str:
        """Convert boolean to Fortran logical."""
        return ".true." if b else ".false."

    def validate_config(self, config: ProcessConfig) -> None:
        """Validate process configuration.

        Args:
            config: Process configuration to validate

        Raises:
            ProcessValidationError: If configuration is invalid
        """
        errors = []

        # Basic validation
        if not config.name:
            errors.append("Process name is required")

        if not config.class_name:
            errors.append("Process class name is required")

        if not config.description:
            errors.append("Process description is required")

        if not config.author:
            errors.append("Author name is required")

        # Name validation
        if not config.name.replace('_', '').replace('-', '').isalnum():
            errors.append("Process name must be alphanumeric (with _ or - allowed)")

        if not config.class_name.replace('_', '').isalnum():
            errors.append("Class name must be alphanumeric (with _ allowed)")

        # Scheme validation
        if not config.schemes:
            errors.append("At least one scheme must be defined")

        scheme_names = [scheme.name for scheme in config.schemes]
        if len(set(scheme_names)) != len(scheme_names):
            errors.append("Scheme names must be unique")

        if config.default_scheme and config.default_scheme not in scheme_names:
            errors.append(f"Default scheme '{config.default_scheme}' not found in schemes")

        # Species validation
        if config.species and len(set(config.species)) != len(config.species):
            errors.append("Species names must be unique")

        # Size bin validation
        if config.has_size_bins and not config.size_bins:
            errors.append("Size bin configuration required when has_size_bins=True")

        if config.size_bins:
            required_keys = ['n_bins', 'type', 'bounds']
            missing_keys = [key for key in required_keys if key not in config.size_bins]
            if missing_keys:
                errors.append(f"Size bin config missing keys: {missing_keys}")

        # Output directory validation
        if config.output_dir:
            output_path = Path(config.output_dir)
            if output_path.exists() and not output_path.is_dir():
                errors.append(f"Output path exists but is not a directory: {config.output_dir}")

        if errors:
            raise ProcessValidationError("\n".join(errors))

    def load_config(self, config_path: Union[str, Path]) -> ProcessConfig:
        """Load and validate process configuration from YAML file.

        Args:
            config_path: Path to YAML configuration file

        Returns:
            Validated ProcessConfig object

        Raises:
            ProcessValidationError: If configuration is invalid
            FileNotFoundError: If config file doesn't exist
            yaml.YAMLError: If YAML parsing fails
        """
        config_path = Path(config_path)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")

        with open(config_path, 'r') as f:
            data = yaml.safe_load(f)

        # Convert schemes to SchemeConfig objects
        schemes = []
        schemes_data = data.get('schemes', {})
        if isinstance(schemes_data, dict):
            # Schemes defined as dict with keys
            for scheme_name, scheme_data in schemes_data.items():
                # Add the scheme name if not already present
                if 'name' not in scheme_data:
                    scheme_data['name'] = scheme_name
                scheme = SchemeConfig(**scheme_data)
                schemes.append(scheme)
        elif isinstance(schemes_data, list):
            # Schemes defined as list
            for scheme_data in schemes_data:
                scheme = SchemeConfig(**scheme_data)
                schemes.append(scheme)

        # Remove schemes from data and create ProcessConfig
        data['schemes'] = schemes
        config = ProcessConfig(**data)

        # Validate configuration
        self.validate_config(config)

        return config

    def generate_process(self, config: ProcessConfig) -> None:
        """Generate complete process implementation.

        Args:
            config: Validated process configuration
        """
        logger.info(f"Generating process: {config.name}")

        # Determine output directory
        if config.output_dir:
            base_output_dir = Path(config.output_dir)
        else:
            # Use the repository root (parent of the script's directory) as base
            script_dir = Path(__file__).resolve().parent
            repo_root = script_dir.parent.parent  # go up from tools/process_generator/
            base_output_dir = repo_root

        # Define directory structure according to CATChem layout
        process_dir = base_output_dir / config.src_base_dir / config.name
        test_dir = base_output_dir / "tests" / "process" / config.name
        docs_dir = base_output_dir / "docs" / "processes" / config.name

        # Create directory structure
        self._create_directory_structure(process_dir, test_dir, docs_dir, config)

        # Generate files
        self._generate_main_interface(process_dir, config)
        self._generate_common_module(process_dir, config)
        self._generate_creator_module(process_dir, config)
        self._generate_schemes(process_dir, config)
        self._generate_cmake_files(process_dir, config)

        if config.generate_tests:
            self._generate_tests(test_dir, config)

        if config.generate_docs:
            self._generate_documentation(docs_dir, config)

        if config.generate_examples:
            self._generate_examples(process_dir, config)

        logger.info(f"Process generation complete: {process_dir}")

    def _create_directory_structure(self, process_dir: Path, test_dir: Path, docs_dir: Path, config: ProcessConfig) -> None:
        """Create the directory structure for the process."""
        logger.info(f"Creating directory structure: {process_dir}")

        # Main process directories
        directories = [
            process_dir,
            process_dir / "schemes",
        ]

        # Test directories
        if config.generate_tests:
            directories.extend([
                test_dir,
                test_dir / "unit",
                test_dir / "integration"
            ])

        # Documentation directories
        if config.generate_docs:
            directories.append(docs_dir)

        # Examples in process directory
        if config.generate_examples:
            directories.append(process_dir / "examples")

        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)

    def _generate_main_interface(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate the main process interface module."""
        logger.info("Generating main process interface")

        template = self.env.get_template('process_interface.F90.j2')
        content = template.render(config=config, timestamp=datetime.now().isoformat())

        filename = f"Process{config.class_name}Interface_Mod.F90"
        output_file = process_dir / filename

        with open(output_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {output_file}")

    def _generate_common_module(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate the common types and utilities module."""
        logger.info("Generating common module")

        template = self.env.get_template('process_common.F90.j2')
        content = template.render(config=config, timestamp=datetime.now().isoformat())

        filename = f"{config.class_name}Common_Mod.F90"
        output_file = process_dir / filename

        with open(output_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {output_file}")

    def _generate_creator_module(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate the process creator module."""
        logger.info("Generating creator module")

        template = self.env.get_template('process_creator.F90.j2')
        content = template.render(config=config, timestamp=datetime.now().isoformat())

        filename = f"{config.class_name}ProcessCreator_Mod.F90"
        output_file = process_dir / filename

        with open(output_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {output_file}")

    def _generate_schemes(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate scheme implementation modules."""
        logger.info("Generating scheme modules")

        from dataclasses import asdict, is_dataclass
        def to_dict(obj):
            if isinstance(obj, dict):
                return obj
            elif is_dataclass(obj):
                return asdict(obj)
            else:
                return obj.__dict__

        # Restore to use the full scheme module template
        template = self.env.get_template('scheme_module.F90.j2')
        schemes_dir = process_dir / "schemes"

        for scheme in config.schemes:
            # Defensive: support both object and dict
            scheme_dict = to_dict(scheme)
            config_dict = to_dict(config)
            logger.info(f"Rendering scheme: {scheme_dict.get('name', '<unknown>')} with class {scheme_dict.get('class_name', '<unknown>')}")
            logger.info(f"SCHEME DICT: {scheme_dict}")
            logger.info(f"CONFIG DICT: {config_dict}")
            try:
                content = template.render(
                    config=config_dict,
                    scheme=scheme_dict,
                    timestamp=datetime.now().isoformat()
                )
                logger.info(f"Template rendered successfully, content length: {len(content)}")
            except Exception as e:
                logger.error(f"Template rendering failed for scheme {scheme_dict.get('name', '<unknown>')}: {e}")
                content = f"! Template rendering failed: {e}\n"

            filename = f"{config_dict['class_name']}Scheme_{scheme_dict['class_name']}_Mod.F90"
            output_file = schemes_dir / filename

            with open(output_file, 'w') as f:
                f.write(content)

            logger.info(f"Generated scheme: {output_file}")

    def _generate_cmake_files(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate CMake configuration files."""
        logger.info("Generating CMake files")

        # Main CMakeLists.txt
        template = self.env.get_template('CMakeLists.txt.j2')
        content = template.render(config=config, timestamp=datetime.now().isoformat())

        cmake_file = process_dir / "CMakeLists.txt"
        with open(cmake_file, 'w') as f:
            f.write(content)

        logger.info(f"Generated: {cmake_file}")

        # Schemes CMakeLists.txt
        schemes_cmake_template = self.env.get_template('schemes_CMakeLists.txt.j2')
        schemes_content = schemes_cmake_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        schemes_cmake_file = process_dir / "schemes" / "CMakeLists.txt"
        with open(schemes_cmake_file, 'w') as f:
            f.write(schemes_content)

        logger.info(f"Generated: {schemes_cmake_file}")

    def _generate_tests(self, test_dir: Path, config: ProcessConfig) -> None:
        """Generate test files in tests/process/<process_name> directory."""
        logger.info(f"Generating test files in: {test_dir}")

        # Unit tests
        unit_template = self.env.get_template('unit_test.F90.j2')
        unit_content = unit_template.render(config=config, timestamp=datetime.now().isoformat())

        unit_file = test_dir / "unit" / f"test_{config.name}_unit.F90"
        with open(unit_file, 'w') as f:
            f.write(unit_content)

        # Integration tests
        integration_template = self.env.get_template('integration_test.F90.j2')
        integration_content = integration_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        integration_file = test_dir / "integration" / f"test_{config.name}_integration.F90"
        with open(integration_file, 'w') as f:
            f.write(integration_content)

        logger.info(f"Generated tests in: {test_dir}")

        # Test CMakeLists.txt
        test_cmake_template = self.env.get_template('test_CMakeLists.txt.j2')
        test_cmake_content = test_cmake_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        test_cmake_file = test_dir / "CMakeLists.txt"
        with open(test_cmake_file, 'w') as f:
            f.write(test_cmake_content)

    def _generate_documentation(self, docs_dir: Path, config: ProcessConfig) -> None:
        """Generate consolidated documentation in docs/processes/<process_name> directory."""
        logger.info(f"Generating documentation in: {docs_dir}")

        # Single comprehensive documentation file
        doc_template = self.env.get_template('process_documentation.md.j2')
        doc_content = doc_template.render(config=config, timestamp=datetime.now().isoformat())

        doc_file = docs_dir / f"{config.name}.md"
        with open(doc_file, 'w') as f:
            f.write(doc_content)

        logger.info(f"Generated documentation in: {docs_dir}")

    def _generate_examples(self, process_dir: Path, config: ProcessConfig) -> None:
        """Generate example files."""
        logger.info("Generating examples")

        examples_dir = process_dir / "examples"

        # Basic usage example
        example_template = self.env.get_template('example_usage.F90.j2')
        example_content = example_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        example_file = examples_dir / f"{config.name}_example.F90"
        with open(example_file, 'w') as f:
            f.write(example_content)

        # Configuration example
        config_template = self.env.get_template('example_config.yaml.j2')
        config_content = config_template.render(
            config=config,
            timestamp=datetime.now().isoformat()
        )

        config_file = examples_dir / f"{config.name}_config.yaml"
        with open(config_file, 'w') as f:
            f.write(config_content)

        logger.info(f"Generated examples in: {examples_dir}")

    def generate_template_config(self, process_type: str = "emission") -> Dict[str, Any]:
        """Generate a template configuration for a given process type.

        Args:
            process_type: Type of process (emission, chemistry, transport, etc.)

        Returns:
            Dictionary containing template configuration
        """
        template_configs = {
            "emission": {
                "name": "my_emission",
                "class_name": "MyEmission",
                "description": "Description of my emission process",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "emission",
                "is_multiphase": False,
                "has_size_bins": False,
                "species": ["species1", "species2"],
                "schemes": [
                    {
                        "name": "simple",
                        "class_name": "Simple",
                        "description": "Simple emission scheme",
                        "author": "Your Name",
                        "required_met_fields": ["temperature", "wind_speed"],
                        "diagnostics": [
                            {
                                "name": "emission_rate",
                                "units": "kg/m2/s",
                                "description": "Emission rate"
                            }
                        ]
                    }
                ],
                "default_scheme": "simple",
                "required_met_fields": ["temperature", "wind_speed"],
                "diagnostics": [
                    {
                        "name": "total_emissions",
                        "units": "kg/m2/s",
                        "description": "Total emission rate"
                    }
                ],
                "generate_tests": True,
                "generate_docs": True,
                "generate_examples": True
            },

            "chemistry": {
                "name": "my_chemistry",
                "class_name": "MyChemistry",
                "description": "Description of my chemistry process",
                "author": "Your Name",
                "version": "1.0.0",
                "process_type": "chemistry",
                "is_multiphase": True,
                "has_size_bins": False,
                "species": ["O3", "NO", "NO2", "OH"],
                "phases": ["gas", "aqueous"],
                "schemes": [
                    {
                        "name": "mechanism1",
                        "class_name": "Mechanism1",
                        "description": "Chemistry mechanism 1",
                        "author": "Your Name",
                        "algorithm_type": "implicit",
                        "required_met_fields": ["temperature", "pressure", "humidity"],
                        "diagnostics": [
                            {
                                "name": "reaction_rates",
                                "units": "molec/cm3/s",
                                "description": "Chemical reaction rates"
                            }
                        ]
                    }
                ],
                "default_scheme": "mechanism1",
                "required_met_fields": ["temperature", "pressure", "humidity"],
                "timestep_dependency": "dependent",
                "parallelization": "column",
                "generate_tests": True,
                "generate_docs": True,
                "generate_examples": True
            }
        }

        return template_configs.get(process_type, template_configs["emission"])


def main():
    """Main entry point for the process generator."""
    parser = argparse.ArgumentParser(
        description="CATChem Process Generator Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate a process from configuration
  python process_generator.py generate --config my_process.yaml

  # Validate a configuration file
  python process_generator.py validate --config my_process.yaml

  # Generate template configuration
  python process_generator.py template --type emission --output emission_template.yaml

  # Generate with custom template directory
  python process_generator.py generate --config my_process.yaml --templates ./my_templates
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # Generate command
    generate_parser = subparsers.add_parser('generate', help='Generate process implementation')
    generate_parser.add_argument('--config', '-c', required=True,
                               help='Path to YAML configuration file')
    generate_parser.add_argument('--templates', '-t',
                               help='Path to template directory')
    generate_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Enable verbose output')

    # Validate command
    validate_parser = subparsers.add_parser('validate', help='Validate configuration file')
    validate_parser.add_argument('--config', '-c', required=True,
                               help='Path to YAML configuration file')
    validate_parser.add_argument('--verbose', '-v', action='store_true',
                               help='Enable verbose output')

    # Template command
    template_parser = subparsers.add_parser('template', help='Generate template configuration')
    template_parser.add_argument('--type', '-t', choices=['emission', 'chemistry', 'transport'],
                               default='emission', help='Type of process template')
    template_parser.add_argument('--output', '-o', required=True,
                               help='Output file for template')

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return 1

    # Set up logging level
    if hasattr(args, 'verbose') and args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        if args.command == 'generate':
            generator = ProcessGenerator(args.templates)
            config = generator.load_config(args.config)
            generator.generate_process(config)

        elif args.command == 'validate':
            generator = ProcessGenerator()
            config = generator.load_config(args.config)
            logger.info(f"Configuration is valid: {args.config}")

        elif args.command == 'template':
            generator = ProcessGenerator()
            template_config = generator.generate_template_config(args.type)

            with open(args.output, 'w') as f:
                yaml.dump(template_config, f, default_flow_style=False, sort_keys=False)

            logger.info(f"Template configuration written to: {args.output}")

    except ProcessValidationError as e:
        logger.error(f"Configuration validation failed: {e}")
        return 1
    except FileNotFoundError as e:
        logger.error(f"File not found: {e}")
        return 1
    except Exception as e:
        logger.error(f"Error: {e}")
        if hasattr(args, 'verbose') and args.verbose:
            import traceback
            traceback.print_exc()
        return 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
