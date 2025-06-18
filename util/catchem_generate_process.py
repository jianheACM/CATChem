#!/usr/bin/env python3
"""
CATChem Process Generator Tool

This tool generates boilerplate code for new atmospheric processes following
the established architecture patterns. It creates process modules, test files,
and updates build systems automatically.

Usage:
    python catchem_generate_process.py --name=NewProcess --type=emission --schemes=scheme1,scheme2
    python catchem_generate_process.py --name=Chemistry --type=multiphase_chemistry --solver=micm --phases=gas,liquid
    python catchem_generate_process.py --config=process_config.yaml

Author: CATChem Development Team
Date: 2025
"""

import argparse
import os
import sys
import logging
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Optional, Union
import json
import yaml
from jinja2 import Environment, FileSystemLoader, select_autoescape

class ProcessGenerator:
    """Main process generator class using Jinja2 templates"""

    def __init__(self, base_path: str = "."):
        self.base_path = Path(base_path)
        self.templates_dir = self.base_path / "util" / "templates"

        # Set up logging
        self.logger = logging.getLogger(self.__class__.__name__)
        self.templates_used = []  # Track templates used during generation

        # Initialize Jinja2 environment
        self.jinja_env = Environment(
            loader=FileSystemLoader(str(self.templates_dir)),
            autoescape=select_autoescape(['html', 'xml']),
            trim_blocks=True,
            lstrip_blocks=True
        )

        self.process_types = {
            'emission': {
                'description': 'Emission process (sources)',
                'common_methods': ['validate_emission_inputs', 'zero_emission_arrays',
                                 'accumulate_emissions', 'apply_emission_scaling'],
                'state_dependencies': ['met_state', 'emis_state'],
                'template': 'process_main.f90.j2'
            },
            'transformation': {
                'description': 'Transformation process (chemistry, physics)',
                'common_methods': ['validate_chem_inputs', 'update_chemistry',
                                 'check_mass_conservation'],
                'state_dependencies': ['chem_state', 'met_state'],
                'template': 'process_main.f90.j2'
            },
            'loss': {
                'description': 'Loss process (deposition, scavenging, decay)',
                'common_methods': ['validate_loss_inputs', 'calculate_loss_rates',
                                 'apply_loss_processes'],
                'state_dependencies': ['chem_state', 'met_state'],
                'template': 'process_main.f90.j2'
            },
            'sink': {
                'description': 'Sink process (deposition, scavenging)',
                'common_methods': ['validate_sink_inputs', 'apply_sink_processes',
                                 'update_deposition_state'],
                'state_dependencies': ['chem_state', 'met_state', 'dep_state'],
                'template': 'process_main.f90.j2'
            },
            'transport': {
                'description': 'Transport process (advection, diffusion)',
                'common_methods': ['validate_transport_inputs', 'apply_transport',
                                 'update_transport_state'],
                'state_dependencies': ['chem_state', 'met_state'],
                'template': 'process_main.f90.j2'
            },
            'multiphase_chemistry': {
                'description': 'Multiphase chemistry process (gas/liquid/solid reactions)',
                'common_methods': ['validate_multiphase_inputs', 'setup_phase_equilibrium',
                                 'calculate_henry_constants', 'update_photolysis_rates',
                                 'solve_chemistry_system', 'partition_species', 'check_phase_balance'],
                'state_dependencies': ['chem_state', 'met_state'],
                'template': 'process_multiphase.f90.j2',
                'requires_multiphase': True,
                'solver_types': ['euler', 'rk4', 'rosenbrock', 'kpp', 'micm'],
                'phase_types': ['gas', 'liquid', 'solid']
            },
            'aqueous_chemistry': {
                'description': 'Aqueous chemistry process (cloud/fog chemistry)',
                'common_methods': ['validate_multiphase_inputs', 'calculate_henry_constants',
                                 'setup_phase_equilibrium', 'solve_chemistry_system',
                                 'partition_species', 'check_phase_balance'],
                'state_dependencies': ['chem_state', 'met_state'],
                'template': 'process_multiphase.f90.j2',
                'requires_multiphase': True,
                'phase_types': ['gas', 'liquid']
            },
            'heterogeneous_chemistry': {
                'description': 'Heterogeneous chemistry process (gas-aerosol reactions)',
                'common_methods': ['validate_multiphase_inputs', 'calculate_uptake_coefficients',
                                 'solve_chemistry_system', 'update_aerosol_surface_area'],
                'state_dependencies': ['chem_state', 'met_state'],
                'template': 'process_multiphase.f90.j2',
                'requires_multiphase': True,
                'phase_types': ['gas', 'solid']
            }
        }

    def load_config_from_yaml(self, config_file: Union[str, Path]) -> Dict:
        """Load process configuration from YAML file"""
        config_path = Path(config_file)
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_file}")

        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        return config

    def generate_process_from_config(self, config: Dict) -> bool:
        """Generate a process from configuration dictionary"""
        return self.generate_process(
            name=config['name'],
            process_type=config['type'],
            schemes=config.get('schemes'),
            multi_phase=config.get('multi_phase', False),
            author=config.get('author', "CATChem Development Team"),
            solver_type=config.get('solver_type', 'euler'),
            phases=config.get('phases'),
            species=config.get('species', []),
            description=config.get('description', '')
        )

    def generate_process(self, name: str, process_type: str, schemes: Optional[List[str]] = None,
                        multi_phase: bool = False, validate_only: bool = False,
                        author: str = "CATChem Development Team", solver_type: str = "euler",
                        phases: Optional[List[str]] = None, species: Optional[List[str]] = None,
                        description: str = "") -> bool:
        """Generate a new process with all required files"""

        # Reset template tracking
        self.templates_used = []

        if validate_only:
            return self._validate_process_config(name, process_type, schemes)

        # Auto-enable multiphase for multiphase process types
        if process_type in ['multiphase_chemistry', 'aqueous_chemistry', 'heterogeneous_chemistry']:
            multi_phase = True
            if phases is None:
                phases = self.process_types[process_type].get('phase_types', ['gas'])

        self.logger.info(f"Starting generation of {name} process (type: {process_type})")
        self.logger.info(f"Author: {author}")
        self.logger.info(f"Schemes: {schemes or 'None'}")
        self.logger.info(f"Multiphase: {multi_phase}")
        if multi_phase:
            self.logger.info(f"Phases: {phases}")
            self.logger.info(f"Solver: {solver_type}")

        try:
            # Create directory structure
            self.logger.info("Creating directory structure...")
            self._create_directory_structure(name, process_type, multi_phase)

            # Generate main process module
            self.logger.info("Generating main process module...")
            self._generate_main_process(name, process_type, schemes, multi_phase, author,
                                      solver_type, phases, species, description)

            # Generate common module if needed
            if schemes:
                self.logger.info("Generating common module...")
                self._generate_common_module(name, schemes, author)

            # Generate scheme modules
            if schemes:
                self.logger.info(f"Generating scheme modules for: {', '.join(schemes)}")
                self._generate_scheme_modules(name, process_type, schemes, author)

            # Generate test files
            self.logger.info("Generating test files...")
            self._generate_test_files(name, process_type, schemes, multi_phase)

            # Update build system
            self.logger.info("Updating build system...")
            self._update_build_system(name, process_type, multi_phase)

            # Generate documentation
            self.logger.info("Generating documentation...")
            self._generate_documentation(name, process_type, schemes)

            # Log template usage summary
            self.logger.info("Template usage summary:")
            for template_info in self.templates_used:
                self.logger.info(f"  {template_info['file']} -> {template_info['template']}")

            self.logger.info(f"Successfully generated {name} process!")
            self.logger.info(f"Files created in: src/process/{name.lower()}/")
            self.logger.info(f"Tests created in: tests/process/{name.lower()}/")

            return True

        except Exception as e:
            self.logger.error(f"Error generating process: {e}")
            return False

    def _create_directory_structure(self, name: str, process_type: str, multi_phase: bool = False):
        """Create the directory structure for the new process"""
        base_dir = self.base_path / "src" / "process" / name.lower()
        base_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        (base_dir / "schemes").mkdir(exist_ok=True)

        # Create multiphase-specific directories
        if multi_phase:
            (base_dir / "phases").mkdir(exist_ok=True)
            (base_dir / "solvers").mkdir(exist_ok=True)
            (base_dir / "equilibrium").mkdir(exist_ok=True)

        # Create test directory
        test_dir = self.base_path / "tests" / "process" / name.lower()
        test_dir.mkdir(parents=True, exist_ok=True)

    def _generate_main_process(self, name: str, process_type: str, schemes: Optional[List[str]],
                             multi_phase: bool, author: str, solver_type: str = "euler",
                             phases: Optional[List[str]] = None, species: Optional[List[str]] = None,
                             description: str = ""):
        """Generate the main process module using Jinja2 templates"""

        template_vars = {
            'process_name': name,
            'process_name_lower': name.lower(),
            'process_name_upper': name.upper(),
            'process_type': process_type,
            'process_description': description or self.process_types[process_type]['description'],
            'author': author,
            'year': datetime.now().year,
            'schemes': schemes or [],
            'multi_phase': multi_phase,
            'solver_type': solver_type,
            'phases': phases or ['gas'],
            'n_phases': len(phases) if phases else 1,
            'has_aqueous': 'liquid' in (phases or []),
            'has_heterogeneous': 'solid' in (phases or []),
            'species': species or [],
            'common_methods': self.process_types[process_type]['common_methods'],
            'state_dependencies': self.process_types[process_type]['state_dependencies']
        }

        # Load and render template
        template_name = self.process_types[process_type]['template']
        self.logger.debug(f"Using template: {template_name}")
        template = self.jinja_env.get_template(template_name)
        process_content = template.render(**template_vars)

        # Write main process file
        output_file = (self.base_path / "src" / "process" / name.lower() /
                      f"{name}Process_Mod.F90")

        # Track template usage
        self.templates_used.append({
            'file': f"{name}Process_Mod.F90",
            'template': template_name,
            'path': str(output_file.relative_to(self.base_path))
        })

        self.logger.debug(f"Generated: {output_file}")

        with open(output_file, 'w') as f:
            f.write(process_content)

    def _generate_common_module(self, name: str, schemes: Optional[List[str]], author: str):
        """Generate the common utility module using Jinja2 templates"""

        template_vars = {
            'process_name': name,
            'process_name_lower': name.lower(),
            'process_name_upper': name.upper(),
            'author': author,
            'year': datetime.now().year,
            'schemes': schemes or []
        }

        template = self.jinja_env.get_template('common_module.f90.j2')
        self.logger.debug(f"Using template: common_module.f90.j2")
        common_content = template.render(**template_vars)

        output_file = (self.base_path / "src" / "process" / name.lower() /
                      f"{name}Common_Mod.F90")

        with open(output_file, 'w') as f:
            f.write(common_content)

        self.logger.info(f"Created: {output_file}")

    def _generate_scheme_modules(self, name: str, process_type: str, schemes: List[str], author: str):
        """Generate scheme modules using Jinja2 templates"""

        schemes_dir = self.base_path / "src" / "process" / name.lower() / "schemes"

        for scheme in schemes:
            template_vars = {
                'scheme_name': scheme,
                'process_name': name,
                'process_type': process_type,
                'author': author,
                'year': datetime.now().year,
                'has_init_routine': False  # Keep schemes simple
            }

            template = self.jinja_env.get_template('scheme_simple.f90.j2')
            self.logger.debug(f"Using template: scheme_simple.f90.j2 for {scheme}")
            scheme_content = template.render(**template_vars)

            output_file = schemes_dir / f"{scheme}Scheme_Mod.F90"

            with open(output_file, 'w') as f:
                f.write(scheme_content)

            self.logger.info(f"Created: {output_file}")

    def _generate_test_files(self, name: str, process_type: str,
                           schemes: Optional[List[str]], multi_phase: bool = False):
        """Generate unit test files using Jinja2 templates"""

        test_dir = self.base_path / "tests" / "process" / name.lower()

        template_vars = {
            'process_name': name,
            'process_name_lower': name.lower(),
            'process_type': process_type,
            'schemes': schemes or [],
            'multi_phase': multi_phase
        }

        # Main process test
        template = self.jinja_env.get_template('test_process.f90.j2')
        self.logger.debug(f"Using template: test_process.f90.j2")
        test_content = template.render(**template_vars)

        test_file = test_dir / f"test_{name.lower()}_process.F90"
        with open(test_file, 'w') as f:
            f.write(test_content)
        self.logger.info(f"Created: {test_file}")

    def _update_build_system(self, name: str, process_type: str, multi_phase: bool = False):
        """Update CMakeLists.txt files"""

        # Create process-specific CMakeLists.txt
        process_cmake = self.base_path / "src" / "process" / name.lower() / "CMakeLists.txt"

        template_vars = {
            'process_name': name,
            'process_name_lower': name.lower(),
            'multi_phase': multi_phase
        }

        template = self.jinja_env.get_template('cmake_process.txt.j2')
        self.logger.debug(f"Using template: cmake_process.txt.j2")
        cmake_content = template.render(**template_vars)

        cmake_file = process_cmake
        with open(process_cmake, 'w') as f:
            f.write(cmake_content)
        self.logger.info(f"Created: {cmake_file}")

    def _generate_documentation(self, name: str, process_type: str, schemes: Optional[List[str]]):
        """Generate documentation files"""

        doc_dir = self.base_path / "docs" / "processes"
        doc_dir.mkdir(parents=True, exist_ok=True)

        template_vars = {
            'process_name': name,
            'process_type': process_type,
            'process_description': self.process_types[process_type]['description'],
            'schemes': schemes or []
        }

        template = self.jinja_env.get_template('process_documentation.md.j2')
        self.logger.debug(f"Using template: process_documentation.md.j2")
        doc_content = template.render(**template_vars)

        doc_file = doc_dir / f"{name.lower()}_process.md"
        with open(doc_file, 'w') as f:
            f.write(doc_content)
        self.logger.info(f"Created: {doc_file}")

    def _validate_process_config(self, name: str, process_type: str, schemes: Optional[List[str]]) -> bool:
        """Validate process configuration"""

        if not name:
            self.logger.error("Process name is required")
            return False

        if process_type not in self.process_types:
            self.logger.error(f"Invalid process type: {process_type}")
            self.logger.error(f"Available types: {list(self.process_types.keys())}")
            return False

        if schemes:
            for scheme in schemes:
                if not scheme.replace('_', '').isalnum():
                    self.logger.error(f"Invalid scheme name: {scheme}")
                    return False

        self.logger.info("Process configuration is valid")
        return True

def main():
    parser = argparse.ArgumentParser(
        description="Generate CATChem atmospheric process boilerplate code",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate dust emission process with multiple schemes
  python catchem_generate_process.py --name=Dust --type=emission --schemes=fengsha,ginoux

  # Generate multiphase chemistry process
  python catchem_generate_process.py --name=Chemistry --type=multiphase_chemistry --solver=micm --phases=gas,liquid

  # Generate from YAML configuration
  python catchem_generate_process.py --config=process_config.yaml

  # Validate configuration only
  python catchem_generate_process.py --name=NewProcess --type=loss --validate-only

Process Types:
  emission       - Emission processes (sources)
  transformation - Transformation processes (chemistry, physics)
  loss           - Loss processes (deposition, scavenging, decay)
  sink           - Sink processes (deposition, scavenging)
  transport      - Transport processes (advection, diffusion)
  multiphase_chemistry    - Multiphase chemistry (gas/liquid/solid)
  aqueous_chemistry       - Aqueous chemistry (cloud/fog)
  heterogeneous_chemistry - Heterogeneous chemistry (gas-aerosol)
        """)

    parser.add_argument('--config',
                       help='YAML configuration file')
    parser.add_argument('--name',
                       help='Process name (e.g., Dust, SeaSalt, Chemistry)')
    parser.add_argument('--type',
                       choices=['emission', 'transformation', 'loss', 'sink', 'transport',
                               'multiphase_chemistry', 'aqueous_chemistry', 'heterogeneous_chemistry'],
                       help='Process type')
    parser.add_argument('--schemes',
                       help='Comma-separated list of scheme names (e.g., fengsha,ginoux)')
    parser.add_argument('--multi-phase', action='store_true',
                       help='Enable multi-phase process support')
    parser.add_argument('--solver', default='euler',
                       choices=['euler', 'rk4', 'rosenbrock', 'kpp', 'micm'],
                       help='Chemistry solver type for multiphase processes')
    parser.add_argument('--phases',
                       help='Comma-separated list of phases (gas,liquid,solid)')
    parser.add_argument('--species',
                       help='Comma-separated list of species names')
    parser.add_argument('--validate-only', action='store_true',
                       help='Only validate configuration, do not generate files')
    parser.add_argument('--author', default='CATChem Development Team',
                       help='Author name for generated files')
    parser.add_argument('--base-path', default='.',
                       help='Base path for CATChem source code')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Enable verbose logging')

    args = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Create generator
    generator = ProcessGenerator(args.base_path)

    # Generate from config file or command line arguments
    if args.config:
        try:
            config = generator.load_config_from_yaml(args.config)
            success = generator.generate_process_from_config(config)
        except Exception as e:
            logging.error(f"Error loading config: {e}")
            sys.exit(1)
    else:
        # Parse command line arguments
        if not args.name or not args.type:
            logging.error("--name and --type are required when not using --config")
            parser.print_help()
            sys.exit(1)

        schemes = []
        if args.schemes:
            schemes = [s.strip().title() for s in args.schemes.split(',')]

        phases = []
        if args.phases:
            phases = [p.strip().lower() for p in args.phases.split(',')]

        species = []
        if args.species:
            species = [s.strip() for s in args.species.split(',')]

        success = generator.generate_process(
            name=args.name,
            process_type=args.type,
            schemes=schemes,
            multi_phase=args.multi_phase,
            validate_only=args.validate_only,
            author=args.author,
            solver_type=args.solver,
            phases=phases,
            species=species
        )

    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
