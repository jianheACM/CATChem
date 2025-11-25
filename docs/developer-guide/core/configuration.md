# Configuration System Developer Guide

CATChem uses a flexible and robust YAML-based configuration system, managed by the `ConfigManager_Mod`. This guide provides a developer-oriented overview of the configuration system, explaining how to work with it to develop, modify, and extend CATChem.

## Overview

The configuration system is designed to be modern, type-safe, and easily extensible. Its key features include:

- **YAML-Based Configuration**: Human-readable and hierarchical configuration files.
- **Type-Safe Access**: A structured `ConfigDataType` in Fortran provides type-safe access to configuration parameters.
- **Schema Validation**: A `ConfigSchemaType` ensures that configuration files have the correct structure and required fields.
- **Modular Loading**: Specialized procedures for loading different parts of the configuration, such as species and emission data.
- **YAML Parsing in Fortran**: The system uses a Fortran-based YAML parser (`yaml_interface_mod`) for seamless integration.

## Architecture

The configuration system is primarily implemented in `src/core/ConfigManager_Mod.F90`.

### `ConfigManagerType`

The `ConfigManagerType` is the central component of the configuration system. It is responsible for loading, validating, and providing access to the configuration.

```fortran
! Located in: src/core/ConfigManager_Mod.F90

type :: ConfigManagerType
  private
  ! Configuration data
  type(yaml_node_t) :: yaml_data
  type(ConfigDataType), public :: config_data
  logical :: is_loaded = .false.

  ! Schema and validation
  type(ConfigSchemaType) :: schema
  logical :: schema_loaded = .false.
contains
  ! Initialization and cleanup
  procedure :: init => config_manager_init
  procedure :: finalize => config_manager_finalize

  ! Configuration loading
  procedure :: load_from_file => config_manager_load_from_file
  procedure :: load_from_string => config_manager_load_from_string

  ! Validation
  procedure :: validate => config_manager_validate

  ! Configuration access
  procedure :: get_string => config_manager_get_string
  procedure :: get_integer => config_manager_get_integer
  procedure :: get_real => config_manager_get_real
  procedure :: get_logical => config_manager_get_logical
  ! ... and more accessors
end type ConfigManagerType
```

### `ConfigDataType`

The `ConfigDataType` is a structured Fortran type that mirrors the structure of the main YAML configuration file. This allows for type-safe, compile-time checked access to configuration parameters.

```fortran
! Located in: src/core/ConfigManager_Mod.F90

type :: ConfigDataType
  ! Configuration categories
  type(RuntimeConfig) :: runtime
  type(FilePathConfig) :: file_paths
  type(ExternalEmisConfig) :: external_emissions
  type(EmissionMappingConfig) :: emission_mapping

  ! ... metadata and procedures ...
end type ConfigDataType
```

Each of the components within `ConfigDataType` (e.g., `RuntimeConfig`, `FilePathConfig`) is also a derived type, corresponding to a section in the YAML file.

## Configuration File Structure

The main configuration is typically provided in a YAML file (e.g., `CATChem_config.yml`). The structure of this file should match the structure of the `ConfigDataType`.

```yaml
# Example CATChem_config.yml

runtime:
  nLevs: 72
  nx: 48
  ny: 48
  nSpecies: 120
  maxSpecies: 200

file_paths:
  Species_File: "species.yml"
  Emission_File: "emissions.yml"
  Input_Directory: "./input"
  Output_Directory: "./output"

external_emissions:
  activate: true
  config_file: "external_emissions.yml"
  global_scale_factor: 1.0

# ... other sections ...
```

## Working with the Configuration System

### Loading a Configuration File

The `load_from_file` procedure is used to load a YAML configuration file. This procedure reads the file, parses the YAML, and populates the `ConfigDataType` structure.

```fortran
! Example of loading a configuration file
subroutine load_configuration(config_manager, filename, rc)
  use ConfigManager_Mod, only: ConfigManagerType
  implicit none

  type(ConfigManagerType), intent(inout) :: config_manager
  character(len=*), intent(in) :: filename
  integer, intent(out) :: rc

  ! Initialize the config manager
  call config_manager%init(rc)
  if (rc /= 0) return

  ! Load the configuration from the file
  call config_manager%load_from_file(filename, rc)
  if (rc /= 0) then
    print *, "Error: Failed to load configuration file: ", trim(filename)
    return
  end if

  ! The configuration is now loaded and accessible
end subroutine load_configuration
```

### Accessing Configuration Parameters

Once the configuration is loaded, you can access parameters in a type-safe manner through the `config_data` member of the `ConfigManagerType` or by using the getter methods.

```fortran
! Example of accessing configuration parameters
subroutine use_configuration(config_manager)
  use ConfigManager_Mod, only: ConfigManagerType
  implicit none

  type(ConfigManagerType), intent(in) :: config_manager
  integer :: n_levels, n_species
  character(len=256) :: species_file
  integer :: rc

  ! Direct access through the structured data type
  n_levels = config_manager%config_data%runtime%nLevs
  species_file = config_manager%config_data%file_paths%Species_File

  print *, "Number of levels: ", n_levels
  print *, "Species file: ", trim(species_file)

  ! Access using getter methods
  call config_manager%get_integer('runtime/nSpecies', n_species, rc)
  if (rc == 0) then
    print *, "Number of species: ", n_species
  end if

end subroutine use_configuration
```

### Schema Validation

The configuration system includes a schema validation mechanism to ensure that the loaded configuration is valid. The `ConfigManagerType` has a default schema, and you can also load a custom schema from a file. The `validate` procedure checks the loaded configuration against the schema.

```fortran
! Example of validating the configuration
subroutine validate_configuration(config_manager, rc)
  use ConfigManager_Mod, only: ConfigManagerType
  implicit none

  type(ConfigManagerType), intent(inout) :: config_manager
  integer, intent(out) :: rc

  ! Validate the loaded configuration against the schema
  call config_manager%validate(rc)
  if (rc /= 0) then
    print *, "Error: Configuration validation failed."
  else
    print *, "Configuration is valid."
  end if

end subroutine validate_configuration
```

### Specialized Loading Procedures

The `ConfigManager_Mod` provides specialized procedures for loading critical parts of the model configuration, such as species and emission data.

#### Loading Species Configuration

The `load_and_init_species` procedure is designed to load a detailed species configuration file and initialize the `ChemStateType` with the loaded data.

```fortran
! Example of loading species data
subroutine load_species(config_manager, chem_state, error_mgr, grid, rc)
  use ConfigManager_Mod, only: ConfigManagerType
  use ChemState_Mod, only: ChemStateType
  use Error_Mod, only: ErrorManagerType
  use GridGeometry_Mod, only: GridGeometryType
  implicit none

  type(ConfigManagerType), intent(inout) :: config_manager
  type(ChemStateType), intent(inout) :: chem_state
  type(ErrorManagerType), pointer, intent(inout) :: error_mgr
  type(GridGeometryType), pointer, intent(in) :: grid
  integer, intent(out) :: rc
  character(len=256) :: species_file

  species_file = config_manager%get_species_file()

  call config_manager%load_and_init_species(species_file, chem_state, error_mgr, grid, rc)
  if (rc /= 0) then
    print *, "Error: Failed to load species configuration."
  end if

end subroutine load_species
```

This procedure reads a YAML file containing a list of species and their properties (molecular weight, type flags, etc.), and populates the `ChemStateType` accordingly.

This powerful and flexible configuration system is a cornerstone of CATChem's modular design, enabling complex model configurations while maintaining clarity and robustness.