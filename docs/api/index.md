# API Reference

Welcome to the CATChem API documentation. This section provides comprehensive reference documentation for all CATChem modules, types, and procedures.

!!! note "Auto-Generated Documentation"
    This API documentation is automatically generated from the source code using [MkDoxy](https://mkdoxy.kubaandrysek.cz/). The documentation is updated with each code change to ensure accuracy.

## 🔗 Browse Complete API Documentation

**[→ Full Auto-Generated API Documentation](../CATChem/)**

The complete API documentation includes all modules, with cross-references, source code links, and detailed documentation for every procedure and type.

## 🎯 Key Modules Quick Access

### Core System
- **[State Management](../CATChem/state__mod_8_f90/)** - Central state container and data management
- **[Column Interface](../CATChem/_column_interface___mod_8_f90/)** - Column virtualization system
- **[Process Manager](../CATChem/_process_manager___mod_8_f90/)** - Process orchestration and lifecycle
- **[Process Interface](../CATChem/_process_interface___mod_8_f90/)** - Base interface for all processes

### Configuration System
- **[Config Manager](../CATChem/_config_manager___mod_8_f90/)** - Configuration file parsing and management
- **[Field Mapping](../CATChem/_field_mapping___mod_8_f90/)** - Input/output field mapping utilities

### Process Implementations

#### Transport Processes
- **[Settling Process](../CATChem/settling_process___mod_8_f90/)** - Gravitational settling with slip correction
- **[Stokes Scheme](../CATChem/_stokesscheme_scheme___mod_8_f90/)** - Modern Stokes settling implementation
- **[YSU Vertical Dispersion](../CATChem/_y_s_u_vertical_dispersion_process___mod_8_f90/)** - Boundary layer mixing

#### Emission Processes
- **[External Emissions](../CATChem/_external_emission_process___mod_8_f90/)** - Anthropogenic and biogenic emissions
- **[Dust Process](../CATChem/_dust_process___mod_8_f90/)** - Mineral dust emission and transport
- **[Sea Salt Process](../CATChem/_sea_salt_process___mod_8_f90/)** - Marine aerosol processes

#### Loss Processes
- **[Dry Deposition](../CATChem/_dry_dep_process___mod_8_f90/)** - Surface deposition processes
- **[Wesely Scheme](../CATChem/_wesely_scheme___mod_8_f90/)** - Dry deposition parametrization

### Diagnostic System
- **[Diagnostic Manager](../CATChem/_diagnostic_manager___mod_8_f90/)** - Diagnostic output management
- **[Diagnostic Interface](../CATChem/_diagnostic_interface___mod_8_f90/)** - Diagnostic data interfaces

### Utilities
- **[Constants](../CATChem/constants_8_f90/)** - Physical and mathematical constants
- **[Utilities](../CATChem/utilities__mod_8_f90/)** - Common utility functions
- **[Error Handling](../CATChem/error__mod_8_f90/)** - Error management system

## Quick Reference

### Key Types

| Type | Module | Description |
|------|--------|-------------|
| `CATChemType` | `CATChemAPI_Mod` | Main API interface |
| `StateContainerType` | `state_mod` | Central data container |
| `ProcessInterface` | `ProcessInterface_Mod` | Base class for processes |
| `ErrorManagerType` | `error_mod` | Error handling |
| `ConfigDataType` | `config_mod` | Configuration management |

### Common Patterns

=== "Process Initialization"

    ```fortran
    use ProcessName_Mod
    type(ProcessNameType) :: process
    call process%init(container, rc)
    ```

=== "Error Handling"

    ```fortran
    use error_mod
    type(ErrorManagerType), pointer :: error_mgr
    error_mgr => container%get_error_manager()
    call error_mgr%push_context('routine_name', 'description')
    ! ... operations ...
    call error_mgr%pop_context()
    ```

=== "Diagnostic Access"

    ```fortran
    use DiagnosticInterface_Mod
    type(DiagnosticFieldType), pointer :: field
    field => diag_mgr%get_field('field_name', rc)
    call field%get_data(data_array, rc)
    ```

## Search Tips

- Use the search box above to find specific procedures or types
- Browse by module for related functionality
- Check the inheritance hierarchy for process types
- Look at usage examples in the source code

## Conventions

### Naming Conventions
- **Modules**: `ModuleName_Mod`
- **Types**: `TypeNameType`
- **Procedures**: `snake_case`
- **Constants**: `UPPER_CASE`

### Return Codes
All procedures use integer return codes following the convention:
- `CC_SUCCESS = 0` - Successful operation
- `CC_FAILURE = -1` - Generic failure
- Specific error codes defined in `error_mod`

### Memory Management
- Pointer associations managed by StateContainer
- No explicit allocation in process modules
- Use `intent(in)` for immutable data, `intent(inout)` for modifications

## Contributing

Found an issue with the documentation? The API docs are generated automatically, so:

1. **Source Code Issues**: Update the source code comments and docstrings
2. **Organization Issues**: Modify the MkDoxy configuration in `mkdocs.yml`
3. **Missing Documentation**: Add Doxygen-style comments to the source code

For details on documentation standards, see the developer guide section on documentation.
