# API Reference

Welcome to the CATChem API documentation. This section provides comprehensive reference documentation for all CATChem modules, types, and procedures.

!!! note "Documentation Structure"
    Our API documentation combines hand-written guides with auto-generated reference material from the source code using [MkDoxy](https://mkdoxy.kubaandrysek.cz/). This ensures both comprehensive coverage and up-to-date accuracy.

## � API Organization

### Core Modules
High-level APIs for the main CATChem systems:

- **[State Management](state-management.md)** - StateContainer, ChemState, MetState data handling
- **[Process Interface](process-interface.md)** - Process development and integration APIs
- **[Column Interface](column-interface.md)** - Column virtualization and 1D processing
- **[Configuration Manager](configuration.md)** - YAML configuration system
- **[Error Handling](error-handling.md)** - Comprehensive error management

### Process Modules
APIs for specific atmospheric process implementations:

- **[Emission Processes](emission-processes.md)** - External emissions, dust, sea salt
- **[Transport Processes](transport-processes.md)** - Settling, vertical mixing, diffusion
- **[Chemistry Processes](chemistry-processes.md)** - Gas-phase and aerosol chemistry
- **[Deposition Processes](deposition-processes.md)** - Dry and wet deposition

### Utility Modules
Supporting APIs and utilities:

- **[Constants](constants.md)** - Physical and mathematical constants
- **[Utilities](utilities.md)** - Common utility functions and tools
- **[Diagnostics](diagnostics.md)** - Diagnostic output and monitoring

## 🔗 Auto-Generated Documentation

**[→ Complete Auto-Generated API Reference](../CATChem/)**

The complete API documentation includes:
- **[Modules Index](../CATChem/modules.md)** - All Fortran modules
- **[Namespaces Index](../CATChem/namespaces.md)** - Module namespaces
- **[Files Index](../CATChem/files.md)** - Source file documentation

### Key Auto-Generated Sections

**Core System Modules:**
- **[State Module](../CATChem/namespacestate__mod.md)** - StateContainer implementation
- **[Process Manager](../CATChem/namespaceprocess__manager__mod.md)** - Process orchestration
- **[Column Interface](../CATChem/namespacecolumn__interface__mod.md)** - Column virtualization
- **[Config Manager](../CATChem/namespaceconfig__manager__mod.md)** - Configuration system

**Process Implementation Modules:**
- **[Emission Processes](../CATChem/namespaceemission__process__mod.md)** - Emission process implementations
- **[Transport Processes](../CATChem/namespacetransport__process__mod.md)** - Transport physics
- **[Chemistry Processes](../CATChem/namespacechemistry__process__mod.md)** - Chemical mechanisms

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
