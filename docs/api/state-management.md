# State Management API

This section covers the core state management APIs in CATChem, providing unified data handling across all atmospheric processes.

## Overview

The state management system provides:
- **StateContainer**: Central data repository for all model state
- **ChemState**: Chemical species concentrations and properties
- **MetState**: Meteorological fields and atmospheric conditions
- **EmisState**: Emission data and accumulation
- **DiagState**: Diagnostic variables and output management

## Core Components

### StateContainer

The `StateContainer` is the central hub for all model data:

```fortran
use State_Mod
type(StateContainerType) :: container

! Initialize with configuration
call container%init(config, rc)

! Access chemical state
call container%get_chem_state(chem_state, rc)

! Access meteorological state
call container%get_met_state(met_state, rc)
```

**Key Features:**
- Unified data access across all processes
- Automatic memory management
- Thread-safe operations
- Efficient column-based access

**Auto-Generated Documentation:** [State Module](../CATChem/namespacestate__mod.md)

### Column Virtualization

Access data efficiently through column interfaces:

```fortran
use ColumnInterface_Mod
type(ColumnType) :: column

! Get column data
call container%get_column(i, j, column, rc)

! Process column data
do k = 1, column%nz
    ! Process atmospheric level k
    column%chem_data(k, species_idx) = new_value
end do

! Update container
call container%update_column(i, j, column, rc)
```

**Performance Benefits:**
- Cache-optimized 1D processing
- Reduced memory allocations
- Natural parallelization

## State Types

### Chemical State (ChemState)

Manages chemical species concentrations:

```fortran
type(ChemStateType) :: chem_state

! Get species concentration
call chem_state%get_species_conc(species_name, concentration, rc)

! Set species concentration
call chem_state%set_species_conc(species_name, new_concentration, rc)

! Get all species data
call chem_state%get_all_species(species_data, rc)
```

### Meteorological State (MetState)

Handles atmospheric conditions:

```fortran
type(MetStateType) :: met_state

! Get temperature profile
call met_state%get_temperature(temperature, rc)

! Get pressure levels
call met_state%get_pressure(pressure, rc)

! Get wind components
call met_state%get_wind(u_wind, v_wind, w_wind, rc)
```

### Emission State (EmisState)

Tracks emission sources and accumulation:

```fortran
type(EmisStateType) :: emis_state

! Set emission rates
call emis_state%set_emission_rates(species_name, rates, rc)

! Accumulate emissions for diagnostics
call emis_state%accumulate_emissions(species_name, increment, rc)

! Get total emissions
call emis_state%get_total_emissions(species_name, total, rc)
```

## Data Access Patterns

### Direct Access

For simple operations:

```fortran
! Get pointer to concentration array
real(fp), pointer :: conc(:,:,:) => null()
call container%get_concentration_ptr(species_idx, conc, rc)

! Modify data directly
conc(i, j, k) = new_value
```

### Column Processing

For optimal performance:

```fortran
! Process all columns
do j = 1, ny
    do i = 1, nx
        call container%get_column(i, j, column, rc)
        call process_column(column)
        call container%update_column(i, j, column, rc)
    end do
end do
```

### Batch Operations

For bulk data handling:

```fortran
! Get multiple species at once
call container%get_species_batch(species_list, concentration_batch, rc)

! Update multiple species
call container%set_species_batch(species_list, new_concentrations, rc)
```

## Error Handling

All state operations include comprehensive error handling:

```fortran
use Error_Mod

! Check return codes
call container%get_species_conc(species_name, conc, rc)
if (rc /= CC_SUCCESS) then
    call error_mgr%report_error(ERROR_DATA_ACCESS, &
                               'Failed to get species concentration', rc)
    return
endif

! Use error context for debugging
call error_mgr%push_context('state_operation', 'Getting chemical data')
call container%operation(data, rc)
call error_mgr%pop_context()
```

## Thread Safety

State operations are designed for parallel processing:

```fortran
! Thread-safe column access
!$OMP PARALLEL DO PRIVATE(column, rc)
do j = 1, ny
    do i = 1, nx
        call container%get_column_safe(i, j, column, rc)
        call process_column(column)
        call container%update_column_safe(i, j, column, rc)
    end do
end do
!$OMP END PARALLEL DO
```

## Memory Management

Automatic memory management with manual control when needed:

```fortran
! Automatic cleanup
call container%finalize(rc)  ! Cleans up all allocated memory

! Manual memory control
call container%allocate_workspace(size, rc)
call container%deallocate_workspace(rc)

! Memory optimization
call container%optimize_memory_layout(rc)
```

## Best Practices

### Performance

1. **Use column processing** for atmospheric calculations
2. **Batch operations** for multiple species
3. **Avoid unnecessary copies** - use pointers when possible
4. **Cache column data** for repeated access

### Safety

1. **Always check return codes**
2. **Use error context** for debugging
3. **Initialize before use**
4. **Clean up resources** in finalize routines

### Debugging

1. **Enable debug mode** for detailed logging
2. **Use error context** to track operation flow
3. **Validate data bounds** in development builds
4. **Monitor memory usage** for large simulations

## See Also

- [Process Interface API](process-interface.md) - How processes interact with state
- [Column Interface API](column-interface.md) - Efficient column-based processing
- [Configuration API](configuration.md) - State initialization and setup
- [Error Handling API](error-handling.md) - Comprehensive error management

---

**Auto-Generated Documentation:** [Complete State Management Reference](../CATChem/namespacestate__mod.md)
