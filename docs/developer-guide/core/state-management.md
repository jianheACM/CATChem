# State Management Developer Guide

CATChem's state management system provides a modern, flexible architecture for handling model state data. This guide describes the state management system, focusing on how to develop, modify, and extend CATChem's state-handling capabilities.

## Overview

The state management system is built around the `StateManagerType`, which acts as a central container for all model state. This approach avoids global variables and provides a structured, maintainable way to manage data. Key features include:

- **`StateManagerType`**: A central container for all model state components.
- **Componentized State**: State is divided into components like `MetStateType` for meteorology and `ChemStateType` for chemistry.
- **Type-Bound Procedures**: State components are managed through modern Fortran type-bound procedures for initialization, cleanup, and validation.
- **Pointer-Based Access**: Efficient, direct access to data is provided via pointers.

## StateManager Architecture

### `StateManagerType`

The `StateManagerType`, defined in `StateManager_Mod.F90`, is the core of the state management system. It encapsulates all other state components.

```fortran
! Located in: src/core/StateManager_Mod.F90

type :: StateManagerType
  private

  ! Core state objects
  type(MetStateType),   allocatable :: met_state   !< Meteorological fields
  type(ChemStateType),  allocatable :: chem_state  !< Chemical species concentrations
  type(ErrorManagerType)            :: error_mgr   !< Error manager

  ! Manager pointers (owned by CATChemCore)
  type(ConfigManagerType), pointer :: config => null()  !< Configuration manager
  type(GridManagerType), pointer :: grid_mgr => null()  !< Grid manager
  type(DiagnosticManagerType), pointer :: diag_mgr => null()  !< Diagnostic manager

  ! ... metadata and procedures ...
contains
  ! Basic lifecycle (called by CATChemCore)
  procedure :: init => manager_init
  procedure :: cleanup => manager_cleanup
  procedure :: finalize => manager_finalize
  procedure :: is_ready => manager_is_ready

  ! State object accessors
  procedure :: get_met_state_ptr => manager_get_met_state_ptr
  procedure :: get_chem_state_ptr => manager_get_chem_state_ptr
  ! ... other accessors ...
end type StateManagerType
```

The `StateManagerType` does not own all the state components directly. Instead, it holds pointers to managers like `ConfigManagerType` and `GridManagerType`, which are managed by a higher-level component (`CATChemCore_Mod`).

### State Component Details

#### `MetStateType` - Meteorological State

The `MetStateType`, defined in `metstate_mod.F90`, manages all meteorological data. This includes a large number of 2D and 3D fields for everything from temperature and pressure to soil properties and cloud fractions.

```fortran
! Located in: src/core/metstate_mod.F90

type, public :: MetStateType
  ! ... numerous meteorological fields ...
  real(fp), allocatable :: T(:,:,:)          !< Temperature [K]
  real(fp), allocatable :: U(:,:,:)          !< E/W component of wind [m s-1]
  real(fp), allocatable :: V(:,:,:)          !< N/S component of wind [m s-1]
  real(fp), allocatable :: QV(:,:,:)         !< Specific Humidity [kg/kg]
  real(fp), allocatable :: PS(:,:)           !< Surface Pressure [Pa]
  ! ... and many more ...

contains
  procedure :: init => metstate_init
  procedure :: cleanup => metstate_cleanup
  procedure :: validate => metstate_validate
  procedure :: get_field_ptr => metstate_get_field_ptr
  procedure :: get_column_ptr_func => metstate_get_column_ptr_func
  ! ... other procedures ...
end type MetStateType
```

#### `ChemStateType` - Chemical State

The `ChemStateType`, defined in `chemstate_mod.F90`, manages chemical species and their concentrations. It contains an array of `SpeciesType` objects, where each object holds the data for a single species.

```fortran
! Located in: src/core/chemstate_mod.F90

type, public :: ChemStateType
  ! ...
  integer :: nSpecies          ! Total Number of Species
  character(len=50), allocatable :: SpeciesNames(:)  ! Species Names
  type(SpeciesType), allocatable :: ChemSpecies(:)
  type(GridGeometryType), pointer :: Grid => null()  ! Pointer to grid geometry

contains
  procedure :: init => chemstate_init
  procedure :: cleanup => chemstate_cleanup
  procedure :: find_species => chemstate_find_species
  procedure :: get_concentration => chemstate_get_concentration
  ! ... other procedures ...
end type ChemStateType
```

## Working with State

### Initialization

The `StateManagerType` and its components are initialized by a higher-level driver. The `init` procedures for each state type handle the allocation of their respective data arrays.

```fortran
! Example of initializing the StateManager and its components
subroutine initialize_model_state(state_manager, config, grid, error_mgr, rc)
  type(StateManagerType), intent(inout) :: state_manager
  type(ConfigManagerType), pointer, intent(in) :: config
  type(GridManagerType), pointer, intent(in) :: grid
  type(ErrorManagerType), pointer, intent(inout) :: error_mgr
  integer, intent(out) :: rc

  type(MetStateType), pointer :: met_state
  type(ChemStateType), pointer :: chem_state
  integer :: nx, ny, nlev, max_species

  ! Initialize the state manager
  call state_manager%init('MyStateManager', rc)
  if (rc /= CC_SUCCESS) return

  ! Set external managers
  call state_manager%set_config(config, rc)
  call state_manager%set_grid_manager(grid, rc)

  ! Get grid dimensions
  call grid%get_dimensions(nx, ny, nlev)

  ! Initialize MetState
  met_state => state_manager%get_met_state_ptr()
  call met_state%init(nx, ny, nlev, error_mgr=error_mgr, rc=rc)
  if (rc /= CC_SUCCESS) return

  ! Initialize ChemState
  max_species = config%get_integer('max_species', default=100)
  chem_state => state_manager%get_chem_state_ptr()
  call chem_state%init(max_species, error_mgr=error_mgr, grid=grid%get_geometry(), rc=rc)
  if (rc /= CC_SUCCESS) return

end subroutine initialize_model_state
```

### Accessing State Data

Access to state data is primarily through pointers. This is an efficient way to work with large data arrays without excessive copying.

#### Accessing Meteorological Data

To access a meteorological field, you first get a pointer to the `MetStateType` object, and then get a pointer to the specific field you need.

```fortran
subroutine example_met_access(state_manager, rc)
  type(StateManagerType), intent(inout) :: state_manager
  integer, intent(out) :: rc

  type(MetStateType), pointer :: met_state
  real(fp), pointer :: temp_ptr(:,:,:)
  real(fp), pointer :: temp_column_ptr(:)
  integer :: i, j

  met_state => state_manager%get_met_state_ptr()
  if (.not. associated(met_state)) then
    rc = CC_FAILURE
    return
  end if

  ! Get a pointer to the entire 3D temperature field
  temp_ptr => met_state%get_field_ptr('T')
  if (associated(temp_ptr)) then
    ! Work with the 3D temperature field
    print *, 'Max temperature:', maxval(temp_ptr)
  end if

  ! Get a pointer to a single vertical column of temperature data
  i = 10; j = 20
  temp_column_ptr => met_state%get_column_ptr_func('T', i, j)
  if (associated(temp_column_ptr)) then
    ! Work with the 1D temperature column
    print *, 'Surface temperature at (10, 20):', temp_column_ptr(1)
  end if

end subroutine example_met_access
```

#### Accessing Chemical Species Data

Accessing chemical species data follows a similar pattern. You get a pointer to the `ChemStateType` object, find the index of the species you are interested in, and then get the concentration data.

```fortran
subroutine example_chem_access(state_manager, rc)
  type(StateManagerType), intent(inout) :: state_manager
  integer, intent(out) :: rc

  type(ChemStateType), pointer :: chem_state
  integer :: o3_index
  real(fp), pointer :: o3_conc_ptr(:,:,:)

  chem_state => state_manager%get_chem_state_ptr()
  if (.not. associated(chem_state)) then
    rc = CC_FAILURE
    return
  end if

  ! Find the index for Ozone
  o3_index = chem_state%find_species('O3')
  if (o3_index > 0) then
    ! Get a pointer to the O3 concentration data
    o3_conc_ptr => chem_state%ChemSpecies(o3_index)%conc
    if (associated(o3_conc_ptr)) then
      ! Work with the 3D O3 concentration field
      o3_conc_ptr = o3_conc_ptr * 1.05  ! Increase O3 by 5%
    end if
  end if

end subroutine example_chem_access
```

### State Validation

Each state component has a `validate` procedure that can be used to check for consistency and physical reasonability of the data.

```fortran
subroutine validate_state(state_manager, error_mgr, rc)
  type(StateManagerType), intent(inout) :: state_manager
  type(ErrorManagerType), pointer, intent(inout) :: error_mgr
  integer, intent(out) :: rc

  type(MetStateType), pointer :: met_state
  type(ChemStateType), pointer :: chem_state

  met_state => state_manager%get_met_state_ptr()
  call met_state%validate(error_mgr, rc)
  if (rc /= CC_SUCCESS) return

  chem_state => state_manager%get_chem_state_ptr()
  call chem_state%validate(error_mgr, rc)
  if (rc /= CC_SUCCESS) return

end subroutine validate_state
```

The `StateValidatorUtilsType` also provides a set of generic validation routines that can be used to check dimensions, bounds, and consistency of arrays.

## Best Practices

- **Use Pointers for Access**: Always use pointers to access large data arrays to avoid unnecessary memory copies.
- **Check for Association**: Before dereferencing a pointer, always check if it is associated using the `associated()` intrinsic function.
- **Validate State**: Use the `validate` procedures to ensure the consistency and correctness of the state data before using it in calculations.
- **Lifecycle Management**: The lifecycle of state components (initialization, cleanup) is managed by a higher-level driver. When developing new processes, you can assume that the state is properly initialized.
- **Error Handling**: All state management procedures return a return code (`rc`). Always check the return code and handle errors appropriately using the `ErrorManagerType`.
