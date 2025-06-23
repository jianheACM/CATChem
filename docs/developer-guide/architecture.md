# Developer Architecture Guide

This guide provides a comprehensive overview of CATChem's software architecture, design principles, and implementation patterns for developers working on the codebase.

## Architectural Overview

CATChem follows a modern, modular architecture designed for:

- **Maintainability**: Clear separation of concerns and well-defined interfaces
- **Extensibility**: Easy addition of new processes and capabilities
- **Performance**: Efficient computation and memory usage patterns
- **Testability**: Comprehensive unit and integration testing support
- **Portability**: Cross-platform compatibility and compiler support

## High-Level Architecture

### System Components

```mermaid
graph TB
    A[User Interface Layer] --> B[API Layer]
    B --> C[Configuration Management]
    B --> D[State Management]
    B --> E[Process Management]

    C --> F[YAML Parser]
    C --> G[Validation Engine]

    D --> H[StateContainer]
    D --> I[Column Virtualization]

    E --> J[Process Registry]
    E --> K[Process Factory]
    E --> L[Process Scheduler]

    J --> M[Chemistry Processes]
    J --> N[Transport Processes]
    J --> O[Emission Processes]
    J --> P[Loss Processes]

    M --> Q[Diagnostic System]
    N --> Q
    O --> Q
    P --> Q

    Q --> R[Output Management]
    Q --> S[Performance Monitoring]
```

### Layer Responsibilities

#### 1. User Interface Layer
- Command-line interface
- Configuration file processing
- Error reporting and logging
- Integration with host models

#### 2. API Layer
- Public interfaces for external integration
- NUOPC, CCPP, and FV3 compatibility
- C/Python binding support
- Version management and compatibility

#### 3. Core System Layer
- Configuration management and validation
- State container and data flow
- Process management and scheduling
- Column virtualization framework

#### 4. Process Layer
- Individual atmospheric processes
- Numerical schemes and algorithms
- Process interdependencies
- Scientific validation and testing

#### 5. Infrastructure Layer
- Diagnostic and monitoring systems
- Performance optimization utilities
- Error handling and recovery
- I/O management and formats

## Design Principles

### 1. Separation of Concerns

Each module has a single, well-defined responsibility:

```fortran
! Configuration is separate from computation
module ConfigManager_Mod
  ! Only handles configuration parsing and validation
end module ConfigManager_Mod

! State management is separate from process logic
module StateInterface_Mod
  ! Only handles data storage and access
end module StateInterface_Mod

! Process logic is separate from infrastructure
module ChemistryProcess_Mod
  ! Only handles chemical calculations
end module ChemistryProcess_Mod
```

### 2. Interface-Based Design

All major components interact through well-defined interfaces:

```fortran
! Abstract process interface
type, abstract :: ProcessInterface_t
contains
  procedure(initialize_interface), deferred :: initialize
  procedure(run_interface), deferred :: run
  procedure(finalize_interface), deferred :: finalize
end type ProcessInterface_t

! Concrete implementations extend the interface
type, extends(ProcessInterface_t) :: ChemistryProcess_t
contains
  procedure :: initialize => chemistry_initialize
  procedure :: run => chemistry_run
  procedure :: finalize => chemistry_finalize
end type ChemistryProcess_t
```

### 3. Dependency Injection

Components receive their dependencies rather than creating them:

```fortran
! Process receives its dependencies at initialization
subroutine chemistry_initialize(this, config, state, diagnostics, rc)
  class(ChemistryProcess_t), intent(inout) :: this
  type(Configuration_t), intent(in) :: config           ! Injected
  type(StateInterface_t), intent(inout) :: state        ! Injected
  type(DiagnosticInterface_t), intent(inout) :: diagnostics ! Injected
  type(ErrorCode_t), intent(out) :: rc
```

### 4. Immutable Configuration

Configuration objects are read-only after initialization:

```fortran
type :: Configuration_t
  private
  type(ConfigurationData_t) :: data
contains
  procedure :: get_value => config_get_value  ! Read-only access
  ! No public setters - configuration is immutable
end type Configuration_t
```

## Core Architecture Components

### Configuration Management

```fortran
module ConfigManager_Mod

  ! Main configuration manager
  type :: ConfigManager_t
    type(ConfigValidator_t) :: validator
    type(ConfigurationTree_t) :: config_tree
    character(len=:), allocatable :: source_file
  contains
    procedure :: load_from_file => cm_load_from_file
    procedure :: validate => cm_validate
    procedure :: get_configuration => cm_get_configuration
  end type ConfigManager_t

  ! Configuration validation
  type :: ConfigValidator_t
    type(SchemaDefinition_t) :: schema
  contains
    procedure :: validate_against_schema => cv_validate_against_schema
    procedure :: check_required_fields => cv_check_required_fields
    procedure :: validate_constraints => cv_validate_constraints
  end type ConfigValidator_t
```

### State Management

```fortran
module StateInterface_Mod

  ! Main state interface
  type :: StateInterface_t
    type(StateContainer_t) :: container
    type(ColumnVirtualization_t) :: column_manager
  contains
    procedure :: get_species => si_get_species
    procedure :: update_species => si_update_species
    procedure :: get_meteorology => si_get_meteorology
    procedure :: update_meteorology => si_update_meteorology
  end type StateInterface_t

  ! State storage container
  type :: StateContainer_t
    type(FieldDictionary_t) :: fields
    type(MetadataDictionary_t) :: metadata
  contains
    procedure :: add_field => sc_add_field
    procedure :: get_field => sc_get_field
    procedure :: update_field => sc_update_field
  end type StateContainer_t
```

### Process Management

```fortran
module ProcessManager_Mod

  ! Process manager coordinates all processes
  type :: ProcessManager_t
    type(ProcessRegistry_t) :: registry
    type(ProcessScheduler_t) :: scheduler
    class(ProcessInterface_t), allocatable :: processes(:)
  contains
    procedure :: initialize => pm_initialize
    procedure :: add_process => pm_add_process
    procedure :: run_processes => pm_run_processes
    procedure :: finalize => pm_finalize
  end type ProcessManager_t

  ! Process scheduling and dependencies
  type :: ProcessScheduler_t
    type(DependencyGraph_t) :: dependencies
    integer, allocatable :: execution_order(:)
  contains
    procedure :: calculate_execution_order => ps_calculate_execution_order
    procedure :: check_dependencies => ps_check_dependencies
  end type ProcessScheduler_t
```

## Column Virtualization Architecture

### Design Philosophy

Column virtualization transforms 3D atmospheric modeling into efficient 1D processing:

```fortran
! Original 3D processing (inefficient)
do i = 1, nx
  do j = 1, ny
    do k = 1, nz
      ! Process point (i,j,k)
      call process_point(data(i,j,k), ...)
    end do
  end do
end do

! Column virtualization (efficient)
do col = 1, num_columns
  ! Process entire column at once
  call process_column(column_data(col), ...)
end do
```

### Column Data Structure

```fortran
type :: ColumnData_t
  ! Spatial information
  real(r8) :: longitude, latitude
  integer :: global_i, global_j
  integer :: num_levels

  ! Atmospheric state
  real(r8), allocatable :: species_concentrations(:,:)  ! (species, levels)
  real(r8), allocatable :: temperature(:)               ! (levels)
  real(r8), allocatable :: pressure(:)                  ! (levels)
  real(r8), allocatable :: density(:)                   ! (levels)

  ! Surface properties
  real(r8) :: surface_pressure
  real(r8) :: surface_temperature
  type(SurfaceProperties_t) :: surface

  ! Process-specific workspace
  type(ProcessWorkspace_t), allocatable :: workspace(:)
```

### Column Processing Framework

```fortran
module ColumnVirtualization_Mod

  type :: ColumnVirtualization_t
    integer :: num_columns
    type(ColumnData_t), allocatable :: columns(:)
    type(ColumnMapper_t) :: mapper
  contains
    procedure :: initialize => cv_initialize
    procedure :: map_from_3d => cv_map_from_3d
    procedure :: map_to_3d => cv_map_to_3d
    procedure :: process_columns => cv_process_columns
  end type ColumnVirtualization_t

  ! Efficient column processing with parallelization
  subroutine cv_process_columns(this, process_kernel, time_step, rc)
    class(ColumnVirtualization_t), intent(inout) :: this
    procedure(column_process_interface) :: process_kernel
    real(r8), intent(in) :: time_step
    type(ErrorCode_t), intent(out) :: rc

    integer :: col

    !$OMP PARALLEL DO PRIVATE(col) SCHEDULE(DYNAMIC)
    do col = 1, this%num_columns
      call process_kernel(this%columns(col), time_step, rc)
      if (rc%is_error()) then
        !$OMP CRITICAL
        call handle_column_error(col, rc)
        !$OMP END CRITICAL
      end if
    end do
    !$OMP END PARALLEL DO
  end subroutine cv_process_columns
```

## Error Handling Architecture

### Error Code System

```fortran
module ErrorHandling_Mod

  ! Comprehensive error information
  type :: ErrorCode_t
    integer :: code = 0
    character(len=:), allocatable :: message
    character(len=:), allocatable :: context(:)
    logical :: is_recoverable = .false.
  contains
    procedure :: set_success => ec_set_success
    procedure :: set_error => ec_set_error
    procedure :: set_warning => ec_set_warning
    procedure :: add_context => ec_add_context
    procedure :: is_error => ec_is_error
    procedure :: is_success => ec_is_success
  end type ErrorCode_t

  ! Error handling patterns
  interface
    subroutine error_handler_interface(error_code, recovery_action)
      import :: ErrorCode_t
      type(ErrorCode_t), intent(in) :: error_code
      character(len=*), intent(out) :: recovery_action
    end subroutine error_handler_interface
  end interface
```

### Error Recovery Strategies

```fortran
! Graceful degradation example
subroutine robust_chemistry_solve(this, column_data, time_step, rc)
  class(ChemistryProcess_t), intent(inout) :: this
  type(ColumnData_t), intent(inout) :: column_data
  real(r8), intent(in) :: time_step
  type(ErrorCode_t), intent(out) :: rc

  type(ErrorCode_t) :: local_rc
  real(r8) :: reduced_timestep

  ! Try normal integration
  call this%solve_chemistry(column_data, time_step, local_rc)

  if (local_rc%is_error() .and. local_rc%is_recoverable()) then
    ! Try with reduced timestep
    reduced_timestep = time_step * 0.1_r8
    call this%solve_chemistry_subcycling(column_data, reduced_timestep, &
                                        int(time_step/reduced_timestep), local_rc)

    if (local_rc%is_success()) then
      call rc%set_warning("Used subcycling for numerical stability")
      return
    end if
  end if

  if (local_rc%is_error()) then
    ! Use fallback scheme
    call this%solve_chemistry_fallback(column_data, time_step, local_rc)
    if (local_rc%is_success()) then
      call rc%set_warning("Used fallback chemistry scheme")
      return
    end if
  end if

  ! If all recovery attempts fail
  call rc%set_error("Chemistry solution failed: " // local_rc%get_message())
end subroutine robust_chemistry_solve
```

## Performance Architecture

### Memory Management

```fortran
module MemoryManager_Mod

  ! Memory pool for efficient allocation
  type :: MemoryPool_t
    real(r8), allocatable :: pool(:)
    integer, allocatable :: allocation_map(:)
    integer :: pool_size, next_free
  contains
    procedure :: allocate_chunk => mp_allocate_chunk
    procedure :: deallocate_chunk => mp_deallocate_chunk
    procedure :: get_statistics => mp_get_statistics
  end type MemoryPool_t

  ! NUMA-aware allocation
  type :: NUMAAllocator_t
    type(MemoryPool_t), allocatable :: pools_per_node(:)
    integer :: num_numa_nodes
  contains
    procedure :: allocate_on_node => na_allocate_on_node
    procedure :: get_local_node => na_get_local_node
  end type NUMAAllocator_t
```

### Computational Kernels

```fortran
! Optimized computational kernels
module ComputationalKernels_Mod

  ! Vectorized operations
  interface vectorized_operation
    module procedure vector_add_r8, vector_multiply_r8, vector_solve_r8
  end interface vectorized_operation

  ! SIMD-optimized vector operations
  pure subroutine vector_add_r8(a, b, result, n)
    integer, intent(in) :: n
    real(r8), intent(in) :: a(n), b(n)
    real(r8), intent(out) :: result(n)

    integer :: i

    !$OMP SIMD
    do i = 1, n
      result(i) = a(i) + b(i)
    end do
    !$OMP END SIMD
  end subroutine vector_add_r8
```

## Testing Architecture

### Test Framework

```fortran
module TestFramework_Mod

  ! Test suite management
  type :: TestSuite_t
    character(len=:), allocatable :: name
    type(TestCase_t), allocatable :: test_cases(:)
    integer :: num_tests, num_passed, num_failed
  contains
    procedure :: add_test => ts_add_test
    procedure :: run_tests => ts_run_tests
    procedure :: generate_report => ts_generate_report
  end type TestSuite_t

  ! Individual test case
  type :: TestCase_t
    character(len=:), allocatable :: name
    procedure(test_procedure_interface), pointer :: test_proc
    logical :: passed = .false.
    character(len=:), allocatable :: failure_message
  contains
    procedure :: run => tc_run
    procedure :: assert_equal => tc_assert_equal
    procedure :: assert_near => tc_assert_near
  end type TestCase_t
```

### Test Categories

```fortran
! Unit tests for individual components
module UnitTests_Mod
  procedure :: test_configuration_loading
  procedure :: test_state_management
  procedure :: test_column_virtualization
  procedure :: test_process_execution
end module UnitTests_Mod

! Integration tests for component interactions
module IntegrationTests_Mod
  procedure :: test_full_process_chain
  procedure :: test_nuopc_integration
  procedure :: test_error_handling
end module IntegrationTests_Mod

! Performance tests for optimization validation
module PerformanceTests_Mod
  procedure :: test_column_processing_speed
  procedure :: test_memory_usage
  procedure :: test_parallel_scalability
end module PerformanceTests_Mod
```

## Documentation Architecture

### Code Documentation Standards

```fortran
!> @brief Brief description of the subroutine
!> @details Detailed description of the algorithm, assumptions, and usage
!> @param[in] input_param Description of input parameter
!> @param[out] output_param Description of output parameter
!> @param[inout] inout_param Description of input/output parameter
!> @return Description of return value
!> @author Developer Name
!> @date Creation date
!> @version Version information
!> @see Related procedures or documentation
!> @warning Important warnings or limitations
!> @note Additional notes or implementation details
subroutine well_documented_procedure(input_param, output_param, inout_param)
  integer, intent(in) :: input_param
  real(r8), intent(out) :: output_param
  type(SomeType_t), intent(inout) :: inout_param
```

### API Documentation

```fortran
!> @defgroup api_core Core API
!> Core functionality for CATChem integration
!> @{

!> @brief Initialize CATChem instance
!> @details This function initializes a CATChem instance with the provided
!> configuration. The instance must be initialized before any other operations.
!> @param[out] instance CATChem instance handle
!> @param[in] config_file Path to configuration file
!> @return Status code (0 = success, non-zero = error)
integer function catchem_initialize(instance, config_file) result(status)
  type(CATChemAPI_t), intent(out) :: instance
  character(len=*), intent(in) :: config_file
end function catchem_initialize

!> @}
```

## Future Architecture Considerations

### Extensibility Points

1. **Plugin Architecture**: Dynamic loading of new processes
2. **GPU Acceleration**: CUDA/OpenCL integration points
3. **Machine Learning**: AI/ML algorithm integration
4. **Cloud Computing**: Containerization and orchestration
5. **Data Formats**: Support for new I/O formats

### Scalability Roadmap

1. **Hybrid Parallelization**: MPI + OpenMP + GPU
2. **Asynchronous Processing**: Overlap computation and communication
3. **Dynamic Load Balancing**: Adaptive work distribution
4. **Memory Hierarchy**: Cache-aware algorithms
5. **Exascale Computing**: Preparation for next-generation HPC

## Related Documentation

- [Build System](build-system.md)
- [Process Architecture](processes/architecture.md)
- [State Management](core/state-management.md)
- [Column Virtualization](../guides/column-virtualization.md)
- [Performance Guide](performance.md)

---

*This architecture guide provides the foundation for understanding CATChem's design and implementation. For specific implementation details, consult the individual component documentation.*
