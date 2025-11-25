# Diagnostic System Developer Guide

CATChem includes a dynamic and extensible diagnostic system that provides comprehensive monitoring, analysis, and output capabilities. This guide describes the architecture of the diagnostic system and explains how to use it to add diagnostics to new or existing processes.

## Overview

The diagnostic system is designed to be process-driven and flexible. Instead of a centralized configuration file that defines all possible diagnostics, each process is responsible for registering the diagnostics it can provide. Key features include:

- **Centralized Management**: The `DiagnosticManagerType` acts as a central hub for all diagnostics.
- **Process-Driven Registration**: Each process registers its own diagnostic fields programmatically.
- **Dynamic Collection**: The diagnostic manager can collect the current values of all registered diagnostics at any point in the simulation.
- **Configurable Output**: The output of diagnostics (e.g., frequency, directory) is configurable through the main CATChem configuration file.

## Architecture

The diagnostic system is primarily implemented in `src/core/DiagnosticManager_Mod.F90` and `src/core/DiagnosticInterface_Mod.F90`.

### `DiagnosticManagerType`

The `DiagnosticManagerType` is the core of the diagnostic system. It manages a collection of `DiagnosticRegistryType` objects, one for each process that provides diagnostics.

```fortran
! Located in: src/core/DiagnosticManager_Mod.F90

type :: DiagnosticManagerType
  private
  ! Error management
  type(ErrorManagerType), pointer :: error_mgr => null()

  ! Process registries management
  type(DiagnosticRegistryType), allocatable :: process_registries(:)
  character(len=64), allocatable :: process_names(:)
  integer :: num_processes = 0

  ! ... configuration and collection state ...
contains
  procedure :: init => diagnostic_manager_init
  procedure :: register_process => diagnostic_manager_register_process
  procedure :: get_process_registry => diagnostic_manager_get_process_registry
  procedure :: collect_all_diagnostics => diagnostic_manager_collect_all
  procedure :: write_output => diagnostic_manager_write_output
  ! ... other procedures ...
end type DiagnosticManagerType
```

### `DiagnosticRegistryType` and `DiagnosticFieldType`

-   **`DiagnosticRegistryType`**: This type, defined in `DiagnosticInterface_Mod.F90`, acts as a container for all the diagnostic fields provided by a single process.
-   **`DiagnosticFieldType`**: Also defined in `DiagnosticInterface_Mod.F90`, this type represents a single diagnostic quantity. It holds metadata such as the field's name, description, and units, as well as a pointer to the actual data.

## Working with the Diagnostic System

### Registering Diagnostics from a Process

To make a process's internal data available as a diagnostic, the process must register itself with the `DiagnosticManager` and then register its diagnostic fields with its own `DiagnosticRegistryType`.

This is typically done during the initialization phase of the process.

```fortran
! Example of a process initializing its diagnostics
subroutine my_process_init(this, diag_manager, rc)
  use DiagnosticManager_Mod, only: DiagnosticManagerType
  use DiagnosticInterface_Mod, only: DiagnosticRegistryType
  implicit none

  class(MyProcessType), intent(inout) :: this
  type(DiagnosticManagerType), intent(inout) :: diag_manager
  integer, intent(out) :: rc

  type(DiagnosticRegistryType), pointer :: registry
  character(len=64), parameter :: process_name = "MyProcess"

  ! Register this process with the diagnostic manager
  call diag_manager%register_process(process_name, rc)
  if (rc /= 0) return

  ! Get the diagnostic registry for this process
  call diag_manager%get_process_registry(process_name, registry, rc)
  if (rc /= 0) return

  ! Register a scalar diagnostic field
  call registry%register_field("my_scalar_diagnostic", &
    description="A scalar value from MyProcess", &
    units="units", &
    data_ptr=this%internal_scalar_value, &
    rc=rc)

  ! Register a 3D array diagnostic field
  call registry%register_field("my_3d_diagnostic", &
    description="A 3D array from MyProcess", &
    units="kg/m3", &
    data_ptr=this%internal_3d_array, &
    rc=rc)

end subroutine my_process_init
```

In this example, the process registers two diagnostic fields: a scalar and a 3D array. It provides pointers to its internal data (`this%internal_scalar_value` and `this%internal_3d_array`). The diagnostic system does not copy this data; it simply holds a pointer to it. This is a "zero-copy" approach that is very efficient.

### Collecting and Writing Diagnostics

The `DiagnosticManager` is responsible for collecting the data from all registered diagnostics and writing it to output files. This is typically done in the main time-stepping loop of the model.

```fortran
! Example of collecting and writing diagnostics in the main loop
subroutine main_time_step_loop(diag_manager, n_steps)
  use DiagnosticManager_Mod, only: DiagnosticManagerType
  implicit none

  type(DiagnosticManagerType), intent(inout) :: diag_manager
  integer, intent(in) :: n_steps
  integer :: i, rc

  do i = 1, n_steps
    ! ... run model processes ...

    ! Collect diagnostics at this timestep
    call diag_manager%collect_all_diagnostics(rc)

    ! Write output if it's the right time
    call diag_manager%write_output(rc)

    ! Advance the diagnostic manager's timestep counter
    call diag_manager%advance_timestep(rc)
  end do

end subroutine main_time_step_loop
```

The `collect_all_diagnostics` procedure iterates through all registered processes and their diagnostic fields, and it would typically store the current values of the diagnostics in an internal buffer. The `write_output` procedure then writes the buffered data to a file (e.g., in NetCDF format).

### Configuration

The behavior of the diagnostic manager is controlled by settings in the main CATChem configuration file.

```yaml
# Example diagnostics configuration in CATChem_config.yml

diagnostics:
  output_enabled: true
  output_prefix: "catchem_diags"
  output_directory: "./output/diagnostics"
  output_frequency: 3600 # in seconds
```

The `DiagnosticManager` reads this configuration during its initialization and configures its output behavior accordingly.

This process-driven, dynamic diagnostic system provides a powerful and flexible way to monitor and analyze the behavior of CATChem and its various components.