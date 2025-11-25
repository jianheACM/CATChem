# Column Virtualization Developer Guide

CATChem's column virtualization system is a key architectural feature that enables efficient processing of atmospheric chemistry. It treats the 3D model domain as a collection of independent 1D vertical columns. This guide describes the column virtualization system and how to work with it as a developer.

## Overview

Column virtualization offers several advantages:

- **Reduced Memory Pressure**: By processing one column at a time, the working memory footprint is significantly reduced.
- **Improved Cache Performance**: Accessing data in a contiguous 1D column improves cache locality.
- **Simplified Physics**: Chemical and physical processes can be written to operate on simple 1D vertical profiles.
- **Efficient Parallelization**: Since columns are independent, they can be processed in parallel, for example using OpenMP.

## Architecture

The column virtualization architecture is centered around three main components defined in `src/core/VirtualColumn_Mod.F90` and `src/core/ColumnInterface_Mod.F90`:

- **`VirtualColumnType`**: A data container for a single vertical column.
- **`ColumnViewType`**: An interface for accessing column data from different grid types.
- **`ColumnProcessorType`**: A manager for processing collections of virtual columns.

### `VirtualColumnType`

The `VirtualColumnType` is the fundamental data structure for column virtualization. It holds all the data for a single vertical column.

A key feature of the `VirtualColumnType` is the `VirtualMetType`, which provides access to meteorological data. The `VirtualMetType` contains pointers to the main meteorological data arrays in `MetStateType`. This "zero-copy" approach is highly efficient as it avoids copying large amounts of data for each column.

```fortran
! Located in: src/core/VirtualColumn_Mod.F90

! Virtual meteorological data container with direct pointers
type :: VirtualMetType
  ! Pointers to meteorological fields, e.g.:
  real(fp), pointer :: T(:)
  real(fp), pointer :: U(:)
  real(fp), pointer :: V(:)
  ! ... and many more
end type VirtualMetType

type :: VirtualColumnType
  ! Meteorological data (pointers managed by VirtualMetType)
  type(VirtualMetType) :: met

  ! Chemical and emission data (copied into arrays for modification)
  real(fp), allocatable :: chem_data(:,:)  ! (nlev, nspec)
  real(fp), allocatable :: emis_data(:,:)  ! (nlev, nspec)

  ! Grid position and metadata
  integer :: grid_i, grid_j
  real(fp) :: lat, lon, area

  ! ... dimensions and status ...
contains
  procedure :: init => virtual_column_init
  procedure :: get_met => virtual_column_get_met
  procedure :: get_chem_field => virtual_column_get_chem_field
  ! ... other procedures ...
end type VirtualColumnType
```

Chemical and emission data, on the other hand, are copied into arrays within the `VirtualColumnType`. This is because these fields are typically modified by the processes, and copying them ensures that modifications to one column do not affect others being processed in parallel.

### Creating and Populating a Virtual Column

A `VirtualColumnType` is created and populated from the main `StateManagerType`. The `StateManagerType` has a `create_virtual_column` procedure that initializes a `VirtualColumnType` for a given grid location `(i, j)`.

```fortran
! Example of creating a virtual column from the StateManager
subroutine create_and_process_column(state_manager, i, j, rc)
  type(StateManagerType), intent(inout) :: state_manager
  integer, intent(in) :: i, j
  integer, intent(out) :: rc

  type(VirtualColumnType) :: column

  ! Create and populate the virtual column for grid cell (i, j)
  call state_manager%create_virtual_column(i, j, column, rc)
  if (rc /= CC_SUCCESS) return

  ! Now the 'column' variable holds the data for the specified column.
  ! The met data is accessible via pointers, e.g., column%met%T(:)
  ! The chemistry data is in column%chem_data(:,:)

  ! ... process the column ...

  ! Apply the changes back to the main state
  call state_manager%apply_virtual_column(column, rc)

  ! Clean up the column
  call column%cleanup()

end subroutine create_and_process_column
```

The `create_virtual_column` procedure performs the following steps:
1.  Initializes a `VirtualColumnType` with the correct dimensions.
2.  Populates the `VirtualMetType` by pointing its components to the appropriate slices of the 3D meteorological fields in `MetStateType`.
3.  Copies the chemical and emission data for the column into the `chem_data` and `emis_data` arrays.

After a process has modified the `chem_data` or `emis_data` in the `VirtualColumnType`, the changes are copied back to the main `ChemStateType` and `EmisStateType` using the `apply_virtual_column` procedure.

### `ColumnProcessorType` and Batch Processing

The `ColumnProcessorType` is designed to manage and process a collection of `VirtualColumnType` objects. This is particularly useful for parallel processing, where you might want to process multiple columns simultaneously across different threads.

```fortran
! Located in: src/core/ColumnInterface_Mod.F90

type :: ColumnProcessorType
  private
  type(VirtualColumnType), allocatable :: columns(:)  !< Array of virtual columns
  integer :: n_columns = 0
contains
  procedure :: init => processor_init
  procedure :: add_column => processor_add_column
  procedure :: process_all => processor_process_all
  ! ... other procedures ...
end type ColumnProcessorType
```

The `process_all` procedure takes a subroutine as an argument, which it then applies to every column in its collection. This subroutine must conform to the `column_process_interface`.

```fortran
! Located in: src/core/ColumnInterface_Mod.F90

abstract interface
  subroutine column_process_interface(column, rc)
    import :: VirtualColumnType
    type(VirtualColumnType), intent(inout) :: column
    integer, intent(out) :: rc
  end subroutine column_process_interface
end interface
```

A typical workflow for parallel processing would be:
1.  Create a `ColumnProcessorType` object.
2.  In a loop over the horizontal grid, create a `VirtualColumnType` for each grid cell and add it to the `ColumnProcessorType`.
3.  Call `process_all` with a pointer to your processing subroutine. The `process_all` procedure can then loop over the columns in parallel (e.g., using OpenMP).

```fortran
subroutine run_chemistry_process_parallel(state_manager, rc)
  type(StateManagerType), intent(inout) :: state_manager
  integer, intent(out) :: rc

  type(ColumnProcessorType) :: processor
  type(VirtualColumnType) :: column
  integer :: i, j, nx, ny, nlev

  call state_manager%get_grid_manager()%get_dimensions(nx, ny, nlev)

  ! Initialize the processor
  call processor%init(nx * ny, rc)

  ! Create virtual columns for the entire domain
  do j = 1, ny
    do i = 1, nx
      call state_manager%create_virtual_column(i, j, column, rc)
      call processor%add_column(column, rc)
    end do
  end do

  ! Process all columns in parallel
  !$OMP PARALLEL DO
  do i = 1, processor%n_columns
    call my_chemistry_subroutine(processor%columns(i), rc)
  end do
  !$OMP END PARALLEL DO

  ! Apply changes back to the main state
  do i = 1, processor%n_columns
    call state_manager%apply_virtual_column(processor%columns(i), rc)
  end do

  call processor%cleanup()

end subroutine run_chemistry_process_parallel
```

## Writing a Column-Based Process

To write a new process that operates on a column, you need to create a subroutine that takes a `VirtualColumnType` as an argument.

```fortran
subroutine my_chemistry_subroutine(column, rc)
  use Precision_Mod, only: fp
  use VirtualColumn_Mod, only: VirtualColumnType
  implicit none

  type(VirtualColumnType), intent(inout) :: column
  integer, intent(out) :: rc

  integer :: k, o3_idx, no2_idx
  real(fp) :: temp_k, o3_conc, no2_conc, reaction_rate

  rc = 0
  o3_idx = 1  ! Placeholder for O3 species index
  no2_idx = 2 ! Placeholder for NO2 species index

  ! Loop over vertical levels in the column
  do k = 1, column%nlev
    ! Access meteorological data via pointers
    temp_k = column%met%T(k)

    ! Access chemistry data from the copied array
    o3_conc = column%chem_data(k, o3_idx)
    no2_conc = column%chem_data(k, no2_idx)

    ! Calculate a simple reaction rate
    reaction_rate = 1.0e-12 * exp(-1400.0 / temp_k) * o3_conc * no2_conc

    ! Update the chemistry data array
    column%chem_data(k, o3_idx) = o3_conc - reaction_rate
    column%chem_data(k, no2_idx) = no2_conc - reaction_rate
  end do

end subroutine my_chemistry_subroutine
```

This modular, column-based approach allows for the development of complex physical and chemical processes that are efficient, scalable, and easy to test and maintain.
