# FV3 Integration

This guide covers integrating CATChem with the FV3 (Finite-Volume Cubed-Sphere) dynamical core used in UFS and GFS models.

## Overview

FV3 integration enables CATChem to work with the FV3 atmospheric model, providing:

- **Cubed-Sphere Grid Support** - Native support for FV3's cubed-sphere grid
- **Efficient Data Exchange** - Optimized interfaces for FV3 data structures
- **Physics Integration** - Seamless integration with FV3 physics packages
- **Performance Optimization** - Optimized for FV3's computational patterns

## FV3 Interface Architecture

### FV3 Chemistry Interface

```fortran
module fv3_catchem_interface
  use fv_arrays_mod, only: fv_atmos_type
  use fv_grid_tools_mod, only: cubed_to_latlon
  use catchem_mod
  use catchem_types
  implicit none

  private

  public :: fv3_catchem_init
  public :: fv3_catchem_run
  public :: fv3_catchem_finalize

contains

  subroutine fv3_catchem_init(Atm, rc)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    integer, intent(out) :: rc

    integer :: n, npx, npy, npz
    type(catchem_config_type) :: config

    rc = 0

    ! Initialize for each atmosphere tile
    do n = 1, size(Atm)
      associate(atm_data => Atm(n))

        ! Get grid dimensions
        npx = atm_data%npx
        npy = atm_data%npy
        npz = atm_data%npz

        ! Setup CATChem configuration for this tile
        call setup_catchem_config_fv3(atm_data, config, rc)
        if (rc /= 0) return

        ! Initialize CATChem for this tile
        call catchem_init_tile(config, n, npx, npy, npz, rc)
        if (rc /= 0) return

      end associate
    end do

  end subroutine fv3_catchem_init

  subroutine fv3_catchem_run(Atm, dt, rc)
    type(fv_atmos_type), intent(inout) :: Atm(:)
    real, intent(in) :: dt
    integer, intent(out) :: rc

    integer :: n

    rc = 0

    ! Process each atmosphere tile
    do n = 1, size(Atm)
      call fv3_catchem_run_tile(Atm(n), n, dt, rc)
      if (rc /= 0) return
    end do

  end subroutine fv3_catchem_run

end module fv3_catchem_interface
```

### Cubed-Sphere Grid Handling

```fortran
subroutine fv3_catchem_run_tile(atm_data, tile_id, dt, rc)
  type(fv_atmos_type), intent(inout) :: atm_data
  integer, intent(in) :: tile_id
  real, intent(in) :: dt
  integer, intent(out) :: rc

  real, allocatable :: catchem_state(:,:,:,:)
  real, pointer :: temp(:,:,:), pres(:,:,:), humid(:,:,:)
  real, pointer :: o3(:,:,:), no2(:,:,:), co(:,:,:)
  integer :: i, j, k, isc, iec, jsc, jec, npz

  rc = 0

  ! Get array bounds for this tile
  isc = atm_data%bd%isc
  iec = atm_data%bd%iec
  jsc = atm_data%bd%jsc
  jec = atm_data%bd%jec
  npz = atm_data%npz

  ! Get meteorological fields from FV3
  temp => atm_data%pt    ! Potential temperature
  pres => atm_data%delp  ! Layer pressure thickness
  humid => atm_data%q    ! Specific humidity

  ! Get chemical tracers from FV3
  o3 => atm_data%q(:,:,:,get_tracer_index('sphum', 'o3'))
  no2 => atm_data%q(:,:,:,get_tracer_index('sphum', 'no2'))
  co => atm_data%q(:,:,:,get_tracer_index('sphum', 'co'))

  ! Convert FV3 data to CATChem format
  call convert_fv3_to_catchem(temp, pres, humid, o3, no2, co, &
                              isc, iec, jsc, jec, npz, &
                              catchem_state, rc)
  if (rc /= 0) return

  ! Process chemistry for this tile
  !$OMP PARALLEL DO PRIVATE(i,j) COLLAPSE(2)
  do j = jsc, jec
    do i = isc, iec
      call catchem_process_column(catchem_state(i,j,:,:), dt, rc)
    end do
  end do
  !$OMP END PARALLEL DO

  ! Convert results back to FV3 format
  call convert_catchem_to_fv3(catchem_state, &
                              o3, no2, co, &
                              isc, iec, jsc, jec, npz, rc)

end subroutine fv3_catchem_run_tile
```

## Grid and Coordinate Transformations

### Cubed-Sphere to Lat-Lon Conversion

```fortran
subroutine convert_cubed_sphere_winds(atm_data, u_latlon, v_latlon)
  type(fv_atmos_type), intent(in) :: atm_data
  real, intent(out) :: u_latlon(:,:,:), v_latlon(:,:,:)

  real :: cosa_u(:,:), sina_u(:,:), cosa_v(:,:), sina_v(:,:)
  real :: u_native(:,:,:), v_native(:,:,:)
  integer :: i, j, k

  ! Get native winds from FV3
  u_native = atm_data%u
  v_native = atm_data%v

  ! Get grid transformation coefficients
  cosa_u = atm_data%gridstruct%cosa_u
  sina_u = atm_data%gridstruct%sina_u
  cosa_v = atm_data%gridstruct%cosa_v
  sina_v = atm_data%gridstruct%sina_v

  ! Transform winds to lat-lon orientation
  do k = 1, size(u_native, 3)
    do j = 1, size(u_native, 2)
      do i = 1, size(u_native, 1)
        u_latlon(i,j,k) = u_native(i,j,k) * cosa_u(i,j) - &
                          v_native(i,j,k) * sina_u(i,j)
        v_latlon(i,j,k) = u_native(i,j,k) * sina_v(i,j) + &
                          v_native(i,j,k) * cosa_v(i,j)
      end do
    end do
  end do

end subroutine convert_cubed_sphere_winds
```

### Pressure Coordinate Conversion

```fortran
subroutine convert_fv3_pressure_coords(atm_data, pressure_levels, rc)
  type(fv_atmos_type), intent(in) :: atm_data
  real, intent(out) :: pressure_levels(:,:,:)
  integer, intent(out) :: rc

  real :: ptop, ak(:), bk(:)
  real :: surface_pressure(:,:)
  integer :: i, j, k, npz

  rc = 0
  npz = atm_data%npz

  ! Get FV3 vertical coordinate parameters
  ptop = atm_data%ak(1)  ! Model top pressure
  ak = atm_data%ak       ! Hybrid coordinate A values
  bk = atm_data%bk       ! Hybrid coordinate B values

  ! Get surface pressure
  surface_pressure = atm_data%ps

  ! Calculate pressure at layer interfaces
  do k = 1, npz + 1
    do j = 1, size(surface_pressure, 2)
      do i = 1, size(surface_pressure, 1)
        pressure_levels(i,j,k) = ak(k) + bk(k) * surface_pressure(i,j)
      end do
    end do
  end do

end subroutine convert_fv3_pressure_coords
```

## Tracer Integration

### FV3 Tracer Registration

```fortran
subroutine register_catchem_tracers_fv3(atm_data, rc)
  type(fv_atmos_type), intent(inout) :: atm_data
  integer, intent(out) :: rc

  integer :: tracer_id
  character(len=32) :: tracer_names(10)
  character(len=64) :: tracer_units(10)
  character(len=128) :: tracer_longnames(10)
  integer :: i, num_tracers

  rc = 0

  ! Define CATChem tracers for FV3
  tracer_names = ['o3    ', 'no    ', 'no2   ', 'co    ', 'so2   ', &
                  'pm25  ', 'pm10  ', 'dust1 ', 'dust2 ', 'dust3 ']

  tracer_units = ['mol/mol', 'mol/mol', 'mol/mol', 'mol/mol', 'mol/mol', &
                  'kg/kg  ', 'kg/kg  ', 'kg/kg  ', 'kg/kg  ', 'kg/kg  ']

  tracer_longnames = ['Ozone                    ', &
                      'Nitric Oxide             ', &
                      'Nitrogen Dioxide         ', &
                      'Carbon Monoxide          ', &
                      'Sulfur Dioxide           ', &
                      'PM 2.5                   ', &
                      'PM 10                    ', &
                      'Dust Size Bin 1          ', &
                      'Dust Size Bin 2          ', &
                      'Dust Size Bin 3          ']

  num_tracers = 10

  ! Register tracers with FV3 tracer manager
  do i = 1, num_tracers
    call register_tracer(atm_data%tracer_manager, &
                        trim(tracer_names(i)), &
                        trim(tracer_units(i)), &
                        trim(tracer_longnames(i)), &
                        tracer_id, rc)
    if (rc /= 0) then
      print *, 'Failed to register tracer: ', trim(tracer_names(i))
      return
    end if

    ! Store tracer ID for later use
    call store_catchem_tracer_id(trim(tracer_names(i)), tracer_id)
  end do

  print *, 'Successfully registered ', num_tracers, ' CATChem tracers with FV3'

end subroutine register_catchem_tracers_fv3
```

### Tracer Transport

```fortran
subroutine apply_fv3_transport_to_catchem_tracers(atm_data, dt, rc)
  type(fv_atmos_type), intent(inout) :: atm_data
  real, intent(in) :: dt
  integer, intent(out) :: rc

  real, pointer :: tracer_data(:,:,:)
  integer :: tracer_id, i
  character(len=32) :: tracer_names(10)

  rc = 0
  tracer_names = ['o3', 'no', 'no2', 'co', 'so2', &
                  'pm25', 'pm10', 'dust1', 'dust2', 'dust3']

  ! Apply FV3 transport to each CATChem tracer
  do i = 1, size(tracer_names)
    ! Get tracer ID
    tracer_id = get_catchem_tracer_id(trim(tracer_names(i)))
    if (tracer_id <= 0) cycle

    ! Get tracer data pointer
    tracer_data => atm_data%q(:,:,:,tracer_id)

    ! Apply FV3 advection (handled by FV3 core)
    ! This is done automatically by FV3 for all tracers in q array

    ! Apply additional transport processes if needed
    call apply_vertical_mixing(atm_data, tracer_data, dt, rc)
    if (rc /= 0) return

    call apply_convective_transport(atm_data, tracer_data, dt, rc)
    if (rc /= 0) return
  end do

end subroutine apply_fv3_transport_to_catchem_tracers
```

## Physics Integration

### Integration with FV3 Physics

```fortran
subroutine integrate_catchem_with_fv3_physics(atm_data, physics_data, dt, rc)
  type(fv_atmos_type), intent(inout) :: atm_data
  type(fv3_physics_type), intent(inout) :: physics_data
  real, intent(in) :: dt
  integer, intent(out) :: rc

  real :: cloud_fraction(:,:,:)
  real :: cloud_water(:,:,:)
  real :: precipitation_rate(:,:)
  real :: boundary_layer_height(:,:)

  rc = 0

  ! Get cloud information from FV3 physics
  cloud_fraction = physics_data%cloud_fraction
  cloud_water = physics_data%cloud_water
  precipitation_rate = physics_data%precip_rate
  boundary_layer_height = physics_data%pbl_height

  ! Use cloud information for aqueous chemistry
  call setup_aqueous_chemistry_from_clouds(cloud_fraction, cloud_water, rc)
  if (rc /= 0) return

  ! Use precipitation for wet deposition
  call apply_wet_deposition(precipitation_rate, dt, rc)
  if (rc /= 0) return

  ! Use boundary layer height for mixing
  call apply_boundary_layer_mixing(boundary_layer_height, dt, rc)
  if (rc /= 0) return

  ! Provide chemistry feedback to radiation
  call update_radiation_from_chemistry(atm_data, physics_data, rc)
  if (rc /= 0) return

end subroutine integrate_catchem_with_fv3_physics
```

## Configuration

### FV3 Configuration (`fv3_catchem_config.yml`)

```yaml
# FV3-specific CATChem configuration
fv3_integration:
  enabled: true

  # Grid configuration
  grid:
    cubed_sphere: true
    tile_count: 6
    grid_type: "C"  # C-grid staggering

  # Tracer configuration
  tracers:
    transport_by_fv3: true
    passive_tracers: [o3, no, no2, co, so2]
    aerosol_tracers: [pm25, pm10, dust1, dust2, dust3]

  # Integration settings
  integration:
    chemistry_timestep: 60.0  # seconds
    transport_timestep: -1    # Use FV3 dynamics timestep
    coupling_method: "operator_splitting"

  # Physics coupling
  physics_coupling:
    aqueous_chemistry: true
    wet_deposition: true
    radiation_feedback: true
    boundary_layer_mixing: true

  # Performance settings
  performance:
    tile_threading: true
    column_parallelization: true
    memory_optimization: true

  # Output configuration
  output:
    native_grid: true       # Output on cubed-sphere grid
    interpolate_to_latlon: false
    diagnostic_level: 1     # 0=minimal, 1=standard, 2=detailed
```

### FV3 Namelist Integration

```fortran
! In FV3 namelist (input.nml)
&gfs_physics_nml
  ! ... other physics options ...

  ! Enable chemistry
  chem_on = .true.
  chem_scheme = 'catchem'
  chem_dt = 60.0

  ! Chemistry tracers
  nchem = 10
  chem_tracer_names = 'o3', 'no', 'no2', 'co', 'so2',
                      'pm25', 'pm10', 'dust1', 'dust2', 'dust3'
/

&fv_core_nml
  ! ... other core options ...

  ! Tracer transport
  kord_tr = 9        ! Tracer transport scheme
  fill = .true.      ! Monotonic transport
  consv_te = 1       ! Total energy conservation

  ! Chemistry-related options
  external_ic = .false.
  mountain = .true.
  warm_start = .false.
/
```

## Testing and Validation

### FV3 Integration Tests

```fortran
program test_fv3_integration
  use fv_arrays_mod
  use fv3_catchem_interface
  implicit none

  type(fv_atmos_type) :: Atm(1)  ! Single tile test
  real :: dt = 60.0
  integer :: rc
  logical :: test_passed = .true.

  ! Setup test atmosphere
  call setup_test_fv3_atmosphere(Atm(1), rc)
  if (rc /= 0) then
    print *, 'Failed to setup test atmosphere'
    test_passed = .false.
  end if

  ! Test initialization
  call fv3_catchem_init(Atm, rc)
  if (rc /= 0) then
    print *, 'FV3 CATChem initialization failed'
    test_passed = .false.
  end if

  ! Test run
  call fv3_catchem_run(Atm, dt, rc)
  if (rc /= 0) then
    print *, 'FV3 CATChem run failed'
    test_passed = .false.
  end if

  ! Test finalization
  call fv3_catchem_finalize(Atm, rc)
  if (rc /= 0) then
    print *, 'FV3 CATChem finalization failed'
    test_passed = .false.
  end if

  if (test_passed) then
    print *, 'FV3 integration tests PASSED'
  else
    print *, 'FV3 integration tests FAILED'
    stop 1
  end if

end program test_fv3_integration
```

## Performance Optimization

### Tile-Level Parallelization

```fortran
subroutine fv3_catchem_run_parallel(Atm, dt, rc)
  type(fv_atmos_type), intent(inout) :: Atm(:)
  real, intent(in) :: dt
  integer, intent(out) :: rc

  integer :: n, num_tiles
  integer :: tile_rc(size(Atm))

  rc = 0
  num_tiles = size(Atm)

  ! Process tiles in parallel
  !$OMP PARALLEL DO PRIVATE(n) SCHEDULE(DYNAMIC)
  do n = 1, num_tiles
    call fv3_catchem_run_tile(Atm(n), n, dt, tile_rc(n))
  end do
  !$OMP END PARALLEL DO

  ! Check for errors
  do n = 1, num_tiles
    if (tile_rc(n) /= 0) then
      print *, 'Error in tile ', n, ': ', tile_rc(n)
      rc = tile_rc(n)
      return
    end if
  end do

end subroutine fv3_catchem_run_parallel
```

## Best Practices

### 1. Memory Management
```fortran
! Use FV3 memory layout for optimal performance
real, pointer :: tracer_ptr(:,:,:) => null()
tracer_ptr => atm_data%q(:,:,:,tracer_id)
! Work directly with pointer, avoid copying
```

### 2. Grid Awareness
```fortran
! Respect FV3 grid staggering and halos
isc = atm_data%bd%isc  ! Start index (computation domain)
iec = atm_data%bd%iec  ! End index (computation domain)
! Don't process halo regions unless necessary
```

### 3. Coordinate Systems
```fortran
! Be aware of FV3's hybrid sigma-pressure coordinates
! Convert to standard pressure when needed for chemistry
call convert_hybrid_to_pressure(atm_data%ak, atm_data%bk, atm_data%ps, pressure)
```

This FV3 integration provides efficient coupling between CATChem and the FV3 dynamical core while respecting FV3's computational patterns and data structures.
