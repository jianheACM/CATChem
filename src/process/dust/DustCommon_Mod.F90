!> \file DustCommon_Mod.F90
!! \brief Common types and utilities for dust process
!!
!! This module defines the configuration and state types used by the
!! dust process and its schemes.
!!
!! Generated on: 2025-07-09T12:43:17.982946
!! Author: Barry Baker
!! Version: 1.0.0

module DustCommon_Mod

   use iso_fortran_env, only: fp => real64
   use ErrorHandler_Mod

   implicit none
   private

   ! Export types
   public :: DustConfig
   public :: DustState
   public :: DustSchemeFENGSHAConfig
   public :: DustSchemeFENGSHAState
   public :: DustSchemeGINOUXConfig
   public :: DustSchemeGINOUXState

   ! Export utility functions
   public :: int_to_string

   !> Main configuration type for dust process
   type :: DustConfig

      ! Process settings
      character(len=32) :: active_scheme = ''
      logical :: is_active = .true.
      real(fp) :: dt_min = 1.0_fp     ! Minimum time step (seconds)
      real(fp) :: dt_max = 3600.0_fp  ! Maximum time step (seconds)

      ! Species configuration
      integer :: n_species = 0
      character(len=32) :: species_names(1)



      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: init => init_dust_config
      procedure, public :: validate => validate_dust_config
      procedure, public :: finalize => finalize_dust_config
      procedure, public :: print_summary => print_dust_config_summary
   end type DustConfig

   !> Main state type for dust process
   type :: DustState

      ! Runtime state
      logical :: is_initialized = .false.
      real(fp) :: current_time = 0.0_fp
      real(fp) :: last_update_time = 0.0_fp
      integer :: n_columns = 0
      integer :: n_levels = 0

      ! Working arrays (allocated during initialization)
      real(fp), allocatable :: emission_rates(:,:)    ! (n_species, n_columns)
      real(fp), allocatable :: total_emissions(:)     ! (n_columns)

      ! Diagnostic arrays
      real(fp), allocatable :: total_dust_emission(:)  ! Total dust emissions for all species

   contains
      procedure, public :: init => init_dust_state
      procedure, public :: allocate_arrays => allocate_dust_state_arrays
      procedure, public :: reset => reset_dust_state
      procedure, public :: finalize => finalize_dust_state
   end type DustState

   !> Configuration type for fengsha scheme
   type :: DustSchemeFENGSHAConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'fengsha'
      character(len=256) :: description = 'Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Scheme parameters
      real(fp) :: alpha = 0.16  ! linear scaling factor
      real(fp) :: beta = 1.0  ! Exponential scaling factor on source parameter
      real(fp) :: drylimit_factor = 1.0  ! Dry Limit factor modifying the Fecan dry limit following Zender 2003
      real(fp) :: drag_option = 1  ! Drag Partition Option: 1 - use input drag, 2 - Darmenova, 3 - Leung 2022, 4 - MB95
      real(fp) :: moist_option = 1 - fecan  ! Moisture parameterization: 1 - Fecan, 2 - shao, 3 - modified shao
      real(fp) :: distribution_option = 1  ! Dust Distribution option: 1 - Kok 2011, 2 - Meng 2022

      ! Required meteorological fields
      integer :: n_required_met_fields = 12
      character(len=32) :: required_met_fields(12)

   contains
      procedure, public :: init => init_fengsha_config
      procedure, public :: validate => validate_fengsha_config
      procedure, public :: finalize => finalize_fengsha_config
   end type DustSchemeFENGSHAConfig

   !> State type for fengsha scheme
   type :: DustSchemeFENGSHAState

      ! Scheme working arrays
      real(fp), allocatable :: work_array_1(:,:)
      real(fp), allocatable :: work_array_2(:,:)

      ! Scheme-specific diagnostic arrays
      real(fp), allocatable :: dust_horizontal_flux(:)  ! Total horizontal flux - Q
      real(fp), allocatable :: dust_moisture_correction(:)  ! Moisture Correction - H
      real(fp), allocatable :: dust_effective_threshold(:)  ! Effective Dust threshold friction velocity: u_thres * H / R

   contains
      procedure, public :: init => init_fengsha_state
      procedure, public :: allocate_arrays => allocate_fengsha_state_arrays
      procedure, public :: reset => reset_fengsha_state
      procedure, public :: finalize => finalize_fengsha_state
   end type DustSchemeFENGSHAState

   !> Configuration type for ginoux scheme
   type :: DustSchemeGINOUXConfig

      ! Scheme metadata
      character(len=64) :: scheme_name = 'ginoux'
      character(len=256) :: description = 'Ginoux dust emission scheme'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Scheme parameters
      real(fp) :: Ch_DU = [0.1, 0.1, 0.1, 0.1, 0.1]  ! Dust tuning coefficient per species 

      ! Required meteorological fields
      integer :: n_required_met_fields = 5
      character(len=32) :: required_met_fields(5)

   contains
      procedure, public :: init => init_ginoux_config
      procedure, public :: validate => validate_ginoux_config
      procedure, public :: finalize => finalize_ginoux_config
   end type DustSchemeGINOUXConfig

   !> State type for ginoux scheme
   type :: DustSchemeGINOUXState

      ! Scheme working arrays
      real(fp), allocatable :: work_array_1(:,:)
      real(fp), allocatable :: work_array_2(:,:)

      ! Scheme-specific diagnostic arrays
      real(fp), allocatable :: dust_horizontal_flux(:)  ! Total horizontal flux - Q
      real(fp), allocatable :: dust_moisture_correction(:)  ! Moisture Correction - H
      real(fp), allocatable :: dust_effective_threshold(:)  ! Effective Dust threshold friction velocity: u_thres * H / R

   contains
      procedure, public :: init => init_ginoux_state
      procedure, public :: allocate_arrays => allocate_ginoux_state_arrays
      procedure, public :: reset => reset_ginoux_state
      procedure, public :: finalize => finalize_ginoux_state
   end type DustSchemeGINOUXState


contains

   !> Initialize dust configuration
   subroutine init_dust_config(this, error_handler)
      class(DustConfig), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! Set default species names



   end subroutine init_dust_config

   !> Validate dust configuration
   subroutine validate_dust_config(this, error_handler)
      class(DustConfig), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      character(len=256) :: error_msg

      ! Validate time step bounds
      if (this%dt_min <= 0.0_fp) then
         call error_handler%set_error(ERROR_CONFIG, &
            "Minimum time step must be positive")
         return
      end if

      if (this%dt_max < this%dt_min) then
         call error_handler%set_error(ERROR_CONFIG, &
            "Maximum time step must be >= minimum time step")
         return
      end if

      ! Validate active scheme
      if (trim(this%active_scheme) /= 'fengsha' .and. &
          trim(this%active_scheme) /= 'ginoux' .and. &
          .true.) then
         write(error_msg, '(A)') "Invalid scheme: " // trim(this%active_scheme)
         call error_handler%set_error(ERROR_CONFIG, error_msg)
         return
      end if

   end subroutine validate_dust_config

   !> Print configuration summary
   subroutine print_dust_config_summary(this)
      class(DustConfig), intent(in) :: this

      write(*, '(A)') "=== Dust Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%active_scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_dust_config_summary

   !> Finalize dust configuration
   subroutine finalize_dust_config(this)
      class(DustConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_dust_config

   !> Initialize dust state
   subroutine init_dust_state(this, n_columns, n_levels, error_handler)
      class(DustState), intent(inout) :: this
      integer, intent(in) :: n_columns
      integer, intent(in) :: n_levels
      type(ErrorHandler), intent(inout) :: error_handler

      this%n_columns = n_columns
      this%n_levels = n_levels
      this%current_time = 0.0_fp
      this%last_update_time = 0.0_fp

      call this%allocate_arrays(error_handler)
      if (error_handler%has_error()) return

      this%is_initialized = .true.

   end subroutine init_dust_state

   !> Allocate state arrays
   subroutine allocate_dust_state_arrays(this, error_handler)
      class(DustState), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      integer :: alloc_stat

      ! Allocate emission arrays
      allocate(this%emission_rates(1, this%n_columns), &
               stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate emission_rates array")
         return
      end if

      allocate(this%total_emissions(this%n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate total_emissions array")
         return
      end if

      ! Allocate diagnostic arrays
      allocate(this%total_dust_emission(this%n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate total_dust_emission array")
         return
      end if

   end subroutine allocate_dust_state_arrays

   !> Reset state arrays to zero
   subroutine reset_dust_state(this)
      class(DustState), intent(inout) :: this

      if (allocated(this%emission_rates)) this%emission_rates = 0.0_fp
      if (allocated(this%total_emissions)) this%total_emissions = 0.0_fp

      if (allocated(this%total_dust_emission)) this%total_dust_emission = 0.0_fp

   end subroutine reset_dust_state

   !> Finalize dust state
   subroutine finalize_dust_state(this)
      class(DustState), intent(inout) :: this

      if (allocated(this%emission_rates)) deallocate(this%emission_rates)
      if (allocated(this%total_emissions)) deallocate(this%total_emissions)

      if (allocated(this%total_dust_emission)) deallocate(this%total_dust_emission)

      this%is_initialized = .false.

   end subroutine finalize_dust_state

   !> Initialize fengsha scheme configuration
   subroutine init_fengsha_config(this, parent_config, error_handler)
      class(DustSchemeFENGSHAConfig), intent(inout) :: this
      type(DustConfig), intent(in) :: parent_config
      type(ErrorHandler), intent(inout) :: error_handler

      ! Set required meteorological fields
      this%required_met_fields(1) = 'IsLand'
      this%required_met_fields(2) = 'USTAR'
      this%required_met_fields(3) = 'LWI'
      this%required_met_fields(4) = 'GVF'
      this%required_met_fields(5) = 'LAI'
      this%required_met_fields(6) = 'FROCEAN'
      this%required_met_fields(7) = 'CLAYFRAC'
      this%required_met_fields(8) = 'SANDFRAC'
      this%required_met_fields(9) = 'FRSNO'
      this%required_met_fields(10) = 'RDRAG'
      this%required_met_fields(11) = 'SSM'
      this%required_met_fields(12) = 'USTAR_THRESHOLD'

      ! Initialize scheme-specific parameters
      ! TODO: Load from configuration file or use defaults

   end subroutine init_fengsha_config

   !> Validate fengsha scheme configuration
   subroutine validate_fengsha_config(this, error_handler)
      class(DustSchemeFENGSHAConfig), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_fengsha_config

   !> Finalize fengsha scheme configuration
   subroutine finalize_fengsha_config(this)
      class(DustSchemeFENGSHAConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_fengsha_config

   !> Initialize fengsha scheme state
   subroutine init_fengsha_state(this, parent_state, error_handler)
      class(DustSchemeFENGSHAState), intent(inout) :: this
      type(DustState), intent(in) :: parent_state
      type(ErrorHandler), intent(inout) :: error_handler

      call this%allocate_arrays(parent_state%n_columns, error_handler)

   end subroutine init_fengsha_state

   !> Allocate fengsha scheme state arrays
   subroutine allocate_fengsha_state_arrays(this, n_columns, error_handler)
      class(DustSchemeFENGSHAState), intent(inout) :: this
      integer, intent(in) :: n_columns
      type(ErrorHandler), intent(inout) :: error_handler

      integer :: alloc_stat

      ! Allocate working arrays
      allocate(this%work_array_1(1, n_columns), &
               stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate fengsha work_array_1")
         return
      end if

      allocate(this%work_array_2(1, n_columns), &
               stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate fengsha work_array_2")
         return
      end if

      ! Allocate diagnostic arrays
      allocate(this%dust_horizontal_flux(n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate fengsha dust_horizontal_flux array")
         return
      end if
      allocate(this%dust_moisture_correction(n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate fengsha dust_moisture_correction array")
         return
      end if
      allocate(this%dust_effective_threshold(n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate fengsha dust_effective_threshold array")
         return
      end if

   end subroutine allocate_fengsha_state_arrays

   !> Reset fengsha scheme state
   subroutine reset_fengsha_state(this)
      class(DustSchemeFENGSHAState), intent(inout) :: this

      if (allocated(this%work_array_1)) this%work_array_1 = 0.0_fp
      if (allocated(this%work_array_2)) this%work_array_2 = 0.0_fp

      if (allocated(this%dust_horizontal_flux)) this%dust_horizontal_flux = 0.0_fp
      if (allocated(this%dust_moisture_correction)) this%dust_moisture_correction = 0.0_fp
      if (allocated(this%dust_effective_threshold)) this%dust_effective_threshold = 0.0_fp

   end subroutine reset_fengsha_state

   !> Finalize fengsha scheme state
   subroutine finalize_fengsha_state(this)
      class(DustSchemeFENGSHAState), intent(inout) :: this

      if (allocated(this%work_array_1)) deallocate(this%work_array_1)
      if (allocated(this%work_array_2)) deallocate(this%work_array_2)

      if (allocated(this%dust_horizontal_flux)) deallocate(this%dust_horizontal_flux)
      if (allocated(this%dust_moisture_correction)) deallocate(this%dust_moisture_correction)
      if (allocated(this%dust_effective_threshold)) deallocate(this%dust_effective_threshold)

   end subroutine finalize_fengsha_state

   !> Initialize ginoux scheme configuration
   subroutine init_ginoux_config(this, parent_config, error_handler)
      class(DustSchemeGINOUXConfig), intent(inout) :: this
      type(DustConfig), intent(in) :: parent_config
      type(ErrorHandler), intent(inout) :: error_handler

      ! Set required meteorological fields
      this%required_met_fields(1) = 'FRLAKE'
      this%required_met_fields(2) = 'GWETTOP'
      this%required_met_fields(3) = 'U10M'
      this%required_met_fields(4) = 'V10M'
      this%required_met_fields(5) = 'SSM'

      ! Initialize scheme-specific parameters
      ! TODO: Load from configuration file or use defaults

   end subroutine init_ginoux_config

   !> Validate ginoux scheme configuration
   subroutine validate_ginoux_config(this, error_handler)
      class(DustSchemeGINOUXConfig), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_ginoux_config

   !> Finalize ginoux scheme configuration
   subroutine finalize_ginoux_config(this)
      class(DustSchemeGINOUXConfig), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_ginoux_config

   !> Initialize ginoux scheme state
   subroutine init_ginoux_state(this, parent_state, error_handler)
      class(DustSchemeGINOUXState), intent(inout) :: this
      type(DustState), intent(in) :: parent_state
      type(ErrorHandler), intent(inout) :: error_handler

      call this%allocate_arrays(parent_state%n_columns, error_handler)

   end subroutine init_ginoux_state

   !> Allocate ginoux scheme state arrays
   subroutine allocate_ginoux_state_arrays(this, n_columns, error_handler)
      class(DustSchemeGINOUXState), intent(inout) :: this
      integer, intent(in) :: n_columns
      type(ErrorHandler), intent(inout) :: error_handler

      integer :: alloc_stat

      ! Allocate working arrays
      allocate(this%work_array_1(1, n_columns), &
               stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate ginoux work_array_1")
         return
      end if

      allocate(this%work_array_2(1, n_columns), &
               stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate ginoux work_array_2")
         return
      end if

      ! Allocate diagnostic arrays
      allocate(this%dust_horizontal_flux(n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate ginoux dust_horizontal_flux array")
         return
      end if
      allocate(this%dust_moisture_correction(n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate ginoux dust_moisture_correction array")
         return
      end if
      allocate(this%dust_effective_threshold(n_columns), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate ginoux dust_effective_threshold array")
         return
      end if

   end subroutine allocate_ginoux_state_arrays

   !> Reset ginoux scheme state
   subroutine reset_ginoux_state(this)
      class(DustSchemeGINOUXState), intent(inout) :: this

      if (allocated(this%work_array_1)) this%work_array_1 = 0.0_fp
      if (allocated(this%work_array_2)) this%work_array_2 = 0.0_fp

      if (allocated(this%dust_horizontal_flux)) this%dust_horizontal_flux = 0.0_fp
      if (allocated(this%dust_moisture_correction)) this%dust_moisture_correction = 0.0_fp
      if (allocated(this%dust_effective_threshold)) this%dust_effective_threshold = 0.0_fp

   end subroutine reset_ginoux_state

   !> Finalize ginoux scheme state
   subroutine finalize_ginoux_state(this)
      class(DustSchemeGINOUXState), intent(inout) :: this

      if (allocated(this%work_array_1)) deallocate(this%work_array_1)
      if (allocated(this%work_array_2)) deallocate(this%work_array_2)

      if (allocated(this%dust_horizontal_flux)) deallocate(this%dust_horizontal_flux)
      if (allocated(this%dust_moisture_correction)) deallocate(this%dust_moisture_correction)
      if (allocated(this%dust_effective_threshold)) deallocate(this%dust_effective_threshold)

   end subroutine finalize_ginoux_state


   !> Convert integer to string (utility function)
   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

end module DustCommon_Mod