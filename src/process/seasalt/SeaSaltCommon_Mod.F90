!> \file SeaSaltCommon_Mod.F90
!! \brief Common types and utilities for seasalt process
!!
!! This module defines the configuration types used by the
!! seasalt process and its schemes.
!!
!! Generated on: 2025-08-06T23:49:33.280423
!! Author: Barry Baker & Wei Li
!! Version: 1.0.0

module SeaSaltCommon_Mod

   use iso_fortran_env, only: fp => real64
   use precision_mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use ConfigManager_Mod, only: ConfigDataType  ! Add ConfigManager integration
   use StateManager_Mod, only: StateManagerType  ! Add StateManager integration

   implicit none
   private

   ! Export types
   public :: SeaSaltProcessConfig  ! New unified process config
   public :: SeaSaltConfig
   public :: SeaSaltSchemeGONG97Config
   public :: SeaSaltSchemeGONG03Config
   public :: SeaSaltSchemeGEOS12Config

   ! Export utility functions
   public :: int_to_string

   !> Unified process configuration type that bridges ConfigManager and process-specific configs
   !! This is the main configuration type that ProcessInterface should use
   type :: SeaSaltProcessConfig
      
      ! Process metadata
      character(len=64) :: process_name = 'seasalt'
      character(len=16) :: process_version = '1.0.0'
      logical :: is_active = .true.
      
      ! Process-specific configuration (delegate to SeaSaltConfig)
      type(SeaSaltConfig) :: seasalt_config
      
      ! Scheme configurations
      type(SeaSaltSchemeGONG97Config) :: gong97_config
      type(SeaSaltSchemeGONG03Config) :: gong03_config
      type(SeaSaltSchemeGEOS12Config) :: geos12_config
      
   contains
      procedure, public :: load_from_config => seasalt_process_load_config
      procedure, public :: load_species_from_chem_state => load_species_from_chem_state
      procedure, public :: validate => seasalt_process_validate
      procedure, public :: finalize => seasalt_process_finalize
      procedure, public :: get_active_scheme_config => get_active_scheme_config
      procedure, private :: load_gong97_config
      procedure, private :: load_gong03_config
      procedure, private :: load_geos12_config
   end type SeaSaltProcessConfig

   !> Main configuration type for seasalt process
   type :: SeaSaltConfig

      ! Process settings
      character(len=32) :: scheme = 'gong97'
      logical :: is_active = .true.
      real(fp) :: dt_min = 1.0_fp     ! Minimum time step (seconds)
      real(fp) :: dt_max = 3600.0_fp  ! Maximum time step (seconds)

      ! Species configuration
      integer :: n_species = 0
      character(len=32), allocatable :: species_names(:)
      integer, allocatable :: species_indices(:)  ! Indices of seasalt species in ChemState

      ! Diagnostic configuration
      logical :: output_diagnostics = .true.
      real(fp) :: diagnostic_frequency = 3600.0_fp  ! Output frequency (seconds)

   contains
      procedure, public :: validate => validate_seasalt_config
      procedure, public :: finalize => finalize_seasalt_config
      procedure, public :: print_summary => print_seasalt_config_summary
   end type SeaSaltConfig

   !> Configuration type for gong97 scheme
   type :: SeaSaltSchemeGONG97Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gong97'
      character(len=256) :: description = 'Gong 1997 sea salt emission scheme'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 2
      character(len=32) :: required_met_fields(2)

   contains
      procedure, public :: validate => validate_gong97_config
      procedure, public :: finalize => finalize_gong97_config
   end type SeaSaltSchemeGONG97Config

   ! gong97 scheme uses local variables only - no persistent state type needed

   !> Configuration type for gong03 scheme
   type :: SeaSaltSchemeGONG03Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'gong03'
      character(len=256) :: description = 'Gong 2003 sea salt emission scheme with improved sub- and super-micron treatment'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor
      logical :: weibull_flag = .false.  ! Apply Weibull distribution for particle size

      ! Required meteorological fields
      integer :: n_required_met_fields = 2
      character(len=32) :: required_met_fields(2)

   contains
      procedure, public :: validate => validate_gong03_config
      procedure, public :: finalize => finalize_gong03_config
   end type SeaSaltSchemeGONG03Config

   ! gong03 scheme uses local variables only - no persistent state type needed

   !> Configuration type for geos12 scheme
   type :: SeaSaltSchemeGEOS12Config

      ! Scheme metadata
      character(len=64) :: scheme_name = 'geos12'
      character(len=256) :: description = 'GEOS-Chem 2012 sea salt emission scheme with observational constraints'
      character(len=64) :: author = 'Barry Baker'
      character(len=16) :: algorithm_type = 'explicit'

      ! Scheme parameters
      real(fp) :: scale_factor = 1.0  ! Emission scale factor

      ! Required meteorological fields
      integer :: n_required_met_fields = 1
      character(len=32) :: required_met_fields(1)

   contains
      procedure, public :: validate => validate_geos12_config
      procedure, public :: finalize => finalize_geos12_config
   end type SeaSaltSchemeGEOS12Config

   ! geos12 scheme uses local variables only - no persistent state type needed


contains

   !> Validate seasalt configuration
   subroutine validate_seasalt_config(this, error_handler)
      class(SeaSaltConfig), intent(inout) :: this
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
      if (trim(this%scheme) /= 'gong97' .and. &
          trim(this%scheme) /= 'gong03' .and. &
          trim(this%scheme) /= 'geos12' .and. &
          .true.) then
         write(error_msg, '(A)') "Invalid scheme: " // trim(this%scheme)
         call error_handler%set_error(ERROR_CONFIG, error_msg)
         return
      end if

   end subroutine validate_seasalt_config

   !> Print configuration summary
   subroutine print_seasalt_config_summary(this)
      class(SeaSaltConfig), intent(in) :: this

      write(*, '(A)') "=== SeaSalt Process Configuration ==="
      write(*, '(A,A)') "  Active scheme: ", trim(this%scheme)
      write(*, '(A,I0)') "  Number of species: ", this%n_species
      write(*, '(A,F0.1,A)') "  Minimum time step: ", this%dt_min, " s"
      write(*, '(A,F0.1,A)') "  Maximum time step: ", this%dt_max, " s"
      write(*, '(A,L1)') "  Output diagnostics: ", this%output_diagnostics
      write(*, '(A)') "============================================="

   end subroutine print_seasalt_config_summary

      !> Finalize seasalt configuration
   subroutine finalize_seasalt_config(this)
      class(SeaSaltConfig), intent(inout) :: this

      ! Deallocate species names array
      if (allocated(this%species_names)) then
         deallocate(this%species_names)
      end if

      ! Deallocate species indices array
      if (allocated(this%species_indices)) then
         deallocate(this%species_indices)
      end if

   end subroutine finalize_seasalt_config

   !> Validate gong97 scheme configuration
   subroutine validate_gong97_config(this, error_handler)
      class(SeaSaltSchemeGONG97Config), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gong97_config

   !> Finalize gong97 scheme configuration
   subroutine finalize_gong97_config(this)
      class(SeaSaltSchemeGONG97Config), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gong97_config


   !> Validate gong03 scheme configuration
   subroutine validate_gong03_config(this, error_handler)
      class(SeaSaltSchemeGONG03Config), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_gong03_config

   !> Finalize gong03 scheme configuration
   subroutine finalize_gong03_config(this)
      class(SeaSaltSchemeGONG03Config), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_gong03_config


   !> Validate geos12 scheme configuration
   subroutine validate_geos12_config(this, error_handler)
      class(SeaSaltSchemeGEOS12Config), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Add scheme-specific validation

   end subroutine validate_geos12_config

   !> Finalize geos12 scheme configuration
   subroutine finalize_geos12_config(this)
      class(SeaSaltSchemeGEOS12Config), intent(inout) :: this

      ! Nothing to deallocate for basic configuration

   end subroutine finalize_geos12_config



   !> Convert integer to string (utility function)
   function int_to_string(int_val) result(str_val)
      integer, intent(in) :: int_val
      character(len=32) :: str_val

      write(str_val, '(I0)') int_val
      str_val = adjustl(str_val)

   end function int_to_string

   !> Load unified process configuration from ConfigManager
   !! This is the main function that ProcessInterface.parse_process_config should call
   !! Process reads its configuration directly from the master YAML via ConfigManager
   subroutine seasalt_process_load_config(this, config_data, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigDataType), intent(in) :: config_data
      type(ErrorManagerType), intent(inout) :: error_handler

      character(len=256) :: scheme_name
      integer :: ierr

      ! Process reads directly from master YAML structure: processes.seasalt.*
      ! ConfigManager provides generic YAML access, process handles its own configuration
      
      ! Load process metadata
      call config_data%get_string("processes/seasalt/name", this%process_name, ierr)
      if (ierr /= CC_SUCCESS) this%process_name = "seasalt"  ! default

      call config_data%get_string("processes/seasalt/version", this%process_version, ierr)
      if (ierr /= CC_SUCCESS) this%process_version = "1.0.0"  ! default

      call config_data%get_logical("processes/seasalt/activate", this%is_active, ierr)
      if (ierr /= CC_SUCCESS) this%is_active = .true.  ! default

      ! Load process-specific configuration directly from master YAML
      call config_data%get_string("processes/seasalt/scheme", this%seasalt_config%scheme, ierr)
      if (ierr /= CC_SUCCESS) then
         call error_handler%set_error(CC_FAILURE, &
            "Missing required 'scheme' in processes/seasalt configuration")
         return
      end if

      ! Species configuration is loaded from ChemState in load_species_from_chem_state
      ! The species come from the master species YAML file (CATChem_species.yml)
      ! and are filtered by is_seasalt property


      ! Load scheme-specific configuration from master YAML
      scheme_name = trim(this%seasalt_config%scheme)
      select case (scheme_name)
      case ('gong97')
         call this%load_gong97_config(config_data, error_handler)
      case ('gong03')
         call this%load_gong03_config(config_data, error_handler)
      case ('geos12')
         call this%load_geos12_config(config_data, error_handler)
      case default
         call error_handler%set_error(CC_FAILURE, &
            "Unknown seasalt scheme: " // trim(scheme_name))
         return
      end select

   end subroutine seasalt_process_load_config

   !> Load species from ChemState based on SeaSaltIndex
   !! This should be called after process configuration is loaded and ChemState is available
   subroutine load_species_from_chem_state(this, chem_state, error_handler)
      use ChemState_Mod, only: ChemStateType
      
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ChemStateType), pointer, intent(in) :: chem_state
      type(ErrorManagerType), intent(inout) :: error_handler
      
      integer :: i
      
      if (.not. associated(chem_state)) then
         call error_handler%set_error(CC_FAILURE, &
            "ChemState not associated in load_species_from_chem_state")
         return
      end if
      
      ! Get the count of seasalt species from ChemState
      this%seasalt_config%n_species = chem_state%nSpeciesSeaSalt
      
      if (this%seasalt_config%n_species <= 0) then
         call error_handler%set_error(CC_FAILURE, &
            "No seasalt species found in ChemState")
         return
      end if
      
      ! Check if SeaSaltIndex is allocated and has correct size
      if (.not. allocated(chem_state%SeaSaltIndex)) then
         call error_handler%set_error(CC_FAILURE, &
            "SeaSaltIndex not allocated in ChemState")
         return
      end if
      
      if (size(chem_state%SeaSaltIndex) < this%seasalt_config%n_species) then
         call error_handler%set_error(CC_FAILURE, &
            "SeaSaltIndex size inconsistent with nSpeciesSeaSalt")
         return
      end if
      
      ! Deallocate existing arrays if allocated
      if (allocated(this%seasalt_config%species_names)) then
         deallocate(this%seasalt_config%species_names)
      end if
      if (allocated(this%seasalt_config%species_indices)) then
         deallocate(this%seasalt_config%species_indices)
      end if
      
      ! Allocate arrays
      allocate(this%seasalt_config%species_names(this%seasalt_config%n_species))
      allocate(this%seasalt_config%species_indices(this%seasalt_config%n_species))
      
      ! Copy species indices directly from ChemState
      this%seasalt_config%species_indices(1:this%seasalt_config%n_species) = &
         chem_state%SeaSaltIndex(1:this%seasalt_config%n_species)
      
      ! Get species names using the indices
      do i = 1, this%seasalt_config%n_species
         if (this%seasalt_config%species_indices(i) > 0 .and. &
             this%seasalt_config%species_indices(i) <= size(chem_state%SpeciesNames)) then
            this%seasalt_config%species_names(i) = &
               trim(chem_state%SpeciesNames(this%seasalt_config%species_indices(i)))
         else
            call error_handler%set_error(CC_FAILURE, &
               "Invalid species index in SeaSaltIndex")
            return
         end if
      end do
      
   end subroutine load_species_from_chem_state

   !> Load gong97 scheme configuration from master YAML
   subroutine load_gong97_config(this, config_data, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigDataType), intent(in) :: config_data
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr

      ! Load scheme parameters directly from processes/seasalt/gong97/* in master YAML
      call config_data%get_real("processes/seasalt/gong97/scale_factor", &
           this%gong97_config%scale_factor, ierr)
      if (ierr /= CC_SUCCESS) this%gong97_config%scale_factor = 1.0
      call config_data%get_logical("processes/seasalt/gong97/weibull_flag", &
           this%gong97_config%weibull_flag, ierr)
      if (ierr /= CC_SUCCESS) this%gong97_config%weibull_flag = .false.


   end subroutine load_gong97_config

   !> Load gong03 scheme configuration from master YAML
   subroutine load_gong03_config(this, config_data, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigDataType), intent(in) :: config_data
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr

      ! Load scheme parameters directly from processes/seasalt/gong03/* in master YAML
      call config_data%get_real("processes/seasalt/gong03/scale_factor", &
           this%gong03_config%scale_factor, ierr)
      if (ierr /= CC_SUCCESS) this%gong03_config%scale_factor = 1.0
      call config_data%get_logical("processes/seasalt/gong03/weibull_flag", &
           this%gong03_config%weibull_flag, ierr)
      if (ierr /= CC_SUCCESS) this%gong03_config%weibull_flag = .false.


   end subroutine load_gong03_config

   !> Load geos12 scheme configuration from master YAML
   subroutine load_geos12_config(this, config_data, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(ConfigDataType), intent(in) :: config_data
      type(ErrorManagerType), intent(inout) :: error_handler

      integer :: ierr

      ! Load scheme parameters directly from processes/seasalt/geos12/* in master YAML
      call config_data%get_real("processes/seasalt/geos12/scale_factor", &
           this%geos12_config%scale_factor, ierr)
      if (ierr /= CC_SUCCESS) this%geos12_config%scale_factor = 1.0


   end subroutine load_geos12_config


   !> Validate unified process configuration
   subroutine seasalt_process_validate(this, state_manager, error_handler)
      class(SeaSaltProcessConfig), intent(inout) :: this
      type(StateManagerType), intent(in) :: state_manager
      type(ErrorManagerType), intent(inout) :: error_handler

      ! Validate main config
      call this%seasalt_config%validate(error_handler)
      if (error_handler%has_error()) return

      ! Validate scheme-specific config
      select case (trim(this%seasalt_config%scheme))
      case ('gong97')
         call this%gong97_config%validate(error_handler)
      case ('gong03')
         call this%gong03_config%validate(error_handler)
      case ('geos12')
         call this%geos12_config%validate(error_handler)
      end select

   end subroutine seasalt_process_validate

   !> Finalize unified process configuration
   subroutine seasalt_process_finalize(this)
      class(SeaSaltProcessConfig), intent(inout) :: this

      call this%seasalt_config%finalize()
      call this%gong97_config%finalize()
      call this%gong03_config%finalize()
      call this%geos12_config%finalize()

   end subroutine seasalt_process_finalize

   !> Get active scheme configuration (polymorphic return)
   function get_active_scheme_config(this) result(scheme_config)
      class(SeaSaltProcessConfig), intent(in) :: this
      class(*), allocatable :: scheme_config

      select case (trim(this%seasalt_config%scheme))
      case ('gong97')
         allocate(scheme_config, source=this%gong97_config)
      case ('gong03')
         allocate(scheme_config, source=this%gong03_config)
      case ('geos12')
         allocate(scheme_config, source=this%geos12_config)
      case default
         ! Return null
      end select

   end function get_active_scheme_config

end module SeaSaltCommon_Mod