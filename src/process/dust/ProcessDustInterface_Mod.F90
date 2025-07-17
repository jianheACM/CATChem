!> \file ProcessDustInterface_Mod.F90
!! \brief Process for computing windblown dust emissions
!!
!! This module provides the main interface for the dust process
!! following the CATChem Process Infrastructure Guide.
!!
!! Generated on: 2025-07-09T12:43:17.965853
!! Author: Barry Baker
!! Version: 1.0.0

module ProcessDustInterface_Mod

   use iso_fortran_env, only: fp => real64, error_unit
   use StateManager_Mod
   use ProcessInterface_Mod
   use DiagnosticManager_Mod
   use ErrorHandler_Mod
   use FieldManager_Mod
   use DustCommon_Mod
   use DustScheme_FENGSHA_Mod
   use DustScheme_GINOUX_Mod

   implicit none
   private

   public :: ProcessDustInterface

   !> Main dust process interface type
   type, extends(ProcessInterface) :: ProcessDustInterface
      private

      ! Process metadata
      character(len=64) :: process_name = 'dust'
      character(len=64) :: process_version = '1.0.0'
      character(len=256) :: process_description = 'Process for computing windblown dust emissions'

      ! Configuration
      type(DustConfig), pointer :: config => null()
      type(DustState), pointer :: state => null()

      ! Scheme selection and configuration
      character(len=32) :: active_scheme = ''
      class(DustSchemeFENGSHAConfig), allocatable :: config_fengsha
      class(DustSchemeFENGSHAState), allocatable :: state_fengsha
      class(DustSchemeGINOUXConfig), allocatable :: config_ginoux
      class(DustSchemeGINOUXState), allocatable :: state_ginoux

      ! Species management
      integer :: n_species = 0
      character(len=32), allocatable :: species_names(:)
      integer, allocatable :: species_indices(:)


      ! Diagnostic indices
      integer :: diag_total_dust_emission_idx = -1
      integer :: diag_fengsha_dust_horizontal_flux_idx = -1
      integer :: diag_fengsha_dust_moisture_correction_idx = -1
      integer :: diag_fengsha_dust_effective_threshold_idx = -1
      integer :: diag_ginoux_dust_horizontal_flux_idx = -1
      integer :: diag_ginoux_dust_moisture_correction_idx = -1
      integer :: diag_ginoux_dust_effective_threshold_idx = -1

      ! Emission tendencies for diagnostic tracking
      real(fp), allocatable :: emission_tendencies(:,:)  ! (n_species, n_columns)

   contains

      ! Required interface methods
      procedure, public :: init => init_dust_process
      procedure, public :: run => run_dust_process
      procedure, public :: finalize => finalize_dust_process

      ! Configuration methods
      procedure, public :: load_config => load_dust_config
      procedure, public :: validate_config => validate_dust_config

      ! Species management methods
      procedure, public :: get_species_list => get_dust_species_list
      procedure, public :: set_species_mapping => set_dust_species_mapping

      ! Diagnostic methods
      procedure, public :: register_diagnostics => register_dust_diagnostics
      procedure, public :: update_diagnostics => update_dust_diagnostics

      ! Required meteorological fields
      procedure, public :: get_required_met_fields => get_dust_required_met_fields

      ! Emission-specific methods
      procedure, public :: apply_emission_tendencies => apply_dust_emission_tendencies

   end type ProcessDustInterface

contains

   !> Initialize the dust process
   subroutine init_dust_process(this, state_manager, config_data, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      character(len=*), intent(in) :: config_data
      type(ErrorHandler), intent(inout) :: error_handler

      character(len=256) :: error_msg
      integer :: i, rc

      ! Initialize parent
      call this%ProcessInterface%init(state_manager, config_data, error_handler)
      if (error_handler%has_error()) return

      ! Load and validate configuration
      call this%load_config(config_data, error_handler)
      if (error_handler%has_error()) return

      call this%validate_config(error_handler)
      if (error_handler%has_error()) return

      ! Set up species mapping
      call this%set_species_mapping(state_manager, error_handler)
      if (error_handler%has_error()) return


      ! Register diagnostics
      call this%register_diagnostics(state_manager, error_handler)
      if (error_handler%has_error()) return

      ! Initialize scheme-specific configurations and states
      if (trim(this%active_scheme) == 'fengsha') then
         allocate(this%config_fengsha)
         allocate(this%state_fengsha)
         call this%config_fengsha%init(this%config, error_handler)
         if (error_handler%has_error()) return
         call this%state_fengsha%init(this%state, error_handler)
         if (error_handler%has_error()) return
      end if
      if (trim(this%active_scheme) == 'ginoux') then
         allocate(this%config_ginoux)
         allocate(this%state_ginoux)
         call this%config_ginoux%init(this%config, error_handler)
         if (error_handler%has_error()) return
         call this%state_ginoux%init(this%state, error_handler)
         if (error_handler%has_error()) return
      end if

      ! Allocate emission tendencies
      allocate(this%emission_tendencies(this%n_species, state_manager%get_n_columns()))
      this%emission_tendencies = 0.0_fp

      this%is_initialized = .true.

   end subroutine init_dust_process

   !> Run the dust process
   subroutine run_dust_process(this, state_manager, dt, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      real(fp), intent(in) :: dt
      type(ErrorHandler), intent(inout) :: error_handler

      character(len=256) :: error_msg
      integer :: i_col, n_columns

      if (.not. this%is_initialized) then
         call error_handler%set_error(ERROR_PROCESS, &
            "Dust process not initialized")
         return
      end if

      n_columns = state_manager%get_n_columns()

      ! Reset emission tendencies
      this%emission_tendencies = 0.0_fp

      ! Process each column
      do i_col = 1, n_columns

         ! Run active scheme
         if (trim(this%active_scheme) == 'fengsha') then
            call this%run_fengsha_scheme(state_manager, i_col, dt, error_handler)
            if (error_handler%has_error()) return
         else if (trim(this%active_scheme) == 'ginoux') then
            call this%run_ginoux_scheme(state_manager, i_col, dt, error_handler)
            if (error_handler%has_error()) return
         else
            call error_handler%set_error(ERROR_PROCESS, &
               "Unknown dust scheme: " // trim(this%active_scheme))
            return
         end if

      end do

      ! Apply emission tendencies to chemical state
      call this%apply_emission_tendencies(state_manager, dt, error_handler)
      if (error_handler%has_error()) return

      ! Update diagnostics
      call this%update_diagnostics(state_manager, error_handler)
      if (error_handler%has_error()) return

   end subroutine run_dust_process

   !> Run the fengsha scheme
   subroutine run_fengsha_scheme(this, state_manager, i_col, dt, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      integer, intent(in) :: i_col
      real(fp), intent(in) :: dt
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Implement fengsha scheme calculation
      ! This is a template - implement actual algorithm here

      ! Example structure:
      ! 1. Get required meteorological fields
      ! 2. Get current chemical concentrations
      ! 3. Calculate process rates/tendencies
      ! 4. Update emission tendencies or chemical state
      ! 5. Update scheme-specific diagnostics

   end subroutine run_fengsha_scheme

   !> Run the ginoux scheme
   subroutine run_ginoux_scheme(this, state_manager, i_col, dt, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      integer, intent(in) :: i_col
      real(fp), intent(in) :: dt
      type(ErrorHandler), intent(inout) :: error_handler

      ! TODO: Implement ginoux scheme calculation
      ! This is a template - implement actual algorithm here

      ! Example structure:
      ! 1. Get required meteorological fields
      ! 2. Get current chemical concentrations
      ! 3. Calculate process rates/tendencies
      ! 4. Update emission tendencies or chemical state
      ! 5. Update scheme-specific diagnostics

   end subroutine run_ginoux_scheme


   !> Load configuration from input data
   subroutine load_dust_config(this, config_data, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      character(len=*), intent(in) :: config_data
      type(ErrorHandler), intent(inout) :: error_handler

      ! Allocate configuration
      allocate(this%config)
      allocate(this%state)

      ! TODO: Parse configuration from config_data (YAML, namelist, etc.)
      ! For now, use defaults
      this%config%active_scheme = ''
      this%active_scheme = this%config%active_scheme

      ! Set species information

   end subroutine load_dust_config

   !> Validate configuration
   subroutine validate_dust_config(this, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      character(len=256) :: error_msg

      ! Validate scheme selection
      if (trim(this%active_scheme) /= 'fengsha' .and. &
          trim(this%active_scheme) /= 'ginoux' .and. &
          .true.) then
         write(error_msg, '(A)') "Invalid dust scheme: " // trim(this%active_scheme)
         call error_handler%set_error(ERROR_CONFIG, error_msg)
         return
      end if

      ! Additional validation as needed

   end subroutine validate_dust_config

   !> Get list of species handled by this process
   function get_dust_species_list(this) result(species_names)
      class(ProcessDustInterface), intent(in) :: this
      character(len=32), allocatable :: species_names(:)

      allocate(species_names(this%n_species))
      species_names = this%species_names

   end function get_dust_species_list

   !> Set species mapping from chemical state
   subroutine set_dust_species_mapping(this, state_manager, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      type(ErrorHandler), intent(inout) :: error_handler

      integer :: i, idx
      character(len=256) :: error_msg

      ! Map species names to indices in chemical state
      do i = 1, this%n_species
         idx = state_manager%get_species_index(this%species_names(i))
         if (idx <= 0) then
            write(error_msg, '(A)') "Species not found in chemical state: " // &
               trim(this%species_names(i))
            call error_handler%set_error(ERROR_SPECIES, error_msg)
            return
         end if
         this%species_indices(i) = idx
      end do

   end subroutine set_dust_species_mapping

   !> Register diagnostic fields
   subroutine register_dust_diagnostics(this, state_manager, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      type(ErrorHandler), intent(inout) :: error_handler

      type(DiagnosticManager), pointer :: diag_mgr

      diag_mgr => state_manager%get_diagnostic_manager()

      ! Register common diagnostics
      this%diag_total_dust_emission_idx = diag_mgr%register_diagnostic( &
         name='dust_total_dust_emission', &
         description='Total dust emissions for all species', &
         units='kg/m2/s', &
         error_handler=error_handler)
      if (error_handler%has_error()) return

      ! Register emission tendency diagnostics for each species
      do i = 1, this%n_species
         write(name, '(A,A,A)') 'dust_emission_', trim(this%species_names(i)), '_tendency'
         write(description, '(A,A)') 'Emission tendency for ', trim(this%species_names(i))
         idx = diag_mgr%register_diagnostic( &
            name=trim(name), &
            description=trim(description), &
            units='kg/m2/s', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
      end do

      ! Register scheme-specific diagnostics
      if (trim(this%active_scheme) == 'fengsha') then
         this%diag_fengsha_dust_horizontal_flux_idx = diag_mgr%register_diagnostic( &
            name='dust_fengsha_dust_horizontal_flux', &
            description='Total horizontal flux - Q', &
            units='kg/m2/s', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
         this%diag_fengsha_dust_moisture_correction_idx = diag_mgr%register_diagnostic( &
            name='dust_fengsha_dust_moisture_correction', &
            description='Moisture Correction - H', &
            units='dimensionless', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
         this%diag_fengsha_dust_effective_threshold_idx = diag_mgr%register_diagnostic( &
            name='dust_fengsha_dust_effective_threshold', &
            description='Effective Dust threshold friction velocity: u_thres * H / R', &
            units='m/s', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
      end if
      if (trim(this%active_scheme) == 'ginoux') then
         this%diag_ginoux_dust_horizontal_flux_idx = diag_mgr%register_diagnostic( &
            name='dust_ginoux_dust_horizontal_flux', &
            description='Total horizontal flux - Q', &
            units='kg/m2/s', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
         this%diag_ginoux_dust_moisture_correction_idx = diag_mgr%register_diagnostic( &
            name='dust_ginoux_dust_moisture_correction', &
            description='Moisture Correction - H', &
            units='dimensionless', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
         this%diag_ginoux_dust_effective_threshold_idx = diag_mgr%register_diagnostic( &
            name='dust_ginoux_dust_effective_threshold', &
            description='Effective Dust threshold friction velocity: u_thres * H / R', &
            units='m/s', &
            error_handler=error_handler)
         if (error_handler%has_error()) return
      end if

   end subroutine register_dust_diagnostics

   !> Update diagnostic fields
   subroutine update_dust_diagnostics(this, state_manager, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      type(ErrorHandler), intent(inout) :: error_handler

      type(DiagnosticManager), pointer :: diag_mgr
      integer :: i

      diag_mgr => state_manager%get_diagnostic_manager()

      ! Update common diagnostics
      if (this%diag_total_dust_emission_idx > 0) then
         ! TODO: Calculate and update total_dust_emission
         ! call diag_mgr%update_diagnostic(this%diag_total_dust_emission_idx, values)
      end if

      ! Update emission tendency diagnostics
      do i = 1, this%n_species
         if (this%diag_emission_tendency_indices(i) > 0) then
            call diag_mgr%update_diagnostic( &
               this%diag_emission_tendency_indices(i), &
               this%emission_tendencies(i, :))
         end if
      end do

      ! Update scheme-specific diagnostics
      if (trim(this%active_scheme) == 'fengsha') then
         if (this%diag_fengsha_dust_horizontal_flux_idx > 0) then
            ! TODO: Calculate and update fengsha dust_horizontal_flux
            ! call diag_mgr%update_diagnostic(this%diag_fengsha_dust_horizontal_flux_idx, values)
         end if
         if (this%diag_fengsha_dust_moisture_correction_idx > 0) then
            ! TODO: Calculate and update fengsha dust_moisture_correction
            ! call diag_mgr%update_diagnostic(this%diag_fengsha_dust_moisture_correction_idx, values)
         end if
         if (this%diag_fengsha_dust_effective_threshold_idx > 0) then
            ! TODO: Calculate and update fengsha dust_effective_threshold
            ! call diag_mgr%update_diagnostic(this%diag_fengsha_dust_effective_threshold_idx, values)
         end if
      end if
      if (trim(this%active_scheme) == 'ginoux') then
         if (this%diag_ginoux_dust_horizontal_flux_idx > 0) then
            ! TODO: Calculate and update ginoux dust_horizontal_flux
            ! call diag_mgr%update_diagnostic(this%diag_ginoux_dust_horizontal_flux_idx, values)
         end if
         if (this%diag_ginoux_dust_moisture_correction_idx > 0) then
            ! TODO: Calculate and update ginoux dust_moisture_correction
            ! call diag_mgr%update_diagnostic(this%diag_ginoux_dust_moisture_correction_idx, values)
         end if
         if (this%diag_ginoux_dust_effective_threshold_idx > 0) then
            ! TODO: Calculate and update ginoux dust_effective_threshold
            ! call diag_mgr%update_diagnostic(this%diag_ginoux_dust_effective_threshold_idx, values)
         end if
      end if

   end subroutine update_dust_diagnostics

   !> Get required meteorological fields
   function get_dust_required_met_fields(this) result(field_names)
      class(ProcessDustInterface), intent(in) :: this
      character(len=32), allocatable :: field_names(:)

      character(len=32) :: common_fields(3)
      character(len=32), allocatable :: scheme_fields(:)
      integer :: n_common, n_scheme, n_total

      ! Common required fields
      common_fields(1) = 'ustar'
      common_fields(2) = 'solar_zenith_angle'
      common_fields(3) = 'leaf_area_index'
      n_common = 3

      ! Scheme-specific fields
      if (trim(this%active_scheme) == 'fengsha') then
         allocate(scheme_fields(12))
         scheme_fields(1) = 'IsLand'
         scheme_fields(2) = 'USTAR'
         scheme_fields(3) = 'LWI'
         scheme_fields(4) = 'GVF'
         scheme_fields(5) = 'LAI'
         scheme_fields(6) = 'FROCEAN'
         scheme_fields(7) = 'CLAYFRAC'
         scheme_fields(8) = 'SANDFRAC'
         scheme_fields(9) = 'FRSNO'
         scheme_fields(10) = 'RDRAG'
         scheme_fields(11) = 'SSM'
         scheme_fields(12) = 'USTAR_THRESHOLD'
         n_scheme = 12
      else if (trim(this%active_scheme) == 'ginoux') then
         allocate(scheme_fields(5))
         scheme_fields(1) = 'FRLAKE'
         scheme_fields(2) = 'GWETTOP'
         scheme_fields(3) = 'U10M'
         scheme_fields(4) = 'V10M'
         scheme_fields(5) = 'SSM'
         n_scheme = 5
      else
         allocate(scheme_fields(0))
         n_scheme = 0
      end if

      ! Combine fields (TODO: remove duplicates)
      n_total = n_common + n_scheme
      allocate(field_names(n_total))
      field_names(1:n_common) = common_fields
      if (n_scheme > 0) then
         field_names(n_common+1:n_total) = scheme_fields
      end if

   end function get_dust_required_met_fields

   !> Apply emission tendencies to chemical state
   subroutine apply_dust_emission_tendencies(this, state_manager, dt, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(StateManagerType), intent(inout) :: state_manager
      real(fp), intent(in) :: dt
      type(ErrorHandler), intent(inout) :: error_handler

      real(fp), pointer :: chem_state(:,:,:)  ! (n_species, n_levels, n_columns)
      integer :: i, i_spec, n_levels, n_columns

      ! Get chemical state
      chem_state => state_manager%get_chem_state()
      n_levels = state_manager%get_n_levels()
      n_columns = state_manager%get_n_columns()

      ! Apply emission tendencies to surface level
      do i = 1, this%n_species
         i_spec = this%species_indices(i)
         chem_state(i_spec, 1, :) = chem_state(i_spec, 1, :) + &
            this%emission_tendencies(i, :) * dt
      end do

   end subroutine apply_dust_emission_tendencies

   !> Finalize the dust process
   subroutine finalize_dust_process(this, error_handler)
      class(ProcessDustInterface), intent(inout) :: this
      type(ErrorHandler), intent(inout) :: error_handler

      ! Deallocate arrays
      if (allocated(this%species_names)) deallocate(this%species_names)
      if (allocated(this%species_indices)) deallocate(this%species_indices)
      if (allocated(this%emission_tendencies)) deallocate(this%emission_tendencies)

      ! Deallocate configuration and state
      if (associated(this%config)) then
         call this%config%finalize()
         deallocate(this%config)
      end if
      if (associated(this%state)) then
         call this%state%finalize()
         deallocate(this%state)
      end if

      ! Finalize scheme configurations and states
      if (allocated(this%config_fengsha)) then
         call this%config_fengsha%finalize()
         deallocate(this%config_fengsha)
      end if
      if (allocated(this%state_fengsha)) then
         call this%state_fengsha%finalize()
         deallocate(this%state_fengsha)
      end if
      if (allocated(this%config_ginoux)) then
         call this%config_ginoux%finalize()
         deallocate(this%config_ginoux)
      end if
      if (allocated(this%state_ginoux)) then
         call this%state_ginoux%finalize()
         deallocate(this%state_ginoux)
      end if

      this%is_initialized = .false.

   end subroutine finalize_dust_process

end module ProcessDustInterface_Mod