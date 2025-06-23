!> \file ExternalEmissionProcess_Mod.F90
!! \brief ExternalEmission process implementation
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! External emission process for applying emission tendencies from driver-loaded data
!!
module ExternalEmissionProcess_Mod
   use precision_mod
   use state_mod, only : StateContainerType
   use error_mod
   use ProcessInterface_Mod
   use DiagnosticInterface_Mod, only: DiagnosticFieldType, DiagnosticRegistryType, &
                                      DIAG_REAL_2D, DIAG_REAL_3D, DIAG_REAL_SCALAR
   use ColumnInterface_Mod, only: VirtualColumnType, ColumnProcessorType
   use EmissionConfigValidator_Mod, only: EmissionConfigValidatorType, ValidationResultType
   use EmisState_Mod, only: EmisStateType
   use ChemState_Mod, only: ChemStateType
   use defaultScheme_Mod

   implicit none
   private

   public :: ExternalEmissionProcessType

   !> ExternalEmission process type extending ColumnProcessInterface for optimal performance
   type, extends(ColumnProcessInterface) :: ExternalEmissionProcessType
      private

      ! Process-specific configuration
      character(len=32) :: selected_scheme = 'default'

      ! External emission configuration (IO handled by driver, process applies tendencies)
      type(EmissionConfigValidatorType) :: config_validator     ! Validates species mapping
      character(len=256) :: config_file_path = ''               ! Path to emission config
      
      ! Dynamic species mapping (configured via YAML)
      integer :: n_emitted_species = 0                          ! Number of species this process emits
      character(len=32), allocatable :: emitted_species(:)      ! Names of emitted species
      
      ! Process-specific emission factors/rates (user-defined based on process needs)
      ! NOTE: External emission data (NetCDF files, etc.) is read by the driver
      ! This process receives processed emission rates and applies them as tendencies
      integer, allocatable :: species_indices(:)                ! Indices in chemical mechanism
      real(fp), allocatable :: species_scale_factors(:)         ! Scale factors for each species
      
      ! Emission processing parameters
      logical :: time_interpolation = .true.                    ! Enable temporal interpolation
      logical :: spatial_interpolation = .true.                 ! Enable spatial interpolation
      real(fp) :: global_scale_factor = 1.0_fp                  ! Global scaling factor

      ! Column processing support - ALWAYS use column virtualization for performance
      logical :: column_processing_enabled = .true.
      logical :: supports_3d_processing = .false.  ! Prefer column processing
      
      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Process-specific data members
      ! ============================================================================
      ! TODO: Add your process-specific data members here
      ! Example data members for external emission processes:
      ! NOTE: External data files are read by the driver, not here
      ! integer, allocatable :: species_indices(:)        ! Indices in chemical state
      ! real(fp), allocatable :: scaling_factors(:)       ! Process-specific scaling
      ! real(fp), allocatable :: temporal_weights(:)      ! Time-varying factors
      ! character(len=32), allocatable :: category_names(:) ! Emission categories
      ! ============================================================================

   contains
      ! Required ProcessInterface methods
      procedure :: init => externalemission_process_init
      procedure :: run => externalemission_process_run
      procedure :: finalize => externalemission_process_finalize
      
      ! Enhanced diagnostic methods
      procedure :: register_diagnostics => externalemission_register_diagnostics
      procedure :: update_diagnostics => externalemission_update_diagnostics
      
      ! Column processing methods
      procedure :: run_column => externalemission_run_column
      procedure :: supports_column_processing => externalemission_supports_column_processing
      
      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Additional methods
      ! ============================================================================
      ! TODO: Add your process-specific methods here
      ! Example methods for external emission processes:
      ! procedure :: validate_emission_config => externalemission_validate_config
      ! procedure :: setup_species_mapping => externalemission_setup_mapping
      ! procedure :: apply_emission_tendencies => externalemission_apply_tendencies
      ! procedure :: get_emission_rates => externalemission_get_rates
      ! ============================================================================
   end type ExternalEmissionProcessType

contains

   !> Initialize ExternalEmission process
   subroutine externalemission_process_init(this, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      type(ConfigDataType), pointer :: config
      character(len=256) :: message

      rc = CC_SUCCESS

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('externalemission_process_init', &
                                  'Initializing ExternalEmission process')

      ! Set process metadata
      this%name = 'ExternalEmission'
      this%version = '1.0'
      this%description = 'External emission process for applying emission tendencies from driver-loaded data'

      ! Get configuration
      config => container%get_config_ptr()
      if (.not. associated(config)) then
         call error_mgr%report_error(ERROR_NOT_INITIALIZED, &
                                    'Configuration not available', rc, &
                                    'externalemission_process_init', &
                                    'Ensure StateContainer is properly initialized')
         call error_mgr%pop_context()
         return
      endif

      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Process-specific initialization
      ! ============================================================================
      ! TODO: Add your process-specific initialization here
      !
      ! For EXTERNAL EMISSION processes:
      ! - External emission data (NetCDF files, etc.) is READ BY THE DRIVER
      ! - This process receives processed emission rates through the StateContainer
      ! - Your process applies these rates as tendencies to species concentrations
      ! - Use EmisState to accumulate emissions from all processes
      !
      ! Example initialization for external emissions:
      ! call this%validate_emission_config(config, rc)
      ! call this%setup_species_mapping(config, rc)
      ! call this%initialize_emission_factors(config, rc)
      ! ============================================================================


      ! Validate inputs
      if (.not. this%validate_emission_inputs(container, rc)) then
         call error_mgr%report_error(ERROR_VALIDATION, &
                                    'ExternalEmission emission process input validation failed', rc)
         call error_mgr%pop_context()
         return
      end if

      ! Register diagnostics with the new diagnostic system
      call this%register_diagnostics(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_DIAGNOSTIC_REGISTRATION, &
                                    'Failed to register ExternalEmission diagnostics', rc)
         call error_mgr%pop_context()
         return
      endif

      this%is_initialized = .true.
      this%is_active = .true.

      ! Register diagnostics after successful initialization
      call this%register_diagnostics(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_DIAGNOSTIC_REGISTRATION, &
                                    'Failed to register diagnostics', rc)
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()

      write(message, '(A,A,A)') 'ExternalEmission process initialized with scheme: ', &
                                trim(this%selected_scheme)
      call error_mgr%report_info(message)
      
      call error_mgr%pop_context()

   end subroutine externalemission_process_init

   !> Run ExternalEmission process
   subroutine externalemission_process_run(this, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS

      if (.not. this%is_ready()) then
         rc = CC_FAILURE
         return
      end if

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('externalemission_process_run', &
                                  'ExternalEmission process execution')

      ! Use column virtualization if supported and enabled
      if (this%supports_column_processing() .and. this%column_processing_enabled) then
         grid_mgr => container%get_grid_manager_ptr()
         if (associated(grid_mgr) .and. grid_mgr%virtualization_enabled) then
            call this%run_column_batch(container, rc)
         else
            call this%run_3d_processing(container, rc)
         endif
      else
         call this%run_3d_processing(container, rc)
      endif

      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_PROCESS_EXECUTION, &
                                    'ExternalEmission process execution failed', rc)
         call error_mgr%pop_context()
         return
      endif

      ! Update diagnostics after successful execution
      call this%update_diagnostics(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_DIAGNOSTIC_UPDATE, &
                                    'Failed to update diagnostics', rc)
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()

   end subroutine externalemission_process_run

   !> Finalize ExternalEmission process
   subroutine externalemission_process_finalize(this, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      this%is_initialized = .false.
      this%is_active = .false.

   end subroutine externalemission_process_finalize

   !> Run column batch processing (optimized for performance)
   subroutine run_column_batch(this, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr
      type(ColumnIteratorType) :: iterator
      type(VirtualColumnType) :: column
      type(ErrorManagerType), pointer :: error_mgr
      integer :: local_rc

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      
      grid_mgr => container%get_grid_manager_ptr()
      call grid_mgr%get_column_iterator(iterator, rc)
      if (rc /= CC_SUCCESS) return

      ! Process all columns using virtualization
      do while (iterator%has_next())
         call iterator%next(column)
         
         call this%run_column(column, container, local_rc)
         if (local_rc /= CC_SUCCESS) then
            rc = local_rc
            call error_mgr%report_error(ERROR_PROCESS_EXECUTION, &
                                       'Column processing failed', rc)
            return
         endif
      end do

   end subroutine run_column_batch

   !> Run 3D processing (fallback when column virtualization not available)
   subroutine run_3d_processing(this, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(EmisStateType), pointer :: emis_state
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()

      ! Get state pointers (no allocation, just access)
      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()
      emis_state => container%get_emis_state_ptr()

      ! Execute scheme-specific calculations on 3D data
      select case (trim(this%selected_scheme))
      case ('default')
         call this%run_default_scheme_3d(met_state, chem_state, emis_state, rc)
      case default
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
                                    'Unknown scheme: ' // trim(this%selected_scheme), rc)
         return
      end select

   end subroutine run_3d_processing

   !> Column processing implementation (required by ColumnProcessInterface)
   subroutine externalemission_run_column(this, column, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Column processing implementation
      ! ============================================================================
      ! TODO: Implement your column-specific processing here
      ! Use column%met_data, column%chem_data, column%emis_data pointers
      ! to access data without additional allocations
      
      ! Example for emission process:
      ! IMPORTANT: External emission data is provided by the driver
      ! Your process receives emission rates through the StateContainer
      ! and applies them as tendencies to species concentrations
      !
      ! Extract meteorological data from column (no allocation - use pointers)
      ! real(fp), pointer :: temperature(:) => column%met_data(:,1)
      ! real(fp), pointer :: pressure(:) => column%met_data(:,2)
      ! real(fp), pointer :: humidity(:) => column%met_data(:,3)
      ! 
      ! Get emission rates for this column (provided by driver via StateContainer)
      ! real(fp), pointer :: emission_rates(:,:) => container%get_emission_rates()
      ! 
      ! Apply emission tendencies to species concentrations
      ! call this%apply_emission_tendencies(column, emission_rates, rc)
      ! if (rc /= CC_SUCCESS) return
      !
      ! Update EmisState for diagnostic accumulation
      ! call container%emis_state%accumulate_emissions(species_idx, column_emissions)
      ! ============================================================================

      rc = CC_SUCCESS

      ! For now, use scheme-specific processing
      select case (trim(this%selected_scheme))
      case ('default')
         call this%run_default_column_scheme(column, container, rc)
      case default
         rc = CC_FAILURE
         return
      end select

   end subroutine externalemission_run_column

   !> Check if this process supports column processing
   function externalemission_supports_column_processing(this) result(supports)
      class(ExternalEmissionProcessType), intent(in) :: this
      logical :: supports

      supports = this%column_processing_enabled

   end function externalemission_supports_column_processing

   !> Register diagnostics for this process
   subroutine externalemission_register_diagnostics(this, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(DiagnosticRegistryType), pointer :: registry
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('externalemission_register_diagnostics', &
                                  'Registering ExternalEmission diagnostics')

      ! Call parent to register process
      call this%ProcessInterface%register_diagnostics(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%pop_context()
         return
      endif

      ! Get diagnostic registry for this process
      diag_mgr => container%get_diagnostic_manager_ptr()
      call diag_mgr%get_process_registry(this%name, registry, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%pop_context()
         return
      endif

      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Register process-specific diagnostics
      ! ============================================================================
      ! TODO: Register your process-specific diagnostic fields here
      
      ! Example diagnostic registrations:
      ! call registry%register_field('surface_emission_flux', &
      !                             'Surface emission flux', &
      !                             'mol m-2 s-1', DIAG_REAL_2D, this%name, rc)
      ! if (rc /= CC_SUCCESS) then
      !    call error_mgr%pop_context()
      !    return
      ! endif
      !
      ! call registry%register_field('total_column_emission', &
      !                             'Total column emission', &
      !                             'mol m-2 s-1', DIAG_REAL_2D, this%name, rc)
      ! ============================================================================

      call error_mgr%pop_context()

   end subroutine externalemission_register_diagnostics

   !> Update diagnostics for this process
   subroutine externalemission_update_diagnostics(this, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(DiagnosticRegistryType), pointer :: registry
      type(DiagnosticFieldType), pointer :: diag_field
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()

      ! Get diagnostic registry for this process
      diag_mgr => container%get_diagnostic_manager_ptr()
      call diag_mgr%get_process_registry(this%name, registry, rc)
      if (rc /= CC_SUCCESS) return

      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Update process-specific diagnostics
      ! ============================================================================
      ! TODO: Update your diagnostic fields with current process data
      
      ! Example diagnostic updates:
      ! diag_field => registry%get_field_ptr('surface_emission_flux', rc)
      ! if (rc == CC_SUCCESS .and. associated(diag_field)) then
      !    call diag_field%update_data(array_2d=this%surface_emission_flux)
      ! endif
      !
      ! diag_field => registry%get_field_ptr('total_column_emission', rc)
      ! if (rc == CC_SUCCESS .and. associated(diag_field)) then
      !    call diag_field%update_data(array_2d=this%total_emission)
      ! endif
      ! ============================================================================

   end subroutine externalemission_update_diagnostics

   !> Run default scheme on 3D data (fallback method)
   subroutine run_default_scheme_3d(this, met_state, chem_state, emis_state, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(MetStateType), intent(in) :: met_state
      type(ChemStateType), intent(inout) :: chem_state
      type(EmisStateType), intent(inout) :: emis_state
      integer, intent(out) :: rc

      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: 3D processing implementation
      ! ============================================================================
      ! TODO: Implement your 3D processing here
      ! Use pointers to access data arrays without additional allocation
      
      ! Access meteorological data with pointers (no allocation)
      ! real(fp), pointer :: temperature(:,:,:) => met_state%temperature
      ! real(fp), pointer :: pressure(:,:,:) => met_state%pressure
      ! real(fp), pointer :: humidity(:,:,:) => met_state%humidity
      ! real(fp), pointer :: wind_u(:,:,:) => met_state%wind_u
      ! real(fp), pointer :: wind_v(:,:,:) => met_state%wind_v
      ! real(fp), pointer :: surface_temperature(:,:) => met_state%surface_temperature
      ! 
      ! Access emission arrays with pointers
      ! real(fp), pointer :: surface_emissions(:,:,:) => emis_state%surface_emissions
      ! 
      ! Call your emission calculation routine
      ! call default_calculate_emission_3d( &
      !    temperature, pressure, humidity, &
      !    wind_u, wind_v, surface_temperature, &
      !    surface_emissions, this%dt, rc)
      ! ============================================================================

      rc = CC_SUCCESS

      ! Placeholder implementation
      write(*,*) 'Running default scheme on 3D data - implement me!'

   end subroutine run_default_scheme_3d

   !> Run default scheme on column data (optimized method)
   subroutine run_default_column_scheme(this, column, container, rc)
      class(ExternalEmissionProcessType), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! ============================================================================
      ! USER CUSTOMIZATION SECTION: Column processing implementation
      ! ============================================================================
      ! TODO: Implement your column processing here
      ! Use column data pointers for efficient access (no allocation needed)
      
      ! Access column meteorological data (pointers to column%met_data)
      ! real(fp), pointer :: temperature(:) => column%met_data(:,1)
      ! real(fp), pointer :: pressure(:) => column%met_data(:,2)
      ! real(fp), pointer :: humidity(:) => column%met_data(:,3)
      ! 
      ! Access surface properties
      ! real(fp) :: surface_temp = column%surface_temperature
      ! real(fp) :: roughness = column%roughness_length
      ! 
      ! Calculate emission for this column
      ! call default_calculate_column_emission( &
      !    column%nz, temperature, pressure, humidity, &
      !    surface_temp, roughness, &
      !    column%emis_data, this%dt, rc)
      ! ============================================================================

      rc = CC_SUCCESS

      ! Placeholder implementation
      write(*,*) 'Running default scheme on column data at (', &
                column%global_i, ',', column%global_j, ') - implement me!'

   end subroutine run_default_column_scheme

end module ExternalEmissionProcess_Mod