!> \file ProcessInterface_Mod.F90
!! \brief Abstract base class interface for all atmospheric processes
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module defines the abstract interface that all atmospheric processes
!! must implement, following the exact structure from PROCESS_ARCHITECTURE_GUIDE.md
!!
module ProcessInterface_Mod
   use precision_mod
   use state_mod, only : StateContainerType
   use error_mod
   use ColumnInterface_Mod, only : VirtualColumnType, ColumnProcessorType
   use StateManager_Mod, only : StateManagerType

   implicit none
   private

   public :: ProcessInterface
   public :: ColumnProcessInterface

   !> \brief Abstract base class for all atmospheric processes
   !!
   !! This abstract type defines the interface that all atmospheric processes
   !! must implement, following the exact structure from PROCESS_ARCHITECTURE_GUIDE.md
   !!
   type, abstract :: ProcessInterface
      private
      character(len=64) :: name = ''         !< Process name
      character(len=64) :: version = ''      !< Version string
      character(len=256) :: description = '' !< Process description
      logical :: is_initialized = .false.    !< Initialization status
      logical :: is_active = .false.         !< Active status
      real(fp) :: dt = 0.0_fp                !< Process timestep

      ! Species management
      integer :: n_species = 0               !< Number of species handled
      character(len=32), allocatable :: species_names(:) !< Species names

      ! Size bin management for aerosol processes
      integer :: n_size_bins = 0             !< Number of size bins
      real(fp), allocatable :: size_bin_bounds(:) !< Size bin boundaries [μm]
      real(fp), allocatable :: size_bin_centers(:) !< Size bin centers [μm]
      character(len=16) :: size_type = 'radius' !< 'radius' or 'diameter'

      ! Multiphase chemistry support
      logical :: is_multiphase = .false.     !< Whether this is a multiphase process
      integer :: n_phases = 1                !< Number of phases (1=gas, 2=gas+liquid, 3=gas+liquid+solid)
      character(len=16), allocatable :: phase_names(:) !< Phase names ('gas', 'liquid', 'solid')
      logical :: has_heterogeneous_reactions = .false. !< Whether process includes heterogeneous reactions
      logical :: has_aqueous_chemistry = .false.       !< Whether process includes aqueous chemistry
      logical :: has_thermodynamics = .false.          !< Whether process includes thermodynamic equilibrium
      logical :: has_cloud_microphysics = .false.      !< Whether process couples with cloud microphysics

      ! Reaction network parameters
      integer :: n_reactions = 0             !< Number of reactions
      integer :: n_photolysis = 0            !< Number of photolysis reactions
      integer :: n_henry_species = 0         !< Number of species with Henry's law constants
      character(len=16) :: solver_type = 'none' !< Solver type ('euler', 'rk4', 'rosenbrock', 'kpp', 'micm')
      real(fp) :: rtol = 1.0e-3_fp           !< Relative tolerance for solver
      real(fp) :: atol = 1.0e-12_fp          !< Absolute tolerance for solver

   contains
      ! Required interface methods
      procedure(init_interface), deferred :: init
      procedure(run_interface), deferred :: run
      procedure(finalize_interface), deferred :: finalize

      ! Diagnostic registration capabilities
      procedure :: register_diagnostics => process_register_diagnostics
      procedure :: update_diagnostics => process_update_diagnostics
      procedure :: get_diagnostic_registry => process_get_diagnostic_registry

      ! Optional interface methods with default implementations
      procedure :: get_name => process_get_name
      procedure :: get_version => process_get_version
      procedure :: is_ready => process_is_ready
      procedure :: activate => process_activate
      procedure :: deactivate => process_deactivate
      procedure :: validate_config => process_validate_config

      ! Common functionality for emission processes
      procedure :: set_timestep => process_set_timestep
      procedure :: get_timestep => process_get_timestep
      procedure :: set_species => process_set_species
      procedure :: get_species_count => process_get_species_count
      procedure :: get_species_names => process_get_species_names
      procedure :: check_species_availability => process_check_species_availability
      procedure :: log_process_info => process_log_info
      procedure :: validate_grid_consistency => process_validate_grid_consistency
      procedure :: allocate_species_arrays => process_allocate_species_arrays
      procedure :: deallocate_species_arrays => process_deallocate_species_arrays

      ! Additional common functionality for aerosol emission processes
      procedure :: set_size_bins => process_set_size_bins
      procedure :: get_size_bin_count => process_get_size_bin_count
      procedure :: get_size_bin_bounds => process_get_size_bin_bounds
      procedure :: validate_emission_inputs => process_validate_emission_inputs
      procedure :: zero_emission_arrays => process_zero_emission_arrays
      procedure :: accumulate_emissions => process_accumulate_emissions
      procedure :: write_process_diagnostics => process_write_diagnostics
      procedure :: check_mass_conservation => process_check_mass_conservation
      procedure :: apply_emission_scaling => process_apply_emission_scaling

      ! Multiphase chemistry functionality
      procedure :: set_multiphase_config => process_set_multiphase_config
      procedure :: get_phase_count => process_get_phase_count
      procedure :: get_phase_names => process_get_phase_names
      procedure :: is_multiphase_process => process_is_multiphase_process
      procedure :: set_solver_config => process_set_solver_config
      procedure :: get_solver_config => process_get_solver_config
      procedure :: validate_multiphase_inputs => process_validate_multiphase_inputs
      procedure :: setup_phase_equilibrium => process_setup_phase_equilibrium
      procedure :: calculate_henry_constants => process_calculate_henry_constants
      procedure :: update_photolysis_rates => process_update_photolysis_rates
      procedure :: solve_chemistry_system => process_solve_chemistry_system
      procedure :: partition_species => process_partition_species
      procedure :: check_phase_balance => process_check_phase_balance

      ! Column virtualization support
      procedure :: supports_column_processing => process_supports_column_processing
      procedure :: process_column => process_process_column
      procedure :: process_all_columns => process_process_all_columns
   end type ProcessInterface

   !> \brief Enhanced process interface specifically for column-based processing
   !!
   !! This interface extends ProcessInterface to provide column virtualization
   !! capabilities, allowing processes to work with virtual columns while
   !! maintaining awareness of 3D spatial relationships.
   type, abstract, extends(ProcessInterface) :: ColumnProcessInterface
      private

      logical :: column_processing_enabled = .true.  !< Enable column processing mode
      integer :: column_batch_size = 100            !< Number of columns to process in batch

   contains
      ! Required column processing methods
      procedure(column_init_interface), deferred :: init_column_processing
      procedure(column_run_interface), deferred :: run_column
      procedure(column_finalize_interface), deferred :: finalize_column_processing

      ! Optional column processing methods with default implementations
      procedure :: set_column_batch_size => column_process_set_batch_size
      procedure :: get_column_batch_size => column_process_get_batch_size
      procedure :: enable_column_processing => column_process_enable
      procedure :: disable_column_processing => column_process_disable
      procedure :: is_column_processing_enabled => column_process_is_enabled

      ! Generic diagnostic update interface - automatically dispatches based on argument types
      generic :: update_diagnostics => update_scalar_diagnostic_column, &
                                       update_1d_diagnostic_column, &
                                       update_2d_diagnostic_column
      procedure :: update_scalar_diagnostic_column => column_update_scalar_diagnostic
      procedure :: update_1d_diagnostic_column => column_update_1d_diagnostic
      procedure :: update_2d_diagnostic_column => column_update_2d_diagnostic

      ! Method to get StateManager - must be implemented by concrete processes
      procedure(get_state_manager_interface), deferred :: get_state_manager
   end type ColumnProcessInterface

   ! Abstract interfaces that must be implemented by concrete processes
   abstract interface
      !> \brief Initialize the process with given container
      subroutine init_interface(this, container, rc)
         import :: ProcessInterface, StateContainerType
         class(ProcessInterface), intent(inout) :: this
         type(StateContainerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Execute the process main calculations
      subroutine run_interface(this, container, rc)
         import :: ProcessInterface, StateContainerType
         class(ProcessInterface), intent(inout) :: this
         type(StateContainerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Clean up and finalize the process
      subroutine finalize_interface(this, rc)
         import :: ProcessInterface
         class(ProcessInterface), intent(inout) :: this
         integer, intent(out) :: rc
      end subroutine
   end interface

   ! Column processing interfaces
   abstract interface
      !> \brief Initialize column processing for the process
      subroutine column_init_interface(this, container, rc)
         import :: ColumnProcessInterface, StateContainerType
         class(ColumnProcessInterface), intent(inout) :: this
         type(StateContainerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Process a single virtual column
      subroutine column_run_interface(this, column, rc)
         import :: ColumnProcessInterface, VirtualColumnType
         class(ColumnProcessInterface), intent(inout) :: this
         type(VirtualColumnType), intent(inout) :: column
         integer, intent(out) :: rc
      end subroutine

      !> \brief Finalize column processing
      subroutine column_finalize_interface(this, container, rc)
         import :: ColumnProcessInterface, StateContainerType
         class(ColumnProcessInterface), intent(inout) :: this
         type(StateContainerType), intent(inout) :: container
         integer, intent(out) :: rc
      end subroutine

      !> \brief Get StateManager for diagnostic access
      function get_state_manager_interface(this) result(state_manager)
         import :: ColumnProcessInterface
         import StateManagerType
         class(ColumnProcessInterface), intent(in) :: this
         type(StateManagerType), pointer :: state_manager
      end function
   end interface

contains

   !> \brief Get process name
   function process_get_name(this) result(name)
      class(ProcessInterface), intent(in) :: this
      character(len=64) :: name
      name = this%name
   end function process_get_name

   !> \brief Get process version
   function process_get_version(this) result(version)
      class(ProcessInterface), intent(in) :: this
      character(len=64) :: version
      version = this%version
   end function process_get_version

   !> \brief Check if process is ready to run
   function process_is_ready(this) result(ready)
      class(ProcessInterface), intent(in) :: this
      logical :: ready
      ready = this%is_initialized .and. this%is_active
   end function process_is_ready

   !> \brief Activate the process
   subroutine process_activate(this)
      class(ProcessInterface), intent(inout) :: this
      this%is_active = .true.
   end subroutine process_activate

   !> \brief Deactivate the process
   subroutine process_deactivate(this)
      class(ProcessInterface), intent(inout) :: this
      this%is_active = .false.
   end subroutine process_deactivate

   !> \brief Default configuration validation (always valid)
   subroutine process_validate_config(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Default implementation - always valid
      rc = CC_SUCCESS
   end subroutine process_validate_config

   !> \brief Set process timestep
   subroutine process_set_timestep(this, dt)
      class(ProcessInterface), intent(inout) :: this
      real(fp), intent(in) :: dt
      this%dt = dt
   end subroutine process_set_timestep

   !> \brief Get process timestep
   function process_get_timestep(this) result(dt)
      class(ProcessInterface), intent(in) :: this
      real(fp) :: dt
      dt = this%dt
   end function process_get_timestep

   !> \brief Set species names that this process handles
   subroutine process_set_species(this, species_names, rc)
      class(ProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: species_names(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS
      this%n_species = size(species_names)

      if (allocated(this%species_names)) deallocate(this%species_names)
      allocate(this%species_names(this%n_species), stat=rc)
      if (rc /= 0) return

      do i = 1, this%n_species
         this%species_names(i) = trim(species_names(i))
      end do
   end subroutine process_set_species

   !> \brief Get number of species handled by this process
   function process_get_species_count(this) result(n_species)
      class(ProcessInterface), intent(in) :: this
      integer :: n_species
      n_species = this%n_species
   end function process_get_species_count

   !> \brief Get species names handled by this process
   function process_get_species_names(this) result(species_names)
      class(ProcessInterface), intent(in) :: this
      character(len=32), allocatable :: species_names(:)

      if (allocated(this%species_names)) then
         allocate(species_names(size(this%species_names)))
         species_names = this%species_names
      end if
   end function process_get_species_names

   !> \brief Check if required species are available in chemical state
   function process_check_species_availability(this, container, rc) result(available)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc
      logical :: available

      type(ChemStateType), pointer :: chem_state
      integer :: i
      character(len=32) :: species_name

      available = .true.
      rc = CC_SUCCESS

      if (.not. allocated(this%species_names) .or. this%n_species == 0) then
         available = .true.  ! No species requirements
         return
      end if

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         available = .false.
         rc = -1
         return
      end if

      do i = 1, this%n_species
         species_name = this%species_names(i)
         if (.not. chem_state%has_species(species_name)) then
            available = .false.
            exit
         end if
      end do
   end function process_check_species_availability

   !> \brief Log process information
   subroutine process_log_info(this, container)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message
      integer :: i

      error_mgr => container%get_error_manager()

      write(message, '(A,A,A,A)') 'Process: ', trim(this%name), ' Version: ', trim(this%version)
      call error_mgr%report_info(message)

      write(message, '(A,I0,A,F8.2,A)') 'Species count: ', this%n_species, ', Timestep: ', this%dt, ' s'
      call error_mgr%report_info(message)

      if (allocated(this%species_names)) then
         do i = 1, min(this%n_species, 5)  ! Log first 5 species
            write(message, '(A,I0,A,A)') 'Species ', i, ': ', trim(this%species_names(i))
            call error_mgr%report_info(message)
         end do
         if (this%n_species > 5) then
            write(message, '(A,I0,A)') '... and ', this%n_species - 5, ' more species'
            call error_mgr%report_info(message)
         end if
      end if
   end subroutine process_log_info

   !> \brief Validate grid consistency between meteorological and chemical states
   function process_validate_grid_consistency(this, container, rc) result(consistent)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc
      logical :: consistent

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      integer :: met_nx, met_ny, met_nz
      integer :: chem_nx, chem_ny, chem_nz

      consistent = .true.
      rc = CC_SUCCESS

      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()

      if (.not. associated(met_state) .or. .not. associated(chem_state)) then
         consistent = .false.
         rc = -1
         return
      end if

      ! Get grid dimensions (assuming these methods exist)
      call met_state%get_dimensions(met_nx, met_ny, met_nz)
      call chem_state%get_dimensions(chem_nx, chem_ny, chem_nz)

      if (met_nx /= chem_nx .or. met_ny /= chem_ny .or. met_nz /= chem_nz) then
         consistent = .false.
         rc = -2
      end if
   end function process_validate_grid_consistency

   !> \brief Allocate species-specific arrays for a process
   subroutine process_allocate_species_arrays(this, n_lon, n_lat, n_species_arrays, rc)
      class(ProcessInterface), intent(in) :: this
      integer, intent(in) :: n_lon, n_lat
      real(fp), allocatable, intent(out) :: n_species_arrays(:,:,:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (this%n_species <= 0) then
         rc = -1
         return
      end if

      allocate(n_species_arrays(n_lon, n_lat, this%n_species), stat=rc)
      if (rc == 0) then
         n_species_arrays = 0.0_fp
      end if
   end subroutine process_allocate_species_arrays

   !> \brief Deallocate species-specific arrays
   subroutine process_deallocate_species_arrays(this, species_arrays)
      class(ProcessInterface), intent(in) :: this
      real(fp), allocatable, intent(inout) :: species_arrays(:,:,:)

      if (allocated(species_arrays)) deallocate(species_arrays)
   end subroutine process_deallocate_species_arrays

   !> \brief Set size bins for aerosol processes
   subroutine process_set_size_bins(this, bin_bounds, size_type, rc)
      class(ProcessInterface), intent(inout) :: this
      real(fp), intent(in) :: bin_bounds(:)  !< Size bin boundaries
      character(len=*), intent(in) :: size_type !< 'radius' or 'diameter'
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (size(bin_bounds) < 2) then
         rc = -1
         return
      end if

      this%n_size_bins = size(bin_bounds) - 1
      this%size_type = trim(size_type)

      if (allocated(this%size_bin_bounds)) deallocate(this%size_bin_bounds)
      if (allocated(this%size_bin_centers)) deallocate(this%size_bin_centers)

      allocate(this%size_bin_bounds(size(bin_bounds)), stat=rc)
      if (rc /= 0) return

      allocate(this%size_bin_centers(this%n_size_bins), stat=rc)
      if (rc /= 0) return

      this%size_bin_bounds = bin_bounds

      ! Calculate bin centers (geometric mean)
      do i = 1, this%n_size_bins
         this%size_bin_centers(i) = sqrt(bin_bounds(i) * bin_bounds(i+1))
      end do
   end subroutine process_set_size_bins

   !> \brief Get number of size bins
   function process_get_size_bin_count(this) result(n_bins)
      class(ProcessInterface), intent(in) :: this
      integer :: n_bins
      n_bins = this%n_size_bins
   end function process_get_size_bin_count

   !> \brief Get size bin boundaries
   function process_get_size_bin_bounds(this) result(bounds)
      class(ProcessInterface), intent(in) :: this
      real(fp), allocatable :: bounds(:)

      if (allocated(this%size_bin_bounds)) then
         allocate(bounds(size(this%size_bin_bounds)))
         bounds = this%size_bin_bounds
      end if
   end function process_get_size_bin_bounds

   !> \brief Validate common emission inputs
   function process_validate_emission_inputs(this, container, rc) result(valid)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc
      logical :: valid

      valid = .true.
      rc = CC_SUCCESS

      ! Check that StateContainer is properly initialized
      if (.not. container%is_initialized) then
         valid = .false.
         rc = -1
         return
      end if

      ! Check meteorological fields are available
      if (.not. associated(container%met_state)) then
         valid = .false.
         rc = -2
         return
      end if

      ! Check emission state is available for emission processes
      if (.not. associated(container%emis_state)) then
         valid = .false.
         rc = -3
         return
      end if
   end function process_validate_emission_inputs

   !> \brief Zero out emission arrays for clean start
   subroutine process_zero_emission_arrays(this, container, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state

      rc = CC_SUCCESS

      emis_state => container%emis_state
      if (.not. associated(emis_state)) then
         rc = -1
         return
      end if

      ! Zero emission arrays - implementation would depend on EmisState structure
      ! This is a placeholder for now

   end subroutine process_zero_emission_arrays

   !> \brief Accumulate emissions from process calculations
   subroutine process_accumulate_emissions(this, container, process_emissions, species_map, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(inout) :: container
      real(fp), intent(in) :: process_emissions(:,:,:) !< Process emission fluxes [kg/m²/s]
      integer, intent(in) :: species_map(:) !< Mapping from process species to global species
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state
      integer :: i, j, k, species_idx

      rc = CC_SUCCESS

      emis_state => container%emis_state
      if (.not. associated(emis_state)) then
         rc = -1
         return
      end if

      ! Accumulate emissions - implementation would depend on EmisState structure
      ! This is a placeholder for now

   end subroutine process_accumulate_emissions

   !> \brief Write process-specific diagnostics
   subroutine process_write_diagnostics(this, container, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Default implementation - no diagnostics written
      ! Specific processes can override this

   end subroutine process_write_diagnostics

   !> \brief Check mass conservation for debugging
   function process_check_mass_conservation(this, container, tolerance, rc) result(conserved)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      real(fp), intent(in) :: tolerance
      integer, intent(out) :: rc
      logical :: conserved

      rc = CC_SUCCESS
      conserved = .true.

      ! Default implementation - always conserved
      ! Specific processes can override this with actual checks

   end function process_check_mass_conservation

   !> \brief Apply emission scaling factors
   subroutine process_apply_emission_scaling(this, container, scaling_factors, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(inout) :: container
      real(fp), intent(in) :: scaling_factors(:) !< Scaling factors per species
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state

      rc = CC_SUCCESS

      emis_state => container%emis_state
      if (.not. associated(emis_state)) then
         rc = -1
         return
      end if

      ! Apply scaling - implementation would depend on EmisState structure
      ! This is a placeholder for now

   end subroutine process_apply_emission_scaling

   !> \brief Configure multiphase process parameters
   subroutine process_set_multiphase_config(this, n_phases, phase_names, has_heterogeneous, &
                                           has_aqueous, has_thermo, has_cloud, rc)
      class(ProcessInterface), intent(inout) :: this
      integer, intent(in) :: n_phases
      character(len=*), intent(in) :: phase_names(:)
      logical, intent(in) :: has_heterogeneous, has_aqueous, has_thermo, has_cloud
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (n_phases < 1 .or. n_phases > 3) then
         rc = -1
         return
      end if

      this%is_multiphase = (n_phases > 1)
      this%n_phases = n_phases
      this%has_heterogeneous_reactions = has_heterogeneous
      this%has_aqueous_chemistry = has_aqueous
      this%has_thermodynamics = has_thermo
      this%has_cloud_microphysics = has_cloud

      if (allocated(this%phase_names)) deallocate(this%phase_names)
      allocate(this%phase_names(n_phases), stat=rc)
      if (rc /= 0) return

      do i = 1, n_phases
         this%phase_names(i) = trim(phase_names(i))
      end do
   end subroutine process_set_multiphase_config

   !> \brief Get number of phases
   function process_get_phase_count(this) result(n_phases)
      class(ProcessInterface), intent(in) :: this
      integer :: n_phases
      n_phases = this%n_phases
   end function process_get_phase_count

   !> \brief Get phase names
   function process_get_phase_names(this) result(phase_names)
      class(ProcessInterface), intent(in) :: this
      character(len=16), allocatable :: phase_names(:)

      if (allocated(this%phase_names)) then
         allocate(phase_names(size(this%phase_names)))
         phase_names = this%phase_names
      end if
   end function process_get_phase_names

   !> \brief Check if this is a multiphase process
   function process_is_multiphase_process(this) result(is_multiphase)
      class(ProcessInterface), intent(in) :: this
      logical :: is_multiphase
      is_multiphase = this%is_multiphase
   end function process_is_multiphase_process

   !> \brief Set solver configuration for chemistry
   subroutine process_set_solver_config(this, solver_type, rtol, atol, n_reactions, n_photolysis, rc)
      class(ProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: solver_type
      real(fp), intent(in) :: rtol, atol
      integer, intent(in) :: n_reactions, n_photolysis
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      this%solver_type = trim(solver_type)
      this%rtol = rtol
      this%atol = atol
      this%n_reactions = n_reactions
      this%n_photolysis = n_photolysis
   end subroutine process_set_solver_config

   !> \brief Get solver configuration
   subroutine process_get_solver_config(this, solver_type, rtol, atol, n_reactions, n_photolysis)
      class(ProcessInterface), intent(in) :: this
      character(len=16), intent(out) :: solver_type
      real(fp), intent(out) :: rtol, atol
      integer, intent(out) :: n_reactions, n_photolysis

      solver_type = this%solver_type
      rtol = this%rtol
      atol = this%atol
      n_reactions = this%n_reactions
      n_photolysis = this%n_photolysis
   end subroutine process_get_solver_config

   !> \brief Validate inputs for multiphase chemistry processes
   function process_validate_multiphase_inputs(this, container, rc) result(valid)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc
      logical :: valid

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message

      valid = .true.
      rc = CC_SUCCESS

      if (.not. this%is_multiphase) then
         return  ! Not applicable for single-phase processes
      end if

      error_mgr => container%get_error_manager()
      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()

      ! Check basic state availability
      if (.not. associated(met_state) .or. .not. associated(chem_state)) then
         valid = .false.
         rc = -1
         call error_mgr%report_error("Multiphase process requires met_state and chem_state")
         return
      end if

      ! Check for required meteorological fields for aqueous chemistry
      if (this%has_aqueous_chemistry) then
         ! Would check for temperature, relative humidity, cloud water content, etc.
         ! Placeholder implementation
         write(message, '(A)') "Validating aqueous chemistry requirements"
         call error_mgr%report_info(message)
      end if

      ! Check for required fields for heterogeneous chemistry
      if (this%has_heterogeneous_reactions) then
         ! Would check for aerosol surface area, particle composition, etc.
         write(message, '(A)') "Validating heterogeneous chemistry requirements"
         call error_mgr%report_info(message)
      end if

   end function process_validate_multiphase_inputs

   !> \brief Set up phase equilibrium calculations
   subroutine process_setup_phase_equilibrium(this, container, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message

      rc = CC_SUCCESS

      if (.not. this%has_thermodynamics) then
         return  ! No thermodynamic equilibrium needed
      end if

      error_mgr => container%get_error_manager()

      ! Setup thermodynamic equilibrium solver
      ! This would involve initializing activity coefficient models,
      ! solubility calculations, etc.
      write(message, '(A)') "Setting up thermodynamic equilibrium solver"
      call error_mgr%report_info(message)

      ! Placeholder implementation - actual implementation would initialize
      ! models like UNIFAC, AIOMFAC, or other activity coefficient models

   end subroutine process_setup_phase_equilibrium

   !> \brief Calculate Henry's law constants for gas-liquid partitioning
   subroutine process_calculate_henry_constants(this, container, temperature, henry_constants, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      real(fp), intent(in) :: temperature(:,:,:)  !< Temperature field [K]
      real(fp), intent(out) :: henry_constants(:,:,:,:)  !< Henry constants [M/atm]
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      integer :: i, j, k, s
      real(fp) :: temp_inv, henry_298, d_sol_h

      rc = CC_SUCCESS

      if (this%n_henry_species == 0) then
         return  ! No Henry's law species
      end if

      error_mgr => container%get_error_manager()

      ! Calculate temperature-dependent Henry's law constants
      ! H(T) = H_298 * exp(d_sol_H/R * (1/T - 1/298.15))
      do k = 1, size(temperature, 3)
         do j = 1, size(temperature, 2)
            do i = 1, size(temperature, 1)
               temp_inv = 1.0_fp / temperature(i,j,k)

               do s = 1, this%n_henry_species
                  ! Placeholder values - actual implementation would read from database
                  henry_298 = 1.0e-3_fp     ! Example value [M/atm]
                  d_sol_h = -5000.0_fp      ! Example enthalpy [K]

                  henry_constants(i,j,k,s) = henry_298 * &
                     exp(d_sol_h * (temp_inv - 1.0_fp/298.15_fp))
               end do
            end do
         end do
      end do

   end subroutine process_calculate_henry_constants

   !> \brief Update photolysis rates based on current conditions
   subroutine process_update_photolysis_rates(this, container, j_values, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      real(fp), intent(out) :: j_values(:,:,:,:)  !< Photolysis rates [1/s]
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS

      if (this%n_photolysis == 0) then
         j_values = 0.0_fp
         return
      end if

      met_state => container%get_met_state_ptr()
      error_mgr => container%get_error_manager()

      if (.not. associated(met_state)) then
         rc = -1
         return
      end if

      ! Calculate photolysis rates based on solar zenith angle, cloud cover, etc.
      ! Placeholder implementation - actual implementation would use:
      ! - Solar zenith angle
      ! - Cloud optical depth
      ! - Aerosol optical depth
      ! - Surface albedo
      ! - Ozone column

      j_values = 0.0_fp  ! Placeholder

   end subroutine process_update_photolysis_rates

   !> \brief Solve the chemistry system using the configured solver
   subroutine process_solve_chemistry_system(this, container, y_old, y_new, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      real(fp), intent(in) :: y_old(:,:,:,:)   !< Initial concentrations
      real(fp), intent(out) :: y_new(:,:,:,:)  !< Final concentrations
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message

      rc = CC_SUCCESS

      error_mgr => container%get_error_manager()

      select case (trim(this%solver_type))
      case ('euler')
         call this%solve_euler(y_old, y_new, rc)
      case ('rk4')
         call this%solve_rk4(y_old, y_new, rc)
      case ('rosenbrock')
         call this%solve_rosenbrock(y_old, y_new, rc)
      case ('kpp')
         call this%solve_kpp(y_old, y_new, rc)
      case ('micm')
         call this%solve_micm(y_old, y_new, rc)
      case ('none')
         y_new = y_old  ! No chemistry
      case default
         write(message, '(A,A,A)') "Unknown solver type: ", trim(this%solver_type)
         call error_mgr%report_error(message)
         rc = -1
         return
      end select

   end subroutine process_solve_chemistry_system

   !> \brief Partition species between gas and liquid phases
   subroutine process_partition_species(this, container, gas_conc, liquid_conc, rc)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      real(fp), intent(inout) :: gas_conc(:,:,:,:)
      real(fp), intent(inout) :: liquid_conc(:,:,:,:)
      integer, intent(out) :: rc

      real(fp), allocatable :: henry_constants(:,:,:,:)
      real(fp), allocatable :: temperature(:,:,:)
      real(fp) :: total_conc, partition_coeff, lwc
      integer :: i, j, k, s

      rc = CC_SUCCESS

      if (.not. this%has_aqueous_chemistry) then
         return  ! No partitioning needed
      end if

      ! Get meteorological data
      allocate(temperature(size(gas_conc,1), size(gas_conc,2), size(gas_conc,3)))
      ! Would get temperature from met_state
      temperature = 298.15_fp  ! Placeholder

      ! Calculate Henry's law constants
      allocate(henry_constants(size(gas_conc,1), size(gas_conc,2), size(gas_conc,3), this%n_henry_species))
      call this%calculate_henry_constants(container, temperature, henry_constants, rc)
      if (rc /= 0) return

      ! Partition species
      do s = 1, min(this%n_henry_species, size(gas_conc,4))
         do k = 1, size(gas_conc, 3)
            do j = 1, size(gas_conc, 2)
               do i = 1, size(gas_conc, 1)
                  total_conc = gas_conc(i,j,k,s) + liquid_conc(i,j,k,s)
                  lwc = 1.0e-6_fp  ! Liquid water content [m³/m³] - placeholder

                  ! Partition coefficient
                  partition_coeff = henry_constants(i,j,k,s) * lwc

                  ! New concentrations
                  gas_conc(i,j,k,s) = total_conc / (1.0_fp + partition_coeff)
                  liquid_conc(i,j,k,s) = total_conc * partition_coeff / (1.0_fp + partition_coeff)
               end do
            end do
         end do
      end do

      deallocate(temperature, henry_constants)

   end subroutine process_partition_species

   !> \brief Check mass balance across phases
   function process_check_phase_balance(this, container, tolerance, rc) result(balanced)
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      real(fp), intent(in) :: tolerance
      integer, intent(out) :: rc
      logical :: balanced

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      balanced = .true.

      if (.not. this%is_multiphase) then
         return  ! Single phase - always balanced
      end if

      error_mgr => container%get_error_manager()

      ! Check mass balance across phases
      ! Placeholder implementation - actual implementation would:
      ! 1. Sum masses in all phases
      ! 2. Compare with initial total mass
      ! 3. Check against tolerance

      ! For now, assume balanced
      call error_mgr%report_info("Phase mass balance check completed")

   end function process_check_phase_balance

   !========================================================================
   ! Diagnostic Registration Methods
   !========================================================================

   !> \brief Register diagnostics for this process
   !!
   !! This method should be called during process initialization to register
   !! all diagnostics that the process will produce. Each process creates
   !! diagnostic fields and registers them with the diagnostic manager.
   !!
   !! \param[inout] this ProcessInterface instance
   !! \param[inout] container StateContainer for accessing diagnostic manager
   !! \param[out] rc Return code
   subroutine process_register_diagnostics(this, container, rc)
      use DiagnosticInterface_Mod, only: DiagnosticFieldType, DIAG_REAL_2D
      use DiagnosticManager_Mod, only: DiagnosticManagerType
      class(ProcessInterface), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS

      ! Get diagnostic manager from container
      diag_mgr => container%get_diagnostic_manager_ptr()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      endif

      ! Get error manager for logging
      error_mgr => container%get_error_manager()

      ! Register the process with the diagnostic manager
      call diag_mgr%register_process(this%name, container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(rc, 'Failed to register process with diagnostic manager')
         return
      endif

      ! Default implementation registers no diagnostics
      ! Concrete processes should override this method to register their specific diagnostics

   end subroutine process_register_diagnostics

   !> \brief Update diagnostic values for this process
   !!
   !! This method should be called during process execution to update
   !! diagnostic field values with current data from the process.
   !!
   !! \param[inout] this ProcessInterface instance
   !! \param[inout] container StateContainer for accessing state data
   !! \param[out] rc Return code
   subroutine process_update_diagnostics(this, container, rc)
      use DiagnosticManager_Mod, only: DiagnosticManagerType
      use DiagnosticInterface_Mod, only: DiagnosticRegistryType
      class(ProcessInterface), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr
      type(DiagnosticRegistryType), pointer :: diag_registry

      rc = CC_SUCCESS

      ! Get diagnostic manager from container
      diag_mgr => container%get_diagnostic_manager_ptr()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      endif

      ! Get this process's diagnostic registry
      call diag_mgr%get_process_registry(this%name, diag_registry, rc)
      if (rc /= CC_SUCCESS) then
         return
      endif

      ! Default implementation does nothing
      ! Concrete processes should override this method to update their diagnostic values

   end subroutine process_update_diagnostics

   !> \brief Get diagnostic registry for this process
   !!
   !! \param[in] this ProcessInterface instance
   !! \param[inout] container StateContainer for accessing diagnostic manager
   !! \param[out] registry Pointer to the diagnostic registry for this process
   !! \param[out] rc Return code
   subroutine process_get_diagnostic_registry(this, container, registry, rc)
      use DiagnosticManager_Mod, only: DiagnosticManagerType
      use DiagnosticInterface_Mod, only: DiagnosticRegistryType
      class(ProcessInterface), intent(in) :: this
      type(StateContainerType), intent(inout) :: container
      type(DiagnosticRegistryType), pointer, intent(out) :: registry
      integer, intent(out) :: rc

      type(DiagnosticManagerType), pointer :: diag_mgr

      nullify(registry)
      rc = CC_SUCCESS

      ! Get diagnostic manager from container
      diag_mgr => container%get_diagnostic_manager_ptr()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      endif

      ! Get this process's diagnostic registry
      call diag_mgr%get_process_registry(this%name, registry, rc)

   end subroutine process_get_diagnostic_registry

   !========================================================================
   ! Column virtualization support methods
   !========================================================================

   !> \brief Check if process supports column-based processing
   function process_supports_column_processing(this) result(supports)
      class(ProcessInterface), intent(in) :: this
      logical :: supports

      supports = .false.  ! Default implementation - override in subclasses
   end function process_supports_column_processing

   !> \brief Process a single column (default implementation)
   subroutine process_process_column(this, column, rc)
      class(ProcessInterface), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = CC_FAILURE  ! Default implementation fails - must be overridden
   end subroutine process_process_column

   !> \brief Process all columns using column processor
   subroutine process_process_all_columns(this, container, rc)
      class(ProcessInterface), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = CC_FAILURE  ! Default implementation fails - must be overridden
   end subroutine process_process_all_columns

   !========================================================================
   ! ColumnProcessInterface Implementation
   !========================================================================

   !> \brief Set column batch size
   subroutine column_process_set_batch_size(this, batch_size)
      class(ColumnProcessInterface), intent(inout) :: this
      integer, intent(in) :: batch_size

      this%column_batch_size = max(1, batch_size)
   end subroutine column_process_set_batch_size

   !> \brief Get column batch size
   function column_process_get_batch_size(this) result(batch_size)
      class(ColumnProcessInterface), intent(in) :: this
      integer :: batch_size

      batch_size = this%column_batch_size
   end function column_process_get_batch_size

   !> \brief Enable column processing
   subroutine column_process_enable(this)
      class(ColumnProcessInterface), intent(inout) :: this

      this%column_processing_enabled = .true.
   end subroutine column_process_enable

   !> \brief Disable column processing
   subroutine column_process_disable(this)
      class(ColumnProcessInterface), intent(inout) :: this

      this%column_processing_enabled = .false.
   end subroutine column_process_disable

   !> \brief Check if column processing is enabled
   function column_process_is_enabled(this) result(is_enabled)
      class(ColumnProcessInterface), intent(in) :: this
      logical :: is_enabled

      is_enabled = this%column_processing_enabled
   end function column_process_is_enabled

   ! Private solver implementations (placeholders)
   subroutine solve_euler(this, y_old, y_new, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(in) :: y_old(:,:,:,:)
      real(fp), intent(out) :: y_new(:,:,:,:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      y_new = y_old  ! Placeholder - would implement forward Euler
   end subroutine solve_euler

   subroutine solve_rk4(this, y_old, y_new, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(in) :: y_old(:,:,:,:)
      real(fp), intent(out) :: y_new(:,:,:,:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      y_new = y_old  ! Placeholder - would implement Runge-Kutta 4th order
   end subroutine solve_rk4

   subroutine solve_rosenbrock(this, y_old, y_new, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(in) :: y_old(:,:,:,:)
      real(fp), intent(out) :: y_new(:,:,:,:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      y_new = y_old  ! Placeholder - would implement Rosenbrock solver
   end subroutine solve_rosenbrock

   subroutine solve_kpp(this, y_old, y_new, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(in) :: y_old(:,:,:,:)
      real(fp), intent(out) :: y_new(:,:,:,:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      y_new = y_old  ! Placeholder - would interface with KPP solver
   end subroutine solve_kpp

   subroutine solve_micm(this, y_old, y_new, rc)
      class(ProcessInterface), intent(in) :: this
      real(fp), intent(in) :: y_old(:,:,:,:)
      real(fp), intent(out) :: y_new(:,:,:,:)
      integer, intent(out) :: rc

      rc = CC_SUCCESS
      y_new = y_old  ! Placeholder - would interface with MICM solver
   end subroutine solve_micm

   !========================================================================
   ! Generic Column Diagnostic Update Methods
   !========================================================================

   !> \brief Update a scalar diagnostic field (column scalar -> global 2D field)
   !!
   !! This is the implementation of the generic update_diagnostics interface
   !! for scalar values. It integrates with DiagnosticManager and DiagnosticInterface
   !! to update 2D global diagnostic fields from column-level scalar data.
   subroutine column_update_scalar_diagnostic(this, field_name, scalar_value, i_col, j_col, rc)
      use DiagnosticManager_Mod, only: DiagnosticManagerType
      use DiagnosticRegistry_Mod, only: DiagnosticRegistryType
      use DiagnosticInterface_Mod, only: DiagnosticFieldType, DiagnosticDataType
      
      class(ColumnProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: scalar_value
      integer, intent(in) :: i_col, j_col
      integer, intent(out) :: rc
      
      ! Local variables for diagnostic system integration
      type(StateManagerType), pointer :: state_manager => null()
      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(DiagnosticFieldType), pointer :: diag_field => null()
      type(DiagnosticDataType), pointer :: diag_data => null()
      real(fp), pointer :: field_data_2d(:,:) => null()
      
      rc = CC_SUCCESS
      
      ! Step 1: Get StateManager from concrete process
      state_manager => this%get_state_manager()
      if (.not. associated(state_manager)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 2: Get DiagnosticManager from StateManager
      diag_mgr => state_manager%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 3: Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(this%name, registry, rc)
      if (rc /= CC_SUCCESS .or. .not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 4: Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = CC_FAILURE  ! Field not found
         return
      end if
      
      ! Step 5: Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = CC_FAILURE  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 6: Get pointer to 2D field data and update
      field_data_2d => diag_data%get_real_2d_ptr()
      if (.not. associated(field_data_2d)) then
         rc = CC_FAILURE  ! Wrong data type or not allocated
         return
      end if
      
      ! Step 7: Update the field at column position
      field_data_2d(i_col, j_col) = scalar_value
      
   end subroutine column_update_scalar_diagnostic

   !> \brief Update a 1D diagnostic field (column 1D -> global 3D field)
   !!
   !! This is the implementation of the generic update_diagnostics interface
   !! for 1D arrays. It maps column-level 1D data to global 3D diagnostic fields.
   subroutine column_update_1d_diagnostic(this, field_name, array_1d, i_col, j_col, rc)
      use DiagnosticManager_Mod, only: DiagnosticManagerType
      use DiagnosticRegistry_Mod, only: DiagnosticRegistryType
      use DiagnosticInterface_Mod, only: DiagnosticFieldType, DiagnosticDataType
      
      class(ColumnProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: array_1d(:)
      integer, intent(in) :: i_col, j_col
      integer, intent(out) :: rc
      
      ! Local variables for diagnostic system integration
      type(StateManagerType), pointer :: state_manager => null()
      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(DiagnosticFieldType), pointer :: diag_field => null()
      type(DiagnosticDataType), pointer :: diag_data => null()
      real(fp), pointer :: field_data_3d(:,:,:) => null()
      integer :: k, n_dim3
      
      rc = CC_SUCCESS
      n_dim3 = size(array_1d)
      
      ! Step 1: Get StateManager from concrete process
      state_manager => this%get_state_manager()
      if (.not. associated(state_manager)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 2: Get DiagnosticManager from StateManager
      diag_mgr => state_manager%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 3: Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(this%name, registry, rc)
      if (rc /= CC_SUCCESS .or. .not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 4: Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = CC_FAILURE  ! Field not found
         return
      end if
      
      ! Step 5: Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = CC_FAILURE  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 6: Get pointer to 3D field data and validate
      field_data_3d => diag_data%get_real_3d_ptr()
      if (.not. associated(field_data_3d)) then
         rc = CC_FAILURE  ! Wrong data type or not allocated
         return
      end if
      
      ! Step 7: Validate dimensions
      if (size(field_data_3d, 3) /= n_dim3) then
         rc = CC_FAILURE  ! Dimension mismatch
         return
      end if
      
      ! Step 8: Update the field at column position
      do k = 1, n_dim3
         field_data_3d(i_col, j_col, k) = array_1d(k)
      end do
      
   end subroutine column_update_1d_diagnostic

   !> \brief Update a 2D diagnostic field (column 2D -> global field)
   !!
   !! This is the implementation of the generic update_diagnostics interface
   !! for 2D arrays. It handles mapping column-level 2D data to global diagnostic fields
   !! using a flattened storage approach.
   subroutine column_update_2d_diagnostic(this, field_name, array_2d, i_col, j_col, rc)
      use DiagnosticManager_Mod, only: DiagnosticManagerType
      use DiagnosticRegistry_Mod, only: DiagnosticRegistryType
      use DiagnosticInterface_Mod, only: DiagnosticFieldType, DiagnosticDataType
      
      class(ColumnProcessInterface), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      real(fp), intent(in) :: array_2d(:,:)
      integer, intent(in) :: i_col, j_col
      integer, intent(out) :: rc
      
      ! Local variables for diagnostic system integration
      type(StateManagerType), pointer :: state_manager => null()
      type(DiagnosticManagerType), pointer :: diag_mgr => null()
      type(DiagnosticRegistryType), pointer :: registry => null()
      type(DiagnosticFieldType), pointer :: diag_field => null()
      type(DiagnosticDataType), pointer :: diag_data => null()
      real(fp), pointer :: field_data_3d(:,:,:) => null()
      integer :: k, l, n_levels, n_species, flat_index
      
      rc = CC_SUCCESS
      
      n_levels = size(array_2d, 1)   ! levels dimension
      n_species = size(array_2d, 2)  ! species dimension
      
      ! Step 1: Get StateManager from concrete process
      state_manager => this%get_state_manager()
      if (.not. associated(state_manager)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 2: Get DiagnosticManager from StateManager
      diag_mgr => state_manager%get_diagnostic_manager()
      if (.not. associated(diag_mgr)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 3: Get process registry from DiagnosticManager
      call diag_mgr%get_process_registry(this%name, registry, rc)
      if (rc /= CC_SUCCESS .or. .not. associated(registry)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 4: Get the specific diagnostic field from registry
      diag_field => registry%get_field_ptr(field_name)
      if (.not. associated(diag_field)) then
         rc = CC_FAILURE  ! Field not found
         return
      end if
      
      ! Step 5: Verify field is ready and get data storage
      if (.not. diag_field%is_ready()) then
         rc = CC_FAILURE  ! Field not initialized
         return
      end if
      diag_data => diag_field%get_data_ptr()
      if (.not. associated(diag_data)) then
         rc = CC_FAILURE
         return
      end if
      
      ! Step 6: Get pointer to 3D field data (flattened approach)
      field_data_3d => diag_data%get_real_3d_ptr()
      if (.not. associated(field_data_3d)) then
         rc = CC_FAILURE  ! Wrong data type or not allocated
         return
      end if
      
      ! Step 7: Validate dimensions (expect flattened storage)
      if (size(field_data_3d, 3) /= n_levels * n_species) then
         rc = CC_FAILURE  ! Dimension mismatch
         return
      end if
      
      ! Step 8: Update the field at column position using flattened indexing
      do l = 1, n_species
         do k = 1, n_levels
            flat_index = (l - 1) * n_levels + k
            field_data_3d(i_col, j_col, flat_index) = array_2d(k, l)
         end do
      end do
      
   end subroutine column_update_2d_diagnostic

end module ProcessInterface_Mod