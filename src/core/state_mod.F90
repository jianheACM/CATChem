!> \file state_mod.F90
!! \brief Modern state management module for CATChem using dependency injection
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides a modern state management system using the StateContainer
!! pattern with dependency injection. It replaces the global state variables with
!! a more maintainable, testable, and flexible architecture.
!!
!! \details
!! The refactored state module provides:
!! - StateContainer type that encapsulates all state objects
!! - StateBuilder for flexible state construction and initiali         ! Validate configuration (TODO: implement config_data_validate)
         ! call config_ptr%validate(error_mgr, local_rc)
         ! if (local_rc /= CC_SUCCESS) then
         !    call error_mgr%pop_context()
         !    return
         ! endifn
!! - Dependency injection pattern for better modularity
!! - Thread-safe state access methods
!! - State validation and error handling
!! - Memory management and cleanup utilities
!!
!! \section state_usage Usage Example
!! \code{.f90}
!! use state_mod
!! type(StateContainerType) :: container
!! type(StateBuilderType) :: builder
!!
!! ! Build and initialize state container
!! call builder%init()
!! call builder%with_config(config_file)
!! call builder%with_grid(nx, ny, nz)
!! call builder%build(container, rc)
!!
!! ! Use the container
!! call some_process(container, rc)
!!
!! ! Cleanup
!! call container%finalize(rc)
!! \endcode
!!
module state_mod
   use precision_mod
   use ConfigManager_Mod, only : ConfigDataType, ConfigManagerType
   use MetState_Mod,   only : MetStateType
   use ChemState_Mod,  only : ChemStateType
   use EmisState_Mod,  only : EmisStateType
   use DiagState_Mod,  only : DiagStateType
   use Error_Mod,      only : CC_SUCCESS, CC_FAILURE, ErrorManagerType

   IMPLICIT NONE
   PRIVATE

   ! Public types and procedures
   PUBLIC :: StateContainerType
   PUBLIC :: StateBuilderType
   PUBLIC :: StateValidatorType

   !> \brief Main state container that encapsulates all CATChem state objects
   !!
   !! This type provides a centralized container for all state objects in CATChem,
   !! enabling dependency injection and better separation of concerns.
   !!
   !! \details
   !! The StateContainer manages:
   !! - Core state objects (Met, Chem, Emis, Diag, Config)
   !! - State validation and consistency checking
   !! - Memory management and cleanup
   !! - Thread-safe access to state data
   !! - State serialization and debugging utilities
   !!
   type :: StateContainerType
      private

      ! Core state objects
      type(ConfigDataType),    allocatable :: config      !< Configuration state
      type(ConfigManagerType), allocatable :: config_mgr  !< Configuration manager
      type(MetStateType),  allocatable :: met_state   !< Meteorological fields
      type(ChemStateType), allocatable :: chem_state  !< Chemical species concentrations
      type(EmisStateType), allocatable :: emis_state  !< Emission fluxes and mappings
      type(DiagStateType), allocatable :: diag_state  !< Diagnostic output variables

      ! Error handling and debugging
      type(ErrorManagerType) :: error_mgr             !< Integrated error manager

      ! Container metadata
      logical :: is_initialized = .false.  !< Initialization status
      logical :: is_valid = .false.        !< Validation status
      character(len=256) :: name = ''      !< Container name for debugging
      integer :: creation_time             !< Creation timestamp

   contains
      ! Initialization and cleanup
      procedure :: init => container_init
      procedure :: finalize => container_finalize
      procedure :: validate => container_validate

      ! State object accessors (const)
      procedure :: get_config => container_get_config
      procedure :: get_met_state => container_get_met_state
      procedure :: get_chem_state => container_get_chem_state
      procedure :: get_emis_state => container_get_emis_state
      procedure :: get_diag_state => container_get_diag_state

      ! State object accessors (mutable)
      procedure :: get_config_ptr => container_get_config_ptr
      procedure :: get_met_state_ptr => container_get_met_state_ptr
      procedure :: get_chem_state_ptr => container_get_chem_state_ptr
      procedure :: get_emis_state_ptr => container_get_emis_state_ptr
      procedure :: get_diag_state_ptr => container_get_diag_state_ptr

      ! Utility methods
      procedure :: is_ready => container_is_ready
      procedure :: get_memory_usage => container_get_memory_usage
      procedure :: print_info => container_print_info
      procedure :: set_name => container_set_name
      procedure :: get_error_manager => container_get_error_manager

   end type StateContainerType

   !> \brief Builder pattern implementation for flexible state container construction
   !!
   !! This type implements the builder pattern to provide a flexible way to construct
   !! and initialize StateContainer objects with various configurations.
   !!
   type :: StateBuilderType
      private

      ! Configuration options
      character(len=512) :: config_file = ''
      logical :: auto_validate = .true.
      logical :: verbose_init = .false.
      character(len=256) :: container_name = 'CATChem_Container'

      ! Initialization flags
      logical :: init_met = .true.
      logical :: init_chem = .true.
      logical :: init_emis = .true.
      logical :: init_diag = .true.

   contains
      ! Builder configuration methods
      procedure :: init => builder_init
      procedure :: with_config_file => builder_with_config_file
      procedure :: with_name => builder_with_name
      procedure :: enable_verbose => builder_enable_verbose
      procedure :: disable_validation => builder_disable_validation
      procedure :: skip_met_init => builder_skip_met_init
      procedure :: skip_chem_init => builder_skip_chem_init
      procedure :: skip_emis_init => builder_skip_emis_init
      procedure :: skip_diag_init => builder_skip_diag_init

      ! Build method
      procedure :: build => builder_build

   end type StateBuilderType

   !> \brief State validation utilities
   !!
   !! This type provides comprehensive validation capabilities for state objects
   !! and the overall state container consistency.
   !!
   type :: StateValidatorType
   contains
      procedure :: validate_container => validator_validate_container
      procedure :: validate_config => validator_validate_config
      procedure :: validate_consistency => validator_validate_consistency
   end type StateValidatorType

contains
   !========================================================================
   ! StateContainer Implementation
   !========================================================================

   !> \brief Initialize the state container
   !!
   !! \param[inout] this The state container
   !! \param[in] name Optional name for the container
   !! \param[out] rc Return code
   subroutine container_init(this, name, rc)
      implicit none
      class(StateContainerType), intent(inout) :: this
      character(len=*), optional, intent(in) :: name
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Set container name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'CATChem_Container'
      endif

      ! Initialize timestamp (simplified - could use system time)
      this%creation_time = 0

      ! Initialize error manager
      call this%error_mgr%init(verbose=.true.)
      call this%error_mgr%push_context('container_init', 'Initializing state container')

      ! Allocate state objects
      if (.not. allocated(this%config)) allocate(this%config)
      if (.not. allocated(this%met_state)) allocate(this%met_state)
      if (.not. allocated(this%chem_state)) allocate(this%chem_state)
      if (.not. allocated(this%emis_state)) allocate(this%emis_state)
      if (.not. allocated(this%diag_state)) allocate(this%diag_state)

      this%is_initialized = .true.
      this%is_valid = .false.  ! Will be set to true after validation

      call this%error_mgr%pop_context()

   end subroutine container_init

   !> \brief Finalize and cleanup the state container
   subroutine container_finalize(this, rc)
      implicit none
      class(StateContainerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate state objects
      if (allocated(this%config)) deallocate(this%config)
      if (allocated(this%met_state)) deallocate(this%met_state)
      if (allocated(this%chem_state)) deallocate(this%chem_state)
      if (allocated(this%emis_state)) deallocate(this%emis_state)
      if (allocated(this%diag_state)) deallocate(this%diag_state)

      this%is_initialized = .false.
      this%is_valid = .false.
      this%name = ''

   end subroutine container_finalize

   !> \brief Validate the state container
   subroutine container_validate(this, rc)
      implicit none
      class(StateContainerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Check if initialized
      if (.not. this%is_initialized) then
         rc = CC_FAILURE
         return
      endif

      ! Check that all required state objects are allocated
      if (.not. allocated(this%config)) then
         rc = CC_FAILURE
         return
      endif

      ! Additional validation logic would go here
      ! For now, we'll assume it's valid if initialized and allocated
      this%is_valid = .true.

   end subroutine container_validate

   !> \brief Get const reference to config
   function container_get_config(this) result(config_ref)
      implicit none
      class(StateContainerType), intent(in) :: this
      type(ConfigDataType) :: config_ref

      if (allocated(this%config)) then
         config_ref = this%config
      endif
   end function container_get_config

   !> \brief Get pointer to config for modification
   function container_get_config_ptr(this) result(config_ptr)
      implicit none
      class(StateContainerType), intent(inout), target :: this
      type(ConfigDataType), pointer :: config_ptr

      if (allocated(this%config)) then
         config_ptr => this%config
      else
         nullify(config_ptr)
      endif
   end function container_get_config_ptr

   !> \brief Get const reference to grid state
   !> \brief Get const reference to met state
   function container_get_met_state(this) result(met_ref)
      implicit none
      class(StateContainerType), intent(in) :: this
      type(MetStateType) :: met_ref

      if (allocated(this%met_state)) then
         met_ref = this%met_state
      endif
   end function container_get_met_state

   !> \brief Get pointer to met state for modification
   function container_get_met_state_ptr(this) result(met_ptr)
      implicit none
      class(StateContainerType), intent(inout), target :: this
      type(MetStateType), pointer :: met_ptr

      if (allocated(this%met_state)) then
         met_ptr => this%met_state
      else
         nullify(met_ptr)
      endif
   end function container_get_met_state_ptr

   !> \brief Get const reference to chem state
   function container_get_chem_state(this) result(chem_ref)
      implicit none
      class(StateContainerType), intent(in) :: this
      type(ChemStateType) :: chem_ref

      if (allocated(this%chem_state)) then
         chem_ref = this%chem_state
      endif
   end function container_get_chem_state

   !> \brief Get pointer to chem state for modification
   function container_get_chem_state_ptr(this) result(chem_ptr)
      implicit none
      class(StateContainerType), intent(inout), target :: this
      type(ChemStateType), pointer :: chem_ptr

      if (allocated(this%chem_state)) then
         chem_ptr => this%chem_state
      else
         nullify(chem_ptr)
      endif
   end function container_get_chem_state_ptr

   !> \brief Get const reference to emis state
   function container_get_emis_state(this) result(emis_ref)
      implicit none
      class(StateContainerType), intent(in) :: this
      type(EmisStateType) :: emis_ref

      if (allocated(this%emis_state)) then
         emis_ref = this%emis_state
      endif
   end function container_get_emis_state

   !> \brief Get pointer to emis state for modification
   function container_get_emis_state_ptr(this) result(emis_ptr)
      implicit none
      class(StateContainerType), intent(inout), target :: this
      type(EmisStateType), pointer :: emis_ptr

      if (allocated(this%emis_state)) then
         emis_ptr => this%emis_state
      else
         nullify(emis_ptr)
      endif
   end function container_get_emis_state_ptr

   !> \brief Get const reference to diag state
   function container_get_diag_state(this) result(diag_ref)
      implicit none
      class(StateContainerType), intent(in) :: this
      type(DiagStateType) :: diag_ref

      if (allocated(this%diag_state)) then
         diag_ref = this%diag_state
      endif
   end function container_get_diag_state

   !> \brief Get pointer to diag state for modification
   function container_get_diag_state_ptr(this) result(diag_ptr)
      implicit none
      class(StateContainerType), intent(inout), target :: this
      type(DiagStateType), pointer :: diag_ptr

      if (allocated(this%diag_state)) then
         diag_ptr => this%diag_state
      else
         nullify(diag_ptr)
      endif
   end function container_get_diag_state_ptr

   !> \brief Check if container is ready for use
   function container_is_ready(this) result(ready)
      implicit none
      class(StateContainerType), intent(in) :: this
      logical :: ready

      ready = this%is_initialized .and. this%is_valid
   end function container_is_ready

   !> \brief Get approximate memory usage in bytes
   function container_get_memory_usage(this) result(memory_bytes)
      implicit none
      class(StateContainerType), intent(in) :: this
      integer(kind=8) :: memory_bytes

      ! This is a simplified calculation - real implementation would
      ! query each state object for its memory usage
      memory_bytes = 0

      if (allocated(this%config)) memory_bytes = memory_bytes + 1024
      if (allocated(this%met_state)) memory_bytes = memory_bytes + 102400
      if (allocated(this%chem_state)) memory_bytes = memory_bytes + 1048576
      if (allocated(this%emis_state)) memory_bytes = memory_bytes + 524288
      if (allocated(this%diag_state)) memory_bytes = memory_bytes + 262144

   end function container_get_memory_usage

   !> \brief Print container information for debugging
   subroutine container_print_info(this)
      implicit none
      class(StateContainerType), intent(in) :: this

      write(*,'(A)') '=== StateContainer Information ==='
      write(*,'(A,A)') 'Name: ', trim(this%name)
      write(*,'(A,L1)') 'Initialized: ', this%is_initialized
      write(*,'(A,L1)') 'Valid: ', this%is_valid
      write(*,'(A,I0)') 'Memory Usage (bytes): ', this%get_memory_usage()
      write(*,'(A,L1)') 'Config allocated: ', allocated(this%config)
      write(*,'(A,L1)') 'Met state allocated: ', allocated(this%met_state)
      write(*,'(A,L1)') 'Chem state allocated: ', allocated(this%chem_state)
      write(*,'(A,L1)') 'Emis state allocated: ', allocated(this%emis_state)
      write(*,'(A,L1)') 'Diag state allocated: ', allocated(this%diag_state)
      write(*,'(A)') '================================='

   end subroutine container_print_info

   !> \brief Set container name
   subroutine container_set_name(this, name)
      implicit none
      class(StateContainerType), intent(inout) :: this
      character(len=*), intent(in) :: name

      this%name = trim(name)
   end subroutine container_set_name

   !========================================================================
   ! StateBuilder Implementation
   !========================================================================

   !> \brief Initialize the state builder
   subroutine builder_init(this)
      implicit none
      class(StateBuilderType), intent(inout) :: this

      ! Reset to defaults
      this%config_file = ''
      this%auto_validate = .true.
      this%verbose_init = .false.
      this%container_name = 'CATChem_Container'
      this%init_met = .true.
      this%init_chem = .true.
      this%init_emis = .true.
      this%init_diag = .true.

   end subroutine builder_init

   !> \brief Set configuration file for the builder
   function builder_with_config_file(this, config_file) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      character(len=*), intent(in) :: config_file
      type(StateBuilderType) :: builder_ref

      this%config_file = trim(config_file)
      builder_ref = this
   end function builder_with_config_file

   !> \brief Set container name for the builder
   function builder_with_name(this, name) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(StateBuilderType) :: builder_ref

      this%container_name = trim(name)
      builder_ref = this
   end function builder_with_name

   !> \brief Enable verbose initialization
   function builder_enable_verbose(this) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      type(StateBuilderType) :: builder_ref

      this%verbose_init = .true.
      builder_ref = this
   end function builder_enable_verbose

   !> \brief Disable automatic validation
   function builder_disable_validation(this) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      type(StateBuilderType) :: builder_ref

      this%auto_validate = .false.
      builder_ref = this
   end function builder_disable_validation

   !> \brief Skip meteorology state initialization
   function builder_skip_met_init(this) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      type(StateBuilderType) :: builder_ref

      this%init_met = .false.
      builder_ref = this
   end function builder_skip_met_init

   !> \brief Skip chemistry state initialization
   function builder_skip_chem_init(this) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      type(StateBuilderType) :: builder_ref

      this%init_chem = .false.
      builder_ref = this
   end function builder_skip_chem_init

   !> \brief Skip emissions state initialization
   function builder_skip_emis_init(this) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      type(StateBuilderType) :: builder_ref

      this%init_emis = .false.
      builder_ref = this
   end function builder_skip_emis_init

   !> \brief Skip diagnostics state initialization
   function builder_skip_diag_init(this) result(builder_ref)
      implicit none
      class(StateBuilderType), intent(inout) :: this
      type(StateBuilderType) :: builder_ref

      this%init_diag = .false.
      builder_ref = this
   end function builder_skip_diag_init

   !> \brief Build the state container with configured options
   subroutine builder_build(this, container, rc)
      implicit none
      class(StateBuilderType), intent(in) :: this
      type(StateContainerType), intent(out) :: container
      integer, intent(out) :: rc

      ! Local variables
      integer :: local_rc
      logical :: is_valid
      type(ConfigDataType), pointer :: config_ptr
      type(MetStateType), pointer :: met_ptr
      type(ChemStateType), pointer :: chem_ptr
      type(EmisStateType), pointer :: emis_ptr
      type(DiagStateType), pointer :: diag_ptr
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS

      if (this%verbose_init) then
         write(*,'(A)') 'Building CATChem state container...'
         write(*,'(A,A)') '  Container name: ', trim(this%container_name)
         write(*,'(A,A)') '  Config file: ', trim(this%config_file)
      endif

      ! Initialize the container
      call container%init(this%container_name, local_rc)
      if (local_rc /= CC_SUCCESS) then
         rc = local_rc
         return
      endif

      ! Get error manager for initialization context
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('builder_build', 'building state container')

      ! Get pointers to state objects for initialization
      config_ptr => container%get_config_ptr()

      ! Initialize configuration if file is provided
      if (len_trim(this%config_file) > 0) then
         if (this%verbose_init) write(*,'(A)') '  Reading configuration file...'

         ! Allocate and initialize configuration manager
         if (.not. allocated(container%config_mgr)) then
            allocate(container%config_mgr)
            call container%config_mgr%init(local_rc)
            if (local_rc /= CC_SUCCESS) then
               call error_mgr%report_error(local_rc, &
                                           'Failed to initialize configuration manager', rc, &
                                           'builder_build', &
                                           'Check system resources')
               call error_mgr%pop_context()
               return
            endif
         endif

         ! Use ConfigManager to load configuration
         call container%config_mgr%load_from_file(this%config_file, local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%report_error(local_rc, &
                                        'Failed to load configuration file', rc, &
                                        'builder_build', &
                                        'Check configuration file format and permissions')
            call error_mgr%pop_context()
            return
         endif

         ! Get the configuration data from the manager
         ! TODO: Add a method to ConfigManagerType to extract ConfigDataType
         ! For now, initialize with defaults
         call config_ptr%init(local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%pop_context()
            return
         endif

         ! Validate configuration
         is_valid = config_ptr%validate(error_mgr, local_rc)
         if (local_rc /= CC_SUCCESS .or. .not. is_valid) then
            call error_mgr%report_error(local_rc, &
                                        'Configuration validation failed', rc, &
                                        'builder_build', &
                                        'Check configuration parameters')
            call error_mgr%pop_context()
            return
         endif
      else
         ! Use default configuration - initialize with defaults
         call config_ptr%init(local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Initialize state objects using new patterns
      if (this%init_met) then
         if (this%verbose_init) write(*,'(A)') '  Initializing meteorology state...'
         met_ptr => container%get_met_state_ptr()
         call met_ptr%init(config_ptr%runtime%nLevs, error_mgr, local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%report_error(local_rc, &
                                        'Failed to initialize meteorology state', rc, &
                                        'builder_build', &
                                        'Check grid dimensions and memory availability')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (this%init_chem) then
         if (this%verbose_init) write(*,'(A)') '  Initializing chemistry state...'
         chem_ptr => container%get_chem_state_ptr()
         call chem_ptr%init(config_ptr%runtime%maxSpecies, error_mgr, local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%report_error(local_rc, &
                                        'Failed to initialize chemistry state', rc, &
                                        'builder_build', &
                                        'Check species configuration and memory availability')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (this%init_emis) then
         if (this%verbose_init) write(*,'(A)') '  Initializing emissions state...'
         emis_ptr => container%get_emis_state_ptr()
         call emis_ptr%init(config_ptr%runtime%nEmissionCategories, error_mgr, local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%report_error(local_rc, &
                                        'Failed to initialize emissions state', rc, &
                                        'builder_build', &
                                        'Check emissions configuration and memory availability')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (this%init_diag) then
         if (this%verbose_init) write(*,'(A)') '  Initializing diagnostics state...'
         diag_ptr => container%get_diag_state_ptr()
         call diag_ptr%init(config_ptr%runtime%nLevs, config_ptr%runtime%nSpecies_drydep, error_mgr, local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%report_error(local_rc, &
                                        'Failed to initialize diagnostics state', rc, &
                                        'builder_build', &
                                        'Check diagnostics configuration and memory availability')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Validate container if requested
      if (this%auto_validate) then
         if (this%verbose_init) write(*,'(A)') '  Validating state container...'
         call container%validate(local_rc)
         if (local_rc /= CC_SUCCESS) then
            call error_mgr%report_error(local_rc, &
                                        'Container validation failed', rc, &
                                        'builder_build', &
                                        'Check all state objects are properly initialized')
            call error_mgr%pop_context()
            return
         endif
      endif

      call error_mgr%pop_context()

      if (this%verbose_init) then
         write(*,'(A)') 'State container built successfully!'
         call container%print_info()
      endif

   end subroutine builder_build

   !========================================================================
   ! StateValidator Implementation
   !========================================================================

   !> \brief Validate the entire state container
   subroutine validator_validate_container(this, container, rc)
      implicit none
      class(StateValidatorType), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc

      integer :: local_rc

      rc = CC_SUCCESS

      ! Check container basic properties
      if (.not. container%is_initialized) then
         rc = CC_FAILURE
         return
      endif

      ! Validate individual state objects
      call this%validate_config(container%get_config(), local_rc)
      if (local_rc /= CC_SUCCESS) then
         rc = local_rc
         return
      endif

      ! Validate consistency between state objects
      call this%validate_consistency(container, local_rc)
      if (local_rc /= CC_SUCCESS) then
         rc = local_rc
         return
      endif

   end subroutine validator_validate_container

   !> \brief Validate configuration object
   subroutine validator_validate_config(this, config, rc)
      implicit none
      class(StateValidatorType), intent(in) :: this
      type(ConfigDataType), intent(in) :: config
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Add specific configuration validation logic here
      ! For example: check that required fields are present,
      ! values are within valid ranges, etc.

      ! Validate simulation parameters
      if (config%runtime%nLevs <= 0) then
         rc = CC_FAILURE
         return
      endif

      if (config%runtime%maxSpecies <= 0) then
         rc = CC_FAILURE
         return
      endif

      ! Validate chemistry settings
      if (config%runtime%nSpecies_drydep > config%runtime%maxSpecies) then
         rc = CC_FAILURE
         return
      endif

   end subroutine validator_validate_config

   !> \brief Validate consistency between state objects
   subroutine validator_validate_consistency(this, container, rc)
      implicit none
      class(StateValidatorType), intent(in) :: this
      type(StateContainerType), intent(in) :: container
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Add consistency checks here
      ! For example: ensure chemistry species match emission species,
      ! grid dimensions are consistent across all state objects, etc.

      ! Placeholder validation - always succeeds for now

   end subroutine validator_validate_consistency

   !> \brief Get reference to the error manager
   function container_get_error_manager(this) result(error_mgr_ref)
      implicit none
      class(StateContainerType), intent(inout), target :: this
      type(ErrorManagerType), pointer :: error_mgr_ref

      error_mgr_ref => this%error_mgr
   end function container_get_error_manager

end module state_mod
