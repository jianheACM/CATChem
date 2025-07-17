!> \file StateManager_Mod.F90
!! \brief Unified state management module for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides unified state management for CATChem, including:
!! - State container and basic management
!! - State validation utilities
!! - State type constants and utilities
!! - Direct field access functions
!!
!! Most complex lifecycle management is delegated to CATChemCore_Mod.
!!
module StateManager_Mod
   use precision_mod, only: fp
   use error_mod, only: CC_SUCCESS, CC_FAILURE, ErrorManagerType
   use ConfigManager_Mod, only: ConfigDataType
   use MetState_Mod, only: MetStateType
   use ChemState_Mod, only: ChemStateType

   implicit none
   private

   ! Public types
   public :: StateManagerType
   public :: StateValidatorUtilsType

   ! Public enumerations and constants
   public :: STATE_TYPE_MET, STATE_TYPE_CHEM, STATE_TYPE_EMIS, STATE_TYPE_DIAG
   public :: STATE_STATUS_UNINITIALIZED, STATE_STATUS_INITIALIZED, STATE_STATUS_VALID, STATE_STATUS_ERROR

   ! Public utility procedures
   public :: get_state_type_name, allocate_met_field

   !=========================================================================
   ! Constants and Enumerations
   !=========================================================================

   !> State type enumeration
   integer, parameter :: STATE_TYPE_MET = 1
   integer, parameter :: STATE_TYPE_CHEM = 2
   integer, parameter :: STATE_TYPE_EMIS = 3
   integer, parameter :: STATE_TYPE_DIAG = 4
   integer, parameter :: STATE_TYPE_CONFIG = 5
   integer, parameter :: STATE_TYPE_GRID = 6

   !> State status enumeration
   integer, parameter :: STATE_STATUS_UNINITIALIZED = 0
   integer, parameter :: STATE_STATUS_INITIALIZED = 1
   integer, parameter :: STATE_STATUS_VALID = 2
   integer, parameter :: STATE_STATUS_ERROR = -1

   !=========================================================================
   ! Types
   !=========================================================================

   !> \brief Simplified state container managed by CATChemCore
   !!
   !! This type provides basic state object management. Complex lifecycle
   !! operations are handled by CATChemCore_Mod.
   !!
   type :: StateManagerType
      private

      ! Core state objects
      type(ConfigDataType), allocatable :: config      !< Configuration state
      type(MetStateType),   allocatable :: met_state   !< Meteorological fields
      type(ChemStateType),  allocatable :: chem_state  !< Chemical species concentrations
      type(ErrorManagerType)            :: error_mgr   !< Error manager

      ! Simple metadata
      logical :: is_initialized = .false.              !< Initialization status
      character(len=256) :: name = ''                  !< Container name

   contains
      ! Basic lifecycle (called by CATChemCore)
      procedure :: init => manager_init
      procedure :: cleanup => manager_cleanup
      procedure :: finalize => manager_finalize
      procedure :: is_ready => manager_is_ready

      ! State object accessors
      procedure :: get_config_ptr => manager_get_config_ptr
      procedure :: get_met_state_ptr => manager_get_met_state_ptr
      procedure :: get_chem_state_ptr => manager_get_chem_state_ptr
      procedure :: get_error_manager => manager_get_error_manager

      ! Utilities
      procedure :: set_name => manager_set_name
      procedure :: print_info => manager_print_info
      procedure :: get_memory_usage => manager_get_memory_usage

   end type StateManagerType

   !> \brief State validation utilities
   !!
   !! This type provides common validation routines that can be used
   !! across different state types.
   !!
   type :: StateValidatorUtilsType

   contains
      procedure :: validate_dimensions => validator_validate_dimensions
      procedure :: validate_bounds => validator_validate_bounds
      procedure :: validate_consistency => validator_validate_consistency
      procedure :: check_nan_values => validator_check_nan_values
      procedure :: check_negative_values => validator_check_negative_values

   end type StateValidatorUtilsType

contains

   !=========================================================================
   ! StateManagerType Implementation
   !=========================================================================

   !> \brief Initialize the state manager (called by CATChemCore)
   subroutine manager_init(this, name, rc)
      class(StateManagerType), intent(inout) :: this
      character(len=*), optional, intent(in) :: name
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Set manager name
      if (present(name)) then
         this%name = trim(name)
      else
         this%name = 'StateManager'
      endif

      ! Allocate state objects
      if (.not. allocated(this%config)) allocate(this%config)
      if (.not. allocated(this%met_state)) allocate(this%met_state)
      if (.not. allocated(this%chem_state)) allocate(this%chem_state)

      this%is_initialized = .true.

   end subroutine manager_init

   !> \brief Clean up the state manager
   subroutine manager_cleanup(this, rc)
      class(StateManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate state objects
      if (allocated(this%config)) deallocate(this%config)
      if (allocated(this%met_state)) deallocate(this%met_state)
      if (allocated(this%chem_state)) deallocate(this%chem_state)

      this%is_initialized = .false.
      this%name = ''

   end subroutine manager_cleanup

   !> \brief Check if manager is ready
   function manager_is_ready(this) result(ready)
      class(StateManagerType), intent(in) :: this
      logical :: ready

      ready = this%is_initialized .and. &
              allocated(this%config) .and. &
              allocated(this%met_state) .and. &
              allocated(this%chem_state)
   end function manager_is_ready

   !> \brief Get pointer to config for modification
   function manager_get_config_ptr(this) result(config_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(ConfigDataType), pointer :: config_ptr

      if (allocated(this%config)) then
         config_ptr => this%config
      else
         nullify(config_ptr)
      endif
   end function manager_get_config_ptr

   !> \brief Get pointer to met state for modification
   function manager_get_met_state_ptr(this) result(met_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(MetStateType), pointer :: met_ptr

      if (allocated(this%met_state)) then
         met_ptr => this%met_state
      else
         nullify(met_ptr)
      endif
   end function manager_get_met_state_ptr

   !> \brief Get pointer to chem state for modification
   function manager_get_chem_state_ptr(this) result(chem_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(ChemStateType), pointer :: chem_ptr

      if (allocated(this%chem_state)) then
         chem_ptr => this%chem_state
      else
         nullify(chem_ptr)
      endif
   end function manager_get_chem_state_ptr

   !> \brief Get pointer to error manager
   function manager_get_error_manager(this) result(error_mgr_ptr)
      class(StateManagerType), intent(inout), target :: this
      type(ErrorManagerType), pointer :: error_mgr_ptr

      error_mgr_ptr => this%error_mgr
   end function manager_get_error_manager

   !> \brief Set manager name
   subroutine manager_set_name(this, name)
      class(StateManagerType), intent(inout) :: this
      character(len=*), intent(in) :: name

      this%name = trim(name)
   end subroutine manager_set_name

   !> \brief Print manager information
   subroutine manager_print_info(this)
      class(StateManagerType), intent(in) :: this

      write(*,'(A)') '=== StateManager Information ==='
      write(*,'(A,A)') 'Name: ', trim(this%name)
      write(*,'(A,L1)') 'Initialized: ', this%is_initialized
      write(*,'(A,L1)') 'Config allocated: ', allocated(this%config)
      write(*,'(A,L1)') 'Met state allocated: ', allocated(this%met_state)
      write(*,'(A,L1)') 'Chem state allocated: ', allocated(this%chem_state)
      write(*,'(A)') '================================='

   end subroutine manager_print_info

   !> \brief Finalize the state manager (alias for cleanup)
   subroutine manager_finalize(this, rc)
      class(StateManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%cleanup(rc)
   end subroutine manager_finalize

   !> \brief Get approximate memory usage in bytes
   function manager_get_memory_usage(this) result(memory_bytes)
      class(StateManagerType), intent(in) :: this
      integer(8) :: memory_bytes

      ! Simplified calculation - real implementation would query each state object
      memory_bytes = 0_8

      if (allocated(this%config)) memory_bytes = memory_bytes + 1024_8
      if (allocated(this%met_state)) memory_bytes = memory_bytes + 102400_8
      if (allocated(this%chem_state)) memory_bytes = memory_bytes + 1048576_8
   end function manager_get_memory_usage

   !=========================================================================
   ! StateValidatorUtilsType Implementation
   !=========================================================================

   !> \brief Validate array dimensions
   subroutine validator_validate_dimensions(this, array_shape, expected_shape, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      integer, intent(in) :: array_shape(:)
      integer, intent(in) :: expected_shape(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (size(array_shape) /= size(expected_shape)) then
         rc = CC_FAILURE
         return
      endif

      do i = 1, size(array_shape)
         if (array_shape(i) /= expected_shape(i)) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_validate_dimensions

   !> \brief Validate value bounds
   subroutine validator_validate_bounds(this, values, min_val, max_val, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      real(fp), intent(in) :: min_val, max_val
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      do i = 1, size(values)
         if (values(i) < min_val .or. values(i) > max_val) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_validate_bounds

   !> \brief Validate consistency between related arrays
   subroutine validator_validate_consistency(this, array1, array2, tolerance, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: array1(:), array2(:)
      real(fp), intent(in) :: tolerance
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      if (size(array1) /= size(array2)) then
         rc = CC_FAILURE
         return
      endif

      do i = 1, size(array1)
         if (abs(array1(i) - array2(i)) > tolerance) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_validate_consistency

   !> \brief Check for NaN values
   subroutine validator_check_nan_values(this, values, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      do i = 1, size(values)
         if (values(i) /= values(i)) then  ! NaN check
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_check_nan_values

   !> \brief Check for negative values where they shouldn't exist
   subroutine validator_check_negative_values(this, values, rc)
      class(StateValidatorUtilsType), intent(in) :: this
      real(fp), intent(in) :: values(:)
      integer, intent(out) :: rc

      integer :: i

      rc = CC_SUCCESS

      do i = 1, size(values)
         if (values(i) < 0.0_fp) then
            rc = CC_FAILURE
            return
         endif
      enddo

   end subroutine validator_check_negative_values

   !=========================================================================
   ! Utility Functions
   !=========================================================================

   !> \brief Get state type name from enumeration
   function get_state_type_name(state_type) result(name)
      integer, intent(in) :: state_type
      character(len=32) :: name

      select case (state_type)
         case (STATE_TYPE_MET)
            name = 'Meteorology'
         case (STATE_TYPE_CHEM)
            name = 'Chemistry'
         case (STATE_TYPE_EMIS)
            name = 'Emissions'
         case (STATE_TYPE_DIAG)
            name = 'Diagnostics'
         case (STATE_TYPE_CONFIG)
            name = 'Configuration'
         case (STATE_TYPE_GRID)
            name = 'Grid'
         case default
            name = 'Unknown'
      end select

   end function get_state_type_name

   !> \brief Allocate a specific field in a MetStateType object
   !!
   !! Direct utility function for field allocation in MetStateType.
   !! For more complex state management, use CATChemCore_Mod.
   !!
   !! \param[inout] met_state   MetStateType object
   !! \param[in]    field_name  Name of the field to allocate (or 'ALL')
   !! \param[out]   rc          Return code (CC_SUCCESS or error code)
   subroutine allocate_met_field(met_state, field_name, rc)
      class(MetStateType), intent(inout) :: met_state
      character(len=*), intent(in) :: field_name
      integer, intent(out) :: rc

      call met_state%allocate_field(field_name, rc)
   end subroutine allocate_met_field

end module StateManager_Mod
