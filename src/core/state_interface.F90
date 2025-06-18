! \file state_interface.F90
!! \brief Abstract interface for all state modules
!!
!! This module defines the standard interface that all state modules
!! (MetState, ChemState, EmisState, DiagState) must implement for consistency.
!!
!! \ingroup core_modules
module StateInterface_Mod
   use Precision_Mod
   use Error_Mod

   implicit none
   private

   !> \brief Abstract interface for state module initialization
   !!
   !! All state modules must implement this interface with consistent
   !! parameter names and error handling patterns.
   abstract interface
      !> \brief Initialize state module
      !!
      !! \param[inout] this The state object to initialize
      !! \param[in] error_mgr Error manager for context and error reporting
      !! \param[out] rc Return code
      subroutine state_init_interface(this, error_mgr, rc)
         import :: ErrorManagerType
         class(*), intent(inout) :: this
         type(ErrorManagerType), pointer, intent(inout) :: error_mgr
         integer, intent(out) :: rc
      end subroutine state_init_interface

      !> \brief Validate state module
      !!
      !! \param[in] this The state object to validate
      !! \param[in] error_mgr Error manager for context and error reporting
      !! \param[out] rc Return code
      subroutine state_validate_interface(this, error_mgr, rc)
         import :: ErrorManagerType
         class(*), intent(in) :: this
         type(ErrorManagerType), pointer, intent(inout) :: error_mgr
         integer, intent(out) :: rc
      end subroutine state_validate_interface

      !> \brief Clean up state module
      !!
      !! \param[inout] this The state object to clean up
      !! \param[in] error_mgr Error manager for context and error reporting
      !! \param[out] rc Return code
      subroutine state_cleanup_interface(this, error_mgr, rc)
         import :: ErrorManagerType
         class(*), intent(inout) :: this
         type(ErrorManagerType), pointer, intent(inout) :: error_mgr
         integer, intent(out) :: rc
      end subroutine state_cleanup_interface

      !> \brief Reset state module
      !!
      !! \param[inout] this The state object to reset
      !! \param[in] error_mgr Error manager for context and error reporting
      !! \param[out] rc Return code
      subroutine state_reset_interface(this, error_mgr, rc)
         import :: ErrorManagerType
         class(*), intent(inout) :: this
         type(ErrorManagerType), pointer, intent(inout) :: error_mgr
         integer, intent(out) :: rc
      end subroutine state_reset_interface

      !> \brief Check if state module is allocated
      !!
      !! \param[in] this The state object to check
      !! \return True if allocated, false otherwise
      logical function state_is_allocated_interface(this)
         class(*), intent(in) :: this
      end function state_is_allocated_interface

      !> \brief Get memory usage of state module
      !!
      !! \param[in] this The state object
      !! \return Memory usage in bytes
      integer(8) function state_get_memory_usage_interface(this)
         class(*), intent(in) :: this
      end function state_get_memory_usage_interface

      !> \brief Print summary of state module
      !!
      !! \param[in] this The state object
      !! \param[in] unit Optional unit to write to (default stdout)
      subroutine state_print_summary_interface(this, unit)
         class(*), intent(in) :: this
         integer, optional, intent(in) :: unit
      end subroutine state_print_summary_interface
   end interface

   !> \brief Standard state module API naming conventions
   !!
   !! All state modules should use these procedure names:
   !! - init: Initialize the state module
   !! - validate: Validate the state module
   !! - cleanup: Clean up allocated memory
   !! - reset: Reset state values to defaults
   !! - is_allocated: Check allocation status
   !! - get_memory_usage: Get memory usage information
   !! - print_summary: Print diagnostic information
   !!
   !! Additional module-specific procedures are allowed but should follow
   !! the pattern: modulename_specific_operation (e.g., chemstate_add_species)

   public :: state_init_interface
   public :: state_validate_interface
   public :: state_cleanup_interface
   public :: state_reset_interface
   public :: state_is_allocated_interface
   public :: state_get_memory_usage_interface
   public :: state_print_summary_interface

end module StateInterface_Mod
