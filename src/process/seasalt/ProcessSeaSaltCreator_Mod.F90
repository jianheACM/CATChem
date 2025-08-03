!> \file ProcessSeaSaltCreator_Mod.F90
!! \brief Creator module for SeaSalt processes
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! This module provides the creator function for SeaSalt processes, following the
!! ProcessCreatorInterface pattern defined in ProcessRegistry_Mod.F90
!!
!! The creator function is responsible for:
!! - Creating (allocating) process instances
!! - Basic initialization (no container required)
!! - Returning polymorphic ProcessInterface instances
!!
!! Note: Full initialization with StateContainer happens later via the init() method.
!!
module ProcessSeaSaltCreator_Mod
   use precision_mod
   use error_mod, only : CC_SUCCESS, CC_FAILURE
   use ProcessInterface_Mod, only : ProcessInterface
   use ProcessSeaSaltInterface_Mod, only : ProcessSeaSaltInterfaceType

   implicit none
   private

   public :: create_seasalt_process

contains

   !> \brief Create a SeaSalt process instance
   !!
   !! This function creates a new SeaSalt process instance following the
   !! ProcessCreatorInterface pattern. The process is allocated but not initialized
   !! with a container - that happens later via the init() method.
   !!
   !! \param[out] process Allocated process instance (polymorphic ProcessInterface)
   !! \param[out] rc Return code (CC_SUCCESS or CC_FAILURE)
   subroutine create_seasalt_process(process, rc)
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ProcessSeaSaltInterfaceType), allocatable :: seasalt_process

      ! Initialize return code
      rc = CC_SUCCESS

      ! Allocate the SeaSalt process interface
      allocate(seasalt_process, stat=rc)
      if (rc /= 0) then
         rc = CC_FAILURE
         return
      endif

      ! Move to the polymorphic output
      ! Note: The process is created but not initialized here.
      ! Initialization with container happens later via the init() method.
      call move_alloc(seasalt_process, process)

   end subroutine create_seasalt_process

end module ProcessSeaSaltCreator_Mod