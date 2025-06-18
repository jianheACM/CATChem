!> \file ProcessManager_Mod.F90
!! \brief High-level process management following the architecture guide
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides high-level management of multiple processes,
!! following the exact structure specified in PROCESS_ARCHITECTURE_GUIDE.md
!!
module ProcessManager_Mod
   use precision_mod
   use state_mod, only : StateContainerType
   use error_mod, only : CC_SUCCESS, CC_FAILURE
   use ProcessInterface_Mod, only : ProcessInterface
   use ProcessFactory_Mod, only : ProcessFactoryType

   implicit none
   private

   public :: ProcessManagerType

   type :: ProcessManagerType
      private
      class(ProcessInterface), allocatable :: processes(:)
      integer :: num_processes = 0
      integer :: max_processes = 50
      type(ProcessFactoryType) :: factory
   contains
      procedure :: init => manager_init
      procedure :: add_process => manager_add_process
      procedure :: run_all => manager_run_all
      procedure :: run_process => manager_run_process
      procedure :: finalize => manager_finalize
      procedure :: list_processes => manager_list_processes
   end type ProcessManagerType

contains

   !> \brief Initialize the process manager
   subroutine manager_init(this, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%factory%init(rc)
      if (rc /= CC_SUCCESS) return

      allocate(this%processes(this%max_processes))
      this%num_processes = 0

      rc = CC_SUCCESS
   end subroutine manager_init

   !> \brief Add a process to the manager
   subroutine manager_add_process(this, process_name, scheme_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: scheme_name
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      class(ProcessInterface), allocatable :: new_process

      if (this%num_processes >= this%max_processes) then
         rc = CC_FAILURE
         return
      endif

      ! Create the process
      call this%factory%create_process(process_name, scheme_name, container, new_process, rc)
      if (rc /= CC_SUCCESS) return

      ! Initialize the process
      call new_process%init(container, rc)
      if (rc /= CC_SUCCESS) return

      ! Add to manager
      this%num_processes = this%num_processes + 1
      this%processes(this%num_processes) = new_process

      rc = CC_SUCCESS
   end subroutine manager_add_process

   !> \brief Run all processes
   subroutine manager_run_all(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = CC_SUCCESS
      do i = 1, this%num_processes
         if (this%processes(i)%is_ready()) then
            call this%processes(i)%run(container, local_rc)
            if (local_rc /= CC_SUCCESS) then
               rc = local_rc
               return
            endif
         endif
      enddo
   end subroutine manager_run_all

   !> \brief Run a specific process by name
   subroutine manager_run_process(this, process_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i

      rc = CC_FAILURE
      do i = 1, this%num_processes
         if (trim(this%processes(i)%get_name()) == trim(process_name)) then
            call this%processes(i)%run(container, rc)
            return
         endif
      enddo
   end subroutine manager_run_process

   !> \brief Finalize all processes
   subroutine manager_finalize(this, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = CC_SUCCESS
      do i = 1, this%num_processes
         call this%processes(i)%finalize(local_rc)
         if (local_rc /= CC_SUCCESS) then
            rc = local_rc
         endif
      enddo

      if (allocated(this%processes)) deallocate(this%processes)
      this%num_processes = 0
   end subroutine manager_finalize

   !> \brief List all processes
   subroutine manager_list_processes(this, process_names, count)
      class(ProcessManagerType), intent(in) :: this
      character(len=64), intent(out) :: process_names(:)
      integer, intent(out) :: count

      integer :: i, max_count

      max_count = min(this%num_processes, size(process_names))
      do i = 1, max_count
         process_names(i) = this%processes(i)%get_name()
      enddo
      count = max_count
   end subroutine manager_list_processes

end module ProcessManager_Mod
