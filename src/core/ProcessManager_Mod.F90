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
   use StateManager_Mod, only : StateManagerType
   use error_mod, only : CC_SUCCESS, CC_FAILURE
   use ProcessInterface_Mod, only : ProcessInterface, ColumnProcessInterface
   use ProcessFactory_Mod, only : ProcessFactoryType
   use GridManager_Mod, only : GridManagerType, ColumnIteratorType
   use ColumnInterface_Mod, only : VirtualColumnType, ColumnProcessorType

   implicit none
   private

   public :: ProcessManagerType

   type :: ProcessManagerType
      private
      class(ProcessInterface), allocatable :: processes(:)
      integer :: num_processes = 0
      integer :: max_processes = 50
      type(ProcessFactoryType) :: factory
      type(ColumnProcessorType) :: column_processor  !< Batch column processor
   contains
      procedure :: init => manager_init
      procedure :: add_process => manager_add_process
      procedure :: run_all => manager_run_all
      procedure :: run_process => manager_run_process
      procedure :: run_column_processes => manager_run_column_processes
      procedure :: run_process_on_columns => manager_run_process_on_columns
      procedure :: finalize => manager_finalize
      procedure :: list_processes => manager_list_processes
      procedure :: get_column_processes => manager_get_column_processes
      procedure :: configure_run_phases => manager_configure_run_phases
      procedure :: run_phase => manager_run_phase
      procedure :: run_all_processes => manager_run_all_processes
      procedure :: set_max_processes => manager_set_max_processes
      procedure :: enable_column_batching => manager_enable_column_batching
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

      ! Initialize column processor
      call this%column_processor%init(rc)
      if (rc /= CC_SUCCESS) return

      rc = CC_SUCCESS
   end subroutine manager_init

   !> \brief Add a process to the manager
   subroutine manager_add_process(this, process_name, scheme_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: scheme_name
      type(StateManagerType), intent(inout) :: container
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

   !> \brief Run all processes using appropriate method (3D or column)
   subroutine manager_run_all(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, local_rc
      type(GridManagerType), pointer :: grid_mgr

      rc = CC_SUCCESS

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = CC_FAILURE
         return
      endif

      do i = 1, this%num_processes
         if (this%processes(i)%is_ready()) then
            ! Check if this is a column process
            select type(proc => this%processes(i))
            class is (ColumnProcessInterface)
               ! Run column-based process
               call this%run_process_on_columns(i, container, rc)
            class default
               ! Run traditional 3D process
               call this%processes(i)%run(container, local_rc)
            end select

            if (local_rc /= CC_SUCCESS) then
               rc = local_rc
               return
            endif
         endif
      enddo
   end subroutine manager_run_all

   !> \brief Run column processes using column virtualization
   subroutine manager_run_column_processes(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, local_rc
      type(GridManagerType), pointer :: grid_mgr
      type(ColumnIteratorType) :: col_iter
      type(VirtualColumnType) :: virtual_col
      integer :: col_idx

      rc = CC_SUCCESS

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = CC_FAILURE
         return
      endif

      ! Initialize column iterator
      call grid_mgr%get_column_iterator(col_iter, rc)
      if (rc /= CC_SUCCESS) return

      ! Process each column
      do while (col_iter%has_next())
         call col_iter%next(col_idx, rc)
         if (rc /= CC_SUCCESS) return

         ! Create virtual column for this index
         call grid_mgr%create_virtual_column(col_idx, virtual_col, rc)
         if (rc /= CC_SUCCESS) return

         ! Extract column data from container
         call virtual_col%extract_from_container(container, rc)
         if (rc /= CC_SUCCESS) return

         ! Run all column processes on this column
         do i = 1, this%num_processes
            select type(proc => this%processes(i))
            class is (ColumnProcessInterface)
               if (proc%is_ready()) then
                  call proc%run_column(virtual_col, local_rc)
                  if (local_rc /= CC_SUCCESS) then
                     rc = local_rc
                     return
                  endif
               endif
            end select
         enddo

         ! Update container with column data
         call virtual_col%update_container(container, rc)
         if (rc /= CC_SUCCESS) return
      enddo
   end subroutine manager_run_column_processes

   !> \brief Run a specific process on all columns
   subroutine manager_run_process_on_columns(this, process_index, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(in) :: process_index
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(GridManagerType), pointer :: grid_mgr
      type(ColumnIteratorType) :: col_iter
      type(VirtualColumnType) :: virtual_col
      integer :: col_idx, local_rc

      rc = CC_SUCCESS

      if (process_index < 1 .or. process_index > this%num_processes) then
         rc = CC_FAILURE
         return
      endif

      ! Get grid manager from container
      grid_mgr => container%get_grid_manager()
      if (.not. associated(grid_mgr)) then
         rc = CC_FAILURE
         return
      endif

      ! Only proceed if this is a column process
      select type(proc => this%processes(process_index))
      class is (ColumnProcessInterface)
         ! Initialize column iterator
         call grid_mgr%get_column_iterator(col_iter, rc)
         if (rc /= CC_SUCCESS) return

         ! Process each column
         do while (col_iter%has_next())
            call col_iter%next(col_idx, rc)
            if (rc /= CC_SUCCESS) return

            ! Create virtual column for this index
            call grid_mgr%create_virtual_column(col_idx, virtual_col, rc)
            if (rc /= CC_SUCCESS) return

            ! Extract column data from container
            call virtual_col%extract_from_container(container, rc)
            if (rc /= CC_SUCCESS) return

            ! Run the process on this column
            call proc%run_column(virtual_col, local_rc)
            if (local_rc /= CC_SUCCESS) then
               rc = local_rc
               return
            endif

            ! Update container with column data
            call virtual_col%update_container(container, rc)
            if (rc /= CC_SUCCESS) return
         enddo
      class default
         ! For non-column processes, run normally
         call this%processes(process_index)%run(container, rc)
      end select
   end subroutine manager_run_process_on_columns

   !> \brief Run a specific process by name
   subroutine manager_run_process(this, process_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i

      rc = CC_FAILURE
      do i = 1, this%num_processes
         if (trim(this%processes(i)%get_name()) == trim(process_name)) then
            ! Check if this is a column process and run appropriately
            select type(proc => this%processes(i))
            class is (ColumnProcessInterface)
               call this%run_process_on_columns(i, container, rc)
            class default
               call this%processes(i)%run(container, rc)
            end select
            return
         endif
      enddo
   end subroutine manager_run_process

   !> \brief Get list of column processes
   subroutine manager_get_column_processes(this, column_indices, count)
      class(ProcessManagerType), intent(in) :: this
      integer, intent(out) :: column_indices(:)
      integer, intent(out) :: count

      integer :: i, max_count

      count = 0
      max_count = min(this%num_processes, size(column_indices))

      do i = 1, this%num_processes
         select type(proc => this%processes(i))
         class is (ColumnProcessInterface)
            count = count + 1
            if (count <= max_count) then
               column_indices(count) = i
            endif
         end select
      enddo
   end subroutine manager_get_column_processes

   !> Configure run phases for multi-phase execution
   subroutine manager_configure_run_phases(this, phase_names, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=64), intent(in) :: phase_names(:)
      integer, intent(out) :: rc

      ! TODO: Implement run phase configuration
      ! This would involve setting up process execution order by phase
      rc = CC_SUCCESS

   end subroutine manager_configure_run_phases

   !> Run a specific run phase
   subroutine manager_run_phase(this, phase_name, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      character(len=*), intent(in) :: phase_name
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = CC_SUCCESS

      ! Run all processes configured for this phase
      do i = 1, this%num_processes
         ! TODO: Check if process belongs to this phase
         ! For now, run all processes
         if (this%processes(i)%is_ready()) then
            select type(proc => this%processes(i))
            class is (ColumnProcessInterface)
               call this%run_process_on_columns(i, container, rc)
            class default
               call this%processes(i)%run(container, local_rc)
            end select

            if (local_rc /= CC_SUCCESS) then
               rc = local_rc
               return
            endif
         endif
      enddo

   end subroutine manager_run_phase

   !> Run all processes (compatibility method)
   subroutine manager_run_all_processes(this, container, rc)
      class(ProcessManagerType), intent(inout) :: this
      type(StateManagerType), intent(inout) :: container
      integer, intent(out) :: rc

      call this%run_all(container, rc)

   end subroutine manager_run_all_processes

   !> Set maximum number of processes
   subroutine manager_set_max_processes(this, max_processes, rc)
      class(ProcessManagerType), intent(inout) :: this
      integer, intent(in) :: max_processes
      integer, intent(out) :: rc

      ! TODO: Reallocate process array if needed
      this%max_processes = max_processes
      rc = CC_SUCCESS

   end subroutine manager_set_max_processes

   !> Enable or disable column batching
   subroutine manager_enable_column_batching(this, enable, rc)
      class(ProcessManagerType), intent(inout) :: this
      logical, intent(in) :: enable
      integer, intent(out) :: rc

      ! TODO: Configure column processor batching mode
      rc = CC_SUCCESS

   end subroutine manager_enable_column_batching

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

      ! Finalize column processor
      call this%column_processor%finalize(local_rc)
      if (local_rc /= CC_SUCCESS) then
         rc = local_rc
      endif

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
