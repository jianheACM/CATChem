!> \file ProcessFactory_Mod.F90
!! \brief Process factory for dynamic process creation following architecture guide
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module implements the ProcessFactory as specified in PROCESS_ARCHITECTURE_GUIDE.md
!!
module ProcessFactory_Mod
   use precision_mod
   use state_mod, only : StateContainerType
   use error_mod
   use ProcessInterface_Mod
   use ProcessRegistry_Mod
   ! Import process creators
   use DustProcessCreator_Mod, only : create_dust_process

   implicit none
   private

   public :: ProcessFactoryType
   public :: create_process

   type :: ProcessFactoryType
      private
      type(ProcessRegistryType) :: registry
   contains
      procedure :: init => factory_init
      procedure :: create_process => factory_create_process
      procedure :: list_available => factory_list_available
      procedure :: register_process => factory_register_process
   end type ProcessFactoryType

contains

   subroutine factory_init(this, rc)
      class(ProcessFactoryType), intent(inout) :: this
      integer, intent(out) :: rc

      call this%registry%init(rc)
      if (rc /= CC_SUCCESS) return

      ! Register built-in processes
      call this%register_builtin_processes(rc)
   end subroutine factory_init

   subroutine factory_create_process(this, process_name, scheme_name, container, process, rc)
      class(ProcessFactoryType), intent(inout) :: this
      character(len=*), intent(in) :: process_name
      character(len=*), intent(in) :: scheme_name
      type(StateContainerType), intent(inout) :: container
      class(ProcessInterface), allocatable, intent(out) :: process
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=128) :: full_process_name

      error_mgr => container%get_error_manager()
      call error_mgr%push_context('factory_create_process', &
                                  'creating process: ' // trim(process_name))

      ! Create full process name including scheme
      if (len_trim(scheme_name) > 0) then
         full_process_name = trim(process_name) // '_' // trim(scheme_name)
      else
         full_process_name = trim(process_name)
      endif

      ! Create process from registry
      call this%registry%create_process(full_process_name, process, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_INVALID_CONFIG, &
                                     'Unknown process: ' // trim(full_process_name), rc, &
                                     'factory_create_process', &
                                     'Check available processes with list_available()')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
   end subroutine factory_create_process

   subroutine factory_list_available(this, process_names, rc)
      class(ProcessFactoryType), intent(in) :: this
      character(len=64), allocatable, intent(out) :: process_names(:)
      integer, intent(out) :: rc

      call this%registry%list_processes(process_names, rc)
   end subroutine factory_list_available

   subroutine factory_register_process(this, name, category, description, creator, rc)
      class(ProcessFactoryType), intent(inout) :: this
      character(len=*), intent(in) :: name, category, description
      procedure(ProcessCreatorInterface) :: creator
      integer, intent(out) :: rc

      call this%registry%register_process(name, category, description, creator, rc)
   end subroutine factory_register_process

   !> \brief Register all built-in processes
   subroutine register_builtin_processes(this, rc)
      class(ProcessFactoryType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Register dust process
      call this%register_process('dust', 'emission', &
                                  'Dust emission process with multiple schemes', &
                                  create_dust_process, rc)
      if (rc /= CC_SUCCESS) return

      ! TODO: Register other processes as they are migrated
      ! call this%register_process('seasalt', 'emission', &
      !                           'Sea salt emission process', &
      !                           create_seasalt_process, rc)
      ! if (rc /= CC_SUCCESS) return

   end subroutine register_builtin_processes

   !> \brief Module-level convenience function
   function create_process(process_name, scheme_name, container, rc) result(process)
      character(len=*), intent(in) :: process_name, scheme_name
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc
      class(ProcessInterface), allocatable :: process

      type(ProcessFactoryType) :: factory

      call factory%init(rc)
      if (rc /= CC_SUCCESS) return

      call factory%create_process(process_name, scheme_name, container, process, rc)
   end function create_process

end module ProcessFactory_Mod
