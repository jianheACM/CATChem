!> \file ProcessSeaSaltInterface_Mod.F90
!! \brief SeaSalt process interface module
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! Sea salt aerosol emission from ocean surfaces via bubble bursting and wave breaking
!!
!! This module provides a clean interface for SeaSalt processes,
!! acting as a coordinator between the framework and scheme implementations.
!! All physics/chemistry calculations are delegated to scheme modules.
!!
module ProcessSeaSaltInterface_Mod
   use precision_mod
   use state_mod, only : StateContainerType
   use error_mod
   use ProcessInterface_Mod
   use DiagnosticInterface_Mod
   use ChemState_Mod
   use Monahan_Mod
   use Smith_Mod
   use Jaegle_Mod
   use Ovadnevaite_Mod

   implicit none
   private

   public :: ProcessSeaSaltInterfaceType

   !> SeaSalt process interface type extending ColumnProcessInterface
   !! Supports both column and full 3D processing modes
   type, extends(ColumnProcessInterface) :: ProcessSeaSaltInterfaceType
      private

      ! Configuration settings
      character(len=32) :: selected_scheme = 'Monahan'
      logical :: scheme_configured = .false.

      ! Diagnostic management
      type(DiagnosticInterface), pointer :: diag_interface => null()
      logical :: diagnostics_registered = .false.

   contains
      ! Required ProcessInterface methods
      procedure :: init => seasalt_interface_init
      procedure :: run => seasalt_interface_run
      procedure :: finalize => seasalt_interface_finalize

      ! Required ColumnProcessInterface methods
      procedure :: init_column_processing => seasalt_init_column_processing
      procedure :: run_column => seasalt_run_column
      procedure :: finalize_column_processing => seasalt_finalize_column_processing

      ! Process-specific configuration methods
      procedure :: set_scheme => seasalt_set_scheme
      procedure :: get_scheme => seasalt_get_scheme
      procedure :: configure_scheme => seasalt_configure_scheme
   end type ProcessSeaSaltInterfaceType

contains

   !> Initialize SeaSalt process interface
   subroutine seasalt_interface_init(this, container, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message

      rc = CC_SUCCESS

      error_mgr => container%get_error_manager()

      ! Set process metadata
      this%name = 'SeaSalt'
      this%version = '1.0'
      this%description = 'Sea salt aerosol emission from ocean surfaces via bubble bursting and wave breaking'

      ! Basic validation
      if (.not. container%is_initialized) then
         call error_mgr%report_error("StateContainer not initialized")
         rc = -1
         return
      end if

      ! Configure the selected scheme
      call this%configure_scheme(container, rc)
      if (rc /= CC_SUCCESS) return

      this%is_initialized = .true.
      this%is_active = .true.

      write(message, '(A,A,A)') 'SeaSalt process interface initialized with scheme: ', &
                                trim(this%selected_scheme)
      call error_mgr%report_info(message)

   end subroutine seasalt_interface_init

   !> Run SeaSalt process interface
   subroutine seasalt_interface_run(this, container, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS

      if (.not. this%is_ready()) then
         rc = -1
         return
      end if

      if (.not. this%scheme_configured) then
         error_mgr => container%get_error_manager()
         call error_mgr%report_error("Scheme not configured")
         rc = -1
         return
      end if

      ! Delegate to the configured scheme
      select case (trim(this%selected_scheme))
      case ('Monahan')
         call monahan_scheme_run(container, rc)
      case ('Smith')
         call smith_scheme_run(container, rc)
      case ('Jaegle')
         call jaegle_scheme_run(container, rc)
      case ('Ovadnevaite')
         call ovadnevaite_scheme_run(container, rc)
      case default
         error_mgr => container%get_error_manager()
         call error_mgr%report_error("Unknown scheme: " // trim(this%selected_scheme))
         rc = -1
         return
      end select

   end subroutine seasalt_interface_run

   !> Finalize SeaSalt process interface
   subroutine seasalt_interface_finalize(this, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      this%is_initialized = .false.
      this%is_active = .false.
      this%scheme_configured = .false.

   end subroutine seasalt_interface_finalize

   !> Set the scheme to use
   subroutine seasalt_set_scheme(this, scheme_name, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      character(len=*), intent(in) :: scheme_name
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Validate scheme name
      select case (trim(scheme_name))
      case ('Monahan')
         this%selected_scheme = trim(scheme_name)
      case ('Smith')
         this%selected_scheme = trim(scheme_name)
      case ('Jaegle')
         this%selected_scheme = trim(scheme_name)
      case ('Ovadnevaite')
         this%selected_scheme = trim(scheme_name)
      case default
         rc = -1
         return
      end select

      this%scheme_configured = .false.  ! Force reconfiguration

   end subroutine seasalt_set_scheme

   !> Get the current scheme name
   function seasalt_get_scheme(this) result(scheme_name)
      class(ProcessSeaSaltInterfaceType), intent(in) :: this
      character(len=32) :: scheme_name

      scheme_name = this%selected_scheme

   end function seasalt_get_scheme

   !> Configure the selected scheme
   subroutine seasalt_configure_scheme(this, container, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Configure the selected scheme
      select case (trim(this%selected_scheme))
      case ('Monahan')
         call monahan_configure(container, rc)
      case ('Smith')
         call smith_configure(container, rc)
      case ('Jaegle')
         call jaegle_configure(container, rc)
      case ('Ovadnevaite')
         call ovadnevaite_configure(container, rc)
      case default
         rc = -1
         return
      end select

      if (rc == CC_SUCCESS) then
         this%scheme_configured = .true.
      end if

   end subroutine seasalt_configure_scheme

   !========================================================================
   ! Column Processing Methods
   !========================================================================

   !> Initialize column processing for SeaSalt
   subroutine seasalt_init_column_processing(this, container, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Column processing initialization
      rc = CC_SUCCESS

      ! Set up any column-specific resources
      ! Default implementation - processes can override if needed

   end subroutine seasalt_init_column_processing

   !> Process a single column for SeaSalt
   subroutine seasalt_run_column(this, column, rc)
      use ColumnInterface_Mod, only : VirtualColumnType
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      type(VirtualColumnType), intent(inout) :: column
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Extract column data and process using the selected scheme
      ! This would call the same scheme as the full 3D version
      ! but operates on column data extracted from VirtualColumnType

      ! TODO: Implement column-specific processing
      ! Example:
      ! - Extract met fields for this column
      ! - Extract species concentrations for this column
      ! - Call the appropriate scheme
      ! - Update column with new values

   end subroutine seasalt_run_column

   !> Finalize column processing for SeaSalt
   subroutine seasalt_finalize_column_processing(this, container, rc)
      class(ProcessSeaSaltInterfaceType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! Column processing cleanup
      rc = CC_SUCCESS

      ! Clean up any column-specific resources
      ! Default implementation - processes can override if needed

   end subroutine seasalt_finalize_column_processing

end module ProcessSeaSaltInterface_Mod