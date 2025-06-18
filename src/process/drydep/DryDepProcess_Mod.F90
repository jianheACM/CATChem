!> \file DryDepProcess_Mod.F90
!! \brief DryDep process implementation
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! Dry deposition process for atmospheric species removal through surface deposition
!!
module DryDepProcess_Mod
   use precision_mod
   use state_mod, only : StateContainerType
   use error_mod
   use ProcessInterface_Mod
   use GOCARTScheme_Mod
   use WeselyScheme_Mod
   use DryDepCommon_Mod

   implicit none
   private

   public :: DryDepProcessType

   !> DryDep process type extending ProcessInterface
   type, extends(ProcessInterface) :: DryDepProcessType
      private

      ! Process-specific configuration
      character(len=32) :: selected_scheme = 'GOCART'

   contains
      ! Required ProcessInterface methods
      procedure :: init => drydep_process_init
      procedure :: run => drydep_process_run
      procedure :: finalize => drydep_process_finalize
   end type DryDepProcessType

contains

   !> Initialize DryDep process
   subroutine drydep_process_init(this, container, rc)
      class(DryDepProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: message

      rc = CC_SUCCESS

      error_mgr => container%get_error_manager()

      ! Set process metadata
      this%name = 'DryDep'
      this%version = '1.0'
      this%description = 'Dry deposition process for atmospheric species removal through surface deposition'


      ! Validate inputs
      ! Basic validation
      if (.not. container%is_initialized) then
         call error_mgr%report_error("StateContainer not initialized")
         rc = -1
         return
      end if

      this%is_initialized = .true.
      this%is_active = .true.

      write(message, '(A,A,A)') 'DryDep process initialized with scheme: ', &
                                trim(this%selected_scheme)
      call error_mgr%report_info(message)

   end subroutine drydep_process_init

   !> Run DryDep process
   subroutine drydep_process_run(this, container, rc)
      class(DryDepProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS

      if (.not. this%is_ready()) then
         rc = -1
         return
      end if

      ! Get state pointers
      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()
      error_mgr => container%get_error_manager()

      ! Execute scheme-specific calculations
      select case (trim(this%selected_scheme))
      case ('GOCART')
         call this%run_gocart_scheme(container, rc)
      case ('Wesely')
         call this%run_wesely_scheme(container, rc)
      case default
         call error_mgr%report_error("Unknown scheme: " // trim(this%selected_scheme))
         rc = -1
         return
      end select

   end subroutine drydep_process_run

   !> Finalize DryDep process
   subroutine drydep_process_finalize(this, rc)
      class(DryDepProcessType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      this%is_initialized = .false.
      this%is_active = .false.

   end subroutine drydep_process_finalize

   !> Run GOCART scheme
   subroutine run_gocart_scheme(this, container, rc)
      class(DryDepProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state

      ! Local arrays for scheme interface
      real(fp), allocatable :: concentrations(:,:,:,:)
      real(fp), allocatable :: temperature(:,:,:)
      real(fp), allocatable :: pressure(:,:,:)
      real(fp), allocatable :: humidity(:,:,:)
      real(fp), allocatable :: loss_rate(:,:,:,:)

      rc = CC_SUCCESS

      ! Get state pointers
      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()

      ! TODO: Extract data from StateContainer to local arrays
      ! This is where you'd get the actual meteorological and chemical data
      ! from the StateContainer and prepare it for the scheme

      ! TODO: Allocate arrays based on grid dimensions
      ! call get_grid_dimensions(nx, ny, nz, n_species)
      ! allocate(temperature(nx, ny, nz), ...)

      ! Call the scheme calculation
      call gocart_calculate( &
         concentrations, temperature, pressure, humidity, loss_rate, &
         this%dt, rc)

      if (rc /= CC_SUCCESS) return

      ! TODO: Update StateContainer with results
      ! This is where you'd take the results from the scheme and update
      ! the appropriate state in the StateContainer

   end subroutine run_gocart_scheme

   !> Run Wesely scheme
   subroutine run_wesely_scheme(this, container, rc)
      class(DryDepProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state

      ! Local arrays for scheme interface
      real(fp), allocatable :: concentrations(:,:,:,:)
      real(fp), allocatable :: temperature(:,:,:)
      real(fp), allocatable :: pressure(:,:,:)
      real(fp), allocatable :: humidity(:,:,:)
      real(fp), allocatable :: loss_rate(:,:,:,:)

      rc = CC_SUCCESS

      ! Get state pointers
      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()

      ! TODO: Extract data from StateContainer to local arrays
      ! This is where you'd get the actual meteorological and chemical data
      ! from the StateContainer and prepare it for the scheme

      ! TODO: Allocate arrays based on grid dimensions
      ! call get_grid_dimensions(nx, ny, nz, n_species)
      ! allocate(temperature(nx, ny, nz), ...)

      ! Call the scheme calculation
      call wesely_calculate( &
         concentrations, temperature, pressure, humidity, loss_rate, &
         this%dt, rc)

      if (rc /= CC_SUCCESS) return

      ! TODO: Update StateContainer with results
      ! This is where you'd take the results from the scheme and update
      ! the appropriate state in the StateContainer

   end subroutine run_wesely_scheme

end module DryDepProcess_Mod