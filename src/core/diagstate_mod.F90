!> \file diagstate_mod.F90
!! \brief Modern diagnostic state management for CATChem
!! \ingroup core_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides diagnostic state management for atmospheric chemistry
!! modeling, including aerosol optical properties, emission fluxes, and
!! deposition diagnostics.
!!
!! \details
!! The diagnostic state tracks key atmospheric chemistry diagnostic variables:
!! - Aerosol optical depth (AOD) at multiple wavelengths
!! - Dust and sea salt emission fluxes
!! - Dry deposition velocities and frequencies
!! - Plume rise heights from different algorithms
!!
!! \section diagstate_usage Usage Example
!! \code{.f90}
!! use diagstate_mod
!! type(DiagStateType) :: diag_state
!! integer :: rc
!!
!! call diag_state%init(nlevs, nspecies_drydep, rc)
!! call diag_state%validate(rc)
!! call diag_state%cleanup(rc)
!! \endcode
!!
module DiagState_Mod
   use precision_mod, only: fp
   use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, &
                        ERROR_MEMORY_ALLOCATION, ERROR_INVALID_INPUT, &
                        ERROR_BOUNDS_CHECK, ERROR_NONE

   implicit none
   private

   public :: DiagStateType

   !> \brief Diagnostic state type for atmospheric chemistry
   !!
   !! This type manages diagnostic variables for atmospheric chemistry modeling
   !! with proper error handling and modern Fortran patterns.
   type :: DiagStateType
      private

      ! Initialization and error management
      logical :: is_initialized = .false.              !< Initialization status
      type(ErrorManagerType) :: error_mgr              !< Integrated error manager

      ! Grid dimensions
      integer :: nlevs = 0                             !< Number of vertical levels
      integer :: nspecies_drydep = 0                   !< Number of dry deposition species

      ! Scalar diagnostic variables
      real(fp) :: dust_total_flux = 0.0_fp             !< Total dust flux [kg m-2 s-1]
      real(fp) :: seasalt_total_flux = 0.0_fp          !< Total sea salt flux [kg m-2 s-1]
      real(fp) :: briggs_plumerise_height = 0.0_fp     !< Briggs plume rise height [m]
      real(fp) :: sofiev_plumerise_height = 0.0_fp     !< Sofiev plume rise height [m]

      ! Arrays for level-dependent diagnostics
      real(fp), allocatable :: AOD550(:)               !< AOD at 550nm [1]
      real(fp), allocatable :: AOD380(:)               !< AOD at 380nm [1]
      real(fp), allocatable :: TOMSAI(:)               !< TOMS Aerosol Index [1]

      ! Dry deposition diagnostics
      real(fp), allocatable :: drydep_velocity(:)      !< Dry deposition velocity [m/s]
      real(fp), allocatable :: drydep_frequency(:)     !< Dry deposition frequency [s-1]

      ! Validation status
      logical :: is_valid = .false.                    !< Validation status

   contains
      ! Enhanced methods with error handling
      procedure :: init => diagstate_init
      procedure :: validate => diagstate_validate
      procedure :: cleanup => diagstate_cleanup
      procedure :: reset => diagstate_reset
      procedure :: is_allocated => diagstate_is_allocated
      procedure :: get_memory_usage => diagstate_get_memory_usage
      procedure :: print_summary => diagstate_print_summary
      procedure :: get_dust_flux => diagstate_get_dust_flux
      procedure :: set_dust_flux => diagstate_set_dust_flux
      procedure :: get_aod550 => diagstate_get_aod550
      procedure :: set_aod550 => diagstate_set_aod550
      procedure :: print_info => diagstate_print_info

   end type DiagStateType

contains

   !> \brief Initialize diagnostic state
   !!
   !! This subroutine initializes the diagnostic state with specified dimensions
   !! and allocates necessary arrays.
   !!
   !! \param[inout] this The diagnostic state object
   !! \param[in] nlevs Number of vertical levels
   !! \param[in] nspecies_drydep Number of species with dry deposition
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   subroutine diagstate_init(this, nlevs, nspecies_drydep, error_mgr, rc)
      use error_mod, only: ErrorManagerType
      implicit none
      class(DiagStateType), intent(inout) :: this
      integer, intent(in) :: nlevs
      integer, intent(in) :: nspecies_drydep
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: alloc_stat

      rc = CC_SUCCESS

      ! Use external error manager
      call error_mgr%push_context('diagstate_init', 'diagstate_mod.F90')

      ! Validate inputs
      if (nlevs <= 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Number of levels must be positive', rc, 'diagstate_init')
         call error_mgr%pop_context()
         return
      endif

      if (nspecies_drydep < 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Number of dry deposition species cannot be negative', rc, 'diagstate_init')
         call error_mgr%pop_context()
         return
      endif

      ! Store dimensions
      this%nlevs = nlevs
      this%nspecies_drydep = nspecies_drydep

      ! Allocate arrays for level-dependent diagnostics
      allocate(this%AOD550(nlevs), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
              'Failed to allocate AOD550 array', rc, 'diagstate_init')
         call error_mgr%pop_context()
         return
      endif
      this%AOD550 = 0.0_fp

      allocate(this%AOD380(nlevs), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
              'Failed to allocate AOD380 array', rc, 'diagstate_init')
         call error_mgr%pop_context()
         return
      endif
      this%AOD380 = 0.0_fp

      allocate(this%TOMSAI(nlevs), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
              'Failed to allocate TOMSAI array', rc, 'diagstate_init')
         call error_mgr%pop_context()
         return
      endif
      this%TOMSAI = 0.0_fp

      ! Allocate dry deposition arrays if needed
      if (nspecies_drydep > 0) then
         allocate(this%drydep_velocity(nspecies_drydep), stat=alloc_stat)
         if (alloc_stat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                 'Failed to allocate drydep_velocity array', rc, 'diagstate_init')
            call error_mgr%pop_context()
            return
         endif
         this%drydep_velocity = 0.0_fp

         allocate(this%drydep_frequency(nspecies_drydep), stat=alloc_stat)
         if (alloc_stat /= 0) then
            call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                 'Failed to allocate drydep_frequency array', rc, 'diagstate_init')
            call error_mgr%pop_context()
            return
         endif
         this%drydep_frequency = 0.0_fp
      endif

      ! Initialize scalar variables to zero
      this%dust_total_flux = 0.0_fp
      this%seasalt_total_flux = 0.0_fp
      this%briggs_plumerise_height = 0.0_fp
      this%sofiev_plumerise_height = 0.0_fp

      ! Mark as initialized
      this%is_initialized = .true.
      this%is_valid = .true.

      call error_mgr%pop_context()

   end subroutine diagstate_init

   !> \brief Validate diagnostic state
   !!
   !! This function validates that the diagnostic state is properly initialized
   !! and all arrays have consistent dimensions.
   !!
   !! \param[in] this The diagnostic state object
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   subroutine diagstate_validate(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType
      implicit none
      class(DiagStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      call error_mgr%push_context('diagstate_validate', 'diagstate_mod.F90')
      rc = CC_SUCCESS

      ! Check initialization
      if (.not. this%is_initialized) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Diagnostic state not initialized', rc, 'diagstate_validate')
         call error_mgr%pop_context()
         return
      endif

      ! Check array allocations and dimensions
      if (.not. allocated(this%AOD550) .or. size(this%AOD550) /= this%nlevs) then
         call error_mgr%report_error(ERROR_BOUNDS_CHECK, &
              'AOD550 array dimension mismatch', rc, 'diagstate_validate')
         call error_mgr%pop_context()
         return
      endif

      if (.not. allocated(this%AOD380) .or. size(this%AOD380) /= this%nlevs) then
         call error_mgr%report_error(ERROR_BOUNDS_CHECK, &
              'AOD380 array dimension mismatch', rc, 'diagstate_validate')
         call error_mgr%pop_context()
         return
      endif

      if (.not. allocated(this%TOMSAI) .or. size(this%TOMSAI) /= this%nlevs) then
         call error_mgr%report_error(ERROR_BOUNDS_CHECK, &
              'TOMSAI array dimension mismatch', rc, 'diagstate_validate')
         call error_mgr%pop_context()
         return
      endif

      ! Check dry deposition arrays if they should exist
      if (this%nspecies_drydep > 0) then
         if (.not. allocated(this%drydep_velocity) .or. size(this%drydep_velocity) /= this%nspecies_drydep) then
            call error_mgr%report_error(ERROR_BOUNDS_CHECK, &
                 'Dry deposition velocity array dimension mismatch', rc, 'diagstate_validate')
            call error_mgr%pop_context()
            return
         endif

         if (.not. allocated(this%drydep_frequency) .or. size(this%drydep_frequency) /= this%nspecies_drydep) then
            call error_mgr%report_error(ERROR_BOUNDS_CHECK, &
                 'Dry deposition frequency array dimension mismatch', rc, 'diagstate_validate')
            call error_mgr%pop_context()
            return
         endif
      endif

      call error_mgr%pop_context()

   end subroutine diagstate_validate

   !> \brief Clean up diagnostic state
   !!
   !! This subroutine deallocates all arrays and resets the diagnostic state.
   !!
   !! \param[inout] this The diagnostic state object
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   subroutine diagstate_cleanup(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType
      implicit none
      class(DiagStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      call error_mgr%push_context('diagstate_cleanup', 'diagstate_mod.F90')

      ! Deallocate arrays
      if (allocated(this%AOD550)) deallocate(this%AOD550)
      if (allocated(this%AOD380)) deallocate(this%AOD380)
      if (allocated(this%TOMSAI)) deallocate(this%TOMSAI)
      if (allocated(this%drydep_velocity)) deallocate(this%drydep_velocity)
      if (allocated(this%drydep_frequency)) deallocate(this%drydep_frequency)

      ! Reset dimensions and status
      this%nlevs = 0
      this%nspecies_drydep = 0
      this%is_initialized = .false.
      this%is_valid = .false.

      ! Reset scalar variables
      this%dust_total_flux = 0.0_fp
      this%seasalt_total_flux = 0.0_fp
      this%briggs_plumerise_height = 0.0_fp
      this%sofiev_plumerise_height = 0.0_fp

      call error_mgr%pop_context()

   end subroutine diagstate_cleanup

   !> \brief Get dust emission flux
   !!
   !! \param[in] this The diagnostic state object
   !! \return Dust emission flux [kg m-2 s-1]
   function diagstate_get_dust_flux(this) result(flux)
      implicit none
      class(DiagStateType), intent(in) :: this
      real(fp) :: flux

      flux = this%dust_total_flux

   end function diagstate_get_dust_flux

   !> \brief Set dust emission flux
   !!
   !! \param[inout] this The diagnostic state object
   !! \param[in] flux Dust emission flux [kg m-2 s-1]
   !! \param[out] rc Return code
   subroutine diagstate_set_dust_flux(this, flux, rc)
      implicit none
      class(DiagStateType), intent(inout) :: this
      real(fp), intent(in) :: flux
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (flux < 0.0_fp) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      this%dust_total_flux = flux

   end subroutine diagstate_set_dust_flux

   !> \brief Get AOD at 550nm for a specific level
   !!
   !! \param[in] this The diagnostic state object
   !! \param[in] level Vertical level index
   !! \param[out] rc Return code
   !! \return AOD at 550nm [dimensionless]
   function diagstate_get_aod550(this, level, rc) result(aod)
      implicit none
      class(DiagStateType), intent(in) :: this
      integer, intent(in) :: level
      integer, intent(out) :: rc
      real(fp) :: aod

      rc = CC_SUCCESS
      aod = 0.0_fp

      if (.not. allocated(this%AOD550)) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      if (level < 1 .or. level > size(this%AOD550)) then
         rc = ERROR_BOUNDS_CHECK
         return
      endif

      aod = this%AOD550(level)

   end function diagstate_get_aod550

   !> \brief Set AOD at 550nm for a specific level
   !!
   !! \param[inout] this The diagnostic state object
   !! \param[in] level Vertical level index
   !! \param[in] aod AOD value [dimensionless]
   !! \param[out] rc Return code
   subroutine diagstate_set_aod550(this, level, aod, rc)
      implicit none
      class(DiagStateType), intent(inout) :: this
      integer, intent(in) :: level
      real(fp), intent(in) :: aod
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (.not. allocated(this%AOD550)) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      if (level < 1 .or. level > size(this%AOD550)) then
         rc = ERROR_BOUNDS_CHECK
         return
      endif

      if (aod < 0.0_fp) then
         rc = ERROR_INVALID_INPUT
         return
      endif

      this%AOD550(level) = aod

   end subroutine diagstate_set_aod550

   !> \brief Reset diagnostic state values
   !!
   !! This subroutine resets diagnostic state values to their defaults
   !! without deallocating arrays.
   !!
   !! \param[inout] this The diagnostic state object
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   subroutine diagstate_reset(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS
      implicit none
      class(DiagStateType), intent(inout) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      call error_mgr%push_context('diagstate_reset', 'diagstate_mod.F90')
      rc = CC_SUCCESS

      ! Reset diagnostic arrays to zero
      if (allocated(this%AOD550)) this%AOD550 = 0.0_fp
      if (allocated(this%AOD380)) this%AOD380 = 0.0_fp
      this%dust_total_flux = 0.0_fp
      this%seasalt_total_flux = 0.0_fp
      if (allocated(this%drydep_velocity)) this%drydep_velocity = 0.0_fp
      if (allocated(this%drydep_frequency)) this%drydep_frequency = 0.0_fp

      call error_mgr%pop_context()
   end subroutine diagstate_reset

   !> \brief Check if diagnostic state is allocated
   !!
   !! \param[in] this The diagnostic state object
   !! \return True if allocated, false otherwise
   logical function diagstate_is_allocated(this)
      implicit none
      class(DiagStateType), intent(in) :: this

      diagstate_is_allocated = allocated(this%AOD550) .or. &
                               allocated(this%AOD380) .or. &
                               allocated(this%TOMSAI) .or. &
                               allocated(this%drydep_velocity) .or. &
                               allocated(this%drydep_frequency)
   end function diagstate_is_allocated

   !> \brief Get memory usage of diagnostic state
   !!
   !! \param[in] this The diagnostic state object
   !! \return Memory usage in bytes
   integer(8) function diagstate_get_memory_usage(this)
      implicit none
      class(DiagStateType), intent(in) :: this

      integer(8) :: total_bytes

      total_bytes = 0_8

      if (allocated(this%AOD550)) total_bytes = total_bytes + size(this%AOD550) * 8_8
      if (allocated(this%AOD380)) total_bytes = total_bytes + size(this%AOD380) * 8_8
      if (allocated(this%TOMSAI)) total_bytes = total_bytes + size(this%TOMSAI) * 8_8
      if (allocated(this%drydep_velocity)) total_bytes = total_bytes + size(this%drydep_velocity) * 8_8
      if (allocated(this%drydep_frequency)) total_bytes = total_bytes + size(this%drydep_frequency) * 8_8

      diagstate_get_memory_usage = total_bytes
   end function diagstate_get_memory_usage

   !> \brief Print summary of diagnostic state
   !!
   !! \param[in] this The diagnostic state object
   !! \param[in] unit Optional unit to write to (default stdout)
   subroutine diagstate_print_summary(this, unit)
      implicit none
      class(DiagStateType), intent(in) :: this
      integer, optional, intent(in) :: unit

      integer :: out_unit

      out_unit = 6  ! Default to stdout
      if (present(unit)) out_unit = unit

      write(out_unit, '(A)') '=== Diagnostic State Summary ==='
      write(out_unit, '(A, I0)') 'Number of levels: ', this%nlevs
      write(out_unit, '(A, I0)') 'Dry deposition species: ', this%nspecies_drydep
      write(out_unit, '(A, L1)') 'Valid: ', this%is_valid
      write(out_unit, '(A, L1)') 'AOD550 allocated: ', allocated(this%AOD550)
      write(out_unit, '(A, L1)') 'AOD380 allocated: ', allocated(this%AOD380)
      write(out_unit, '(A, L1)') 'TOMSAI allocated: ', allocated(this%TOMSAI)
      write(out_unit, '(A, L1)') 'Dry dep velocity allocated: ', allocated(this%drydep_velocity)
      write(out_unit, '(A, L1)') 'Dry dep frequency allocated: ', allocated(this%drydep_frequency)
      write(out_unit, '(A, I0, A)') 'Memory usage: ', this%get_memory_usage(), ' bytes'
      write(out_unit, '(A)') '================================='
   end subroutine diagstate_print_summary

   !> \brief Print diagnostic state information
   !!
   !! This subroutine prints detailed information about the diagnostic state
   !! for debugging and diagnostics.
   !!
   !! \param[in] this The diagnostic state object
   subroutine diagstate_print_info(this)
      implicit none
      class(DiagStateType), intent(in) :: this

      write(*, '(A)') '=== Diagnostic State Information ==='
      write(*, '(A,L1)') 'Initialized: ', this%is_initialized
      write(*, '(A,L1)') 'Valid: ', this%is_valid
      write(*, '(A,I0)') 'Number of levels: ', this%nlevs
      write(*, '(A,I0)') 'Number of dry deposition species: ', this%nspecies_drydep
      write(*, '(A,F12.6)') 'Dust total flux [kg m-2 s-1]: ', this%dust_total_flux
      write(*, '(A,F12.6)') 'Sea salt total flux [kg m-2 s-1]: ', this%seasalt_total_flux
      write(*, '(A,F12.3)') 'Briggs plume rise height [m]: ', this%briggs_plumerise_height
      write(*, '(A,F12.3)') 'Sofiev plume rise height [m]: ', this%sofiev_plumerise_height

      if (allocated(this%AOD550)) then
         write(*, '(A,I0)') 'AOD550 array size: ', size(this%AOD550)
      else
         write(*, '(A)') 'AOD550: Not allocated'
      endif

      if (allocated(this%drydep_velocity)) then
         write(*, '(A,I0)') 'Dry deposition velocity array size: ', size(this%drydep_velocity)
      else
         write(*, '(A)') 'Dry deposition velocity: Not allocated'
      endif

      write(*, '(A)') '=================================='

   end subroutine diagstate_print_info

end module DiagState_Mod
