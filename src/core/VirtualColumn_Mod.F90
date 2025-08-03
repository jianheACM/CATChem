!> \file VirtualColumn_Mod.F90
!! \brief Standalone module for VirtualColumnType - simple data container
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! \details This module provides a simple virtual column data container
!! that is independent of StateManager to avoid circular dependencies.
!! It holds column data and metadata, but delegates all grid operations
!! to higher-level modules.

module VirtualColumn_Mod
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   implicit none
   private

   public :: VirtualColumnType

   !> \brief Virtual column data container for process-level column virtualization
   !! \details
   !! This is a simple data container that holds:
   !! - Column data arrays (met, chem, emissions)
   !! - Grid position and metadata
   !! - Basic accessor methods
   !!
   !! It does NOT know about StateManager to avoid circular dependencies.
   !! Grid operations are handled by ColumnInterface_Mod or StateManager_Mod.
   type :: VirtualColumnType
      private

      ! Column data arrays
      real(fp), allocatable :: met_data(:)        !< Meteorology column (nlev)
      real(fp), allocatable :: chem_data(:,:)     !< Chemistry column (nlev, nspec_chem)
      real(fp), allocatable :: emis_data(:,:)     !< Emissions column (nlev, nspec_emis)

      ! Dimensions
      integer :: nlev = 0
      integer :: nspec_chem = 0
      integer :: nspec_emis = 0

      ! Grid position and metadata
      real(fp) :: lat = 0.0_fp
      real(fp) :: lon = 0.0_fp
      real(fp) :: area = 0.0_fp
      integer :: grid_i = 1
      integer :: grid_j = 1

      ! State
      logical :: is_valid = .false.

   contains
      procedure :: init => virtual_column_init
      procedure :: get_met_field => virtual_column_get_met_field
      procedure :: get_chem_field => virtual_column_get_chem_field
      procedure :: get_emis_field => virtual_column_get_emis_field
      procedure :: set_met_field => virtual_column_set_met_field
      procedure :: set_chem_field => virtual_column_set_chem_field
      procedure :: set_emis_field => virtual_column_set_emis_field
      procedure :: get_position => virtual_column_get_position
      procedure :: get_metadata => virtual_column_get_metadata
      procedure :: get_dimensions => virtual_column_get_dimensions
      procedure :: is_initialized => virtual_column_is_initialized
      procedure :: cleanup => virtual_column_cleanup
   end type VirtualColumnType

contains

   !> \brief Initialize virtual column with dimensions and grid position
   !! \details Simple initialization - no StateManager dependency
   subroutine virtual_column_init(this, nlev, nspec_chem, nspec_emis, &
                                  grid_i, grid_j, lat, lon, area, rc)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: nlev, nspec_chem, nspec_emis
      integer, intent(in) :: grid_i, grid_j
      real(fp), intent(in) :: lat, lon, area
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Store dimensions
      this%nlev = nlev
      this%nspec_chem = nspec_chem
      this%nspec_emis = nspec_emis

      ! Store position/metadata
      this%grid_i = grid_i
      this%grid_j = grid_j
      this%lat = lat
      this%lon = lon
      this%area = area

      ! Allocate data arrays
      if (nlev > 0) then
         allocate(this%met_data(nlev), stat=rc)
         if (rc /= 0) then
            rc = CC_FAILURE
            return
         endif
         this%met_data = 0.0_fp
      endif

      if (nlev > 0 .and. nspec_chem > 0) then
         allocate(this%chem_data(nlev, nspec_chem), stat=rc)
         if (rc /= 0) then
            rc = CC_FAILURE
            return
         endif
         this%chem_data = 0.0_fp
      endif

      if (nlev > 0 .and. nspec_emis > 0) then
         allocate(this%emis_data(nlev, nspec_emis), stat=rc)
         if (rc /= 0) then
            rc = CC_FAILURE
            return
         endif
         this%emis_data = 0.0_fp
      endif

      this%is_valid = .true.
      rc = CC_SUCCESS
   end subroutine virtual_column_init

   !> \brief Get meteorological field value at level k
   function virtual_column_get_met_field(this, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(in) :: k
      real(fp) :: value

      if (allocated(this%met_data) .and. k >= 1 .and. k <= this%nlev) then
         value = this%met_data(k)
      else
         value = 0.0_fp
      endif
   end function virtual_column_get_met_field

   !> \brief Get chemical species concentration at level k
   function virtual_column_get_chem_field(this, ispec, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(in) :: ispec, k
      real(fp) :: value

      if (allocated(this%chem_data) .and. &
          k >= 1 .and. k <= this%nlev .and. &
          ispec >= 1 .and. ispec <= this%nspec_chem) then
         value = this%chem_data(k, ispec)
      else
         value = 0.0_fp
      endif
   end function virtual_column_get_chem_field

   !> \brief Get emission flux at level k
   function virtual_column_get_emis_field(this, ispec, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(in) :: ispec, k
      real(fp) :: value

      if (allocated(this%emis_data) .and. &
          k >= 1 .and. k <= this%nlev .and. &
          ispec >= 1 .and. ispec <= this%nspec_emis) then
         value = this%emis_data(k, ispec)
      else
         value = 0.0_fp
      endif
   end function virtual_column_get_emis_field

   !> \brief Set meteorological field value at level k
   subroutine virtual_column_set_met_field(this, k, value)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: k
      real(fp), intent(in) :: value

      if (allocated(this%met_data) .and. k >= 1 .and. k <= this%nlev) then
         this%met_data(k) = value
      endif
   end subroutine virtual_column_set_met_field

   !> \brief Set chemical species concentration at level k
   subroutine virtual_column_set_chem_field(this, ispec, k, value)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: ispec, k
      real(fp), intent(in) :: value

      if (allocated(this%chem_data) .and. &
          k >= 1 .and. k <= this%nlev .and. &
          ispec >= 1 .and. ispec <= this%nspec_chem) then
         this%chem_data(k, ispec) = value
      endif
   end subroutine virtual_column_set_chem_field

   !> \brief Set emission flux at level k
   subroutine virtual_column_set_emis_field(this, ispec, k, value)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(in) :: ispec, k
      real(fp), intent(in) :: value

      if (allocated(this%emis_data) .and. &
          k >= 1 .and. k <= this%nlev .and. &
          ispec >= 1 .and. ispec <= this%nspec_emis) then
         this%emis_data(k, ispec) = value
      endif
   end subroutine virtual_column_set_emis_field

   !> \brief Get grid position
   subroutine virtual_column_get_position(this, grid_i, grid_j)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(out) :: grid_i, grid_j

      grid_i = this%grid_i
      grid_j = this%grid_j
   end subroutine virtual_column_get_position

   !> \brief Get column metadata
   subroutine virtual_column_get_metadata(this, lat, lon, area)
      class(VirtualColumnType), intent(in) :: this
      real(fp), intent(out) :: lat, lon, area

      lat = this%lat
      lon = this%lon
      area = this%area
   end subroutine virtual_column_get_metadata

   !> \brief Get column dimensions
   subroutine virtual_column_get_dimensions(this, nlev, nspec_chem, nspec_emis)
      class(VirtualColumnType), intent(in) :: this
      integer, intent(out) :: nlev, nspec_chem, nspec_emis

      nlev = this%nlev
      nspec_chem = this%nspec_chem
      nspec_emis = this%nspec_emis
   end subroutine virtual_column_get_dimensions

   !> \brief Check if column is initialized
   function virtual_column_is_initialized(this) result(initialized)
      class(VirtualColumnType), intent(in) :: this
      logical :: initialized

      initialized = this%is_valid
   end function virtual_column_is_initialized

   !> \brief Clean up virtual column
   subroutine virtual_column_cleanup(this)
      class(VirtualColumnType), intent(inout) :: this

      if (allocated(this%met_data)) deallocate(this%met_data)
      if (allocated(this%chem_data)) deallocate(this%chem_data)
      if (allocated(this%emis_data)) deallocate(this%emis_data)

      this%nlev = 0
      this%nspec_chem = 0
      this%nspec_emis = 0
      this%is_valid = .false.
   end subroutine virtual_column_cleanup

end module VirtualColumn_Mod
