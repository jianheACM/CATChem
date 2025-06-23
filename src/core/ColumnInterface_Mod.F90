!========================================================================
!! Virtual Column Interface Module
!!
!! Provides a column-based interface to multi-dimensional atmospheric data.
!! This allows column-based processes to work transparently with both 1D
!! (column) and multi-dimensional (2D/3D) host models.
!!
!! Key features:
!! - Zero-copy access via pointers
!! - OpenMP/OpenACC friendly
!! - Support for 1D, 2D, and 3D grids
!!
!! @author CATChem Development Team
!! @date 2025
!! @version 1.0
!!
!! This module provides a column-based interface to multi-dimensional atmospheric data,
!! allowing column-based processes to work transparently with both 1D (column) and
!! multi-dimensional (2D/3D) host models.
!!
!! @details
!! - Zero-copy access via pointers
!! - OpenMP/OpenACC friendly
!! - Support for 1D, 2D, and 3D grids
!!
module ColumnInterface_Mod
   use Precision_Mod, only: fp
   use Error_Mod, only: CC_SUCCESS, CC_FAILURE
   use State_Mod, only: StateContainerType
   use ChemState_Mod, only: ChemStateType
   use MetState_Mod, only: MetStateType
   use EmisState_Mod, only: EmisStateType

   implicit none
   private

   public :: ColumnViewType
   public :: create_column_view
   public :: VirtualColumnType
   public :: ColumnProcessorType

   !========================================================================
   !! Column View Type
   !!
   !! Provides column-based access to multi-dimensional state data.
   !! Supports both direct column models and virtual columns from 3D data.
   !========================================================================

   !> \brief Column view abstraction for multi-dimensional state data.
   !! \details
   !! The ColumnViewType provides a unified interface for accessing and manipulating
   !! a single atmospheric column within a multi-dimensional (1D/2D/3D) model grid.
   !! It enables column-based processes to operate transparently on both column models
   !! and virtual columns extracted from 2D/3D data, supporting zero-copy pointer access
   !! for performance and compatibility with OpenMP/OpenACC parallelism.
   !!
   !! **Members:**
   !! - grid_type: Grid type (1=column, 2=2D, 3=3D)
   !! - nx, ny, nlev: Grid dimensions
   !! - current_i, current_j: Current column indices
   !! - container: Pointer to the state container
   !! - met_column: Pointer to meteorology column (nlev)
   !! - chem_column: Pointer to chemistry column (nlev, nspecies)
   !! - emis_column: Pointer to emissions column (nlev, nspecies)
   !!
   !! **Usage:**
   !! 1. Initialize with `init()`
   !! 2. Set column position with `set_column_position()`
   !! 3. Access or modify column data via pointers
   !!
   type :: ColumnViewType
      private

      ! Grid information
      integer :: grid_type  !< Grid type: 1=column, 2=2D, 3=3D
      integer :: nx         !< Number of grid points in x-direction
      integer :: ny         !< Number of grid points in y-direction
      integer :: nlev       !< Number of vertical levels
      integer :: current_i  !< Current column i-index
      integer :: current_j  !< Current column j-index

      ! State container reference
      type(StateContainerType), pointer :: container => null() !< Pointer to state container

      ! Column data pointers (for current column)
      real(fp), pointer :: met_column(:) => null()            !< Pointer to meteorology column (nlev)
      real(fp), pointer :: chem_column(:,:) => null()         !< Pointer to chemistry column (nlev, nspecies)
      real(fp), pointer :: emis_column(:,:) => null()         !< Pointer to emissions column (nlev, nspecies)

   contains
      procedure :: init => column_view_init
      procedure :: set_column_position => column_view_set_position
      procedure :: get_met_column_ptr => column_view_get_met_ptr
      procedure :: get_chem_column_ptr => column_view_get_chem_ptr
      procedure :: get_emis_column_ptr => column_view_get_emis_ptr
      procedure :: apply_column_tendency => column_view_apply_tendency
      procedure :: cleanup => column_view_cleanup

   end type ColumnViewType

   !> \brief Virtual column abstraction for process-level column virtualization
   !!
   !! This type provides a completely virtualized column interface where processes
   !! see only 1D column data, but the underlying grid can be 1D, 2D, or 3D.
   !! The virtual column handles all the complexity of 3D grid management internally.
   type :: VirtualColumnType
      private

      ! Column data (always 1D from process perspective)
      real(fp), pointer :: met_data(:) => null()      !< 1D meteorology data
      real(fp), pointer :: chem_data(:,:) => null()   !< 1D chemistry data (nlev, nspec)
      real(fp), pointer :: emis_data(:,:) => null()   !< 1D emission data (nlev, nspec)

      ! Column metadata
      integer :: nlev = 0              !< Number of vertical levels
      integer :: nspec_chem = 0        !< Number of chemical species
      integer :: nspec_emis = 0        !< Number of emission species
      real(fp) :: lat = 0.0_fp         !< Column latitude
      real(fp) :: lon = 0.0_fp         !< Column longitude
      real(fp) :: area = 0.0_fp        !< Column area [m²]

      ! Grid position (hidden from processes)
      integer :: grid_i = 1            !< Grid i-index
      integer :: grid_j = 1            !< Grid j-index
      type(StateContainerType), pointer :: container => null()

      logical :: is_valid = .false.    !< Column validity flag

   contains
      procedure :: init => virtual_column_init
      procedure :: update_from_grid => virtual_column_update_from_grid
      procedure :: apply_to_grid => virtual_column_apply_to_grid
      procedure :: get_met_field => virtual_column_get_met_field
      procedure :: get_chem_field => virtual_column_get_chem_field
      procedure :: get_emis_field => virtual_column_get_emis_field
      procedure :: set_met_field => virtual_column_set_met_field
      procedure :: set_chem_field => virtual_column_set_chem_field
      procedure :: set_emis_field => virtual_column_set_emis_field
      procedure :: get_metadata => virtual_column_get_metadata
      procedure :: cleanup => virtual_column_cleanup
   end type VirtualColumnType

   !> \brief Column processor for managing multiple virtual columns
   !!
   !! This type manages a collection of virtual columns and provides
   !! batch processing capabilities while maintaining the column virtualization.
   type :: ColumnProcessorType
      private

      type(VirtualColumnType), allocatable :: columns(:)  !< Array of virtual columns
      integer :: n_columns = 0                            !< Number of columns
      integer :: current_column = 1                       !< Current column index
      type(StateContainerType), pointer :: container => null()

   contains
      procedure :: init => processor_init
      procedure :: add_column => processor_add_column
      procedure :: process_all => processor_process_all
      procedure :: get_column => processor_get_column
      procedure :: get_current_column => processor_get_current_column
      procedure :: next_column => processor_next_column
      procedure :: reset => processor_reset
      procedure :: cleanup => processor_cleanup
   end type ColumnProcessorType

contains

   !========================================================================
   !! Initialize column view for a state container
   !========================================================================
   subroutine column_view_init(this, container, rc)
      class(ColumnViewType), intent(inout) :: this
      type(StateContainerType), intent(in), target :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state

      rc = CC_SUCCESS

      ! Store container reference
      this%container => container

      ! Get grid dimensions from MetState
      met_state => container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         print *, 'ERROR: MetState not available for column view'
         rc = CC_FAILURE
         return
      endif

      ! Determine grid type and dimensions
      call met_state%get_dimensions(this%nx, this%ny, this%nlev)

      if (this%nx == 1 .and. this%ny == 1) then
         this%grid_type = 1  ! Pure column model
      else if (this%ny == 1) then
         this%grid_type = 2  ! 2D (x-z) model
      else
         this%grid_type = 3  ! 3D (x-y-z) model
      endif

      ! Initialize position
      this%current_i = 1
      this%current_j = 1

   end subroutine column_view_init

   !========================================================================
   !! Set current column position (for multi-dimensional grids)
   !========================================================================
   subroutine column_view_set_position(this, i, j, rc)
      class(ColumnViewType), intent(inout) :: this
      integer, intent(in) :: i, j
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Validate position
      if (i < 1 .or. i > this%nx .or. j < 1 .or. j > this%ny) then
         print *, 'ERROR: Invalid column position:', i, j
         rc = CC_FAILURE
         return
      endif

      this%current_i = i
      this%current_j = j

      ! Update pointers to current column - TODO: implement if needed
      rc = CC_SUCCESS

   end subroutine column_view_set_position

   !========================================================================
   !! Get meteorological column pointer
   !========================================================================
   function column_view_get_met_ptr(this, field_name) result(column_ptr)
      class(ColumnViewType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      real(fp), pointer :: column_ptr(:)

      type(MetStateType), pointer :: met_state

      column_ptr => null()

      if (.not. associated(this%container)) return

      met_state => this%container%get_met_state_ptr()
      if (.not. associated(met_state)) return

      ! For current column model, MetState fields are already 1D
      column_ptr => met_state%get_field_ptr(field_name)

   end function column_view_get_met_ptr

   !========================================================================
   !! Get chemistry column pointer for a species
   !========================================================================
   function column_view_get_chem_ptr(this, species_name) result(column_ptr)
      class(ColumnViewType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      real(fp), pointer :: column_ptr(:)

      type(ChemStateType), pointer :: chem_state

      column_ptr => null()

      if (.not. associated(this%container)) return

      chem_state => this%container%get_chem_state_ptr()
      if (.not. associated(chem_state)) return

      ! For current column model, ChemState species are already 1D
      column_ptr => chem_state%get_species_ptr(species_name)

   end function column_view_get_chem_ptr

   !========================================================================
   !! Get emissions column pointer for a species
   !========================================================================
   function column_view_get_emis_ptr(this, species_name) result(column_ptr)
      class(ColumnViewType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      real(fp), pointer :: column_ptr(:)

      type(EmisStateType), pointer :: emis_state

      column_ptr => null()

      if (.not. associated(this%container)) return

      emis_state => this%container%get_emis_state_ptr()
      if (.not. associated(emis_state)) return

      ! For current column model, EmisState species are already 1D
      column_ptr => emis_state%get_emission_ptr(species_name)

   end function column_view_get_emis_ptr

   !========================================================================
   !! Apply tendency to current column and accumulate in emissions
   !========================================================================
   subroutine column_view_apply_tendency(this, species_name, tendency_column, &
                                        dt, layer_height, rc)
      class(ColumnViewType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      real(fp), intent(in) :: tendency_column(:)  ! Tendency flux (kg/m²/s)
      real(fp), intent(in) :: dt
      real(fp), intent(in), optional :: layer_height
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      real(fp), pointer :: chem_column(:), emis_column(:)
      real(fp) :: height
      integer :: k

      rc = CC_SUCCESS

      height = 50.0_fp  ! Default surface layer height
      if (present(layer_height)) height = layer_height

      ! Get pointers to current column data
      chem_column => this%get_chem_column_ptr(species_name)
      emis_column => this%get_emis_column_ptr(species_name)

      if (.not. associated(chem_column) .or. .not. associated(emis_column)) then
         print *, 'ERROR: Could not get column pointers for ', trim(species_name)
         rc = CC_FAILURE
         return
      endif

      ! Get chemistry state for species metadata
      chem_state => this%container%get_chem_state_ptr()

      ! Apply tendency to each vertical level
      do k = 1, this%nlev
         ! Apply tendency directly to chemistry column
         chem_column(k) = chem_column(k) + tendency_column(k) * dt

         ! Accumulate in emissions for diagnostics
         emis_column(k) = emis_column(k) + tendency_column(k)
      end do

   end subroutine column_view_apply_tendency

   !========================================================================
   !! Cleanup column view
   !========================================================================
   subroutine column_view_cleanup(this)
      class(ColumnViewType), intent(inout) :: this

      ! Nullify pointers
      this%container => null()
      this%met_column => null()
      this%chem_column => null()
      this%emis_column => null()

   end subroutine column_view_cleanup

   !========================================================================
   !! Convenience function to create column view
   !========================================================================
   function create_column_view(container) result(col_view)
      type(StateContainerType), intent(in), target :: container
      type(ColumnViewType) :: col_view

      integer :: rc

      call col_view%init(container, rc)
      if (rc /= CC_SUCCESS) then
         print *, 'WARNING: Failed to initialize column view'
      endif

   end function create_column_view

   !========================================================================
   ! VirtualColumnType Implementation
   !========================================================================

   !> \brief Initialize virtual column
   subroutine virtual_column_init(this, container, grid_i, grid_j, rc)
      class(VirtualColumnType), intent(inout) :: this
      type(StateContainerType), target, intent(in) :: container
      integer, intent(in) :: grid_i, grid_j
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(EmisStateType), pointer :: emis_state

      rc = CC_SUCCESS

      this%container => container
      this%grid_i = grid_i
      this%grid_j = grid_j

      ! Get state pointers
      met_state => container%get_met_state_ptr()
      chem_state => container%get_chem_state_ptr()
      emis_state => container%get_emis_state_ptr()

      if (.not. associated(met_state)) then
         rc = CC_FAILURE
         return
      endif

      ! Get dimensions from met state
      call met_state%get_dimensions(nx_dummy, ny_dummy, this%nlev)

      ! Set up data pointers (this would extract column data from 3D grids)
      call this%update_from_grid(rc)

      this%is_valid = (rc == CC_SUCCESS)

   end subroutine virtual_column_init

   !> \brief Update virtual column data from 3D grid
   subroutine virtual_column_update_from_grid(this, rc)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      integer :: nx_dummy, ny_dummy

      rc = CC_SUCCESS

      if (.not. associated(this%container)) then
         rc = CC_FAILURE
         return
      endif

      met_state => this%container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         rc = CC_FAILURE
         return
      endif

      ! Extract column data from 3D grid
      call met_state%get_dimensions(nx_dummy, ny_dummy, this%nlev)

      ! Get column pointers for meteorology
      this%met_data => met_state%get_column_ptr('temperature', this%grid_i, this%grid_j)

      ! Similar for chemistry and emissions...
      ! this%chem_data => chem_state%get_column_ptr(this%grid_i, this%grid_j)
      ! this%emis_data => emis_state%get_column_ptr(this%grid_i, this%grid_j)

   end subroutine virtual_column_update_from_grid

   !> \brief Apply virtual column changes back to 3D grid
   subroutine virtual_column_apply_to_grid(this, rc)
      class(VirtualColumnType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Apply any changes made to the virtual column back to the 3D grid
      ! This handles the reverse mapping from 1D column to 3D grid

      ! Implementation would copy modified data back to appropriate grid locations

   end subroutine virtual_column_apply_to_grid

   !> \brief Get meteorological field value at level k
   function virtual_column_get_met_field(this, field_name, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: k
      real(fp) :: value

      ! For now, return from met_data array
      if (associated(this%met_data) .and. k >= 1 .and. k <= this%nlev) then
         value = this%met_data(k)
      else
         value = 0.0_fp
      endif

   end function virtual_column_get_met_field

   !> \brief Get chemical species concentration at level k
   function virtual_column_get_chem_field(this, species_name, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(in) :: k
      real(fp) :: value

      ! Implementation would look up species by name and return concentration
      value = 0.0_fp  ! Placeholder

   end function virtual_column_get_chem_field

   !> \brief Get emission flux at level k
   function virtual_column_get_emis_field(this, species_name, k) result(value)
      class(VirtualColumnType), intent(in) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(in) :: k
      real(fp) :: value

      ! Implementation would look up species by name and return emission flux
      value = 0.0_fp  ! Placeholder

   end function virtual_column_get_emis_field

   !> \brief Set meteorological field value at level k
   subroutine virtual_column_set_met_field(this, field_name, k, value)
      class(VirtualColumnType), intent(inout) :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in) :: k
      real(fp), intent(in) :: value

      if (associated(this%met_data) .and. k >= 1 .and. k <= this%nlev) then
         this%met_data(k) = value
      endif

   end subroutine virtual_column_set_met_field

   !> \brief Set chemical species concentration at level k
   subroutine virtual_column_set_chem_field(this, species_name, k, value)
      class(VirtualColumnType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(in) :: k
      real(fp), intent(in) :: value

      ! Implementation would look up species by name and set concentration

   end subroutine virtual_column_set_chem_field

   !> \brief Set emission flux at level k
   subroutine virtual_column_set_emis_field(this, species_name, k, value)
      class(VirtualColumnType), intent(inout) :: this
      character(len=*), intent(in) :: species_name
      integer, intent(in) :: k
      real(fp), intent(in) :: value

      ! Implementation would look up species by name and set emission flux

   end subroutine virtual_column_set_emis_field

   !> \brief Get column metadata
   subroutine virtual_column_get_metadata(this, lat, lon, area, nlev)
      class(VirtualColumnType), intent(in) :: this
      real(fp), intent(out) :: lat, lon, area
      integer, intent(out) :: nlev

      lat = this%lat
      lon = this%lon
      area = this%area
      nlev = this%nlev

   end subroutine virtual_column_get_metadata

   !> \brief Clean up virtual column
   subroutine virtual_column_cleanup(this)
      class(VirtualColumnType), intent(inout) :: this

      this%met_data => null()
      this%chem_data => null()
      this%emis_data => null()
      this%container => null()
      this%is_valid = .false.

   end subroutine virtual_column_cleanup

   !========================================================================
   ! ColumnProcessorType Implementation
   !========================================================================

   !> \brief Initialize column processor
   subroutine processor_init(this, container, max_columns, rc)
      class(ColumnProcessorType), intent(inout) :: this
      type(StateContainerType), target, intent(in) :: container
      integer, intent(in) :: max_columns
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      this%container => container
      allocate(this%columns(max_columns))
      this%n_columns = 0
      this%current_column = 1

   end subroutine processor_init

   !> \brief Add a column to the processor
   subroutine processor_add_column(this, grid_i, grid_j, rc)
      class(ColumnProcessorType), intent(inout) :: this
      integer, intent(in) :: grid_i, grid_j
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (this%n_columns >= size(this%columns)) then
         rc = CC_FAILURE
         return
      endif

      this%n_columns = this%n_columns + 1
      call this%columns(this%n_columns)%init(this%container, grid_i, grid_j, rc)

   end subroutine processor_add_column

   !> \brief Process all columns with a given procedure
   subroutine processor_process_all(this, process_proc, rc)
      class(ColumnProcessorType), intent(inout) :: this
      procedure(column_process_interface) :: process_proc
      integer, intent(out) :: rc

      integer :: i, local_rc

      rc = CC_SUCCESS

      do i = 1, this%n_columns
         call process_proc(this%columns(i), local_rc)
         if (local_rc /= CC_SUCCESS) then
            rc = local_rc
            ! Continue processing other columns or exit based on error policy
         endif
      enddo

   end subroutine processor_process_all

   !> \brief Get specific column by index
   function processor_get_column(this, index) result(column)
      class(ColumnProcessorType), intent(in) :: this
      integer, intent(in) :: index
      type(VirtualColumnType) :: column

      if (index >= 1 .and. index <= this%n_columns) then
         column = this%columns(index)
      endif

   end function processor_get_column

   !> \brief Get current column
   function processor_get_current_column(this) result(column)
      class(ColumnProcessorType), intent(in) :: this
      type(VirtualColumnType) :: column

      column = this%get_column(this%current_column)

   end function processor_get_current_column

   !> \brief Move to next column
   subroutine processor_next_column(this, rc)
      class(ColumnProcessorType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      if (this%current_column < this%n_columns) then
         this%current_column = this%current_column + 1
      else
         rc = CC_FAILURE  ! No more columns
      endif

   end subroutine processor_next_column

   !> \brief Reset to first column
   subroutine processor_reset(this)
      class(ColumnProcessorType), intent(inout) :: this

      this%current_column = 1

   end subroutine processor_reset

   !> \brief Clean up column processor
   subroutine processor_cleanup(this)
      class(ColumnProcessorType), intent(inout) :: this

      integer :: i

      if (allocated(this%columns)) then
         do i = 1, this%n_columns
            call this%columns(i)%cleanup()
         enddo
         deallocate(this%columns)
      endif

      this%container => null()
      this%n_columns = 0

   end subroutine processor_cleanup

   !> \brief Abstract interface for column processing procedures
   abstract interface
      subroutine column_process_interface(column, rc)
         import :: VirtualColumnType
         type(VirtualColumnType), intent(inout) :: column
         integer, intent(out) :: rc
      end subroutine column_process_interface
   end interface

end module ColumnInterface_Mod
