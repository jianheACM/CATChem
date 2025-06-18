!> \file DustProcess_Mod.F90
!! \brief Dust emission process implementation using StateContainer and ProcessInterface
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module implements dust emission processes that extend the ProcessInterface
!! and integrate with the StateContainer architecture. It provides a framework
!! for multiple dust emission schemes (AFWA, Fengsha, GOCART) while leveraging
!! common utilities from DustCommon_Mod.
!!
!! The process handles:
!! - Configuration management through ConfigDataType
!! - State management through StateContainer components
!! - Error handling through ErrorManager
!! - Species mapping and emission calculations
!! - Grid consistency validation
!!
!! \see DustCommon_Mod for low-level utility functions
!! \see ProcessInterface_Mod for base class functionality
!! \see StateContainer_Mod for state management
!!
module DustProcess_Mod

   use precision_mod
   use ProcessInterface_Mod, only: ProcessInterface
   use state_mod, only: StateContainerType, ConfigDataType, MetStateType, &
                       ChemStateType, EmisStateType, DiagStateType
   use DustCommon_Mod
   use error_mod

   implicit none
   private

   public :: DustProcessType
   public :: DUST_SCHEME_AFWA, DUST_SCHEME_FENGSHA, DUST_SCHEME_GOCART

   ! Dust scheme identifiers
   integer, parameter :: DUST_SCHEME_AFWA = 1
   integer, parameter :: DUST_SCHEME_FENGSHA = 2
   integer, parameter :: DUST_SCHEME_GOCART = 3

   !> \brief Dust emission process implementation
   !!
   !! This type extends ProcessInterface to provide dust emission calculations
   !! integrated with the StateContainer architecture. It supports multiple
   !! dust emission schemes and manages all state through the container.
   !!
   type, extends(ProcessInterface) :: DustProcessType
      private

      ! Process-specific configuration
      integer :: dust_scheme = DUST_SCHEME_AFWA    !< Selected dust scheme
      real(fp) :: min_wind_speed = 3.0_fp          !< Minimum wind speed for emission [m/s]
      real(fp) :: max_soil_moisture = 0.5_fp       !< Maximum soil moisture for emission [-]
      logical :: use_soil_moisture = .true.        !< Enable soil moisture correction
      logical :: use_roughness_correction = .true. !< Enable roughness correction

      ! Grid information
      integer :: n_lon = 0                         !< Number of longitude points
      integer :: n_lat = 0                         !< Number of latitude points
      integer :: n_lev = 0                         !< Number of vertical levels

      ! Species mapping
      integer, allocatable :: species_indices(:)   !< Indices of dust species in chem_state

      ! Diagnostic arrays
      real(fp), allocatable :: friction_velocity(:,:)      !< Surface friction velocity [m/s]
      real(fp), allocatable :: threshold_velocity(:,:)     !< Threshold velocity [m/s]
      real(fp), allocatable :: soil_moisture_factor(:,:)   !< Soil moisture attenuation [-]
      real(fp), allocatable :: total_dust_flux(:,:)        !< Total dust emission [kg/m²/s]

   contains
      ! Required ProcessInterface methods
      procedure :: init => dust_process_init
      procedure :: run => dust_process_run
      procedure :: finalize => dust_process_finalize

      ! Process-specific methods
      procedure :: set_dust_scheme => dust_process_set_scheme
      procedure :: get_dust_scheme => dust_process_get_scheme
      procedure :: calculate_emissions => dust_process_calculate_emissions
      procedure :: apply_soil_moisture_correction => dust_process_soil_moisture
      procedure :: apply_roughness_correction => dust_process_roughness
      procedure :: map_species_to_chemical_state => dust_process_map_species
      procedure :: validate_dust_config => dust_process_validate_config

      ! Scheme-specific implementations
      procedure, private :: run_afwa_scheme => dust_process_run_afwa
      procedure, private :: run_fengsha_scheme => dust_process_run_fengsha
      procedure, private :: run_gocart_scheme => dust_process_run_gocart
   end type DustProcessType

contains

   !> \brief Initialize dust emission process
   !!
   !! Sets up the dust process by reading configuration, validating parameters,
   !! allocating arrays, and mapping species to the chemical state.
   !!
   !! \param[inout] this DustProcessType instance
   !! \param[inout] container StateContainer with all state objects
   !! \param[out] rc Return code
   !!
   subroutine dust_process_init(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ConfigDataType), pointer :: config
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: dust_species_list
      character(len=32) :: dust_species(10)  ! Support up to 10 dust species
      integer :: n_dust_species, i
      logical :: species_available

      ! Initialize
      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_init', 'Initializing dust emission process')

      ! Set process metadata
      this%name = 'DustEmission'
      this%version = '2.0'

      ! Get configuration
      config => container%get_config_ptr()
      if (.not. associated(config)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Configuration not available', rc)
         return
      end if

      ! Read dust-specific configuration
      call config%get_parameter('dust_scheme', this%dust_scheme, default=DUST_SCHEME_AFWA)
      call config%get_parameter('dust_min_wind_speed', this%min_wind_speed, default=3.0_fp)
      call config%get_parameter('dust_max_soil_moisture', this%max_soil_moisture, default=0.5_fp)
      call config%get_parameter('dust_use_soil_moisture', this%use_soil_moisture, default=.true.)
      call config%get_parameter('dust_use_roughness', this%use_roughness_correction, default=.true.)
      call config%get_parameter('dust_species_list', dust_species_list, default='DST01,DST02,DST03,DST04,DST05')
      call config%get_parameter('process_timestep', this%dt, default=3600.0_fp)

      ! Parse dust species list
      call parse_species_list(dust_species_list, dust_species, n_dust_species)
      call this%set_species(dust_species(1:n_dust_species), rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_CONFIGURATION, 'Failed to set dust species', rc)
         return
      end if

      ! Get grid dimensions from meteorological state
      met_state => container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Meteorological state not available', rc)
         return
      end if

      call met_state%get_dimensions(this%n_lon, this%n_lat, this%n_lev)

      ! Validate grid consistency
      if (.not. this%validate_grid_consistency(container, rc)) then
         call error_mgr%report_error(ERROR_GRID_MISMATCH, 'Grid dimensions inconsistent', rc)
         return
      end if

      ! Check species availability
      species_available = this%check_species_availability(container, rc)
      if (.not. species_available) then
         call error_mgr%report_error(ERROR_SPECIES_NOT_FOUND, 'Required dust species not available', rc)
         return
      end if

      ! Map species to chemical state indices
      call this%map_species_to_chemical_state(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_SPECIES_MAPPING, 'Failed to map dust species', rc)
         return
      end if

      ! Allocate diagnostic arrays
      allocate(this%friction_velocity(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%threshold_velocity(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%soil_moisture_factor(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%total_dust_flux(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      ! Initialize arrays
      this%friction_velocity = 0.0_fp
      this%threshold_velocity = 0.0_fp
      this%soil_moisture_factor = 1.0_fp
      this%total_dust_flux = 0.0_fp

      ! Validate process configuration
      call this%validate_dust_config(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_CONFIGURATION, 'Invalid dust configuration', rc)
         return
      end if

      ! Mark as initialized
      this%is_initialized = .true.

      ! Log process information
      call this%log_process_info(container)

      call error_mgr%pop_context()
      return

999   continue
      call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate dust arrays', rc)
      call error_mgr%pop_context()

   end subroutine dust_process_init

   !> \brief Execute dust emission calculations
   !!
   !! Performs the main dust emission calculation based on the selected scheme,
   !! updates emission and diagnostic states, and applies corrections.
   !!
   !! \param[inout] this DustProcessType instance
   !! \param[inout] container StateContainer with all state objects
   !! \param[out] rc Return code
   !!
   subroutine dust_process_run(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_run', 'Executing dust emission calculations')

      ! Check if process is ready
      if (.not. this%is_ready()) then
         call error_mgr%report_error(ERROR_NOT_INITIALIZED, 'Dust process not ready', rc)
         return
      end if

      ! Reset diagnostic arrays
      this%friction_velocity = 0.0_fp
      this%threshold_velocity = 0.0_fp
      this%soil_moisture_factor = 1.0_fp
      this%total_dust_flux = 0.0_fp

      ! Execute scheme-specific calculations
      select case (this%dust_scheme)
      case (DUST_SCHEME_AFWA)
         call this%run_afwa_scheme(container, rc)
      case (DUST_SCHEME_FENGSHA)
         call this%run_fengsha_scheme(container, rc)
      case (DUST_SCHEME_GOCART)
         call this%run_gocart_scheme(container, rc)
      case default
         call error_mgr%report_error(ERROR_INVALID_OPTION, 'Unknown dust scheme', rc)
         return
      end select

      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_COMPUTATION, 'Dust scheme calculation failed', rc)
         return
      end if

      ! Update diagnostic state with computed values
      call this%update_diagnostic_state(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_STATE_UPDATE, 'Failed to update diagnostic state', rc)
         return
      end if

      call error_mgr%pop_context()

   end subroutine dust_process_run

   !> \brief Finalize dust emission process
   !!
   !! Cleans up allocated memory and finalizes the dust process.
   !!
   !! \param[inout] this DustProcessType instance
   !! \param[out] rc Return code
   !!
   subroutine dust_process_finalize(this, rc)
      class(DustProcessType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate arrays
      if (allocated(this%species_indices)) deallocate(this%species_indices)
      if (allocated(this%friction_velocity)) deallocate(this%friction_velocity)
      if (allocated(this%threshold_velocity)) deallocate(this%threshold_velocity)
      if (allocated(this%soil_moisture_factor)) deallocate(this%soil_moisture_factor)
      if (allocated(this%total_dust_flux)) deallocate(this%total_dust_flux)

      ! Deallocate species arrays from parent class
      call this%deallocate_species_arrays(this%species_arrays)

      ! Reset state
      this%is_initialized = .false.
      this%is_active = .false.
      this%n_species = 0

   end subroutine dust_process_finalize

   !> \brief Set dust emission scheme
   subroutine dust_process_set_scheme(this, scheme)
      class(DustProcessType), intent(inout) :: this
      integer, intent(in) :: scheme
      this%dust_scheme = scheme
   end subroutine dust_process_set_scheme

   !> \brief Get dust emission scheme
   function dust_process_get_scheme(this) result(scheme)
      class(DustProcessType), intent(in) :: this
      integer :: scheme
      scheme = this%dust_scheme
   end function dust_process_get_scheme

   !> \brief Calculate dust emissions using selected scheme
   !!
   !! This is a high-level interface that calls the appropriate scheme-specific
   !! calculation method based on the configured dust scheme.
   !!
   subroutine dust_process_calculate_emissions(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! This method delegates to the scheme-specific run method
      call this%run(container, rc)

   end subroutine dust_process_calculate_emissions

   !> \brief Apply soil moisture correction to dust emissions
   !!
   !! Uses the utilities from DustCommon_Mod to calculate soil moisture
   !! correction factors based on the selected parameterization.
   !!
   subroutine dust_process_soil_moisture(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ErrorManagerType), pointer :: error_mgr
      real(fp), pointer :: soil_moisture(:,:), clay_frac(:,:), sand_frac(:,:)
      integer :: i, j

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_soil_moisture', 'Applying soil moisture correction')

      if (.not. this%use_soil_moisture) then
         this%soil_moisture_factor = 1.0_fp
         call error_mgr%pop_context()
         return
      end if

      met_state => container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Met state not available', rc)
         return
      end if

      ! Get soil properties (assuming these methods exist in MetState)
      soil_moisture => met_state%get_field_ptr('soil_moisture')
      clay_frac => met_state%get_field_ptr('clay_fraction')
      sand_frac => met_state%get_field_ptr('sand_fraction')

      if (.not. (associated(soil_moisture) .and. associated(clay_frac) .and. associated(sand_frac))) then
         call error_mgr%report_error(ERROR_FIELD_NOT_FOUND, 'Required soil fields not available', rc)
         return
      end if

      ! Calculate soil moisture correction using Fecan parameterization
      do j = 1, this%n_lat
         do i = 1, this%n_lon
            if (soil_moisture(i,j) <= this%max_soil_moisture) then
               call Fecan_SoilMoisture(clay_frac(i,j), sand_frac(i,j), soil_moisture(i,j), &
                                     this%soil_moisture_factor(i,j), error_mgr, rc)
               if (rc /= CC_SUCCESS) then
                  call error_mgr%report_error(ERROR_COMPUTATION, 'Fecan soil moisture calculation failed', rc)
                  return
               end if
            else
               this%soil_moisture_factor(i,j) = 0.0_fp  ! Suppress emission for wet soils
            end if
         end do
      end do

      call error_mgr%pop_context()

   end subroutine dust_process_soil_moisture

   !> \brief Apply surface roughness correction
   !!
   !! Uses utilities from DustCommon_Mod to apply drag partition correction
   !! based on surface roughness elements.
   !!
   subroutine dust_process_roughness(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ErrorManagerType), pointer :: error_mgr
      real(fp), pointer :: roughness(:,:)
      real(fp) :: drag_partition
      integer :: i, j

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_roughness', 'Applying roughness correction')

      if (.not. this%use_roughness_correction) then
         call error_mgr%pop_context()
         return
      end if

      met_state => container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Met state not available', rc)
         return
      end if

      roughness => met_state%get_field_ptr('surface_roughness')
      if (.not. associated(roughness)) then
         call error_mgr%report_error(ERROR_FIELD_NOT_FOUND, 'Surface roughness not available', rc)
         return
      end if

      ! Apply drag partition correction
      do j = 1, this%n_lat
         do i = 1, this%n_lon
            call MB95_DragPartition(roughness(i,j), drag_partition, error_mgr, rc)
            if (rc /= CC_SUCCESS) then
               call error_mgr%report_error(ERROR_COMPUTATION, 'Drag partition calculation failed', rc)
               return
            end if

            ! Apply correction to friction velocity
            this%friction_velocity(i,j) = this%friction_velocity(i,j) * drag_partition
         end do
      end do

      call error_mgr%pop_context()

   end subroutine dust_process_roughness

   !> \brief Map dust species names to chemical state indices
   !!
   !! Creates a mapping between the dust species handled by this process
   !! and their indices in the chemical state arrays.
   !!
   subroutine dust_process_map_species(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr
      integer :: i, species_idx

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_map_species', 'Mapping dust species')

      chem_state => container%get_chem_state_ptr()
      if (.not. associated(chem_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Chemical state not available', rc)
         return
      end if

      if (allocated(this%species_indices)) deallocate(this%species_indices)
      allocate(this%species_indices(this%n_species), stat=rc)
      if (rc /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate species indices', rc)
         return
      end if

      ! Map each species name to its index in chemical state
      do i = 1, this%n_species
         species_idx = chem_state%get_species_index(this%species_names(i))
         if (species_idx <= 0) then
            call error_mgr%report_error(ERROR_SPECIES_NOT_FOUND, &
                 'Species not found: ' // trim(this%species_names(i)), rc)
            return
         end if
         this%species_indices(i) = species_idx
      end do

      call error_mgr%pop_context()

   end subroutine dust_process_map_species

   !> \brief Validate dust-specific configuration parameters
   subroutine dust_process_validate_config(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_validate_config', 'Validating dust configuration')

      ! Validate dust scheme
      if (this%dust_scheme < DUST_SCHEME_AFWA .or. this%dust_scheme > DUST_SCHEME_GOCART) then
         call error_mgr%report_error(ERROR_INVALID_OPTION, 'Invalid dust scheme', rc)
         return
      end if

      ! Validate physical parameters
      if (this%min_wind_speed < 0.0_fp .or. this%min_wind_speed > 10.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid minimum wind speed', rc)
         return
      end if

      if (this%max_soil_moisture < 0.0_fp .or. this%max_soil_moisture > 1.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid maximum soil moisture', rc)
         return
      end if

      if (this%dt <= 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid process timestep', rc)
         return
      end if

      call error_mgr%pop_context()

   end subroutine dust_process_validate_config

   !> \brief Run AFWA dust emission scheme
   !!
   !! Implements the Air Force Weather Agency (AFWA) dust emission scheme
   !! using utilities from DustCommon_Mod.
   !!
   subroutine dust_process_run_afwa(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(EmisStateType), pointer :: emis_state
      type(ErrorManagerType), pointer :: error_mgr
      real(fp), pointer :: wind_speed(:,:), clay_frac(:,:), sand_frac(:,:)
      real(fp) :: dust_flux, horizontal_flux
      integer :: i, j, k

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_run_afwa', 'Running AFWA dust scheme')

      ! Get state pointers
      met_state => container%get_met_state_ptr()
      emis_state => container%get_emis_state_ptr()

      if (.not. (associated(met_state) .and. associated(emis_state))) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Required states not available', rc)
         return
      end if

      ! Get meteorological fields
      wind_speed => met_state%get_field_ptr('wind_speed_10m')
      clay_frac => met_state%get_field_ptr('clay_fraction')
      sand_frac => met_state%get_field_ptr('sand_fraction')

      if (.not. (associated(wind_speed) .and. associated(clay_frac) .and. associated(sand_frac))) then
         call error_mgr%report_error(ERROR_FIELD_NOT_FOUND, 'Required met fields not available', rc)
         return
      end if

      ! Apply soil moisture correction
      call this%apply_soil_moisture_correction(container, rc)
      if (rc /= CC_SUCCESS) return

      ! Calculate emissions for each grid cell
      do j = 1, this%n_lat
         do i = 1, this%n_lon

            ! Skip if wind speed too low
            if (wind_speed(i,j) < this%min_wind_speed) cycle

            ! Calculate friction velocity (simplified - would use more sophisticated method)
            this%friction_velocity(i,j) = wind_speed(i,j) * 0.04_fp  ! Rough approximation

            ! Calculate threshold velocity using MB97 parameterization
            call MB97_threshold_velocity(100.0_fp, 2650.0_fp, this%threshold_velocity(i,j), error_mgr, rc)
            if (rc /= CC_SUCCESS) return

            ! Skip if below threshold
            if (this%friction_velocity(i,j) < this%threshold_velocity(i,j)) cycle

            ! Calculate horizontal flux using Kawamura parameterization
            call Kawamura_HorizFlux(this%friction_velocity(i,j), this%threshold_velocity(i,j), &
                                   horizontal_flux, error_mgr, rc)
            if (rc /= CC_SUCCESS) return

            ! Calculate dust flux (simplified size distribution)
            dust_flux = horizontal_flux * this%soil_moisture_factor(i,j) * 1.0e-6_fp

            this%total_dust_flux(i,j) = dust_flux

            ! Distribute flux among dust species (simplified)
            do k = 1, this%n_species
               call emis_state%add_emission(this%species_indices(k), i, j, 1, &
                                          dust_flux / real(this%n_species, fp))
            end do

         end do
      end do

      ! Apply roughness correction
      call this%apply_roughness_correction(container, rc)
      if (rc /= CC_SUCCESS) return

      call error_mgr%pop_context()

   end subroutine dust_process_run_afwa

   !> \brief Run Fengsha dust emission scheme
   !!
   !! Placeholder for Fengsha dust emission scheme implementation.
   !!
   subroutine dust_process_run_fengsha(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_run_fengsha', 'Running Fengsha dust scheme')

      ! TODO: Implement Fengsha scheme
      call error_mgr%report_info('Fengsha scheme not yet implemented')

      call error_mgr%pop_context()

   end subroutine dust_process_run_fengsha

   !> \brief Run GOCART dust emission scheme
   !!
   !! Placeholder for GOCART dust emission scheme implementation.
   !!
   subroutine dust_process_run_gocart(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('dust_process_run_gocart', 'Running GOCART dust scheme')

      ! TODO: Implement GOCART scheme
      call error_mgr%report_info('GOCART scheme not yet implemented')

      call error_mgr%pop_context()

   end subroutine dust_process_run_gocart

   !> \brief Update diagnostic state with computed dust values
   subroutine update_diagnostic_state(this, container, rc)
      class(DustProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(DiagStateType), pointer :: diag_state
      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('update_diagnostic_state', 'Updating diagnostic state')

      diag_state => container%get_diag_state_ptr()
      if (.not. associated(diag_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Diagnostic state not available', rc)
         return
      end if

      ! Update diagnostic fields (assuming these methods exist)
      call diag_state%set_field('dust_friction_velocity', this%friction_velocity)
      call diag_state%set_field('dust_threshold_velocity', this%threshold_velocity)
      call diag_state%set_field('dust_soil_moisture_factor', this%soil_moisture_factor)
      call diag_state%set_field('dust_total_flux', this%total_dust_flux)

      call error_mgr%pop_context()

   end subroutine update_diagnostic_state

   !> \brief Parse comma-separated species list into array
   subroutine parse_species_list(species_list, species_array, n_species)
      character(len=*), intent(in) :: species_list
      character(len=32), intent(out) :: species_array(:)
      integer, intent(out) :: n_species

      integer :: i, start, comma_pos
      character(len=len(species_list)) :: remaining

      n_species = 0
      remaining = trim(species_list)
      start = 1

      do while (len_trim(remaining) > 0 .and. n_species < size(species_array))
         comma_pos = index(remaining, ',')

         if (comma_pos > 0) then
            n_species = n_species + 1
            species_array(n_species) = trim(remaining(1:comma_pos-1))
            remaining = remaining(comma_pos+1:)
         else
            n_species = n_species + 1
            species_array(n_species) = trim(remaining)
            exit
         end if
      end do

   end subroutine parse_species_list

end module DustProcess_Mod
