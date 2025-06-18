!> \file SeaSaltProcess_Mod.F90
!! \brief Sea salt emission process implementation using StateContainer and ProcessInterface
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module implements sea salt emission processes that extend the ProcessInterface
!! and integrate with the StateContainer architecture. It provides a framework
!! for multiple sea salt emission schemes (Monahan, Gong, SSLT) while leveraging
!! common utilities from SeaSaltCommon_Mod.
!!
!! The process handles:
!! - Configuration management through ConfigDataType
!! - State management through StateContainer components
!! - Error handling through ErrorManager
!! - Species mapping and emission calculations
!! - Ocean mask and surface conditions
!!
!! \see SeaSaltCommon_Mod for low-level utility functions
!! \see ProcessInterface_Mod for base class functionality
!! \see StateContainer_Mod for state management
!!
module SeaSaltProcess_Mod

   use precision_mod
   use ProcessInterface_Mod, only: ProcessInterface
   use state_mod, only: StateContainerType, ConfigDataType, MetStateType, &
                       ChemStateType, EmisStateType, DiagStateType
   use SeaSaltCommon_Mod
   use error_mod

   implicit none
   private

   public :: SeaSaltProcessType
   public :: SEASALT_SCHEME_MONAHAN, SEASALT_SCHEME_GONG, SEASALT_SCHEME_SSLT

   ! Sea salt scheme identifiers
   integer, parameter :: SEASALT_SCHEME_MONAHAN = 1
   integer, parameter :: SEASALT_SCHEME_GONG = 2
   integer, parameter :: SEASALT_SCHEME_SSLT = 3

   !> \brief Sea salt emission process implementation
   !!
   !! This type extends ProcessInterface to provide sea salt emission calculations
   !! integrated with the StateContainer architecture. It supports multiple
   !! sea salt emission schemes and manages all state through the container.
   !!
   type, extends(ProcessInterface) :: SeaSaltProcessType
      private

      ! Process-specific configuration
      integer :: seasalt_scheme = SEASALT_SCHEME_MONAHAN  !< Selected sea salt scheme
      real(fp) :: min_wind_speed = 2.0_fp                !< Minimum wind speed for emission [m/s]
      real(fp) :: max_wind_speed = 50.0_fp               !< Maximum wind speed for emission [m/s]
      real(fp) :: min_sst = 271.15_fp                    !< Minimum sea surface temperature [K]
      logical :: use_sst_correction = .true.             !< Enable SST correction
      logical :: use_whitecap_correction = .true.        !< Enable whitecap parameterization

      ! Particle size configuration
      real(fp) :: min_particle_radius = 0.1_fp           !< Minimum particle radius [μm]
      real(fp) :: max_particle_radius = 10.0_fp          !< Maximum particle radius [μm]
      integer :: n_size_bins = 5                         !< Number of size bins

      ! Grid information
      integer :: n_lon = 0                               !< Number of longitude points
      integer :: n_lat = 0                               !< Number of latitude points
      integer :: n_lev = 0                               !< Number of vertical levels

      ! Species mapping
      integer, allocatable :: species_indices(:)         !< Indices of sea salt species in chem_state

      ! Size bin configuration
      real(fp), allocatable :: size_bin_centers(:)       !< Particle radius bin centers [μm]
      real(fp), allocatable :: size_bin_widths(:)        !< Particle radius bin widths [μm]

      ! Diagnostic arrays
      real(fp), allocatable :: wind_speed_10m(:,:)       !< 10m wind speed [m/s]
      real(fp), allocatable :: whitecap_coverage(:,:)    !< Whitecap coverage fraction [-]
      real(fp), allocatable :: sea_surface_temp(:,:)     !< Sea surface temperature [K]
      real(fp), allocatable :: total_seasalt_flux(:,:)   !< Total sea salt emission [kg/m²/s]
      logical, allocatable :: ocean_mask(:,:)            !< Ocean/land mask

   contains
      ! Required ProcessInterface methods
      procedure :: init => seasalt_process_init
      procedure :: run => seasalt_process_run
      procedure :: finalize => seasalt_process_finalize

      ! Process-specific methods
      procedure :: set_seasalt_scheme => seasalt_process_set_scheme
      procedure :: get_seasalt_scheme => seasalt_process_get_scheme
      procedure :: calculate_emissions => seasalt_process_calculate_emissions
      procedure :: calculate_whitecap_coverage => seasalt_process_whitecap
      procedure :: apply_sst_correction => seasalt_process_sst_correction
      procedure :: setup_size_bins => seasalt_process_setup_size_bins
      procedure :: map_species_to_chemical_state => seasalt_process_map_species
      procedure :: validate_seasalt_config => seasalt_process_validate_config

      ! Scheme-specific implementations
      procedure, private :: run_monahan_scheme => seasalt_process_run_monahan
      procedure, private :: run_gong_scheme => seasalt_process_run_gong
      procedure, private :: run_sslt_scheme => seasalt_process_run_sslt
   end type SeaSaltProcessType

contains

   !> \brief Initialize sea salt emission process
   !!
   !! Sets up the sea salt process by reading configuration, validating parameters,
   !! allocating arrays, setting up size bins, and mapping species.
   !!
   !! \param[inout] this SeaSaltProcessType instance
   !! \param[inout] container StateContainer with all state objects
   !! \param[out] rc Return code
   !!
   subroutine seasalt_process_init(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ConfigDataType), pointer :: config
      type(MetStateType), pointer :: met_state
      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr
      character(len=256) :: seasalt_species_list
      character(len=32) :: seasalt_species(10)  ! Support up to 10 sea salt species
      integer :: n_seasalt_species, i
      logical :: species_available

      ! Initialize
      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_init', 'Initializing sea salt emission process')

      ! Set process metadata
      this%name = 'SeaSaltEmission'
      this%version = '2.0'

      ! Get configuration
      config => container%get_config_ptr()
      if (.not. associated(config)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Configuration not available', rc)
         return
      end if

      ! Read sea salt-specific configuration
      call config%get_parameter('seasalt_scheme', this%seasalt_scheme, default=SEASALT_SCHEME_MONAHAN)
      call config%get_parameter('seasalt_min_wind_speed', this%min_wind_speed, default=2.0_fp)
      call config%get_parameter('seasalt_max_wind_speed', this%max_wind_speed, default=50.0_fp)
      call config%get_parameter('seasalt_min_sst', this%min_sst, default=271.15_fp)
      call config%get_parameter('seasalt_use_sst_correction', this%use_sst_correction, default=.true.)
      call config%get_parameter('seasalt_use_whitecap', this%use_whitecap_correction, default=.true.)
      call config%get_parameter('seasalt_min_radius', this%min_particle_radius, default=0.1_fp)
      call config%get_parameter('seasalt_max_radius', this%max_particle_radius, default=10.0_fp)
      call config%get_parameter('seasalt_n_size_bins', this%n_size_bins, default=5)
      call config%get_parameter('seasalt_species_list', seasalt_species_list, &
                               default='SSLT01,SSLT02,SSLT03,SSLT04,SSLT05')
      call config%get_parameter('process_timestep', this%dt, default=3600.0_fp)

      ! Parse sea salt species list
      call parse_species_list(seasalt_species_list, seasalt_species, n_seasalt_species)
      call this%set_species(seasalt_species(1:n_seasalt_species), rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_CONFIGURATION, 'Failed to set sea salt species', rc)
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
         call error_mgr%report_error(ERROR_SPECIES_NOT_FOUND, 'Required sea salt species not available', rc)
         return
      end if

      ! Setup size bins
      call this%setup_size_bins(rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_CONFIGURATION, 'Failed to setup size bins', rc)
         return
      end if

      ! Map species to chemical state indices
      call this%map_species_to_chemical_state(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_SPECIES_MAPPING, 'Failed to map sea salt species', rc)
         return
      end if

      ! Allocate diagnostic arrays
      allocate(this%wind_speed_10m(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%whitecap_coverage(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%sea_surface_temp(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%total_seasalt_flux(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      allocate(this%ocean_mask(this%n_lon, this%n_lat), stat=rc)
      if (rc /= 0) goto 999

      ! Initialize arrays
      this%wind_speed_10m = 0.0_fp
      this%whitecap_coverage = 0.0_fp
      this%sea_surface_temp = 273.15_fp
      this%total_seasalt_flux = 0.0_fp
      this%ocean_mask = .false.

      ! Setup ocean mask (simplified - assumes all points are ocean for now)
      this%ocean_mask = .true.

      ! Validate process configuration
      call this%validate_seasalt_config(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_CONFIGURATION, 'Invalid sea salt configuration', rc)
         return
      end if

      ! Mark as initialized
      this%is_initialized = .true.

      ! Log process information
      call this%log_process_info(container)

      call error_mgr%pop_context()
      return

999   continue
      call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate sea salt arrays', rc)
      call error_mgr%pop_context()

   end subroutine seasalt_process_init

   !> \brief Execute sea salt emission calculations
   !!
   !! Performs the main sea salt emission calculation based on the selected scheme,
   !! updates emission and diagnostic states.
   !!
   !! \param[inout] this SeaSaltProcessType instance
   !! \param[inout] container StateContainer with all state objects
   !! \param[out] rc Return code
   !!
   subroutine seasalt_process_run(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_run', 'Executing sea salt emission calculations')

      ! Check if process is ready
      if (.not. this%is_ready()) then
         call error_mgr%report_error(ERROR_NOT_INITIALIZED, 'Sea salt process not ready', rc)
         return
      end if

      ! Reset diagnostic arrays
      this%wind_speed_10m = 0.0_fp
      this%whitecap_coverage = 0.0_fp
      this%total_seasalt_flux = 0.0_fp

      ! Execute scheme-specific calculations
      select case (this%seasalt_scheme)
      case (SEASALT_SCHEME_MONAHAN)
         call this%run_monahan_scheme(container, rc)
      case (SEASALT_SCHEME_GONG)
         call this%run_gong_scheme(container, rc)
      case (SEASALT_SCHEME_SSLT)
         call this%run_sslt_scheme(container, rc)
      case default
         call error_mgr%report_error(ERROR_INVALID_OPTION, 'Unknown sea salt scheme', rc)
         return
      end select

      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_COMPUTATION, 'Sea salt scheme calculation failed', rc)
         return
      end if

      ! Update diagnostic state with computed values
      call this%update_diagnostic_state(container, rc)
      if (rc /= CC_SUCCESS) then
         call error_mgr%report_error(ERROR_STATE_UPDATE, 'Failed to update diagnostic state', rc)
         return
      end if

      call error_mgr%pop_context()

   end subroutine seasalt_process_run

   !> \brief Finalize sea salt emission process
   !!
   !! Cleans up allocated memory and finalizes the sea salt process.
   !!
   !! \param[inout] this SeaSaltProcessType instance
   !! \param[out] rc Return code
   !!
   subroutine seasalt_process_finalize(this, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate arrays
      if (allocated(this%species_indices)) deallocate(this%species_indices)
      if (allocated(this%size_bin_centers)) deallocate(this%size_bin_centers)
      if (allocated(this%size_bin_widths)) deallocate(this%size_bin_widths)
      if (allocated(this%wind_speed_10m)) deallocate(this%wind_speed_10m)
      if (allocated(this%whitecap_coverage)) deallocate(this%whitecap_coverage)
      if (allocated(this%sea_surface_temp)) deallocate(this%sea_surface_temp)
      if (allocated(this%total_seasalt_flux)) deallocate(this%total_seasalt_flux)
      if (allocated(this%ocean_mask)) deallocate(this%ocean_mask)

      ! Reset state
      this%is_initialized = .false.
      this%is_active = .false.
      this%n_species = 0

   end subroutine seasalt_process_finalize

   !> \brief Set sea salt emission scheme
   subroutine seasalt_process_set_scheme(this, scheme)
      class(SeaSaltProcessType), intent(inout) :: this
      integer, intent(in) :: scheme
      this%seasalt_scheme = scheme
   end subroutine seasalt_process_set_scheme

   !> \brief Get sea salt emission scheme
   function seasalt_process_get_scheme(this) result(scheme)
      class(SeaSaltProcessType), intent(in) :: this
      integer :: scheme
      scheme = this%seasalt_scheme
   end function seasalt_process_get_scheme

   !> \brief Calculate sea salt emissions using selected scheme
   subroutine seasalt_process_calculate_emissions(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      ! This method delegates to the scheme-specific run method
      call this%run(container, rc)

   end subroutine seasalt_process_calculate_emissions

   !> \brief Calculate whitecap coverage
   !!
   !! Uses utilities from SeaSaltCommon_Mod to calculate whitecap coverage
   !! based on wind speed and optionally sea surface temperature.
   !!
   subroutine seasalt_process_whitecap(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(MetStateType), pointer :: met_state
      type(ErrorManagerType), pointer :: error_mgr
      real(fp), pointer :: wind_speed(:,:), sst(:,:)
      integer :: i, j

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_whitecap', 'Calculating whitecap coverage')

      met_state => container%get_met_state_ptr()
      if (.not. associated(met_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Met state not available', rc)
         return
      end if

      ! Get meteorological fields
      wind_speed => met_state%get_field_ptr('wind_speed_10m')
      sst => met_state%get_field_ptr('sea_surface_temperature')

      if (.not. associated(wind_speed)) then
         call error_mgr%report_error(ERROR_FIELD_NOT_FOUND, 'Wind speed not available', rc)
         return
      end if

      ! Calculate whitecap coverage for each grid cell
      do j = 1, this%n_lat
         do i = 1, this%n_lon

            ! Skip land points
            if (.not. this%ocean_mask(i,j)) cycle

            this%wind_speed_10m(i,j) = wind_speed(i,j)

            if (associated(sst) .and. this%use_sst_correction) then
               this%sea_surface_temp(i,j) = sst(i,j)
               call Whitecap_Coverage_Salisbury(wind_speed(i,j), sst(i,j), &
                                              this%whitecap_coverage(i,j), error_mgr, rc)
            else
               call Whitecap_Coverage_Monahan(wind_speed(i,j), &
                                            this%whitecap_coverage(i,j), error_mgr, rc)
            end if

            if (rc /= CC_SUCCESS) then
               call error_mgr%report_error(ERROR_COMPUTATION, 'Whitecap calculation failed', rc)
               return
            end if
         end do
      end do

      call error_mgr%pop_context()

   end subroutine seasalt_process_whitecap

   !> \brief Apply sea surface temperature correction
   subroutine seasalt_process_sst_correction(this, container, base_flux, corrected_flux, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      real(fp), intent(in) :: base_flux
      real(fp), intent(out) :: corrected_flux
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()

      if (.not. this%use_sst_correction) then
         corrected_flux = base_flux
         return
      end if

      ! Use average SST for correction (simplified)
      call Sea_Surface_Temperature_Correction(base_flux, 288.15_fp, corrected_flux, error_mgr, rc)

   end subroutine seasalt_process_sst_correction

   !> \brief Setup particle size bins for sea salt emissions
   subroutine seasalt_process_setup_size_bins(this, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      integer, intent(out) :: rc

      real(fp) :: log_min, log_max, log_width
      integer :: i

      rc = CC_SUCCESS

      if (allocated(this%size_bin_centers)) deallocate(this%size_bin_centers)
      if (allocated(this%size_bin_widths)) deallocate(this%size_bin_widths)

      allocate(this%size_bin_centers(this%n_size_bins), stat=rc)
      if (rc /= 0) return

      allocate(this%size_bin_widths(this%n_size_bins), stat=rc)
      if (rc /= 0) return

      ! Setup logarithmically-spaced size bins
      log_min = log10(this%min_particle_radius)
      log_max = log10(this%max_particle_radius)
      log_width = (log_max - log_min) / real(this%n_size_bins, fp)

      do i = 1, this%n_size_bins
         this%size_bin_centers(i) = 10.0_fp**(log_min + (real(i,fp) - 0.5_fp) * log_width)
         this%size_bin_widths(i) = 10.0_fp**(log_min + real(i,fp) * log_width) - &
                                  10.0_fp**(log_min + real(i-1,fp) * log_width)
      end do

   end subroutine seasalt_process_setup_size_bins

   !> \brief Map sea salt species names to chemical state indices
   subroutine seasalt_process_map_species(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ChemStateType), pointer :: chem_state
      type(ErrorManagerType), pointer :: error_mgr
      integer :: i, species_idx

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_map_species', 'Mapping sea salt species')

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

   end subroutine seasalt_process_map_species

   !> \brief Validate sea salt-specific configuration parameters
   subroutine seasalt_process_validate_config(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_validate_config', 'Validating sea salt configuration')

      ! Validate sea salt scheme
      if (this%seasalt_scheme < SEASALT_SCHEME_MONAHAN .or. this%seasalt_scheme > SEASALT_SCHEME_SSLT) then
         call error_mgr%report_error(ERROR_INVALID_OPTION, 'Invalid sea salt scheme', rc)
         return
      end if

      ! Validate physical parameters
      if (this%min_wind_speed < 0.0_fp .or. this%min_wind_speed > 20.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid minimum wind speed', rc)
         return
      end if

      if (this%max_wind_speed < this%min_wind_speed .or. this%max_wind_speed > 100.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid maximum wind speed', rc)
         return
      end if

      if (this%min_sst < 200.0_fp .or. this%min_sst > 350.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid minimum SST', rc)
         return
      end if

      if (this%min_particle_radius <= 0.0_fp .or. this%min_particle_radius >= this%max_particle_radius) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid particle radius range', rc)
         return
      end if

      if (this%n_size_bins < 1 .or. this%n_size_bins > 20) then
         call error_mgr%report_error(ERROR_INVALID_RANGE, 'Invalid number of size bins', rc)
         return
      end if

      call error_mgr%pop_context()

   end subroutine seasalt_process_validate_config

   !> \brief Run Monahan sea salt emission scheme
   !!
   !! Implements the Monahan et al. (1986) sea salt emission scheme
   !! using utilities from SeaSaltCommon_Mod.
   !!
   subroutine seasalt_process_run_monahan(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(EmisStateType), pointer :: emis_state
      type(ErrorManagerType), pointer :: error_mgr
      real(fp) :: seasalt_flux, total_flux
      integer :: i, j, k, ibin

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_run_monahan', 'Running Monahan sea salt scheme')

      ! Get emission state
      emis_state => container%get_emis_state_ptr()
      if (.not. associated(emis_state)) then
         call error_mgr%report_error(ERROR_NULL_POINTER, 'Emission state not available', rc)
         return
      end if

      ! Calculate whitecap coverage
      call this%calculate_whitecap_coverage(container, rc)
      if (rc /= CC_SUCCESS) return

      ! Calculate emissions for each grid cell
      do j = 1, this%n_lat
         do i = 1, this%n_lon

            ! Skip land points
            if (.not. this%ocean_mask(i,j)) cycle

            ! Skip if wind speed too low or too high
            if (this%wind_speed_10m(i,j) < this%min_wind_speed .or. &
                this%wind_speed_10m(i,j) > this%max_wind_speed) cycle

            ! Skip if SST too low
            if (this%sea_surface_temp(i,j) < this%min_sst) cycle

            total_flux = 0.0_fp

            ! Calculate emissions for each size bin
            do ibin = 1, this%n_size_bins
               call Monahan_SeaSaltFlux(this%size_bin_centers(ibin), &
                                       this%wind_speed_10m(i,j), &
                                       this%whitecap_coverage(i,j), &
                                       seasalt_flux, error_mgr, rc)
               if (rc /= CC_SUCCESS) then
                  call error_mgr%report_error(ERROR_COMPUTATION, 'Monahan flux calculation failed', rc)
                  return
               end if

               ! Convert from number flux to mass flux (simplified)
               seasalt_flux = seasalt_flux * this%size_bin_widths(ibin) * 1.0e-12_fp  ! Convert to kg/m²/s

               total_flux = total_flux + seasalt_flux

               ! Distribute among species (simplified - equal distribution)
               do k = 1, this%n_species
                  call emis_state%add_emission(this%species_indices(k), i, j, 1, &
                                             seasalt_flux / real(this%n_species, fp))
               end do
            end do

            this%total_seasalt_flux(i,j) = total_flux

         end do
      end do

      call error_mgr%pop_context()

   end subroutine seasalt_process_run_monahan

   !> \brief Run Gong sea salt emission scheme
   !!
   !! Placeholder for Gong sea salt emission scheme implementation.
   !!
   subroutine seasalt_process_run_gong(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_run_gong', 'Running Gong sea salt scheme')

      ! TODO: Implement Gong scheme using Gong_SeaSaltFlux
      call error_mgr%report_info('Gong scheme not yet fully implemented')

      ! For now, use Monahan as placeholder
      call this%run_monahan_scheme(container, rc)

      call error_mgr%pop_context()

   end subroutine seasalt_process_run_gong

   !> \brief Run SSLT sea salt emission scheme
   !!
   !! Placeholder for SSLT sea salt emission scheme implementation.
   !!
   subroutine seasalt_process_run_sslt(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
      type(StateContainerType), intent(inout) :: container
      integer, intent(out) :: rc

      type(ErrorManagerType), pointer :: error_mgr

      rc = CC_SUCCESS
      error_mgr => container%get_error_manager()
      call error_mgr%push_context('seasalt_process_run_sslt', 'Running SSLT sea salt scheme')

      ! TODO: Implement SSLT scheme
      call error_mgr%report_info('SSLT scheme not yet implemented')

      call error_mgr%pop_context()

   end subroutine seasalt_process_run_sslt

   !> \brief Update diagnostic state with computed sea salt values
   subroutine update_diagnostic_state(this, container, rc)
      class(SeaSaltProcessType), intent(inout) :: this
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
      call diag_state%set_field('seasalt_wind_speed_10m', this%wind_speed_10m)
      call diag_state%set_field('seasalt_whitecap_coverage', this%whitecap_coverage)
      call diag_state%set_field('seasalt_sea_surface_temp', this%sea_surface_temp)
      call diag_state%set_field('seasalt_total_flux', this%total_seasalt_flux)

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

end module SeaSaltProcess_Mod
