!> \file DustProcessCreator_Mod.F90
!! \brief Factory for creating dust process instances
!!
!! This module provides the factory functions for creating dust
!! process instances following the CATChem Process Factory pattern.
!!
!! Generated on: 2025-07-09T12:43:17.993393
!! Author: Barry Baker
!! Version: 1.0.0

module DustProcessCreator_Mod

   use iso_fortran_env, only: fp => real64
   use ProcessInterface_Mod
   use ProcessDustInterface_Mod
   use ErrorHandler_Mod

   implicit none
   private

   public :: create_dust_process
   public :: get_dust_process_info

   ! Process metadata constants
   character(len=*), parameter :: PROCESS_NAME = 'dust'
   character(len=*), parameter :: PROCESS_VERSION = '1.0.0'
   character(len=*), parameter :: PROCESS_DESCRIPTION = 'Process for computing windblown dust emissions'
   character(len=*), parameter :: PROCESS_AUTHOR = 'Barry Baker'

contains

   !> Create a new dust process instance
   !!
   !! This factory function creates and returns a new instance of the
   !! dust process. The process is not initialized - the caller
   !! must call the init() method with appropriate configuration.
   !!
   !! @param[out] process     Allocated process instance
   !! @param[inout] error_handler Error handling object
   function create_dust_process(error_handler) result(process)
      type(ErrorHandler), intent(inout) :: error_handler
      class(ProcessInterface), allocatable :: process

      type(ProcessDustInterface), allocatable :: dust_process
      integer :: alloc_stat

      ! Allocate the process instance
      allocate(dust_process, stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_handler%set_error(ERROR_MEMORY, &
            "Failed to allocate dust process instance")
         return
      end if

      ! Move to polymorphic variable
      call move_alloc(dust_process, process)

   end function create_dust_process

   !> Get information about the dust process
   !!
   !! This function returns metadata about the dust process
   !! including name, version, description, and capabilities.
   !!
   !! @param[out] info Process information structure
   subroutine get_dust_process_info(info)
      type(ProcessInfo), intent(out) :: info

      ! Basic information
      info%name = PROCESS_NAME
      info%version = PROCESS_VERSION
      info%description = PROCESS_DESCRIPTION
      info%author = PROCESS_AUTHOR

      ! Process characteristics
      info%process_type = 'emission'
      info%is_multiphase = .false.
      info%has_size_bins = .false.
      info%supports_vectorization = .true.

      ! Integration characteristics
      info%timestep_dependency = 'independent'
      info%parallelization = 'column'
      info%memory_requirements = 'low'

      ! Species information
      info%n_species = 0



      ! Available schemes
      info%n_schemes = 2
      allocate(info%scheme_names(info%n_schemes))
      allocate(info%scheme_descriptions(info%n_schemes))
      info%scheme_names(1) = 'fengsha'
      info%scheme_descriptions(1) = 'Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS'
      info%scheme_names(2) = 'ginoux'
      info%scheme_descriptions(2) = 'Ginoux dust emission scheme'
      info%default_scheme = ''

      ! Required meteorological fields
      info%n_required_met_fields = 3
      allocate(info%required_met_fields(info%n_required_met_fields))
      info%required_met_fields(1) = 'ustar'
      info%required_met_fields(2) = 'solar_zenith_angle'
      info%required_met_fields(3) = 'leaf_area_index'

      ! Optional meteorological fields
      info%n_optional_met_fields = 0

      ! Required chemical fields
      info%n_required_chem_fields = 0

      ! Diagnostic information
      info%n_diagnostics = 1
      allocate(info%diagnostic_names(info%n_diagnostics))
      allocate(info%diagnostic_descriptions(info%n_diagnostics))
      allocate(info%diagnostic_units(info%n_diagnostics))
      info%diagnostic_names(1) = 'total_dust_emission'
      info%diagnostic_descriptions(1) = 'Total dust emissions for all species'
      info%diagnostic_units(1) = 'kg/m2/s'

   end subroutine get_dust_process_info

   !> Check if dust process is compatible with host model
   !!
   !! This function checks whether the dust process is compatible
   !! with the host model configuration and requirements.
   !!
   !! @param[in] host_config Host model configuration
   !! @param[inout] error_handler Error handling object
   !! @return .true. if compatible, .false. otherwise
   function is_dust_compatible(host_config, error_handler) result(is_compatible)
      type(HostModelConfig), intent(in) :: host_config
      type(ErrorHandler), intent(inout) :: error_handler
      logical :: is_compatible

      is_compatible = .true.

      ! Check required species

      ! Check required meteorological fields
      if (.not. host_config%has_met_field('ustar')) then
         call error_handler%set_warning(WARNING_COMPATIBILITY, &
            "Required meteorological field not available: ustar")
         is_compatible = .false.
      end if
      if (.not. host_config%has_met_field('solar_zenith_angle')) then
         call error_handler%set_warning(WARNING_COMPATIBILITY, &
            "Required meteorological field not available: solar_zenith_angle")
         is_compatible = .false.
      end if
      if (.not. host_config%has_met_field('leaf_area_index')) then
         call error_handler%set_warning(WARNING_COMPATIBILITY, &
            "Required meteorological field not available: leaf_area_index")
         is_compatible = .false.
      end if



   end function is_dust_compatible

   !> Validate dust process configuration
   !!
   !! This function validates the process configuration against the
   !! process requirements and constraints.
   !!
   !! @param[in] config_data Configuration data (YAML, namelist, etc.)
   !! @param[inout] error_handler Error handling object
   !! @return .true. if valid, .false. otherwise
   function validate_dust_config_data(config_data, error_handler) result(is_valid)
      character(len=*), intent(in) :: config_data
      type(ErrorHandler), intent(inout) :: error_handler
      logical :: is_valid

      ! TODO: Implement configuration validation
      ! This should parse the config_data and validate:
      ! - Scheme selection is valid
      ! - Parameters are within acceptable ranges
      ! - Required fields are present
      ! - No conflicting options

      is_valid = .true.

      ! Placeholder validation
      if (len_trim(config_data) == 0) then
         call error_handler%set_error(ERROR_CONFIG, &
            "Empty configuration data provided")
         is_valid = .false.
      end if

   end function validate_dust_config_data

   !> Get default configuration for dust process
   !!
   !! This function returns a default configuration string that can be
   !! used to initialize the dust process with reasonable defaults.
   !!
   !! @param[out] config_data Default configuration string
   subroutine get_dust_default_config(config_data)
      character(len=*), intent(out) :: config_data

      ! Return default YAML configuration
      config_data = &
         '# Default dust process configuration' // new_line('A') // &
         'process:' // new_line('A') // &
         '  name: "dust"' // new_line('A') // &
         '  version: "1.0.0"' // new_line('A') // &
         '  active_scheme: ""' // new_line('A') // &
         '  is_active: true' // new_line('A') // &
         '' // new_line('A') // &
         '# Scheme configuration' // new_line('A') // &
         'schemes:' // new_line('A') // &
         '  fengsha:' // new_line('A') // &
         '    description: "Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      alpha: 0.16' // new_line('A') // &
         '      beta: 1.0' // new_line('A') // &
         '      drylimit_factor: 1.0' // new_line('A') // &
         '      drag_option: 1' // new_line('A') // &
         '      moist_option: 1 - fecan' // new_line('A') // &
         '      distribution_option: 1' // new_line('A') // &
         '' // new_line('A') // &
         '  ginoux:' // new_line('A') // &
         '    description: "Ginoux dust emission scheme"' // new_line('A') // &
         '    algorithm_type: "explicit"' // new_line('A') // &
         '    parameters:' // new_line('A') // &
         '      Ch_DU: [0.1, 0.1, 0.1, 0.1, 0.1]' // new_line('A') // &
         '' // new_line('A') // &
         '# Diagnostic configuration' // new_line('A') // &
         'diagnostics:' // new_line('A') // &
         '  output_frequency: 3600.0  # seconds' // new_line('A') // &
         '  output_diagnostics: true'

   end subroutine get_dust_default_config

end module DustProcessCreator_Mod