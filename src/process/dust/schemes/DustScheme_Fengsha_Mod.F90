!> \file DustScheme_FENGSHA_Mod.F90
!! \brief Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS
!!
!! Pure science kernel for fengsha scheme in dust process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_fengsha (search for "TODO")
!! 2. Add scheme-specific helper subroutines as needed
!! 3. Update physical constants for your scheme
!! 4. Customize the environmental response functions
!!
!! INFRASTRUCTURE RESPONSIBILITIES (handled by host model):
!! - Parameter initialization and validation
!! - Input array validation and error handling
!! - Memory management and array allocation
!! - Integration with host model time stepping
!!
!! Generated on: 2025-08-03T14:41:50.766057
!! Author: Barry Baker
!! Reference: Zhang et al. 2022
module DustScheme_FENGSHA_Mod

   use iso_fortran_env, only: fp => real64

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_fengsha
   public :: fengsha_params_t

   ! Physical constants (modify as needed for your scheme)
   real(fp), parameter :: R_GAS = 8.314_fp           ! Universal gas constant [J/mol/K]
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor


   !> Science parameters for fengsha scheme
   !! Host model is responsible for initializing and validating these
   type :: fengsha_params_t
      real(fp) :: alpha  ! linear scaling factor
      real(fp) :: beta  ! Exponential scaling factor on source parameter
      real(fp) :: drylimit_factor  ! Dry Limit factor modifying the Fecan dry limit following Zender 2003
      real(fp) :: drag_option  ! Drag Partition Option: 1 - use input drag, 2 - Darmenova, 3 - Leung 2022, 4 - MB95
      real(fp) :: moist_option  ! Moisture parameterization: 1 - Fecan, 2 - shao, 3 - modified shao
      real(fp) :: distribution_option  ! Dust Distribution option: 1 - Kok 2011, 2 - Meng 2022
   end type fengsha_params_t


contains

   !> Pure science computation for fengsha scheme
   !!
   !! This is a pure computational kernel implementing Fengsha Dust emission scheme developed at NOAA ARL for use at NOAA NWS.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  {'name': 'IsLand', 'description': 'Land mask', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'is_land'}    {'name': 'IsLand', 'description': 'Land mask', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'is_land'} field [appropriate units]
   !! @param[in]  {'name': 'USTAR', 'description': 'Friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'friction_velocity'}    {'name': 'USTAR', 'description': 'Friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'friction_velocity'} field [appropriate units]
   !! @param[in]  {'name': 'LWI', 'description': 'Land-water index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'land_water_index'}    {'name': 'LWI', 'description': 'Land-water index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'land_water_index'} field [appropriate units]
   !! @param[in]  {'name': 'GVF', 'description': 'Green vegetation fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'green_veg_fraction'}    {'name': 'GVF', 'description': 'Green vegetation fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'green_veg_fraction'} field [appropriate units]
   !! @param[in]  {'name': 'LAI', 'description': 'Leaf area index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'leaf_area_index'}    {'name': 'LAI', 'description': 'Leaf area index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'leaf_area_index'} field [appropriate units]
   !! @param[in]  {'name': 'FROCEAN', 'description': 'Ocean fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'ocean_fraction'}    {'name': 'FROCEAN', 'description': 'Ocean fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'ocean_fraction'} field [appropriate units]
   !! @param[in]  {'name': 'CLAYFRAC', 'description': 'Clay fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'clay_fraction'}    {'name': 'CLAYFRAC', 'description': 'Clay fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'clay_fraction'} field [appropriate units]
   !! @param[in]  {'name': 'SANDFRAC', 'description': 'Sand fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'sand_fraction'}    {'name': 'SANDFRAC', 'description': 'Sand fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'sand_fraction'} field [appropriate units]
   !! @param[in]  {'name': 'FRSNO', 'description': 'Snow fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'snow_fraction'}    {'name': 'FRSNO', 'description': 'Snow fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'snow_fraction'} field [appropriate units]
   !! @param[in]  {'name': 'RDRAG', 'description': 'Drag partition', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'drag_partition'}    {'name': 'RDRAG', 'description': 'Drag partition', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'drag_partition'} field [appropriate units]
   !! @param[in]  {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'}    {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'} field [appropriate units]
   !! @param[in]  {'name': 'USTAR_THRESHOLD', 'description': 'Threshold friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'threshold_ustar'}    {'name': 'USTAR_THRESHOLD', 'description': 'Threshold friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'threshold_ustar'} field [appropriate units]
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[out] emission_flux  Emission fluxes [kg/m²/s] (num_layers, num_species)
   pure subroutine compute_fengsha( &
      num_layers, &
      num_species, &
      params, &
      {'name': 'IsLand', 'description': 'Land mask', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'is_land'}, &      {'name': 'USTAR', 'description': 'Friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'friction_velocity'}, &      {'name': 'LWI', 'description': 'Land-water index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'land_water_index'}, &      {'name': 'GVF', 'description': 'Green vegetation fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'green_veg_fraction'}, &      {'name': 'LAI', 'description': 'Leaf area index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'leaf_area_index'}, &      {'name': 'FROCEAN', 'description': 'Ocean fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'ocean_fraction'}, &      {'name': 'CLAYFRAC', 'description': 'Clay fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'clay_fraction'}, &      {'name': 'SANDFRAC', 'description': 'Sand fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'sand_fraction'}, &      {'name': 'FRSNO', 'description': 'Snow fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'snow_fraction'}, &      {'name': 'RDRAG', 'description': 'Drag partition', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'drag_partition'}, &      {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'}, &      {'name': 'USTAR_THRESHOLD', 'description': 'Threshold friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'threshold_ustar'}, &
      species_conc, &
      emission_flux &
   )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(fengsha_params_t), intent(in) :: params
      real(fp), intent(in) :: {'name': 'IsLand', 'description': 'Land mask', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'is_land'}(num_layers)
      real(fp), intent(in) :: {'name': 'USTAR', 'description': 'Friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'friction_velocity'}(num_layers)
      real(fp), intent(in) :: {'name': 'LWI', 'description': 'Land-water index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'land_water_index'}(num_layers)
      real(fp), intent(in) :: {'name': 'GVF', 'description': 'Green vegetation fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'green_veg_fraction'}(num_layers)
      real(fp), intent(in) :: {'name': 'LAI', 'description': 'Leaf area index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'leaf_area_index'}(num_layers)
      real(fp), intent(in) :: {'name': 'FROCEAN', 'description': 'Ocean fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'ocean_fraction'}(num_layers)
      real(fp), intent(in) :: {'name': 'CLAYFRAC', 'description': 'Clay fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'clay_fraction'}(num_layers)
      real(fp), intent(in) :: {'name': 'SANDFRAC', 'description': 'Sand fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'sand_fraction'}(num_layers)
      real(fp), intent(in) :: {'name': 'FRSNO', 'description': 'Snow fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'snow_fraction'}(num_layers)
      real(fp), intent(in) :: {'name': 'RDRAG', 'description': 'Drag partition', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'drag_partition'}(num_layers)
      real(fp), intent(in) :: {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'}(num_layers)
      real(fp), intent(in) :: {'name': 'USTAR_THRESHOLD', 'description': 'Threshold friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'threshold_ustar'}(num_layers)
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(out) :: emission_flux(num_layers, num_species)

      ! Local variables
      integer :: k, species_idx
      real(fp) :: base_emission_factor
      real(fp) :: temperature_factor
      real(fp) :: light_factor
      real(fp) :: combined_factor

      ! Initialize output (pure subroutines must initialize all outputs)
      emission_flux = 0.0_fp

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! TODO: Replace this generic implementation with your scheme's algorithm
         ! This is a placeholder that demonstrates the expected structure

         ! Initialize environmental factors
         temperature_factor = 1.0_fp
         light_factor = 1.0_fp

         ! Apply scheme-specific environmental responses
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'IsLand', 'description': 'Land mask', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'is_land'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'USTAR', 'description': 'Friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'friction_velocity'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'LWI', 'description': 'Land-water index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'land_water_index'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'GVF', 'description': 'Green vegetation fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'green_veg_fraction'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'LAI', 'description': 'Leaf area index', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'leaf_area_index'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'FROCEAN', 'description': 'Ocean fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'ocean_fraction'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'CLAYFRAC', 'description': 'Clay fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'clay_fraction'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'SANDFRAC', 'description': 'Sand fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'sand_fraction'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'FRSNO', 'description': 'Snow fraction', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'snow_fraction'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'RDRAG', 'description': 'Drag partition', 'units': 'dimensionless', 'dimensions': 'scalar', 'variable_name': 'drag_partition'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'SSM', 'description': 'Surface soil moisture', 'units': 'm3/m3', 'dimensions': 'scalar', 'variable_name': 'soil_moisture'} affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how {'name': 'USTAR_THRESHOLD', 'description': 'Threshold friction velocity', 'units': 'm/s', 'dimensions': 'scalar', 'variable_name': 'threshold_ustar'} affects your emissions

         ! Combine environmental factors
         combined_factor = temperature_factor * light_factor

         ! Apply to each species
         do species_idx = 1, num_species
            ! Base emission factor (customize this for species-specific emissions)
            base_emission_factor = DEFAULT_SCALING

            ! Compute emission flux using your scheme's formula
            ! This is a simple example - replace with your actual algorithm
            emission_flux(k, species_idx) = base_emission_factor * combined_factor * &
                                          species_conc(k, species_idx)

            ! Ensure non-negative emissions
            emission_flux(k, species_idx) = max(0.0_fp, emission_flux(k, species_idx))
         end do

      end do

   end subroutine compute_fengsha

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   ! Example helper subroutine template:
   !
   ! !> Compute custom environmental response for fengsha scheme
   ! pure function compute_temperature_response_fengsha(temperature, params) result(factor)
   !    real(fp), intent(in) :: temperature     ! Temperature [K]
   !    type(fengsha_params_t), intent(in) :: params
   !    real(fp) :: factor
   !
   !    ! Add your temperature response algorithm here
   !    factor = exp(params%temperature_dependency * (temperature - T_STANDARD) / T_STANDARD)
   ! end function compute_temperature_response_fengsha

end module DustScheme_FENGSHA_Mod