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
!! Generated on: 2025-08-05T10:07:05.820923
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
   real(fp), parameter :: PI = 3.14159265359_fp     ! Pi constant

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
   !! @param[in]  is_land    Land mask [dimensionless]
   !! @param[in]  friction_velocity    Friction velocity [m/s]
   !! @param[in]  land_water_index    Land-water index [dimensionless]
   !! @param[in]  green_veg_fraction    Green vegetation fraction [dimensionless]
   !! @param[in]  leaf_area_index    Leaf area index [dimensionless]
   !! @param[in]  ocean_fraction    Ocean fraction [dimensionless]
   !! @param[in]  clay_fraction    Clay fraction [dimensionless]
   !! @param[in]  sand_fraction    Sand fraction [dimensionless]
   !! @param[in]  snow_fraction    Snow fraction [dimensionless]
   !! @param[in]  drag_partition    Drag partition [dimensionless]
   !! @param[in]  soil_moisture    Surface soil moisture [m3/m3]
   !! @param[in]  threshold_ustar    Threshold friction velocity [m/s]
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[out] emission_flux  Emission fluxes [kg/m²/s] (num_layers, num_species)
   pure subroutine compute_fengsha( &
      num_layers, &
      num_species, &
      params, &
      is_land, &      friction_velocity, &      land_water_index, &      green_veg_fraction, &      leaf_area_index, &      ocean_fraction, &      clay_fraction, &      sand_fraction, &      snow_fraction, &      drag_partition, &      soil_moisture, &      threshold_ustar, &
      species_conc, &
      emission_flux &
   )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(fengsha_params_t), intent(in) :: params
      real(fp), intent(in) :: is_land  ! Land mask [dimensionless]
      real(fp), intent(in) :: friction_velocity  ! Friction velocity [m/s]
      real(fp), intent(in) :: land_water_index  ! Land-water index [dimensionless]
      real(fp), intent(in) :: green_veg_fraction  ! Green vegetation fraction [dimensionless]
      real(fp), intent(in) :: leaf_area_index  ! Leaf area index [dimensionless]
      real(fp), intent(in) :: ocean_fraction  ! Ocean fraction [dimensionless]
      real(fp), intent(in) :: clay_fraction  ! Clay fraction [dimensionless]
      real(fp), intent(in) :: sand_fraction  ! Sand fraction [dimensionless]
      real(fp), intent(in) :: snow_fraction  ! Snow fraction [dimensionless]
      real(fp), intent(in) :: drag_partition  ! Drag partition [dimensionless]
      real(fp), intent(in) :: soil_moisture  ! Surface soil moisture [m3/m3]
      real(fp), intent(in) :: threshold_ustar  ! Threshold friction velocity [m/s]
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(out) :: emission_flux(num_layers, num_species)

      ! Local variables
      integer :: k, species_idx
      real(fp) :: base_emission_factor
      real(fp) :: environmental_factor
      real(fp) :: species_factor

      ! Initialize output (pure subroutines must initialize all outputs)
      emission_flux = 0.0_fp

      ! Main computation loop - CUSTOMIZE THIS SECTION FOR YOUR SCHEME
      do k = 1, num_layers

         ! TODO: Replace this generic implementation with your scheme's algorithm
         ! This is a placeholder that demonstrates the expected structure

         ! Initialize environmental factors
         environmental_factor = 1.0_fp

         ! Apply scheme-specific environmental responses based on meteorological fields
         ! Land mask usage (customize for your scheme)
         ! TODO: Consider how IsLand affects your emissions
         ! Field available: is_land [dimensionless]
         ! Wind/friction velocity response (customize for your scheme)
         ! TODO: Implement wind-dependent emission formula
         environmental_factor = environmental_factor * &
                               max(0.0_fp, friction_velocity)
         ! Land-water index usage (customize for your scheme)
         ! TODO: Consider how LWI affects your emissions
         ! Field available: land_water_index [dimensionless]
         ! Green vegetation fraction usage (customize for your scheme)
         ! TODO: Consider how GVF affects your emissions
         ! Field available: green_veg_fraction [dimensionless]
         ! Vegetation response (customize for your scheme)
         environmental_factor = environmental_factor * &
                               (1.0_fp - min(leaf_area_index, 5.0_fp) / 5.0_fp)
         ! Ocean fraction usage (customize for your scheme)
         ! TODO: Consider how FROCEAN affects your emissions
         ! Field available: ocean_fraction [dimensionless]
         ! Clay fraction usage (customize for your scheme)
         ! TODO: Consider how CLAYFRAC affects your emissions
         ! Field available: clay_fraction [dimensionless]
         ! Sand fraction usage (customize for your scheme)
         ! TODO: Consider how SANDFRAC affects your emissions
         ! Field available: sand_fraction [dimensionless]
         ! Snow fraction usage (customize for your scheme)
         ! TODO: Consider how FRSNO affects your emissions
         ! Field available: snow_fraction [dimensionless]
         ! Drag partition usage (customize for your scheme)
         ! TODO: Consider how RDRAG affects your emissions
         ! Field available: drag_partition [dimensionless]
         ! Soil moisture response (customize for your scheme)
         ! TODO: Implement moisture inhibition
         environmental_factor = environmental_factor * &
                               max(0.1_fp, 1.0_fp - soil_moisture)
         ! Wind/friction velocity response (customize for your scheme)
         ! TODO: Implement wind-dependent emission formula
         environmental_factor = environmental_factor * &
                               max(0.0_fp, threshold_ustar)

         ! Apply to each species
         do species_idx = 1, num_species
            ! Base emission factor (customize this for species-specific emissions)
            base_emission_factor = DEFAULT_SCALING

            ! Species-specific factor (customize based on species properties)
            species_factor = 1.0_fp  ! TODO: Add species-specific scaling

            ! Compute emission flux using your scheme's formula
            ! This is a simple example - replace with your actual algorithm
            emission_flux(k, species_idx) = base_emission_factor * &
                                          environmental_factor * &
                                          species_factor * &
                                          (1.0_fp + species_conc(k, species_idx))

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

   !> Example helper function for environmental response
   pure function compute_fengsha(met_value, reference_value) result(factor)
      real(fp), intent(in) :: met_value       ! Meteorological value
      real(fp), intent(in) :: reference_value ! Reference value
      real(fp) :: factor

      ! Simple exponential response - customize for your scheme
      factor = exp((met_value - reference_value) / reference_value)
      factor = max(0.0_fp, min(10.0_fp, factor))  ! Reasonable bounds
   end function compute_fengsha

   !> Example helper function for species-specific scaling
   pure function compute_species_scaling_fengsha(species_idx, params) result(scaling)
      integer, intent(in) :: species_idx
      type(fengsha_params_t), intent(in) :: params
      real(fp) :: scaling

      ! Species-specific scaling - customize for your scheme
      select case (species_idx)
      case (1)
         scaling = 1.0_fp    ! First species baseline
      case (2:3)
         scaling = 0.5_fp    ! Reduced emission for species 2-3
      case default
         scaling = 0.1_fp    ! Low emission for other species
      end select

   end function compute_species_scaling_fengsha

end module DustScheme_FENGSHA_Mod