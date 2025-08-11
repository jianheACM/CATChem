!> \file SeaSaltScheme_GONG97_Mod.F90
!! \brief Gong 1997 sea salt emission scheme
!!
!! Pure science kernel for gong97 scheme in seasalt process.
!! This module contains ONLY the computational algorithm with NO infrastructure dependencies.
!! Uses only basic Fortran types for maximum portability and reusability.
!!
!! SCIENCE CUSTOMIZATION GUIDE:
!! 1. Modify the algorithm in compute_gong97 (search for "TODO")
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
!! Generated on: 2025-08-11T15:55:19.557843
!! Author: Barry Baker
!! Reference: Gong et al. [1997]
module SeaSaltScheme_GONG97_Mod

   use iso_fortran_env, only: fp => real64

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_gong97
   public :: gong97_params_t

   ! Physical constants (modify as needed for your scheme)
   real(fp), parameter :: R_GAS = 8.314_fp           ! Universal gas constant [J/mol/K]
   real(fp), parameter :: T_STANDARD = 303.15_fp    ! Standard reference temperature [K]
   real(fp), parameter :: DEFAULT_SCALING = 1.0e-9_fp ! Default emission scaling factor
   real(fp), parameter :: PI = 3.14159265359_fp     ! Pi constant

   !> Science parameters for gong97 scheme
   !! Host model is responsible for initializing and validating these
   type :: gong97_params_t
      real(fp) :: scale_factor  ! Emission scale factor
      real(fp) :: weibull_flag  ! Apply Weibull distribution for particle size
   end type gong97_params_t

contains

   !> Pure science computation for gong97 scheme
   !!
   !! This is a pure computational kernel implementing Gong 1997 sea salt emission scheme.
   !! NO error checking, validation, or infrastructure concerns.
   !! Host model must ensure all inputs are valid before calling.
   !!
   !! @param[in]  num_layers     Number of vertical layers
   !! @param[in]  num_species    Number of chemical species
   !! @param[in]  params         Scheme parameters (pre-validated by host)
   !! @param[in]  frocean    FROCEAN field [appropriate units]
   !! @param[in]  frseaice    FRSEAICE field [appropriate units]
   !! @param[in]  sst    SST field [appropriate units]
   !! @param[in]  u10m    U10M field [appropriate units]
   !! @param[in]  v10m    V10M field [appropriate units]
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[out] emission_flux  Emission fluxes [kg/m²/s] (num_layers, num_species)
   pure subroutine compute_gong97( &
      num_layers, &
      num_species, &
      params, &
      frocean, &      frseaice, &      sst, &      u10m, &      v10m, &
      species_conc, &
      emission_flux &
   )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(gong97_params_t), intent(in) :: params
      real(fp), intent(in) :: frocean(num_layers)
      real(fp), intent(in) :: frseaice(num_layers)
      real(fp), intent(in) :: sst(num_layers)
      real(fp), intent(in) :: u10m(num_layers)
      real(fp), intent(in) :: v10m(num_layers)
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
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how FROCEAN affects your emissions
         ! environmental_factor = environmental_factor * some_function(frocean(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how FRSEAICE affects your emissions
         ! environmental_factor = environmental_factor * some_function(frseaice(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how SST affects your emissions
         ! environmental_factor = environmental_factor * some_function(sst(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how U10M affects your emissions
         ! environmental_factor = environmental_factor * some_function(u10m(k))
         ! Generic field usage (customize for your scheme)
         ! TODO: Consider how V10M affects your emissions
         ! environmental_factor = environmental_factor * some_function(v10m(k))

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

   end subroutine compute_gong97

   ! =======================================================================
   ! SCHEME-SPECIFIC HELPER SUBROUTINES
   ! =======================================================================
   ! Add your custom scientific algorithms here as pure functions/subroutines
   ! Examples: environmental response functions, species-specific calculations, etc.

   !> Example helper function for environmental response
   pure function compute_gong97(met_value, reference_value) result(factor)
      real(fp), intent(in) :: met_value       ! Meteorological value
      real(fp), intent(in) :: reference_value ! Reference value
      real(fp) :: factor

      ! Simple exponential response - customize for your scheme
      factor = exp((met_value - reference_value) / reference_value)
      factor = max(0.0_fp, min(10.0_fp, factor))  ! Reasonable bounds
   end function compute_gong97

   !> Example helper function for species-specific scaling
   pure function compute_species_scaling_gong97(species_idx, params) result(scaling)
      integer, intent(in) :: species_idx
      type(gong97_params_t), intent(in) :: params
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

   end function compute_species_scaling_gong97

end module SeaSaltScheme_GONG97_Mod