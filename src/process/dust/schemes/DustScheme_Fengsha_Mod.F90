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
!! Generated on: 2025-07-09T12:43:17.999787
!! Author: Barry Baker
!! Reference: Zhang et al. 2022
module DustScheme_FENGSHA_Mod

   use iso_fortran_env, only: fp => real64
   use constants, only: g0

   implicit none
   private

   ! Public interface - pure science only
   public :: compute_fengsha
   public :: fengsha_params_t


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
   !! @param[in]  IsLand    IsLand field [appropriate units]
   !! @param[in]  USTAR    USTAR field [appropriate units]
   !! @param[in]  LWI    LWI field [appropriate units]
   !! @param[in]  GVF    GVF field [appropriate units]
   !! @param[in]  LAI    LAI field [appropriate units]
   !! @param[in]  FROCEAN    FROCEAN field [appropriate units]
   !! @param[in]  CLAYFRAC    CLAYFRAC field [appropriate units]
   !! @param[in]  SANDFRAC    SANDFRAC field [appropriate units]
   !! @param[in]  FRSNO    FRSNO field [appropriate units]
   !! @param[in]  RDRAG    RDRAG field [appropriate units]
   !! @param[in]  SSM    SSM field [appropriate units]
   !! @param[in]  USTAR_THRESHOLD    USTAR_THRESHOLD field [appropriate units]
   !! @param[in]  species_conc   Species concentrations [mol/mol] (num_layers, num_species)
   !! @param[out] emission_flux  Emission fluxes [kg/m²/s] (num_layers, num_species)
   pure subroutine compute_fengsha(num_layers, num_species, params, IsLand, USTAR, &
      LWI, GVF, LAI, FROCEAN, CLAYFRAC, SANDFRAC, FRSNO, RDRAG, SSM, USTAR_THRESHOLD, &
      species_conc, emission_flux, horizontal_flux, vertical_flux, effective_threshold)
   )

      ! Arguments
      integer, intent(in) :: num_layers
      integer, intent(in) :: num_species
      type(fengsha_params_t), intent(in) :: params
      real(fp), intent(in) :: IsLand(num_layers)
      real(fp), intent(in) :: USTAR(num_layers)
      real(fp), intent(in) :: LWI(num_layers)
      real(fp), intent(in) :: GVF(num_layers)
      real(fp), intent(in) :: LAI(num_layers)
      real(fp), intent(in) :: FROCEAN(num_layers)
      real(fp), intent(in) :: CLAYFRAC(num_layers)
      real(fp), intent(in) :: SANDFRAC(num_layers)
      real(fp), intent(in) :: FRSNO(num_layers)
      real(fp), intent(in) :: RDRAG(num_layers)
      real(fp), intent(in) :: SSM(num_layers)
      real(fp), intent(in) :: USTAR_THRESHOLD(num_layers)
      real(fp), intent(in) :: species_conc(num_layers, num_species)
      real(fp), intent(out) :: emission_flux(num_layers, num_species)
      real(fp), intent(out) :: horizontal_flux(num_layers, num_species)
      real(fp), intent(out) :: vertical_flux(num_layers, num_species)
      real(fp), intent(out) :: effective_threshold(num_layers, num_species)

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
         ! Consider how IsLand affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how USTAR affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how LWI affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how GVF affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how LAI affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how FROCEAN affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how CLAYFRAC affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how SANDFRAC affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how FRSNO affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how RDRAG affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how SSM affects your emissions
         ! Generic field usage (customize for your scheme)
         ! Consider how USTAR_THRESHOLD affects your emissions

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