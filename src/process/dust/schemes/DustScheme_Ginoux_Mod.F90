!> \file DustScheme_Ginoux_Mod.F90
!! \brief Modern Ginoux dust emission scheme implementation - StateContainer independent
!! \ingroup dust_schemes
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! \details
!! The Ginoux scheme calculates dust emissions based on surface wind speed,
!! soil moisture, and surface characteristics. It provides size-resolved
!! dust emission fluxes for integration with atmospheric transport models.
!!
!! \section ginoux_reference Reference
!! Ginoux, P., "Sources and distributions of dust aerosols simulated with the GOCART model",
!! Journal of Geophysical Research, vol. 106, no. D17, pp. 20,255–20,273, 2001.
!! doi:10.1029/2000JD000053.

!!
!! \author Barry baker
!! \date 05/2024
!! \ingroup catchem_dust_process
module DustScheme_Ginoux_Mod
   use precision_mod, only: fp, ZERO
   use constants, only: g0
   use DustCommon_Mod
   implicit none
   private

   public :: run_ginoux_scheme

contains

   !> \brief Execute Ginoux dust emission scheme (single-column)
   !!
   !! This is a placeholder implementation for the Ginoux dust emission scheme
   !! for a single column. Currently returns zero emissions.
   !!
   subroutine run_ginoux_scheme(n_dust_species, &
                                ! Meteorological inputs (scalars)
                                wind_friction_velocity, &
                                wind_speed_10m, &
                                surface_temperature, &
                                air_density, &
                                soil_moisture_top, &
                                surface_roughness, &
                                ! Surface properties (scalars)
                                clay_fraction, &
                                sand_fraction, &
                                soil_type, &
                                ocean_fraction, &
                                land_ice_fraction, &
                                snow_fraction, &
                                sediment_supply_map, &
                                roughness_drag, &
                                ! Dust size distribution (1D arrays)
                                effective_radius, &
                                lower_radius, &
                                upper_radius, &
                                ! Configuration parameters
                                alpha_scale, &
                                beta_scale, &
                                moist_opt, &
                                drag_opt, &
                                horiz_flux_opt, &
                                ! Outputs (1D array for size bins)
                                total_emission, &
                                emission_per_species, &
                                rc)
      implicit none

      ! Inputs
      integer, intent(in) :: n_dust_species

      ! Meteorological fields (scalars - single column)
      real(fp), intent(in) :: wind_friction_velocity   ! [m/s]
      real(fp), intent(in) :: wind_speed_10m           ! [m/s]
      real(fp), intent(in) :: surface_temperature      ! [K]
      real(fp), intent(in) :: air_density              ! [kg/m³]
      real(fp), intent(in) :: soil_moisture_top        ! [m³/m³]
      real(fp), intent(in) :: surface_roughness        ! [m]

      ! Surface property fields (scalars - single column)
      real(fp), intent(in) :: clay_fraction            ! [0-1]
      real(fp), intent(in) :: sand_fraction            ! [0-1]
      integer, intent(in)  :: soil_type                ! [1-19]
      real(fp), intent(in) :: ocean_fraction           ! [0-1]
      real(fp), intent(in) :: land_ice_fraction        ! [0-1]
      real(fp), intent(in) :: snow_fraction            ! [0-1]
      real(fp), intent(in) :: sediment_supply_map      ! [0-1]
      real(fp), intent(in) :: roughness_drag           ! [0-1]

      ! Size distribution [n_dust_species]
      real(fp), intent(in) :: effective_radius(:)      ! [m]
      real(fp), intent(in) :: lower_radius(:)          ! [m]
      real(fp), intent(in) :: upper_radius(:)          ! [m]

      ! Configuration parameters
      real(fp), intent(in) :: alpha_scale, beta_scale
      integer, intent(in)  :: moist_opt, drag_opt, horiz_flux_opt

      ! Outputs (scalars and 1D array for size bins)
      real(fp), intent(out) :: total_emission              ! [μg/m²/s]
      real(fp), intent(out) :: emission_per_species(:)     ! [μg/m²/s] [n_dust_species]
      integer, intent(out)  :: rc

      ! Local variables for Ginoux scheme computation
      integer :: n
      logical :: do_dust
      real(fp) :: u_threshold, topographic_factor, bare_soil_factor
      real(fp) :: erodibility_factor, wind_factor, emission_strength
      real(fp) :: distribution(n_dust_species)

      ! Ginoux scheme constants
      real(fp), parameter :: wind_threshold = 6.5_fp    ! Threshold wind speed [m/s]
      real(fp), parameter :: clay_limit = 0.2_fp        ! Clay fraction limit
      real(fp), parameter :: sand_optimal = 0.3_fp      ! Optimal sand fraction

      ! Initialize
      rc = 0
      total_emission = ZERO
      emission_per_species = ZERO

      ! Check if dust emission should occur
      do_dust = .true.

      ! Don't do dust over water, ice, or snow
      if (ocean_fraction > 0.1_fp .or. land_ice_fraction > 0.1_fp .or. snow_fraction > 0.1_fp) then
         do_dust = .false.
      endif

      ! Don't do dust over frozen soil
      if (surface_temperature <= 273.15_fp) do_dust = .false.

      ! Don't do dust over inappropriate soil types
      if (soil_type == 15 .or. soil_type == 16 .or. soil_type == 18) then
         do_dust = .false.
      endif

      if (.not. do_dust) return

      ! Ginoux scheme: Compute topographic source function
      ! Based on sediment supply map and surface characteristics
      topographic_factor = sediment_supply_map * alpha_scale

      ! Compute bare soil factor (inverse of clay content up to a limit)
      if (clay_fraction < clay_limit) then
         bare_soil_factor = (1.0_fp - clay_fraction / clay_limit)
      else
         bare_soil_factor = 0.0_fp  ! Too much clay inhibits emission
      endif

      ! Compute erodibility factor based on sand content
      ! Ginoux scheme favors moderate sand content
      if (sand_fraction < sand_optimal) then
         erodibility_factor = sand_fraction / sand_optimal
      else
         erodibility_factor = (1.0_fp - sand_fraction) / (1.0_fp - sand_optimal)
      endif

      ! Compute wind factor (Ginoux uses a different formulation than Fengsha)
      if (wind_speed_10m > wind_threshold) then
         wind_factor = (wind_speed_10m - wind_threshold)**2 * wind_speed_10m
      else
         wind_factor = 0.0_fp
      endif

      ! Apply soil moisture reduction (simple exponential)
      if (soil_moisture_top > 0.0_fp) then
         wind_factor = wind_factor * exp(-beta_scale * soil_moisture_top)
      endif

      ! Compute emission strength [μg/m²/s]
      emission_strength = topographic_factor * bare_soil_factor * erodibility_factor * &
                         wind_factor * air_density * 1.0e6_fp  ! Convert to μg

      ! Get size distribution for dust bins
      call KokDistribution(effective_radius, lower_radius, upper_radius, distribution)

      ! Total emission
      total_emission = emission_strength

      ! Distribute among size bins
      do n = 1, n_dust_species
         emission_per_species(n) = distribution(n) * total_emission
      enddo

   end subroutine run_ginoux_scheme

end module DustScheme_Ginoux_Mod



