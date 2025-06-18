!> \file DustScheme_Fengsha_Mod.F90
!! \brief Modern Fengsha dust emission scheme implementation - StateContainer independent
!! \ingroup dust_schemes
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module implements the FENGSHA dust emission scheme using modern
!! architecture but independent of StateContainer to ease science developer burden.
!! Based on the original ccpr_scheme_fengsha_mod.F90 but simplified for modern use.
!!
!! \module DustScheme_Fengsha_Mod
!!!>
module DustScheme_Fengsha_Mod
   use precision_mod, only: fp, ZERO
   use constants, only: g0
   use DustCommon_Mod
   implicit none
   private

   public :: run_fengsha_scheme

contains

   !> \brief Execute Fengsha dust emission scheme (single-column)
   !!
   !! This is the main interface for the FENGSHA dust emission scheme for a single column.
   !! Science developers can call this directly with scalar/1D variables without
   !! worrying about StateContainer complexity.
   !!
   !! Based on: Zhang, L., Montuoro, R., McKeen, S. A., Baker, B., et al. (2022):
   !! Development and evaluation of the Aerosol Forecast Member in NCEP's GEFS-Aerosols v1,
   !! Geosci. Model Dev., 15, 5337–5369, https://doi.org/10.5194/gmd-15-5337-2022
   !!
   subroutine run_fengsha_scheme(n_dust_species, &
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

      ! Local variables for computation (all scalars for single column)
      integer :: n
      logical :: do_dust
      real(fp) :: u_friction, u_threshold, skin_temp, soil_moisture
      real(fp) :: clay_frac, sand_frac, air_dens, ssm, rdrag
      real(fp) :: frac_ocean, frac_ice, frac_snow, z0
      real(fp) :: H, R, SEP, horiz_flux, h_to_v_ratio
      real(fp) :: fengsha_scaling, alpha_grav
      real(fp) :: distribution(n_dust_species)

      ! Constants
      real(fp), parameter :: clay_threshold = 0.2_fp
      real(fp), parameter :: kvhmax = 2.0e-4_fp      ! Max vertical-to-horizontal flux ratio
      real(fp), parameter :: z0s = 0.00003_fp        ! Smooth roughness length [m]

      ! Initialize
      rc = 0
      total_emission = ZERO
      emission_per_species = ZERO
      alpha_grav = alpha_scale / g0

      ! Extract single-column values
      u_friction = wind_friction_velocity
      skin_temp = surface_temperature
      soil_moisture = soil_moisture_top
      clay_frac = clay_fraction
      sand_frac = sand_fraction
      air_dens = air_density
      ssm = sediment_supply_map
      rdrag = roughness_drag
      frac_ocean = ocean_fraction
      frac_ice = land_ice_fraction
      frac_snow = snow_fraction
      z0 = surface_roughness

      ! Check if dust emission should occur
      do_dust = .true.

      ! Don't do dust over bedrock, lava, or permanent ice
      if (soil_type == 15 .or. soil_type == 16 .or. soil_type == 18) then
         do_dust = .false.
      endif

      ! Check for valid sediment supply and drag partition
      if (ssm < 0.15_fp .or. ssm > 1.0_fp) do_dust = .false.
      if (rdrag < 0.15_fp .or. rdrag > 1.0_fp) do_dust = .false.

      ! Don't do dust over frozen soil
      if (skin_temp <= 273.15_fp) do_dust = .false.

      if (.not. do_dust) return

      ! Compute soil moisture attenuation factor
      select case (moist_opt)
      case (1)
         call Fecan_SoilMoisture(clay_frac, sand_frac, soil_moisture, H)
      case (2)
         call Shao_SoilMoisture(soil_moisture, H)
      case default
         H = 1.0_fp
      end select

      ! Get size distribution for dust bins
      call KokDistribution(effective_radius, lower_radius, upper_radius, distribution)

      ! Compute vertical-to-horizontal mass flux ratio
      if (clay_frac > clay_threshold) then
         h_to_v_ratio = kvhmax
      else
         h_to_v_ratio = 10.0_fp**(13.4_fp * clay_frac - 6.0_fp)
      endif

      ! Compute soil erosion potential (simple version)
      call Soil_Erosion_Potential(clay_frac, sand_frac, (1.0_fp - clay_frac - sand_frac), SEP)

      ! Compute drag partition factor
      select case (drag_opt)
      case (1)
         R = rdrag  ! Use input drag partition
      case (2)
         call MB95_DragPartition(z0, z0s, R)
      case default
         R = 1.0_fp
      end select

      ! Compute threshold friction velocity (simplified)
      u_threshold = 0.2_fp  ! Simplified - could use MB97_threshold_velocity

      ! Compute horizontal mass flux
      select case (horiz_flux_opt)
      case (1)
         call Kawamura_HorizFlux(u_friction, u_threshold, R, H, horiz_flux)
      case (2)
         call Draxler_HorizFlux(u_friction, u_threshold, R, H, horiz_flux)
      case default
         horiz_flux = 0.0_fp
      end select

      ! Compute Fengsha scaling factor
      fengsha_scaling = (1.0_fp - frac_ice) * (1.0_fp - frac_snow) * &
                        (1.0_fp - frac_ocean) * alpha_grav * &
                        (ssm ** beta_scale) * air_dens * 1.0e9_fp

      ! Compute total emission [μg/m²/s]
      total_emission = fengsha_scaling * horiz_flux * h_to_v_ratio

      ! Distribute among size bins
      do n = 1, n_dust_species
         emission_per_species(n) = distribution(n) * total_emission
      enddo

   end subroutine run_fengsha_scheme

end module DustScheme_Fengsha_Mod
