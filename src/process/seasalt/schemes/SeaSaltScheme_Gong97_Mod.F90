!> \file ccpr_scheme_gong97_mod.F90
!! \brief Gong 1997 sea salt emission scheme
!! \ingroup catchem_seasalt_process
!!
!! \author Barry Baker
!! \date 05/2024
!!
!! This module implements the Gong 1997 sea salt emission scheme for calculating
!! sea salt aerosol emissions in the CATChem atmospheric chemistry model.
!!
!! \details
!! The Gong 1997 scheme calculates sea salt emissions based on surface wind speed
!! and provides size-resolved emission fluxes. It is one of the foundational
!! parameterizations for sea salt emission modeling.
!!
!! \section gong97_reference Reference
!! Gong, S. L., L. A. Barrie, and J.-P. Blanchet (1997), Modeling sea-salt aerosols in the atmosphere:
!! 1. Model development, J. Geophys. Res., 102(D3), 3805–3818, doi:10.1029/96JD02953.
!!
MODULE SeaSaltScheme_Gong97_Mod

   USE Precision_Mod
   USE constants
   USE SeaSaltCommon_Mod

   IMPLICIT NONE
   PRIVATE

   ! Public procedures
   PUBLIC :: seasalt_scheme_gong97

CONTAINS

   !> Gong 1997 sea salt emission scheme
   !!
   !! Computes sea salt mass and number emissions for a single column
   !! using the Gong 1997 parameterization.
   !!
   !! \param frocean Ocean fraction [0-1]
   !! \param frseaice Sea ice fraction [0-1]
   !! \param u10m 10-m u-wind [m/s]
   !! \param v10m 10-m v-wind [m/s]
   !! \param sst Sea surface temperature [K]
   !! \param weibull_flag Apply Weibull distribution
   !! \param scale_factor Emission scale factor [dimensionless]
   !! \param upper_radius Upper bin radius [um]
   !! \param lower_radius Lower bin radius [um]
   !! \param effective_radius Effective bin radius [um]
   !! \param density Sea salt density [kg/m³]
   !! \param mass_emission Output mass emission [ug/m²/s]
   !! \param number_emission Output number emission [#/m²/s]
   !! \param total_mass_emission Total mass emission [ug/m²/s]
   !! \param total_number_emission Total number emission [#/m²/s]
   !! \param rc Return code
   SUBROUTINE seasalt_scheme_gong97( &
      frocean, frseaice, u10m, v10m, sst, weibull_flag, scale_factor, &
      upper_radius, lower_radius, effective_radius, density, &
      mass_emission, number_emission, &
      total_mass_emission, total_number_emission, rc)

      IMPLICIT NONE

      ! Arguments
      REAL(fp), INTENT(IN) :: frocean, frseaice, u10m, v10m, sst
      LOGICAL, INTENT(IN) :: weibull_flag
      REAL(fp), INTENT(IN) :: scale_factor
      REAL(fp), INTENT(IN) :: upper_radius(:), lower_radius(:)
      REAL(fp), INTENT(IN) :: effective_radius(:), density(:)
      REAL(fp), INTENT(OUT) :: mass_emission(:), number_emission(:)
      REAL(fp), INTENT(OUT) :: total_mass_emission, total_number_emission
      INTEGER, INTENT(OUT) :: rc

      ! Local variables
      INTEGER :: n_bins, i_bin, i_sub
      INTEGER, PARAMETER :: n_sub = 10  ! Number of sub-bins
      REAL(fp), PARAMETER :: r80_factor = 1.65_fp  ! RH=80% to dry radius factor
      REAL(fp) :: open_ocean_fraction
      REAL(fp) :: wind_speed_10m, effective_wind
      REAL(fp) :: dry_radius, wet_radius, bin_width, sub_bin_width
      REAL(fp) :: particle_volume, particle_mass
      REAL(fp) :: emission_rate, number_flux, mass_flux
      REAL(fp) :: whitecap_frac

      ! Gong 1997 parameters
      REAL(fp), PARAMETER :: a_factor = 3.0_fp
      REAL(fp), PARAMETER :: b_factor = 0.38_fp
      REAL(fp), PARAMETER :: r_power = 1.05_fp
      REAL(fp), PARAMETER :: exp_power = 1.19_fp
      REAL(fp), PARAMETER :: w_power = 3.41_fp

      ! Weibull parameters
      REAL(fp), PARAMETER :: weibull_shape = 2.0_fp
      REAL(fp), PARAMETER :: weibull_scale = 1.2_fp

      ! Initialize
      rc = 0
      n_bins = SIZE(upper_radius)

      ! Initialize output arrays
      mass_emission(:) = 0.0_fp
      number_emission(:) = 0.0_fp
      total_mass_emission = 0.0_fp
      total_number_emission = 0.0_fp

      ! Check if we have ocean conditions
      open_ocean_fraction = frocean * (1.0_fp - frseaice)
      IF (open_ocean_fraction <= 0.0_fp) RETURN

      ! Calculate 10-m wind speed
      wind_speed_10m = compute_wind_speed_10m(u10m, v10m)
      IF (wind_speed_10m <= 0.0_fp) RETURN

      ! Apply Weibull distribution if requested
      IF (weibull_flag) THEN
         effective_wind = weibull_distribution(wind_speed_10m, weibull_shape, weibull_scale)
      ELSE
         effective_wind = wind_speed_10m
      END IF

      ! Calculate whitecap fraction
      whitecap_frac = whitecap_fraction(effective_wind, 1)  ! Use Monahan parameterization

      ! Loop over size bins
      DO i_bin = 1, n_bins

         ! Calculate bin width (wet radius)
         bin_width = (upper_radius(i_bin) - lower_radius(i_bin)) * r80_factor
         sub_bin_width = bin_width / REAL(n_sub, fp)

         ! Loop over sub-bins for integration
         DO i_sub = 1, n_sub

            ! Calculate wet radius for this sub-bin
            wet_radius = (lower_radius(i_bin) + &
                         (REAL(i_sub, fp) - 0.5_fp) * sub_bin_width / r80_factor) * r80_factor

            ! Calculate dry radius
            dry_radius = wet_radius / r80_factor

            ! Skip if radius is too small
            IF (dry_radius <= 0.01_fp) CYCLE

            ! Calculate particle volume and mass
            particle_volume = (4.0_fp/3.0_fp) * PI * (dry_radius * 1.0e-6_fp)**3  ! m³
            particle_mass = density(i_bin) * particle_volume * 1.0e12_fp  ! ug

            ! Calculate emission rate using Gong 1997 parameterization
            emission_rate = seasalt_emission_gong( &
               wet_radius, sub_bin_width, effective_wind, scale_factor, &
               a_factor, b_factor, r_power, exp_power, w_power)

            ! Apply whitecap fraction and ocean fraction
            emission_rate = emission_rate * whitecap_frac * open_ocean_fraction

            ! Convert to number and mass fluxes
            number_flux = emission_rate  ! #/m²/s
            mass_flux = emission_rate * particle_mass  ! ug/m²/s

            ! Add to bin totals
            number_emission(i_bin) = number_emission(i_bin) + number_flux
            mass_emission(i_bin) = mass_emission(i_bin) + mass_flux

         END DO

         ! Add to total emissions
         total_number_emission = total_number_emission + number_emission(i_bin)
         total_mass_emission = total_mass_emission + mass_emission(i_bin)

      END DO

   END SUBROUTINE seasalt_scheme_gong97

END MODULE SeaSaltScheme_Gong97_Mod
