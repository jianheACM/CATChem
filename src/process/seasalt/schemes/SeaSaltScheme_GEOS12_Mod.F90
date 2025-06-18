!> \file SeaSaltScheme_GEOS12_Mod.F90
!! \brief GEOS-Chem 2012 sea salt emission scheme
!! \ingroup catchem_seasalt_process
!!
!! \author Barry Baker
!! \date 05/2024
!!
!! This module implements the GEOS-Chem 2012 sea salt emission scheme for calculating
!! sea salt aerosol emissions in the CATChem atmospheric chemistry model.
!!
!! \details
!! The GEOS-Chem 2012 scheme provides an improved parameterization for sea salt
!! emissions based on extensive observational constraints and updated understanding
!! of sea salt emission processes.
!!
!! \section geos12_reference Reference
!! Jaeglé, L., Quinn, P. K., Bates, T. S., Alexander, B., and Lin, J.-T.: Global distribution of sea salt aerosols:
!! new constraints from in situ and remote sensing observations, Atmos. Chem. Phys., 11, 3137–3157,
!! https://doi.org/10.5194/acp-11-3137-2011, 2011.
!!
MODULE SeaSaltScheme_GEOS12_Mod

   USE Precision_Mod
   USE constants
   USE SeaSaltCommon_Mod

   IMPLICIT NONE
   PRIVATE

   ! Public procedures
   PUBLIC :: seasalt_scheme_geos12

CONTAINS

   !> GEOS-12 sea salt emission scheme
   !!
   !! Computes sea salt mass and number emissions for a single column
   !! using the GEOS-12 parameterization.
   !!
   !! \param frocean Ocean fraction [0-1]
   !! \param frseaice Sea ice fraction [0-1]
   !! \param ustar Friction velocity [m/s]
   !! \param sst Sea surface temperature [K]
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
   SUBROUTINE seasalt_scheme_geos12( &
      frocean, frseaice, ustar, sst, scale_factor, &
      upper_radius, lower_radius, effective_radius, density, &
      mass_emission, number_emission, &
      total_mass_emission, total_number_emission, rc)

      IMPLICIT NONE

      ! Arguments
      REAL(fp), INTENT(IN) :: frocean, frseaice, ustar, sst, scale_factor
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
      REAL(fp) :: dry_radius, wet_radius, bin_width, sub_bin_width
      REAL(fp) :: particle_volume, particle_mass
      REAL(fp) :: emission_rate, number_flux, mass_flux
      REAL(fp) :: sst_correction
      REAL(fp) :: u_star_cm  ! Friction velocity in cm/s

      ! Scheme parameters
      REAL(fp), PARAMETER :: a_factor = 3.0_fp
      REAL(fp), PARAMETER :: b_factor = 0.5_fp
      REAL(fp), PARAMETER :: r_power = 1.05_fp
      REAL(fp), PARAMETER :: exp_power = 1.19_fp
      REAL(fp), PARAMETER :: u_power = 3.41_fp

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
      IF (open_ocean_fraction <= 0.0_fp .OR. ustar <= 0.0_fp) RETURN

      ! Convert friction velocity to cm/s
      u_star_cm = ustar * 100.0_fp

      ! Apply SST correction
      sst_correction = jeagle_sst_correction(sst)

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

            ! Calculate emission rate using GEOS-12 parameterization
            emission_rate = seasalt_emission_gong( &
               wet_radius, sub_bin_width, u_star_cm, scale_factor, &
               a_factor, b_factor, r_power, exp_power, u_power)

            ! Apply corrections
            emission_rate = emission_rate * sst_correction * open_ocean_fraction

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

   END SUBROUTINE seasalt_scheme_geos12

END MODULE SeaSaltScheme_GEOS12_Mod
