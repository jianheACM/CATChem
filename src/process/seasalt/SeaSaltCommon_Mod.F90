!> \file SeaSaltCommon_Mod.F90
!! \brief Common utilities for sea salt emission calculations
!! \ingroup process_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 2.0
!!
!! This module provides pure utility functions for sea salt emission calculations.
!! All functions are independent of StateContainer and can be used by different
!! sea salt emission schemes. State management is handled at the process level
!! through the StateContainer.
!!
!! The module contains implementations of various sea salt emission parameterizations
!! from the literature, including Monahan et al. (1986), Gong (2003), and others.
!!
!! \see Monahan, E. C., D. E. Spiel, and K. L. Davidson (1986), A model of marine
!!      aerosol generation via whitecaps and wave disruption, in Oceanic Whitecaps,
!!      edited by E. C. Monahan and G. Mac Niocaill, pp. 167–174, D. Reidel, Norwell, Mass.
!! \see Gong, S. L. (2003), A parameterization of sea-salt aerosol source function
!!      for sub- and super-micron particles, Global Biogeochem. Cycles, 17(4), 1097
!!
module SeaSaltCommon_Mod

   use precision_mod, only: fp
   use error_mod, only: ErrorManagerType, ERROR_INVALID_INPUT, ERROR_COMPUTATION, CC_SUCCESS

   implicit none
   private

   ! Public procedures - all utility functions for sea salt calculations
   public :: Monahan_SeaSaltFlux
   public :: Gong_SeaSaltFlux
   public :: Whitecap_Coverage_Monahan
   public :: Whitecap_Coverage_Salisbury
   public :: Particle_SizeDistribution_Monahan
   public :: Particle_SizeDistribution_Gong
   public :: Wind_Speed_Correction_10m
   public :: Sea_Surface_Temperature_Correction

contains

   !> \brief Calculate sea salt flux using Monahan et al. (1986) parameterization
   !!
   !! Monahan, E. C., D. E. Spiel, and K. L. Davidson (1986), A model of marine
   !! aerosol generation via whitecaps and wave disruption, in Oceanic Whitecaps,
   !! edited by E. C. Monahan and G. Mac Niocaill, pp. 167–174, D. Reidel, Norwell, Mass.
   !!
   !! \param particle_radius Particle radius at 80% RH [μm]
   !! \param wind_speed_10m Wind speed at 10m height [m/s]
   !! \param whitecap_coverage Whitecap coverage fraction [-]
   !! \param sea_salt_flux Calculated sea salt flux [particles/m²/s/μm]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Monahan_SeaSaltFlux(particle_radius, wind_speed_10m, whitecap_coverage, &
                                  sea_salt_flux, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: particle_radius       !< Particle radius at 80% RH [μm]
      real(fp), intent(in) :: wind_speed_10m        !< Wind speed at 10m [m/s]
      real(fp), intent(in) :: whitecap_coverage     !< Whitecap coverage fraction [-]
      real(fp), intent(out) :: sea_salt_flux        !< Sea salt flux [particles/m²/s/μm]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp), parameter :: A_MONAHAN = 1.373_fp   !< Monahan coefficient A
      real(fp), parameter :: B_MONAHAN = 3.0_fp     !< Monahan exponent B
      real(fp) :: log_radius                        !< Log of particle radius

      call error_mgr%push_context('Monahan_SeaSaltFlux', 'Computing Monahan sea salt flux')

      ! Initialize
      sea_salt_flux = 0.0_fp
      rc = CC_SUCCESS

      ! Input validation
      if (particle_radius <= 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Particle radius must be positive', rc)
         return
      endif

      if (wind_speed_10m < 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Wind speed must be non-negative', rc)
         return
      endif

      if (whitecap_coverage < 0.0_fp .or. whitecap_coverage > 1.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Whitecap coverage must be between 0 and 1', rc)
         return
      endif

      ! Calculate flux using Monahan parameterization
      ! dF/dr = A * r^(-B) * (1 + 0.057*r^1.05) * 10^(1.19*exp(-B²))
      ! where r is particle radius in μm

      log_radius = log10(particle_radius)

      sea_salt_flux = A_MONAHAN * (particle_radius**(-B_MONAHAN)) * &
                     (1.0_fp + 0.057_fp * (particle_radius**1.05_fp)) * &
                     (10.0_fp**(1.19_fp * exp(-((log_radius)**2)))) * &
                     whitecap_coverage

      call error_mgr%pop_context()

   end subroutine Monahan_SeaSaltFlux

   !> \brief Calculate sea salt flux using Gong (2003) parameterization
   !!
   !! Gong, S. L. (2003), A parameterization of sea-salt aerosol source function
   !! for sub- and super-micron particles, Global Biogeochem. Cycles, 17(4), 1097
   !!
   !! \param particle_radius Particle radius at 80% RH [μm]
   !! \param wind_speed_10m Wind speed at 10m height [m/s]
   !! \param whitecap_coverage Whitecap coverage fraction [-]
   !! \param sea_salt_flux Calculated sea salt flux [particles/m²/s/μm]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Gong_SeaSaltFlux(particle_radius, wind_speed_10m, whitecap_coverage, &
                               sea_salt_flux, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: particle_radius       !< Particle radius at 80% RH [μm]
      real(fp), intent(in) :: wind_speed_10m        !< Wind speed at 10m [m/s]
      real(fp), intent(in) :: whitecap_coverage     !< Whitecap coverage fraction [-]
      real(fp), intent(out) :: sea_salt_flux        !< Sea salt flux [particles/m²/s/μm]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables - Gong (2003) parameters
      real(fp), parameter :: A_GONG = 4.7e6_fp      !< Gong coefficient A
      real(fp), parameter :: ALPHA = -3.0_fp        !< Size distribution exponent
      real(fp), parameter :: R0 = 2.1_fp            !< Reference radius [μm]
      real(fp), parameter :: SIGMA = 0.7_fp         !< Standard deviation
      real(fp) :: size_factor                       !< Size-dependent factor

      call error_mgr%push_context('Gong_SeaSaltFlux', 'Computing Gong sea salt flux')

      ! Initialize
      sea_salt_flux = 0.0_fp
      rc = CC_SUCCESS

      ! Input validation
      if (particle_radius <= 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Particle radius must be positive', rc)
         return
      endif

      if (wind_speed_10m < 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Wind speed must be non-negative', rc)
         return
      endif

      if (whitecap_coverage < 0.0_fp .or. whitecap_coverage > 1.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Whitecap coverage must be between 0 and 1', rc)
         return
      endif

      ! Calculate flux using Gong parameterization
      ! Modified size distribution for both sub- and super-micron particles

      if (particle_radius <= 2.0_fp) then
         ! Sub-micron parameterization
         size_factor = (particle_radius**ALPHA) * &
                      exp(-0.5_fp * ((log(particle_radius) - log(R0))/SIGMA)**2)
      else
         ! Super-micron parameterization (enhanced production)
         size_factor = (particle_radius**ALPHA) * &
                      exp(-0.5_fp * ((log(particle_radius) - log(R0))/(SIGMA*1.5_fp))**2) * &
                      (1.0_fp + 0.2_fp * particle_radius)
      endif

      sea_salt_flux = A_GONG * size_factor * whitecap_coverage

      call error_mgr%pop_context()

   end subroutine Gong_SeaSaltFlux

   !> \brief Calculate whitecap coverage using Monahan and O'Muircheartaigh (1980)
   !!
   !! Monahan, E. C. and I. O'Muircheartaigh (1980), Optimal power-law description
   !! of oceanic whitecap coverage dependence on wind speed, J. Phys. Oceanogr., 10, 2094-2099
   !!
   !! \param wind_speed_10m Wind speed at 10m height [m/s]
   !! \param whitecap_coverage Calculated whitecap coverage fraction [-]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Whitecap_Coverage_Monahan(wind_speed_10m, whitecap_coverage, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: wind_speed_10m        !< Wind speed at 10m [m/s]
      real(fp), intent(out) :: whitecap_coverage    !< Whitecap coverage fraction [-]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp), parameter :: U_THRESHOLD = 3.84_fp  !< Threshold wind speed [m/s]
      real(fp), parameter :: WC_COEFF = 3.84e-6_fp  !< Whitecap coefficient

      call error_mgr%push_context('Whitecap_Coverage_Monahan', 'Computing Monahan whitecap coverage')

      ! Initialize
      whitecap_coverage = 0.0_fp
      rc = CC_SUCCESS

      ! Input validation
      if (wind_speed_10m < 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Wind speed must be non-negative', rc)
         return
      endif

      ! Calculate whitecap coverage: W = 3.84e-6 * U10^3.41
      ! Only for wind speeds above threshold
      if (wind_speed_10m > U_THRESHOLD) then
         whitecap_coverage = WC_COEFF * (wind_speed_10m**3.41_fp)

         ! Cap at physically reasonable maximum
         whitecap_coverage = min(whitecap_coverage, 0.5_fp)
      endif

      call error_mgr%pop_context()

   end subroutine Whitecap_Coverage_Monahan

   !> \brief Calculate whitecap coverage using Salisbury et al. (2013)
   !!
   !! Salisbury, D. J., et al. (2013), A global climatology of whitecap fraction
   !! based on satellite observations, J. Geophys. Res. Oceans, 118, 1824–1831
   !!
   !! \param wind_speed_10m Wind speed at 10m height [m/s]
   !! \param sea_surface_temp Sea surface temperature [K]
   !! \param whitecap_coverage Calculated whitecap coverage fraction [-]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Whitecap_Coverage_Salisbury(wind_speed_10m, sea_surface_temp, &
                                          whitecap_coverage, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: wind_speed_10m        !< Wind speed at 10m [m/s]
      real(fp), intent(in) :: sea_surface_temp      !< Sea surface temperature [K]
      real(fp), intent(out) :: whitecap_coverage    !< Whitecap coverage fraction [-]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp), parameter :: U_THRESHOLD = 3.0_fp   !< Threshold wind speed [m/s]
      real(fp), parameter :: T_REF = 273.15_fp      !< Reference temperature [K]
      real(fp) :: temp_factor                       !< Temperature correction factor

      call error_mgr%push_context('Whitecap_Coverage_Salisbury', 'Computing Salisbury whitecap coverage')

      ! Initialize
      whitecap_coverage = 0.0_fp
      rc = CC_SUCCESS

      ! Input validation
      if (wind_speed_10m < 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Wind speed must be non-negative', rc)
         return
      endif

      if (sea_surface_temp < 200.0_fp .or. sea_surface_temp > 350.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Sea surface temperature out of reasonable range', rc)
         return
      endif

      ! Calculate whitecap coverage with temperature dependence
      if (wind_speed_10m > U_THRESHOLD) then
         temp_factor = 1.0_fp + 0.01_fp * (sea_surface_temp - T_REF)
         whitecap_coverage = 3.97e-6_fp * (wind_speed_10m**3.52_fp) * temp_factor

         ! Cap at physically reasonable maximum
         whitecap_coverage = max(0.0_fp, min(whitecap_coverage, 0.5_fp))
      endif

      call error_mgr%pop_context()

   end subroutine Whitecap_Coverage_Salisbury

   !> \brief Calculate particle size distribution parameters for Monahan scheme
   !!
   !! \param particle_radius Particle radius [μm]
   !! \param distribution_value Size distribution value [particles/μm]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Particle_SizeDistribution_Monahan(particle_radius, distribution_value, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: particle_radius       !< Particle radius [μm]
      real(fp), intent(out) :: distribution_value   !< Distribution value [particles/μm]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp) :: log_r                             !< Log of radius

      call error_mgr%push_context('Particle_SizeDistribution_Monahan', &
                                  'Computing Monahan size distribution')

      ! Initialize
      distribution_value = 0.0_fp
      rc = CC_SUCCESS

      ! Input validation
      if (particle_radius <= 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Particle radius must be positive', rc)
         return
      endif

      ! Monahan size distribution: dN/dr ∝ r^(-3) * (1 + 0.057*r^1.05) * exp(-B²)
      log_r = log10(particle_radius)

      distribution_value = (particle_radius**(-3.0_fp)) * &
                          (1.0_fp + 0.057_fp * (particle_radius**1.05_fp)) * &
                          exp(-((log_r)**2))

      call error_mgr%pop_context()

   end subroutine Particle_SizeDistribution_Monahan

   !> \brief Calculate particle size distribution parameters for Gong scheme
   !!
   !! \param particle_radius Particle radius [μm]
   !! \param distribution_value Size distribution value [particles/μm]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Particle_SizeDistribution_Gong(particle_radius, distribution_value, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: particle_radius       !< Particle radius [μm]
      real(fp), intent(out) :: distribution_value   !< Distribution value [particles/μm]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp), parameter :: R0 = 2.1_fp            !< Reference radius [μm]
      real(fp), parameter :: SIGMA = 0.7_fp         !< Standard deviation

      call error_mgr%push_context('Particle_SizeDistribution_Gong', &
                                  'Computing Gong size distribution')

      ! Initialize
      distribution_value = 0.0_fp
      rc = CC_SUCCESS

      ! Input validation
      if (particle_radius <= 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Particle radius must be positive', rc)
         return
      endif

      ! Gong size distribution with bi-modal characteristics
      distribution_value = (particle_radius**(-3.0_fp)) * &
                           exp(-0.5_fp * ((log(particle_radius) - log(R0))/SIGMA)**2)

      call error_mgr%pop_context()

   end subroutine Particle_SizeDistribution_Gong

   !> \brief Correct wind speed to 10m height using neutral stability
   !!
   !! \param wind_speed Wind speed at reference height [m/s]
   !! \param reference_height Reference height [m]
   !! \param wind_speed_10m Corrected wind speed at 10m [m/s]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Wind_Speed_Correction_10m(wind_speed, reference_height, wind_speed_10m, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: wind_speed            !< Wind speed at reference height [m/s]
      real(fp), intent(in) :: reference_height      !< Reference height [m]
      real(fp), intent(out) :: wind_speed_10m       !< Wind speed at 10m [m/s]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp), parameter :: Z0_OCEAN = 1.5e-4_fp   !< Ocean roughness length [m]
      real(fp), parameter :: TARGET_HEIGHT = 10.0_fp !< Target height [m]
      real(fp) :: log_correction                     !< Logarithmic correction factor

      call error_mgr%push_context('Wind_Speed_Correction_10m', 'Correcting wind speed to 10m')

      ! Initialize
      wind_speed_10m = wind_speed
      rc = CC_SUCCESS

      ! Input validation
      if (wind_speed < 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Wind speed must be non-negative', rc)
         return
      endif

      if (reference_height <= Z0_OCEAN) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Reference height must be greater than roughness length', rc)
         return
      endif

      ! Apply logarithmic wind profile correction (neutral stability)
      if (abs(reference_height - TARGET_HEIGHT) > 0.1_fp) then
         log_correction = log((TARGET_HEIGHT - Z0_OCEAN) / Z0_OCEAN) / &
                         log((reference_height - Z0_OCEAN) / Z0_OCEAN)
         wind_speed_10m = wind_speed * log_correction
      endif

      call error_mgr%pop_context()

   end subroutine Wind_Speed_Correction_10m

   !> \brief Apply sea surface temperature correction to emission flux
   !!
   !! \param base_flux Base emission flux [particles/m²/s]
   !! \param sea_surface_temp Sea surface temperature [K]
   !! \param corrected_flux Temperature-corrected emission flux [particles/m²/s]
   !! \param error_mgr Error manager for error handling
   !! \param rc Return code
   !!
   !! \ingroup catchem_seasalt_process
   subroutine Sea_Surface_Temperature_Correction(base_flux, sea_surface_temp, &
                                                 corrected_flux, error_mgr, rc)
      implicit none

      ! Parameters
      real(fp), intent(in) :: base_flux             !< Base emission flux [particles/m²/s]
      real(fp), intent(in) :: sea_surface_temp      !< Sea surface temperature [K]
      real(fp), intent(out) :: corrected_flux       !< Corrected emission flux [particles/m²/s]
      type(ErrorManagerType), intent(inout) :: error_mgr !< Error manager
      integer, intent(out) :: rc                    !< Return code

      ! Local Variables
      real(fp), parameter :: T_REF = 288.15_fp      !< Reference temperature [K]
      real(fp), parameter :: TEMP_COEFF = 0.02_fp   !< Temperature coefficient [1/K]
      real(fp) :: temp_factor                       !< Temperature correction factor

      call error_mgr%push_context('Sea_Surface_Temperature_Correction', &
                                  'Applying SST correction to flux')

      ! Initialize
      corrected_flux = base_flux
      rc = CC_SUCCESS

      ! Input validation
      if (base_flux < 0.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Base flux must be non-negative', rc)
         return
      endif

      if (sea_surface_temp < 200.0_fp .or. sea_surface_temp > 350.0_fp) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
              'Sea surface temperature out of reasonable range', rc)
         return
      endif

      ! Apply temperature correction: F_corrected = F_base * (1 + α*(T - T_ref))
      temp_factor = 1.0_fp + TEMP_COEFF * (sea_surface_temp - T_REF)
      corrected_flux = base_flux * max(0.1_fp, temp_factor)  ! Prevent negative corrections

      call error_mgr%pop_context()

   end subroutine Sea_Surface_Temperature_Correction

end module SeaSaltCommon_Mod
