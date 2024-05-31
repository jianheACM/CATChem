! Jian.He@noaa.gov, 05/2023
! Move to parameters folder for CATChem

module catchem_constants

  implicit none

  public

  integer, parameter :: kind_chem = 8

  real(kind=kind_chem),parameter:: con_g      =9.80665e+0_kind_chem                !< gravity (\f$m/s^{2}\f$)
  real(kind=kind_chem),parameter:: con_rd     =2.8705e+2_kind_chem                 !< gas constant air (\f$J/kg/K\f$)
  real(kind=kind_chem),parameter:: con_rv     =4.6150e+2_kind_chem                 !< gas constant H2O (\f$J/kg/K\f$)
  real(kind=kind_chem),parameter:: con_cp     =1.0046e+3_kind_chem                 !< spec heat air at p (\f$J/kg/K\f$)
  real(kind=kind_chem),parameter:: con_cv     =7.1760e+2_kind_chem                 !< spec heat air at v (\f$J/kg/K\f$)
  real(kind=kind_chem),parameter:: con_pi     =4.0d0*atan(1.0d0)                   !< pi

  real(kind=kind_chem),parameter :: epsilc     = 1.e-30_kind_chem   

  real(kind=kind_chem), public, parameter :: RADIUS = 6.3712e+6_kind_chem           !< Radius of the Earth [m]
  real(kind=kind_chem), public, parameter :: PI_8   = 3.1415926535897931_kind_chem  !< Ratio of circle circumference to diameter [N/A]
  real(kind=kind_chem), public, parameter :: PI     = 3.1415926535897931_kind_chem  !< Ratio of circle circumference to diameter [N/A] (REAL(KIND=8))
  real(kind=kind_chem), public, parameter :: OMEGA  = 7.2921e-5_kind_chem   !< Rotation rate of the Earth [1/s]
  real(kind=kind_chem), public, parameter :: GRAV   = 9.80665_kind_chem     !< Acceleration due to gravity [m/s^2]
  real(kind=kind_chem), public, parameter :: GRAV_8 = 9.80665_kind_chem     !< Acceleration due to gravity [m/s^2] (REAL(KIND=8))
  real(kind=kind_chem), public, parameter :: RDGAS  = 287.05_kind_chem      !< Gas constant for dry air [J/kg/deg]
  real(kind=kind_chem), public, parameter :: RVGAS  = 461.50_kind_chem      !< Gas constant for water vapor [J/kg/deg]
! Extra:
  real(kind=kind_chem), public, parameter :: HLV      = 2.5e6_kind_chem     !< Latent heat of evaporation [J/kg]
  real(kind=kind_chem), public, parameter :: HLF      = 3.3358e5_kind_chem  !< Latent heat of fusion [J/kg]
  real(kind=kind_chem), public, parameter :: con_cliq = 4.1855e+3_kind_chem !< spec heat H2O liq [J/kg/K]
  real(kind=kind_chem), public, parameter :: con_csol = 2.1060e+3_kind_chem !< spec heat H2O ice [J/kg/K]
  real(kind=kind_chem), public, parameter :: CP_AIR = 1004.6_kind_chem      !< Specific heat capacity of dry air at constant pressure [J/kg/deg]
  real(kind=kind_chem), public, parameter :: KAPPA  = RDGAS/CP_AIR        !< RDGAS / CP_AIR [dimensionless]
  real(kind=kind_chem), public, parameter :: TFREEZE = 273.15_kind_chem     !< Freezing temperature of fresh water [K]

  real(kind=kind_chem), public, parameter :: STEFAN  = 5.6734e-8_kind_chem !< Stefan-Boltzmann constant [W/m^2/deg^4]

  real(kind=kind_chem), public, parameter :: CP_VAPOR = 4.0_kind_chem*RVGAS      !< Specific heat capacity of water vapor at constant pressure [J/kg/deg]
  real(kind=kind_chem), public, parameter :: CP_OCEAN = 3989.24495292815_kind_chem !< Specific heat capacity taken from McDougall (2002) 
                                                               !! "Potential
                                                               !Enthalpy ..."
                                                               ![J/kg/deg]
  real(kind=kind_chem), public, parameter :: RHO0    = 1.035e3_kind_chem  !< Average density of sea water [kg/m^3]
  real(kind=kind_chem), public, parameter :: RHO0R   = 1.0_kind_chem/RHO0 !< Reciprocal of average density of sea water [m^3/kg]
  real(kind=kind_chem), public, parameter :: RHO_CP  = RHO0*CP_OCEAN    !< (kg/m^3)*(cal/kg/deg C)(joules/cal) = (joules/m^3/deg C) [J/m^3/deg]

  real(kind=kind_chem), public, parameter :: ES0 = 1.0_kind_chem        !< Humidity factor. Controls the humidity content of the atmosphere through
                                                    !! the Saturation Vapour
                                                    !Pressure expression when
                                                    !using DO_SIMPLE.
                                                    ![dimensionless]
  real(kind=kind_chem), public, parameter :: DENS_H2O = 1000._kind_chem !< Density of liquid water [kg/m^3]
  real(kind=kind_chem), public, parameter :: HLS = HLV + HLF          !< Latent heat of sublimation [J/kg]

  real(kind=kind_chem), public, parameter :: WTMAIR   = 2.896440E+01_kind_chem   !< Molecular weight of air [AMU]
  real(kind=kind_chem), public, parameter :: WTMH2O   = WTMAIR*(RDGAS/RVGAS)   !< Molecular weight of water [AMU]
  real(kind=kind_chem), public, parameter :: WTMOZONE =  47.99820_kind_chem      !< Molecular weight of ozone [AMU]
  real(kind=kind_chem), public, parameter :: WTMC     =  12.00000_kind_chem      !< Molecular weight of carbon [AMU]
  real(kind=kind_chem), public, parameter :: WTMCO2   =  44.00995_kind_chem      !< Molecular weight of carbon dioxide [AMU]
  real(kind=kind_chem), public, parameter :: WTMCH4   =  16.0425_kind_chem       !< Molecular weight of methane [AMU]
  real(kind=kind_chem), public, parameter :: WTMO2    =  31.9988_kind_chem       !< Molecular weight of molecular oxygen [AMU]
  real(kind=kind_chem), public, parameter :: WTMCFC11 = 137.3681_kind_chem       !< Molecular weight of CFC-11 (CCl3F) [AMU]
  real(kind=kind_chem), public, parameter :: WTMCFC12 = 120.9135_kind_chem       !< Molecular weight of CFC-21 (CCl2F2) [AMU]
  real(kind=kind_chem), public, parameter :: WTMN     =  14.0067_kind_chem       !< Molecular weight of Nitrogen [AMU]
  real(kind=kind_chem), public, parameter :: DIFFAC   = 1.660000E+00_kind_chem   !< Diffusivity factor [dimensionless]
  real(kind=kind_chem), public, parameter :: AVOGNO   = 6.023000E+23_kind_chem   !< Avogadro's number [atoms/mole]
  real(kind=kind_chem), public, parameter :: PSTD     = 1.013250E+06_kind_chem   !< Mean sea level pressure [dynes/cm^2]
  real(kind=kind_chem), public, parameter :: PSTD_MKS = 101325.0_kind_chem       !< Mean sea level pressure [N/m^2]

  real(kind=kind_chem), public, parameter :: SECONDS_PER_DAY    = 8.640000E+04_kind_chem !< Seconds in a day [s]
  real(kind=kind_chem), public, parameter :: SECONDS_PER_HOUR   = 3600._kind_chem        !< Seconds in an hour [s]
  real(kind=kind_chem), public, parameter :: SECONDS_PER_MINUTE = 60._kind_chem          !< Seconds in a minute [s]
  real(kind=kind_chem), public, parameter :: RAD_TO_DEG         = 180._kind_chem/PI      !< Degrees per radian [deg/rad]
  real(kind=kind_chem), public, parameter :: DEG_TO_RAD         = PI/180._kind_chem      !< Radians per degree [rad/deg]
  real(kind=kind_chem), public, parameter :: RADIAN             = RAD_TO_DEG           !< Equal to RAD_TO_DEG for backward compatability. [rad/deg]
  real(kind=kind_chem), public, parameter :: ALOGMIN            = -50.0_kind_chem        !< Minimum value allowed as argument to log function [N/A]
  real(kind=kind_chem), public, parameter :: EPSLN              = 1.0e-40_kind_chem      !< A small number to prevent divide by zero exceptions [N/A]

  real(kind=kind_chem), public, parameter :: RADCON = ((1.0E+02*GRAV)/(1.0E+04*CP_AIR))*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                                      !heating rate in degrees per day [deg sec/(cm day)]
  real(kind=kind_chem), public, parameter :: RADCON_MKS  = (GRAV/CP_AIR)*SECONDS_PER_DAY !< Factor used to convert flux divergence to
                                                                       !heating rate in degrees per day [deg sec/(cm day)]
  real(kind=kind_chem), public, parameter :: O2MIXRAT    = 2.0953E-01_kind_chem !< Mixing ratio of molecular oxygen in air [dimensionless]
  real(kind=kind_chem), public, parameter :: RHOAIR      = 1.292269_kind_chem   !< Reference atmospheric density [kg/m^3]
  real(kind=kind_chem), public, parameter :: VONKARM     = 0.40_kind_chem       !< Von Karman constant [dimensionless]
  real(kind=kind_chem), public, parameter :: C2DBARS     = 1.e-4_kind_chem      !< Converts rho*g*z (in mks) to dbars: 1dbar = 10^4 (kg/m^3)(m/s^2)m [dbars]
  real(kind=kind_chem), public, parameter :: KELVIN      = 273.15_kind_chem     !< Degrees Kelvin at zero Celsius [K]

end module catchem_constants
