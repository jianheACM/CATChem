!This code calculates tracer tendencies due to trop/strat chemistry
!Code is adopted from GFDL-AM4 model and restructured/updated
!by Jian He (Jian.He@noaa.gov), 09/2023

module gfdl_am4chem_mod
!-----------------------------------------------------------------------

use         tropchem_types_mod, only : tropchem_opt, tropchem_diag, &
                                       tropchem_types_init, missing_value
use         gfdl_time_utls_mod, only : time_type, &
                                       get_date, &
                                       set_date, &
                                       set_time, &
                                       days_in_year, &
                                       real_to_time_type, &
                                       time_type_to_real, &
                                       operator(+), operator(-), &
                                       time_interp
use          catchem_constants, only : grav, rdgas, WTMAIR, WTMH2O, AVOGNO, &
                                       PI, DEG_TO_RAD, SECONDS_PER_DAY, kind_chem
use              mo_chemdr_mod, only : chemdr, chemdr_init
use              mo_setsox_mod, only : setsox_init
use             mo_chemini_mod, only : chemini
use             M_TRACNAME_MOD, only : tracnam
use                MO_GRID_MOD, only : pcnstm1
use              CHEM_MODS_MOD, only : phtcnt, gascnt
use   strat_chem_utilities_mod, only : strat_chem_utilities_init, &
                                       strat_chem_dcly_dt, &
                                       strat_chem_dcly_dt_time_vary, &
                                       strat_chem_dcly_dt_endts, &
                                       strat_chem_get_aerosol, &
                                       psc_type, &
                                       strat_chem_get_h2so4, &
                                       strat_chem_get_psc, &
                                       strat_chem_destroy_psc, &
                                       strat_chem_psc_sediment, &
                                       strat_chem_get_extra_h2o
use           mo_chem_utls_mod, only : get_spc_ndx, get_rxt_ndx, &
                                       get_tracer_ndx, &
                                       get_solar_flux_by_band, &
                                       NO_TRACER
use        gfdl_astronomy_mod,  only : astronomy_init, astronomy_end, &
                                       diurnal_solar, universal_time
use gfdl_sat_vapor_pres_mod, only : sat_vapor_pres_init
use atmos_radon_mod,       only : atmos_radon_sourcesink,   &
                                  atmos_radon_init,         &
                                  atmos_radon_end
use atmos_regional_tracer_driver_mod,only : regional_tracer_driver, &
                                            regional_tracer_driver_init
use atmos_age_tracer_mod,  only : atmos_age_tracer_init, atmos_age_tracer, &
                                  atmos_age_tracer_end
use cloud_chem, only: CLOUD_CHEM_PH_LEGACY, CLOUD_CHEM_PH_BISECTION, &
                      CLOUD_CHEM_PH_CUBIC, CLOUD_CHEM_F1P,&
                      CLOUD_CHEM_F1P_BUG, CLOUD_CHEM_F1P_BUG2, CLOUD_CHEM_LEGACY
use aerosol_thermodynamics, only: AERO_ISORROPIA, AERO_LEGACY, NO_AERO
use mo_usrrxt_mod, only: HET_CHEM_LEGACY, HET_CHEM_J1M
use mo_errmsg,             only : errmsg

implicit none

private

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public  gfdl_am4chem_driver, gfdl_am4chem_init,  &
        gfdl_am4chem_time_vary, gfdl_am4chem_endts, gfdl_am4chem_end

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!     ...  declare type that will store the field infomation for the
!          emission file
!-----------------------------------------------------------------------
type,public :: field_init_type
   character(len=64), pointer :: field_names(:)
end type field_init_type

!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
integer, parameter :: maxinv = 100
real               :: relaxed_dt = SECONDS_PER_DAY*10.,     & ! relaxation timescale (sec) for the upper boundary values
                      relaxed_dt_lbc = SECONDS_PER_DAY*10., & ! relaxation timescale (sec) for the lower boundary values
                      ub_pres = 100.e2,               & ! pressure (Pa) above which to apply chemical upper boundary conditions
                      lb_pres = 950.e2                  ! pressure (Pa) below which to apply chemical lower boundary conditions
character(len=64)  :: file_sulfate = 'sulfate.nc',    & ! NetCDF file for sulfate concentrations
                      file_conc = 'conc_all.nc',      & ! NetCDF file for tracer concentrations (initial and fixed)
                      file_ub = 'ub_vals.nc'            ! NetCDF file for chemical upper boundary conditions
character(len=64)  :: file_dry = 'depvel.nc',         & ! NetCDF file for dry deposition velocities
                      file_jval_lut = 'jvals.v5',     & ! ascii file for photolysis rate lookup table
                      file_jval_lut_min = ''            ! ascii file for photolysis rate LUT (for solar min)
character(len=10), dimension(maxinv) :: inv_list =''    ! list of invariant (fixed) tracers
!real               :: lght_no_prd_factor = 1.           ! lightning NOx scale factor
!logical            :: normalize_lght_no_prd_area = .false. ! normalize lightning NOx production by grid cell area
!real               :: min_land_frac_lght = -999.        ! minimum land fraction for lightning NOx calculation
real               :: strat_chem_age_factor = 1.        ! scale factor for age of air
real               :: strat_chem_dclydt_factor = 1.     ! scale factor for dcly/dt
logical            :: do_tropchem = .false.             ! Do tropospheric chemistry?
logical            :: use_tdep_jvals = .false.          ! Use explicit temperature dependence for photolysis rates
real               :: o3_column_top = 10.               ! O3 column above model top (DU)
real               :: jno_scale_factor = 1.             ! scale factor for NO photolysis rate (jNO)
logical            :: repartition_water_tracers = .false. ! Allow PSC scheme to act on total water (vapor+condensed)
logical            :: allow_negative_cosz = .false.     ! Allow negative values for cosine of solar zenith angle
logical            :: allow_psc_settling_type1 = .false.! Allow Type-I (NAT) PSCs to settle
logical            :: allow_psc_settling_type2 = .false.! Allow Type-II (ice) PSCs to settle
logical            :: force_cly_conservation = .false.  ! Force chemical conservation of Cly
logical            :: rescale_cly_components = .false.  ! Rescale individual Cly components to total Cly VMR
logical            :: set_min_h2o_strat = .false.       ! Don't allow total water concentration in the stratosphere to fall below 2*CH4_trop
character(len=64)  :: ch4_filename = 'ch4_gblannualdata'! Methane timeseries filename
real               :: ch4_scale_factor = 1.             ! Methane scale factor to convert to VMR (mol/mol)
character(len=64)  :: cfc_lbc_filename = 'chemlbf'      ! Input file for CFC lower boundary conditions
logical            :: time_varying_cfc_lbc = .true.     ! Allow time variation of CFC lower boundary conditions
integer, dimension(6) :: cfc_lbc_dataset_entry = (/ 1, 1, 1, 0, 0, 0 /) ! Entry date for CFC lower boundary condition file
integer            :: verbose = 3                       ! level of diagnostic output
logical            :: retain_cm3_bugs = .false.         ! retain bugs present in code used in CM3
logical            :: do_fastjx_photo = .false.         ! use fastjx routine ?
character(len=32)  :: clouds_in_fastjx = 'lsc_only'     ! nature of clouds seen in fastjx calculation; may currently be 'none' or 'lsc_only' (default)
logical            :: check_convergence = .false.       ! if T, non-converged chem tendencies will not be used
real               :: e90_tropopause_vmr = 9.e-8        ! e90 tropopause concentration
logical            :: time_varying_solarflux = .false.  ! allow sloar cycle on fastjx v7.1

! namelist to fix solar flux bug
! if set to true then solar flux will vary with time
logical :: solar_flux_bugfix = .false.  ! reproduce original behavior


!JianHe: maybe we can pass co2 from rad to here
!co2
real*8             :: co2_fixed_value   = 400. ! JianHe
!character(len=64)  :: co2_filename = 'co2_gblannualdata'
character(len=64)  :: co2_filename = 'no_file'
real               :: co2_scale_factor = 1.e-6  ! ppm input
real               :: co2_fixed_year   = -999

character(len=64)  :: cloud_chem_pH_solver     = 'bisection'
character(len=64)  :: cloud_chem_type          = 'f1p_bug2'
logical            :: het_chem_fine_aerosol_only = .false.
real               :: min_lwc_for_cloud_chem     = 1.e-8
real               :: frac_dust_incloud          = 0
real               :: frac_aerosol_incloud       = 1
real               :: cloud_pH                   = -999  !<0 do not force
real               :: max_rh_aerosol             = 9999      !max rh used for aerosol thermo (to make sure no filter)
logical            :: limit_no3                  = .true.   !for isorropia/stratosphere

character(len=64)  :: aerosol_thermo_method = 'legacy'               ! other choice isorropia
character(len=64)  :: het_chem_type         = 'legacy'
real               :: gN2O5                 = 0.1
real               :: gNO2                  = 1e-4
real               :: gSO2                  = 0.
real               :: gSO2_dust             = 0.
real               :: gNH3                  = 0.05
real               :: gHNO3_dust            = 0.
real               :: gNO3_dust             = -999.
real               :: gN2O5_dust            = -999.
integer            :: gHNO3_dust_dynamic    = 0
real               :: gH2SO4_dust           = 0.
logical            :: do_h2so4_nucleation   = .false.
logical            :: cloud_ho2_h2o2        = .true.
real               :: gNO3                  = 0.1
real               :: gHO2                  = 1.

character(len=128) :: sim_data_filename = 'sim.dat'      ! Input file for chemistry pre-processor

character(len=64)  :: gso2_dynamic          = 'none'

type(tropchem_diag),  save :: trop_diag
type(tropchem_opt),   save :: trop_option

namelist /gfdl_am4chem_nml/    &
    relaxed_dt, relaxed_dt_lbc, ub_pres, lb_pres, file_sulfate, file_conc, file_ub, & 
    file_dry, inv_list, &
    strat_chem_age_factor, strat_chem_dclydt_factor, do_tropchem, use_tdep_jvals, file_jval_lut, &
    file_jval_lut_min, o3_column_top, jno_scale_factor, repartition_water_tracers, allow_negative_cosz, &
    allow_psc_settling_type1, allow_psc_settling_type2, force_cly_conservation, rescale_cly_components, & 
    set_min_h2o_strat, ch4_filename, ch4_scale_factor, co2_fixed_value, co2_fixed_year, co2_filename, &
    co2_scale_factor, cfc_lbc_filename, time_varying_cfc_lbc, cfc_lbc_dataset_entry, verbose, &
    retain_cm3_bugs, do_fastjx_photo, clouds_in_fastjx, check_convergence, solar_flux_bugfix, &
    e90_tropopause_vmr, aerosol_thermo_method, het_chem_type, gn2o5,gno2,gno3,gso2,gnh3, &
    ghno3_dust,gh2so4_dust,gho2,ghno3_dust_dynamic,gso2_dust,gn2o5_dust,gno3_dust, do_h2so4_nucleation, &
    cloud_chem_pH_solver, cloud_chem_type, min_lwc_for_cloud_chem,het_chem_fine_aerosol_only, cloud_pH, &
    frac_dust_incloud, frac_aerosol_incloud, max_rh_aerosol, limit_no3, cloud_ho2_h2o2, sim_data_filename,&
    time_varying_solarflux, gso2_dynamic

public           &
    relaxed_dt, relaxed_dt_lbc, ub_pres, lb_pres, file_sulfate, file_conc, file_ub, &
    file_dry, inv_list, &
    strat_chem_age_factor, strat_chem_dclydt_factor, do_tropchem, use_tdep_jvals, file_jval_lut, &
    file_jval_lut_min, o3_column_top, jno_scale_factor, repartition_water_tracers, allow_negative_cosz, &
    allow_psc_settling_type1, allow_psc_settling_type2, force_cly_conservation, rescale_cly_components, & 
    set_min_h2o_strat, ch4_filename, ch4_scale_factor, co2_fixed_value, co2_fixed_year, co2_filename, &
    co2_scale_factor, cfc_lbc_filename, time_varying_cfc_lbc, cfc_lbc_dataset_entry, verbose, &
    retain_cm3_bugs, do_fastjx_photo, clouds_in_fastjx, check_convergence, solar_flux_bugfix, &
    e90_tropopause_vmr, aerosol_thermo_method, het_chem_type, gn2o5,gno2,gno3,gso2,gnh3, &
    ghno3_dust,gh2so4_dust,gho2,ghno3_dust_dynamic,gso2_dust,gn2o5_dust,gno3_dust, do_h2so4_nucleation, &
    cloud_chem_pH_solver, cloud_chem_type, min_lwc_for_cloud_chem,het_chem_fine_aerosol_only, cloud_pH, &
    frac_dust_incloud, frac_aerosol_incloud, max_rh_aerosol, limit_no3, cloud_ho2_h2o2, sim_data_filename,&
    time_varying_solarflux,gso2_dynamic

integer                     :: nco2 = 0
character(len=7), parameter :: module_name = 'tracers'
real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4,     & !conversion factor (cm2/m2)
                   twopi      = 2.*PI

real, parameter :: mw_so4     = 96e-3 !kg/mol

logical, dimension(pcnstm1) :: has_ubc = .false., &
                               has_lbc = .false., &
                               fixed_lbc_time = .false.
type(time_type), dimension(pcnstm1) :: lbc_entry
logical, dimension(pcnstm1) :: has_airc = .false.
character(len=64),dimension(pcnstm1) :: ub_names, airc_names
real, parameter :: small = 1.e-50
integer :: sphum_ndx=0, cl_ndx=0, clo_ndx=0, hcl_ndx=0, hocl_ndx=0, clono2_ndx=0, &
           cl2o2_ndx=0, cl2_ndx=0, clno2_ndx=0, br_ndx=0, bro_ndx=0, hbr_ndx=0, &
           hobr_ndx=0, brono2_ndx=0, brcl_ndx=0, &
           hno3_ndx=0, o3_ndx=0, &
           no_ndx=0, no2_ndx=0, no3_ndx=0, n_ndx=0, n2o5_ndx=0, ho2no2_ndx=0, &
           pan_ndx=0, onit_ndx=0, mpan_ndx=0, isopno3_ndx=0, onitr_ndx=0, &
           extinct_ndx=0, noy_ndx=0, cly_ndx=0, bry_ndx=0, ch4_ndx=0, &
           dms_ndx=0, so4_ndx=0, co_ndx=0, n2o_ndx=0, lch4_ndx=0, oh_ndx=0
!JianHe: For additional NOy component
integer :: nh4no3_ndx=0, isn1_ndx=0, ino2_ndx=0, isnooa_ndx=0, inpn_ndx=0, &
           isopnb_ndx=0, isopnbo2_ndx=0, macrn_ndx=0, macrno2_ndx=0, &
           mvkn_ndx=0, r4n1_ndx=0, r4n2_ndx=0

integer :: o3s_ndx=0
integer :: o3s_e90_ndx=0
integer :: e90_ndx=0
integer :: nSOA      =0
integer :: nOH       =0
integer :: nC4H10    =0
integer :: ncodirect =0
integer :: ne90 =0
integer :: nsulfate  =0
integer :: nISOP     =0
integer :: nage      =0
integer :: naoanh    =0

logical :: do_interactive_h2o = .false.         ! Include chemical sources/sinks of water vapor?
real, parameter :: solarflux_min = 1.09082, &   ! solar minimum flux (band 18) [W/m2]
                   solarflux_max = 1.14694      ! solar maximum flux (band 18) [W/m2]
integer :: num_solar_bands = 18 !JianHe: hardcoded based on esfsw_2015


logical :: prevent_flux_through_ice = .true., step_update_tracer = .true.
                               ! when true, tracers will only be fluxed
                               ! through the non-ice-covered portions of
                               ! ocean grid boxes


!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer :: id_sul, id_temp, id_dclydt, id_dbrydt, id_dclydt_chem, &
           id_psc_sat, id_psc_nat, id_psc_ice, id_volc_aer, &
           id_imp_slv_nonconv, id_srf_o3, id_coszen, id_h2o_chem
integer :: inqa, inql, inqi !index of the three water species(nqa, nql, nqi)
integer :: age_ndx ! index of age tracer
logical :: module_is_initialized=.false.
logical :: use_lsc_in_fastjx

integer :: jo2_ndx, jo1d_ndx, jno2_ndx

integer, dimension(pcnstm1) :: indices, id_prod, id_loss, id_chem_tend, &
                               id_ub, id_lb
!new diagnostics (f1p)
integer, dimension(pcnstm1) :: id_prod_mol, id_loss_mol
integer :: id_pso4_h2o2,id_pso4_o3,id_ghno3_d,id_phno3_d(5), id_phno3_g_d, id_pso4_d(5), &
           id_pso4_g_d, id_gso2, id_aerosol_pH, id_cloud_pH, id_cloud_pHw, id_cld_amt_chem

logical :: has_ts_avg = .true.   ! currently reading in from monthly mean files.
integer, dimension(phtcnt)  :: id_jval
integer, dimension(gascnt)  :: id_rate_const
integer :: id_prodox, id_lossox ! for production and loss of ox.(jmao,1/1/2011)

type :: lb_type
   real, dimension(:), pointer :: gas_value
   type(time_type), dimension(:), pointer :: gas_time
end type lb_type
type(lb_type), dimension(pcnstm1) :: lb

type :: co2_type
   logical                                :: use_fix_value
   real                                   :: fixed_value
   real,dimension(:), pointer             :: gas_value
   type(time_type), dimension(:), pointer :: gas_time
   logical                                :: use_fix_time
   type(time_type)                        :: fixed_entry
end type co2_type
type(co2_type) :: co2_t
!type(interpolate_type), save :: drydep_data_default
integer :: clock_id,ndiag

integer, dimension(:), pointer :: nradon
integer, dimension(:), pointer :: nconvect

!++van
!real, allocatable    :: nb_N_Ox(:)
!--van

!---- version number ---------------------------------------------------
character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
!-----------------------------------------------------------------------

contains


!#######################################################################

! <SUBROUTINE NAME="gfdl_am4chem_driver">
!   <OVERVIEW>
!     Tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates the sources and sinks of tracers
!     due to tropospheric chemistry. It is called from atmos_tracer_driver.
!   </DESCRIPTION>
!   <TEMPLATE>
!     call gfdl_am4chem_driver (lon, lat, land, ocn_flx_fraction, pwt, r, chem_dt,           &
!                           Time, phalf, pfull, t, is, ie, js, je, dt, &
!                           z_half, z_full, q, tsurf, albedo, coszen,  &
!                           area, w10m, flux_sw_down_vis_dir, flux_sw_down_vis_dif, &
!                           half_day, &
!                           Time_next, rdiag,  do_nh3_atm_ocean_exchange, kbot)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="land" TYPE="real" DIM="(:,:)">
!     Land fraction
!   </IN>
!   <IN NAME="ocn_flx_fraction" TYPE="real" DIM="(:,:)">
!     Fraction of cell through which ocean flux is allowed
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     Pressure weighting (air mass) for each layer (kg/m2)
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="time_type">
!     Model time
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model full levels (Pa)
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>
!   <IN NAME="ie, je" TYPE="integer">
!     Local domain end indices
!   </IN>
!   <IN NAME="dt" TYPE="real">
!     Model physics timestep (s)
!   </IN>
!   <IN NAME="z_half" TYPE="real" DIM="(:,:,:)">
!     Height at model half levels (m)
!   </IN>
!   <IN NAME="z_full" TYPE="real" DIM="(:,:,:)">
!     Height at model full levels (m)
!   </IN>
!   <IN NAME="q" TYPE="real" DIM="(:,:,:)">
!     Specific humidity (kg/kg)
!   </IN>
!   <IN NAME="tsurf" TYPE="real" DIM="(:,:)">
!     Surface temperature (K)
!   </IN>
!   <IN NAME="albedo" TYPE="real" DIM="(:,:)">
!     Surface albedo
!   </IN>
!   <IN NAME="coszen" TYPE="real" DIM="(:,:)">
!     Cosine of the solar zenith angle
!   </IN>
!   <IN NAME="area" TYPE="real" DIM="(:,:)">
!     Grid box area (m^2)
!   </IN>
!   <IN NAME="w10m" TYPE="real" DIM="(:,:)">
!     Windspeed at 10m (m/s)
!   </IN>
!   <IN NAME="half_day" TYPE="real" DIM="(:,:)">
!     Half-day length  (dimensionless; 0 to pi)
!   </IN>
!   <IN NAME="do_nh3_atm_ocean_exchange," TYPE="logical">
!     Allow interactive atm-ocn exchange of NH3?
!   </IN>
!   <OUT NAME="chem_dt" TYPE="real" DIM="(:,:,:,:)">
!     Tracer tendencies from tropospheric chemistry (VMR/s)
!   </OUT>
!   <INOUT NAME="rdiag" TYPE="real" DIM="(:,:,:,:)">
!     Diagnostic tracer mixing ratios (tropchem tracers in VMR),
!     updated on output
!   </INOUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

subroutine gfdl_am4chem_driver(me,master,Time_init,Time,Time_next,dt_time,tracer_names, &
                            lon, lat, land, ocn_flx_fraction, &
                            r, chem_dt, phalf, pfull, t, is, ie, js, je, dt, & 
                            z_half, z_full, q, cldfrac, dfdage_interp, tsurf, albedo, &
                            icoszen,solcon, area, w10m, rdiag, jvals_out, &
                            prodox,lossox,prodlch4,lossch4,prodoh,lossoh)


!-----------------------------------------------------------------------
   integer, intent (in) :: me
   integer, intent (in) :: master

   character(len=32),   intent(in)                :: tracer_names(:)
   real, intent(in),    dimension(:,:)            :: lon, lat  !radian
   real, intent(in),    dimension(:,:)            :: land    ! land fraction
   real, intent(in),    dimension(:,:)            :: ocn_flx_fraction ! grid box fraction over which DMS flux from ocean occurs
   real, intent(in),    dimension(:,:,:,:)        :: r       ! prognostic tracer array 
   real, intent(out),   dimension(:,:,:,:)        :: chem_dt ! tendency array for all met+chem tracers (prog+diag)
   type(time_type), intent(in)                    :: Time_init,Time, Time_next,dt_time
   integer, intent(in)                            :: is, ie, js, je
   real, intent(in),    dimension(:,:,:)          :: phalf,pfull,t
   real, intent(in)                               :: dt      ! timestep (s)
   real, intent(in),    dimension(:,:,:)          :: z_half  ! height in meters at half levels
   real, intent(in),    dimension(:,:,:)          :: z_full  ! height in meters at full levels
   real, intent(in),    dimension(:,:,:)          :: q       ! specific humidity at current time step (kg/kg)
   real, intent(in),    dimension(:,:,:)          :: cldfrac ! total cloud fration
   real, intent(in),    dimension(:,:,:,:)        :: dfdage_interp
   real, intent(in),    dimension(:,:)            :: tsurf   ! surface temperature (K)
   real, intent(in),    dimension(:,:)            :: albedo  ! surface albedo
   real, intent(in),    dimension(:,:)            :: icoszen  
   real, intent(in)                               :: solcon !, rrsun
   real, intent(in),    dimension(:,:)            :: area    ! grid box area (m^2)
   real, intent(in),    dimension(:,:)            :: w10m    ! wind speed at 10m (m/s)
   real, intent(inout), dimension(:,:,:,:)        :: rdiag   ! diagnostic tracer concentrations
   real, intent(inout), dimension(:,:,:,:)        :: jvals_out ! jo2,j1d,jno2
   real, intent(inout), dimension(:,:,:)          :: prodox,lossox,prodlch4,lossch4,prodoh,lossoh
!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
   real, dimension(size(r,1),size(r,2),size(r,3)) :: rtnd, pwt
   real, dimension(size(r,1),size(r,2),size(r,3)) :: sulfate_data
!  real, dimension(size(r,1),size(r,2),size(r,3)) :: ub_temp,rno
   real, dimension(size(r,1),size(r,2),size(r,3),maxinv) :: inv_data
   real, dimension(size(r,1),size(r,2),size(r,3)) :: age, cly0, cly, cly_ratio, &
                                                     bry, dclydt, dbrydt, noy, &
                                                     extinct, strat_aerosol
   real, dimension(size(r,1),size(r,2),size(r,3),3) :: psc_vmr_save, dpsc_vmr
   real, dimension(size(r,1),size(r,2)) :: tsfcair, flux_sw_down_vis
   integer :: i,j,k,n,kb,id,jd,kd,ninv,nt,ntp, idx, index1, index2
!  integer :: nno,nno2
   integer :: inv_index
   integer :: plonl, nnn
   logical :: used
   real :: scale_factor, frac, ico2
   real,  dimension(size(r,1),size(r,3)) :: pdel, h2so4, h2o_temp, qlocal, cloud_water, co2_2d, frac_liq, s_down
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1)  :: r_temp, r_in, r_ub
   real, dimension(size(r,1),size(r,2),size(r,3)) :: tend_tmp, extra_h2o
   real, dimension(pcnstm1) :: r_lb
   real, dimension(size(land,1), size(land,2)) :: oro ! 0 and 1 rep. of land
   real, dimension(size(r,1),size(r,2)) :: coszen_local, fracday_local, half_day_local
   real, dimension(size(r,1),size(r,2)) :: coszen
   real :: rrsun_local, rrsun
   integer :: year,month,day,hour,minute,second
   integer :: jday
   real, dimension(size(r,1),size(r,2),size(r,3),pcnstm1) :: prod, loss
!   real, dimension(size(r,1),size(r,2),size(r,3)):: prodox, lossox
   real, dimension(size(r,1),size(r,2),size(r,3),phtcnt) :: jvals
   real, dimension(size(r,1),size(r,2),size(r,3),gascnt) :: rate_constants
   real, dimension(size(r,1),size(r,2),size(r,3)) :: imp_slv_nonconv
   real, dimension(size(r,1),size(r,2),size(r,3)):: e90_vmr
   real, dimension(size(r,1),size(r,2),size(r,3)):: dz
   real, dimension(size(r,1),size(r,2),size(r,3),size(r,4)):: reg_dtend
   real :: solar_phase
   real :: solflxband(num_solar_bands)
   type(psc_type) :: psc
   type(time_type) :: lbc_Time
   !f1p
   !trop diag arrays
   real, dimension(size(r,1),size(r,2),size(r,3),trop_diag%nb_diag) :: trop_diag_array
   type(time_type) :: co2_time
   character(len=32) :: tracer_name, noytracer

!-----------------------------------------------------------------------

!<ERROR MSG="gfdl_am4chem_init must be called first." STATUS=".true.">
!   Tropchem_driver_init needs to be called before tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call errmsg ('Tropchem_driver','gfdl_am4chem_init must be called first.', .true.)

   ntp = size(r,4)
   plonl = size(r,1)

!initialize diagnostic array
   trop_diag_array(:,:,:,:) = 0.
   where(land(:,:) >= 0.5)
      oro(:,:) = 1.
   elsewhere
      oro(:,:) = 0.
   endwhere

   id=size(r,1); jd=size(r,2); kd=size(r,3)

   ninv=0
   do n = 1, size(inv_list)
      if(inv_list(n) /= '') then
         ninv = ninv + 1
      else
         exit
      end if
   end do

   tsfcair(:,:) = t(:,:,kd)
   dz(:,:,:) = z_half(:,:,1:kd) - z_half(:,:,2:(kd+1))

   !JianHe: instead of reading offline so4, we use simulated so4
   sulfate_data(:,:,:) = r(:,:,:,indices(so4_ndx))  ! vmr

!------------------------------------------------------------------------
! Get air mass in layer (in kg/m2), equal to dP/g
! Get air density in layer (in kg/m3), equal to dP/(g*dz)
! JianHe: for non-hydrostatic, we do not use below eqn to calculate rho
!------------------------------------------------------------------------
      do k=1,kd
         pwt(:,:,k)=(phalf(:,:,k+1)-phalf(:,:,k))/grav
         !rho(:,:,k) = pwt(:,:,k)/(z_half(:,:,k) - z_half(:,:,k+1))
      enddo

!-----------------------------------------------------------------------
!--------- Get current date
!-----------------------------------------------------------------------
!      call get_date(Time, year, month, day, hour, minute, second)

!-----------------------------------------------------------------------
!     ... read in the concentrations of "invariant" (i.e., prescribed)
!         species
!-----------------------------------------------------------------------
!JianHe: no inv_list
!   do n = 1,ninv
!      call interpolator( conc, Time, phalf, inv_data(:,:,:,n), &
!                         trim(inv_list(n)), is, js)
!      inv_index = get_tracer_ndx( tracer_names, trim(inv_list(n)) ) - ntp
!      rdiag(:,:,:,inv_index) = inv_data(:,:,:,n)
!   end do

   chem_dt(:,:,:,:) =0.

!-----------------------------------------------------------------------
!     ... assign concentrations of prognostic (r) and diagnostic (rdiag)
!         species to r_temp
!-----------------------------------------------------------------------
   do n = 1,pcnstm1
      if(indices(n) <= ntp) then
         !JianHe: indices(n) is based on the entire tracer array 
         r_temp(:,:,:,n) = r(:,:,:,indices(n))  
      else
         r_temp(:,:,:,n) = rdiag(:,:,:,indices(n)-ntp)
      end if
   end do

!------------------------------------------------------------------------
! Compute radon source-sink tendency
!------------------------------------------------------------------------
    do nnn = 1, size(nradon(:))
     if (nradon(nnn) > 0) then
       if (nradon(nnn) > ntp) call errmsg ('gfdl_am4chem_driver', &
                            'Number of tracers .lt. number for radon', .true.)
         call atmos_radon_sourcesink (lon,lat,land,pwt,r(:,:,:,nradon(nnn)),  &
                                 rtnd)
       chem_dt(:,:,:,nradon(nnn))=chem_dt(:,:,:,nradon(nnn))+rtnd(:,:,:)
    endif

   end do

!-----------------------------------------------------------------------
!     ... convert to H2O VMR
!-----------------------------------------------------------------------
   if (sphum_ndx > 0) then
      r_temp(:,:,:,sphum_ndx) = r_temp(:,:,:,sphum_ndx) * WTMAIR / WTMH2O
   end if

!-----------------------------------------------------------------------
!     ... convert volcanic aerosol extinction into aerosol surface area
!-----------------------------------------------------------------------
   if (extinct_ndx > 0 .and. extinct_ndx <= ntp) then
      extinct(:,:,:) = r(:,:,:,extinct_ndx)
   else if (extinct_ndx > ntp) then
      extinct(:,:,:) = rdiag(:,:,:,extinct_ndx-ntp)
   else
      extinct(:,:,:) = 0.
   end if
   call strat_chem_get_aerosol( extinct, strat_aerosol )

!-----------------------------------------------------------------------
!     ... get age of air
!-----------------------------------------------------------------------
   if(age_ndx > 0 .and. age_ndx <= ntp) then
!JianHe: include age sourcesink here:
      call atmos_age_tracer( lon, lat, pwt,  &
                             r(:,:,:,age_ndx),  &
                             rtnd)
      chem_dt(:,:,:,age_ndx)=chem_dt(:,:,:,age_ndx)+rtnd(:,:,:)

      age(:,:,:) = r(:,:,:,age_ndx)
   else
      age(:,:,:) = 0.
   end if

!-----------------------------------------------------------------------
!     ... Chemical families
!-----------------------------------------------------------------------
   cly0(:,:,:) = 0.
   if (cl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl_ndx)
   end if
   if (clo_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clo_ndx)
   end if
   if (hcl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,hcl_ndx)
   end if
   if (hocl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,hocl_ndx)
   end if
   if (clono2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (cl2o2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl2o2_ndx)*2
   end if
   if (cl2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,cl2_ndx)*2
   end if
   if (clno2_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,clno2_ndx)
   end if
   if (brcl_ndx>0) then
      cly0(:,:,:) = cly0(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if

!-----------------------------------------------------------------------
!     ... cosine of solar zenith angle
!-----------------------------------------------------------------------
!JianHe: coszen is very different from those calculated in ccpp/physcis
   if (allow_negative_cosz) then
      call diurnal_solar( lat, lon, Time, coszen_local, fracday_local, &
                          rrsun_local, dt_time, &
                          allow_negative_cosz=.true. )
   else
      call diurnal_solar( lat, lon, Time, coszen_local, fracday_local, &
                          rrsun_local, dt_time, &
                          half_day_out=half_day_local )
   end if
   
!JianHe: we use coszen calculated from physics routine or here????
   coszen(:,:) = coszen_local(:,:) 
   rrsun = rrsun_local
   !coszen(:,:) = icoszen(:,:)

   r_temp(:,:,:,:) = MAX(r_temp(:,:,:,:),small)

!set CO2
   if (nco2 == NO_TRACER) then
      if (co2_t%use_fix_value) then
         co2_2d(:,:) = co2_t%fixed_value
      else
         if (co2_t%use_fix_time) then
            co2_time = co2_t%fixed_entry
         else
            co2_time = Time
         end if
         call time_interp( co2_time, co2_t%gas_time(:), frac, index1, index2 )
         ico2 = co2_t%gas_value(index1) + &
              frac*( co2_t%gas_value(index2) - co2_t%gas_value(index1) )
         co2_2d(:,:) = ico2
      end if
   end if

   do j = 1,jd
      do k = 1,kd
         pdel(:,k) = phalf(:,j,k+1) - phalf(:,j,k)
      end do
      qlocal(:,:) = q(:,j,:)
      if (nco2>0) then
         co2_2d(:,:) = r(:,j,:,nco2)
      end if

!-----------------------------------------------------------------------
!     ... get stratospheric h2so4
!-----------------------------------------------------------------------
      call strat_chem_get_h2so4( pfull(:,j,:), age(:,j,:), h2so4 )

!-----------------------------------------------------------------------
!     ... compute PSC amounts
!-----------------------------------------------------------------------
      if (sphum_ndx>0) then
         h2o_temp(:,:) = r_temp(:,j,:,sphum_ndx)
      else
         h2o_temp(:,:) = qlocal(:,:) * WTMAIR/WTMH2O
      end if
      cloud_water(:,:) = MAX(r(:,j,:,inql)+r(:,j,:,inqi),0.)
      where ( cloud_water(:,:) .gt. 0 )
          frac_liq(:,:) =  MAX( r(:,j,:,inql) , 0.)/cloud_water
      elsewhere
      frac_liq(:,:) = 1.
      endwhere

      if (repartition_water_tracers) then
         h2o_temp(:,:) = h2o_temp(:,:) + cloud_water(:,:) * WTMAIR/WTMH2O
      end if
      if (set_min_h2o_strat) then
         call strat_chem_get_extra_h2o( h2o_temp, age(:,j,:), r_temp(:,j,:,ch4_ndx), Time, extra_h2o(:,j,:) )
         h2o_temp(:,:) = h2o_temp(:,:) + extra_h2o(:,j,:)
      end if

      call strat_chem_get_psc( t(:,j,:), pfull(:,j,:), &
                               r_temp(:,j,:,hno3_ndx), h2o_temp(:,:), &
                               h2so4, strat_aerosol(:,j,:), psc, psc_vmr_out=psc_vmr_save(:,j,:,:) )

      if (repartition_water_tracers) then
         cloud_water(:,:) = MAX(0.,cloud_water(:,:) - psc_vmr_save(:,j,:,3)*WTMH2O/WTMAIR) ! reduce cloud_water by amount of type-II PSC
         h2o_temp(:,:) = h2o_temp(:,:) - cloud_water(:,:) * WTMAIR/WTMH2O                  ! remaining water is present as vapor
      end if

      if (sphum_ndx>0) then
         r_temp(:,j,:,sphum_ndx) = h2o_temp(:,:)
      end if
      qlocal(:,:) = h2o_temp(:,:) * WTMH2O/WTMAIR
      r_in(:,j,:,:) = r_temp(:,j,:,:)

!NEED TO BE MERGED
!!!!!!!
!!!!!!!
!      where( cloud_liq(:,:)+cloud_ice(:,:) .gt. small )
!         s_down(:,:)    = cloud_water(:,:)/(cloud_liq(:,:)+cloud_ice(:,:))
!         cloud_liq(:,:) = cloud_liq(:,:)*s_down(:,:)
!         cloud_ice(:,:) = cloud_ice(:,:)*s_down(:,:)
!      end where
!!!!!!
!!!!!!

!-----------------------------------------------------------------------
!     ... get solar cycle phase (use radiation band #18)
!-----------------------------------------------------------------------
      if (solar_flux_bugfix) then
         call get_solar_flux_by_band(Time,solflxband)
      else
         call get_solar_flux_by_band(Time,solflxband, ref=.true.)
      endif
      solar_phase = solflxband(num_solar_bands)
      solar_phase = (solar_phase-solarflux_min)/(solarflux_max-solarflux_min)

      if (me==master) then
        print *, 'solar phase', solflxband(num_solar_bands), solar_phase
        print *, 'solar data', solcon, coszen(10,j), rrsun
      endif
!-----------------------------------------------------------------------
!     ... get e90 concentrations
!-----------------------------------------------------------------------

   e90_ndx = get_tracer_ndx( tracer_names,'e90' )
   if (e90_ndx > 0) then
      e90_vmr(:,j,:) = r(:,j,:,e90_ndx)
   else
      e90_vmr(:,j,:) = 0.
   end if

!-----------------------------------------------------------------------
!     ... call chemistry driver
!-----------------------------------------------------------------------
      call chemdr(me, master, &
                  r_temp(:,j,:,:),             & ! species volume mixing ratios (VMR)
                  r(:,j,:,:),                  &
                  tracer_names,                &
                  phalf(:,j,:),                & ! pressure at boundaries (Pa)
                  pwt(:,j,:) ,                 & ! column air density (Kg/m2)
                  j,                           & ! j
                  Time_next,                   & ! time
                  lat(:,j),                    & ! latitude
                  lon(:,j),                    & ! longitude
                  dt,                          & ! timestep in seconds
                  phalf(:,j,SIZE(phalf,3)),    & ! surface press ( pascals )
                  phalf(:,j,1),                & ! model top pressure (pascals)
                  pfull(:,j,:),                & ! midpoint press ( pascals )
                  pdel,                        & ! delta press across midpoints
                  z_full(:,j,:),               & ! height at midpoints ( m )
                  z_half(:,j,:),               & ! height at interfaces ( m )
                  !MAX(r(:,j,:,inqa),0.),       & ! cloud fraction
                  MAX(cldfrac(:,j,:),0.),      & ! cloud fraction
                  cloud_water(:,:),            & ! total cloud water (kg/kg)
                  frac_liq(:,:),               & ! fraction of liquid water
                  t(:,j,:),                    & ! temperature
                  inv_data(:,j,:,:),           & ! invariant species
                  qlocal(:,:),                 & ! specific humidity ( kg/kg )
                  albedo(:,j),                 & ! surface albedo
                  coszen(:,j),                 & ! cosine of solar zenith angle
!                  rrsun,                       & ! earth-sun distance factor
!                  coszen_local(:,j),                 & ! cosine of solar zenith angle, wiggle
                  rrsun_local,                       & ! earth-sun distance factor
                  prod(:,j,:,:),               & ! chemical production rate
                  loss(:,j,:,:),               & ! chemical loss rate
                  jvals(:,j,:,:),              & ! photolysis rates (s^-1)
                  rate_constants(:,j,:,:),     & ! kinetic rxn rate constants (cm^3 molec^-1 s^-1 for 2nd order)
                  sulfate_data(:,j,:),         & ! sulfate aerosol
                  psc,                         & ! polar stratospheric clouds (PSCs)
                  do_interactive_h2o,          & ! include h2o sources/sinks?
                  solar_phase,                 & ! solar cycle phase (1=max, 0=min)
                  imp_slv_nonconv(:,j,:),      & ! flag for non-convergence of implicit solver
                  plonl,                       & ! number of longitudes
                  prodox(:,j,:),               & ! production of ox(jmao,1/1/2011)
                  lossox(:,j,:),               & ! loss of ox(jmao,1/1/2011)
                  e90_vmr(:,j,:),              & ! e90 concentrations
                  e90_tropopause_vmr,          & ! e90 tropopause threshold
                  co2_2d,                      &
                  trop_diag_array(:,j,:,:),    &
                  trop_option,                 &
                  trop_diag)

      !JianHe: diag out
      jvals_out(:,j,:,1) = jvals(:,j,:,jo2_ndx)
      jvals_out(:,j,:,2) = jvals(:,j,:,jo1d_ndx)
      jvals_out(:,j,:,3) = jvals(:,j,:,jno2_ndx)

      call strat_chem_destroy_psc( psc )
   end do

   r_temp(:,:,:,:) = MAX( r_temp(:,:,:,:), small )
   if (allow_psc_settling_type1 .or. allow_psc_settling_type2) then
      call strat_chem_psc_sediment( psc_vmr_save, pfull, dt, dpsc_vmr )
      if (.not. allow_psc_settling_type1) dpsc_vmr(:,:,:,2) = 0.
      if (.not. allow_psc_settling_type2) dpsc_vmr(:,:,:,3) = 0.
   end if

!-----------------------------------------------------------------------
!     ... output diagnostics
!-----------------------------------------------------------------------
   prodlch4 = prod(:,:,:,lch4_ndx)  ! vmr/s
   lossch4 = loss(:,:,:,ch4_ndx)  ! vmr/s

   prodoh = prod(:,:,:,oh_ndx)
   lossoh = loss(:,:,:,oh_ndx)

   do n = 1,pcnstm1

!-----------------------------------------------------------------------
!     ... compute tendency
!-----------------------------------------------------------------------
      tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - r_in(:,:,:,n) )/dt
      if(indices(n) <= ntp) then
!-----------------------------------------------------------------------
!     ... prognostic species
!-----------------------------------------------------------------------
!        tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - MAX(r(:,:,:,indices(n))*scale_factor,small) )/dt
!         chem_dt(:,:,:,indices(n)) = airc_emis(:,:,:,n) + emis_source(:,:,:,n) + tend_tmp(:,:,:)
!JianHe: we do not include emis here
         chem_dt(:,:,:,indices(n)) = tend_tmp(:,:,:)  
      else
!-----------------------------------------------------------------------
!     ... diagnostic species
!-----------------------------------------------------------------------
!        tend_tmp(:,:,:) = ( r_temp(:,:,:,n) - MAX(rdiag(:,:,:,indices(n)-ntp)*scale_factor,small) )/dt
         rdiag(:,:,:,indices(n)-ntp) = r_temp(:,:,:,n)
      end if
!-----------------------------------------------------------------------
!     ... output diagnostic tendency
!-----------------------------------------------------------------------
!      if(id_chem_tend(n)>0) then
!         used = send_data( id_chem_tend(n), tend_tmp(:,:,:), Time_next, is_in=is,js_in=js)
!      end if

!-----------------------------------------------------------------------
!     ... apply upper boundary condition
!-----------------------------------------------------------------------
      if(has_ubc(n)) then
         !call interpolator(ub(n), Time, phalf, r_ub(:,:,:,n), trim(ub_names(n)), is, js)
         !if(id_ub(n)>0) then
         !   used = send_data(id_ub(n), r_ub(:,:,:,n), Time_next, is_in=is, js_in=js)
         !end if
         where (pfull(:,:,:) < ub_pres)
            chem_dt(:,:,:,indices(n)) = (r_ub(:,:,:,n) - r(:,:,:,indices(n))) / relaxed_dt
         endwhere
      end if

!-----------------------------------------------------------------------
!     ... apply lower boundary condition
!-----------------------------------------------------------------------
      if(has_lbc(n)) then
         if (fixed_lbc_time(n)) then
            lbc_Time = lbc_entry(n)
         else
            lbc_Time = Time
         end if
         call time_interp( lbc_Time, lb(n)%gas_time(:), frac, index1, index2 )
         r_lb(n) = lb(n)%gas_value(index1) + frac*( lb(n)%gas_value(index2) - lb(n)%gas_value(index1) )
         !if(id_lb(n)>0) then
         !   used = send_data(id_lb(n), r_lb(n), Time_next)
         !end if
         if (me == master) then
            print *, 'lower bound n2o: ', r_lb(n)
            print *, 'model n2o: ', r(10,1,60,indices(n)), lat(10,1), lon(10,1)
         end if
         where (pfull(:,:,:) > lb_pres)
            chem_dt(:,:,:,indices(n)) = (r_lb(n) - r(:,:,:,indices(n))) / relaxed_dt_lbc
         endwhere
      end if

   end do

!-----------------------------------------------------------------------
!     ... ox budget diagnostics
!-----------------------------------------------------------------------
!   if(id_prodox>0) then
!      used = send_data(id_prodox, prodox(:,:,:), Time_next, is_in=is, js_in=js)
!   end if
!   if(id_lossox>0) then
!      used = send_data(id_lossox, lossox(:,:,:), Time_next, is_in=is, js_in=js)
!   end if

!-----------------------------------------------------------------------
!     ... Chemical families (Cly)
!-----------------------------------------------------------------------
   cly(:,:,:) = 0.
   if (cl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl_ndx)
   end if
   if (clo_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clo_ndx)
   end if
   if (hcl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,hcl_ndx)
   end if
   if (hocl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,hocl_ndx)
   end if
   if (clono2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clono2_ndx)
   end if
   if (cl2o2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl2o2_ndx)*2
   end if
   if (cl2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,cl2_ndx)*2
   end if
   if (clno2_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,clno2_ndx)
   end if
   if (brcl_ndx>0) then
      cly(:,:,:) = cly(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if

!-----------------------------------------------------------------------
!     ... Cly chemical tendency diagnostic
!-----------------------------------------------------------------------
!   if (id_dclydt_chem>0) then
!      used = send_data(id_dclydt_chem, (cly(:,:,:)-cly0(:,:,:))/dt, Time_next, is_in=is, js_in=js)
!   end if

!-----------------------------------------------------------------------
!     ... Cly conservation
!-----------------------------------------------------------------------
   if (force_cly_conservation .or. rescale_cly_components) then
      if (rescale_cly_components) then
         cly_ratio(:,:,:) = r(:,:,:,cly_ndx) / MAX( cly(:,:,:), small )
         cly(:,:,:) = r(:,:,:,cly_ndx)
      else if (force_cly_conservation) then
         cly_ratio(:,:,:) = cly0(:,:,:) / MAX( cly(:,:,:), small )
         cly(:,:,:) = cly0(:,:,:)
      end if
      if (cl_ndx>0) then
         r_temp(:,:,:,cl_ndx) = r_temp(:,:,:,cl_ndx) * cly_ratio(:,:,:)
      end if
      if (clo_ndx>0) then
         r_temp(:,:,:,clo_ndx) = r_temp(:,:,:,clo_ndx) * cly_ratio(:,:,:)
      end if
      if (hcl_ndx>0) then
         r_temp(:,:,:,hcl_ndx) = r_temp(:,:,:,hcl_ndx) * cly_ratio(:,:,:)
      end if
      if (hocl_ndx>0) then
         r_temp(:,:,:,hocl_ndx) = r_temp(:,:,:,hocl_ndx) * cly_ratio(:,:,:)
      end if
      if (clono2_ndx>0) then
         r_temp(:,:,:,clono2_ndx) = r_temp(:,:,:,clono2_ndx) * cly_ratio(:,:,:)
      end if
      if (cl2o2_ndx>0) then
         r_temp(:,:,:,cl2o2_ndx) = r_temp(:,:,:,cl2o2_ndx) * cly_ratio(:,:,:)
      end if
      if (cl2_ndx>0) then
         r_temp(:,:,:,cl2_ndx) = r_temp(:,:,:,cl2_ndx) * cly_ratio(:,:,:)
      end if
      if (clno2_ndx>0) then
         r_temp(:,:,:,clno2_ndx) = r_temp(:,:,:,clno2_ndx) * cly_ratio(:,:,:)
      end if
      if (brcl_ndx>0) then
         r_temp(:,:,:,brcl_ndx) = r_temp(:,:,:,brcl_ndx) * cly_ratio(:,:,:)
      end if
   end if

!-----------------------------------------------------------------------
!     ... Chemical families (Bry, NOy)
!-----------------------------------------------------------------------
   bry(:,:,:) = 0.
   if (br_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,br_ndx)
   end if
   if (bro_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,bro_ndx)
   end if
   if (hbr_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,hbr_ndx)
   end if
   if (hobr_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,hobr_ndx)
   end if
   if (brono2_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,brono2_ndx)
   end if
   if (brcl_ndx>0) then
      bry(:,:,:) = bry(:,:,:) + r_temp(:,:,:,brcl_ndx)
   end if
   
!++van
   noy(:,:,:) = 0.
! Loop over total number of atmospheric tracers (nt), not just solver tracers (pcnstm1)
!   call get_number_tracers(tracer_names, num_tracers=nt)
!   do n = 1,nt
!      if ( nb_N_Ox(n) .gt. 0.) then
!        tracer_name = tracer_name(n)
!        if (tracer_name .eq. 'brono2') then
!            noytracer = 'BrONO2' 
!        else if (tracer_name .eq. 'clono2') then
!            noytracer = 'ClONO2'
!        else
!            noytracer =  uppercase(tracer_name)
!        end if
!        idx = get_spc_ndx(noytracer)
!        noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,idx)*nb_N_ox(n)
!      end if
!   end do
!--van

! JianHe: NOy: no, no2, hno3, ho2no2, n2o5, pan, mpan, nh4no3, clono2, brono2,
! isn1, ino2, isnooa, inpn, isopnb, isopnbo2, macrn, macrno2, mvkn,
! r4n1, r4n2

   noy(:,:,:) = noy(:,:,:) + r_temp(:,:,:,no_ndx) + &
                r_temp(:,:,:,no2_ndx) + &
                r_temp(:,:,:,hno3_ndx) + &
                r_temp(:,:,:,ho2no2_ndx) + &
                r_temp(:,:,:,n2o5_ndx)*2. + &
                r_temp(:,:,:,pan_ndx) + &
                r_temp(:,:,:,mpan_ndx) + &
                r_temp(:,:,:,nh4no3_ndx) + &
                r_temp(:,:,:,clono2_ndx) + &
                r_temp(:,:,:,brono2_ndx) + &
                r_temp(:,:,:,isn1_ndx) + &
                r_temp(:,:,:,ino2_ndx) + &
                r_temp(:,:,:,isnooa_ndx) + &
                r_temp(:,:,:,inpn_ndx) + &
                r_temp(:,:,:,isopnb_ndx) + &
                r_temp(:,:,:,isopnbo2_ndx) + &
                r_temp(:,:,:,macrn_ndx) + &
                r_temp(:,:,:,macrno2_ndx) + &
                r_temp(:,:,:,mvkn_ndx) + &
                r_temp(:,:,:,r4n1_ndx) + &
                r_temp(:,:,:,r4n2_ndx)

!-----------------------------------------------------------------------
!     ... stratospheric Cly and Bry source
!-----------------------------------------------------------------------
   if(age_ndx > 0 .and. age_ndx <= ntp) then
      call strat_chem_dcly_dt(me,master,Time, phalf, is, js, age, cly, bry, dfdage_interp,dclydt, dbrydt)
      do k = 1,kd
         where( coszen(:,:) > 0. )
            dclydt(:,:,k) = 2*dclydt(:,:,k)
            dbrydt(:,:,k) = 2*dbrydt(:,:,k)
         elsewhere
            dclydt(:,:,k) = 0.
            dbrydt(:,:,k) = 0.
         end where
      end do
   else
      dclydt(:,:,:) = 0.
      dbrydt(:,:,:) = 0.
   end if
   
!   if(me==master) then
!     do i = 1,id
!       do j = 1,jd
!         write(*,*) 'solar data: ', i,j,coszen(i,j),coszen_local(i,j),rrsun,rrsun_local
!         !write(*,*) 'dclydt: ', dclydt(i,j,13), dbrydt(i,j,13)
!         !write(*,*) 'chem_dt: ', chem_dt(i,j,13,indices(cl_ndx)),chem_dt(i,j,13,indices(br_ndx))
!       enddo
!     enddo
!   endif

   if (cl_ndx>0) then
      chem_dt(:,:,:,indices(cl_ndx)) = chem_dt(:,:,:,indices(cl_ndx)) + dclydt(:,:,:)
      !used = send_data(id_dclydt, dclydt, Time_next, is_in=is, js_in=js)
   end if
   if (br_ndx>0) then
      chem_dt(:,:,:,indices(br_ndx)) = chem_dt(:,:,:,indices(br_ndx)) + dbrydt(:,:,:)
      !used = send_data(id_dbrydt, dbrydt, Time_next, is_in=is, js_in=js)
   end if

!-----------------------------------------------------------------------
!     ... Set diagnostic tracers for chemical families
!-----------------------------------------------------------------------
   if (noy_ndx > ntp) then
      rdiag(:,:,:,noy_ndx-ntp) = noy(:,:,:)
   end if
   if (cly_ndx > ntp) then
      rdiag(:,:,:,cly_ndx-ntp) = cly(:,:,:) + dclydt(:,:,:)*dt
   else if (cly_ndx > 0) then
      chem_dt(:,:,:,cly_ndx) = dclydt(:,:,:)
   end if
   if (bry_ndx > ntp) then
      rdiag(:,:,:,bry_ndx-ntp) = bry(:,:,:) + dbrydt(:,:,:)*dt
   end if

!-----------------------------------------------------------------------
!     ... convert H2O VMR tendency to specific humidity tendency
!-----------------------------------------------------------------------
   if (sphum_ndx > 0) then
      n = indices(sphum_ndx)
      !if (me == master) then
      ! print *, 'Check H2O index: ', sphum_ndx, n
      !end if
      chem_dt(:,:,:,n) = chem_dt(:,:,:,n) * WTMH2O / WTMAIR
!     chem_dt(:,:,:,n) = 0.
   end if

!   call mpp_clock_end(clock_id)

!------------------------------------------------------------------------
! Regional tracer driver
!------------------------------------------------------------------------
   if (ne90 > 0) then
      if (ne90 > ntp) call errmsg ('Tracer_driver', &
                         'Number of tracers < number for e90', .true.)
! r is array for all tracers, reg_dtend is tendency arrary for prognostic tracer 
      call regional_tracer_driver( lon, lat, pwt, r, reg_dtend, &
                                   phalf, is, js)
      chem_dt(:,:,:,1:ntp) = chem_dt(:,:,:,1:ntp) + reg_dtend(:,:,:,1:ntp)
   endif

!-----------------------------------------------------------------------


end subroutine gfdl_am4chem_driver
!</SUBROUTINE>

!#######################################################################

! <FUNCTION NAME="gfdl_am4chem_init">
!   <OVERVIEW>
!     Initializes the tropospheric chemistry driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the tropospheric chemistry module.
!     It is called from atmos_tracer_driver_init.
!     Data sets are read in for dry deposition, upper boundary conditions,
!     and emissions. Off-line sulfate concentrations are also read in for
!     use in calculating heterogeneous reaction rates (if SO4 is not
!     included as a tracer).
!   </DESCRIPTION>
!   <TEMPLATE>
!     Ltropchem = gfdl_am4chem_init( r, mask, axes, Time, &
!                                       lonb_mod, latb_mod, phalf, &
!                                       drydep_data )
!   </TEMPLATE>
!   <IN NAME="mask" TYPE="real, optional" DIM="(:,:,:)">
!      optional mask that designates which grid points
!      are above (1) or below (0) the ground
!   </IN>
!   <IN NAME="axes" TYPE="integer" DIM="(4)">
!     The axes relating to the tracer array
!   </IN>
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="lonb_mod" TYPE="real" DIM="(:,:)">
!     The longitude corners for the local domain.
!   </IN>
!   <IN NAME="latb_mod" TYPE="real" DIM="(:,:)">
!     The latitude corners for the local domain.
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <OUT NAME="drydep_data" TYPE="interpolate_type" DIM="(:)">
!     Tracer dry deposition velocities
!   </OUT>
!   <OUT NAME="Ltropchem" TYPE="logical">
!     Do tropospheric chemistry? (Output as function value)
!   </OUT>
!   <INOUT NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (tropchem tracers in VMR)
!   </INOUT>

 subroutine gfdl_am4chem_init (me, master, nlunit, input_nml_file, &
                             logunit, fn_nml, tracer_names, &
                             ntche, ndchm, Time_init, Time, restart)
!-----------------------------------------------------------------------
!
!   r    = tracer fields dimensioned as (nlon,nlat,nlev,ntrace)
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   implicit none

   integer, intent (in) :: me
   integer, intent (in) :: master
   type(time_type) :: Time_init,Time
   integer, intent (in) :: nlunit
   integer, intent (in) :: logunit
   integer, intent (in) :: ntche, ndchm
   character (len = 64), intent (in) :: fn_nml
   character (len = *),  intent (in) :: input_nml_file(:)
   character(len=32), intent(in) :: tracer_names(:)
   logical,          intent(in)  :: restart

   integer :: ios
   logical :: exists

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
      integer :: n, outunit, ix
!<f1p
      real    :: tracer_mw, sum_N_ox
      character(len=64), parameter    :: sub_name = 'gfdl_am4chem_init'

   real    :: small_value

   integer :: flag_file, flag_spec, flag_fixed
   integer :: i, nt, ntp, ntd
   integer :: ierr, io
   character(len=64) :: nc_file,filename,specname
   character(len=256) :: control=''
   character(len=64) :: name=''
   !type(interpolate_type) :: init_conc
   character(len=64),dimension(pcnstm1) :: &
                                           conc_files = '', &
                                           ub_files = '', &
                                           lb_files = '', &
                                           conc_names
   logical :: tracer_initialized

   integer :: unit
   character(len=16) :: fld
   character(len=32) :: tracer_name

   integer :: flb, series_length, year, diy
   real :: input_time
   real :: scale_factor, extra_seconds, fixed_year
   type(time_type) :: Year_t

!-----------------------------------------------------------------------
!     ... read namelist
!-----------------------------------------------------------------------

#ifdef INTERNAL_FILE_NML
   if (me == master) then
       write (*, *) "Read internal nml"
   endif
   read (input_nml_file, nml = gfdl_am4chem_nml)
#else
   inquire (file = trim (fn_nml), exist = exists)
   if (.not. exists) then
       write (*, *) 'gfdlam4_chem :: namelist file: ', trim (fn_nml), ' does not exist'
       stop
   else
       open (unit = nlunit, file = fn_nml, action = 'read' , status = 'old', iostat = ios)
   endif
   read (nlunit, nml = gfdl_am4chem_nml)
   close (nlunit)
#endif

   ! write version number and namelist to log file
   if (me == master) then
       write (*, *) " ================================================================== "
       write (*, *) "gfdl_am4chem_mod"
       write (*, nml = gfdl_am4chem_nml)
       verbose = verbose + 1
   endif

!----- initialize the number of tracers
    !JianHe:
    !nt, total # of tracers
    !ntp, total prog tracers (met+chem)
    !ntd, total diag tracers (chem only for now)
      nt = size(tracer_names)
      ntp = ntche   
      ntd = ndchm

    if (me == master) then
       write (*, *) "total # of tracers", nt
       write (*, *) "total # of prog tracers", ntp
       write (*, *) "total # of diag tracers", ntd
    endif
!-----------------------------------------------------------------------
!
!  When initializing additional tracers, the user needs to make changes
!
!-----------------------------------------------------------------------

      if (module_is_initialized) return

!-------------------------------------------------------------------------
!     ... Make sure input value for clouds_in_fastjx is a valid option.
!-------------------------------------------------------------------------
   if (trim(clouds_in_fastjx) == 'none') then
     use_lsc_in_fastjx = .false.
   else if (trim(clouds_in_fastjx) == 'lsc_only') then
     use_lsc_in_fastjx = .true.
   else
     call errmsg ('gfdl_am4chem_init', &
                     ' invalid string for clouds_in_fastjx', .true.)
   endif


if ( (gNO3_dust .gt. 0. .or. gN2O5_dust .gt. 0.) .and. .not. het_chem_fine_aerosol_only ) then
   call errmsg ('gfdl_am4chem_init', 'uptake on dust + all aerosol => double counting', .true. )
end if

if ( trim(cloud_chem_type) == 'legacy' ) then
   trop_option%cloud_chem = CLOUD_CHEM_LEGACY
   if(me == master) write(*,*) 'legacy_cloud'
elseif ( trim(cloud_chem_type) == 'f1p' ) then
   trop_option%cloud_chem = CLOUD_CHEM_F1P
elseif ( trim(cloud_chem_type) == 'f1p_bug' ) then
   trop_option%cloud_chem = CLOUD_CHEM_F1P_BUG
elseif ( trim(cloud_chem_type) == 'f1p_bug2' ) then
   trop_option%cloud_chem = CLOUD_CHEM_F1P_BUG2
end if

!cloud chem pH solver


if ( trim(cloud_chem_pH_solver) == 'am3' ) then
   trop_option%cloud_chem_pH_solver  = CLOUD_CHEM_PH_LEGACY
   if(me == master) write(*,*) 'legacypH'
elseif ( trim(cloud_chem_pH_solver) == 'bisection' ) then
   trop_option%cloud_chem_pH_solver   = CLOUD_CHEM_PH_BISECTION
elseif ( trim(cloud_chem_pH_solver) == 'cubic' ) then
   trop_option%cloud_chem_pH_solver   = CLOUD_CHEM_PH_CUBIC
else
   call errmsg ('gfdl_am4chem_init', 'undefined cloud chem', .true. )
end if

if ( cloud_pH .lt. 0 ) then
   trop_option%cloud_H = -999
else
   trop_option%cloud_H = 10**(-cloud_pH)
end if

   if(me == master) write(*,*) 'cloud_H',trop_option%cloud_H


if ( trim(het_chem_type) == 'legacy' ) then
   trop_option%het_chem = HET_CHEM_LEGACY
   if(me == master) write(*,*) 'Using legacy heterogeneous chemistry'
elseif ( trim(het_chem_type) == 'j1m' ) then
   trop_option%het_chem = HET_CHEM_J1M
   if(me == master) write(*,*) 'Using new (J1M) heterogeneous chemistry'
end if

! Chemical pre-processor input filename
trop_option%sim_data_flsp = sim_data_filename
if(me == master) write(*,*) 'Chemical pre-processor input filename: ', TRIM(sim_data_filename)

!gammas to be added when het chem is working
trop_option%gN2O5                    = gN2O5
if(me == master)    write(*,*)     "gN2O5:",trop_option%gN2O5
trop_option%gNO3                     = gNO3
if(me == master)    write(*,*)     "gNO3:",trop_option%gNO3
trop_option%gNO2                     = gNO2
if(me == master)    write(*,*)     "gNO2:",trop_option%gNO2
trop_option%gHO2                     = gHO2
if(me == master)    write(*,*)     "gHO2:",trop_option%gHO2
trop_option%gNH3                     = gNH3
if(me == master)    write(*,*)     "gNH3:",trop_option%gNH3
trop_option%retain_cm3_bugs = retain_cm3_bugs
trop_option%do_fastjx_photo = do_fastjx_photo
trop_option%min_lwc_for_cloud_chem = min_lwc_for_cloud_chem
trop_optioN%check_convergence = check_convergence
trop_optioN%use_lsc_in_fastjx = use_lsc_in_fastjx
trop_option%het_chem_fine_aerosol_only = het_chem_fine_aerosol_only
trop_option%cloud_ho2_h2o2 = cloud_ho2_h2o2
trop_option%max_rh_aerosol = max_rh_aerosol
trop_option%limit_no3      = limit_no3
trop_option%frac_aerosol_incloud = frac_aerosol_incloud
trop_option%time_varying_solarflux = time_varying_solarflux

trop_option%gSO2                     = gSO2
if(me == master) write(*,*) 'gSO2: ',trop_option%gSO2
if (trim(gso2_dynamic).eq.'none') then
   trop_option%gSO2_dynamic             = -1
else if (trim(gso2_dynamic).eq.'wang2014') then
   trop_option%gSO2_dynamic             = 1   
   !http://onlinelibrary.wiley.com/doi/10.1002/2013JD021426/full  
else if (trim(gso2_dynamic).eq.'zheng2015') then
   trop_option%gSO2_dynamic             = 2
!   http://www.atmos-chem-phys.net/15/2031/2015/
end if
if(me == master) write(*,*) 'gso2_dynamic case:',trop_option%gSO2_dynamic

!aerosol thermo
if    ( trim(aerosol_thermo_method)   == 'legacy' ) then
   trop_option%aerosol_thermo = AERO_LEGACY
   if(me == master) write(*,*) 'legacy no3'
elseif ( trim(aerosol_thermo_method)   == 'isorropia' ) then
   trop_option%aerosol_thermo = AERO_ISORROPIA
elseif ( trim(aerosol_thermo_method)   == 'no_thermo' ) then
   trop_option%aerosol_thermo = NO_AERO
else
   call errmsg ('gfdl_am4chem_init', 'undefined aerosol thermo', .true. )
end if

!---------------------------------------------------------------------
!  make sure that astronomy_mod has been initialized (if radiation
!  not being called in this run, it will not have previously been
!  initialized).
!---------------------------------------------------------------------
   if(me == master) write(*,*) 'Initialize astronomy'
   call astronomy_init

!   Initialize lookup tables needed for sat_vapor_pres 
  if(me == master) write(*,*) 'Initialize sat_vapor_pres'
    call sat_vapor_pres_init

!----- set initial value of radon ------------
  if(me == master) write(*,*) 'Initialize radon_tracer'
    call atmos_radon_init(me, master, tracer_names, nradon)

!----- initialize the age tracer ------------
   if(me == master) write(*,*) 'Initialize atmos_age_tracer'
   call atmos_age_tracer_init( tracer_names, nage)

!-----------------------------------------------------------------------
!     ... Initialize chemistry driver
!-----------------------------------------------------------------------

   if(me == master) write(*,*) 'Initialize chemistry driver ', trop_option

   call chemini( file_jval_lut, file_jval_lut_min, use_tdep_jvals, &
                 tracer_names, me, master, &
                 o3_column_top, jno_scale_factor, verbose,   &
                 retain_cm3_bugs, do_fastjx_photo, trop_option)

   if(me == master) write(*,*) 'Finish initializing chemistry driver'
!-----------------------------------------------------------------------
!     ... set initial value of indices
!-----------------------------------------------------------------------
   indices(:) = 0
   do i=1,pcnstm1
      !JianHe: this is based on the entire tracer array (met+chem)
      n = get_tracer_ndx(tracer_names, tracnam(i))
      if (trim(tracnam(i)) == 'H2O') then
         if (n <= 0) then
            if (me == master) then
              write(*,*) 'No H2O tracer find in arrary, use sphum instead'
            end if
            n = get_tracer_ndx(tracer_names, 'sphum')
         end if
         sphum_ndx = i
         do_interactive_h2o = .true.
      end if
      if (n >0) then
         indices(i) = n
         if (indices(i) > 0 .and. me == master) then
            write(*,30) tracnam(i),indices(i)
            write(logunit,30) trim(tracnam(i)),indices(i)
         end if
      else
!<ERROR MSG="Tropospheric chemistry tracer not found in field table" STATUS="WARNING">
!   A tropospheric chemistry tracer was not included in the field table
!</ERROR>
         call errmsg ('gfdl_am4chem_init', trim(tracnam(i)) // ' is not found', .true.)
      end if
   end do
30 format (A,' was initialized as tracer number ',i3)

   nco2 = get_tracer_ndx(tracer_names, 'co2' )
   if (nco2 > 0 .and. me == master) then
      call errmsg ('gfdl_am4chem_driver', 'CO2 is active',.true.)
   endif
   cl_ndx     = get_spc_ndx('Cl')
   clo_ndx    = get_spc_ndx('ClO')
   hcl_ndx    = get_spc_ndx('HCl')
   hocl_ndx   = get_spc_ndx('HOCl')
   clono2_ndx = get_spc_ndx('ClONO2')
   cl2o2_ndx  = get_spc_ndx('Cl2O2')
   cl2_ndx    = get_spc_ndx('Cl2')
   clno2_ndx  = get_spc_ndx('ClNO2')
   br_ndx     = get_spc_ndx('Br')
   bro_ndx    = get_spc_ndx('BrO')
   hbr_ndx    = get_spc_ndx('HBr')
   hobr_ndx   = get_spc_ndx('HOBr')
   brono2_ndx = get_spc_ndx('BrONO2')
   brcl_ndx   = get_spc_ndx('BrCl')
   hno3_ndx   = get_spc_ndx('HNO3')
   no_ndx     = get_spc_ndx('NO')
   no2_ndx    = get_spc_ndx('NO2')
   no3_ndx    = get_spc_ndx('NO3')
   n_ndx      = get_spc_ndx('N')
   n2o5_ndx   = get_spc_ndx('N2O5')
   ho2no2_ndx = get_spc_ndx('HO2NO2')
   pan_ndx    = get_spc_ndx('PAN')
   onit_ndx   = get_spc_ndx('ONIT')
   mpan_ndx   = get_spc_ndx('MPAN')
   isopno3_ndx= get_spc_ndx('ISOPNO3')
   onitr_ndx  = get_spc_ndx('ONITR')
   o3_ndx     = get_spc_ndx('O3')
   ch4_ndx    = get_spc_ndx('CH4')
   lch4_ndx   = get_spc_ndx('LCH4')
   dms_ndx    = get_spc_ndx('DMS')
   so4_ndx    = get_spc_ndx('SO4')
   co_ndx     = get_spc_ndx('CO')
   n2o_ndx    = get_spc_ndx('N2O')
   oh_ndx     = get_spc_ndx('OH')

   o3s_ndx     = get_spc_ndx('O3S')
   o3s_e90_ndx = get_spc_ndx('O3S_E90')
   e90_ndx     = get_tracer_ndx(tracer_names,'e90')
   extinct_ndx = get_tracer_ndx(tracer_names, 'Extinction')
   noy_ndx     = get_tracer_ndx(tracer_names, 'NOy')
   cly_ndx     = get_tracer_ndx(tracer_names, 'Cly')
   bry_ndx     = get_tracer_ndx(tracer_names, 'Bry')

!JianHe: For additional NOy component
   nh4no3_ndx  = get_spc_ndx('NH4NO3')
   isn1_ndx    = get_spc_ndx('ISN1')
   ino2_ndx    = get_spc_ndx('INO2')
   isnooa_ndx  = get_spc_ndx('ISNOOA')
   inpn_ndx    = get_spc_ndx('INPN')
   isopnb_ndx  = get_spc_ndx('ISOPNB')
   isopnbo2_ndx= get_spc_ndx('ISOPNBO2')
   macrn_ndx   = get_spc_ndx('MACRN')
   macrno2_ndx = get_spc_ndx('MACRNO2')
   mvkn_ndx    = get_spc_ndx('MVKN')
   r4n1_ndx    = get_spc_ndx('R4N1')
   r4n2_ndx    = get_spc_ndx('R4N2')

   jo2_ndx     = get_rxt_ndx( 'jo2' )
   jno2_ndx    = get_rxt_ndx( 'jno2' )
   jo1d_ndx    = get_rxt_ndx( 'jo1d' )

!-----------------------------------------------------------------------
!     ... Check Cly settings
!-----------------------------------------------------------------------
   if (rescale_cly_components) then
      !if (cly_ndx == NO_TRACER .or. .not. check_if_prognostic(tracer_names,cly_ndx)) then
      if (cly_ndx == NO_TRACER ) then
         call errmsg ('gfdl_am4chem_init', &
                          'rescale_cly_components=T requires Cly to be registered as a prognostic tracer', .true.)
      end if
      if (force_cly_conservation) then
         call errmsg ('gfdl_am4chem_init', &
                          'rescale_cly_components=T incompatible with force_cly_conservation=T setting', .true.)
      end if
   end if

   do i = 1,pcnstm1
!-----------------------------------------------------------------------
!     ... Upper boundary condition
!-----------------------------------------------------------------------
!      if( query_method('upper_bound', tracer_names,indices(i),name,control) ) then
!         if( trim(name)=='file' ) then
!            flag_file = parse(control, 'file',filename)
!            flag_spec = parse(control, 'name',specname)
!
!            if( flag_file > 0 .and. trim(filename) /= trim(file_ub) ) then
!               ub_files(i) = trim(filename)
!               call interpolator_init(ub(i), trim(filename), lonb_mod, latb_mod, &
!                       data_out_of_bounds=(/CONSTANT/),          &
!                       vert_interp=(/INTERP_WEIGHTED_P/))
!            else
!               ub_files(i) = trim(file_ub)
!               ub(i) = ub_default
!            end if
!            if(flag_spec > 0) then
!               ub_names(i) = trim(specname)
!            else
!               ub_names(i) = trim(lowercase(tracnam(i)))
!            end if
!
!            has_ubc(i) = .true.
!
!         end if
!      end if
 
!-----------------------------------------------------------------------
!     ... Lower boundary condition, for n2o....
!-----------------------------------------------------------------------
      lbc_entry(i) = Time

      flb = 66
      if ( i == n2o_ndx ) then
          has_lbc(i) = .true.
          lb_files(i) = 'INPUT/n2o_gblannualdata'
          scale_factor = 1.e-9  ! ppb to vmr
          inquire(file=trim(lb_files(i)),exist=exists)
          if (exists) then
            open(unit=flb,file=trim(lb_files(i)),action='read', &
                 status='old',iostat=ios)
                 read(flb, FMT='(i12)') series_length
                 allocate( lb(i)%gas_value(series_length), &
                           lb(i)%gas_time(series_length) )
!---------------------------------------------------------------------
!    convert the time stamps of the series to time_type variables.
!---------------------------------------------------------------------
                  do n = 1,series_length
                     read (flb, FMT = '(2f12.4)') input_time, lb(i)%gas_value(n)
                     year = INT(input_time)
                     Year_t = set_date(year,1,1,0,0,0)
                     diy = days_in_year (Year_t)
                     extra_seconds = (input_time - year)*diy*SECONDS_PER_DAY
                     lb(i)%gas_time(n) = Year_t + set_time(NINT(extra_seconds), 0)
                  end do
                  lb(i)%gas_value(:) = lb(i)%gas_value(:) * scale_factor
                  close( flb )

           else
               call errmsg ('gfdl_am4chem_init', &
                                   'Failed to find input file '//trim(lb_files(i)), .true.)
           end if
      end if

!-----------------------------------------------------------------------
!     ... Initial conditions
!-----------------------------------------------------------------------
!!      tracer_initialized = .false.
!      if (restart) then
!          tracer_initialized = .true.
!      else
!          !JianHe: this should be done in _prep
!      end if
!
!      if(.not. tracer_initialized) then
!         !JianHe: preprocessed ic
!        tracer_initialized = .true.
!      end if

   end do

!move CO2 input out of the loop of "do i = 1,pcnstm1", 2016-07-25
!fp
!CO2
         co2_t%use_fix_value  = .true.
         co2_t%fixed_value    = co2_fixed_value

!-----------------------------------------------------------------------
!     ... Print out settings for tracer
!-----------------------------------------------------------------------
   if( me == master ) then
      write(logunit,*) '---------------------------------------------------------------------------------------'
      do i = 1,pcnstm1
         write(logunit,*) 'The tracname index is ',i
         write(logunit,*) 'The tracname is ',tracnam(i)
         write(logunit,*) '---------------------------------------------------------------------------------------'
      end do
   end if


!-----------------------------------------------------------------------
!     ... Get the index number for the cloud variables
!-----------------------------------------------------------------------
   !JianHe: we do not have cld_amt in tracer
   !inqa = get_tracer_ndx(tracer_names,'cld_amt') ! cloud fraction
   inql = get_tracer_ndx(tracer_names,'liq_wat') ! cloud liquid specific humidity
   inqi = get_tracer_ndx(tracer_names,'ice_wat') ! cloud ice water specific humidity

   if( me == master ) then
     print *, 'cloud water indexes: ', inql, inqi
   end if

   age_ndx = get_tracer_ndx(tracer_names,'age')  ! age tracer

   nSOA      = get_tracer_ndx(tracer_names,'SOA')
   nOH       = get_tracer_ndx(tracer_names,'oh')
   nC4H10    = get_tracer_ndx(tracer_names,'c4h10')
   nISOP     = get_tracer_ndx(tracer_names,'isop')

! Check for presence of OH and C4H10 (diagnostic) tracers
! If not present set index to 1 so interface calls do not fail,
! but FATAL error will be issued by atmos_sulfate_init,
! atmos_carbon_aerosol_init (if do_dynamic_bc or do_dynamic_om), and
! atmos_SOA_init (if use_interactive_tracers)
      if (nOH == NO_TRACER) nOH = 1
      if (nC4H10 == NO_TRACER) nC4H10 = 1
      ncodirect = get_tracer_ndx(tracer_names,'codirect')
      ne90      = get_tracer_ndx(tracer_names,'e90')
      naoanh    = get_tracer_ndx(tracer_names,'aoanh')

!-----------------------------------------------------------------------
!     ... Initializations for stratospheric chemistry
!-----------------------------------------------------------------------
!++lwh
   if (set_min_h2o_strat) then
      if (ch4_ndx>0) then
         if (.not. has_lbc(ch4_ndx)) then
            call errmsg ('Tropchem_driver','set_min_h2o_strat=T, but LBC not set for CH4', .true.)
         end if
      else
         call errmsg ('Tropchem_driver','set_min_h2o_strat=T, but CH4 not included in chemistry solver', .true.)
      end if
   end if
   call strat_chem_utilities_init( Time_init, me, master, &
                                   strat_chem_age_factor, strat_chem_dclydt_factor, &
                                   set_min_h2o_strat, ch4_filename, ch4_scale_factor, &
                                   fixed_lbc_time(ch4_ndx), lbc_entry(ch4_ndx), &
                                   cfc_lbc_filename, time_varying_cfc_lbc, cfc_lbc_dataset_entry )
!--lwh

!-----------------------------------------------------------------------
!     ... initialize time_interp
!-----------------------------------------------------------------------
   !call time_interp_init


!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
!   clock_id = mpp_clock_id('Chemistry')
   call setsox_init(trop_option,tracer_names)
   call chemdr_init(trop_option,tracer_names)

!initialize diag array
   if ( trop_option%aerosol_thermo == AERO_ISORROPIA ) then
      small_value = 1.e-25
   else
      small_value = 1.e-20
   end if

   call tropchem_types_init(trop_diag,small_value)


! regional tracer driver
      if (ncodirect > 0 .or. ne90 > 0) then
        call regional_tracer_driver_init (me, master, tracer_names)
      endif


   module_is_initialized = .true.


!-----------------------------------------------------------------------

end subroutine gfdl_am4chem_init
!</FUNCTION>

!#####################################################################

subroutine gfdl_am4chem_time_vary (Time)

type(time_type), intent(in) :: Time

      integer :: yr, mo,day, hr,min, sec
      integer :: n

      call strat_chem_dcly_dt_time_vary (Time)

end subroutine gfdl_am4chem_time_vary




!#####################################################################

subroutine gfdl_am4chem_endts


      integer :: n


      call strat_chem_dcly_dt_endts



end subroutine gfdl_am4chem_endts


!######################################################################

subroutine gfdl_am4chem_end

!-----------------------------------------------------------------------
!     ... initialize mpp clock id
!-----------------------------------------------------------------------
   
!---------------------------------------------------------------------
!    verify that the module is initialized.
!---------------------------------------------------------------------
      if ( .not. module_is_initialized) then
        call errmsg ('gfdl_am4chem_end',  &
              'module has not been initialized',.true.)
      endif

      write (*,'(/,(a))') 'Exiting gfdl_am4chem_mod, have a nice day ...'

      call astronomy_end
      call atmos_radon_end
      call atmos_age_tracer_end

!   deallocate(nb_N_Ox)
   module_is_initialized = .false.
   
    
!-----------------------------------------------------------------------

end subroutine gfdl_am4chem_end

!############################################################################
end module gfdl_am4chem_mod
