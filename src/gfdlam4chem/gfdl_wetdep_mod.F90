module gfdl_wetdep_mod

  ! <CONTACT EMAIL="William.Cooke@noaa.gov">
  !   William Cooke
  ! </CONTACT>

  ! <REVIEWER EMAIL="Bruce.Wyman@noaa.gov">
  !   Bruce Wyman
  ! </REVIEWER>

  ! <HISTORY SRC="http://www.gfdl.noaa.gov/fms-cgi-bin/cvsweb.cgi/FMS/"/>

  ! <OVERVIEW>
  !     This code provides some utility routines for atmospheric tracers in the FMS framework.
  ! </OVERVIEW>
  ! <DESCRIPTION>
  !    This module gives utility routines which can be used to provide
  !    consistent removal mechanisms for atmospheric tracers.
  !
  !    In particular it provides schemes for wet and dry deposiiton that
  !    can be easily utilized.
  !
  ! </DESCRIPTION>


  ! --->h1g, add a scale factor for aerosol wet deposition, 2014-04-10
  use  catchem_constants, only : GRAV, &     ! acceleration due to gravity [m/s2]
       RDGAS, &    ! gas constant for dry air [J/kg/deg]
       vonkarm, &
       PI, &
       DENS_H2O, & ! Water density [kg/m3]
       WTMH2O, &   ! Water molecular weight [g/mole]
       WTMAIR, &   ! Air molecular weight [g/mole]
       AVOGNO      ! Avogadro's number
  use  gfdl_astronomy_mod,  only : astronomy_init, astronomy_end, &
                                   diurnal_solar, universal_time

  use  gfdl_time_utls_mod, only : time_type, &
                                       get_date, &
                                       set_date, &
                                       set_time, &
                                       days_in_year, &
                                       real_to_time_type, &
                                       time_type_to_real, &
                                       operator(+), operator(-), &
                                       time_interp, lowercase
  use  mo_chem_utls_mod, only : get_spc_ndx, get_rxt_ndx, &
                                       get_tracer_ndx, &
                                       get_solar_flux_by_band, &
                                       NO_TRACER

  use mo_errmsg,         only : errmsg

  implicit none
  private
  !-----------------------------------------------------------------------
  !----- interfaces -------

  public  gfdl_wet_deposition,    &
       gfdl_wetdep_end, &
       gfdl_wetdep_init

  !---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

  logical :: module_is_initialized = .FALSE.

  character(len=7), parameter :: mod_name = 'tracers'
  integer, parameter :: max_tracers = 300

  real, parameter :: T_homogeneous = 233.15
  !-----------------------------------------------------------------------
  !--- identification numbers for  diagnostic fields and axes ----
  integer :: &
       id_tracer_wdep_ls(max_tracers),   id_tracer_wdep_cv(max_tracers),  &
       id_tracer_wdep_lsin(max_tracers), id_tracer_wdep_cvin(max_tracers),&
       id_tracer_wdep_lsbc(max_tracers), id_tracer_wdep_cvbc(max_tracers),&
       id_tracer_reevap_ls(max_tracers), id_tracer_reevap_cv(max_tracers),&
       id_tracer_wdep_ls_3d(max_tracers)

  !cmip6 (f1p)
  !

  integer :: id_w10m, id_delm
  integer :: id_u_star, id_b_star, id_rough_mom, id_z_pbl,  &
       id_mo_length_inv, id_vds
  !character(len=32),  dimension(max_tracers) :: tracer_names     = ' '
  character(len=32),  dimension(max_tracers) :: tracer_units     = ' '
  character(len=128), dimension(max_tracers) :: tracer_longnames = ' '
  character(len=32),  dimension(max_tracers) :: tracer_wdep_names     = ' '
  character(len=32),  dimension(max_tracers) :: tracer_wdep_units     = ' '
  character(len=128), dimension(max_tracers) :: tracer_wdep_longnames = ' '
  !type(cmip_diag_id_type) :: ID_so2_reevap_ls
  !----------------parameter values for the diagnostic units--------------
  real, parameter :: mw_air = WTMAIR/1000.  ! Convert from [g/mole] to [kg/mole]
  real, parameter :: mw_h2o = WTMH2O/1000.  ! Convert from [g/mole] to [kg/mole]
  real, parameter :: mw_so4 = 96./1000.     ! Convert from [g/mole] to [kg/mole]
  real, parameter :: twopi = 2*PI

  type wetdep_type
     character (len=500) :: scheme, text_in_scheme, control
     real  :: Henry_constant
     real  :: Henry_variable
     real  :: frac_in_cloud
     real  :: frac_in_cloud_snow
     real  :: frac_in_cloud_snow_homogeneous
     real  :: alpha_r
     real  :: alpha_s
     logical :: Lwetdep, Lgas, Laerosol, Lice, so2_so4_evap, is_so2
  end type wetdep_type

  type(wetdep_type), dimension(:), allocatable :: Wetdep


  ! --->h1g, add a scale factor for aerosol wet deposition, 2014-04-10
  real ::                scale_aerosol_wetdep =1.0
  real ::                scale_aerosol_wetdep_snow =1.0
  character(len=64)  :: file_dry = 'depvel.nc'  ! NetCDF file for dry deposition velocities
  logical :: drydep_exp = .true.
  real :: T_snow_dep = 263.15
  real :: kbs_val   = 0.5 ! surface conductance of rough sea (m/s)
  !namelist /wetdep_nml/  scale_aerosol_wetdep,  scale_aerosol_wetdep_snow, file_dry, drydep_exp, T_snow_dep, &
  !                       kbs_val
  ! <---h1g,
contains

  !
  ! ######################################################################
  !
  !<SUBROUTINE NAME="gfdl_wetdep_init">
  !<OVERVIEW>
  ! This is a routine to create and register the dry and wet deposition
  ! fields of the tracers.
  !</OVERVIEW>
  !<DESCRIPTION>
  !  This routine creates diagnostic names for dry and wet deposition fields of the tracers.
  !  It takes the tracer name and appends "ddep" for the dry deposition field and "wdep" for
  !  the wet deposition field. This names can then be entered in the diag_table for
  !  diagnostic output of the tracer dry and wet deposition. The module name associated with
  !  these fields in "tracers". The units of the deposition fields are assumed to be kg/m2/s.
  !</DESCRIPTION>
  !<TEMPLATE>
  ! call gfdl_wetdep_init(lonb,latb, mass_axes, Time)
  !</TEMPLATE>
  !   <IN NAME="lonb" TYPE="real" DIM="(:,:)">
  !     The longitude corners for the local domain.
  !   </IN>
  !   <IN NAME="latb" TYPE="real" DIM="(:,:)">
  !     The latitude corners for the local domain.
  !   </IN>
  !   <IN NAME="mass_axes" TYPE="integer" DIM="(3)">
  !     The axes relating to the tracer array.
  !   </IN>
  !   <IN NAME="Time" TYPE="type(time_type)">
  !     Model time.
  !   </IN>

  subroutine gfdl_wetdep_init(me, master, tracer_names)

    ! Routine to initialize the tracer identification numbers.
    ! This registers the 2D fields for the wet and dry deposition.
    integer, intent (in) :: me
    integer, intent (in) :: master
    character(len=32), intent(in) :: tracer_names(:)

    integer :: ntrace
    character(len=20) :: units =''
    !
    integer :: n, logunit
    character(len=128) :: name

    logical  :: flag

    ! --->h1g, add a scale factor for aerosol wet deposition, 2014-04-10
    !   local variables:
    integer   :: unit, io, ierr
    ! <---h1g,

!    do n = 1, max_tracers
!       !write ( tracer_names(n),     100 ) n
!       write ( tracer_longnames(n), 102 ) n
!       tracer_units(n) = 'none'
!    enddo
!100 format ('tr',i3.3)
!102 format ('tracer ',i3.3)

    !call get_number_tracers(MODEL_ATMOS, num_tracers= ntrace)
    ntrace = size(tracer_names)

    if (ntrace > 0) then
       allocate (Wetdep(ntrace))
    endif
    do n = 1, ntrace
       !--- set tracer tendency names where tracer names have changed ---

       !call get_tracer_names(MODEL_ATMOS,n,tracer_names(n),tracer_longnames(n),tracer_units(n))
       name = tracer_names(n)
       tracer_longnames(n) = tracer_names(n)
       tracer_units(n)='vmr'

       if ((name=='sphum') .or. (name=='liq_wat') .or. (name=='ice_wat') &
          .or. (name=='rainwat') .or. (name=='snowwat')  .or. (name=='graupel') &
          .or. (name=='ice_nc') .or. (name=='rain_nc') .or. (name=='o3mr') &
          .or. (name=='rainwat') .or. (name=='soa') .or. (name=='sulf') &
          .or. (name=='pp25') .or. (name=='bc1') .or. (name=='bc2') &
          .or. (name=='oc1') .or. (name=='oc2') .or. (name=='dust1') &
          .or. (name=='dust2') .or. (name=='dust3') .or. (name=='dust4') &
          .or. (name=='dust5') .or. (name=='seas1') .or. (name=='seas2') &
          .or. (name=='seas3') .or. (name=='seas4') .or. (name=='seas5') &
          .or. (name=='pp10')) then
          tracer_units(n)='mmr'
       endif
       if (name=='cld_amt') then
          tracer_units(n)='1'
       endif
       if ((name=='age') .or. (name=='aoanh')) then
          tracer_units(n)='years'
       endif
       if (name=='lch4') then
          tracer_units(n)='ppm/s'
       endif

!       write (name,100) n
!       if (trim(tracer_names(n)) /= name) then
!          tracer_wdep_names(n) = trim(tracer_names(n)) //'_wdep'
!       endif
!       write (name,102) n
!       if (trim(tracer_longnames(n)) /= name) then
!          tracer_wdep_longnames(n) = &
!               trim(tracer_longnames(n)) // ' wet deposition for tracers'
!       endif

!       select case (trim(tracer_units(n)))
!       case ('mmr')
!          units = 'kg/m2/s'
!       case ('kg/kg')
!          units = 'kg/m2/s'
!       case ('vmr')
!          units = 'mole/m2/s'
!       case ('mol/mol')
!          units = 'mole/m2/s'
!       case ('mole/mole')
!          units = 'mole/m2/s'
!       case default
!          units = trim(tracer_units(n))//' kg/(m2 s)'
!          call errmsg('gfdl_wetdep_init',&
!               ' Wet dep units set to '//trim(units)//' in gfdl_wetdep for '//trim(tracer_names(n)),&
!               .true.)
!       end select

       !-----------------------------------------------------------------------
       !    read namelist.
       !-----------------------------------------------------------------------
!JianHe: hardcoded for now
!#ifdef INTERNAL_FILE_NML
!       read (input_nml_file, nml=wetdep_nml, iostat=io)
!       ierr = check_nml_error(io,'wetdep_nml')
!#else
!       if ( file_exist('input.nml')) then
!          unit =  open_namelist_file ( )
!          ierr=1; do while (ierr /= 0)
!          read  (unit, nml=wetdep_nml, iostat=io, end=10)
!          ierr = check_nml_error(io,'wetdep_nml')
!       end do
!10     call close_file (unit)
!    endif
!#endif

    !JianHe: hardcoded here, but would be flexible to have in field table or
    !external file in the future.
       
 !Default
    Wetdep(n)%control = "henry"
    Wetdep(n)%scheme = "None"
    Wetdep(n)%Henry_constant = 0.
    Wetdep(n)%Henry_variable = 0.
    Wetdep(n)%frac_in_cloud = 0.
    Wetdep(n)%frac_in_cloud_snow = 0.
    Wetdep(n)%alpha_r = 0.
    Wetdep(n)%alpha_s = 0.
    Wetdep(n)%Lwetdep = .false.
    Wetdep(n)%Lgas = .false.
    Wetdep(n)%Laerosol = .false.
    Wetdep(n)%Lice = .true.
    Wetdep(n)%frac_in_cloud_snow_homogeneous = 0.
    Wetdep(n)%so2_so4_evap = .false.

    !aerosol
    if (name=='soa') then
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.3*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=0.
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if ((name=='so4') .or. (name=='sulf')) then
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.3*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=0.
      Wetdep(n)%frac_in_cloud_snow_homogeneous = 0.1
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if ((name=='dust1') .or. (name=='dust2') .or. (name=='dust3') &
       .or. (name=='dust4') .or. (name=='dust5')) then
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.2*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=0.05*scale_aerosol_wetdep_snow
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if (name=='seas1') then
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.2*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=Wetdep(n)%frac_in_cloud
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if (name=='seas2') then
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.3*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=Wetdep(n)%frac_in_cloud
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if ((name=='seas3') .or. (name=='seas4') .or. (name=='seas5')) then
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.5*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=Wetdep(n)%frac_in_cloud
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if ((name=='bc1') .or. (name=='oc1')) then  ! bcphob, omphob
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=Wetdep(n)%frac_in_cloud
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if (name=='bc2') then  ! bcphil
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.2*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=0.
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if (name=='oc2') then  ! omphil
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.3*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=0.
      Wetdep(n)%frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if ((name=='nh4no3') .or. (name=='nh4')) then  
      Wetdep(n)%scheme = "aerosol_below"
      Wetdep(n)%frac_in_cloud = 0.3*scale_aerosol_wetdep
      Wetdep(n)%frac_in_cloud_snow=0.*scale_aerosol_wetdep_snow
      Wetdep(n)%frac_in_cloud_snow_homogeneous = 0.1
      Wetdep(n)%alpha_r = 0.001*scale_aerosol_wetdep
      Wetdep(n)%alpha_s = 0.001*scale_aerosol_wetdep_snow
      Wetdep(n)%Laerosol = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    ! gas
    if (name=='ch2o') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=6.22e-2
      Wetdep(n)%Henry_variable=6460.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='hno3') then
      Wetdep(n)%scheme = "henry_below"
      Wetdep(n)%Henry_constant=3.19e6
      Wetdep(n)%Henry_variable=8700.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if (name=='ho2no2') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=1.97e-1
      Wetdep(n)%Henry_variable=0.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='ch3ooh') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=2.24e-3
      Wetdep(n)%Henry_variable=5610.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='h2o2') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=7.35e-1
      Wetdep(n)%Henry_variable=6620.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='ch3cho') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=1.38e-4
      Wetdep(n)%Henry_variable=5600.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='mvk') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=2.07e-4
      Wetdep(n)%Henry_variable=7800.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='macr') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=4.24e-5
      Wetdep(n)%Henry_variable=5300.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if ((name=='ch3oh') .or. (name=='c2h5oh')) then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=2.24e-3
      Wetdep(n)%Henry_variable=5610.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='glyald') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=4.05e-1
      Wetdep(n)%Henry_variable=4600.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='hyac') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=6.22e-2
      Wetdep(n)%Henry_variable=6460.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if ((name=='hcl') .or. (name=='hbr') & 
        .or. (name=='isopooh')) then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=3.19e6
      Wetdep(n)%Henry_variable=8700.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='so2') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=1.49e-2
      Wetdep(n)%Henry_variable=5080.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%so2_so4_evap = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='dms') then
      Wetdep(n)%scheme = "henry_below"
      Wetdep(n)%Henry_constant=4.74e-6
      Wetdep(n)%Henry_variable=3100.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
    end if

    if (name=='nh3') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=3.3d1
      Wetdep(n)%Henry_variable=4.1d3
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='atooh') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=1.14e-1
      Wetdep(n)%Henry_variable=6300.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if ((name=='glyx') .or. (name=='mgly')) then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=3.6e0
      Wetdep(n)%Henry_variable=7200.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if (name=='iepox') then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=8.3e-1
      Wetdep(n)%Henry_variable=7400.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if ((name=='isopnb') .or. (name=='macrn') &
       .or. (name=='mvkn')) then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=1.7e-1
      Wetdep(n)%Henry_variable=9200.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    if ((name=='r4n1') .or. (name=='r4n2')) then
      Wetdep(n)%scheme = "henry_below_noice"
      Wetdep(n)%Henry_constant=2.1e0
      Wetdep(n)%Henry_variable=8700.
      Wetdep(n)%Lgas = .true.
      Wetdep(n)%Lwetdep = .true.
      Wetdep(n)%Lice = .false.
    end if

    Wetdep(n)%text_in_scheme = Wetdep(n)%scheme
    if (name=='so2') then
      Wetdep(n)%text_in_scheme = 'henry_below_noice_so2'
    end if

!    if (me.eq.master) then
!       if (Wetdep(n)%Lwetdep) then
!          write(*,*) 'name: ',trim(tracer_names(n))
!          write(*,*) 'scheme: ',trim(Wetdep(n)%scheme)
!          write(*,*) 'H:',Wetdep(n)%Henry_constant,Wetdep(n)%Henry_variable
!          write(*,*) 'frac_in_cloud', Wetdep(n)%frac_in_cloud
!          write(*,*) 'frac_in_cloud_snow', Wetdep(n)%frac_in_cloud_snow
!          write(*,*) 'frac_in_cloud_snow_homogeneous', Wetdep(n)%frac_in_cloud_snow_homogeneous
!          write(*,*) 'so2_so4_evap', Wetdep(n)%so2_so4_evap
!          write(*,*) 'alpha_r,alpha_s', Wetdep(n)%alpha_r, Wetdep(n)%alpha_s
!       end if
!    end if


    if ( lowercase(trim(tracer_names(n))) .eq. "so2" .or. lowercase(trim(tracer_names(n))) .eq. "simpleso2" ) then
       Wetdep(n)%is_so2 = .true.
    else
       Wetdep(n)%is_so2 = .false.
    end if

 enddo

 module_is_initialized = .TRUE.

end subroutine gfdl_wetdep_init


!####################################################################

!<SUBROUTINE NAME = "wet_deposition">
!<TEMPLATE>
!CALL wet_deposition( n, T, pfull, phalf, zfull, zhalf, &
!                     rain, snow, qdt, cloud, rain3d, snow3d, &
!                     tracer, tracer_dt, Time, cloud_param, is, js, dt )
!</TEMPLATE>
subroutine gfdl_wet_deposition( me, master, tracer_names,&
    n, T, pfull, phalf, zfull, zhalf, &
    rain, snow, qdt, cloud, cloud_frac, &
    f_snow_berg, rain3d, snow3d, &
    tracer, tracer_dt, cloud_param, &
    is, js, dt, sum_wdep_out, so2_so4_out )
  !
  !<OVERVIEW>
  ! Routine to calculate the fraction of tracer removed by wet deposition
  !</OVERVIEW>
  !
  !<IN NAME="n" TYPE="integer">
  !   Tracer number
  !</IN>
  !<IN NAME="is, js" TYPE="integer">
  !   start indices for array (computational indices)
  !</IN>
  !<IN NAME="T" TYPE="real" DIM="(:,:,:)">
  !   Temperature
  !</IN>
  !<IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
  !   Full level pressure field (Pa)
  !</IN>
  !<IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
  !   Half level pressure field (Pa)
  !</IN>
  !<IN NAME="zfull" TYPE="real" DIM="(:,:,:)">
  !   Full level height field (m)
  !</IN>
  !<IN NAME="zhalf" TYPE="real" DIM="(:,:,:)">
  !   Half level height field (m)
  !</IN>
  !<IN NAME="rain" TYPE="real" DIM="(:,:)">
  !   Precipitation in the form of rain
  !</IN>
  !<IN NAME="snow" TYPE="real" DIM="(:,:)">
  !   Precipitation in the form of snow
  !</IN>
  !<IN NAME="qdt" TYPE="real" DIM="(:,:,:)">
  !   The tendency of the specific humidity (+ condenstate) due to the cloud parametrization (kg/kg/s)
  !</IN>
  !<IN NAME="cloud" TYPE="real" DIM="(:,:,:)">
  !   Cloud amount (liquid + ice) (kg/kg)
  !</IN>
  !<IN NAME="cloud_frac" TYPE="real" DIM="(:,:,:)">
  !   Cloud area fraction
  !</IN>
  !<IN NAME="rain3d" TYPE="real" DIM="(:,:,:)">
  !   Precipitation in the form of rain (kg/m2/s)
  !</IN>
  !<IN NAME="snow3d" TYPE="real" DIM="(:,:,:)">
  !   Precipitation in the form of snow (kg/m2/s)
  !</IN>
  !<IN NAME="tracer" TYPE="real" DIM="(:,:,:)">
  !   The tracer field
  !</IN>
  !<IN NAME="Time" TYPE="type(time_type)">
  !   The time structure for submitting wet deposition as a diagnostic
  !</IN>
  !<IN NAME="cloud_param" TYPE="character">
  !   Is this a convective (convect) or large scale (lscale) cloud parametrization?
  !</IN>
  !<IN NAME="dt" TYPE="real">
  !   The model timestep (in seconds)
  !</IN>
  !<OUT NAME="tracer_dt" TYPE="real" DIM="(:,:,:)">
  !   The tendency of the tracer field due to wet deposition
  !</OUT>
  !<DESCRIPTION>
  ! Schemes allowed here are
  !
  ! 1) Deposition removed in the same fractional amount as the modeled precipitation rate is to
  !    a standardized precipitation rate.
  !    Basically this scheme assumes that a fractional area of the gridbox is affected by
  !    precipitation and that this precipitation rate is due to a cloud of standardized cloud
  !    liquid water content. Removal is constant throughout the column where precipitation is occuring.
  !
  ! 2) Removal according to Henry's Law. This law states that the ratio of the concentation in
  !    cloud water and the partial pressure in the interstitial air is a constant. If tracer
  !    is in VMR, the units for Henry's constant are mole/L/Pa (normally it is mole/L/atm).
  !    Parameters for a large number of species can be found at
  !    http://www.mpch-mainz.mpg.de/~sander/res/henry.html
  !
  ! 3) Aerosol removal, using specified in-cloud tracer fraction

  ! 4) Similar as 3) with some lwh modifications
  !
  ! To utilize this section of code add one of the following lines as
  ! a method for the tracer of interest in the field table.
  !<PRE>
  ! "wet_deposition","henry","henry=XXX, dependence=YYY"
  !     where XXX is the Henry's constant for the tracer in question
  !       and YYY is the temperature dependence of the Henry's Law constant.
  !
  ! "wet_deposition","fraction","lslwc=XXX, convlwc=YYY"
  !     where XXX is the liquid water content of a standard large scale cloud
  !       and YYY is the liquid water content of a standard convective cloud.
  !</PRE>

  !</DESCRIPTION>

  !-----------------------------------------------------------------------
  !     ... dummy arguments
  !-----------------------------------------------------------------------
 integer, intent (in) :: me
 integer, intent (in) :: master
 character(len=32), intent(in) :: tracer_names(:)

 integer,          intent(in)                     :: n, is, js
 real,             intent(in),  dimension(:,:,:)  :: T, pfull,phalf, zfull, zhalf, qdt, cloud, tracer
 real,             intent(in),  dimension(:,:,:)  :: cloud_frac
 real,             intent(in),  dimension(:,:,:)  :: f_snow_berg
 ! snow production by Bergeron process
 real,             intent(in),  dimension(:,:)    :: rain, snow
 character(len=*), intent(in)                     :: cloud_param
 real,             intent(out), dimension(:,:,:)  :: tracer_dt
 real,             intent(in)                     :: dt
 real,             intent(in),  dimension(:,:,:)  :: rain3d, snow3d
 real,             intent(out),  dimension(:,:),   optional :: sum_wdep_out
 real,             intent(out),  dimension(:,:,:), optional :: so2_so4_out

 !-----------------------------------------------------------------------
 !     ... local variables
 !-----------------------------------------------------------------------
 real, dimension(size(T,1),size(T,2),size(pfull,3)) :: &
      Htemp, xliq, n_air, rho_air, pwt, zdel, precip3d, scav_fact3d, &
      precip3ds, precip3dr
 real, dimension(size(T,1),size(T,2)) :: &
      temp_factor, scav_factor, washout, sum_wdep, &
      w_h2o, K1, K2, beta, f_a,  scav_factor_s,  &
      wdep_in, wdep_bc, fluxr,fluxs, tracer_flux
 real, dimension(size(T,1),size(T,2),size(pfull,3)) :: &
      in_temp, bc_temp, dt_temp, reevap_fraction, reevap_diag
 integer, dimension(size(T,1),size(T,2)) :: &
      ktopcd, kendcd
 real, dimension(size(rain3d,1),size(rain3d,2),size(rain3d,3)) :: rainsnow3d

 integer :: i, j, k, kk, id, jd, kd, flaglw

 real, dimension(size(T,1),size(T,2),size(T,3)) :: conc

 real :: conc_rain, conc_rain_total, conc_sat

 real, parameter ::  DENS_SNOW = 500.    ! Snow density [kg/m3]
 real, parameter ::  RGAS      = 8.3143  ! ideal gas constant Pa m3/mol/K

 real    :: &
      Henry_constant, Henry_variable, &
      clwc, wash, premin, prenow, hwtop, &
      diag_scale

 real*8, parameter :: &
      inv298p15 = 1./298.15, &     ! 1/K
      kboltz = 1.38E-23,         & ! J/K
      rain_diam  = 1.89e-3,     &  ! mean diameter of rain drop (m)
      rain_vterm = 7.48,        &  ! rain drop terminal velocity (m/s)
      vk_air = 6.18e-6,         &  ! kinematic viscosity of air (m^2/s)
      d_g = 1.12e-5,            &  ! diffusive coefficient (m^2/s)
      geo_fac = 6.,             &  ! geometry factor (surface area/volume = geo_fac/diameter)
      cm3_2_m3 = 1.e-6             ! m3/cm3

 real :: &
      k_g,                       & ! mass transfer coefficient (m/s)
      stay,                      & ! fraction
      fall_time                    ! fall time through layer (s)

 real :: f_a0, scav_factor0, sa_drop0, fgas0
 real :: frac_in_cloud, frac_in_cloud_snow, frac_in_cloud_snow_homogeneous, frac_int, ph
 real , parameter :: &
      R_r = 0.001, &               ! radius of cloud-droplets for rain
      R_s = 0.001, &               ! radius of cloud-droplets for snow
      frac_int_gas = 1.0,   &
      frac_int_aerosol= 0.5

 real :: alpha_r, alpha_s

 logical :: &
      used, &
      Lwetdep, Lgas, Laerosol, Lice
 character(len=500) :: &
      tracer_name, control, scheme, units, &
      text_in_scheme

 !-----------------------------------------------------------------------

 ktopcd = 0
 kendcd = 0

 tracer_dt   = 0.
 wdep_in     = 0.
 wdep_bc     = 0.
 beta        = 0.
 reevap_fraction = 0.
 reevap_diag = 0.
 tracer_flux = 0.

 sum_wdep = 0.

 id = size(T,1)
 jd = size(T,2)
 kd = size(T,3)

 !call get_tracer_names(MODEL_ATMOS,n,tracer_name, units = units)
 tracer_name = tracer_names(n)
 units = tracer_units(n) 

    if (me.eq.master) then
      if (Wetdep(n)%Lwetdep) then
       if ((tracer_name=='so2') .or. (tracer_name=='nh3') .or. &
          (tracer_name=='so4'))  then
          write(*,*) 'name: ',trim(tracer_names(n))
          write(*,*) 'scheme: ',trim(Wetdep(n)%scheme)
          write(*,*) 'H:',Wetdep(n)%Henry_constant,Wetdep(n)%Henry_variable
          write(*,*) 'frac_in_cloud', Wetdep(n)%frac_in_cloud
          write(*,*) 'frac_in_cloud_snow', Wetdep(n)%frac_in_cloud_snow
          write(*,*) 'frac_in_cloud_snow_homogeneous', Wetdep(n)%frac_in_cloud_snow_homogeneous
          write(*,*) 'so2_so4_evap', Wetdep(n)%so2_so4_evap
          write(*,*) 'alpha_r,alpha_s', Wetdep(n)%alpha_r, Wetdep(n)%alpha_s
       end if
     end if
    end if

 if ( .not. Wetdep(n)%Lwetdep) return
 text_in_scheme = Wetdep(n)%text_in_scheme
 control = Wetdep(n)%control
 scheme = Wetdep(n)%scheme
 Henry_constant = Wetdep(n)%Henry_constant
 Henry_variable = Wetdep(n)%Henry_variable
 frac_in_cloud = Wetdep(n)%frac_in_cloud
 frac_in_cloud_snow = Wetdep(n)%frac_in_cloud_snow
 frac_in_cloud_snow_homogeneous = Wetdep(n)%frac_in_cloud_snow_homogeneous
 alpha_r = Wetdep(n)%alpha_r
 alpha_s = Wetdep(n)%alpha_s
 Lwetdep = Wetdep(n)%Lwetdep
 Lgas    = Wetdep(n)%Lgas
 Laerosol = Wetdep(n)%Laerosol
 Lice    = Wetdep(n)%Lice
 
 rho_air(:,:,:) = pfull(:,:,:) / ( T(:,:,:)*RDGAS ) ! kg/m3
 !   Lice = .not. (scheme=='henry_noice' .or. scheme=='henry_below_noice' .or. &
 !                 scheme=='aerosol_noice' .or. scheme=='aerosol_below_noice' )
 if (Lice) then 
    rainsnow3d(:,:,:) = rain3d(:,:,:) + snow3d(:,:,:)
 else
    rainsnow3d(:,:,:) = rain3d(:,:,:)
 end if
 do k=1,kd
    precip3d(:,:,k) = rainsnow3d(:,:,k+1)-rainsnow3d(:,:,k)
    precip3dr(:,:,k) = rain3d(:,:,k+1)-rain3d(:,:,k)
    pwt(:,:,k)  = ( phalf(:,:,k+1) - phalf(:,:,k) )/GRAV ! kg/m2
    zdel(:,:,k) = zhalf(:,:,k) - zhalf(:,:,k+1) ! m
    !zdel(:,:,k) = pwt(:,:,k)/rhoa(:,:,k)     ! JianHe: follow GOCART?
 end do
 if (Lice) then
    do k=1,kd
       precip3ds(:,:,k) = snow3d(:,:,k+1)-snow3d(:,:,k) 
    end do
 else
    do k=1,kd
       precip3ds(:,:,k) = 0.
    end do
 endif
 !
 !+++ pag:
 !
 if(lowercase(scheme)=='wdep_gas' .or. lowercase(scheme)=='wdep_aerosol') then
    ph = 5.0
    !  cloud liquid water content
    xliq(:,:,:) = 0.5E-3                           !default value
    if(trim(cloud_param) .eq. 'convect') then
       xliq(:,:,:) = 1.0e-3
    elseif(trim(cloud_param) .eq. 'lscale') then
       xliq(:,:,:) = 0.5e-3
    endif

    !
    ! tracer == gas
    !
    if(lowercase(scheme)=="wdep_gas") THEN
       frac_int = frac_int_gas
       do k = 1, kd
          do j = 1, jd
             do i = 1, id
                Htemp(i,j,k)=henry_constant &
                     *exp(-Henry_variable*(1./298.-1./T(i,j,k)))
                K1(i,j)=1.2e-2*exp(-2010*(1/298.-1/T(i,j,k)))
                K2(i,j)=6.6e-8*exp(-1510*(1/298.-1/T(i,j,k)))
                HTemp(i,j,k)=Htemp(i,j,k)*(1 + K1(i,j) &
                     /10.**(-ph) + K1(i,j)*K2(i,j)/(10.**(-ph))**2)
                f_a(i,j) = Htemp(i,j,k)/101.325*RGAS &
                     *T(i,j,k)*xliq(i,j,k)*rho_air(i,j,k)/DENS_H2O
                scav_fact3d(i,j,k)=f_a(i,j)/(1.+f_a(i,j))
             enddo
          enddo
       enddo
    elseif(lowercase(scheme)=="wdep_aerosol") THEN
       frac_int = frac_int_aerosol
       scav_fact3d(:,:,:)=frac_in_cloud
    else
       print *,' Aerosol number =',n,' tracer_name=',tracer_name,' scheme=',text_in_scheme
       print *, 'Please check "am2p12.ft'
       call errmsg('wet_deposition', 'Tracer is neither aerosol NOR gas.', .true. )
    endif
    !
    !in cloud scavenging
    !
    do k=1,kd
       do j = 1, jd
          do i = 1, id
             beta(i,j)  = MAX( 0.0, precip3d(i,j,k)/pwt(i,j,k)/xliq(i,j,k))
             in_temp(i,j,k) = (exp(-beta(i,j)*scav_fact3d(i,j,k)*dt)-1.0)
             if ( tracer(i,j,k) .gt. 0.) then
                wdep_in(i,j)=wdep_in(i,j) &
                     - in_temp(i,j,k)*tracer(i,j,k)*pwt(i,j,k)
                tracer_dt(i,j,k) = tracer_dt(i,j,k) &
                     - in_temp(i,j,k)*tracer(i,j,k)/dt
             endif
             !
             !--reevaporation
             !    calculation of fracion of aerosols to-be
             !    reevaporated to the atmosphere:

             beta(i,j)=precip3d(i,j,k)
             if (beta(i,j) < 0.) then
                beta(i,j) = beta(i,j)/rainsnow3d(i,j,k)
             endif
             if (rainsnow3d(i,j,k+1) == 0. ) then
                !--reevaporation total
                beta(i,j)=MIN(MAX(0.,-beta(i,j)),1.)
             else
                beta(i,j)=MIN(MAX(0.,-beta(i,j))*frac_int,1.)
             end if
             ! reevporating to atmosphere
             reevap_diag(i,j,k)=beta(i,j)*wdep_in(i,j)
             wdep_in(i,j) = wdep_in(i,j)*(1.-beta(i,j))
             tracer_dt(i,j,k) = tracer_dt(i,j,k) &
                  - reevap_diag(i,j,k)/pwt(i,j,k)/dt
          enddo
       enddo
    enddo
    ! Below cloud scavenging
    do k=1,kd
       do j = 1, jd
          do i = 1, id
             fluxs(i,j) = (snow3d(i,j,k+1)+snow3d(i,j,k))/2.0
             fluxr(i,j) = (rain3d(i,j,k+1)+rain3d(i,j,k))/2.0
             bc_temp(i,j,k) = 3./4.*dt* &
                  (fluxr(i,j)*alpha_r/R_r/DENS_H2O &
                  + fluxs(i,j)*alpha_s/R_s/DENS_SNOW)
             if ( tracer(i,j,k) .gt. 0. ) then
                wdep_bc(i,j)=wdep_bc(i,j) &
                     + bc_temp(i,j,k)*tracer(i,j,k)*pwt(i,j,k)
                tracer_dt(i,j,k) = tracer_dt(i,j,k) &
                     + bc_temp(i,j,k)*tracer(i,j,k)/dt
             endif
          enddo
       enddo
    enddo
    !
    !  end wdep_gas or wdep_aerosol
    !
 else

    ! Calculate fraction of precipitation reevaporated in layer
    do k=1,kd
       where( rainsnow3d(:,:,k) > 0. .and. precip3d(:,:,k) < 0. )
          reevap_fraction(:,:,k) = &
               -precip3d(:,:,k) / (rainsnow3d(:,:,k)) ! fraction
       end where
       ! Assume that the tracer reevaporation fraction is 50% of the precip
       ! reevaporation fraction, except when fraction = 100%
       !      where( reevap_fraction(:,:,k) < 1. )
       !         reevap_fraction(:,:,k) = 0.5*reevap_fraction(:,:,k)
       !      end where
    end do
    !  cloud liquid water content
    !   xliq = 0.5E-3                           !default value
    !   if(trim(cloud_param) .eq. 'convect') then
    !      xliq = 1.0e-3
    !   elseif(trim(cloud_param) .eq. 'lscale') then
    !      xliq = 0.5e-3
    !   endif

    ! Lgas = lowercase(scheme)=='henry' .or. lowercase(scheme)=='henry_below' .or. &
    !        lowercase(scheme)=='henry_noice' .or. lowercase(scheme)=='henry_below_noice'
    ! Laerosol = lowercase(scheme)=='aerosol' .or. lowercase(scheme)=='aerosol_below' .or. &
    !            lowercase(scheme)=='aerosol_noice' .or. lowercase(scheme)=='aerosol_below_noice'
    ! Assume that the aerosol reevaporation fraction is 50% of the precip
    ! reevaporation fraction, except when fraction = 100%
    if( Lgas ) then
       frac_int = frac_int_gas
    elseif( Laerosol ) then
       frac_int = frac_int_aerosol
    else
       frac_int = 1.
    end if

    if( Lgas .or. Laerosol ) then
       ! units = VMR
       !
       ! Henry_constant (mole/L/Pa) = [X](aq) / Px(g)
       ! where [X](aq) is the concentration of tracer X in precipitation (mole/L)
       !       Px(g) is the partial pressure of the tracer in the air (Pa)
       !
       ! VMR (total) = VMR (gas) + VMR (aq)
       !             = VMR (gas) + [X] * L
       !
       ! where L = cloud liquid amount (kg H2O/mole air)
       !
       ! Using Henry's Law, [X] = H * Px = H * VMR(gas) * Pfull
       !
       ! So, VMR (total) =  VMR(gas) * [ 1 + H * Pfull * L ]
       !
       ! VMR(gas) = VMR(total) / [1 + H * Pfull * L]
       !
       ! [X] = H * Pfull * VMR(total) / [ 1 + H * Pfull * L]
       !
       ! Following Giorgi and Chameides, JGR, 90(D5), 1985, the first-order loss
       ! rate constant (s^-1) of X due to wet deposition equals:
       !
       ! k = W_X / n_X
       !
       ! where W_x = the loss rate (molec/cm3/s), and n_X = the number density (molec/cm3)
       !
       ! W_X = [X] * W_H2O / (55 mole/L)
       ! n_x = VMR(total) * n_air (molec/cm3) = VMR(total) * P/(kT) * 1E-6 m3/cm3
       !
       ! where P = atmospheric pressure (Pa)
       !       k = Boltzmann's constant = 1.38E-23 J/K
       !       T = temperature (K)
       !       W_H2O = removal rate of water (molec/cm3/s)
       !
       !             [X] * W_H2O / 55
       ! So, k = ------------------------------
       !         VMR(total) * P/(kT) * 1E-6
       !
       !         W_H2O    H * VMR(total) * P / [ 1 + H * P *L ]
       !       = ----- * ---------------------------------------
       !          55          VMR(total) * P/(kT) * 1E-6
       !
       !         W_H2O     H * kT * 1E6
       !       = ----- *  -------------
       !          55      1 + H * P * L
       !
       !         W_H2O     1     1     H * P * L
       !       = ----- * ----- * - * -------------
       !          55     n_air   L   1 + H * P * L
       !
       ! where W_H2O = precip3d (kg/m2/s) * (AVOGNO/mw_h2o) (molec/kg) / zdel (m) * 1E-6 m3/cm3
       !
       if( (Lgas .and. Henry_constant > 0) .or. Laerosol ) then
          in_temp(:,:,:) = 0.
          bc_temp(:,:,:) = 0.
          do k=1,kd
             ! Calculate the temperature dependent Henry's Law constant
             scav_factor(:,:) = 0.0
             xliq(:,:,k)  = MAX( cloud(:,:,k) * mw_air, 0. ) ! (kg H2O)/(mole air)
             n_air(:,:,k) = pfull(:,:,k) / (kboltz*T(:,:,k)) * cm3_2_m3 ! molec/cm3
             if (Lgas) then
                temp_factor(:,:) = 1/T(:,:,k)-inv298p15
                Htemp(:,:,k) = Henry_constant * &
                     exp( Henry_variable*temp_factor )
                f_a(:,:) = Htemp(:,:,k) * pfull(:,:,k) * xliq(:,:,k) ! / cloud_frac
                scav_factor(:,:) = f_a(:,:) / ( 1.+f_a(:,:) )
                scav_factor_s(:,:) = scav_factor(:,:)
             else if (Laerosol) then
                scav_factor(:,:) = frac_in_cloud
                if (frac_in_cloud == frac_in_cloud_snow) then
                   scav_factor_s(:,:) = scav_factor(:,:)
                else
                   scav_factor_s(:,:) = f_snow_berg(:,:,k)*frac_in_cloud_snow +&
                        (1.-f_snow_berg(:,:,k))*frac_in_cloud
                   if (frac_in_cloud_snow_homogeneous .ne. frac_in_cloud) then
                      where (T(:,:,k).lt.T_homogeneous)
                         scav_factor_s(:,:) = f_snow_berg(:,:,k)*frac_in_cloud_snow +&
                              (1.-f_snow_berg(:,:,k))*frac_in_cloud_snow_homogeneous
                      end where
                   end if
                endif
             end if
             !        where (precip3d(:,:,k) > 0.0)
             where (precip3d(:,:,k) > 0. .and. xliq(:,:,k) > 0.)
                w_h2o(:,:) = precip3d(:,:,k) * (AVOGNO/mw_h2o) / zdel(:,:,k) * cm3_2_m3 ! molec/cm3/s
                beta(:,:) = w_h2o(:,:) * mw_h2o  / (n_air(:,:,k) * xliq(:,:,k))
                where (precip3ds(:,:,k) > 0.0 .and. precip3dr(:,:,k) > 0.0)
                   where (scav_factor(:,:) /= scav_factor_s(:,:))
                      scav_factor(:,:) = ( scav_factor(:,:)*precip3dr(:,:,k) +  &
                           scav_factor_s(:,:) * precip3ds(:,:,k)) / precip3d(:,:,k)
                   end where
                elsewhere
                   where (precip3ds(:,:,k) > 0.0)
                      scav_factor(:,:) = scav_factor_s(:,:)
                   endwhere
                endwhere
                in_temp(:,:,k) = beta(:,:) * scav_factor(:,:) ! 1/s
             endwhere
          enddo
          !-----------------------------------------------------------------
          ! Below-cloud wet scavenging
          !-----------------------------------------------------------------
          if( lowercase(scheme)=='henry_below' .or. lowercase(scheme)=='henry_below_noice') then
             k_g = d_g/rain_diam * &
                  ( 2. + 0.6 * sqrt( rain_diam*rain_vterm/vk_air ) * (vk_air/d_g)**(1./3.) )
             conc(:,:,:) = tracer(:,:,:) * n_air(:,:,:) / cm3_2_m3 ! Convert from VMR to molec/m3
             do kk = 1,kd
                do j = 1,jd
                   do i = 1,id
                      stay = 1.
                      if( precip3d(i,j,kk) > 0. ) then
                         conc_rain_total = 0.
                         stay = zfull(i,j,kk) / (rain_vterm * dt)
                         stay = min( stay, 1. )
                         do k = kk,kd
                            f_a0 = Htemp(i,j,k) * pfull(i,j,k) * xliq(i,j,kk) * n_air(i,j,kk)/n_air(i,j,k)
                            scav_factor0 = f_a0 / ( 1.+f_a0 )
                            conc_sat = conc(i,j,k) * scav_factor0 ! molec/m3 <== (xeqca1)
                            sa_drop0 = geo_fac / rain_diam * xliq(i,j,kk) * n_air(i,j,kk) / &
                                 ( DENS_H2O * AVOGNO * cm3_2_m3 ) ! (m2 H2O) / (m3 air)
                            fgas0 = conc(i,j,k) * k_g ! molec/m2/s
                            fall_time = zdel(i,j,k) / rain_vterm ! sec
                            conc_rain = fgas0 * sa_drop0 * fall_time ! molec/m3 <== (xca1)
                            conc_rain_total = conc_rain_total + conc_rain ! molec/m3 <== (all1)
                            if ( conc_rain_total < conc_sat ) then
                               conc(i,j,k) = max( conc(i,j,k)-conc_rain, 0. )
                            end if
                         end do
                         conc(i,j,kk) = conc(i,j,kk) / n_air(i,j,kk) * cm3_2_m3 ! Convert to VMR
                         conc(i,j,kk) = tracer(i,j,kk) - conc(i,j,kk)
                         if ( conc(i,j,kk) /= 0. .and. tracer(i,j,kk) /= 0. ) then
                            fall_time = zdel(i,j,kk)/rain_vterm
                            bc_temp(i,j,kk) = bc_temp(i,j,kk) + &
                                 conc(i,j,kk) / (tracer(i,j,kk) * fall_time) * stay ! 1/s
                         end if
                      end if
                   end do
                end do
             end do

          else if ( lowercase(scheme) == 'aerosol_below' .or. lowercase(scheme) == 'aerosol_below_noice') then

             do k=1,kd
                fluxs = (snow3d(:,:,k+1)+snow3d(:,:,k))*0.5
                fluxr = (rain3d(:,:,k+1)+rain3d(:,:,k))*0.5
                bc_temp(:,:,k) = 0.75 * &
                     (fluxr(:,:)*alpha_r/R_r/DENS_H2O + &
                     fluxs(:,:)*alpha_s/R_s/DENS_SNOW)
             end do

          end if


          do k = 1,kd
             wdep_in(:,:) = wdep_in(:,:) - &
                  in_temp(:,:,k)*tracer(:,:,k)*pwt(:,:,k)*  &
                  cloud_frac(:,:,k)
             wdep_bc(:,:) = wdep_bc(:,:) - &
                  bc_temp(:,:,k)*tracer(:,:,k)*pwt(:,:,k)
          enddo
          dt_temp(:,:,:) = 1. - exp( -bc_temp(:,:,:)*dt ) & ! fractional loss/timestep
               * ( cloud_frac(:,:,:)*exp( -in_temp(:,:,:)*dt ) + (1-cloud_frac(:,:,:)) )
          tracer_dt(:,:,:) = dt_temp(:,:,:) / dt !+ve loss frequency (1/sec)
       endif

    else if(lowercase(scheme)=='fraction') then
       tracer_dt = 0.0
       !-----------------------------------------------------------------------
       !
       !     Compute areal fractions experiencing wet deposition:
       !
       !     Set minimum precipitation rate below which no wet removal
       !     occurs to 0.01 cm/day ie 1.16e-6 mm/sec (kg/m2/s)
       premin=1.16e-6
       !
       !     Large scale cloud liquid water content (kg/m3)
       !     and below cloud washout efficiency (cm-1):
       !flaglw =parse(control,'lslwc',clwc)
       !if (flaglw == 0 ) clwc=0.5e-3

       !JianHe: temporary fix, not used
       if(trim(cloud_param) .eq. 'lscale') clwc=0.5e-3
       wash=1.0
       !
       !     When convective adjustment occurs, use convective cloud liquid water content:
       !
       if(trim(cloud_param) .eq. 'convect') then
          !flaglw = parse(control,'convlwc',clwc)
          !if (flaglw == 0) clwc=2.0e-3
          clwc=2.0e-3
          wash=0.3
       end if
       !
       do j=1,size(rain,2)
          do i=1,size(rain,1)
             tracer_dt(i,j,:)=0.0
             washout(i,j)=0.0
             prenow = rain(i,j) + snow(i,j)
             if(prenow .gt. premin) then
                !
                ! Assume that the top of the cloud is where the highest model level
                ! specific humidity is reduced. And the the bottom of the cloud is the
                ! lowest model level where specific humidity is reduced.
                !
                ktopcd(i,j) = 0
                do k = kd,1,-1
                   if (qdt(i,j,k) < 0.0 ) ktopcd(i,j) = k
                enddo
                kendcd(i,j) = 0
                do k = 1,kd
                   if (qdt(i,j,k) < 0.0 ) kendcd(i,j) = k
                enddo
                !
                !     Thickness of precipitating cloud deck:
                !
                if(ktopcd(i,j).gt.1) then
                   hwtop = 0.0
                   do k=ktopcd(i,j),kendcd(i,j)
                      hwtop=hwtop+(phalf(i,j,k+1)-phalf(i,j,k))*rdgas*T(i,j,k)/grav/pfull(i,j,k)
                   enddo
                   do k=ktopcd(i,j),kendcd(i,j)
                      !     Areal fraction affected by precip clouds (max = 0.5):
                      tracer_dt(i,j,k)=prenow/(clwc*hwtop)
                   end do
                endif

                washout(i,j)=prenow*wash
             endif
          end do
       end do
    endif

    ! Now multiply by the tracer mixing ratio to get the actual tendency.
    tracer_dt(:,:,:) = MIN( MAX(tracer_dt(:,:,:), 0.0E+00), 0.5/dt)
    where (tracer > 0.)
       tracer_dt = tracer_dt*tracer
    else where
       tracer_dt = 0.
    end where

    !++lwh
    !
    ! Re-evaporation
    !
    do k = 1,kd
       where (reevap_fraction(:,:,k) > 0.)
          reevap_diag(:,:,k) = reevap_fraction(:,:,k) * tracer_flux(:,:)
          ! tracer reevaporation fraction is reduced from precip reevaporation,
          ! except when complete reevaporation occurs
          where( reevap_fraction(:,:,k) < 1. )
             reevap_diag(:,:,k) = reevap_diag(:,:,k) * frac_int
          end where
          tracer_dt(:,:,k) = tracer_dt(:,:,k) - reevap_diag(:,:,k) / pwt(:,:,k)
       end where
       tracer_flux(:,:) = tracer_flux(:,:) + tracer_dt(:,:,k)*pwt(:,:,k)
    end do
    !--lwh
    !


    if ( present(so2_so4_out) )then
       so2_so4_out = 0.
       if ( wetdep(n)%is_so2 ) then
          if ( wetdep(n)%so2_so4_evap ) then
             so2_so4_out = reevap_diag / pwt
          end if
       end if
    end if


 endif ! End branching pag/lwh
 !
 ! Output diagnostics in kg/m2/s (if MMR) or mole/m2/s (if VMR)
 if(trim(units) .eq. 'mmr') then
    diag_scale = 1.
 elseif(trim(units) .eq. 'vmr') then
    diag_scale = mw_air ! kg/mole
 else
    diag_scale = 1.
    !if (me == master) then
    !write(*,*) ' Tracer number =',n,' tracer_name=',tracer_name
    !write(*,*) ' scheme=',text_in_scheme
    !write(*,*) ' control=',control
    !write(*,*) ' scheme=',scheme
    !write(*,*) 'Please check field table'
    !write(*,*) 'tracers units =',trim(units),'it should be either  mmr or vmr!'
    !end if
    !  <ERROR MSG="Unsupported tracer units" STATUS="FATAL">
    !     Tracer units must be either VMR or MMR
    !  </ERROR>
    !call errmsg('wet_deposition', 'Unsupported tracer units.', .true. )
 endif

!if (me == master) then
!  write(*,*) ' Tracer number =',n,' tracer_name=',tracer_name
!  write(*,*) 'tracers units =',trim(units),'it should be either  mmr or vmr!'
!endif  

 ! Column integral of wet deposition
 sum_wdep = 0.
 do k=1,kd
    sum_wdep = sum_wdep + tracer_dt(:,:,k)*pwt(:,:,k)/diag_scale
 end do

 if (present (sum_wdep_out))  sum_wdep_out = -sum_wdep

end subroutine gfdl_wet_deposition
!</SUBROUTINE>
!
subroutine wet_deposition_0D( Henry_constant, Henry_variable, &
                              frac_in_cloud, alpha_r, alpha_s, &
                              T, p0, p1, rho_air, &
                              cloud, rain, snow, &
                              tracer, Lgas, Laerosol, Lice, &
                              delta_tracer )
implicit none
!      
!<OVERVIEW>
! Routine to calculate the fraction of tracer removed by wet deposition
!</OVERVIEW>
!
!<IN NAME="T" TYPE="real">
!   Temperature (K)
!</IN>
!<IN NAME="p0" TYPE="real">
!   Pressure (Pa) at layer closer to surface
!</IN>
!<IN NAME="p1" TYPE="real">
!   Pressure (Pa) at layer farther from surface
!</IN>
!<IN NAME="rho_air" TYPE="real">
!   Air density (kg/m3)
!</IN>
!<IN NAME="cloud" TYPE="real">
!   Cloud amount (liquid+ice) (kg/kg)
!</IN>
!<IN NAME="rain" TYPE="real">
!   Precipitation increment (rain) (kg/m3)
!</IN>
!<IN NAME="snow" TYPE="real">
!   Precipitation increment (snow) (kg/m3)
!</IN>
!<IN NAME="tracer" TYPE="real">
!   The tracer field (tracer units)
!</IN>
!<IN NAME="Lgas" TYPE="logical">
!   Is tracer a gas?
!</IN>
!<IN NAME="Laerosol" TYPE="logical">
!   Is tracer an aerosol?
!</IN>
!<IN NAME="Lice" TYPE="logical">
!   Is tracer removed by snow (or only by rain)?
!</IN>
!<OUT NAME="delta_tracer" TYPE="real">
!   The change (increment) of the tracer field due to wet deposition (tracer
!   units)
!/OUT>
!<DESCRIPTION>
! Schemes allowed here are:
!
! 1) Removal according to Henry's Law. This law states that the ratio of the
! concentation in 
!    cloud water and the partial pressure in the interstitial air is a constant.
!    In this 
!    instance, the units for Henry's constant are kg/L/Pa (normally it is
!    M/L/Pa)
!    Parameters for a large number of species can be found at
!    http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!
! 2) Aerosol removal, using specified in-cloud tracer fraction

! To utilize this section of code add one of the following lines as 
! a method for the tracer of interest in the field table.
!<PRE>
! "wet_deposition","henry","henry=XXX, dependence=YYY"
! "wet_deposition","henry_below","henry=XXX, dependence=YYY"
!     where XXX is the Henry's constant for the tracer in question
!       and YYY is the temperature dependence of the Henry's Law constant.
!
! "wet_deposition","aerosol","frac_incloud=XXX"
! "wet_deposition","aerosol_below","frac_incloud=XXX"
!     where XXX is the in-cloud fraction of the aerosol tracer
!</PRE>

!</DESCRIPTION>

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
real,             intent(in)                     :: Henry_constant, Henry_variable, &
                                                    frac_in_cloud, alpha_r, alpha_s
real,             intent(in)                     :: T, p0, p1, rho_air
real,             intent(in)                     :: cloud, rain, snow
real,             intent(in)                     :: tracer
logical,          intent(in)                     :: Lgas, Laerosol, Lice
real,             intent(out)                    :: delta_tracer

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
real :: &
      Htemp, xliq, n_air, pwt, pmid, precip
real :: &
      temp_factor, scav_factor, &
      w_h2o, beta, f_a, in_temp
real, parameter :: &
      GRAV = 9.80,              &  ! acceleration due to gravity [m/s2]
!      RDGAS = 287.04,           &  ! gas constant for dry air [J/kg/deg]
      AVOGNO = 6.023000E+23,    &  ! Avogadro's number
      inv298p15 = 1./298.15,    &  ! 1/K
      cm3_2_m3 = 1.e-6             ! m3/cm3
real, parameter :: mw_air = 28.96440E-03 ! molar mass of air (kg/mole)
real, parameter :: mw_h2o = 18.0E-03  ! molar mass of H2O (kg/mole)

!-----------------------------------------------------------------------

delta_tracer = 0.

pmid = 0.5 * (p0+p1)         ! Pa
pwt     = ( p0 - p1 )/GRAV   ! kg/m2

if( Lgas .or. Laerosol ) then
!++lwh
! units = VMR
!
! Henry_constant (mole/L/Pa) = [X](aq) / Px(g) 
! where [X](aq) is the concentration of tracer X in precipitation (mole/L)
!       Px(g) is the partial pressure of the tracer in the air (Pa)
!
! VMR (total) = VMR (gas) + VMR (aq)
!             = VMR (gas) + [X] * L
!
! where L = cloud liquid amount (kg H2O/mole air)
!
! Using Henry's Law, [X] = H * Px = H * VMR(gas) * Pfull
!
! So, VMR (total) =  VMR(gas) * [ 1 + H * Pfull * L ]
! 
! VMR(gas) = VMR(total) / [1 + H * Pfull * L]
!
! [X] = H * Pfull * VMR(total) / [ 1 + H * Pfull * L]
!
! Following Giorgi and Chameides, JGR, 90(D5), 1985, the first-order loss
! rate constant (s^-1) of X due to wet deposition equals:
!
! k = W_X / n_X
!
! where W_x = the loss rate (molec/cm3/s), and n_X = the number density
! (molec/cm3)
! 
! W_X = [X] * W_H2O / (55 mole/L)
! n_x = VMR(total) * n_air (molec/cm3) = VMR(total) * P/(kT) * 1E-6 m3/cm3
! 
! where P = atmospheric pressure (Pa)
!       k = Boltzmann's constant = 1.38E-23 J/K
!       T = temperature (K)
!       W_H2O = removal increment of water (molec/cm3)
! 
!             [X] * W_H2O / 55         
! So, k = ------------------------------
!         VMR(total) * P/(kT) * 1E-6
! 
!         W_H2O    H * VMR(total) * P / [ 1 + H * P *L ]
!       = ----- * ---------------------------------------
!          55          VMR(total) * P/(kT) * 1E-6
! 
!         W_H2O     H * kT * 1E6
!       = ----- *  -------------    
!          55      1 + H * P * L 
!
!         W_H2O     1     1     H * P * L
!       = ----- * ----- * - * -------------
!          55     n_air   L   1 + H * P * L
!
! where W_H2O = precip (kg/m3) * (AVOGNO/mw_h2o) (molec/kg) * 1E-6 m3/cm3
!
   if( (Lgas .and. Henry_constant > 0.) .or. Laerosol ) then
      if (Lice) then
         precip = rain+snow
      else
         precip = rain
      end if
      in_temp = 0.

      scav_factor = 0.0
      xliq = MAX( cloud * mw_air, 0. ) ! (kg H2O)/(mole air)
      n_air = rho_air * (AVOGNO/mw_air) * cm3_2_m3 ! molec/cm3
      if (Lgas) then
! Calculate the temperature dependent Henry's Law constant
         temp_factor = 1./T-inv298p15
         Htemp = Henry_constant * exp( Henry_variable*temp_factor )
         f_a = Htemp * pmid * xliq
         scav_factor = f_a / ( 1.+f_a )
      else if (Laerosol) then
         scav_factor = frac_in_cloud
      end if
      if (precip > 0. .and. xliq > 0.) then
         w_h2o = precip * (AVOGNO/mw_h2o) * cm3_2_m3 ! molec/cm3
         beta = w_h2o * mw_h2o  / (n_air * xliq)   ! fraction of condensed water removed
         beta = MAX(MIN(beta,1.),0.)
         in_temp = beta * scav_factor              ! fraction of tracer removed
      end if

!     wdep_in = - in_temp*tracer*pwt
!     dt_temp = 1. - exp( -in_temp*dt ) ! fractional loss/timestep
!     tracer_dt = dt_temp / dt !+ve loss frequency (1/sec)
      delta_tracer = in_temp ! fraction of tracer removed
!--lwh
   endif

end if

! Now multiply by the tracer mixing ratio to get the actual tendency.
! tracer_dt = MIN( MAX(tracer_dt, 0.0E+00), 0.5/dt)
if (tracer > 0.) then
   delta_tracer = delta_tracer*tracer
else
   delta_tracer = 0.
end if

! Output diagnostics in kg/m2/s (if MMR) or mole/m2/s (if VMR)
! if(trim(units) .eq. 'mmr') then
!    diag_scale = 1.
! else if(trim(units) .eq. 'vmr') then
!    diag_scale = mw_air ! kg/mole
! else
!    write(*,*) ' Tracer number =',n,' tracer_name=',tracer_name
!    write(*,*) ' scheme=',text_in_scheme
!    write(*,*) ' control=',control
!    write(*,*) ' scheme=',scheme
!    write(*,*) 'Please check field table'
!    write(*,*) 'tracers units =',trim(units),'it should be either  mmr or vmr!'
!  <ERROR MSG="Unsupported tracer units" STATUS="FATAL">
!     Tracer units must be either VMR or MMR
!  </ERROR>
!    call error_mesg('wet_deposition', 'Unsupported tracer units.', FATAL )
! end if

! if(trim(cloud_param) == 'donner') then
!    if (id_tracer_wdep_donin(n) > 0 ) then
!        used = send_data ( id_tracer_wdep_donin(n), wdep_in/diag_scale, Time,
!        is_in=is, js_in=js)
!    endif
!    if(id_tracer_wdep_donin_dt(n) > 0) then
!       used = send_data ( id_tracer_wdep_donin_dt(n), in_temp, Time, is_in=is,
!       js_in=js, ks_in=1)
!    endif
!    if(id_tracer_wdep_don_dt(n) > 0) then
!       used = send_data ( id_tracer_wdep_don_dt(n), dt_temp/dt, Time, is_in=is,
!       js_in=js, ks_in=1)
!    endif
! endif

end subroutine wet_deposition_0D
!</SUBROUTINE>
!
!#######################################################################
!
!
!<SUBROUTINE NAME="tracer_utilities_end">
!<OVERVIEW>
!  The destructor routine for the tracer utilities module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits.
! </DESCRIPTION>

subroutine gfdl_wetdep_end

 module_is_initialized = .FALSE.

end subroutine gfdl_wetdep_end
!</SUBROUTINE>

end module gfdl_wetdep_mod

