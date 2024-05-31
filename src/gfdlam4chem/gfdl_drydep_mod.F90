module gfdl_drydep_mod

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

  !JianHe: read prescribed depvel from external file (drydep_data)
  integer, parameter :: ngas_drydep = 22
   character(len=32), dimension(ngas_drydep), save :: depvel_name = &      !
      (/ "co", "ch2o", "o3", "no", "no2", "hno3", "hno4", "n2o5", "ch4", "ch3ooh", &
         "h2o2", "pan", "pmn", "ch3coch3", "glyc", "hac", "rip", "so2", "nh3", "hobr", &
         "hbr", "brno3" /)

  !----- interfaces -------

  public    &
       gfdl_dry_deposition,    &
       gfdl_drydep_init,gfdl_drydep_end

  !---- version number -----
  character(len=128) :: version = '$Id$'
  character(len=128) :: tagname = '$Name$'

  logical :: module_is_initialized = .FALSE.

  character(len=7), parameter :: mod_name = 'tracers'
  integer, parameter :: max_tracers = 300

  real, parameter :: T_homogeneous = 233.15
  !-----------------------------------------------------------------------
  !--- identification numbers for  diagnostic fields and axes ----
  integer :: id_tracer_ddep(max_tracers), id_tracer_dvel(max_tracers)

  !cmip6 (f1p)
  !


  integer :: id_tracer_ddep_cmip(max_tracers)
  integer :: id_w10m, id_delm
  integer :: id_u_star, id_b_star, id_rough_mom, id_z_pbl,  &
       id_mo_length_inv, id_vds
  !character(len=32),  dimension(max_tracers) :: tracer_names     = ' '
  character(len=32),  dimension(max_tracers) :: tracer_units     = ' '
  character(len=128), dimension(max_tracers) :: tracer_longnames = ' '
  character(len=32),  dimension(max_tracers) :: tracer_ddep_names     = ' '
  character(len=32),  dimension(max_tracers) :: tracer_dvel_names     = ' '
  character(len=32),  dimension(max_tracers) :: tracer_ddep_units     = ' '
  character(len=32),  dimension(max_tracers) :: tracer_dvel_units     = ' '
  character(len=128), dimension(max_tracers) :: tracer_ddep_longnames = ' '
  character(len=128), dimension(max_tracers) :: tracer_dvel_longnames = ' '
  !----------------parameter values for the diagnostic units--------------
  real, parameter :: mw_air = WTMAIR/1000.  ! Convert from [g/mole] to [kg/mole]
  real, parameter :: mw_h2o = WTMH2O/1000.  ! Convert from [g/mole] to [kg/mole]
  real, parameter :: mw_so4 = 96./1000.     ! Convert from [g/mole] to [kg/mole]
  real, parameter :: twopi = 2*PI

  type drydep_type
     character (len=500) :: scheme, name, control
     real  :: land_dry_dep_vel
     real  :: sea_dry_dep_vel
     real  :: ice_dry_dep_vel
     real  :: snow_dry_dep_vel
     real  :: vegn_dry_dep_vel
     logical :: land_does_drydep ! if true, then land model handles dry deposition ,
     ! over land surfaces and therefore this module should scale down the deposition
     ! it calculates by 1 - fraction of land
     logical :: Ldrydep
     real    :: surfr,snowr,landr,sear
  end type drydep_type

  type(drydep_type), dimension(:), allocatable :: Drydep


  ! --->h1g, add a scale factor for aerosol wet deposition, 2014-04-10
  real ::                scale_aerosol_wetdep =1.0
  real ::                scale_aerosol_wetdep_snow =1.0
  character(len=64)  :: file_dry = 'depvel.nc'  ! NetCDF file for dry deposition velocities
  logical :: drydep_exp = .true.
  real :: T_snow_dep = 263.15
  real :: kbs_val   = 0.5 ! surface conductance of rough sea (m/s)
  ! <---h1g,
contains

  !
  ! ######################################################################
  !
  !<SUBROUTINE NAME="gfdl_drydep_init">
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
  ! call gfdl_drydep_init(lonb,latb, mass_axes, Time)
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

  subroutine gfdl_drydep_init(me,master,tracer_names)

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

!---------------------------------------------------------------------
!  make sure that astronomy_mod has been initialized (if radiation
!  not being called in this run, it will not have previously been
!  initialized).
!---------------------------------------------------------------------
    call astronomy_init

!    do n = 1, max_tracers
!       write ( tracer_names(n),     100 ) n
!       write ( tracer_longnames(n), 102 ) n
!       tracer_units(n) = 'none'
!    enddo
!100 format ('tr',i3.3)
!102 format ('tracer ',i3.3)

    !call get_number_tracers(MODEL_ATMOS, num_tracers= ntrace)
    ntrace = size(tracer_names)

    if (ntrace > 0) then
       allocate (Drydep(ntrace))
    endif

    if (me == master) then
      print *, "Total number of tracers passed to drydep: ", ntrace
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
!          tracer_ddep_names(n) = trim(tracer_names(n)) //'_ddep'
!          tracer_dvel_names(n) = trim(tracer_names(n)) //'_dvel'
!       endif
!       write (name,102) n
!       if (trim(tracer_longnames(n)) /= name) then
!          tracer_ddep_longnames(n) = &
!               trim(tracer_longnames(n)) // ' dry deposition for tracers'
!          tracer_dvel_longnames(n) = &
!               trim(tracer_longnames(n)) // ' dry deposition velocity for tracers'
!       endif

       select case (trim(tracer_units(n)))
       case ('mmr')
          units = 'kg/m2/s'
       case ('kg/kg')
          units = 'kg/m2/s'
       case ('vmr')
          units = 'mole/m2/s'
       case ('mol/mol')
          units = 'mole/m2/s'
       case ('mole/mole')
          units = 'mole/m2/s'
       case default
          units = trim(tracer_units(n))//' kg/(m2 s)'
          call errmsg('gfdl_drydep_init',&
               ' Dry dep units set to '//trim(units)//' in gfdl_drydep for '//trim(tracer_names(n)),&
               .true.)
       end select

       !JianHe: hardcoded here, but would be flexible to have in field table or
       !external file in the future.
       Drydep(n)%Ldrydep = .False.
       Drydep(n)%name = tracer_names(n)
       Drydep(n)%land_does_drydep = .False.
       Drydep(n)%scheme = "None"
       Drydep(n)%land_dry_dep_vel = 0.
       Drydep(n)%sea_dry_dep_vel = 0.

       Drydep(n)%surfr = 500.
       Drydep(n)%landr = 500.
       Drydep(n)%snowr = 500.
       Drydep(n)%sear = 500.

       if ((name=='soa') .or. (name=='so4') .or. (name=='sulf') .or. (name=='dust1') &
          .or. (name=='seas1') .or. (name=='bc1') .or. (name=='bc2') &
          .or. (name=='oc1') .or. (name=='oc2') .or. (name=='nh4')) then
           Drydep(n)%Ldrydep = .True.
           !Drydep(n)%control
           Drydep(n)%scheme = "williams_wind_driven"
           Drydep(n)%surfr = 200.
           Drydep(n)%landr = 200.
           Drydep(n)%snowr = 1200.
           Drydep(n)%sear = 200.
       end if
    
       if ((name=='dust2') .or. (name=='dust3') .or. (name=='dust4') &
          .or. (name=='dust5') .or. (name=='seas2') .or. (name=='nh4no3')) then
           Drydep(n)%Ldrydep = .True.
           !Drydep(n)%control
           Drydep(n)%scheme = "williams_wind_driven"
           Drydep(n)%surfr = 100.
           Drydep(n)%landr = 100.
           Drydep(n)%snowr = 1200.
           Drydep(n)%sear = 100.
       end if

       if ((name=='seas3') .or. (name=='seas4') .or. (name=='seas5')) then
           Drydep(n)%Ldrydep = .True.
           !Drydep(n)%control
           Drydep(n)%scheme = "williams_wind_driven"
           Drydep(n)%surfr = 50.
           Drydep(n)%landr = 50.
           Drydep(n)%snowr = 1200.
           Drydep(n)%sear = 50.
       end if

       
       if ((name=='co') .or. (name=='ch2o') .or. (name=='o3') &
          .or. (name=='no') .or. (name=='no2') .or. (name=='hno3') &
          .or. (name=='ho2no2') .or. (name=='n2o5') & ! .or. (name=='ch4') & ! exclude ch4 for now
          .or. (name=='ch3ooh') .or. (name=='h2o2') .or. (name=='pan') &
          .or. (name=='mpan') .or. (name=='ch3coch3') .or. (name=='ch3oh') &
          .or. (name=='c2h5oh') .or. (name=='glyald') .or. (name=='hyac') &
          .or. (name=='isopooh') .or. (name=='h2') .or. (name=='so2') &
          .or. (name=='nh3') .or. (name=='hcl') .or. (name=='hobr') &
          .or. (name=='hbr') .or. (name=='bronor') .or. (name=='atooh') &
          .or. (name=='glyx') .or. (name=='iepox') .or. (name=='isopnb') &
          .or. (name=='macrn') .or. (name=='mgly') .or. (name=='mvkn') &
          .or. (name=='r4n1') .or. (name=='r4n2') .or. (name=='o3s') &
          .or. (name=='o3s_e90')) then
           Drydep(n)%Ldrydep = .True.
           !Drydep(n)%control
           Drydep(n)%scheme = "file"
        end if

      if (name=='dms') then
          Drydep(n)%Ldrydep = .True.
          Drydep(n)%scheme = "fixed"
          Drydep(n)%land_dry_dep_vel=0.11e-2
          Drydep(n)%sea_dry_dep_vel = 0.09e-2
      end if

    ! check that the corresonding land tracer is present in the land if the
    ! dry deposition is done on the land side
    if (Drydep(n)%land_does_drydep) then
          call errmsg('gfdl_drydep_init', &
               'land_does_drydep is not working!', .true.)
    endif

    if (me.eq.master) then
       write(*,*) 'name: ',trim(tracer_names(n))
       write(*,*) 'scheme: ',trim(Drydep(n)%scheme)
       write(*,*) 'land does drydep: ',Drydep(n)%land_does_drydep
       write(*,*) 'land dvel: ',Drydep(n)%land_dry_dep_vel
       write(*,*) 'sea dvel: ',Drydep(n)%sea_dry_dep_vel
       write(*,*) 'surfr,landr,snowr,sear: ',Drydep(n)%surfr,Drydep(n)%landr,Drydep(n)%snowr,Drydep(n)%sear
    end if

 enddo

! call write_version_number (version, tagname)

! if ( me == master ) then
!    logunit=101
!    call write_namelist_values (logunit,ntrace)
! endif

 module_is_initialized = .TRUE.

end subroutine gfdl_drydep_init


!####################################################################


!</SUBROUTINE>
!
!#######################################################################
!
subroutine write_namelist_values (unit, ntrace)
 integer, intent(in) :: unit, ntrace
 integer :: n

 write (unit,10)
 do n = 1, ntrace
    write (unit,11) trim(tracer_ddep_names(n)),     &
         trim(tracer_ddep_longnames(n)), &
         trim(tracer_ddep_units(n))
    write (unit,11) trim(tracer_dvel_names(n)),     &
         trim(tracer_dvel_longnames(n)), &
         'm/s'
 enddo

10 format (' &TRACER_DIAGNOSTICS_NML', &
      /,'    TRACER:  names  longnames  (units)')
11 format (a16,2x,a,2x,'(',a,')')

end subroutine write_namelist_values

!
!#######################################################################
!
!<SUBROUTINE NAME = "dry_deposition">
subroutine gfdl_dry_deposition( me, master, tracer_names, &
    n, is, js, u, v, T, pwt, pfull, dz, &
    u_star, landfrac, frac_open_sea,dsinku, dt, tracer,  &
    lat, lon, Time, drydep_data, drydep_vel, con_atm)
  ! When formulation of dry deposition is resolved perhaps use the following?
  !                           landfr, seaice_cn, snow_area, &
  !                           vegn_cover, vegn_lai, &
  !                           b_star, z_pbl, rough_mom )
  !
  !<OVERVIEW>
  ! Routine to calculate the fraction of tracer to be removed by dry
  ! deposition.
  !</OVERVIEW>
  !<DESCRIPTION>
  ! There are three types of dry deposition coded.
  !
  ! 1) Wind driven derived dry deposition velocity.
  !
  ! 2) Fixed dry deposition velocity.
  !
  ! 3) Dry deposition velocities read in from input file

  ! There are an addition three types of dry deposition coded
  ! but presently commented out.
  !
  ! 4) Wind driven derived dry deposition velocity, surface dependent.
  !
  ! 5) Wind driven derived dry deposition velocity, surface and boundary
  !    layer stability dependent.
  !
  ! 6) Fixed dry deposition velocity, difference between land, snow-covered
  !    land, sea and ice-covered sea.
  !
  ! The theory behind the wind driven dry deposition velocity calculation
  ! assumes that the deposition can be modeled as a parallel resistance type
  ! problem.
  !
  !  Total resistance to HNO3-type dry deposition,
  !<PRE>       R = Ra + Rb
  !  resisa = aerodynamic resistance
  !  resisb = surface resistance (laminar layer + uptake)
  !         = 5/u*  [s/cm]        for neutral stability
  !      Vd = 1/R
  !</PRE>
  ! For the fixed dry deposition velocity, there is no change in the
  ! deposition velocity but the variation of the depth of the surface
  ! layer implies that there is variation in the amount deposited.
  !
  ! To utilize this section of code add one of the following lines as
  ! a method for the tracer of interest in the field table.
  !<PRE>
  ! "dry_deposition","wind_driven","surfr=XXX"
  !     where XXX is the total resistance defined above.
  !
  ! "dry_deposition","fixed","land=XXX, sea=YYY"
  !     where XXX is the dry deposition velocity (m/s) over land
  !       and YYY is the dry deposition velocity (m/s) over sea.
  !
  ! "dry_deposition","file","FILENAME.NC"
  !     where FILENAME.NC is the NetCDF file name.
  !</PRE>
  !</DESCRIPTION>
  !<TEMPLATE>
  ! call dry_deposition( n, is, js, u, v, T, pwt, pfull, dz,
  !                      u_star, landfrac, dsinku, tracer, Time, drydep_data)
  !</TEMPLATE>
  !
  !  <IN NAME="n" TYPE="integer">
  !    The tracer number.
  !  </IN>
  !  <IN NAME="is, js" TYPE="integer">
  !    Start indices for array (computational indices).
  !  </IN>
  !  <IN NAME="u" TYPE="real" DIM="(:,:)">
  !    U wind field.
  !  </IN>
  !  <IN NAME="v" TYPE="real" DIM="(:,:)">
  !    V wind field.
  !  </IN>
  !  <IN NAME="T" TYPE="real" DIM="(:,:)">
  !    Temperature.
  !  </IN>
  !  <IN NAME="pwt" TYPE="real" DIM="(:,:)">
  !     Pressure differential of half levels.
  !  </IN>
  !  <IN NAME="pfull" TYPE="real" DIM="(:,:)">
  !     Full pressure levels.
  !  </IN>
  !  <IN NAME="u_star" TYPE="real" DIM="(:,:)">
  !     Friction velocity.
  !  </IN>
  !  <IN NAME="lon" TYPE="real" DIM="(:,:)">
  !     Longitude.
  !  </IN>
  !  <IN NAME="landfrac" TYPE="logical">
  !     Fraction of land in a grid cell.
  !  </IN>
  !  <INOUT NAME="drydep_data" TYPE="interpolate_type">
  !     Dry deposition data interpolated from input file.
  !  </INOUT>
  !
  !  <OUT NAME="dsinku" TYPE="real" DIM="(:,:)">
  !    The amount of tracer in the surface layer which is dry deposited per second.
  !  </OUT>
  !
 integer, intent (in) :: me
 integer, intent (in) :: master
 character(len=32), intent(in) :: tracer_names(:)

 integer, intent(in)                 :: n, is, js
 real, intent(in), dimension(:,:)    :: u, v, T, pwt, pfull, u_star, tracer, dz
 real, intent(in), dimension(:,:)    :: lon, lat
 real, intent(in), dimension(:,:)    :: landfrac,frac_open_sea
 real, intent(in), dimension(:,:), optional    :: con_atm
 ! When formulation of dry deposition is resolved perhaps use the following?
 !real, intent(in), dimension(:,:)    :: landfr, z_pbl, b_star, rough_mom
 !real, intent(in), dimension(:,:)    :: seaice_cn, snow_area, vegn_cover,  &
 !                                       vegn_lai
 type(time_type), intent(in)         :: Time !, Time_next
 real, intent(in), dimension(:,:,:)  :: drydep_data
 real, intent(in)                   :: dt
 real, intent(out), dimension(:,:)   :: dsinku
 real, intent(out), dimension(:,:)   :: drydep_vel

 real,dimension(size(u,1),size(u,2))   :: hwindv,frictv,resisa,ka,kss,kbs,km,vd_ocean,A,B,alpha,landr2
 !real,dimension(size(u,1),size(u,2))   :: mo_length_inv, vds, rs, k1, k2
 integer :: i,j, flagsr, id, jd
 real    :: land_dry_dep_vel, sea_dry_dep_vel, ice_dry_dep_vel,  &
      snow_dry_dep_vel, vegn_dry_dep_vel,   &
      surfr, sear,  snowr, vegnr, landr
 real    :: diag_scale
 real    :: factor_tmp, gmt, dv_on, dv_off, dayfrac, vd_night, vd_day, loc_angle
 logical :: used, diurnal
 integer :: flag_species, flag_diurnal
 character(len=10) ::units,names
 character(len=500) :: name,control,scheme, speciesname,dummy
 real, dimension(size(tracer,1),size(tracer,2)) :: coszen, fracday, half_day
 real :: rrsun

  call diurnal_solar( lat, lon, Time, coszen, fracday, &
                      rrsun, dt_time=real_to_time_type(dt), &
                      half_day_out=half_day )


 ! Default zero
 dsinku = 0.0
 if (.not. Drydep(n)%Ldrydep) return
 name =Drydep(n)%name
 control = Drydep(n)%control
 scheme = Drydep(n)%scheme
 land_dry_dep_vel = Drydep(n)%land_dry_dep_vel
 sea_dry_dep_vel = Drydep(n)%sea_dry_dep_vel
 !ice_dry_dep_vel = Drydep(n)%ice_dry_dep_vel
 !snow_dry_dep_vel = Drydep(n)%snow_dry_dep_vel
 !vegn_dry_dep_vel = Drydep(n)%vegn_dry_dep_vel

 ! delta z = dp/(rho * grav)
 ! delta z = RT/g*dp/p    pwt = dp/g
 !dz(:,:) = pwt(:,:)*rdgas*T(:,:)/pfull(:,:)
 id=size(pfull,1); jd=size(pfull,2)

 surfr=Drydep(n)%surfr
 snowr=Drydep(n)%snowr
 landr=Drydep(n)%landr
 sear=Drydep(n)%sear
 

!    if (me.eq.master) then
!     if (name == "so4") then
!       write(*,*) 'name: ',trim(tracer_names(n))
!       write(*,*) 'scheme: ',trim(Drydep(n)%scheme)
!       write(*,*) 'land does drydep: ',Drydep(n)%land_does_drydep
!       write(*,*) 'land dvel: ',Drydep(n)%land_dry_dep_vel
!       write(*,*) 'sea dvel: ',Drydep(n)%sea_dry_dep_vel
!       write(*,*) 'surfr,landr,snowr,sear: ',Drydep(n)%surfr,Drydep(n)%landr,Drydep(n)%snowr,Drydep(n)%sear
!     end if
!   end if

 select case(lowercase(scheme))

 
 case ('williams_wind_driven')

    where(T.lt.T_snow_dep)
       landr2=snowr
    elsewhere
       landr2=landr
    endwhere

    frictv=u_star
    where (frictv .lt. 0.1) frictv=0.1

    hwindv=sqrt(u**2+v**2)

    if (present(con_atm)) then
       ka = con_atm
       resisa = 1/max(con_atm,1.e-25)
    else
       resisa=hwindv/(u_star*u_star)
       ka=1./max(resisa,1.e-25)
    end if

    !f1p
    !alpha = max(min(1.7e-6*hwindv**3.75,1.),0.) !WU 1979
    !alpha = max(min(1.e-6*u_star**3,0.) !Wu 1988, variations of whitecap coverage with wind stress and water temperature
!   alpha = max(min(2.81e-5*(hwindv-3.87)**2.76,1.),0.) !Observations of whitecap coverage and the relation to wind stress, wave slope, and turbulent dissipation, Schwendeman and Thomson, 2016, JGR ocean
!Observations of whitecap coverage and the relation to wind stress, wave slope, and turbulent dissipation, Schwendeman and Thomson, 2016, JGR ocean
    where (hwindv .gt. 3.87)
       alpha = max(min(2.81e-5*(hwindv-3.87)**2.76,1.),0.)
    elsewhere
       alpha = 0.
    endwhere

    kss   = frictv/sear
    kbs   = kbs_val !set to very high value
    km    = hwindv !lateral transport

    A = km*ka+(1.-alpha)*ka*alpha*(ka+kbs)
    B = km*((1.-alpha)*(ka+kss)+alpha*(ka+kbs))+(1.-alpha)*(ka+kss)*alpha*(ka+kbs)
    
    vd_ocean = A/B*((1.-alpha)*kss &
         + km * alpha * kbs/(km + alpha*(ka+kbs)) &
         + alpha * kbs * alpha * ka / (km+alpha*(ka+kbs)))

    drydep_vel(:,:) = (1./(landr2/frictv + resisa)) * (1.-frac_open_sea) &
         +     frac_open_sea * vd_ocean

    dsinku = drydep_vel(:,:)/dz(:,:)

 case('wind_driven')
    ! Calculate horizontal wind velocity and aerodynamic resistance:
    !   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2),
    !   and  u*=sqrt(momentum flux)  is friction velocity.
    !
    !****  Compute dry sinks (loss frequency, need modification when
    !****    different vdep values are to be used for species)
    hwindv=sqrt(u**2+v**2)
    frictv=u_star
    resisa=hwindv/(u_star*u_star)
    where (frictv .lt. 0.1) frictv=0.1
    drydep_vel(:,:) = (1./(surfr/frictv + resisa))
    dsinku = drydep_vel(:,:)/dz(:,:)

    !    case('sfc_dependent_wind_driven')
    !! Calculate horizontal wind velocity and aerodynamic resistance:
    !!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2),
    !!   and  u*=sqrt(momentum flux)  is friction velocity.
    !!
    !!****  Compute dry sinks (loss frequency, need modification when
    !!****    different vdep values are to be used for species)
    !        flagsr=parse(control,'surfr',surfr)
    !        if(flagsr == 0) surfr=500.
    !
    !        flagsr=parse(control,'sear',sear)
    !        if(flagsr == 0) sear=surfr
    !
    !        flagsr=parse(control,'icer',icer)
    !        if(flagsr == 0) icer=surfr
    !
    !        flagsr=parse(control,'snowr',snowr)
    !        if(flagsr == 0) snowr=surfr
    !
    !        flagsr=parse(control,'vegnr',vegnr)
    !        if(flagsr == 0) vegnr=surfr
    !
    !        hwindv=sqrt(u**2+v**2)
    !        frictv=u_star
    !        resisa=hwindv/(u_star*u_star)
    !        where (frictv .lt. 0.1) frictv=0.1
    !        drydep_vel(:,:) = (1./(surfr/frictv + resisa))*  &
    !                                         (landfr(:,:) - snow_area(:,:)) + &
    !                          (1./(snowr/frictv + resisa))*snow_area(:,:)  +&
    !                          (1./(sear/frictv + resisa))*   &
    !                                 (1. - landfr(:,:) - seaice_cn(:,:))   + &
    !                          (1./(icer/frictv + resisa))*seaice_cn(:,:)
    !        dsinku(:,:) = drydep_vel(:,:) / dz(:,:)
    !
    !    case('sfc_BL_dependent_wind_driven')
    !! Calculate horizontal wind velocity and aerodynamic resistance:
    !!   where xxfm=(u*/u) is drag coefficient, Ra=u/(u*^2),
    !!   and  u*=sqrt(momentum flux)  is friction velocity.
    !!
    !!****  Compute dry sinks (loss frequency, need modification when
    !!****    different vdep values are to be used for species)
    !        flagsr=parse(control,'surfr',surfr)
    !        if(flagsr == 0) surfr=500.
    !
    !        flagsr=parse(control,'sear',sear)
    !        if(flagsr == 0) sear=surfr
    !
    !        flagsr=parse(control,'icer',icer)
    !        if(flagsr == 0) icer=surfr
    !
    !        flagsr=parse(control,'snowr',snowr)
    !        if(flagsr == 0) snowr=surfr
    !
    !        flagsr=parse(control,'vegnr',vegnr)
    !        if(flagsr == 0) vegnr=surfr
    !
    !        hwindv=sqrt(u**2+v**2)
    !        frictv=u_star
    !        resisa=hwindv/(u_star*u_star)
    !        where (frictv .lt. 0.1) frictv=0.1
    !        mo_length_inv = - vonkarm * b_star/(frictv*frictv)
    !        where (rough_mom > 0.005)
    !           k1 =  0.001222 * log10(rough_mom) + 0.003906
    !        elsewhere
    !           k1 =  0.001222 * log10(0.005) + 0.003906
    !        endwhere
    !        where(mo_length_inv < 0)
    !           k2 = 0.0009 * ( - z_pbl * mo_length_inv)**(2.0/3.0)
    !        elsewhere
    !           k2 = 0.0
    !        endwhere
    !        vds = frictv * (k1 + k2)
    !        rs  = 1.0 / vds
    !        drydep_vel(:,:) = (1./(rs + resisa))*      &
    !                                   (landfr(:,:) - snow_area(:,:))  +  &
    !                          (1./(snowr/frictv + resisa))*snow_area(:,:)  + &
    !                          (1./(sear/frictv + resisa))*   &
    !                                (1. - landfr(:,:) - seaice_cn(:,:))  + &
    !                          (1./(icer/frictv + resisa))*seaice_cn(:,:)
    !        dsinku(:,:) = drydep_vel(:,:) / dz(:,:)

 case('fixed')
    !JianHe: only for DMS
    ! For the moment let's try to calculate the delta-z of the bottom
    ! layer and using a simple dry deposition velocity times the
    ! timestep, idt, calculate the fraction of the lowest layer which
    ! deposits.
    where (landfrac(:,:)> 0.5 )
       ! dry dep value over the land surface
       drydep_vel(:,:) = land_dry_dep_vel
    elsewhere
       ! dry dep value over the sea surface
       drydep_vel(:,:) = sea_dry_dep_vel
    endwhere
    dsinku(:,:) = drydep_vel(:,:) / dz(:,:)

    !    case('sfc_dependent_fixed')
    !      drydep_vel(:,:) = land_dry_dep_vel*    &
    !                                     (landfr(:,:) - snow_area(:,:))   + &
    !                        snow_dry_dep_vel*snow_area(:,:)   +  &
    !                        sea_dry_dep_vel*  &
    !                                  (1. - landfr(:,:) - seaice_cn(:,:))  + &
    !                        ice_dry_dep_vel*seaice_cn(:,:)
    !      dsinku(:,:) = drydep_vel(:,:) / dz(:,:)

 case('file')
    !flag_species = parse(control,'name',speciesname)
    !if(flag_species>0) then
    !   name = trim(speciesname)
    !else
    !   call get_tracer_names(MODEL_ATMOS,n,name)
    !endif
    
    !flag_diurnal = parse(control,'diurnal',dummy)
    !diurnal = (flag_diurnal > 0)

    name = tracer_names(n)
    if ((tracer_names(n)=='o3') .or. (tracer_names(n)=='o3s') &
        .or. (tracer_names(n)=='o3s_e90')) then
      diurnal = .True.
    else
      diurnal = .False.
    endif
 
    !JianHe: hard-coded here
    if ((tracer_names(n)=="co") .or. &
        (tracer_names(n)=="h2")) then
       drydep_vel(:,:) = drydep_data(:,:,1)
    else if (tracer_names(n) == "ch2o") then
       drydep_vel(:,:) = drydep_data(:,:,2)
    else if ((tracer_names(n)=='o3') .or. &
             (tracer_names(n)=='o3s') .or. &
             (tracer_names(n)=='o3s_e90')) then
       drydep_vel(:,:) = drydep_data(:,:,3)
    else if (tracer_names(n)=='no') then
       drydep_vel(:,:) = drydep_data(:,:,4)
    else if (tracer_names(n)=='no2') then
       drydep_vel(:,:) = drydep_data(:,:,5)
    else if ((tracer_names(n)=='hno3') .or. &
             (tracer_names(n)=='hcl') .or. & 
             (tracer_names(n)=='isopnb') .or. &
             (tracer_names(n)=='macrn') .or. &
             (tracer_names(n)=='mvkn') .or. &
             (tracer_names(n)=='r4n1') .or. &
             (tracer_names(n)=='r4n2')) then
       drydep_vel(:,:) = drydep_data(:,:,6)
    else if (tracer_names(n)=='ho2no2') then
       drydep_vel(:,:) = drydep_data(:,:,7)
    else if (tracer_names(n)=='n2o5') then
       drydep_vel(:,:) = drydep_data(:,:,8)
    else if (tracer_names(n)=='ch4') then
       drydep_vel(:,:) = drydep_data(:,:,9)
    else if ((tracer_names(n)=='ch3ooh') .or. &
             (tracer_names(n)=='ch3oh') .or. &
             (tracer_names(n)=='c2h5oh') .or. &
             (tracer_names(n)==' atooh')) then
       drydep_vel(:,:) = drydep_data(:,:,10)
    else if ((tracer_names(n)=='h2o2') .or. &
             (tracer_names(n)=='glyx') .or. &
             (tracer_names(n)=='iepox') .or. &
             (tracer_names(n)=='mgly')) then
       drydep_vel(:,:) = drydep_data(:,:,11)
    else if (tracer_names(n)=='pan') then
       drydep_vel(:,:) = drydep_data(:,:,12)
    else if (tracer_names(n)=='mpan') then
       drydep_vel(:,:) = drydep_data(:,:,13)
    else if (tracer_names(n)=='ch3coch3') then
       drydep_vel(:,:) = drydep_data(:,:,14)
    else if (tracer_names(n)=='glyald') then
       drydep_vel(:,:) = drydep_data(:,:,15)
    else if (tracer_names(n)=='hyac') then
       drydep_vel(:,:) = drydep_data(:,:,16)
    else if (tracer_names(n)=='isopooh') then
       drydep_vel(:,:) = drydep_data(:,:,17)
    else if (tracer_names(n)=='so2') then
       drydep_vel(:,:) = drydep_data(:,:,18)
    else if (tracer_names(n)=='nh3') then
       drydep_vel(:,:) = drydep_data(:,:,19)
    else if (tracer_names(n)=='hobr') then
       drydep_vel(:,:) = drydep_data(:,:,20)
    else if (tracer_names(n)=='hbr') then
       drydep_vel(:,:) = drydep_data(:,:,21)
    else if (tracer_names(n)=='brono2') then
       drydep_vel(:,:) = drydep_data(:,:,22)
    else
       drydep_vel(:,:) = 0.
    end if

    if (diurnal) then
       do j = 1,jd
          do i = 1,id
             ! half_day is between 0 and pi, so dv_off btwn 0 to pi, dv_on btwn -pi and 0
             dv_off = MIN( 1.2*half_day(i,j), PI )
             dv_on = -dv_off
             dayfrac = dv_off/PI
             ! apply the mean dep vel during polar day or polar night (or nearby)
             if (dv_off > 0 .and. dv_off < PI  ) then
                vd_night = MIN(0.001, 0.5*drydep_vel(i,j))
                vd_day = ( drydep_vel(i,j)-vd_night*(1.-dayfrac) ) / dayfrac
                gmt = universal_time(Time)
                loc_angle = gmt + lon(i,j) - PI
                if (loc_angle >= PI) loc_angle = loc_angle - twopi
                if (loc_angle < -PI) loc_angle = loc_angle + twopi
                if( loc_angle >= dv_off .or. loc_angle <= dv_on ) then
                   drydep_vel(i,j) = vd_night
                else
                   factor_tmp = loc_angle - dv_on
                   factor_tmp = factor_tmp / MAX(2*dv_off,1.e-6)
                   drydep_vel(i,j) = 0.5*PI*sin(factor_tmp*PI)*(vd_day-vd_night) + vd_night
                end if
             end if

          !   if (me == master) then
          !     if (drydep_vel(i,j) < 0.) then
          !       print *, "Negative drydep_vel for ", tracer_names(n)
          !     end if
          !   end if
          end do
       end do
    end if !(diurnal)

    dsinku(:,:) = drydep_vel(:,:) / dz(:,:)
 case('default')
    drydep_vel(:,:) = 0.
 end select

 if (Drydep(n)%land_does_drydep) then
    ! land handles dry deposition, so we need to scale the calculated values of
    ! sink by the fraction of the non-land in the grid cell
    dsinku = dsinku*(1-landfrac)
 endif
 dsinku(:,:) = MAX(dsinku(:,:), 0.0E+00)
 if ( drydep_exp ) then
    where(tracer>0)
       dsinku=tracer*(1. - exp(-dsinku*dt))/dt
    elsewhere
       dsinku=0.0
    endwhere
 else
    where(tracer>0)
       dsinku=dsinku*tracer
    elsewhere
       dsinku=0.0
    endwhere
 end if

end subroutine gfdl_dry_deposition
!</SUBROUTINE>
!
!#######################################################################
!
!######################################################################
!<SUBROUTINE NAME="tracer_utilities_end">
!<OVERVIEW>
!  The destructor routine for the tracer utilities module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits.
! </DESCRIPTION>

subroutine gfdl_drydep_end


 module_is_initialized = .FALSE.

end subroutine gfdl_drydep_end
!</SUBROUTINE>


end module gfdl_drydep_mod

