module atmos_regional_tracer_driver_mod
!
! <CONTACT EMAIL="Larry.Horowitz@noaa.gov">
!   Larry W. Horowitz
! </CONTACT>

! <OVERVIEW>
!     This code calculates tendencies for regional chemical tracers
! </OVERVIEW>

! <DESCRIPTION>
!
! This code calculates emissions and chemical production/loss of regional
! chemical tracers tracers.
!
! Add more chemical tracers in extra regions like: central Africa, Russian Asia
! and southeastern Asia ---- by yyf, Jan 27, 2009
!
! </DESCRIPTION>


!-----------------------------------------------------------------------

use mo_chem_utls_mod, only : get_spc_ndx, get_tracer_ndx, &
                             NO_TRACER
use catchem_constants,     only : WTMAIR, AVOGNO,       &
                                  DEG_TO_RAD
use mo_errmsg,             only : errmsg

implicit none

private

!-----------------------------------------------------------------------
!     ... interfaces
!-----------------------------------------------------------------------
public  regional_tracer_driver, regional_tracer_driver_init

!-----------------------------------------------------------------------
!     ... namelist
!-----------------------------------------------------------------------
!  set a maximum number of tracers to assign the array
!  and update the actual number (ntracers) from the namelist
!  ntracers is read in regional_tracer_driver_init and 
!  must be publicly accessiable in the entire module
!-----------------------------------------------------------------------
!JianHe: hardcoded here
integer, parameter :: ntracers = 3
integer, parameter :: max_ntracers = 50
real, dimension(ntracers) :: loss_freq = (/1.286008e-07, 0., 2.314815e-07/)
real, parameter           :: loss_freq_all = 0.
real, parameter           :: co_yield_from_avoc = 0.7, &
                                 co_yield_from_bvoc = 0.4, &
                                 co_yield_from_ch4  = 0.86,&
                                 tch4_vmr = 0.
real, parameter :: sec_per_day    = 86400., &
                   age_relax_time = 0.04166, & ! timescale for relaxation to zero (days)
                   k_relax_aoanh  = 1./(age_relax_time*sec_per_day), & ! (1/sec)
                   k_relax_nh50   = 1./(age_relax_time*sec_per_day), & ! (1/sec)
                   days_per_year  = 365.25, &
                   k_aging        = 1./(days_per_year*sec_per_day), & ! increase age at 1 yr/yr (convert to yr/sec)
                   lat30=30.*DEG_TO_RAD, lat50=50.*DEG_TO_RAD, & ! tagged latitude (deg) for aoanh/nh50 tracers
                   nh50_fixed_val = 100.e-9 ! prescribed VMR of NH50 tracer (in surface layer, 30-50N)
character(len=16), dimension(ntracers) :: regtrc_names = (/"e90", "aoanh", "nh50"/)
 
!namelist /regional_tracer_driver_nml/     &
!                               ntracers,  &
!                               loss_freq, &
!                               loss_freq_all, &
!                               regtrc_names, &
!                               co_yield_from_avoc, &
!                               co_yield_from_bvoc, &
!                               co_yield_from_ch4, &
!                               tch4_vmr

!-----------------------------------------------------------------------
!     ...  declare type that will store the field infomation for the 
!          emission file
!-----------------------------------------------------------------------

type,public :: field_init_type
   character(len=64), pointer :: field_names(:)
end type field_init_type



character(len=7), parameter :: module_name = 'tracers'
real, parameter :: g_to_kg    = 1.e-3,    & !conversion factor (kg/g)
                   m2_to_cm2  = 1.e4        !conversion factor (cm2/m2)
real, parameter :: emis_cons = WTMAIR * g_to_kg * m2_to_cm2 / AVOGNO
logical, dimension(max_ntracers) :: Lemis = .false.
     
integer, dimension(max_ntracers) :: tracer_indices = 0

logical :: module_is_initialized=.false.

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer, dimension(max_ntracers) :: id_prod, id_loss, id_chem_tend, id_emiss

!-----------------------------------------------------------------------
!     ... tracer numbers
!-----------------------------------------------------------------------
integer :: id_tch4=0, id_avoc=0, id_bvoc=0, &
           id_cofromch4=0, id_cofromavoc=0, id_cofrombvoc=0, &
           id_aoanh=0, id_nh50=0

!---- version number ---------------------------------------------------
character(len=128), parameter :: version     = ''
character(len=128), parameter :: tagname     = ''
!-----------------------------------------------------------------------

contains


!#######################################################################

! <SUBROUTINE NAME="regional_tracer_driver">
!   <OVERVIEW>
!     Regional tracer driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine calculates the sources and sinks of regional tracers
!     for "CO" and "SA"
!   </DESCRIPTION>
!   <TEMPLATE>
!     call regional_tracer_driver( lon, lat, pwt, r, chem_dt, &
!                                  Time, phalf, is, js, kbot)
!   </TEMPLATE>
!   <IN NAME="lon" TYPE="real" DIM="(:,:)">
!     The longitudes for the local domain.
!   </IN>
!   <IN NAME="lat" TYPE="real" DIM="(:,:)">
!     The latitudes for the local domain.
!   </IN>
!   <IN NAME="pwt" TYPE="real" DIM="(:,:,:)">
!     Pressure weighting (air mass) for each layer (kg/m2)
!   </IN>
!   <IN NAME="r" TYPE="real" DIM="(:,:,:,:)">
!     Tracer mixing ratios (regional_tracer tracers in VMR)
!   </IN>
!   <IN NAME="Time, Time_next" TYPE="time_type">
!     Model time
!   </IN>
!   <IN NAME="phalf" TYPE="real" DIM="(:,:,:)">
!     Pressure on the model half levels (Pa)
!   </IN>
!   <IN NAME="is, js" TYPE="integer">
!     Local domain start indices
!   </IN>
!   <OUT NAME="chem_dt" TYPE="real" DIM="(:,:,:,:)">
!     Tracer tendencies from tropospheric chemistry (VMR/s)
!   </OUT>
!   <IN NAME="kbot" TYPE="integer, optional" DIM="(:,:)">
!     Integer array describing which model layer intercepts the surface.
!   </IN>

subroutine regional_tracer_driver( lon, lat, pwt, r, chem_dt, &
                                   phalf, is, js, kbot)

!-----------------------------------------------------------------------
   real, intent(in),    dimension(:,:)            :: lon, lat
   real, intent(in),    dimension(:,:,:)          :: pwt
   real, intent(in),    dimension(:,:,:,:)        :: r
   real, intent(out),   dimension(:,:,:,:)        :: chem_dt
   integer, intent(in)                            :: is,js
   real, intent(in),    dimension(:,:,:)          :: phalf
   integer, intent(in),  dimension(:,:), optional :: kbot
!-----------------------------------------------------------------------
   real, dimension(size(r,1),size(r,2)) :: emis
   real, dimension(size(r,1),size(r,2),ntracers) :: emis_save
   integer :: i,j,k,n,kb,id,jd,kd,trind
   integer :: nt, ntp
   logical :: used
   real, dimension(size(r,1),size(r,2),size(r,3),ntracers) :: emis_source
   character(len=64) :: name
   real, dimension(size(r,1),size(r,2),size(r,3),ntracers) :: prod, loss

   integer :: omp_get_num_threads
!-----------------------------------------------------------------------

!<ERROR MSG="regional_tracer_driver_init must be called first." STATUS="FATAL">
!   Regional_tracer_driver_init needs to be called before regional_tracer_driver.
!</ERROR>
   if (.not. module_is_initialized)  &
      call errmsg ('Regional_tracer_driver','regional_tracer_driver_init must be called first.', .true.)

   id=size(r,1); jd=size(r,2); kd=size(r,3)
   nt = size(r,4); ntp = size(chem_dt,4)
 
   emis_source(:,:,:,:) = 0.0  
! JianHe: source is done in other place, we only calculate loss here
! we only have emissions for e90
   chem_dt(:,:,:,:) = 0.0

   do n = 1, ntracers
      trind = tracer_indices(n)

!-----------------------------------------------------------------------
!     ... calculate chemical losses
!-----------------------------------------------------------------------
      if (trind > 0) then
         if( n==id_tch4 .and. tch4_vmr>0.) then
            loss(:,:,:,n) = loss_freq(n) * tch4_vmr
         else
            loss(:,:,:,n) = loss_freq(n) * r(:,:,:,trind)
         end if
      else
         loss(:,:,:,n) = 0.
      end if

   end do

!-----------------------------------------------------------------------
!     ... calculate chemical production
!-----------------------------------------------------------------------
   prod(:,:,:,:) = 0.
   if (id_tch4>0 .and. id_cofromch4>0) then
      prod(:,:,:,id_cofromch4)  = co_yield_from_ch4  * loss(:,:,:,id_tch4)
   end if
   if (id_avoc>0 .and. id_cofromavoc>0) then
      prod(:,:,:,id_cofromavoc) = co_yield_from_avoc * loss(:,:,:,id_avoc)
   end if
   if (id_bvoc>0 .and. id_cofrombvoc>0) then
      prod(:,:,:,id_cofrombvoc) = co_yield_from_bvoc * loss(:,:,:,id_bvoc)
   end if

!-----------------------------------------------------------------------
!     ... modify prod and loss for NH CMIP tracers
!-----------------------------------------------------------------------
   if (id_aoanh > 0) then
      trind = tracer_indices(id_aoanh)
      if (trind > 0) then
         prod(:,:,:,id_aoanh) = k_aging
         where (lat(:,:)>=lat30 .and. lat(:,:)<=lat50)
            loss(:,:,kd,id_aoanh) = k_relax_aoanh * r(:,:,kd,trind)
            prod(:,:,kd,id_aoanh) = 0.
         endwhere
      end if
   end if
   if (id_nh50 > 0) then
      trind = tracer_indices(id_nh50)
      if (trind > 0) then
         where (lat(:,:)>=lat30 .and. lat(:,:)<=lat50)
            prod(:,:,kd,id_nh50) = k_relax_nh50 * (nh50_fixed_val-r(:,:,kd,trind))
         endwhere
      end if
   end if


   do n = 1, ntracers
      trind = tracer_indices(n)

!-----------------------------------------------------------------------
!     ... compute tendency
!-----------------------------------------------------------------------
      if (trind>0 .and. trind<=ntp) then
         chem_dt(:,:,:,trind) = emis_source(:,:,:,n) + prod(:,:,:,n) - loss(:,:,:,n)
      end if

   end do     
   
!-----------------------------------------------------------------------
    
end subroutine regional_tracer_driver
!</SUBROUTINE>

!#######################################################################

! <SUBROUTINE NAME="regional_tracer_driver_init">
!   <OVERVIEW>
!     Initializes the regional tracer driver.
!   </OVERVIEW>
!   <DESCRIPTION>
!     This subroutine initializes the regional tracer module.
!     It is called from atmos_tracer_driver_init.
!     Data sets are read in for surface emissions.
!   </DESCRIPTION>
!   <TEMPLATE>
!       call regional_tracer_driver_init (lonb_mod, latb_mod, axes, Time, mask)
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

subroutine regional_tracer_driver_init( me, master, tracer_names)

!-----------------------------------------------------------------------
!
!   mask = optional mask (0. or 1.) that designates which grid points
!          are above (=1.) or below (=0.) the ground dimensioned as
!          (nlon,nlat,nlev).
!
!-----------------------------------------------------------------------
   integer, intent (in) :: me
   integer, intent (in) :: master

   character(len=32),   intent(in)                :: tracer_names(:)

   integer :: n, trind
   character(len=64) :: trname
   character(len=64) :: diag_name
   integer :: ierr, io
   integer :: unit

   integer :: omp_get_num_threads


!-----------------------------------------------------------------------
!
   if (module_is_initialized) return

!-----------------------------------------------------------------------
!     ... write version number
!-----------------------------------------------------------------------
!   call write_version_number(version, tagname)
    
!-----------------------------------------------------------------------
!     ... set initial value of indices
!-----------------------------------------------------------------------
   tracer_indices(:) = 0
   do n=1,ntracers
      trind = get_tracer_ndx(tracer_names, regtrc_names(n))
      if (trind >0) then
         tracer_indices(n) = trind
         if ( me == master) write(*,30) regtrc_names(n),tracer_indices(n)
      else
!<ERROR MSG="Tropospheric chemistry tracer not found in field table" STATUS="WARNING">
!   A tropospheric chemistry tracer was not included in the field table
!</ERROR>
         call errmsg ('regional_tracer_driver_init', trim(regtrc_names(n)) // ' is not found', .true.)
      end if
   end do
30 format (A,' was initialized as tracer number ',i3)

!-----------------------------------------------------------------------
!     ... Keep track of tracer indices
!-----------------------------------------------------------------------
   do n=1,ntracers
      trname = regtrc_names(n)
      select case (trname)
         case('tch4')
            id_tch4 = n
         case('avoc')
            id_avoc = n
         case('bvoc')
            id_bvoc = n
         case('cofromch4')
            id_cofromch4 = n
         case('cofromavoc')
            id_cofromavoc = n
         case('cofrombvoc')
            id_cofrombvoc = n
         case('aoanh')
            id_aoanh = n
         case('nh50')
            id_nh50 = n
      end select
   end do

   module_is_initialized = .true.
      
      
!-----------------------------------------------------------------------
      
end subroutine regional_tracer_driver_init
!</SUBROUTINE>



!############################################################################
end module atmos_regional_tracer_driver_mod
