module atmos_tropopause_mod
! <CONTACT EMAIL="Larry.Horowitz@noaa.gov">
!   Larry.Horowitz
! </CONTACT>

! <OVERVIEW>
!   Tropopause diagnostics
! </OVERVIEW>

! <DESCRIPTION>
!   Tropopause diagnostics
! </DESCRIPTION>

use gfdl_time_utls_mod, only : time_type

use mo_errmsg,             only : errmsg

implicit none

private

!-----------------------------------------------------------------------
!----- interfaces -------

public  atmos_tropopause, &
        atmos_tropopause_init

!-----------------------------------------------------------------------
!----------- module data -------------------
!-----------------------------------------------------------------------

character(len=48), parameter :: module_name = 'tracers'

logical :: module_is_initialized = .FALSE.
integer :: logunit

real, parameter    :: ztrop_low = 5.e3   ! lowest tropopause level allowed (m)
real, parameter    :: ztrop_high = 20.e3 ! highest tropopause level allowed (m)
real, parameter    :: max_dtdz   = 2.e-3  ! max dt/dz for tropopause level (K/m)

!-----------------------------------------------------------------------
!     ... identification numbers for diagnostic fields
!-----------------------------------------------------------------------
integer :: id_ptp, id_tatp, id_ztp

!---------------------------------------------------------------------
!-------- namelist  ---------
!-----------------------------------------------------------------------

logical  :: do_tropopause_diagnostics = .true.

!namelist /atmos_tropopause_nml/  &
!          do_tropopause_diagnostics

!---- version number -----
character(len=128) :: version = '$$'
character(len=128) :: tagname = '$$'

!-----------------------------------------------------------------------

contains

!#######################################################################

!<SUBROUTINE NAME ="atmos_tropopause">
!<OVERVIEW>
!  A subroutine to calculate tropopause diagnostics.
!
! do_co2_restore   = logical to turn co2_restore on/off: default = .false.
! restore_co2_dvmr = partial pressure of co2 to which to restore  (mol/mol)
! restore_klimit   = atmospheric level to which to restore starting from top
! restore_tscale   = timescale in seconds with which to restore
!
!</OVERVIEW>
!<DESCRIPTION>
! A routine to calculate tropopause diagnostics.
!</DESCRIPTION>
!<TEMPLATE>
!call atmos_tropopause (Time, Time_next, t, pfull)
!</TEMPLATE>
!
!   <IN NAME="Time" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="Time_next" TYPE="type(time_type)">
!     Model time.
!   </IN>
!   <IN NAME="t" TYPE="real" DIM="(:,:,:)">
!     Temperature.
!   </IN>
!   <IN NAME="pfull" TYPE="real" DIM="(:,:,:)">
!     Pressures on the model full levels.
!   </IN>
!   <IN NAME="z_full" TYPE="real" DIM="(:,:,:)">
!     Height of the model full levels.
!   </IN>

subroutine atmos_tropopause(is, ie, js, je, Time, Time_next, t, pfull, z_full, &
                            ptp, tropopause_ind)

   integer, intent(in)                   :: is, ie, js, je
   type (time_type), intent(in)          :: Time, Time_next
   real,    intent(in), dimension(:,:,:) :: t            ! K
   real,    intent(in), dimension(:,:,:) :: pfull        ! Pa
   real,    intent(in), dimension(:,:,:) :: z_full       ! m
   integer, intent(out), dimension(:,:)  :: tropopause_ind ! 1
   real,    intent(inout), dimension(:,:)  :: ptp

!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

integer   :: i,j,k,id,jd,kd
real      :: dtemp

real, dimension(size(t,1),size(t,2)) :: tatp, ztp

!-----------------------------------------------------------------------

    if (.not. module_is_initialized)  &
       call errmsg ('Atmos_tropopause','atmos_tropopause_init must be called first.', .true.)

    id=size(t,1); jd=size(t,2); kd=size(t,3)

    do j=1,jd
    do i=1,id

       do k = kd-1,2,-1
          if (z_full(i,j,k) < ztrop_low ) then
              cycle
          else if( z_full(i,j,k) > ztrop_high ) then
              tropopause_ind(i,j)    = k
              exit
          end if
          dtemp = t(i,j,k) - t(i,j,k-1)
          if( dtemp < max_dtdz*(z_full(i,j,k-1) - z_full(i,j,k)) ) then
             tropopause_ind(i,j)    = k
             exit
          end if
       end do

       ptp(i,j) = pfull(i,j,tropopause_ind(i,j))
       tatp(i,j) = t(i,j,tropopause_ind(i,j))
       ztp(i,j) = z_full(i,j,tropopause_ind(i,j))

    enddo
    enddo


end subroutine atmos_tropopause
!</SUBROUTINE >



!#######################################################################

!<SUBROUTINE NAME ="atmos_tropopause_init">

!<OVERVIEW>
! Subroutine to initialize the tropopause diagnostics module.
!</OVERVIEW>

 subroutine atmos_tropopause_init (me, master)

!
!-----------------------------------------------------------------------
!     arguments
!-----------------------------------------------------------------------
!
   integer, intent (in) :: me
   integer, intent (in) :: master

!type(time_type),       intent(in)                   :: Time
!
!-----------------------------------------------------------------------
!     local variables
!         unit       io unit number used to read namelist file
!         ierr       error code
!         io         error status returned from io operation
!-----------------------------------------------------------------------
!
integer :: ierr, unit, io
!
!-----------------------------------------------------------------------
!     local parameters
!-----------------------------------------------------------------------
!

    if (module_is_initialized) return

    !call write_version_number (version, tagname)

    if(me==master) then
        write(*,*) 'atmos_tropopause_init is done'
    end if
    module_is_initialized = .TRUE.


!-----------------------------------------------------------------------

end subroutine atmos_tropopause_init
!</SUBROUTINE>

end module atmos_tropopause_mod
