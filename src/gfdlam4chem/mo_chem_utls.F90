      module mo_chem_utls_mod

implicit none
      private
      public :: adjh2o, inti_mr_xform, mmr2vmr, vmr2mmr, negtrc, &
                get_spc_ndx, get_het_ndx, get_extfrc_ndx, &
                has_drydep, has_srfems, get_rxt_ndx, get_grp_ndx, &
                get_grp_mem_ndx, chem_utls_init, &
                get_solar_flux_by_band

      integer, parameter :: NO_TRACER = -99

      public :: get_tracer_ndx, NO_TRACER, sjl_fillz

!     save

      integer :: ox_ndx, o3_ndx, o1d_ndx, o_ndx
      logical :: do_ox

character(len=128), parameter :: version     = '$Id$'
character(len=128), parameter :: tagname     = '$Name$'
logical                       :: module_is_initialized = .false.

      contains

      subroutine chem_utls_init (retain_cm3_bugs)
!-----------------------------------------------------------------------
!     ... Initialize the chem utils module
!-----------------------------------------------------------------------

      implicit none

      logical, intent(in) :: retain_cm3_bugs

!      ox_ndx = get_spc_ndx( 'OX' )
! for ox budget (jmao,1/7/2011)
   if (retain_cm3_bugs) then
      ox_ndx = get_spc_ndx( 'OX' )
   else
      ox_ndx = get_spc_ndx( 'O3' )
   endif
      if( ox_ndx > 0 ) then
         o3_ndx  = get_grp_mem_ndx( 'O3' )
         o1d_ndx = get_grp_mem_ndx( 'O1D' )
         o_ndx   = get_grp_mem_ndx( 'O' )
         do_ox   = o3_ndx > 0 .and. o1d_ndx > 0 .and. o_ndx > 0
      else
         o3_ndx  = 1
         o1d_ndx = 1
         o_ndx   = 1
         do_ox = .false.
      end if

      end subroutine chem_utls_init

      subroutine adjh2o( h2o, sh, mbar, vmr, do_interactive_h2o, plonl )
!-----------------------------------------------------------------------
!     ... transform water vapor from mass to volumetric mixing ratio
!-----------------------------------------------------------------------


      implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in) :: plonl
      real, dimension(:,:,:), intent(in)  :: vmr         ! xported species vmr
      real, dimension(:,:),   intent(in)  :: sh          ! specific humidity ( mmr )
      real, dimension(:,:),   intent(in)  :: mbar        ! atmos mean mass
      logical,                intent(in)  :: do_interactive_h2o ! include h2o sources/sinks?
      real, dimension(:,:),   intent(out) :: h2o         ! water vapor vmr

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      real, parameter :: mh2o = 1. /18.01528

      integer ::   k, ndx_ch4
      real    ::   t_value(plonl)
      integer ::   plev

      plev = SIZE(vmr,2)
!-----------------------------------------------------------------------
!       ... if not using interactive water vapor, adjust model
!           water vapor in stratosphere for source from CH4 oxidation
!-----------------------------------------------------------------------
      ndx_ch4 = get_spc_ndx( 'CH4' )
!++lwh
         do k = 1,plev
            h2o(:,k)   = mbar(:,k) * sh(:plonl,k) * mh2o
            if( .not. do_interactive_h2o .and. ndx_ch4 > 0 ) then
               t_value(:) = 6.e-6 - 2.*vmr(:,k,ndx_ch4)
!              where( t_value(:) > h2o(:,k) )
!                 h2o(:,k) = t_value(:)
!              end where
               h2o(:,k) = MAX(h2o(:,k),t_value(:))
            end if
         end do
!--lwh

      end subroutine adjh2o      

      subroutine inti_mr_xform( sh, mbar, plonl )
!-----------------------------------------------------------------
!       ... initialize mean atmospheric "wet" mass
!-----------------------------------------------------------------


      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)  :: sh(:,:)     ! specific humidity (kg/kg)
      real, intent(out) :: mbar(:,:)   ! mean wet atm mass ( amu )

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
      real, parameter :: dry_mass = 28.966    ! amu
      real, parameter :: mfac = 1. / .622

      integer :: k
      integer :: plev
      
      plev = size(sh,2)

      do k = 1,plev
         mbar(:,k) = dry_mass
      end do

      end subroutine inti_mr_xform

      subroutine mmr2vmr( vmr, mmr, mbar, plonl )
!-----------------------------------------------------------------
!       ... xfrom from mass to volume mixing ratio
!-----------------------------------------------------------------

      use chem_mods_mod, only : adv_mass
      use mo_grid_mod,   only : pcnstm1, pcnst

      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)  :: mbar(:,:)
      real, intent(in)  :: mmr(:,:,:)
      real, intent(out) :: vmr(:,:,:)

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
      integer :: k, m
      integer :: plev
      
      plev = size(mbar,2)

      do m = 1,pcnstm1
         if( adv_mass(m) /= 0. ) then
            do k = 1,plev
               vmr(:,k,m) = mbar(:,k) * mmr(:,k,m) / adv_mass(m)
            end do
         end if
      end do

      end subroutine mmr2vmr

      subroutine vmr2mmr( vmr, mmr, nas, grp_ratios, mbar, plonl )
!-----------------------------------------------------------------
!       ... xfrom from mass to volume mixing ratio
!-----------------------------------------------------------------

      use chem_mods_mod, only : adv_mass, nadv_mass, grpcnt
      use mo_grid_mod,   only : pcnstm1, pcnst

      implicit none

!-----------------------------------------------------------------
!       ... dummy args
!-----------------------------------------------------------------
      integer, intent(in) :: plonl
      real, intent(in)  :: mbar(:,:)
      real, intent(in)  :: vmr(:,:,:)
      real, intent(out) :: mmr(:,:,:)
      real, intent(in)  :: grp_ratios(:,:,:)
      real, intent(out) :: nas(:,:,:)

!-----------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------
      integer :: k, m
      integer :: plev
      real    :: grp_mass(plonl)            ! weighted group mass

      plev = size(mbar,2)

!-----------------------------------------------------------------
!       ... the non-group species
!-----------------------------------------------------------------
      do m = 1,pcnstm1
         if( adv_mass(m) /= 0. ) then
            do k = 1,plev
               mmr(:,k,m) = adv_mass(m) * vmr(:,k,m) / mbar(:,k)
            end do
         end if
      end do
!-----------------------------------------------------------------
!       ... the "group" species
!-----------------------------------------------------------------
      if( do_ox ) then
         do k = 1,plev
            grp_mass(:)     = grp_ratios(:,k,o3_ndx) * nadv_mass(o3_ndx) &
                              + grp_ratios(:,k,o_ndx) * nadv_mass(o_ndx) &
                              + grp_ratios(:,k,o1d_ndx) * nadv_mass(o1d_ndx)      
            mmr(:,k,ox_ndx)  = grp_mass(:) * vmr(:,k,ox_ndx) / mbar(:,k)
            grp_mass(:)     = mmr(:,k,ox_ndx) / grp_mass(:)
            nas(:,k,o3_ndx)  = nadv_mass(o3_ndx) * grp_ratios(:,k,o3_ndx) * grp_mass(:)
            nas(:,k,o_ndx)   = nadv_mass(o_ndx) * grp_ratios(:,k,o_ndx) * grp_mass(:)
            nas(:,k,o1d_ndx) = nadv_mass(o1d_ndx) * grp_ratios(:,k,o1d_ndx) * grp_mass(:)
         end do
      end if

      end subroutine vmr2mmr

      subroutine negtrc( lat, header, fld, plonl )
!-----------------------------------------------------------------------
!       ... check for negative constituent values and
!           replace with zero value
!-----------------------------------------------------------------------

      use mo_grid_mod,    only : pcnstm1
      use m_tracname_mod, only : tracnam

      implicit none

!-----------------------------------------------------------------------
!       ... dummy arguments
!-----------------------------------------------------------------------
      integer, intent(in)          :: lat                      ! current latitude
      integer, intent(in)          :: plonl
      character(len=*), intent(in) :: header                   ! caller tag
      real, intent(inout)          :: fld(:,:,:)               ! field to check

!-----------------------------------------------------------------------
!       ... local variables
!-----------------------------------------------------------------------
      integer :: m
      integer :: nneg                       ! flag counter

      do m  = 1,pcnstm1
         nneg = count( fld(:,:,m) < 0. )
         if( nneg > 0 ) then
            where( fld(:,:,m) < 0. )
               fld(:,:,m) = 0.
            endwhere
!           if( pdiags%negtrc ) then
!              worst     = minval( fld(:,:,m) )
!              windex(:) = minloc( fld(:,:,m) )
!              iw        = windex(1)
!              kw        = windex(2)
!           end if
         end if
!        if( pdiags%negtrc .and. nneg > 0 ) then
!           write(*,*) header(:len(header)), tracnam(m), ' has ',nneg,' neg values'
!           write(*,*) ' worst =',worst,' @ long = ',iw,' lat = ',lat,' eta = ',kw
!        end if
      end do

      end subroutine negtrc

      integer function get_spc_ndx( spc_name )
!-----------------------------------------------------------------------
!     ... return overall species index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : pcnstm1
      use m_tracname_mod, only : tracnam

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_spc_ndx = -1
      do m = 1,pcnstm1
         if( trim( spc_name ) == trim( tracnam(m) ) ) then
            get_spc_ndx = m
            exit
         end if
      end do

      end function get_spc_ndx

      integer function get_grp_ndx( grp_name )
!-----------------------------------------------------------------------
!     ... return group index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : ngrp, grp_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: grp_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_grp_ndx = -1
      do m = 1,ngrp
         if( trim( grp_name ) == trim( grp_lst(m) ) ) then
            get_grp_ndx = m
            exit
         end if
      end do

      end function get_grp_ndx

      integer function get_grp_mem_ndx( mem_name )
!-----------------------------------------------------------------------
!     ... return group member index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : grpcnt
      use m_tracname_mod, only : natsnam

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: mem_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_grp_mem_ndx = -1
      if( grpcnt > 0 ) then
         do m = 1,max(1,grpcnt)
            if( trim( mem_name ) == trim( natsnam(m) ) ) then
               get_grp_mem_ndx = m
               exit
            end if
         end do
      end if

      end function get_grp_mem_ndx

      integer function get_het_ndx( het_name )
!-----------------------------------------------------------------------
!     ... return overall het process index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : hetcnt, het_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: het_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_het_ndx = -1
      do m = 1,max(1,hetcnt)
         if( trim( het_name ) == trim( het_lst(m) ) ) then
            get_het_ndx = m
            exit
         end if
      end do

      end function get_het_ndx

      integer function get_extfrc_ndx( frc_name )
!-----------------------------------------------------------------------
!     ... return overall external frcing index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : extcnt, extfrc_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: frc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_extfrc_ndx = -1
      if( extcnt > 0 ) then
         do m = 1,max(1,extcnt)
            if( trim( frc_name ) == trim( extfrc_lst(m) ) ) then
               get_extfrc_ndx = m
               exit
            end if
         end do
      end if

      end function get_extfrc_ndx

      integer function get_rxt_ndx( rxt_alias )
!-----------------------------------------------------------------------
!     ... return overall external frcing index associated with spc_name
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : rxt_alias_cnt, rxt_alias_lst, rxt_alias_map

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: rxt_alias

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      get_rxt_ndx = -1
      do m = 1,rxt_alias_cnt
         if( trim( rxt_alias ) == trim( rxt_alias_lst(m) ) ) then
            get_rxt_ndx = rxt_alias_map(m)
            exit
         end if
      end do

      end function get_rxt_ndx

      logical function has_drydep( spc_name )
!-----------------------------------------------------------------------
!     ... return logical for species dry deposition
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : drydep_cnt, drydep_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      has_drydep = .false.
      do m = 1,drydep_cnt
         if( trim( spc_name ) == trim( drydep_lst(m) ) ) then
            has_drydep = .true.
            exit
         end if
      end do

      end function has_drydep

      logical function has_srfems( spc_name )
!-----------------------------------------------------------------------
!     ... return logical for species surface emission
!-----------------------------------------------------------------------

      use chem_mods_mod,  only : srfems_cnt, srfems_lst

      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
      character(len=*), intent(in) :: spc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
      integer :: m

      has_srfems = .false.
      do m = 1,srfems_cnt
         if( trim( spc_name ) == trim( srfems_lst(m) ) ) then
            has_srfems = .true.
            exit
         end if
      end do

      end function has_srfems

      integer function get_tracer_ndx( tracer_names, trc_name )
!-----------------------------------------------------------------------
!     ... return overall tracer index associated with tracer_name
!-----------------------------------------------------------------------

      use gfdl_time_utls_mod, only : lowercase
      implicit none

!-----------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------
    character(len=32), intent(in) :: tracer_names(:)
    character(len=*),  intent(in) :: trc_name

!-----------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------
    integer :: i

    get_tracer_ndx = NO_TRACER

    do i=1, size(tracer_names)
       if (trim(lowercase(trc_name)) == trim(lowercase(tracer_names(i)))) then
           get_tracer_ndx = i
           exit
       endif
    enddo

    end function get_tracer_ndx

! ######################################################################
!

subroutine GET_RH (T,Q,P,RH,mask)
  !***********************************************************************
  !  SUBROUTINE GET_RH
  !  PURPOSE
  !     VECTOR COMPUTATION OF RELATIVE HUMIDITY
  !  DESCRIPTION OF PARAMETERS
  !     T        TEMPERATURE VECTOR (DEG K)
  !     P        PRESSURE VECTOR (Pa)
  !     Q        SPECIFIC HUMIDITY (kg/kg)
  !     RH       RELATIVE HUMIDITY
  !     MASK     EARTH SURFACE BELOW GROUND (in pressure coordinates)
  !
  !***********************************************************************
  !

 Real, parameter :: ONE    = 1.
 Real, parameter :: ZP622  = 0.622
 Real, parameter :: Z1P0S1 = 1.00001
 Real, parameter :: Z1P622 = 1.622
 Real, parameter :: Z138P9 = 138.90001
 Real, parameter :: Z198P9 = 198.99999
 Real, parameter :: Z200   = 200.0
 Real, parameter :: Z337P9 = 337.9

 real, intent(in), dimension(:,:,:) :: T !temp at curr time step [ deg k ]
 real, intent(in), dimension(:,:,:) :: P !pressure at full levels [ Pa ]
 real, intent(in), dimension(:,:,:) :: Q !specific humidity at current time step [ kg / kg ]
 real, intent(in), dimension(:,:,:), optional :: mask !
 real, intent(out), dimension(:,:,:) :: RH !relative humidity [0-1]

 ! Dynamic Work Space
 ! ------------------
 real :: A1622
 real, dimension(size(q,1),size(q,2),size(q,3)) :: e1, e2, tq, qs
 integer, dimension(size(q,1),size(q,2),size(q,3)) :: i1, i2
 integer :: i,j,k

 !
 real, dimension(67) :: EST1
 data EST1/       0.31195E-02, 0.36135E-02, 0.41800E-02, &
      0.48227E-02, 0.55571E-02, 0.63934E-02, 0.73433E-02, &
      0.84286E-02, 0.96407E-02, 0.11014E-01, 0.12582E-01, &
      0.14353E-01, 0.16341E-01, 0.18574E-01, 0.21095E-01, &
      0.23926E-01, 0.27096E-01, 0.30652E-01, 0.34629E-01, &
      0.39073E-01, 0.44028E-01, 0.49546E-01, 0.55691E-01, &
      0.62508E-01, 0.70077E-01, 0.78700E-01, 0.88128E-01, &
      0.98477E-01, 0.10983E+00, 0.12233E+00, 0.13608E+00, &
      0.15121E+00, 0.16784E+00, 0.18615E+00, 0.20627E+00, &
      0.22837E+00, 0.25263E+00, 0.27923E+00, 0.30838E+00, &
      0.34030E+00, 0.37520E+00, 0.41334E+00, 0.45497E+00, &
      0.50037E+00, 0.54984E+00, 0.60369E+00, 0.66225E+00, &
      0.72589E+00, 0.79497E+00, 0.86991E+00, 0.95113E+00, &
      0.10391E+01, 0.11343E+01, 0.12372E+01, 0.13484E+01, &
      0.14684E+01, 0.15979E+01, 0.17375E+01, 0.18879E+01, &
      0.20499E+01, 0.22241E+01, 0.24113E+01, 0.26126E+01, &
      0.28286E+01, 0.30604E+01, 0.33091E+01, 0.35755E+01/
 !
 real, dimension(72)  :: EST2
 data EST2/ &
      0.38608E+01, 0.41663E+01, 0.44930E+01, 0.48423E+01, &
      0.52155E+01, 0.56140E+01, 0.60394E+01, 0.64930E+01, &
      0.69767E+01, 0.74919E+01, 0.80406E+01, 0.86246E+01, &
      0.92457E+01, 0.99061E+01, 0.10608E+02, 0.11353E+02, &
      0.12144E+02, 0.12983E+02, 0.13873E+02, 0.14816E+02, &
      0.15815E+02, 0.16872E+02, 0.17992E+02, 0.19176E+02, &
      0.20428E+02, 0.21750E+02, 0.23148E+02, 0.24623E+02, &
      0.26180E+02, 0.27822E+02, 0.29553E+02, 0.31378E+02, &
      0.33300E+02, 0.35324E+02, 0.37454E+02, 0.39696E+02, &
      0.42053E+02, 0.44531E+02, 0.47134E+02, 0.49869E+02, &
      0.52741E+02, 0.55754E+02, 0.58916E+02, 0.62232E+02, &
      0.65708E+02, 0.69351E+02, 0.73168E+02, 0.77164E+02, &
      0.81348E+02, 0.85725E+02, 0.90305E+02, 0.95094E+02, &
      0.10010E+03, 0.10533E+03, 0.11080E+03, 0.11650E+03, &
      0.12246E+03, 0.12868E+03, 0.13517E+03, 0.14193E+03, &
      0.14899E+03, 0.15634E+03, 0.16400E+03, 0.17199E+03, &
      0.18030E+03, 0.18895E+03, 0.19796E+03, 0.20733E+03, &
      0.21708E+03, 0.22722E+03, 0.23776E+03, 0.24871E+03/
 !
 real, dimension(139) :: EST
 EQUIVALENCE (EST(1)  , EST1(1)), (EST(68),EST2(1))
 !***********************************************************************
 !
 A1622   = ONE  / Z1P622
 TQ = T - Z198P9
 I1(:,:,:) = 1
 I2(:,:,:) = 1
 where ( T < Z200 ) TQ = Z1P0S1
 where ( T > Z337P9 ) TQ = Z138P9
 IF (present(mask)) THEN
    where ( mask > 0. )
       I1 = int(TQ)
       I2 = I1 + 1
    end where
 else
    I1 = int(TQ)
    I2 = I1 + 1
 endif
 do i=1,size(q,1)
    do j=1,size(q,2)
       do k=1,size(q,3)
          E1(i,j,k) =  EST( I1(i,j,k) )
          E2(i,j,k) =  EST( I2(i,j,k) )
       enddo
    enddo
 enddo
 QS(:,:,:) = TQ(:,:,:) - float(I1(:,:,:))
 QS(:,:,:) = E1(:,:,:) + QS(:,:,:) * ( E2(:,:,:)-E1(:,:,:) )
 E1(:,:,:) = (0.01 * P(:,:,:)) * A1622
 where ( E1 < QS ) QS = E1
 if (present(mask)) then
    where ( mask > 0. )  QS = ZP622 * QS / ( P * 0.01)
 else
    QS(:,:,:) = ZP622 * QS(:,:,:) / ( P(:,:,:) * 0.01)
 endif
 RH(:,:,:) = Q(:,:,:)/QS(:,:,:)

end subroutine GET_RH
! ######################################################################
!
subroutine get_cldf(ps, pfull, rh, cldf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This subroutine estimates the cloud fraction "cldf" for
!!! each grid box using an empirical function of the relative
!!! humidity in that grid box, following Sundqvist et al., Mon. Weather Rev.,
!!! v117, 164101657, 1989:
!!!
!!!             cldf = 1 - sqrt[ 1 - (RH - RH0)/(1 - RH0) ]
!!!
!!! where RH is the relative humidity and RH0 is the threshold relative
!!! humidity for condensation specified as a function of pressure based
!!! on Xu and Krueger, Mon. Weather Rev., v119, 342-367, 1991.
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real, intent(in),  dimension(:,:)   :: ps
 real, intent(in),  dimension(:,:,:) :: pfull
 real, intent(in),  dimension(:,:,:) :: rh
 real, intent(out), dimension(:,:,:) :: cldf
 real, parameter :: zrt = 0.6
 real, parameter :: zrs = 0.99
 integer           :: i,j,k, id, jd, kd
 real              :: p, r, r0, b0
 id=size(pfull,1); jd=size(pfull,2); kd=size(pfull,3)

 do k = 1, kd
    do j = 1, jd
       do i = 1, id
          P = pfull(i,j,k)
          R = RH(i,j,k)
          R0 = ZRT + (ZRS-ZRT) * exp(1.-(PS(i,j)/P)**2.5)
          B0 = (R-R0) / (1.-R0)
          if (R .lt.R0) B0 = 0.
          if (B0.gt.1.) B0 = 1.
          CLDF(i,j,k) = 1.-sqrt(1.-B0)
       end do
    end do
 end do

end subroutine get_cldf

!######################################################################
! !IROUTINE: sjl_fillz --- Fill from neighbors below and above
!
! !INTERFACE:
subroutine sjl_fillz(im, km, nq, q, dp)

 implicit none

 ! !INPUT PARAMETERS:
 integer,  intent(in):: im                ! No. of longitudes
 integer,  intent(in):: km                ! No. of levels
 integer,  intent(in):: nq                ! Total number of tracers

 real, intent(in)::  dp(im,km)       ! pressure thickness
 ! !INPUT/OUTPUT PARAMETERS:
 real, intent(inout) :: q(im,km,nq)   ! tracer mixing ratio

 ! !DESCRIPTION:
 !   Check for "bad" data and fill from east and west neighbors
 !
 ! !BUGS:
 !   Currently this routine only performs the east-west fill algorithm.
 !   This is because the N-S fill is very hard to do in a reproducible
 !   fashion when the problem is decomposed by latitudes.
 !
 ! !REVISION HISTORY:
 !   00.04.01   Lin        Creation
 !
 !EOP
 !-----------------------------------------------------------------------
 !BOC
 !
 ! !LOCAL VARIABLES:
 integer i, k, ic
 real qup, qly, dup

 do ic=1,nq
    ! Top layer
    do i=1,im
       if( q(i,1,ic) < 0.) then
          q(i,2,ic) = q(i,2,ic) + q(i,1,ic)*dp(i,1)/dp(i,2)
          q(i,1,ic) = 0.
       endif
    enddo

    ! Interior
    do k=2,km-1
       do i=1,im
          if( q(i,k,ic) < 0. ) then
             ! Borrow from above
             qup =  q(i,k-1,ic)*dp(i,k-1)
             qly = -q(i,k  ,ic)*dp(i,k  )
             dup =  min( 0.75*qly, qup )        !borrow no more than 75% from top
             q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
             ! Borrow from below: q(i,k,ic) is still negative at this stage
             q(i,k+1,ic) = q(i,k+1,ic) + (dup-qly)/dp(i,k+1)
             q(i,k  ,ic) = 0.
          endif
       enddo
    enddo

    ! Bottom layer
    k = km
    do i=1,im
       if( q(i,k,ic) < 0.) then
          ! Borrow from above
          qup =  q(i,k-1,ic)*dp(i,k-1)
          qly = -q(i,k  ,ic)*dp(i,k  )
          dup =  min( qly, qup )
          q(i,k-1,ic) = q(i,k-1,ic) - dup/dp(i,k-1)
          q(i,k,ic) = 0.
       endif
    enddo
 enddo
end subroutine sjl_fillz

!=======================================================================

subroutine get_solar_flux_by_band (Time,solarflux, ref)

      use mo_fastjx_mod, only : get_solar_data
      use gfdl_time_utls_mod, only : time_type, get_date

      implicit none

real,    intent(out)          :: solarflux(:)
logical, intent(in), optional :: ref
type(time_type), intent(in)   :: Time
!---------------------------------------------------------------------
integer, parameter :: nbands = 18
     
real,    dimension(:), allocatable  :: solflxbandref, solflxband
real  :: solar_constant
integer :: yr, mo, day, hr, minute, sec
integer :: nband, iounit

allocate ( solflxbandref  (nbands) )
allocate ( solflxband (nbands) )

iounit = 66

      if (present(ref)) then  ! no time-varying solar
         if (ref) then
            open(iounit, file='INPUT/esf_sw_input_data', form='formatted', status='OLD')
            read(iounit,'(12f10.4)') ( solflxbandref(nband), nband=1,nbands )
            solarflux = solflxbandref
            close(iounit)
            return
         endif
      else ! time-varying solar
         call get_date(Time,yr,mo,day,hr,minute,sec) 
         call get_solar_data(yr, mo, solar_constant, solflxband)
         solarflux = solflxband
      endif

!---------------------------------------------------------------------
end subroutine get_solar_flux_by_band

      end module mo_chem_utls_mod
