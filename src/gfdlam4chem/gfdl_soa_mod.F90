module gfdl_soa_mod
! <DESCRIPTION>
!   This module is an implementation of Secondary organic aerosols (SOA)
!   from anthropogenic activities, and is based on Tie et al. (JGR, 2003).
!   The only souce of SOA is due to the oxydation of C4H10 by OH.
!   The concentrations of these 2 gas species are read as input.
! </DESCRIPTION>
! <WARNING>
!  To save space only the actual month of input files are kept in memory. 
!  This implies that the "gfdl_soa_init" should be executed at the begining 
!  of each month. In other words, the script should not run more than 1 month
!  without a restart.
! </WARNING>
! <CONTACT EMAIL="Paul.Ginoux@noaa.gov">
!   Paul Ginoux
! </CONTACT>
!-----------------------------------------------------------------------

use           mo_chem_utls_mod, only : get_tracer_ndx, &
                                       NO_TRACER
use          catchem_constants, only : PI, GRAV, RDGAS, WTMAIR, AVOGNO
use      mo_errmsg,             only : errmsg

implicit none

private
!-----------------------------------------------------------------------
!----- interfaces -------
!
public  gfdl_soa_init, gfdl_soa_end, gfdl_soa_chem, &
        gfdl_soa_endts

!--- Arrays to help calculate tracer sources/sinks ---

character(len=6), parameter :: module_name = 'tracer'

integer :: nSOA    = 0  ! tracer number for Secondary Organic Aerosol
integer :: nOH     = 0  ! tracer number for OH
integer :: nC4H10  = 0  ! tracer number for C4H10
integer :: nISOP   = 0  ! tracer number for isoprene
integer :: nTERP   = 0  ! tracer number for terpenes (currently monoterpenes)
integer :: nO3     = 0  ! tracer number for ozone
integer :: nNO3    = 0  ! tracer number for NO3
!--- identification numbers for  diagnostic fields and axes ----
integer ::   id_OH_conc            = 0
integer ::   id_C4H10_conc         = 0
integer ::   id_ISOP_conc          = 0
integer ::   id_TERP_conc          = 0
integer ::   id_SOA_chem           = 0
integer ::   id_SOA_chem_col       = 0
integer ::   id_chepsoa            = 0 ! cmip field
integer ::   id_chepasoa           = 0 ! cmip field
integer ::   id_SOA_isoprene       = 0
integer ::   id_SOA_terpene        = 0
integer ::   id_SOA_biogenic       = 0

real, parameter       :: wtm_C = 12.
real, parameter       :: wtm_C4H10 = 58.
integer, parameter    :: carbon_per_isoprene = 5
integer, parameter    :: carbon_per_terpene = 10
real, parameter       :: cm2_per_m2 = 1.e4
real, parameter       :: kg_per_g = 1.e-3
real, parameter       :: om_oc_ratio = 1.5
! Factors to convert from molecules(BVOC)/cm2/s to kg(OM)/m2/s
real, parameter       :: isoprene_factor = cm2_per_m2 / AVOGNO * wtm_C * kg_per_g * carbon_per_isoprene * om_oc_ratio
real, parameter       :: terpene_factor = cm2_per_m2 / AVOGNO * wtm_C * kg_per_g * carbon_per_terpene * om_oc_ratio
integer, parameter    :: TS_CONSTANT=1, TS_FIXED=2, TS_VARYING=3
real, parameter       :: boltz = 1.38044e-16      ! Boltzmann's Constant (erg/K)

!-----------------------------------------------------------------------
!----------- namelist -------------------
!-----------------------------------------------------------------------

character(len=32)  :: gas_conc_filename = 'gas_conc_3D.nc'
character(len=32)  :: isoprene_filename = 'biogenic_emis.nc'
character(len=32)  :: terpene_filename  = 'biogenic_emis.nc'
character(len=32), dimension(2) :: gas_conc_name
data gas_conc_name/'OH','C4H10'/
character(len=32), dimension(1) :: isoprene_input_name = (/ 'isoprene' /)
character(len=32), dimension(1) :: terpene_input_name = (/ 'terpenes' /)
integer, dimension(6) :: isoprene_dataset_entry  = (/ 1, 1, 1, 0, 0, 0 /)
integer, dimension(6) :: terpene_dataset_entry   = (/ 1, 1, 1, 0, 0, 0 /)
character(len=80)     :: isoprene_source = 'guenther'
character(len=80)     :: terpene_source = 'guenther'
character(len=80)     :: isoprene_time_dependency_type = 'constant'
character(len=80)     :: terpene_time_dependency_type  = 'constant'
real                  :: isoprene_SOA_yield = 0.1
real                  :: terpene_SOA_yield  = 0.1

real                  :: ISOP_OH_SOA_yield = 0.05
real                  :: TERP_OH_SOA_yield  = 0.05
real                  :: ISOP_NO3_SOA_yield = 0.05
real                  :: TERP_NO3_SOA_yield  = 0.05


!logical               :: use_interactive_tracers = .true., &
!                         use_interactive_BVOC_emis = .false.

logical :: module_is_initialized=.FALSE.
logical :: used

!---- version number -----
character(len=128) :: version = '$Id$'
character(len=128) :: tagname = '$Name$'
!-----------------------------------------------------------------------

contains


!#######################################################################

!<SUBROUTINE NAME="gfdl_soa_init">
!<OVERVIEW>
! The constructor routine for the soa module.
!</OVERVIEW>
subroutine gfdl_soa_init (me, master, tracer_names, &
                          isoprene_SOA_yield,terpene_SOA_yield,&
                          use_interactive_BVOC_emis,mask)
!-----------------------------------------------------------------------
   integer, intent (in) :: me
   integer, intent (in) :: master
   character(len=32), intent(in) :: tracer_names(:)
   real, intent(in) :: isoprene_SOA_yield,terpene_SOA_yield
   logical , intent(in) :: use_interactive_BVOC_emis
   real, intent(in), dimension(:,:,:), optional        :: mask
   character(len=7), parameter :: mod_name = 'tracers'
!
!-----------------------------------------------------------------------
!
      integer  unit,io,ierr, logunit, outunit
      character(len=3) :: SOA_tracer
      integer :: ios
      logical :: exists

!
      data SOA_tracer/'SOA'/

!
      if (module_is_initialized) return

!----- set initial value of soa ------------

      nSOA = get_tracer_ndx(tracer_names,'SOA')
      if (nSOA > 0) then
         if (nSOA > 0 .and. me == master) &
                 write (*,30) SOA_tracer,nsoa
      endif

  30   format (A,' was initialized as tracer number ',i2)

!----- check for other required tracers ------------

      nOH = get_tracer_ndx(tracer_names,'OH')
      if (nOH > 0 .and. me == master) then
         print *, "OH was initialized as tracer number ", nOH
      endif
      nC4H10 = get_tracer_ndx(tracer_names,'C4H10')
      if (nC4H10 > 0 .and. me == master) then
         print *, "C4H10 was initialized as tracer number ", nC4H10
      endif
      nISOP = get_tracer_ndx(tracer_names,'ISOP')
      if (nISOP > 0 .and. me == master) then
         print *, "ISOP was initialized as tracer number ", nISOP
      end if
      nTERP = get_tracer_ndx(tracer_names,'C10H16')
      if (nTERP > 0 .and. me == master) then
         print *, "C10H16 was initialized as tracer number ", nTERP
      end if
      nO3 = get_tracer_ndx(tracer_names,'O3')
      if (nO3 > 0 .and. me == master) then
         print *, "O3 was initialized as tracer number ", nO3
      end if
      nNO3 = get_tracer_ndx(tracer_names,'NO3')
      if (nNO3 > 0 .and. me == master) then
         print *, "NO3 was initialized as tracer number ", nNO3
      end if
  40  format (A,' was initialized as tracer number ',i2)

         if (me == master) then
            write (*,*) 'gfdl_soa_mod: Using interactive tracers'
         end if
         if (nOH==NO_TRACER .or. nC4H10==NO_TRACER &
             .or. nOH <= 0 .or. nC4H10 <= 0 ) &
           call errmsg ('gfdl_soa_mod', &
            'OH and C4H10 tracers must be present if use_interactive_tracers=T', .true.)

      if (use_interactive_BVOC_emis) then
         if (me == master) then
            write (*,*) 'gfdl_soa_mod: Using interactive BVOC emissions'
         endif
      endif
     
      module_is_initialized = .TRUE.

!-----------------------------------------------------------------------
end subroutine gfdl_soa_init



!#####################################################################

subroutine gfdl_soa_endts             


!      call unset_interpolator_time_flag (gas_conc_interp)
!      call unset_interpolator_time_flag (isoprene_interp)
!      call unset_interpolator_time_flag (terpene_interp)


end subroutine gfdl_soa_endts



!#####################################################################

!</SUBROUTINE>

!#######################################################################
!<SUBROUTINE NAME="gfdl_soa_end">
!<OVERVIEW>
!  The destructor routine for the soa module.
!</OVERVIEW>
! <DESCRIPTION>
! This subroutine writes the version name to logfile and exits. 
! </DESCRIPTION>
!<TEMPLATE>
! call gfdl_soa_end
!</TEMPLATE>
 subroutine gfdl_soa_end

!      call interpolator_end (gas_conc_interp)
!      call interpolator_end (isoprene_interp)
!      call interpolator_end (terpene_interp)
      module_is_initialized = .FALSE.

 end subroutine gfdl_soa_end

!</SUBROUTINE>
!-----------------------------------------------------------------------
      SUBROUTINE gfdl_soa_chem(pwt,temp,pfull, phalf, dt,     &
                          SOA, OH, C4H10, xbvoc, &
                          use_interactive_BVOC_emis, &
                          SOA_dt, is,ie,js,je,kbot)

! ****************************************************************************
      real, intent(in),    dimension(:,:,:)          :: pwt
      real, intent(in),    dimension(:,:,:)          :: temp,pfull,phalf
      real, intent(in)                               :: dt
      real, intent(in),    dimension(:,:,:)          :: SOA, OH, C4H10
      real, intent(in),    dimension(:,:,:)          :: xbvoc
      logical, intent(in)                            :: use_interactive_BVOC_emis
      real, intent(out),   dimension(:,:,:)          :: SOA_dt 
      integer, intent(in),  dimension(:,:), optional :: kbot
      integer, intent(in)                            :: is,ie,js,je
! Working vectors
      real, dimension(size(SOA,1),size(SOA,2),size(SOA,3)) :: &
               SOA_chem, OH_conc, C4H10_conc
      real, dimension(size(SOA,1),size(SOA,2)) :: &
               SOA_prod, &
               xu, dayl, h, hl, hc, hred, fac_OH, fact_OH, &
	       isoprene_emis, terpene_emis

      real, parameter                            :: c4h10_SOA_yield = 0.1
      real, parameter                            :: A0 = 0.006918
      real, parameter                            :: A1 = 0.399912
      real, parameter                            :: A2 = 0.006758
      real, parameter                            :: A3 = 0.002697
      real, parameter                            :: B1 = 0.070257
      real, parameter                            :: B2 = 0.000907
      real, parameter                            :: B3 = 0.000148
      real, parameter                            :: A_C4H10_OH = 1.55E-11
      real, parameter                            :: B_C4H10_OH = 540.

      real, parameter                            :: A_ISOP_OH  = 3.1E-11
      real, parameter                            :: B_ISOP_OH  = 350.
      real, parameter                            :: A_TERP_OH  = 1.2E-11
      real, parameter                            :: B_TERP_OH  = 440.
      real, parameter                            :: A_ISOP_O3  = 1.0E-14
      real, parameter                            :: B_ISOP_O3  = -1970.
      real, parameter                            :: A_TERP_O3  = 5.3E-16
      real, parameter                            :: B_TERP_O3  = -530.
      real, parameter                            :: A_ISOP_NO3 = 3.3E-12
      real, parameter                            :: B_ISOP_NO3 = -450.
      real, parameter                            :: A_TERP_NO3 = 1.2E-12
      real, parameter                            :: B_TERP_NO3 = 490.
   
      real                                       :: decl, hd, x
      integer                                    :: i,j,k,id,jd,kd,logunit
      integer                                    :: istep, nstep
      real                                       :: air_dens

      
! Local grid sizes
      id=size(SOA,1); jd=size(SOA,2); kd=size(SOA,3)

      isoprene_emis(:,:)    = 0.0
      terpene_emis(:,:)    = 0.0
!----------------------------------------------------------------------
!    SOA_dt initially contains chemical production (pseudo-emission added later)
!----------------------------------------------------------------------
         do i=1,id
         do j=1,jd
         do k=1,kd
            air_dens = 10.*pfull(i,j,k)/(boltz*temp(i,j,k)) ! molec/cm3
            SOA_dt(i,j,k) = A_C4H10_OH * exp( -B_C4H10_OH/temp(i,j,k) ) * c4h10_SOA_yield &
                          * C4H10(i,j,k)*OH(i,j,k)*air_dens
         enddo
         enddo
         enddo

           SOA_dt(:,:,kd) = SOA_dt(:,:,kd) +  (isoprene_emis(:,:) + terpene_emis(:,:))/ pwt(:,:,kd) 

      SOA_chem(:,:,:)=SOA_dt(:,:,:)*pwt(:,:,:)

!----------------------------------------------------------------------
!    Pseudo-emission of SOA from biogenic VOCs,
!    ... add to SOA_dt (in lowest model layer)
!----------------------------------------------------------------------
      if (.not. use_interactive_BVOC_emis) then
         isoprene_emis(:,:)    = 0.0
         terpene_emis(:,:)    = 0.0
        !JianHe: placeholder, currently read offline biogenic
         isoprene_emis(:,:) = xbvoc(:,:,1) !isop, molec/cm2/s
         terpene_emis(:,:)  = xbvoc(:,:,2) !terp
      else
        !JianHe: placeholder, include xbvoc in the future
        ! isoprene_emis(:,:) = xbvoc(:,:,ind_xbvoc_ISOP)
        ! terpene_emis(:,:)  = xbvoc(:,:,ind_xbvoc_TERP)
               
      endif

!----------------------------------------------------------------------
!    Scale from BVOC emissions to SOA pseudo-emission
!----------------------------------------------------------------------
      isoprene_emis(:,:) = isoprene_emis(:,:) * isoprene_factor * isoprene_SOA_yield
      terpene_emis(:,:) = terpene_emis(:,:) * terpene_factor * terpene_SOA_yield

      SOA_dt(:,:,kd) = SOA_dt(:,:,kd) &
                     + (isoprene_emis(:,:) + terpene_emis(:,:)) / pwt(:,:,kd)

! column production of SOA 


      SOA_prod = 0.
      do k=1,kd
        SOA_prod = SOA_prod +  SOA_chem(:,:,k)
      end do

end subroutine gfdl_soa_chem


end module gfdl_soa_mod
