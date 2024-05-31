!>\file catchem_drydep_wrapper.F90
!! This file is GSDChem dry deposition wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Kate.Zhang@noaa.gov 04/2023
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov

 module catchem_drydep_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use catchem_constants, only : WTMAIR, WTMH2O, epsilc
   use drydep_wesely_mod
   use drydep_gocart_mod
   use dep_vertmx_mod
   use gocart_diag_mod

   use gfdl_time_utls_mod, only : time_type, &
                                  get_date, &
                                  set_date, &
                                  set_time, &
                                  days_in_year, &
                                  real_to_time_type, &
                                  time_type_to_real, &
                                  operator(+), operator(-), &
                                  time_interp
   use gfdl_drydep_mod

   implicit none

   private

   public :: catchem_drydep_wrapper_init, catchem_drydep_wrapper_run, catchem_drydep_wrapper_finalize

   !JianHe: read prescribed depvel from external file (depvel_in)
   integer, parameter :: ngas_drydep = 22
   character(len=32), dimension(ngas_drydep), save :: depvel_name = &      !
      (/ "co", "ch2o", "o3", "no", "no2", "hno3", "hno4", "n2o5", "ch4", "ch3ooh", &
         "h2o2", "pan", "pmn", "ch3coch3", "glyc", "hac", "rip", "so2", "nh3", "hobr", &
         "hbr", "brno3" /)
   
   type(time_type) :: Time_init,Time,Time_next,dt_time

   logical :: is_initialized = .false.

contains

!>\brief The subroutine initializes the GFDL AM4 drydep.
!!
!> \section arg_table_catchem_drydep_wrapper_init Argument Table
!! \htmlinclude catchem_drydep_wrapper_init.html
!!
      subroutine catchem_drydep_wrapper_init(me, master, tracer_names, &
                      gaschem_opt, do_am4chem, errmsg, errflg)

       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master
       integer, intent(in)  :: gaschem_opt
       logical, intent(in)  :: do_am4chem

       character(len=32), intent(in) :: tracer_names(:)
       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (is_initialized) return

       if (gaschem_opt == 1 .or. do_am4chem) then
          if (me == master) then
            print *, 'Start GFDL AM4 dry depostion for both gas and aerosols'
          endif

          call gfdl_drydep_init(me, master, tracer_names)

       end if

       is_initialized = .true.
      end subroutine catchem_drydep_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_drydep_wrapper_finalize Argument Table
!!
      subroutine catchem_drydep_wrapper_finalize()
      end subroutine catchem_drydep_wrapper_finalize

!> \defgroup catchem_group CATChem drydep wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_drydep_wrapper CATChem drydep wrapper Module  
!> \ingroup catchem_drydep_group
!! This is the CATChem drydep wrapper Module
!! \section arg_table_catchem_drydep_wrapper_run Argument Table
!! \htmlinclude catchem_drydep_wrapper_run.html
!!
!> @{
    subroutine catchem_drydep_wrapper_run(me, master, tracer_names, &
                   im, kte, kme, ktau, dt, land,      &
                   landfrac,oceanfrac, fice, ugrs, vgrs, &
                   ustar, rlat, rlon, tskin, julian, rainc_cpl, hf2d, pb2d,   &
                   pr3d, ph3d, phl3d, prl3d, tk3d, spechum, exch, depvel_in,   &
                   vegtype, sigmaf, jdate, idat, dswsfc, zorl, snow_cplchm,      & 
                   ntrac,ntchs,ntchm,ndchm,ntso2,ntsulf,ntDMS,ntmsa,ntpp25,         &
                   ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,     &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,            &
                   ntsoa, ntage, ntaoanh, ntextinct, &
                   ntchmdiag,gq0,qgrs,drydep,ddep,chem_conv_tr_in,                 &
                   chem_opt,gaschem_opt,do_am4chem,errmsg,errflg)

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    character(len=32), intent(in) :: tracer_names(:)

    integer,        intent(in) :: im,kte,kme,ktau,jdate(8),idat(8)
    integer,        intent(in) :: ntchs,ntchm,ndchm
    integer,        intent(in) :: ntrac,ntchmdiag,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    integer,        intent(in) :: ntage, ntaoanh,ntsoa,ntextinct
    real(kind_phys),intent(in) :: dt,julian

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land, vegtype     
    real(kind_phys), dimension(im), intent(in) :: ustar,                  &
                rlat,rlon, tskin, rainc_cpl,                              &
                hf2d, pb2d, sigmaf, dswsfc, zorl, snow_cplchm, &
                landfrac, oceanfrac, fice
    real(kind_phys), optional, intent(in) :: depvel_in(:,:)  !JianHe
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d, spechum, exch
    real(kind_phys), dimension(im,kte), intent(in) :: ugrs,vgrs
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: drydep
    real(kind_phys), optional, intent(inout) :: ddep(:,:)
    integer,         intent(in) :: chem_conv_tr_in, chem_opt, gaschem_opt
    logical,         intent(in) :: do_am4chem
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,        &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, zmid, exch_h, pwt, &
                     u_phy, v_phy

    real(kind_phys), dimension(ims:im, jms:jme) :: ust, tsk,              &
                     xland, xlat, xlong, rcav, hfx, pbl, frac_open_sea, &
                     rlat_in,rlon_in

!>- chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: dry_fall
    real(kind_phys), dimension(im, 1, ntchmdiag, 4) :: trdf

    integer :: ide, ime, ite, kde, julday

    real(kind_phys), dimension(ims:im, jms:jme) :: vegfrac, rmol, gsw, znt
    real(kind_phys), dimension(ims:im, jms:jme) :: snowh  
    integer,         dimension(ims:im, jms:jme) :: ivgtyp

    integer :: current_month

    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: ac3, ahno3, anh3, asulf, cor3, h2oai, h2oaj, nu3
    real(kind_phys), dimension(ims:im, jms:jme) :: dep_vel_o3, e_co

    real(kind_phys) :: gmt
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

!>-- local variables
    integer :: i, j, jp, k, kp, n, nv, nt
    integer :: chem_conv_tr
    ! dry deposition velocity
    real(kind_phys), dimension(ims:im, jts:jte, num_chem ) ::   ddvel, dsinku
    real(kind_phys), dimension(ims:im, jms:jme, 1:ngas_drydep) :: depvel_am4

    real(kind_phys) :: delz_at_w, cdt, factor

    ! turbulent transport
    real(kind_phys), dimension(kts:kte ) :: pblst,zz
    real(kind_phys), dimension(kts:kte+1):: ekmfull,zzfull
    real(kind_phys), dimension(kts:kte ) :: dryrho_1d


    errmsg = ''
    errflg = 0

    chem_conv_tr      = chem_conv_tr_in

    h2oai = 0.
    h2oaj = 0.
    nu3   = 0.
    ac3   = 0.
    cor3  = 0.
    asulf = 0.
    ahno3 = 0.
    anh3  = 0.
    e_co  = 0.
    dep_vel_o3 = 0.
    ddvel(:,:,:) = 0.
    dsinku(:,:,:) = 0.

    gmt = real(idat(5))
    julday = real(julian)                                       

    current_month=jdate(2)

    !JianHe: Model current time
    Time = set_date(jdate(1),jdate(2),jdate(3), &
                        jdate(5),jdate(6),jdate(7))

    dry_fall = 0._kind_phys

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
    enddo

!!!

!>- get ready for chemistry run
    call catchem_prep_drydep(                                          &
        gaschem_opt, do_am4chem, &
        ustar,land,rlat,rlon,tskin,                                     &
        landfrac,oceanfrac, fice,  &
        ugrs, vgrs, &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,exch,      &
        vegtype,sigmaf,dswsfc,zorl,snow_cplchm,hf2d,pb2d,                  &
        ust,tsk,xland,frac_open_sea,xlat,xlong,                           &
        rlat_in,rlon_in, &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,pwt,                               &
        u_phy, v_phy, &
        t8w,exch_h,z_at_w,zmid,                                         &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,          &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntchs,ntchm,ndchm,ntrac,gq0,num_chem, num_moist,                &
        ppm2ugkg,moist,chem, &
        ivgtyp,vegfrac,rmol,gsw,znt,hfx,pbl,snowh,                      &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    !>-- compute dry deposition
 
    !if (gaschem_opt == 1 .or. do_am4chem) then
#ifdef AM4_CHEM
     ! JianHe: we need initialize in the wrapper file
     ! instead of in lower subroutines to avoid issue
     ! for diagnostic output (required by debugger turned on)
     ddep = 0._kind_phys

     ! read external drydep velocity  
     do j=jts,jte
        do i=its,ite
          depvel_am4(i,j,1:ngas_drydep)=depvel_in(i,1:ngas_drydep)
        enddo
     enddo

      ! AM4 drydep for prog chem tracer only
       do n = 1,ntchm  ! for gas+aerosol, tracer index in chem array
          nt = ntchs+n-1   ! tracer index in tracer array
          call gfdl_dry_deposition(me, master, tracer_names, &
                     nt, its, jts, u_phy(:,kts,:), v_phy(:,kts,:), &
                     t_phy(:,kts,:), pwt(:,kts,:), p_phy(:,kts,:), &
                     dz8w(:,kts,:), ust, &
                     xland, frac_open_sea, dsinku(:,:,n), dt, &
                     chem(:,kts,:,n), rlat_in, rlon_in, Time, &
                     depvel_am4(:,:,:),ddvel(:,:,n))
       end do

       !
       !   Compute dry deposition according to NGAC
       !
       cdt = real(dt, kind=kind_phys)
       do j = jts,jte
         do i = its,ite
           do nv = 1, ntchm
             factor = 1._kind_phys - exp(-ddvel(i,j,nv)*cdt/dz8w(i,kts,j))
             !JianHe: ug/m2/s for aerosol, ppm/m2/s for gas
             dry_fall(i,j,nv) = max(0.0, factor * chem(i,kts,j,nv)) & 
                         * (p8w(i,kts,j)-p8w(i,kts+1,j))/g/dt

             
          end do
         end do
       end do

      ! JianHe:output dry dep for o3, nh3, and so4
      !note so4 is in ppm in the tracer array
      ddep(:,1) = dry_fall(:,1,p_o3)*1.e-6*48./WTMAIR ! kg/m2/s
      ddep(:,2) = dry_fall(:,1,p_nh3)*1.e-6*17./WTMAIR ! kg/m2/s
      ddep(:,3) = dry_fall(:,1,p_so4)*1.e-6*96./WTMAIR ! kg/m2/s

#else
      do j = jts,jte
        do i = its,ite
          delz_at_w = z_at_w(i,kts+1,j) - z_at_w(i,kts,j)

        !JianHe: placeholder, in the future, we will do
        !gas and aerosol drydep seperately
        !We will not based on chem_opt, but on gas/aero schemes

        IF( chem_opt /= GOCART_SIMPLE ) THEN
          ! wesely for gases 
       !   call wesely_driver(current_month,julday, &
       !       t_phy(i,kts,j),moist(i,kts,j,:),p8w(i,kts,j),     &
       !       rcav(i,j),p_phy(i,kts,j),ddvel(i,j,:),     &
       !       ivgtyp(i,j), chem_opt,&
       !       tsk(i,j),gsw(i,j),vegfrac(i,j),rmol(i,j),        &
       !       ust(i,j),znt(i,j),delz_at_w,snowh(i,j))
        ENDIF

        IF (( chem_opt == GOCART_SIMPLE ) .or.            &
              ( chem_opt == GOCARTRACM_KPP)  .or.            &
              ( chem_opt == 316)  .or.            &
              ( chem_opt == 317)  .or.            &
!             ( chem_opt == 502)  .or.            &
              (chem_opt == 304          )) then
          ! this does aerosol species (dust,seas, bc,oc,sulf) for gocart only
          call gocart_drydep_driver(  &
              p8w(i,kts,j),rho_phy(i,kts,j),dz8w(i,kts,j), &
              ddvel(i,j,:),xland(i,j),hfx(i,j),ivgtyp(i,j), &
              tsk(i,j),pbl(i,j),ust(i,j),znt(i,j))
        ELSE if (chem_opt == 501 ) then
! for caesium .1cm/s
!
          ddvel(i,j,:)=.001

        ELSE if (chem_opt == 108 ) then
!!       call soa_vbs_depdriver (ust,t_phy,                    &
!!               moist,p8w,rmol,znt,pbl,           &
!!               alt,p_phy,chem,rho_phy,dz8w,                    &
!!               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,      &
!!               aer_res,ddvel(:,:,numgas+1:num_chem),           &
!!               num_chem-numgas,                                &
!!               ids,ide, jds,jde, kds,kde,                      &
!!               ims,ime, jms,jme, kms,kme,                      &
!!               its,ite, jts,jte, kts,kte                       )
! limit aerosol ddvels to <= 0.5 m/s
! drydep routines occasionally produce unrealistically-large particle
! diameter leading to unrealistically-large sedimentation velocity 
          ddvel(i,j,numgas+1:num_chem) = min( 0.50, ddvel(i,j,numgas+1:num_chem))
        ELSE
          !Set dry deposition velocity to zero when using the
          !chemistry tracer mode.
          ddvel(i,j,:) = 0.
        END IF

        !
        !   Compute dry deposition according to NGAC
        !
        cdt = real(dt, kind=kind_phys)
        do nv = 1, num_chem
          factor = 1._kind_phys - exp(-ddvel(i,j,nv)*cdt/dz8w(i,kts,j))
          dry_fall(i,j,nv) = max(0.0, factor * chem(i,kts,j,nv)) & !ug/m2/s
                         * (p8w(i,kts,j)-p8w(i,kts+1,j))/g/dt
        end do
      end do
    end do

#endif 
! gaschem

        !
        !   This will be called later from subgrd_transport_driver.F !!!!!!!!
        !
    do j=jts,jte
      do i=its,ite
        !if(p_dust_1.gt.1)dep_vel_o3(i,j)=ddvel(i,j,p_dust_1)
        if (p_o3 > 0) dep_vel_o3(i,j)=ddvel(i,j,p_o3)
        pblst=0.

        !
        !-- start with vertical mixing
        !
        do k=kts,kte+1
           zzfull(k)=z_at_w(i,k,j)-z_at_w(i,kts,j)
        enddo

        if (chem_conv_tr == CTRA_OPT_NONE) then
          ekmfull = 0.
        else
          ekmfull(kts)=0.
          do k=kts+1,kte
           ekmfull(k)=max(1.e-6,exch_h(i,k,j))
          enddo
          ekmfull(kte+1)=0.
        end if
 
        do k=kts,kte
           zz(k)=zmid(i,k,j)-z_at_w(i,kts,j)
        enddo
!   vertical mixing routine (including deposition)
!   need to be careful here with that dumm tracer in spot 1
!   do not need lho,lho2
!   (03-may-2006 rce - calc dryrho_1d and pass it to vertmx)
!
!     if(p_o3.gt.1)dep_vel_o3(i,j)=ddvel(i,j,p_o3)
!        do nv=1,num_chem-0
      !JianHe: should only for prog chem tracer?
         do nv=1,ntchm
           do k=kts,kte
              pblst(k)=max(epsilc,chem(i,k,j,nv))
              dryrho_1d(k) = 1./rri(i,k,j)
           enddo

           mix_select: SELECT CASE(chem_opt)
           CASE (RADM2SORG_AQ, RACMSORG_AQ, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ)
!           if(.not.is_aerosol(nv))then ! mix gases not aerosol
               call vertmx(dt,pblst,ekmfull,dryrho_1d, &
                           zzfull,zz,ddvel(i,j,nv),kts,kte)

!           endif

           CASE DEFAULT
               call vertmx(dt,pblst,ekmfull,dryrho_1d, &
                           zzfull,zz,ddvel(i,j,nv),kts,kte)

           END SELECT mix_select

           do k=kts,kte
              chem(i,k,j,nv)=max(epsilc,pblst(k))
           enddo
        enddo 
      enddo
    enddo

    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       !JianHe: update here
       gq0(i,k,ntchs:ntrac) = ppm2ugkg(1:num_chem)* max(epsilc,chem(i,k,1,1:num_chem))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntchs:ntrac) =gq0(i,k,ntchs:ntrac)
     enddo
    enddo

    ! -- output dry deposition
    call gocart_diag_store(2, dry_fall, trdf)

    drydep (:,:)=trdf(:,1,:,2)

!
   end subroutine catchem_drydep_wrapper_run
!> @}
  subroutine catchem_prep_drydep(                                     &
        gaschem_opt, do_am4chem, &
        ustar,land,rlat,rlon,ts2d,                                     &
        landfrac, oceanfrac, fice, &
        ugrs, vgrs, &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,exch,            &
        vegtype,sigmaf,dswsfc,zorl,snow_cplchm,hf2d,pb2d,                 &
        ust,tsk,xland,frac_open_sea,xlat,xlong,                        &
        rlat_in, rlon_in, &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,pwt,                              &
        u_phy, v_phy, &
        t8w,exch_h,z_at_w,zmid,                                        &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,ntss1,ntss2,ntss3,ntss4,ntss5,         &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntchs,ntchm,ndchm,ntrac,gq0,num_chem, num_moist,               &
        ppm2ugkg,moist,chem, &
        ivgtyp,vegfrac,rmol,gsw,znt,hfx,pbl,snowh,                     &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !FV3 input variables
    integer, intent(in)  :: gaschem_opt, ntchs
    integer, dimension(ims:ime), intent(in) :: land, vegtype
    integer, intent(in) :: ntchm,ndchm
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer, intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         ustar, rlat, rlon, ts2d, sigmaf, dswsfc, zorl, snow_cplchm, hf2d, pb2d
    real(kind=kind_phys), dimension(ims:ime), intent(in) :: landfrac, oceanfrac,fice
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: phl3d,tk3d,prl3d,spechum,exch
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: ugrs,vgrs
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0
    logical, intent(in) :: do_am4chem

    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(1:num_chem), intent(in) :: ppm2ugkg

    
    integer,dimension(ims:ime, jms:jme), intent(out) :: ivgtyp
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, p_phy, rho_phy, dz8w, p8w, t8w, zmid,         &
         exch_h, pwt, u_phy, v_phy
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         ust, tsk, xland, xlat, xlong, vegfrac, rmol, gsw, znt, hfx,    &
         pbl, snowh, frac_open_sea, rlat_in, rlon_in
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w

    ! -- local variables
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: delp
    integer i,ip,j,jp,k,kp,kk,kkp,nv,l,ll,n,nt


    ! -- initialize output arrays
    ivgtyp         = 0._kind_phys
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    u_phy          = 0._kind_phys 
    v_phy          = 0._kind_phys 
    t8w            = 0._kind_phys
    zmid           = 0._kind_phys
    exch_h         = 0._kind_phys
    ust            = 0._kind_phys
    tsk            = 0._kind_phys
    xland          = 0._kind_phys
    frac_open_sea  = 0._kind_phys ! ocean fraction exclu. ice
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    vegfrac        = 0._kind_phys
    rmol           = 0._kind_phys
    gsw            = 0._kind_phys
    znt            = 0._kind_phys
    hfx            = 0._kind_phys
    pbl            = 0._kind_phys
    snowh          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    rlat_in        = 0._kind_phys
    rlon_in        = 0._kind_phys

    pwt            = 0._kind_phys ! Pressure weighting (air mass) for each layer (kg/m2)
    delp           = 0._kind_phys

    do i=its,ite
     tsk  (i,1)=ts2d (i)
     ust  (i,1)=ustar(i)
     !xland(i,1)=real(land(i))
     xland(i,1)=landfrac(i)
     frac_open_sea(i,1) = min(max(0., oceanfrac(i)*(1.-fice(i))),oceanfrac(i))
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     rlat_in(i,1)=rlat(i)
     rlon_in(i,1)=rlon(i)
     gsw  (i,1)=dswsfc(i)
     znt  (i,1)=zorl(i)*0.01
     hfx  (i,1)=hf2d(i)
     pbl  (i,1)=pb2d(i)
     snowh(i,1)=snow_cplchm(i)*0.001
     ivgtyp (i,1)=vegtype(i)
     vegfrac(i,1)=sigmaf (i)
    enddo
   
    rmol=0.

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
          z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
          u_phy(i,k,j) = ugrs(ip,kp)
          v_phy(i,k,j) = vgrs(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          p8w(i,k,j)=pr3d(ip,kp)
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte+1
        kk=min(k,kte)
        kkp = kk - kts + 1
        do i=its,ite
          ip = i - its + 1
          dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
          t_phy(i,k,j)=tk3d(ip,kkp)
          p_phy(i,k,j)=prl3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          moist(i,k,j,:)=0.
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          if (t_phy(i,k,j) > 265.) then
            moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
            moist(i,k,j,3)=0.
            if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
          else
            moist(i,k,j,2)=0.
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
            if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
          endif
          !--
          zmid(i,k,j)=phl3d(ip,kkp)/g
        enddo
      enddo
    enddo

    ! -- the imported atmospheric heat diffusivity is only available up to kte-1
    do j=jts,jte
      jp = j - jts + 1
      do k=kts,kte-1
        kkp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
          exch_h(i,k,j)=exch(ip,kkp)
        enddo
      enddo
    enddo

    do j=jts,jte
      do k=2,kte
        do i=its,ite
          t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
        enddo
      enddo
    enddo

    ! -- only used in phtolysis....
    do j=jts,jte
      do i=its,ite
        t8w(i,1,j)=t_phy(i,1,j)
        t8w(i,kte+1,j)=t_phy(i,kte,j)
      enddo
    enddo

    do k=kts,kte
       delp(:,k,:)=p8w(:,k,:)-p8w(:,k+1,:)
       pwt(:,k,:)=delp(:,k,:)/g
    enddo

    do k=kms,kte
     do i=ims,ime
       !JianHe: 08/2023: can we use some kind of index loop for this?
         do nt = 1, ntchm+ndchm
           ! start with so2, we need to make sure the order of tracers
           ! in chem array consistent with field table
           ! chem array stores all chem tracers (prog+diag)
           ! num_chem shoud be ntchm + ndchm
           chem(i,k,jts,nt) = max(epsilc,gq0(i,k,ntchs+nt-1)/ppm2ugkg(nt))
         end do
     enddo
    enddo

  end subroutine catchem_prep_drydep

!> @}
  end module catchem_drydep_wrapper
