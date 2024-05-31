!>\file catchem_soa_wrapper.F90
!! This file is CATChem soa wrapper with CCPP coupling to FV3
!! Jian.He@noaa.gov 11/2023
!! Revision History:
!! 11/2023, Create wrapper for SOA

 module catchem_soa_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use catchem_constants, only: AVOGNO, epsilc
   use gfdl_soa_mod, only : gfdl_soa_init, gfdl_soa_chem

   implicit none

   private

   public :: catchem_soa_wrapper_init, catchem_soa_wrapper_run, catchem_soa_wrapper_finalize

   logical :: is_initialized = .false.
   logical :: use_interactive_BVOC_emis 

contains

!> \brief Brief description of the subroutine
!!
!> \section arg_table_catchem_soa_wrapper_init Argument Table
!! \htmlinclude catchem_soa_wrapper_init.html
!!
      subroutine catchem_soa_wrapper_init(me, master, soa_opt, tracer_names, &
                                          isoprene_SOA_yield, &
                                          terpene_SOA_yield, &
                                          use_interactive_BVOC_emis,errmsg, errflg)

       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master
       integer, intent (in) :: soa_opt
       character(len=32), intent(in) :: tracer_names(:)
       real(kind_phys), intent (in) :: isoprene_SOA_yield,terpene_SOA_yield
       logical, intent (in) :: use_interactive_BVOC_emis

       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

       integer :: ios
       logical :: exists

       
       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (is_initialized) return

       if (soa_opt >= 1) then
          if (me == master) then
            print *, 'Initialize SOA: ', isoprene_SOA_yield, &
               terpene_SOA_yield, use_interactive_BVOC_emis
          endif
       end if

       if (soa_opt==1) then
         call gfdl_soa_init(me, master, tracer_names, &
              isoprene_SOA_yield, terpene_SOA_yield, use_interactive_BVOC_emis )
       endif

       is_initialized = .true.

      end subroutine catchem_soa_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_soa_wrapper_finalize Argument Table
!!
      subroutine catchem_soa_wrapper_finalize()
      end subroutine catchem_soa_wrapper_finalize

!> \defgroup catchem_group CATChem soa wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_soa_wrapper CATChem soa wrapper Module  
!> \ingroup catchem_soa_group
!! This is the CATChem soa wrapper Module
!! \section arg_table_catchem_soa_wrapper_run Argument Table
!! \htmlinclude catchem_soa_wrapper_run.html
!!
!>\section catchem_soa_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_soa_wrapper_run(me, master,ntextinct, &
                   im, kte, kme, dt, garea, rlat, rlon,   &
                   landfrac,pr3d, ph3d,phl3d, prl3d, tk3d,  spechum,     &
                   bioem, ntrac, ntsoa, ntoh, ntc4h10, &
                   ntisop,ntc10h16, gq0, qgrs, soa_opt,        &
                   errmsg,errflg)

    implicit none

    integer,        intent(in) :: me, master,ntextinct

    integer,        intent(in) :: im,kte,kme
    integer,        intent(in) :: ntrac,ntsoa,ntoh,ntc4h10,ntisop,ntc10h16
    real(kind_phys),intent(in) :: dt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im), intent(in) :: landfrac
    real(kind_phys), dimension(im), intent(in) :: garea, rlat,rlon
    real(kind_phys), dimension(im,kme),intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte),intent(in) :: phl3d, prl3d, tk3d,  spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,4), intent(in) :: bioem     

    integer,           intent(in) :: soa_opt
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,        &
                     dz8w, p8w, rho_phy, z_phy,p_phy, z_at_w, pwt

    !JianHe: AM4 lev is in 3rd dimension 
    real(kind_phys), dimension(1:im,jms:jme,1:kte) :: t_am4, z_am4, p_am4
    real(kind_phys), dimension(1:im,jms:jme,1:kme) :: z8w_am4, p8w_am4

    real(kind_phys), dimension(ims:im, jms:jme) :: xland, dxy, xlat, xlong
    real(kind_phys), dimension(ims:im, jms:jme) :: ilat
    real(kind_phys), dimension(ims:im,jms:jme, 2) :: xbvoc
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    integer :: ide, ime, ite, kde
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte) :: pwt_am4
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte) :: soa_am4, oh_am4, c4h10_am4, soa_dt

    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg


!>-- local variables
    integer :: i, j, jp, k, kp, n
    real(kind_phys) :: factor2

    errmsg = ''
    errflg = 0

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

!>- get ready for chemistry run
if (soa_opt > 0) then

    call catchem_soa_prep(                                         &
        landfrac,garea,rlat,rlon,                           &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,                   &
        bioem,xbvoc,xland,dxy,xlat,xlong,                   &
        rri,t_phy,rho_phy,z_phy,p_phy,z_at_w,dz8w,p8w,pwt,                    &
        ntsoa,ntoh,ntc4h10,ntisop,ntc10h16, &
        ntrac,gq0,num_chem,ppm2ugkg,chem,                    &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


#ifdef AM4_CHEM
    if (soa_opt == SOA_OPT_AM4) then
      do j=jts,jte
        do i=its,ite
         !JianHe: AM4 needs vertical level from model top to bottom
          do k = kts, kte
            kp = kte-k+1
            pwt_am4(i,j,kp) = pwt(i,k,j)
            t_am4(i,j,kp) = t_phy(i,k,j)  !K
            p_am4(i,j,kp) = p_phy(i,k,j)  !Pa
            z_am4(i,j,kp) = z_phy(i,k,j)
            soa_am4(i,j,kp) = chem(i,k,j,p_soa)*1.e-9 ! ug/kg to kg/kg
            oh_am4(i,j,kp) = chem(i,k,j,p_oh)*1.e-6  ! ppm to mol/mol
            c4h10_am4(i,j,kp) = chem(i,k,j,p_c4h10)*1.e-6  ! ppm to mol/mol
          enddo
          do k = kts, kte+1
            kp = kte+1-k+1
            z8w_am4(i,j,kp) = z_at_w(i,k,j)
            p8w_am4(i,j,kp) = p8w(i,k,j)
          end do

          ilat(i,j) = xlat(i,j)*pi/180.  ! need radiance input
          !
        enddo
      enddo  
           
      ! -- GFDL AM4 soa scheme 
            
      call gfdl_soa_chem(pwt_am4, t_am4, p_am4, &
                    p8w_am4, dt, &
                    soa_am4,oh_am4,c4h10_am4, &
                    xbvoc,use_interactive_BVOC_emis, &
                    soa_dt, its,ite,jts,jte)           

      soa_am4(its:ite,jts:jte,kts:kte) = soa_am4(its:ite,jts:jte,kts:kte) + &
                        soa_dt(its:ite,jts:jte,kts:kte)*dt
      do j=jts,jte
        do i=its,ite
          do k = kts, kte
            kp = kte-k+1
            gq0(i,k,ntsoa) = soa_am4(i,j,kp)*1.e9  ! reverse, kg/kg to ug/kg
          enddo
        enddo
      enddo
    endif

    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntsoa) = max(epsilc,gq0(i,k,ntsoa))
       qgrs(i,k,ntsoa )=gq0(i,k,ntsoa  )
     enddo
    enddo
#endif

endif
!
   end subroutine catchem_soa_wrapper_run
!> @}
  subroutine catchem_soa_prep(                                    &
        landfrac,garea,rlat,rlon,                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,                  &
        bioem,xbvoc,xland,dxy,xlat,xlong,              &
        rri,t_phy,rho_phy,  &
        z_phy,p_phy,z_at_w,dz8w,p8w,pwt,                    &
        ntsoa,ntoh,ntc4h10,ntisop,ntc10h16, &
        ntrac,gq0,num_chem,ppm2ugkg,chem,                       &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)


    !FV3 input variables
    real(kind=kind_phys), dimension(ims:ime), intent(in) :: landfrac
    integer, intent(in) :: ntrac
    integer,        intent(in) :: ntsoa,ntoh,ntc4h10,ntisop,ntc10h16
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         garea, rlat, rlon
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,spechum
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0
    real(kind_phys), dimension(ims:ime,4), intent(in) :: bioem

    !GSD Chem variables
    integer,intent(in) ::  num_chem
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(1:num_chem), intent(in) :: ppm2ugkg

    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, rho_phy, z_phy, p_phy, z_at_w, dz8w,p8w, pwt
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         xland, dxy, xlat, xlong
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem
    real(kind_phys), dimension(ims:ime, jms:jme, 2), intent(out) :: xbvoc

    ! -- local variables
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: delp
    integer i,ip,j,jp,k,kp,kk,kkp,l,ll,n

    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    z_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    z_at_w         = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    xland          = 0._kind_phys
    dxy            = 0._kind_phys
    chem           = 0._kind_phys
    xbvoc          = 0._kind_phys

    pwt  = 0._kind_phys  ! Pressure weighting (air mass) for each layer (kg/m2)
    delp           = 0._kind_phys

    do i=its,ite
     dxy  (i,1)=garea(i)
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     xland(i,1)=landfrac(i)
    enddo
   
    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
         ip = i - its + 1
         z_at_w(i,kts,j)=max(0.,ph3d(ip,1)/g)  !model interface
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
          z_phy(i,k,j)=max(0.,phl3d(ip,kp)/g)  !model midpoint
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
          !--
        enddo
      enddo
    enddo

    do k=kts,kte
       delp(:,k,:)=p8w(:,k,:)-p8w(:,k+1,:)
       pwt(:,k,:)=delp(:,k,:)/g
    enddo

#ifdef AM4_CHEM
    do k=kms,kte
     do i=ims,ime
       chem(i,k,jts,p_soa   )=max(epsilc,gq0(i,k,ntsoa  )/ppm2ugkg(p_soa))
       chem(i,k,jts,p_oh   )=max(epsilc,gq0(i,k,ntoh  )/ppm2ugkg(p_oh))
       chem(i,k,jts,p_c4h10   )=max(epsilc,gq0(i,k,ntc4h10  )/ppm2ugkg(p_c4h10))
       chem(i,k,jts,p_isop   )=max(epsilc,gq0(i,k,ntisop  )/ppm2ugkg(p_isop))
       chem(i,k,jts,p_c10h16   )=max(epsilc,gq0(i,k,ntc10h16)/ppm2ugkg(p_c10h16))
     enddo
    enddo

     do i=its,ite
       xbvoc(i,1,1) = bioem(i,1)*AVOGNO*1.e-10/3600. ! isop, mol/km2/hr to molec/cm2/s
       xbvoc(i,1,2) = bioem(i,2)*AVOGNO*1.e-10/3600. ! c10h16
     enddo
#endif

  end subroutine catchem_soa_prep
!> @}
  end module catchem_soa_wrapper

