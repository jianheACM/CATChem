!>\file catchem_lnox_wrapper.F90
!! This file is CATChem lightning nox emission wrapper with CCPP coupling to FV3
!! Jian.He@noaa.gov 10/2023
!! Revision History:
!! 10/2023, Create wrapper for LNOx

 module catchem_lnox_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use catchem_constants, only: WTMAIR, AVOGNO,epsilc
   use gfdl_lnox_mod, only : gfdl_lnox_init, gfdl_lnox_hook

   implicit none

   private

   public :: catchem_lnox_wrapper_init, catchem_lnox_wrapper_run, catchem_lnox_wrapper_finalize

   logical :: is_initialized = .false.

contains

!> \brief Brief description of the subroutine
!!
!> \section arg_table_catchem_lnox_wrapper_init Argument Table
!! \htmlinclude catchem_lnox_wrapper_init.html
!!
      subroutine catchem_lnox_wrapper_init(me, master, lnox_opt, lght_no_prd_factor, &
                                           min_land_frac_lght, &
                                           normalize_lght_no_prd_area,errmsg, errflg)

       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master
       integer, intent (in) :: lnox_opt
       real(kind_phys), intent (in) :: lght_no_prd_factor,min_land_frac_lght
       logical, intent (in) :: normalize_lght_no_prd_area

       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

       integer :: ios
       logical :: exists

       
       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (is_initialized) return

       if (lnox_opt > 0) then
          if (me == master) then
            print *, 'Initialize lightning NOx: ', lght_no_prd_factor, &
              normalize_lght_no_prd_area, min_land_frac_lght
          endif

       call gfdl_lnox_init(me, master, lght_no_prd_factor, normalize_lght_no_prd_area, min_land_frac_lght )

       is_initialized = .true.

       end if

      end subroutine catchem_lnox_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_lnox_wrapper_finalize Argument Table
!!
      subroutine catchem_lnox_wrapper_finalize()
      end subroutine catchem_lnox_wrapper_finalize

!> \defgroup catchem_group CATChem lnox wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_lnox_wrapper CATChem lnox wrapper Module  
!> \ingroup catchem_lnox_group
!! This is the CATChem lnox wrapper Module
!! \section arg_table_catchem_lnox_wrapper_run Argument Table
!! \htmlinclude catchem_lnox_wrapper_run.html
!!
!>\section catchem_lnox_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_lnox_wrapper_run(me, master, &
                   im, kte, kme, dt, garea, rlat, rlon,   &
                   landfrac, ktop, kbot,                  &
                   pr3d, ph3d,phl3d, prl3d, tk3d,  spechum,     &
                   bioem, ntrac, ntno, gq0, qgrs, lnox_opt,        &
                   errmsg,errflg)

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master
    integer,        intent(in) :: im,kte,kme
    integer,        intent(in) :: ntrac,ntno
    integer, dimension(im), intent(in) :: kbot,ktop
    real(kind_phys),intent(in) :: dt

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im), intent(in) :: landfrac
    real(kind_phys), dimension(im), intent(in) :: garea, rlat,rlon
    real(kind_phys), dimension(im,kme),intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte),intent(in) :: phl3d, prl3d, tk3d,  spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), optional, intent(inout) :: bioem(:,:)     

    integer,           intent(in) :: lnox_opt
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,        &
                     dz8w, p8w, rho_phy, z_phy,p_phy, z_at_w, pwt

    !JianHe: AM4 lev is in 3rd dimension 
    real(kind_phys), dimension(1:im,jms:jme,1:kte) :: t_am4, z_am4
    real(kind_phys), dimension(1:im,jms:jme,1:kme) :: z8w_am4, p8w_am4

    real(kind_phys), dimension(ims:im, jms:jme) :: xland, dxy, xlat, xlong
    integer, dimension(ims:im, jms:jme) :: cldtop ! cloud top level index
    integer, dimension(ims:im, jms:jme) :: cldbot ! cloud bottom level index
    real(kind_phys), dimension(ims:im, jms:jme) :: ilat
    real(kind_phys), dimension(ims:im, jms:jme, kte) :: zm     ! geopot height above surface at midpoints (m)
    real(kind_phys), dimension(ims:im, jms:jme, kme) :: zint   ! geopot height above surface at interfaces (m)

    real(kind_phys), dimension(ims:im, jms:jme, 1:kte) :: prod_no ! production of NOx (molec cm^-3 s^-1)
    real(kind_phys), dimension(ims:im, jms:jme) :: prod_no_col ! production of NOx (molec cm^-2 s^-1)
    
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    integer :: ide, ime, ite, kde
    integer, dimension(ims:im, jms:jme)  :: iktop,ikbot

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
if (lnox_opt > 0) then
    call catchem_lnox_prep(                                         &
        landfrac,garea,rlat,rlon,                           &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,                   &
        ktop,kbot,xland,dxy,xlat,xlong,                   &
        rri,t_phy,rho_phy,z_phy,p_phy,z_at_w,dz8w,p8w,pwt,                 &
        ntno,ntrac,gq0,num_chem,ppm2ugkg,                              &
        chem,cldtop,cldbot,                                   &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)



    !if (lnox_opt == LNOx_OPT_AM4) then
#ifdef AM4_CHEM
      do j=jts,jte
        do i=its,ite
         !JianHe: AM4 needs vertical level from model top to bottom
          do k = kts, kte
            kp = kte-k+1
            t_am4(i,j,kp) = t_phy(i,k,j)  !K
            z_am4(i,j,kp) = z_phy(i,k,j)
          enddo
          do k = kts, kte+1
            kp = kte+1-k+1
            z8w_am4(i,j,kp) = z_at_w(i,k,j)
          end do

          iktop(i,j) = kte - cldtop(i,j) + 1
          ikbot(i,j) = kte - cldbot(i,j) + 1
 
          if (me==master) then
            write(*,*) 'catchem_lnox: ', cldtop(i,j),cldbot(i,j),iktop(i,j),ikbot(i,j)
          endif
          ilat(i,j) = xlat(i,j)*pi/180.  ! need radiance input
          !
        enddo
      enddo  
   
      ! -- GFDL AM4 lnox scheme 
            
      call gfdl_lnox_hook(me, master, &
                    iktop, ikbot, xland, &
                    z_am4, z8w_am4, t_am4, &
                    prod_no,prod_no_col,dxy,ilat, &
                    its,jts )
   
      bioem(:,4)=prod_no_col(:,1)/AVOGNO*1.e10*3600.  !molec cm^-2 s^-1 to mol/km2/hr

      do j=jts,jte
        do i=its,ite
          do k = kts, kte
            kp = kte-k+1
            factor2=4.828e-4*dt*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_no)=chem(i,k,j,p_no)+ &
                             prod_no(i,j,kp)/AVOGNO*3600.*1.e15*dz8w(i,k,j)*1.e-3* & ! molec cm^-3 s^-1 to mol/km2/hr
                             factor2
          enddo
        enddo
      enddo

    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntno  )=ppm2ugkg(p_no   ) * max(epsilc,chem(i,k,1,p_no)) 
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntno )=gq0(i,k,ntno  )
     enddo
    enddo
#endif

end if

!
   end subroutine catchem_lnox_wrapper_run
!> @}
  subroutine catchem_lnox_prep(                                    &
        landfrac,garea,rlat,rlon,                                     &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,                  &
        ktop,kbot,xland,dxy,xlat,xlong,              &
        rri,t_phy,rho_phy,  &
        z_phy,p_phy,z_at_w,dz8w,p8w,pwt,                   &
        ntno,ntrac,gq0,num_chem,ppm2ugkg,                             &
        chem,cldtop,cldbot,                                      &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)


    !FV3 input variables
    real(kind=kind_phys), dimension(ims:ime), intent(in) :: landfrac
    integer, dimension(ims:ime), intent(in) :: ktop,kbot
    integer, intent(in) :: ntrac
    integer,        intent(in) :: ntno
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::                & 
         garea, rlat, rlon
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::     &
         pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,spechum
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


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
    integer, dimension(ims:ime, jms:jme),          intent(out) ::   cldtop, cldbot
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem

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
    cldtop         = 0
    cldbot         = 0

    pwt  = 0._kind_phys  ! Pressure weighting (air mass) for each layer (kg/m2)
    delp           = 0._kind_phys

    do i=its,ite
     dxy  (i,1)=garea(i)
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     xland(i,1)=landfrac(i)
     cldtop(i,1)=ktop(i)
     cldbot(i,1)=kbot(i)
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
       chem(i,k,jts,p_no   )=max(epsilc,gq0(i,k,ntno  )/ppm2ugkg(p_no))
     enddo
    enddo
#endif 

  end subroutine catchem_lnox_prep
!> @}
  end module catchem_lnox_wrapper
