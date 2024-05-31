!>\file catchem_anthropogenic_wrapper.F90
!! This file is GSDChem anthropogenic emission wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 07/2020
!! Kate.Zhang@noaa.gov 04/2023
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov
!! 05/2024, Include emissions for AM4 tracers

 module catchem_anthropogenic_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use catchem_constants, only: epsilc
   use mo_vinterp,      only : vinterp

   implicit none

   private

   public :: catchem_anthropogenic_wrapper_init, catchem_anthropogenic_wrapper_run, catchem_anthropogenic_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
      subroutine catchem_anthropogenic_wrapper_init()
      end subroutine catchem_anthropogenic_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_anthropogenic_wrapper_finalize Argument Table
!!
      subroutine catchem_anthropogenic_wrapper_finalize()
      end subroutine catchem_anthropogenic_wrapper_finalize

!> \defgroup catchem_anthropogenic_group CATChem anthro wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_anthropogenic_wrapper CATChem anthro  wrapper Module  
!> \ingroup catchem_anthropogenic_group
!! This is the CATChem anthro  wrapper Module
!! \section arg_table_catchem_anthropogenic_wrapper_run Argument Table
!! \htmlinclude catchem_anthropogenic_wrapper_run.html
!!
!>\section catchem_anthropogenic_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_anthropogenic_wrapper_run(me, master, &
                   im, kte, kme, ktau, dt,               &
                   jdate, garea, rlat, rlon,    &
                   pr3d, ph3d,phl3d, prl3d, tk3d, spechum,emi_in,                       &
                   ntrac,ntso2,ntsulf,ntpp25,ntbc1,ntoc1,ntpp10,                        &
                   ntextinct, &
                   ntisop,ntno,ntno2,ntco,ntc2h4,ntc2h6,ntc3h6,ntc3h8,                 &
                   ntnh3,ntch2o,ntch4,ntc4h10,ntc10h16,ntch3oh,ntc2h5oh,               &
                   ntch3coch3,nth2,nte90,gaschem_opt,do_am4chem, &
                   ntchs, ntche,chemic_in,chem_in_opt, &
                   nvar_emi,nvar_chemic, &
                   gq0,qgrs,abem,bioem,antem,chem_opt_in,kemit_in,pert_scale_anthro,                &
                   emis_amp_anthro,do_sppt_emis,sppt_wts,errmsg,errflg)

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master

    integer,        intent(in) :: im,kte,kme,ktau, jdate(8)
    integer,        intent(in) :: ntrac,ntchs,ntche
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf
    integer,        intent(in) :: ntisop,ntno,ntno2,ntco,ntc2h4,ntc2h6,ntc3h6,ntc3h8,  &
                                  ntnh3,ntch2o,ntch4,ntc4h10,ntc10h16,ntch3oh,ntc2h5oh,&
                                  ntch3coch3,nth2,nte90,ntextinct
    integer,        intent(in) :: nvar_emi,nvar_chemic
    real(kind_phys),intent(in) :: dt, emis_amp_anthro, pert_scale_anthro
    real(kind_phys),optional, intent(in) :: sppt_wts(:,:)
    logical,        intent(in) :: do_sppt_emis

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    !JianHe, we may need more tracers in the emissions input
    real(kind_phys), dimension(im, nvar_emi), intent(in) :: emi_in
    real(kind_phys), dimension(im), intent(in) :: garea, rlat,rlon
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d, spechum
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,12), intent(inout) :: abem
    integer,         intent(in) :: chem_opt_in, kemit_in, gaschem_opt,chem_in_opt
    logical,         intent(in) :: do_am4chem
    real(kind_phys), optional, intent(inout) :: bioem(:,:)
    real(kind_phys), optional, intent(inout) :: antem(:,:) !4 species for now
    real(kind_phys), optional, intent(in) :: chemic_in(:,:,:)  !JianHe

    character(len=*),intent(out) :: errmsg
    integer,         intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,       &
                     p_phy, z_at_w, dz8w, p8w, rho_phy

    real(kind_phys), dimension(ims:im, jms:jme) :: xlat, xlong, dxy

!>- chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    integer :: ide, ime, ite, kde
    real(kind_phys), dimension(ims:im, kms:kemit, jms:jme, 1:num_emis_ant) :: emis_ant
    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg
    real(kind_phys), parameter :: ugkg = 1.e-09_kind_phys !lzhang

    !JianHe: currently only for isop + c10h16
    real(kind_phys), dimension(ims:im, kms:kemit_in, jms:jme, 1:num_emis_bio) :: emis_bio
    integer, parameter :: nvl_am4ic = 49 ! JianHe:number of input levels from AM4 chemic

    integer :: i, j, jp, k, kp, n, nt
    real(kind_phys) :: random_factor(ims:im,jms:jme)
  

    errmsg = ''
    errflg = 0

    chem_opt          = chem_opt_in
    kemit             = kemit_in

    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    if(do_sppt_emis) then
      random_factor(:,jms) = pert_scale_anthro*max(min(1+(sppt_wts(:,kme/2)-1)*emis_amp_anthro,2.0),0.0)
    else
      random_factor = 1.0
    endif

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

!>- get ready for chemistry run
    ! -- anthropogenic emission
    call catchem_prep_anthropogenic(                                   &
        me, master, &
        ktau,dtstep,                                                    &
        jdate,garea,rlat,rlon, &
        xlat,xlong,dxy,  &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,emi_in,                      &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,z_at_w,                        &
        ntisop,ntno,ntno2,ntco,ntc2h4,ntc2h6,ntc3h6,ntc3h8,             &
        ntnh3,ntch2o,ntch4,ntc4h10,ntc10h16,ntch3oh,ntc2h5oh,           &
        ntch3coch3,nth2,nte90,gaschem_opt,do_am4chem,chem_opt,  &
        ntso2,ntsulf,ntpp25,ntbc1,ntoc1,ntpp10,ntrac,gq0,               &
        num_chem, num_ebu_in,num_emis_ant,emis_ant,                     &
        num_emis_bio,emis_bio, &
        ppm2ugkg,chem,random_factor,              &
        chem_in_opt, ntchs, ntche, &
        nvar_emi, chemic_in, nvl_am4ic, nvar_chemic,    &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       if (gaschem_opt == 1 .or. do_am4chem) then  !JianHe: for AM4
         do n = 1, num_chem
           nt = n+ntchs-1
           gq0(i,k,nt) = ppm2ugkg(n) * max(epsilc,chem(i,k,1,n))
         enddo
       else
         gq0(i,k,ntso2  )=ppm2ugkg(p_so2   ) * max(epsilc,chem(i,k,1,p_so2))
         gq0(i,k,ntsulf )=ppm2ugkg(p_sulf  ) * max(epsilc,chem(i,k,1,p_sulf))
         gq0(i,k,ntpp25 )=ppm2ugkg(p_p25   ) * max(epsilc,chem(i,k,1,p_p25))
         gq0(i,k,ntbc1  )=ppm2ugkg(p_bc1   ) * max(epsilc,chem(i,k,1,p_bc1))
         gq0(i,k,ntoc1  )=ppm2ugkg(p_oc1   ) * max(epsilc,chem(i,k,1,p_oc1))
         gq0(i,k,ntpp10 )=ppm2ugkg(p_p10   ) * max(epsilc,chem(i,k,1,p_p10))
       end if
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       qgrs(i,k,ntso2  )=gq0(i,k,ntso2  )
       qgrs(i,k,ntsulf )=gq0(i,k,ntsulf )
       qgrs(i,k,ntpp25 )=gq0(i,k,ntpp25 )
       qgrs(i,k,ntbc1  )=gq0(i,k,ntbc1  )
       qgrs(i,k,ntoc1  )=gq0(i,k,ntoc1  )
       qgrs(i,k,ntpp10 )=gq0(i,k,ntpp10 )
       if (gaschem_opt == 1 .or. do_am4chem) then !JianHe: For AM4
         do n = ntchs,ntche  ! only prog chem
           qgrs(i,k,n) = gq0(i,k,n)
         enddo
       end if
     enddo
    enddo

    abem(:,1)=ugkg*emis_ant(:,kts,1,p_e_bc )
    abem(:,2)=ugkg*emis_ant(:,kts,1,p_e_oc )
    abem(:,3)=ugkg*emis_ant(:,kts,1,p_e_so2)

#ifdef AM4_CHEM
    bioem(:,2)=emis_bio(:,kts,1,p_ebio_isop)
    bioem(:,3)=emis_bio(:,kts,1,p_ebio_c10h16)

    antem(:,1)=emis_ant(:,kts,1,p_e_co)
    antem(:,2)=emis_ant(:,kts,1,p_e_ch4)
    antem(:,3)=emis_ant(:,kts,1,p_e_nh3)
    antem(:,4)=emis_ant(:,kts,1,p_e_so2)
#endif

!
   end subroutine catchem_anthropogenic_wrapper_run
!> @}

   subroutine catchem_prep_anthropogenic(                               &
        me, master, &
        ktau,dtstep,                                                     &
        jdate,garea,rlat,rlon, &
        xlat,xlong,dxy,  &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,emi_in,                       &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,z_at_w,                         &
        ntisop,ntno,ntno2,ntco,ntc2h4,ntc2h6,ntc3h6,ntc3h8,             &
        ntnh3,ntch2o,ntch4,ntc4h10,ntc10h16,ntch3oh,ntc2h5oh,           &
        ntch3coch3,nth2,nte90,gaschem_opt,do_am4chem,chem_opt,           &
        ntso2,ntsulf,ntpp25,ntbc1,ntoc1,ntpp10,ntrac,gq0,                &
        num_chem, num_ebu_in,num_emis_ant,emis_ant,                      &
        num_emis_bio,emis_bio,ppm2ugkg,chem,random_factor,    &
        chem_in_opt, ntchs, ntche, &
        nvar_emi,chemic_in,nvl_am4ic, nvar_chemic,     &
        ids,ide, jds,jde, kds,kde,                                       &
        ims,ime, jms,jme, kms,kme,                                       &
        its,ite, jts,jte, kts,kte)

    integer, intent (in) :: me
    integer, intent (in) :: master

    !Chem input configuration
    integer, intent(in) :: ktau,jdate(8)
    real(kind=kind_phys), intent(in) :: dtstep

    !Stochastic physics variables
    real(kind_phys), intent(in) :: random_factor(ims:ime,jms:jme)

    !FV3 input variables
    integer, intent(in) :: ntrac,ntchs, ntche
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer, intent(in) :: ntsulf
    integer, intent(in) :: ntisop,ntno,ntno2,ntco,ntc2h4,ntc2h6,ntc3h6,ntc3h8,  &
                           ntnh3,ntch2o,ntch4,ntc4h10,ntc10h16,ntch3oh,ntc2h5oh,&
                           ntch3coch3,nth2,nte90
    real(kind=kind_phys), dimension(ims:ime), intent(in) :: garea, rlat, rlon
    real(kind=kind_phys), dimension(ims:ime, nvar_emi),   intent(in) :: emi_in
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: phl3d,tk3d,prl3d,spechum
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0

    integer,           intent(in) :: gaschem_opt,chem_opt,chem_in_opt
    logical,           intent(in) :: do_am4chem


    !chemical IC
    integer,intent(in) :: nvl_am4ic, nvar_chemic, nvar_emi
    real(kind=kind_phys), optional, intent(in) :: chemic_in(:,:,:)

    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_ebu_in,num_emis_ant
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(1:num_chem), intent(in) :: ppm2ugkg

    real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, 1:num_emis_ant), intent(inout) :: emis_ant

    integer,intent(in) ::  num_emis_bio  ! offline biogenic emis
    real(kind_phys), dimension(ims:ime, kms:kemit, jms:jme, 1:num_emis_bio), intent(inout) :: emis_bio

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, p_phy, rho_phy, dz8w, p8w
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) :: xlat, xlong, dxy
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, jms:jme, 1:num_emis_ant) :: emiss_ab
    real(kind_phys), dimension(ims:ime, jms:jme, 1:num_emis_bio) :: emiss_bio
    real(kind_phys), parameter :: frac_so2_ant = 1.0_kind_phys  ! antropogenic so2 fraction

!>- volcanic stuff
    integer ::ko,k_final,k_initial,kl,kk4,curr_hours,curr_secs,curr_day,curr_mth,curr_yr
    integer :: ivolcano,num_emis_voll
    real(kind_phys) :: x1,ashz_above_vent,mindist,currdist
    real(kind_phys), DIMENSION (kms:kme) :: vert_mass_dist
    real(kind_phys) :: eh,maxth,lon_vol,lat_vol,gmm,erup_beg,erup_dt,erup_end
    real(kind_phys), dimension(ims:ime, kms:kme,   jms:jme, num_emis_vol) :: emis_vol

    ! -- volcano ashes parameters
    real(kind_phys), dimension(7) :: h = (/ (240., i = 1, 6), huge(1.0) /)
    real(kind_phys), dimension(6) :: emiss_ash_table = (/  5834.,  3834., 5834.,  3334.,  3334.,  2334. /)
    real(kind_phys), dimension(6) :: eh_table        = (/ 3.11e5, 3.87e4, 3.11e5, 2.17e4, 2.17e4, 4.93e3 /)
    real(kind_phys), parameter :: percen_mass_umbrel = 0.75
    real(kind_phys), parameter :: base_umbrel        = 0.25    ! fraction
    real(kind_phys), parameter :: base_umbrel2       = 1.0     ! evenly distribution

    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour,igbox,jgbox,nt
    real(kind_phys), DIMENSION (ims:ime,jms:jme) :: so2_mass,emiss_ash_dtt
    real(kind_phys), dimension(ims:ime, jms:jme) :: emiss_ash_mass
    real(kind_phys), dimension(ims:ime, jms:jme) :: emiss_ash_height
    real(kind_phys), dimension(ims:ime, jms:jme) :: emiss_ash_dt
    real(kind_phys) ::  factor,factor2

    !--chemical IC array
    real(kind_phys), dimension(ims:ime, jms:jme, 1:nvl_am4ic, 1:nvar_chemic) :: chemic
    real(kind_phys), dimension(ims:ime, kms:kte, jms:jme, 1:nvar_chemic) :: chemic_interp
    real(kind_phys), dimension(1:nvl_am4ic) :: p_am4ic
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_mb

    p_am4ic = (/ 997.948596287477, 992.786752544069, 985.401223968179, &
    975.648851831372, 963.308651425801, 948.051493845347, 929.509176539644, &
    907.245072729549, 880.739861311195, 849.397132370557, 812.629155924972, &
    769.874744124755, 720.696844719115, 664.880933508635, 602.621201211844, &
    535.059070791703, 464.941771304376, 396.36295331148, 333.155174318448,  &
    277.364765336414, 229.263219871856, 188.254870627357, 153.542046894013, &
    124.359351881068, 99.9990654267351, 79.8128107819713, 63.2117171006669, &
    49.6656638058622, 38.7016977185877, 29.9017286729427, 22.8996206774882, &
    17.3778029011962, 13.0635264807004, 9.72489105230714, 7.16675888193994, &
    5.2266642666333, 3.77081383122845, 2.6902576168577, 1.89729594414021,   &
    1.32216980330045, 0.910066120368144, 0.618455975354221,   &
    0.414767919596273, 0.274388498599214, 0.178752860685895,  &
    0.11361803458164, 0.0686551589322532, 0.0380102383882785, &
    0.0171052512649617 /)

!    print*,'hli into volc'
    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    dxy            = 0._kind_phys
    chemic_interp  = 0._kind_phys !
    chemic         = 0._kind_phys

    emiss_ash_dtt  = 0._kind_phys
    num_emis_voll  =0._kind_phys
    so2_mass       = 0._kind_phys
    vert_mass_dist = 0._kind_phys

    lon_vol=-10000
    lat_vol=-10000

    curr_yr  = jdate(1)
    curr_mth = jdate(2)
    curr_day = jdate(3)
    ! -- initialize local arrays
    !idate=20220115  
    if (curr_yr==2022.and. curr_mth==1.and. curr_day==15) then
    erup_beg=3*3600.
    erup_dt=18000. !second
    erup_end=erup_beg+erup_dt
    gmm=float(ktau)*dtstep
    num_emis_voll =4
    !-- locatoin of volcano
    lon_vol=-175.38+360
    lat_vol=-20.57
    endif

    mindist=1.e9
    igbox=-1
    jgbox=-1

    ! -- sanity check for volcanic emissions
    if (num_emis_voll > 0) then
      select case (chem_opt)
        case (316)
          jmax = 10
        case (317, 502)
          jmax = 4
        case (CHEM_OPT_GOCART)
          jmax = 4
        case default
          jmax = num_emis_voll
      end select
      !if (num_emis_voll /= jmax) then
      !  call chem_rc_set(CHEM_RC_FAILURE, &
      !    msg="Inconsistent volcanic ash settings", &
      !    file=__FILE__, line=__LINE__, rc=rc)
      !  return
      !end if
    end if

    ! -- initialize fire emissions
    if (ktau <= 1) then
      emis_ant = 0.
      emis_bio = 0.
      emis_vol = 0.
    end if
    
    do i=its,ite
     dxy  (i,1)=garea(i)
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
    enddo

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

    p_mb = 0.01*p_phy ! Pa to mb

!lzhang

    ! -- volcanics
   if (num_emis_voll >0) then
    do j=jts,jte
      do i=its,ite
        if(abs(lon_vol-xlong(i,j)) .gt. 2. .or.&
           abs(lat_vol-xlat(i,j)) .gt. 2. ) cycle

           currdist=sqrt((lon_vol-xlong(i,j))**2+&
                         (lat_vol-xlat(i,j))**2)

        if(currdist<mindist) then
           mindist=currdist
           igbox=i ; jgbox=j
        endif
      enddo
    enddo

    if (igbox .ne. -1 .and. jgbox .ne.-1) then
       emiss_ash_dtt(igbox,jgbox)     = 3600.*5
       !emiss_ash_height(igbox,jgbox)  = 18000.
       emiss_ash_height(igbox,jgbox)  = 12000.
    endif
    endif !num_emis_voll 

    ! -- anthropagenic
    emiss_ab  = 0.   ! background
    emiss_bio = 0.   ! biogenic
    do j=jts,jte
     do i=its,ite
      emiss_ab(i,j,p_e_bc)   =emi_in(i,1)*random_factor(i,j)
      emiss_ab(i,j,p_e_oc)   =emi_in(i,2)*random_factor(i,j)
      emiss_ab(i,j,p_e_sulf) =emi_in(i,3)*random_factor(i,j)
      emiss_ab(i,j,p_e_pm_25)=emi_in(i,4)*random_factor(i,j)
      emiss_ab(i,j,p_e_so2)  =emi_in(i,5)*random_factor(i,j)
      emiss_ab(i,j,p_e_pm_10)=emi_in(i,6)*random_factor(i,j)
     enddo
    enddo

    factor=0.
    k=kts
    if (p_bc2 > 1) then
      do j=jts,jte
        do i=its,ite
          emis_ant(i,k,j,p_e_bc)=emiss_ab(i,j,p_e_bc)
          emis_ant(i,k,j,p_e_oc)=emiss_ab(i,j,p_e_oc)
          emis_ant(i,k,j,p_e_sulf)=emiss_ab(i,j,p_e_sulf)
          emis_ant(i,k,j,p_e_so2)=frac_so2_ant * emiss_ab(i,j,p_e_so2)
          emis_ant(i,k,j,p_e_dms)= 0. !emiss_ab(j,p_e_dms)
          emis_ant(i,k,j,p_e_pm_25)=emiss_ab(i,j,p_e_pm_25)
          emis_ant(i,k,j,p_e_pm_10)=emiss_ab(i,j,p_e_pm_10)
        enddo
      enddo
    endif
  
 
#ifdef AM4_CHEM
    ! -- anthropagenic
    do j=jts,jte
     do i=its,ite
          !JH: need update later, in emi_in, more tracers
          emiss_ab(i,j,p_e_no) =emi_in(i,11)*random_factor(i,j)
          emiss_ab(i,j,p_e_no2) =emi_in(i,12)*random_factor(i,j)
          emiss_ab(i,j,p_e_nh3) =emi_in(i,13)*random_factor(i,j)
          emiss_ab(i,j,p_e_co) =emi_in(i,14)*random_factor(i,j)
          emiss_ab(i,j,p_e_ch4) =emi_in(i,15)*random_factor(i,j)
          emiss_ab(i,j,p_e_ch2o) =emi_in(i,16)*random_factor(i,j)
          emiss_ab(i,j,p_e_c2h4) =emi_in(i,17)*random_factor(i,j)
          emiss_ab(i,j,p_e_c2h6) =emi_in(i,18)*random_factor(i,j)
          emiss_ab(i,j,p_e_c3h6) =emi_in(i,19)*random_factor(i,j)
          emiss_ab(i,j,p_e_c3h8) =emi_in(i,20)*random_factor(i,j)
          emiss_ab(i,j,p_e_c4h10) =emi_in(i,21)*random_factor(i,j)
          emiss_ab(i,j,p_e_isop) =emi_in(i,22)*random_factor(i,j)
          emiss_ab(i,j,p_e_c10h16) =emi_in(i,23)*random_factor(i,j)
          emiss_ab(i,j,p_e_ch3oh) =emi_in(i,24)*random_factor(i,j)
          emiss_ab(i,j,p_e_c2h5oh) =emi_in(i,25)*random_factor(i,j)
          emiss_ab(i,j,p_e_ch3coch3) =emi_in(i,26)*random_factor(i,j)
          emiss_ab(i,j,p_e_h2) =emi_in(i,27)*random_factor(i,j)
          emiss_ab(i,j,p_e_e90) =emi_in(i,28)*random_factor(i,j)
          !biogenic
          emiss_bio(i,j,p_ebio_isop) =emi_in(i,29)*random_factor(i,j)
          emiss_bio(i,j,p_ebio_c10h16) =emi_in(i,30)*random_factor(i,j)
     enddo
    enddo

    k=kts
    do j=jts,jte
      do i=its,ite
            emis_ant(i,k,j,p_e_no) =emiss_ab(i,j,p_e_no)
            emis_ant(i,k,j,p_e_no2) =emiss_ab(i,j,p_e_no2)
            emis_ant(i,k,j,p_e_nh3) =emiss_ab(i,j,p_e_nh3)
            emis_ant(i,k,j,p_e_co) =emiss_ab(i,j,p_e_co)
            emis_ant(i,k,j,p_e_ch4) =emiss_ab(i,j,p_e_ch4)
            emis_ant(i,k,j,p_e_ch2o) =emiss_ab(i,j,p_e_ch2o)
            emis_ant(i,k,j,p_e_c2h4) =emiss_ab(i,j,p_e_c2h4)
            emis_ant(i,k,j,p_e_c2h6) =emiss_ab(i,j,p_e_c2h6)
            emis_ant(i,k,j,p_e_c3h6) =emiss_ab(i,j,p_e_c3h6)
            emis_ant(i,k,j,p_e_c3h8) =emiss_ab(i,j,p_e_c3h8)
            emis_ant(i,k,j,p_e_c4h10) =emiss_ab(i,j,p_e_c4h10)
            emis_ant(i,k,j,p_e_isop) =emiss_ab(i,j,p_e_isop)
            emis_ant(i,k,j,p_e_c10h16) =emiss_ab(i,j,p_e_c10h16)
            emis_ant(i,k,j,p_e_ch3oh) =emiss_ab(i,j,p_e_ch3oh)
            emis_ant(i,k,j,p_e_c2h5oh) =emiss_ab(i,j,p_e_c2h5oh)
            emis_ant(i,k,j,p_e_ch3coch3) =emiss_ab(i,j,p_e_ch3coch3)
            emis_ant(i,k,j,p_e_h2) =emiss_ab(i,j,p_e_h2)
            emis_ant(i,k,j,p_e_e90) =emiss_ab(i,j,p_e_e90)

            emis_bio(i,k,j,p_ebio_isop) =emiss_bio(i,j,p_ebio_isop)
            emis_bio(i,k,j,p_ebio_c10h16) =emiss_bio(i,j,p_ebio_c10h16)
        enddo
      enddo
#endif


    do k=kms,kte
     do i=ims,ime
       !JianHe: use num_chem loop here to read from tracer arrary to chem array
       !num_chem = ntchm + ndchm (all chem prog+diag tracers in tracer array)
       !tracer index starting from p_so2 (ntso2=ntchs)
       !in the future, we do not need specify tracer index in catchem_config. 
       !For now, it should be in consistent order in chem array and tracer array
       if (gaschem_opt == 1 .or. do_am4chem) then
         do n = 1,num_chem
            nt = n+ntchs-1  ! tracer index in tracer array
            chem(i,k,jts,n) = max(epsilc,gq0(i,k,nt)/ppm2ugkg(n))
         enddo
       else
         chem(i,k,jts,p_so2   )=max(epsilc,gq0(i,k,ntso2  )/ppm2ugkg(p_so2))
         chem(i,k,jts,p_sulf  )=max(epsilc,gq0(i,k,ntsulf )/ppm2ugkg(p_sulf))
         chem(i,k,jts,p_p25   )=max(epsilc,gq0(i,k,ntpp25 )/ppm2ugkg(p_p25))
         chem(i,k,jts,p_bc1   )=max(epsilc,gq0(i,k,ntbc1  )/ppm2ugkg(p_bc1))
         chem(i,k,jts,p_oc1   )=max(epsilc,gq0(i,k,ntoc1  )/ppm2ugkg(p_oc1))
         chem(i,k,jts,p_p10   )=max(epsilc,gq0(i,k,ntpp10 )/ppm2ugkg(p_p10))
       end if
     enddo
    enddo

    !JianHe: let's include chemical ICs here before adding emissions
    if(ktau.le.1)then
      if (chem_in_opt == 1 ) then
#ifdef AM4_CHEM
        !if (gaschem_opt == 1 .or. do_am4chem) then
          do i=its,ite
            do k=1,nvl_am4ic
              do nt=1,nvar_chemic
              !if (me == master) then
              !   write (*, *) "Read chemic: ", nt, chemic_in(i,k,nt)
              !endif
                chemic(i,1,k,nt)=max(epsilc,chemic_in(i,k,nt))
              enddo
            enddo
          enddo

          do i=its,ite
            do j=jts,jte
              do k=kts,kte
                call vinterp(chemic(i,j,1:nvl_am4ic,1:nvar_chemic),nvar_chemic,nvl_am4ic, &
                    p_am4ic,p_mb(i,k,j), chemic_interp(i,k,j,1:nvar_chemic))
              enddo
              chem(i,kts:kte,j,p_co) = chemic_interp(i,kts:kte,j,1)
              chem(i,kts:kte,j,p_o3) = chemic_interp(i,kts:kte,j,2)
              chem(i,kts:kte,j,p_n2o) = chemic_interp(i,kts:kte,j,3)
              chem(i,kts:kte,j,p_no) = chemic_interp(i,kts:kte,j,4)
              chem(i,kts:kte,j,p_no2) = chemic_interp(i,kts:kte,j,5)
              chem(i,kts:kte,j,p_hno3) = chemic_interp(i,kts:kte,j,6)
              chem(i,kts:kte,j,p_n2o5) = chemic_interp(i,kts:kte,j,7)
              chem(i,kts:kte,j,p_ch4) = chemic_interp(i,kts:kte,j,8)
              chem(i,kts:kte,j,p_pan) = chemic_interp(i,kts:kte,j,9)
              chem(i,kts:kte,j,p_c2h6) = chemic_interp(i,kts:kte,j,10)
              chem(i,kts:kte,j,p_ch3coch3) = chemic_interp(i,kts:kte,j,11)
              chem(i,kts:kte,j,p_hcl) = chemic_interp(i,kts:kte,j,12)
              chem(i,kts:kte,j,p_hocl) = chemic_interp(i,kts:kte,j,13)
              chem(i,kts:kte,j,p_clono2) = chemic_interp(i,kts:kte,j,14)
              chem(i,kts:kte,j,p_clo) = chemic_interp(i,kts:kte,j,15)
              chem(i,kts:kte,j,p_hobr) = chemic_interp(i,kts:kte,j,16)
              chem(i,kts:kte,j,p_hbr) = chemic_interp(i,kts:kte,j,17)
              chem(i,kts:kte,j,p_brono2) = chemic_interp(i,kts:kte,j,18)
              chem(i,kts:kte,j,p_bro) = chemic_interp(i,kts:kte,j,19)
              chem(i,kts:kte,j,p_age) = chemic_interp(i,kts:kte,j,20)
              chem(i,kts:kte,j,p_cl) = chemic_interp(i,kts:kte,j,21)
              chem(i,kts:kte,j,p_cl2) = chemic_interp(i,kts:kte,j,22)
              chem(i,kts:kte,j,p_cl2o2) = chemic_interp(i,kts:kte,j,23)
              chem(i,kts:kte,j,p_cl2o2) = chemic_interp(i,kts:kte,j,23)
              chem(i,kts:kte,j,p_br) = chemic_interp(i,kts:kte,j,24)
              chem(i,kts:kte,j,p_brcl) = chemic_interp(i,kts:kte,j,25)
              chem(i,kts:kte,j,p_extinction) = chemic_interp(i,kts:kte,j,26)
            enddo
          enddo
        ! endif !(gaschem_opt)

#else 
      !placeholer
#endif
    
      else !chem_in_opt = 0
          do i=its,ite
            do j=jts,jte
              chem(i,kts:kte,j,ntchs:ntche) = max(epsilc, chem(i,kts:kte,j,ntchs:ntche))
            enddo
          enddo
      endif !(chem_in_opt)
    endif !(ktau<=1)

    !
    ! -- gocart background fields only if gocart is called
    !
!   emis_ant=0.
    nv=1
    k=kts
    factor2=0.
    factor=0.
    if (p_bc2 > 1)then
      if (chem_opt == CHEM_OPT_GOCART) then
        do j=jts,jte
          do i=its,ite
            factor=dtstep*rri(i,k,j)/dz8w(i,k,j)
            factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
            chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+emis_ant(i,k,j,p_e_bc)*factor
            chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+emis_ant(i,k,j,p_e_oc)*factor
            chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+emis_ant(i,k,j,p_e_pm_25)*factor
            chem(i,k,j,p_p10)=chem(i,k,j,p_p10)+emis_ant(i,k,j,p_e_pm_10)*factor
            chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+emis_ant(i,k,j,p_e_sulf)*factor
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+emis_ant(i,k,j,p_e_so2)*factor2
          enddo
        enddo
      endif

     !JianHe: convert aerosol to ug/kg and gas to ppm
     !if (gaschem_opt == 1 .or. do_am4chem) then
#ifdef AM4_CHEM
        do j=jts,jte
          do i=its,ite
            factor=dtstep*rri(i,k,j)/dz8w(i,k,j)
            factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))

            chem(i,k,j,p_bc1)=chem(i,k,j,p_bc1)+emis_ant(i,k,j,p_e_bc)*factor    !ug/kg
            chem(i,k,j,p_oc1)=chem(i,k,j,p_oc1)+emis_ant(i,k,j,p_e_oc)*factor
            chem(i,k,j,p_p25)=chem(i,k,j,p_p25)+emis_ant(i,k,j,p_e_pm_25)*factor
            chem(i,k,j,p_p10)=chem(i,k,j,p_p10)+emis_ant(i,k,j,p_e_pm_10)*factor
            chem(i,k,j,p_sulf)=chem(i,k,j,p_sulf)+emis_ant(i,k,j,p_e_sulf)*factor
            chem(i,k,j,p_so2)=chem(i,k,j,p_so2)+emis_ant(i,k,j,p_e_so2)*factor2
            chem(i,k,j,p_no) =chem(i,k,j,p_no) +emis_ant(i,k,j,p_e_no )*factor2
            chem(i,k,j,p_no2) =chem(i,k,j,p_no2) +emis_ant(i,k,j,p_e_no2)*factor2
            chem(i,k,j,p_nh3) =chem(i,k,j,p_nh3) +emis_ant(i,k,j,p_e_nh3)*factor2
            chem(i,k,j,p_co) =chem(i,k,j,p_co) +emis_ant(i,k,j,p_e_co )*factor2
            chem(i,k,j,p_ch4) =chem(i,k,j,p_ch4) +emis_ant(i,k,j,p_e_ch4)*factor2
            chem(i,k,j,p_ch2o) =chem(i,k,j,p_ch2o) +emis_ant(i,k,j,p_e_ch2o)*factor2
            chem(i,k,j,p_c2h4) =chem(i,k,j,p_c2h4) +emis_ant(i,k,j,p_e_c2h4)*factor2
            chem(i,k,j,p_c2h6) =chem(i,k,j,p_c2h6) +emis_ant(i,k,j,p_e_c2h6)*factor2
            chem(i,k,j,p_c3h6) =chem(i,k,j,p_c3h6) +emis_ant(i,k,j,p_e_c3h6)*factor2
            chem(i,k,j,p_c3h8) =chem(i,k,j,p_c3h8) +emis_ant(i,k,j,p_e_c3h8)*factor2
            chem(i,k,j,p_c4h10) =chem(i,k,j,p_c4h10) +emis_ant(i,k,j,p_e_c4h10)*factor2
            chem(i,k,j,p_isop) =chem(i,k,j,p_isop) +emis_ant(i,k,j,p_e_isop)*factor2 + &
                                emis_bio(i,k,j,p_ebio_isop)*factor2
            chem(i,k,j,p_c10h16) =chem(i,k,j,p_c10h16)+emis_ant(i,k,j,p_e_c10h16 )*factor2 + &
                                  emis_bio(i,k,j,p_ebio_c10h16)*factor2
            chem(i,k,j,p_ch3oh) =chem(i,k,j,p_ch3oh) +emis_ant(i,k,j,p_e_ch3oh)*factor2
            chem(i,k,j,p_c2h5oh) =chem(i,k,j,p_c2h5oh)+emis_ant(i,k,j,p_e_c2h5oh )*factor2
            chem(i,k,j,p_ch3coch3) =chem(i,k,j,p_ch3coch3)+emis_ant(i,k,j,p_e_ch3coch3 )*factor2
            chem(i,k,j,p_h2) =chem(i,k,j,p_h2) +emis_ant(i,k,j,p_e_h2 )*factor2
            chem(i,k,j,p_e90) =chem(i,k,j,p_e90) +emis_ant(i,k,j,p_e_e90)*factor2
          enddo
        enddo
#endif

    else if (p_tr2 > 1)then    !co2 here
      do j=jts,jte
        do i=its,ite
!           factor2=dtstep*rri(i,k,j)/dz8w(i,k,j)
          !factor2=4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
          !chem(i,k,j,p_tr1)=chem(i,k,j,p_tr1)+emis_ant(i,k,j,p_e_tr1)*factor2
          !chem(i,k,j,p_tr2)=chem(i,k,j,p_tr2)+emis_ant(i,k,j,p_e_tr2)*factor2
        enddo
      enddo
    else if ((p_tr2 > 1) .and. (p_bc2 > 1))then
      !call chem_rc_set(CHEM_RC_FAILURE, msg="Inconsistent options detected.", &
      !  file=__FILE__, line=__LINE__, rc=rc)
      return
    endif

    curr_secs=ktau*ifix(dtstep)
    curr_hours=curr_secs/3600
!
!     do volcanoes if avaiable
!
!     if(chem_opt == 502 ) then
!       do j=jts,jte
!       if(emiss_ash_dt(j).le.0)CYCLE
!       emiss_ash_mass(j)=0.
!       emiss_ash_height(j)=0.
!       enddo
!
! default
!

    do j=jts,jte
      jp = j - jts + 1
      do i=its,ite
        ip = i - its + 1
        if (emiss_ash_dtt(i,j) > 0.) then
          !lzhang so2_mass(i,j)=1.5e4*3600.*1.e9/64./area(ip,jp)
         !so2_mass(i,j)=6.4e7*1.e9/64./area(ip,jp)
          !so2 unit should be mol/km2/hour
          so2_mass(i,j)=1.24e7*1.e9/64./garea(i) ! 1.24e7 is kg (5 hour total is 62 kilotons), garea is m2
          eh=2600.*(emiss_ash_height(i,j)*.0005)**4.1494
          emiss_ash_mass(i,j)=eh*1.e9/garea(i) !eh is kg
         !emiss_ash_mass(i,j)=eh*1.e9/area(ip,jp)
        end if
      enddo
    enddo

    ! -- real-time application, keeping eruption constant
!
!    if (ktau <= 2) then
      ! -- volcanic emissions
      if (num_emis_voll > 0) then

        emis_vol = 0._kind_phys
  !      if(curr_hours.eq.h1 .or. curr_hours.eq.h2 .or. curr_hours.eq.h3 &
  !         .or. curr_hours.eq.h4 .or. curr_hours.eq.h5 .or. curr_hours.eq.h6
  !         .or. h1.gt.239)then
  !         .or. curr_hours.eq.0)then
  !         if(chem_opt == 316 .or. chem_opt == 317 .or. chem_opt == 502) then

        do j=jts,jte
          do i=its,ite
            if (emiss_ash_dtt(i,j)     <= 0) cycle
            if (emiss_ash_height(i,j) <= 0) cycle
            ashz_above_vent=emiss_ash_height(i,j) +z_at_w(i,kts,j)
            do k=kte-1,kts,-1
              if (z_at_w(i,k,j) < ashz_above_vent)then
                k_final=k+1
                exit
              endif !inner
            enddo
            vert_mass_dist=0.
  !              k_initial=int((k_final+k_initial)*0.5)

            ! -- parabolic vertical distribution between k_initial and k_final
            kk4 = k_final-k_initial+2
            do ko=1,kk4-1
              kl=ko+k_initial-1
              vert_mass_dist(kl) = 6.*percen_mass_umbrel*float(ko)/float(kk4)**2 * (1. - float(ko)/float(kk4))
            enddo
            if (sum(vert_mass_dist(kts:kte)) /= percen_mass_umbrel) then
              x1= ( percen_mass_umbrel- sum(vert_mass_dist(kts:kte)))/float(k_final-k_initial+1)
              do ko=k_initial,k_final
                vert_mass_dist(ko) = vert_mass_dist(ko)+ x1 !- values between 0 and 1.
              enddo
                     !pause
            endif !inner
            !k_final > 0 .and. k_initial >

            ! -- linear detrainment from vent to base of umbrella
            do ko=1,k_initial-1
              vert_mass_dist(ko)=float(ko)/float(k_initial-1)
            enddo
            x1=sum(vert_mass_dist(1:k_initial-1))

            do ko=1,k_initial-1
              vert_mass_dist(ko)=(1.-percen_mass_umbrel)*vert_mass_dist(ko)/x1
            enddo

            select case (chem_opt)
              case (CHEM_OPT_GOCART, CHEM_OPT_GFDL_AM4)
                ! -- if applied to gocart we only need finest ash bins, we use
                ! the coarse one for so2
                do ko=1,k_final
                  emis_vol(i,ko,j,p_e_vash1)=vert_mass_dist(ko)*so2_mass(i,j)
                  emis_vol(i,ko,j,p_e_vash2)=.08*vert_mass_dist(ko)*emiss_ash_mass(i,j)
                  emis_vol(i,ko,j,p_e_vash3)=.05*vert_mass_dist(ko)*emiss_ash_mass(i,j)
                  emis_vol(i,ko,j,p_e_vash4)=.035*vert_mass_dist(ko)*emiss_ash_mass(i,j)
                enddo
              case default
                ! -- no default action
            end select

            do ko=k_final+1,kte
              emis_vol(i,ko,j,p_e_vash1)=0.
              emis_vol(i,ko,j,p_e_vash2)=0.
              emis_vol(i,ko,j,p_e_vash3)=0.
              emis_vol(i,ko,j,p_e_vash4)=0.
            enddo

          enddo
        enddo
      end if  !num_emis_vol

    ! -- add volcanic emissions
    if (num_emis_voll > 0) then

      select case (chem_opt)
        case (CHEM_OPT_GOCART, CHEM_OPT_GFDL_AM4)
          ! -- for gocart only lump ash into p25 and p10
          !if (num_emis_voll /= 4) then
          !  call chem_rc_set(CHEM_RC_FAILURE, msg="num_emis_vol must be 4", &
          !    file=__FILE__, line=__LINE__, rc=rc)
          !  return
          !end if
! now we got volcanoc emissions, they need to be added to chem array
!
!       write(message,'(" Do volcanic emissions ")') 
!       CALL WRF_MESSAGE (message)

      do j=jts,jte
      do i=its,ite
          ivolcano = 0


      if (gmm >=erup_beg .and. gmm<=erup_end) then
          ivolcano = 1
      endif

                   do k=kts,kte
                if (emiss_ash_dtt(i,j) <= 0.) cycle
                factor=float(ivolcano)*4.828e-4*dtstep*rri(i,k,j)/(60.*dz8w(i,k,j))
                factor2=float(ivolcano)*dtstep*rri(i,k,j)/dz8w(i,k,j)
                chem(i,k,j,p_p25)=chem(i,k,j,p_p25)                          &
                                 +emis_vol(i,k,j,p_e_vash4)*factor2
                chem(i,k,j,p_so2)=chem(i,k,j,p_so2)                          &
                                 +emis_vol(i,k,j,p_e_vash1)*factor
                chem(i,k,j,p_p10)=chem(i,k,j,p_p10)                          &
  !                              +.5* emis_vol(i,k,j,p_e_vash4)*factor2      &
                                 +1.* emis_vol(i,k,j,p_e_vash3)*factor2      &
                                 +.5* emis_vol(i,k,j,p_e_vash2)*factor2
              enddo
            enddo
          enddo
        case default
          ! -- volcanic emissions not included by default
      end select
    end if


  end subroutine catchem_prep_anthropogenic

!> @}
  end module catchem_anthropogenic_wrapper
