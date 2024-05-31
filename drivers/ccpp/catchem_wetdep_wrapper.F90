!>\file catchem_wetdep_wrapper.F90
!! This file is GSDChem large-scale wet deposition wrapper with CCPP coupling to FV3
!! Haiqin.Li@noaa.gov 06/2020
!! Kate.Zhang@noaa.gov 02/2023
!! Revision History:
!! 05/2023, Restructure for CATChem, Jian.He@noaa.gov
!! Kate.Zhang@noaa.gov 11/2023
!! Jian.He@noaa.gov: include AM4 wetdep, 05/2024

 module catchem_wetdep_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use catchem_constants,  only : WTMAIR, epsilc
   use wetdep_ls_mod
   use gocart_aerosols_mod
   use dust_data_mod
   use gocart_diag_mod
   use gfdl_wetdep_mod

   implicit none

   private

   public :: catchem_wetdep_wrapper_init, catchem_wetdep_wrapper_run, catchem_wetdep_wrapper_finalize

   logical :: is_initialized = .false.
   logical :: do_cv_wetdep = .false.
   real(kind_phys),parameter :: rhow = 1.0e3, rhor = 1.0e3, rhos = 1.0e2, rhog = 4.0e2

contains

!> \brief Brief description of the subroutine
!!
!> \section arg_table_catchem_wetdep_wrapper_init Argument Table
!! \htmlinclude catchem_wetdep_wrapper_init.html
!!
      subroutine catchem_wetdep_wrapper_init(me, master, tracer_names, &
                                           gaschem_opt, errmsg, errflg)

       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master
       integer, intent(in)  :: gaschem_opt

       character(len=32), intent(in) :: tracer_names(:)
       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (is_initialized) return

       if (gaschem_opt == 1) then
          if (me == master) then
            print *, 'Start GFDL AM4 wet depostion for both gas and aerosols'
          endif

          call gfdl_wetdep_init(me, master, tracer_names)

          is_initialized = .true.
       end if

      end subroutine catchem_wetdep_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_wetdep_wrapper_finalize Argument Table
!!
      subroutine catchem_wetdep_wrapper_finalize()
      end subroutine catchem_wetdep_wrapper_finalize

!> \defgroup catchem_group CATChem wetdep wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_wetdep_wrapper CATChem wetdep wrapper Module  
!> \ingroup catchem_wetdep_group
!! This is the CATChem wetdep wrapper Module
!! \section arg_table_catchem_wetdep_wrapper_run Argument Table
!! \htmlinclude catchem_wetdep_wrapper_run.html
!!
!>\section catchem_wetdep_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_wetdep_wrapper_run(me,master, &
                   im, kte, kme, ktau, dt,       &
                   imp_physics, imp_physics_gfdl, imp_physics_thompson, &
                   rain_cplchm, rainc_cpl,rlat,                         &
                   pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum,  &
                   w, dqdt, lmk, cld_frac, cnvc, & 
                   ntrac,ntchs, ntchm, ndchm, ntchmdiag,                &
                   ntcw, ntiw, pfi_lsan, pfl_lsan, &
                   ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                     &
                   ntbc1,ntbc2,ntoc1,ntoc2,                             &
                   ntss1,ntss2,ntss3,ntss4,ntss5,                       &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,      &
                   ntsoa, ntage, ntaoanh,tracer_names, &
                   gq0,qgrs,wetdpl,wdep,wetdep_ls_opt_in,gaschem_opt,      &
                   do_am4chem, errmsg,errflg)

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master

    integer,        intent(in) :: im,kte,kme,ktau,ntchmdiag, &
                                  ntcw,ntiw,ntchs,ntchm,ndchm
    integer,        intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    integer,        intent(in) :: ntage, ntaoanh,ntsoa
    real(kind_phys),intent(in) :: dt
    logical,        intent(in) :: do_am4chem

    character(len=32), intent(in) :: tracer_names(:)

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    real(kind_phys), dimension(im),     intent(in) :: rain_cplchm,rainc_cpl,rlat
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: phl3d, prl3d, tk3d,        &
                us3d, vs3d, spechum, w, dqdt
    integer,           intent(in) :: lmk
    real(kind_phys), dimension(im,lmk), intent(in) ::  cld_frac
    real(kind_phys), dimension(im,kte), intent(in) ::  &
                                           cnvc,pfi_lsan, pfl_lsan
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    real(kind_phys), dimension(im,ntchmdiag), intent(inout) :: wetdpl
    real(kind_phys), optional, intent(inout) :: wdep(:,:)
    integer,           intent(in) :: wetdep_ls_opt_in,gaschem_opt
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy, u_phy, v_phy,       &
                     p_phy, z_at_w, dz8w, p8w, t8w, rho_phy, vvel, dqdti
    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: z_phy, sphum,cldfrac

    real(kind_phys), dimension(ims:im, jms:jme) :: rcav, rnav,xlat
    real(kind_phys), dimension(ims:im, jms:jme) :: surfrain, surfsnow

!>- vapor & chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem )  :: var_rmv

    !JianHe: AM4 var
    real(kind_phys), dimension(ims:im, jms:jme, kts:kte) :: so2_so4_evap, so2_so4
    real(kind_phys), dimension(ims:im, jms:jme, kts:kte) :: wetdeptnd
    real(kind_phys), dimension(ims:im, jms:jme, 1:num_chem ) :: ls_wetdep, cv_wetdep
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte,1:num_chem ) :: trac_am4
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte,1:ntrac) :: trac_dt
    real(kind_phys), dimension(1:im,jms:jme,1:kte) :: t_am4,    &
                    p_am4, z_am4, sphum_am4, cldfrac_am4, dqdt_am4, &
                    cldamt_am4, drain, dsnow, pdel, f_snow_berg
    real(kind_phys), dimension(1:im,jms:jme,1:kme) :: z8w_am4, p8w_am4, &
                                           rain3d_am4, snow3d_am4 ! need in phalf level

    integer :: ide, ime, ite, kde

    real(kind_phys) :: dtstep, qls,qcv,delz
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    ! -- output tracers
    real(kind_phys), dimension(im, 1, ntchmdiag, 4) :: trdf


!>-- local variables
    integer :: i, j, jp, k, kp, n, nt
  

    errmsg = ''
    errflg = 0

    if (wetdep_ls_opt_in > 0) then
    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- initialize large-sacle wet depostion
!    if (ktau==1) then
!     call dep_wet_ls_init()
!    endif

    ! -- set control flags

    ! -- compute accumulated large-scale and convective rainfall since last call
    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

    ! -- compute incremental convective and large-scale rainfall
    do i=its,ite
     rcav(i,1)=max(rainc_cpl(i)*1000.              , 0.) ! meter to mm
     rnav(i,1)=max((rain_cplchm(i)-rainc_cpl(i))*1000., 0.) ! meter to mm
    enddo

!!!
!>- get ready for chemistry run
    call catchem_prep_wetdep(ktau,dtstep,                               &
        imp_physics, imp_physics_gfdl, imp_physics_thompson,            &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w, dqdt,           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,z_phy,dz8w,p8w,lmk,cld_frac,cnvc, &
        t8w,dqdti,cldfrac,z_at_w,vvel,rlat,xlat,                                &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,                                        &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                  &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                 &
        ntrac,ntchm,ndchm,ntchs,gq0,num_chem, num_moist,                &
        ppm2ugkg,moist,chem,gaschem_opt,                                &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

    ! -- ls wet deposition
    var_rmv(:,:,:)=0._kind_phys

    select case (wetdep_ls_opt_in)
      case (WDLS_OPT_GSD)
        do j=jts,jte
          do i=its,ite
            call wetdep_ls(dt,chem(i,:,j,:),rnav(i,j),moist(i,:,j,:),   &
                        rho_phy(i,:,j),var_rmv(i,j,:),xlat(i,j),        &
                        p_qc,p_qi,dz8w(i,:,j),vvel(i,:,j),              &
                        kms,kme,kts,kte)
          enddo
        end do

      case (WDLS_OPT_NGAC)
        do j=jts,jte
          do i=its,ite
            call WetRemovalGOCART(kts,kte, 1,1, dt,   &
                               var_rmv(i,j,:),chem(i,:,j,:),    &
                               p_phy(i,:,j),t_phy(i,:,j),       &
                               rho_phy(i,:,j),dqdti(i,:,j),     &
                               rcav(i,j),rnav(i,j),             &
                               kms,kme)
          enddo
        enddo

         !if (chem_rc_check(localrc, msg="Failure in NGAC wet removal scheme", &
         !  file=__FILE__, line=__LINE__, rc=rc)) return
#ifdef AM4_CHEM
    ! JianHe: we need initialize in the wrapper file 
    ! instead of in lower subroutines to avoid issue 
    ! for diagnostic output (required by debugger turned on)
    wdep = 0._kind_phys

      case (WDLS_OPT_AM4)
        do j = jts, jte
          do i = its, ite
           !JianHe: AM4 needs vertical level from model top to bottom
            rain3d_am4(i,j,1) = 0. ! we set 0 at the model top for now
            snow3d_am4(i,j,1) = 0. ! input are accumulated precip
            surfrain(i,j) = pfl_lsan(i,1)   ! actually not used in the new scheme
            surfsnow(i,j) = pfi_lsan(i,1)   ! actually not used in the new scheme

            do k = kts, kte
              kp = kte-k+1
              !trac_am4(i,j,kp,1:num_chem)=chem(i,k,j,1:num_chem)
              t_am4(i,j,kp) = t_phy(i,k,j)  !K
              p_am4(i,j,kp) = p_phy(i,k,j)  !Pa
              z_am4(i,j,kp) = z_phy(i,k,j)  !m
              !cldfrac_am4(i,j,kp) = cldfrac(i,k,j)
              !cldamt_am4(i,j,kp) = gq0(i,k,ntcw)+gq0(i,k,ntiw) ! kg/kg
              rain3d_am4(i,j,kp+1) = pfl_lsan(i,k) ! kg m-2 s-1  !JianHe: is this the correct variablie?
              snow3d_am4(i,j,kp+1) = pfi_lsan(i,k) ! AM4 needs var in phalf levels    

              do n = 1, num_chem
              ! We need to make sure the order of tracers
              ! consistent with field table
              ! trac_prog array stores prog met+chem tracers only
              ! gas vmr, aerosol mmr
                nt = n+ntchs-1  ! tracer index in the tracer array
                if ((nt == ntage) .or. (nt == ntaoanh) ) then
                  trac_am4(i,j,kp,n) = chem(i,k,j,n)
                else if ((nt==ntdust1) .or. (nt==ntdust2) .or. (nt==ntdust3) &
                  .or. (nt==ntdust4) .or. (nt==ntdust5) .or. (nt==ntss1) &
                  .or. (nt==ntss2) .or. (nt==ntss3) .or. (nt==ntss4) &
                  .or. (nt==ntss5) .or. (nt==ntoc1) .or. (nt==ntoc2) &
                  .or. (nt==ntbc1) .or. (nt==ntbc2) .or. (nt==ntsoa) &
                  .or. (nt==ntpp25) .or. (nt==ntpp10) .or. (nt==ntsulf)) then
                  trac_am4(i,j,kp,n) = chem(i,k,j,n)*1.e-9  ! ug/kg to kg/kg
                else
                  trac_am4(i,j,kp,n) = chem(i,k,j,n)*1.e-6  ! ppm to mol/mol
                endif
              end do

            end do

            do k = kts, kte+1
              kp = kte+1-k+1
              z8w_am4(i,j,kp) = z_at_w(i,k,j)
              p8w_am4(i,j,kp) = p8w(i,k,j)
            end do

            drain(i,j,1:kte)=rain3d_am4(i,j,2:kte+1)-rain3d_am4(i,j,1:kte) ! kg/m2/s
            dsnow(i,j,1:kte)=snow3d_am4(i,j,2:kte+1)-snow3d_am4(i,j,1:kte) ! kg/m2/s
            pdel(i,j,1:kte)=(p8w_am4(i,j,2:kte+1)-p8w_am4(i,j,1:kte))/g    ! kg/m2

            do k = kts, kte
              kp = kte-k+1
              !not used in the new scheme
              !not used in the new scheme
              dqdt_am4(i,j,kp) = (drain(i,j,kp)+dsnow(i,j,kp))/pdel(i,j,kp)! kg/kg/s, may include wv, liquid, ice
              qls = dqdt_am4(i,j,kp)*rho_phy(i,k,j)                        ! kg/m3/s
              if (qls > 0.) then
                !JianHe: cldfrac is not calculated from MP
                !fraction of grid box covered by precipitating clouds.
                !following WetRemovalGOCART
                ! F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
                !may need a better way to calculate cloud fraction for wetdep
                 cldfrac_am4(i,j,kp) = 1./(1.+1.*1.0e-4*5.0e-4/qls)
              else
                 cldfrac_am4(i,j,kp) = 0.
              endif
              !Need add cloud amount at end of timestep.
              !WetRemovalGOCART has its own way to calculate cldliq and cldice
              cldamt_am4(i,j,kp) = gq0(i,k,ntcw)+gq0(i,k,ntiw) ! kg/kg
              cldamt_am4(i,j,kp) = cldamt_am4(i,j,kp) + dqdt_am4(i,j,kp)*dt

            end do
          end do
        end do

        !----------------------------------------------------------------------
        !    initialize field to hold so2--> so4 tendency returned from
        !    subroutine wet_deposition.
        !----------------------------------------------------------------------
        so2_so4_evap(:,:,:) = 0.
        trac_dt(:,:,:,:) = 0.
        !----------------------------------------------------------------------
        !    loop over each tracer, skipping the cloud tracers, calling the
        !    wet deposition routine to compute the amount of that tracer removed 
        !    from the column in each layer (Tend_mp%wetdeptnd(i,j,k)) and in the
        !    column (Removal_mp%ls_wetdep(i,j,n)).
        !----------------------------------------------------------------------
        do n = 1,ntchm  ! for gas+aerosol, tracer index in chem array
          nt = ntchs+n-1   ! tracer index in tracer array

          wetdeptnd(:,:,:) = 0.0
          f_snow_berg(:,:,:) = 0.0 ! Thompson scheme does not address Bergeron
                                   ! Process, set to 0 for now; 
                                   ! could be calculated in morrison scheme 

          !Do large-scale first
          call gfdl_wet_deposition (        &
              me, master, tracer_names,&
              nt, t_am4, p_am4, p8w_am4,   &
              z_am4, z8w_am4, surfrain,  &
              surfsnow, dqdt_am4, cldamt_am4, &
              cldfrac_am4, f_snow_berg, rain3d_am4,  &
              snow3d_am4, trac_am4(:,:,:,n),  &
              wetdeptnd, 'lscale', its, jts, dt,  &
              sum_wdep_out=ls_wetdep(:,:,n),  &
              so2_so4_out=so2_so4(:,:,:))

!-----------------------------------------------------------------------
!    add the wet deposition tendency for the tracer to the accumulated
!    total tracer tendency at each level.
!-----------------------------------------------------------------------
          ! ug/kg/s for aerosol and ppm/s for gas
          ! output in positive
          trac_dt (:,:,:,n) = wetdeptnd(:,:,:)

          if (tracer_names(nt) .eq. 'so2') then
            so2_so4_evap = so2_so4
          end if
        enddo

!-----------------------------------------------------------------------
!    correct so2 and so4 tendency. so2 is converted to so4.
!-----------------------------------------------------------------------
        !JianHe: this is already handled in wet_deposition sub?
        trac_dt (:,:,:,p_so2) = trac_dt (:,:,:,p_so2) + so2_so4_evap
        trac_dt (:,:,:,p_so4) = trac_dt (:,:,:,p_so4) - so2_so4_evap

      if(do_cv_wetdep) then
        !For convective removal, need fix
!------------------------------------------------------------------------
!    initialize fields needed for call to wet deposition routine.
!------------------------------------------------------------------------
     ! f_snow_berg = 0.
     ! C2ls_mp%wet_data = 0.0
     ! C2ls_mp%cloud_frac = 0.1
     ! C2ls_mp%cloud_wet = 1.e-3

        do n = 1,ntchm  ! for gas+aerosol, tracer index in chem array
          nt = ntchs+n-1   ! tracer index in tracer array

          wetdeptnd(:,:,:) = 0.0
          f_snow_berg(:,:,:) = 0.0 ! Thompson scheme does not address Bergeron
                                   ! Process, set to 0 for now; 
                                   ! could be calculated in morrison scheme 
          cldfrac_am4(:,:,:) = 0.1   !Based on AM4
          cldamt_am4(:,:,:) = 1.e-3
          surfrain(:,:) = rcav/1000.*rho_phy(:,kts,:)/dt    ! kg/m2/s
          surfsnow(:,:) = 0.
          !we need convective precip in kg/m2/s
          rain3d_am4(:,:,:) = 0.  ! kg/m2/s, we need convective precip in kg/m2/s
          snow3d_am4(:,:,:) = 0.

          do k = kts, kte
            kp = kte-k+1
            !not used in the new scheme
            dqdt_am4(i,j,kp) = dqdti(i,k,j)! kg/kg/s, may include wv, liquid, ice
            qcv = -dqdt_am4(i,j,kp)*rho_phy(i,k,j)                        ! kg/m3/s
            if (qcv > 0.) then
              !JianHe: cldfrac is not calculated from MP
              !fraction of grid box covered by precipitating clouds.
              !following WetRemovalGOCART
              ! F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
              !may need a better way to calculate cloud fraction for wetdep
               cldfrac_am4(i,j,kp) = 0.3/(1.+0.3*1.5e-3*2.0e-3/qcv)
            else
               cldfrac_am4(i,j,kp) = 0.
            endif

            delz = abs(z8w_am4(i,j,kp+1)-z8w_am4(i,j,kp))
            drain(i,j,kp)=qcv*delz    ! kg/m2/s
          end do

          do k = 2,kte+1
            rain3d_am4(:,:,k) = rain3d_am4(:,:,k-1)+drain(i,j,k)
          enddo

          call gfdl_wet_deposition (        &
              me, master, tracer_names,&
              nt, t_am4, p_am4, p8w_am4,   &
              z_am4, z8w_am4, surfrain,  &
              surfsnow, dqdt_am4, cldamt_am4, &
              cldfrac_am4, f_snow_berg, rain3d_am4,  &
              snow3d_am4, trac_am4(:,:,:,n),  &
              wetdeptnd, 'convect', its, jts, dt)

          trac_dt (:,:,:,n) = trac_dt (:,:,:,n) + wetdeptnd(:,:,:)

        end do
     end if

        do j = jts, jte
          do i = its, ite
           !JianHe: AM4 needs vertical level from model top to bottom
            do k = kts, kte
              kp = kte-k+1
              !chem(i,k,j,1:num_chem) = trac_am4(i,j,kp,1:num_chem) - &
              !                         trac_dt (i,j,kp,1:num_chem)*dt
              do n = 1, num_chem
                nt = n+ntchs-1  ! tracer index in the tracer array
                if ((nt == ntage) .or. (nt == ntaoanh) ) then
                  chem(i,k,j,n) = trac_am4(i,j,kp,n) - trac_dt(i,j,kp,n)*dt
                else if ((nt==ntdust1) .or. (nt==ntdust2) .or. (nt==ntdust3) &
                  .or. (nt==ntdust4) .or. (nt==ntdust5) .or. (nt==ntss1) &
                  .or. (nt==ntss2) .or. (nt==ntss3) .or. (nt==ntss4) &
                  .or. (nt==ntss5) .or. (nt==ntoc1) .or. (nt==ntoc2) &
                  .or. (nt==ntbc1) .or. (nt==ntbc2) .or. (nt==ntsoa) &
                  .or. (nt==ntpp25) .or. (nt==ntpp10) .or. (nt==ntsulf)) then
                  chem(i,k,j,n) = trac_am4(i,j,kp,n) - trac_dt(i,j,kp,n)*dt
                  chem(i,k,j,n) = chem(i,k,j,n)*1.e9  ! kg/kg to ug/kg
                else
                  chem(i,k,j,n) = trac_am4(i,j,kp,n) - trac_dt(i,j,kp,n)*dt
                  chem(i,k,j,n) = chem(i,k,j,n)*1.e6  ! mol/mol to ppm
                endif
              end do

            enddo
            !JianHe: need attension on the conversion
            !ls_wetdep: kg/m2/s for aerosol and mol/m2/s for gas
            !kg/m2/s (if MMR) or mole/m2/s (if VMR)
            var_rmv(i,j,1:num_chem) = abs(ls_wetdep(i,j,1:num_chem))
            wdep(i,1) = var_rmv(i,j,p_so2)*1.e-3*64. ! kg/m2/s
            wdep(i,2) = var_rmv(i,j,p_nh3)*1.e-3*17. ! kg/m2/s
            wdep(i,3) = var_rmv(i,j,p_so4)*1.e-3*96. ! kg/m2/s
          enddo
        enddo
#endif
      case default
        ! -- no further option implemented
        errmsg = 'Logic error in catchem_wetdep_wrapper_run: invalid wdls_opt'
        errflg = 1
        return
    end select


    ! -- put chem stuff back into tracer array
    do k=kts,kte
     do i=its,ite
       gq0(i,k,ntchs:ntrac) = ppm2ugkg(1:num_chem)* max(epsilc,chem(i,k,1,1:num_chem))
       qgrs(i,k,ntchs:ntrac) = gq0(i,k,ntchs:ntrac)
     enddo
    enddo

    ! -- output large-scale wet deposition
    call gocart_diag_store(3, var_rmv, trdf)

     wetdpl (:,:)=trdf(:,1,:,3)
    end if ! wetdep_ls_opt_in 
!
   end subroutine catchem_wetdep_wrapper_run
!> @}

  subroutine catchem_prep_wetdep(ktau,dtstep,                          &
        imp_physics, imp_physics_gfdl, imp_physics_thompson,           &
        pr3d,ph3d,phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,dqdt,           &
        rri,t_phy,u_phy,v_phy,p_phy,rho_phy,z_phy,dz8w,p8w,lmk,cld_frac,cnvc,&
        t8w,dqdti,cldfrac,z_at_w,vvel,rlat,xlat,                               &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,                                       &
        ntss1,ntss2,ntss3,ntss4,ntss5,                                 &
        ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ntpp10,                &
        ntrac,ntchm,ndchm,ntchs,gq0,num_chem, num_moist,               &
        ppm2ugkg,moist,chem,gaschem_opt,                               &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    integer, intent(in) :: ktau
    real(kind=kind_phys), intent(in) :: dtstep
    integer, intent(in) :: gaschem_opt

    !FV3 input variables
    integer, intent(in) :: ntrac,ntss1,ntss2,ntss3,ntss4,ntss5
    integer, intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson
    integer, intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    integer, intent(in) :: ntchm,ndchm,ntchs
    real(kind=kind_phys), dimension(ims:ime), intent(in) ::rlat
    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         phl3d,tk3d,prl3d,us3d,vs3d,spechum,w,dqdt
    integer,        intent(in) :: lmk
    real(kind=kind_phys), dimension(ims:ime,lmk) :: cld_frac
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::  cnvc
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !GSD Chem variables
    integer,intent(in) ::  num_chem, num_moist
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(1:num_chem), intent(in) :: ppm2ugkg
    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, u_phy, v_phy, p_phy, rho_phy, z_phy, dz8w, p8w, t8w, vvel, dqdti, cldfrac
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem

    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) :: xlat
    ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
    real(kind_phys) ::  factor,factor2,pu,pl,aln,pwant
    real(kind_phys) ::  xhour,xmin,xlonn,xtime,real_time
    real(kind_phys), DIMENSION (1,1) :: sza,cosszax
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour,nt

    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    u_phy          = 0._kind_phys
    v_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    z_phy          = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    vvel           = 0._kind_phys
    dqdti          = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    xlat           = 0._kind_phys
    cldfrac        = 0._kind_phys

    do i=its,ite
    xlat (i,1)=rlat(i)*180./pi
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
          ! Geopotential height in m2 s-2 to height in m
          z_phy(i,k,j) = phl3d(ip,kp)/g
        enddo
      enddo
    enddo

    do j=jts,jte
      jp = j - jts + 1
      do k=kts,lmk
        kp = k - kts + 1
        do i=its,ite
          ip = i - its + 1
           cldfrac(i,k,j) = cld_frac(ip,kp)
           cldfrac(i,k,j) = max(0.0,cldfrac(i,k,j))
           cldfrac(i,k,j) = min(cldfrac(i,k,j),1.)
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
          u_phy(i,k,j)=us3d(ip,kkp)
          dqdti(i,k,j)=dqdt(ip,kkp)
          v_phy(i,k,j)=vs3d(ip,kkp)
          rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
          rri(i,k,j)=1./rho_phy(i,k,j)
          vvel(i,k,j)=-w(ip,kkp)*rri(i,k,j)/g 
          moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
          moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
          if (moist(i,k,j,2) < 1.e-30) moist(i,k,j,2)=0.
            if (imp_physics==imp_physics_thompson) then
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldi)
            else if (imp_physics==imp_physics_gfdl) then
            moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldii)
            endif
          if(moist(i,k,j,3) < 1.e-30) moist(i,k,j,3)=0.
          !--
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


  end subroutine catchem_prep_wetdep
!> @}
  end module catchem_wetdep_wrapper
