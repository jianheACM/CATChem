!>\file catchem_gfdlam4chem_wrapper.F90
!! This file is AM4Chem  wrapper with CCPP coupling to FV3
!! 08/2023, Include GFDL-AM4 gas-phase chemistry, Jian.He@noaa.gov

 module catchem_gfdlam4chem_wrapper

   use physcons,        only : g => con_g, pi => con_pi
   use machine ,        only : kind_phys
   use catchem_config
   use catchem_constants, only : AVOGNO,WTMAIR,epsilc
   use gocart_aerosols_mod
   use dust_data_mod
   use gfdl_am4chem_mod
   use mo_vinterp,  only : vinterp
   use gfdl_time_utls_mod, only : time_type, &
                                  get_date, &
                                  set_date, &
                                  set_time, &
                                  days_in_year, &
                                  real_to_time_type, &
                                  time_type_to_real, &
                                  operator(+), operator(-), &
                                  time_interp
   use mo_chem_utls_mod, only : sjl_fillz

   implicit none

   private

   public :: catchem_gfdlam4chem_wrapper_init, catchem_gfdlam4chem_wrapper_run, catchem_gfdlam4chem_wrapper_finalize

   logical :: is_initialized = .false.
   logical :: do_positive_adjust = .false.

   type(time_type) :: Time_init,Time,Time_next,dt_time

contains

! -----------------------------------------------------------------------
! CCPP entry points for gfdl am4 chemistry
! -----------------------------------------------------------------------

!>\brief The subroutine initializes the GFDL AM4 chemistry.
!!
!> \section arg_table_catchem_gfdlam4chem_wrapper_init Argument Table
!! \htmlinclude catchem_gfdlam4chem_wrapper_init.html
!!
      subroutine catchem_gfdlam4chem_wrapper_init(me, master, nlunit, input_nml_file, logunit, fn_nml, &
                      gaschem_opt, do_am4chem, phot_opt, idat, jdat, dt, &
                                                  tracer_names, restart, ntche, ndchm, errmsg, errflg)
       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master
       integer, intent (in) :: nlunit
       integer, intent (in) :: logunit
       integer, intent (in) :: idat(8), jdat(8)
       integer, intent (in) :: ntche, ndchm
       real(kind_phys),  intent(in)  :: dt
       character(len=*), intent(in)  :: fn_nml
       character(len=*), intent(in)  :: input_nml_file(:)
       integer,          intent(in)  :: gaschem_opt
       integer,          intent(in)  :: phot_opt
       logical,          intent(in)  :: do_am4chem
       logical,          intent(in)  :: restart 
       character(len=32), intent(in) :: tracer_names(:)
       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg
    
       ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (is_initialized) return

       if (do_am4chem) then

          if ((phot_opt .lt. 1) .or. (gaschem_opt .lt. 1)) then
          write(errmsg,'(*(a))') 'Namelist option gfdl_am4chem does not match choice in suite definition file'
          errflg = 1
          return
          end if
          if (me == master) then
            print *, 'Start GFDL AM4 gas-phase chemistry'
          endif

          ! Model initial time
          Time_init = set_date(idat(1),idat(2),idat(3), &
                        idat(5),idat(6),idat(7))
          if (me == master) then
               write (*, *) "catchem_gfdlam4chem_wrapper_init-Time_init: ",idat(1),idat(2),idat(3), &
                        idat(5),idat(6),idat(7)
          endif


          ! Model current time
          Time = set_date(jdat(1),jdat(2),jdat(3), &
                        jdat(5),jdat(6),jdat(7))

          if (me == master) then
              write (*, *) "catchem_gfdlam4chem_wrapper_init-Time_current: ",jdat(1),jdat(2),jdat(3), &
                        jdat(5),jdat(6),jdat(7)
          endif

          call gfdl_am4chem_init(me, master, nlunit, input_nml_file, logunit, fn_nml, &
            tracer_names, ntche, ndchm, Time_init,Time, restart)

       is_initialized = .true.

       end if

      end subroutine catchem_gfdlam4chem_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_gfdlam4chem_wrapper_finalize Argument Table
!!
      subroutine catchem_gfdlam4chem_wrapper_finalize(me, master,errmsg, errflg)

       implicit none

       integer, intent (in) :: me
       integer, intent (in) :: master

       character(len=*), intent(out) :: errmsg
       integer,          intent(out) :: errflg

   ! Initialize CCPP error handling variables
       errmsg = ''
       errflg = 0

       if (me==master) then
          write(*,*) 'Finishing catchem_gfdlam4chem_wrapper!'
       endif

       if (.not.is_initialized) return

       is_initialized = .false.

      end subroutine catchem_gfdlam4chem_wrapper_finalize

!> \defgroup catchem_group CATChem AM4-gocart wrapper Module
!! This is the Configurable ATmospheric Chemistry (CATChem)
!>\defgroup catchem_gfdlam4chem_wrapper CATChem AM4-gocart wrapper Module  
!> \ingroup catchem_gfdlam4chem_group
!! This is the CATChem AM4-gocart wrapper Module
!! \section arg_table_catchem_gfdlam4chem_wrapper_run Argument Table
!! \htmlinclude catchem_gfdlam4chem_wrapper_run.html
!!
!>\section catchem_gfdlam4chem_wrapper CATChem Scheme General Algorithm
!> @{
    subroutine catchem_gfdlam4chem_wrapper_run(me, master,im, kte, kme, &
                   ktau, dt, garea,   &
                   ntrac, ntchm, ndchm, ntchs, ntche, ntqv, &
                   land, landfrac,oceanfrac, fice, lmk, cld_frac, rlat, rlon, &
                   tskin, julian, xcosz,solcon,&
                   pr3d,ph3d,phl3d,prl3d, tk3d, spechum,&
                   u10m, v10m,idat, jdat, jval, emi2_in, bioem, &
                   ox_prod, ox_loss, lch4_prod, ch4_loss, &
                   oh_prod, oh_loss, &
                   dfdage_in, sfalb, &
                   ntso2, ntsulf, ntDMS, ntmsa, ntpp25,                     &
                   ntdust1,ntdust2,ntdust3,ntdust4,ntdust5, &
                   ntss1,ntss2,ntss3,ntss4,ntss5, ntsoa, ntso4, &
                   ntbc1, ntbc2, ntoc1, ntoc2, ntpp10, ntage, ntaoanh, &
                   ntextinct,ntoz, &
                   gaschem_opt, do_am4chem, tracer_names,   &               
                   chem_in_opt,chem_opt,phot_opt,                    &
                   aer_ra_frq_in,        &
                   gq0, qgrs, tile_num, errmsg, errflg) 

    implicit none

    integer, intent (in) :: me
    integer, intent (in) :: master

    integer,        intent(in) :: im,kte,kme,ktau,tile_num
    integer,        intent(in) :: idat(8), jdat(8)
    integer,        intent(in) :: ntrac, ntchm, ndchm, ntchs, ntche, ntqv
    integer,        intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer,        intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    integer,        intent(in) :: ntdust1,ntdust2,ntdust3,ntdust4,ntdust5
    integer,        intent(in) :: ntss1,ntss2,ntss3,ntss4,ntss5
    integer,        intent(in) :: ntoz, ntage, ntsoa, ntaoanh, ntextinct,ntso4
    character(len=32), dimension(ntrac), intent(in) :: tracer_names
    logical,        intent(in)  :: do_am4chem

    real(kind_phys),intent(in) :: dt,julian,solcon

    integer, parameter :: ids=1,jds=1,jde=1, kds=1
    integer, parameter :: ims=1,jms=1,jme=1, kms=1
    integer, parameter :: its=1,jts=1,jte=1, kts=1

    integer, dimension(im), intent(in) :: land
    real(kind_phys), dimension(im), intent(in) :: u10m, v10m, landfrac, oceanfrac, fice, &
                garea, rlat,rlon, tskin, sfalb, xcosz

    real(kind_phys), dimension(im,kte,3), intent(inout) :: jval  !jo2,jo1d,jno2
    real(kind_phys), dimension(im,kte), intent(inout) :: ox_prod, ox_loss, lch4_prod, ch4_loss, &
                                                         oh_prod, oh_loss
    real(kind_phys), dimension(im,64, 3), intent(in) :: emi2_in  !JianHe: vertical
    real(kind_phys), dimension(im,3), intent(in) :: bioem
    real(kind_phys), dimension(im,72, 8), intent(in) :: dfdage_in  !JianHe
    real(kind_phys), dimension(im,kme), intent(in) :: ph3d, pr3d
    real(kind_phys), dimension(im,kte), intent(in) :: prl3d,phl3d,tk3d, spechum
    integer,           intent(in) :: lmk
    real(kind_phys), dimension(im,lmk), intent(in) ::  cld_frac
    real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0, qgrs
    integer,           intent(in) :: chem_in_opt
    integer,           intent(in) :: chem_opt, phot_opt, aer_ra_frq_in
    integer,           intent(in) :: gaschem_opt
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: rri, t_phy,        &
                     p_phy, z_at_w, z_phy, dz8w, p8w, t8w, rho_phy, pwt, &
                     sphum,cldfrac

    real(kind_phys), dimension(ims:im, jms:jme) :: xlat, xlong, dxy, xland, &
                                                   frocean,fraci,tsk, &
                                                   frac_open_sea,w10m, &
                                                   albedo,icoszen, &
                                                   rlat_in,rlon_in

!>- chemistry variables
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist 
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem
    real(kind_phys), dimension(ims:im, kms:kte, jms:jme, 1:ntche) :: gq0_prog

!   JianHe: AM4 lev is in 3rd dimension (e.g., tracer(:,j,k,n)
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte,1:ntche) :: trac_prog ! conc array for all prog tracers, met+chem
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte,1:ndchm) :: chem_diag ! conc array for diag chem tracers
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte,1:ntrac) :: trac_dt   ! conc tendency array for all tracers, prog+diag
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte,3) :: jvals_out
    real(kind_phys), dimension(ims:im,jms:jme,kts:kte) :: prodox,lossox, &
                                                          prodlch4,lossch4, &
                                                          prodoh,lossoh
    real(kind_phys), dimension(1:im,jms:jme,1:kte) :: t_am4,    &
                    p_am4, z_am4, sphum_am4, cldfrac_am4
    real(kind_phys), dimension(1:im,jms:jme,1:kme) :: z8w_am4, p8w_am4

    integer :: ide, ime, ite, kde

    logical, parameter :: readrestart = .false.
    integer, parameter :: nvl_dfdage = 72 ! number of input levels from AM4 dfdage
 
    integer, parameter :: nspecies_age=8
 
    real(kind_phys), dimension(ims:im, jms:jme, kms:kte, 1:nspecies_age) :: dfdage_interp,dfdage_am4
    real(kind_phys), dimension(ims:im, kms:kme, jms:jme) :: pm10, pm2_5_dry, pm2_5_dry_ec

    real(kind_phys) :: dtstep
    real(kind_phys), dimension(1:num_chem) :: ppm2ugkg

    ! -- output tracers
    real(kind_phys), dimension(ims:im, jms:jme, 1:kme) :: p10, pm25!, ebu_oc


!>-- local variables
    logical :: call_radiation, call_am4
    integer :: i, j, jp, k, kp, n, nt, nv
    integer :: call_chem, call_phot
    real(kind_phys) :: curr_secs, esfact
    integer :: iyear, imon, iday, ihour, imin, isec
    

    errmsg = ''
    errflg = 0

if (do_am4chem) then

    call_chem=0
    call_phot=0
    curr_secs = ktau * dt

    !JianHe: we set to be 0 for now
    pm2_5_dry = 0.

! Model initial time
   Time_init = set_date(idat(1),idat(2),idat(3), &
                        idat(5),idat(6),idat(7))
   if (me == master) then
       write (*, *) "catchem_gfdlam4chem_wrapper_run-Time_init: ",im,tile_num,idat(1),idat(2),idat(3), &
                        idat(5),idat(6),idat(7)
   endif


! Model current time
   Time = set_date(jdat(1),jdat(2),jdat(3), &
                        jdat(5),jdat(6),jdat(7))

   if (me == master) then
       write (*, *) "catchem_gfdlam4chem_wrapper_run-Time_current: ",jdat(1),jdat(2),jdat(3), &
                        jdat(5),jdat(6),jdat(7)
   endif

! Model time step
   dt_time = real_to_time_type(dt)

! Model next time
   Time_next = Time + dt_time

   if (me == master) then
       call get_date( Time_next, iyear, imon, iday, ihour, imin, isec )
       write (*, *) "catchem_gfdlam4chem_wrapper_run-Time_next: ",iyear, imon, iday, ihour, imin, isec
   endif

   jvals_out = 0._kind_phys
   prodox = 0._kind_phys
   lossox = 0._kind_phys
   prodlch4 = 0._kind_phys
   lossch4 = 0._kind_phys
   prodoh = 0._kind_phys
   lossoh = 0._kind_phys

   trac_dt = 0._kind_phys 
    ! -- set domain
    ide=im 
    ime=im
    ite=im
    kde=kte

    ! -- volume to mass fraction conversion table (ppm -> ug/kg)
    ppm2ugkg         = 1._kind_phys
   !ppm2ugkg(p_so2 ) = 1.e+03_kind_phys * mw_so2_aer / mwdry
    !chem(p_sulf) is in ppm, gq(ntsulf) is in ug/kg
    ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

    ! -- set control flags
    call_radiation   = (mod(int(curr_secs), max(1, 60*aer_ra_frq_in)) == 0) .or. (ktau == 1)
    call_am4         = (mod(ktau, call_chemistry) == 0) .or. (ktau == 1)

    if ((mod(ktau,call_chemistry)==0).or.(ktau==1)) then
       call_chem=1
    end if

    if (call_radiation .and. phot_opt >= 1) then
       call_phot=1
    end if

    if (ktau > 1) then
      dtstep = call_chemistry * dt
    else
      dtstep = dt
    end if

!!!

!>- get ready for chemistry run
    call catchem_gfdlam4chem_prep(                                             &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                      &
        garea,land,rlat,rlon,u10m,v10m,                    &
        landfrac,oceanfrac, fice, tskin, sfalb, &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,lmk,cld_frac,                 &
        sphum,cldfrac,emi2_in,bioem,dfdage_in,                               &
        xlat,xlong,dxy,xland,frac_open_sea,tsk,w10m,                    &
        rlat_in,rlon_in, &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,                  &
        t8w,z_at_w, z_phy, pwt,albedo,icoszen,                  &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                                &
        ntbc1,ntbc2,ntoc1,ntoc2,ntpp10,ntoz,                      &
        chem_opt, gaschem_opt, tracer_names, ntchm, ntchs,  &
        ntche,ndchm,ntrac,gq0,                                 &
        num_chem, num_moist,                                 &
        call_am4,nvl_dfdage,                    &
        dfdage_interp,ppm2ugkg,                                         &
        moist,chem, gq0_prog,                               &
        ids,ide, jds,jde, kds,kde,                                      &
        ims,ime, jms,jme, kms,kme,                                      &
        its,ite, jts,jte, kts,kte)

    if (call_chem == 1) then
      !JianHe: trop+strat chem
      !We need rearange tracer array for AM4
        do j = jts, jte
          do i = its, ite
           !JianHe: AM4 needs vertical level from model top to bottom
            do k = kts, kte
              kp = kte-k+1
              do nt = 1, ntche
              ! We need to make sure the order of tracers
              ! consistent with field table
              ! trac_prog array stores prog met+chem tracers only
              ! gas vmr, aerosol mmr
                if ( (nt < ntchs) .or. (nt == ntage) .or. (nt == ntaoanh) ) then
                  trac_prog(i,j,kp,nt) = gq0_prog(i,k,j,nt)
                else if ((nt==ntdust1) .or. (nt==ntdust2) .or. (nt==ntdust3) &
                  .or. (nt==ntdust4) .or. (nt==ntdust5) .or. (nt==ntss1) &
                  .or. (nt==ntss2) .or. (nt==ntss3) .or. (nt==ntss4) &
                  .or. (nt==ntss5) .or. (nt==ntoc1) .or. (nt==ntoc2) &
                  .or. (nt==ntbc1) .or. (nt==ntbc2) .or. (nt==ntsoa) &
                  .or. (nt==ntpp25) .or. (nt==ntpp10) .or. (nt==ntsulf)) then 
                  trac_prog(i,j,kp,nt) = gq0_prog(i,k,j,nt)*1.e-9  ! ug/kg to kg/kg
                else
                  trac_prog(i,j,kp,nt) = gq0_prog(i,k,j,nt)*1.e-6  ! ppm to mol/mol
                endif
              end do
              do nt = 1, ndchm
              ! start with so2, we need to make sure the order of tracers
              ! consistent with field table
              ! chem_diag array stores diag chem tracers only, all are gases for
              ! now
                nv = ntche+nt   ! index in tracer array
                if (nv==ntextinct) then
                  ! In AM4, this is SW extinction (band 4) for vocanoes 
                  ! We currently read AM4 output for this tracer
                  ! could be linked to diag extinction in the future
                  chem_diag(i,j,kp,nt) = chem(i,k,j,ntchm+nt)    ! /m
                else
                  chem_diag(i,j,kp,nt) = chem(i,k,j,ntchm+nt)*1.e-6  ! ppm to mol/mol
                endif
              end do
      
              t_am4(i,j,kp) = t_phy(i,k,j)  !K
              p_am4(i,j,kp) = p_phy(i,k,j)  !Pa
              z_am4(i,j,kp) = z_phy(i,k,j)
              sphum_am4(i,j,kp) = sphum(i,k,j)
              cldfrac_am4(i,j,kp) = cldfrac(i,k,j)
              dfdage_am4(i,j,kp,1:nspecies_age) = dfdage_interp(i,j,k,1:nspecies_age)
            end do
            do k = kts, kte+1
              kp = kte+1-k+1
              z8w_am4(i,j,kp) = z_at_w(i,k,j)
              p8w_am4(i,j,kp) = p8w(i,k,j)
            end do
          end do
        end do

        call gfdl_am4chem_driver( me,master,Time_init,Time,Time_next,dt_time,&
             tracer_names, rlon_in, rlat_in, xland,  &
             frac_open_sea, trac_prog, trac_dt, p8w_am4,p_am4, &
             t_am4, its, ite, jts, jte,  &
             dtstep, z8w_am4, z_am4, sphum_am4, &
             cldfrac_am4,dfdage_am4, &
             tsk,albedo,icoszen,solcon, &
             dxy, w10m, chem_diag, jvals_out, &
             prodox,lossox,prodlch4,lossch4,prodoh,lossoh)


        do nt = ntchs, ntche
          trac_prog(its:ite,jts:jte,kts:kte,nt) = trac_prog(its:ite,jts:jte,kts:kte,nt) + &
                                   trac_dt(its:ite,jts:jte,kts:kte,nt)*dtstep 
        enddo                 

        do j = jts, jte
          do i = its, ite
            do k = kts, kte
              kp = kte-k+1 
              ! prognostic tracer with tendency
              ! JianHe: we do not adjust sphum for now
              do nt = ntchs, ntche
                if ((nt == ntage) .or. (nt == ntaoanh)) then
                    gq0(i,k,nt) = trac_prog(i,j,kp,nt)
                else if ((nt==ntdust1) .or. (nt==ntdust2) .or. (nt==ntdust3) &
                  .or. (nt==ntdust4) .or. (nt==ntdust5) .or. (nt==ntss1) &
                  .or. (nt==ntss2) .or. (nt==ntss3) .or. (nt==ntss4) &
                  .or. (nt==ntss5) .or. (nt==ntoc1) .or. (nt==ntoc2) &
                  .or. (nt==ntbc1) .or. (nt==ntbc2) .or. (nt==ntsoa) &
                  .or. (nt==ntpp25) .or. (nt==ntpp10) .or. (nt==ntsulf)) then
                  gq0(i,k,nt) = trac_prog(i,j,kp,nt)*1.e9  ! reverse, kg/kg to ug/kg
                  gq0(i,k,nt) =  max(epsilc,gq0(i,k,nt))
                else
                  gq0(i,k,nt) = trac_prog(i,j,kp,nt)*1.e6  ! reverse, vmr to ppm
                  gq0(i,k,nt) =  max(epsilc,gq0(i,k,nt))
                endif
              end do
              ! diagnostic tracer
              do nt = 1,ndchm
                nv = ntche+nt
                if (nv==ntextinct) then
                  gq0(i,k,nv) = chem_diag(i,j,kp,nt)
                else
                  gq0(i,k,nv) = chem_diag(i,j,kp,nt)*1.e6 !reverse, vmr to ppm
                endif
                gq0(i,k,nv) = max(epsilc,gq0(i,k,nv))
              end do

              ! update chem array
              chem(i,k,j,1:num_chem) = gq0(i,k,ntchs:ntrac)/ppm2ugkg(1:num_chem)

              jval(i,k,:) = jvals_out(i,j,kp,:)
              ox_prod(i,k) = prodox(i,j,kp)
              ox_loss(i,k) = lossox(i,j,kp)
              lch4_prod(i,k) = prodlch4(i,j,kp)
              ch4_loss(i,k) = lossch4(i,j,kp)
              oh_prod(i,k) = prodoh(i,j,kp)
              oh_loss(i,k) = lossoh(i,j,kp)

            end do
          end do
        end do

      call gocart_aerosols_driver(ktau,dtstep,t_phy,moist,              &
           chem,rho_phy,dz8w,p8w,dxy,g,                              &
           chem_opt,num_chem,num_moist,                                 &
           ids,ide, jds,jde, kds,kde,                                   &
           ims,ime, jms,jme, kms,kme,                                   &
           its,ite, jts,jte, kts,kte                        )
    endif  ! call_chem

    call sum_pm_gocart (                                                &
         rri, chem,pm2_5_dry, pm2_5_dry_ec, pm10,                       &
         num_chem,chem_opt,                                             &
         ids,ide, jds,jde, kds,kde,                                     &
         ims,ime, jms,jme, kms,kme,                                     &
         its,ite, jts,jte, kts,kte)

    ! -- pm25 and pm10 for output , not for tracer options
    do j = jts, jte
      do k = kts, kte
        do i = its, ite
          pm25  (i,j,k) = pm2_5_dry(i,k,j)
          p10   (i,j,k) = pm10     (i,k,j)
        end do
      end do
    end do

    ! -- put chem stuff back into tracer array
    do k=kts,kte
     kp = kte-k+1
     do i=its,ite
        gq0(i,k,ntchs:ntrac) = ppm2ugkg(1:num_chem)* max(epsilc,chem(i,k,1,1:num_chem))
     enddo
    enddo

    do k=kts,kte
     do i=its,ite
       !JianHe: 08/2023: can we use some kind of index loop for this?
        qgrs(i,k,1:ntrac)=gq0(i,k,1:ntrac)
     enddo
    enddo

end if  ! do_am4chem
!
   end subroutine catchem_gfdlam4chem_wrapper_run
!> @}

   subroutine catchem_gfdlam4chem_prep(                                        &
        readrestart,chem_in_opt,ktau,dtstep,xcosz,                 &
        garea,land,rlat,rlon,u10m,v10m,                   &
        landfrac, oceanfrac, fice, ts2d, sfalb,    &
        pr3d,ph3d,phl3d,tk3d,prl3d,spechum,lmk,cld_frac,                &
        sphum,cldfrac,emi2_in,bioem,dfdage_in,                               &
        xlat,xlong,dxy,xland,frac_open_sea,tsk,w10m,           &
        rlat_in, rlon_in, &
        rri,t_phy,p_phy,rho_phy,dz8w,p8w,                 &
        t8w,z_at_w,z_phy,pwt,albedo,icoszen,                  &
        ntso2,ntsulf,ntDMS,ntmsa,ntpp25,                               &
        ntbc1,ntbc2,ntoc1,ntoc2,ntpp10,ntoz,                             &
        chem_opt, gaschem_opt, tracer_names, ntchm, ntchs,  &
        ntche,ndchm, ntrac,gq0,                              &
        num_chem, num_moist,                                &
        call_am4, nvl_dfdage,               &
        dfdage_interp,ppm2ugkg,                                        &
        moist,chem,gq0_prog,                              &
        ids,ide, jds,jde, kds,kde,                                     &
        ims,ime, jms,jme, kms,kme,                                     &
        its,ite, jts,jte, kts,kte)

    !Chem input configuration
    logical, intent(in) :: readrestart
    integer, intent(in) :: chem_in_opt, ktau
    real(kind=kind_phys), intent(in) :: dtstep

    !FV3 input variables
    integer, dimension(ims:ime), intent(in) :: land
    integer, intent(in) :: ntrac,ntchm,ntchs,ntche,ndchm
    integer, intent(in) :: ntso2,ntpp25,ntbc1,ntoc1,ntpp10
    integer, intent(in) :: ntsulf,ntbc2,ntoc2,ntDMS,ntmsa
    integer, intent(in) :: ntoz
    integer, intent(in) :: chem_opt, gaschem_opt,lmk
    character(len=32), dimension(ntrac), intent(in) :: tracer_names

    real(kind=kind_phys), dimension(ims:ime), intent(in) :: landfrac, oceanfrac,fice, &
         u10m, v10m, garea, rlat, rlon, ts2d, xcosz, sfalb
    real(kind=kind_phys), dimension(ims:ime, 64, 3),   intent(in) :: emi2_in
    real(kind_phys), dimension(ims:ime,3), intent(in) :: bioem
    real(kind=kind_phys), dimension(ims:ime, 72, 8),   intent(in) :: dfdage_in

    real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d,ph3d
    real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) ::       &
         tk3d,prl3d,spechum,phl3d
    real(kind=kind_phys), dimension(ims:ime, lmk), intent(in) :: cld_frac
    real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0


    !Chem variables
    integer,intent(in) ::  num_chem, num_moist, nvl_dfdage
    logical,intent(in) ::  call_am4
    integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
                           ims,ime, jms,jme, kms,kme,                      &
                           its,ite, jts,jte, kts,kte

    real(kind_phys), dimension(1:num_chem), intent(in) :: ppm2ugkg
    
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::              & 
         rri, t_phy, p_phy, rho_phy, dz8w, p8w, t8w, sphum, cldfrac
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) ::  z_at_w, pwt, & 
                                                                           z_phy

    real(kind_phys), dimension(ims:ime, jms:jme),          intent(out) ::              &
         xlat, xlong, dxy, xland, frac_open_sea, w10m,  tsk, &
         albedo, icoszen, rlat_in, rlon_in
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_moist), intent(out) :: moist
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, 1:num_chem),  intent(out) :: chem
    real(kind_phys), dimension(ims:ime, kms:kte, jms:jme, 1:ntche), intent(out) :: gq0_prog
    real(kind_phys), dimension(ims:ime, jms:jme, kms:kte, 8), intent(out) :: dfdage_interp

    real(kind_phys), dimension(nvl_dfdage) :: p_dfdage

    ! -- local variables
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: delp
    real(kind_phys), dimension(ims:ime) :: chem_adjust

    real(kind_phys) ::  pu,pl,aln,pwant
    real(kind_phys), DIMENSION (1,1) :: sza,cosszax
    integer i,ip,j,jp,k,kp,kk,kkp,nv,jmax,jmaxi,l,ll,n,ndystep,ixhour,nt

    integer, parameter :: nspecies_age=8
    integer, parameter :: nspecies_ic=25
    integer, parameter :: nspecies_bvoc=2
    real(kind_phys), dimension(ims:ime, jms:jme, 1:nvl_dfdage, 1:nspecies_age) :: dfdage
    real(kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_mb

    p_dfdage = (/ 996.1072, 987.3705, 976.5037, 963.1111, 946.795, 927.1469, 903.7763, &
    876.2767, 844.2427, 807.2709, 764.9445, 719.5322, 673.5149, 627.1146, &
    580.5856, 534.2095, 488.279, 443.0946, 398.9813, 356.2639, 315.2542, &
    276.2535, 239.5307, 205.2897, 173.6734, 144.9057, 119.4669, 97.04092, &
    77.65289, 61.8385, 49.24479, 39.21585, 31.22935, 24.86934, 19.80458, &
    15.77128, 12.55939, 10.00161, 7.964731, 6.342675, 5.050959, 4.022306, &
    3.203144, 2.550808, 2.031324, 1.617635, 1.288196, 1.025849, 0.8169296, &
    0.6505579, 0.5180687, 0.4125615, 0.3285413, 0.2616323, 0.2083496, &
    0.1659182, 0.1321282, 0.1052196, 0.08379117, 0.06672669, 0.05313747, & 
    0.04231574, 0.03369796, 0.02683523, 0.02137008, 0.01701797, 0.01355218, &
    0.0107922, 0.008594343, 0.006844047, 0.004711564, 0.002703662 /)


    ! -- initialize output arrays
    rri            = 0._kind_phys
    t_phy          = 0._kind_phys
    p_phy          = 0._kind_phys
    rho_phy        = 0._kind_phys
    dz8w           = 0._kind_phys
    p8w            = 0._kind_phys
    t8w            = 0._kind_phys
    sphum          = 0._kind_phys
    cldfrac        = 0._kind_phys
    xland          = 0._kind_phys
    xlat           = 0._kind_phys
    xlong          = 0._kind_phys
    dxy            = 0._kind_phys
    moist          = 0._kind_phys  
    chem           = 0._kind_phys
    z_at_w         = 0._kind_phys
    z_phy          = 0._kind_phys
    tsk            = 0._kind_phys ! Surface skin temperature (K) 
    w10m           = 0._kind_phys ! 10-m wind speed (m/s)
    frac_open_sea  = 0._kind_phys ! ocean fraction exclu. ice
    albedo         = 0._kind_phys ! Surface albedo
    icoszen        = 0._kind_phys 
    dfdage_interp  = 0._kind_phys 
    gq0_prog       = 0._kind_phys
    rlat_in        = 0._kind_phys
    rlon_in        = 0._kind_phys

    pwt  = 0._kind_phys  ! Pressure weighting (air mass) for each layer (kg/m2)
    delp           = 0._kind_phys

    do i=its,ite
     xland(i,1)=landfrac(i)
     dxy  (i,1)=garea(i)
     xlat (i,1)=rlat(i)*180./pi
     xlong(i,1)=rlon(i)*180./pi
     rlat_in(i,1)=rlat(i)
     rlon_in(i,1)=rlon(i)
     tsk  (i,1)=ts2d (i)
     w10m (i,1)=sqrt(u10m(i)*u10m(i)+v10m(i)*v10m(i))
     frac_open_sea(i,1) = min(max(0., oceanfrac(i)*(1.-fice(i))),oceanfrac(i))
     albedo(i,1)=sfalb(i)
     icoszen(i,1)=xcosz(i)
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
          sphum(i,k,j) = spechum(ip,kp)
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

    p_mb = 0.01*p_phy ! Pa to mb
   
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
         do nt = 1,ntchs-1  ! prognostic metvar
           gq0_prog(i,k,jts,nt) = gq0(i,k,nt)
         end do
     enddo
    enddo

     do i=its,ite
       do j=jts,jte
         do k=kms,kte
           do nt = ntchs,ntche   ! prognostic chemvar
             gq0_prog(i,k,j,nt) = max(epsilc,chem(i,k,j,nt-ntchs+1)*ppm2ugkg(nt-ntchs+1))
           end do
         enddo
       enddo
     enddo

     do i=its,ite
       do k=1,nvl_dfdage
          dfdage(i,1,k,1:nspecies_age)=dfdage_in(i,k,1:nspecies_age)
       enddo
     enddo


    !
    if (call_am4 .and. (gaschem_opt == 1)) then
      do j=jts,jte
        do i=its,ite
          do k=kts,kte
              call vinterp(dfdage(i,j,1:nvl_dfdage,1:nspecies_age),nspecies_age,nvl_dfdage,p_dfdage,p_mb(i,k,j),dfdage_interp(i,j,k,1:nspecies_age))
          enddo
        enddo
      enddo
    endif  

  end subroutine catchem_gfdlam4chem_prep

!> @}
  end module catchem_gfdlam4chem_wrapper

