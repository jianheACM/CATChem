!> \file catchem_wrapper_utils.F90
!> \brief CATCHEM-CCPP interface utilities module
!>
!> \details
!> Provides wrapper utilities for interfacing CATCHEM chemistry model with CCPP
!> framework. Handles data transformation and management between host model and
!> CATCHEM chemistry calculations.
!>
!> \author Barry Baker
!>
!> \date 11/2024
!>
!> \ingroup catchem_ccpp_group
!!!>
module catchem_wrapper_utils

    use CATChem, only: ConfigType, MetStateType, ChemStateType, &
                      EmisStateType, DiagStateType, cc_read_config, &
                      cc_allocate_metstate, cc_allocate_chemstate, &
                      cc_allocate_diagstate, cc_allocate_emisstate, &
                      cc_checkdeallocate, cc_checkallocate, CC_SUCCESS
    use catchem_types, only: catchem_container_type

    implicit none
    private

    public :: read_catchem_config
    public :: allocate_catchem_container
    public :: transform_ccpp_to_catchem

contains

    !> Reads and initializes the CATChem configuration from a specified file
    !!
    !! This subroutine reads configuration settings and initializes the states
    !! for meteorology, emissions, and chemistry using the CATChem interface.
    !!
    !! \param[inout] config         Configuration type containing CATChem settings
    !! \param[inout] catchem_states Grid state container for CATChem components
    !! \param[in]    config_file    Path to the CATChem configuration file
    !! \param[out]   errflg         Error flag (0=success, non-zero=failure)
    !! \param[out]   errmsg         Error message if errflg is non-zero
    !!
    !! \ingroup catchem_ccpp_group
    !!!>
    subroutine read_catchem_config(config, catchem_states, config_file, errflg, errmsg)
        type(ConfigType), intent(inout) :: config
        type(catchem_container_type), intent(inout) :: catchem_states
        character(len=*), intent(in) :: config_file
        integer, intent(out) :: errflg
        character(len=*), intent(out) :: errmsg

        call cc_read_config(config, &
                          catchem_states%MetState(1), &
                          catchem_states%EmisState(1), &
                          catchem_states%ChemState(1), &
                          errflg, &
                          config_file)
        if (errflg /= 0) then
            errmsg = 'Error reading CATChem configuration file: '//trim(config_file)
        end if
    end subroutine read_catchem_config

    !> Allocates memory for the entire catchem container
    !!
    !! \param[in]    config         CATChem configuration
    !! \param[inout] catchem_states Container for all CATChem states
    !! \param[in]    im            Number of horizontal points
    !! \param[in]    kme           Number of vertical levels
    !! \param[in]    nsoil         Number of soil levels
    !! \param[out]   errflg        Error flag (0=success, non-zero=failure)
    !! \param[out]   errmsg        Error message if errflg is non-zero
    !!
    !! \ingroup catchem_ccpp_group
    !!!>
    subroutine allocate_catchem_container(config, catchem_states, im, kme, nsoil, nLandTypes, errflg, errmsg)
        type(ConfigType), intent(in) :: config
        type(catchem_container_type), intent(inout) :: catchem_states
        integer, intent(in) :: im
        integer, intent(in) :: kme
        integer, intent(in) :: nsoil
        integer, intent(in) :: nLandTypes
        integer, intent(out) :: errflg
        character(len=*), intent(out) :: errmsg

        integer :: i

        ! Initialize
        errflg = 0
        errmsg = ''

        ! Allocate State arrays in the container
        if (.not. allocated(catchem_states%MetState)) then
            allocate(catchem_states%MetState(im), stat=errflg)
            if (check_allocation_error('catchem_states%MetState', errflg, errmsg)) return

            ! Initialize each MetState using CC API
            do i = 1, im
                call CC_Allocate_MetState(catchem_states%MetState(i), errflg)
                if (errflg /= CC_SUCCESS) then
                    errmsg = 'Error in CC_Allocate_MetState'
                    return
                endif
            end do
        endif

        if (.not. allocated(catchem_states%ChemState)) then
            allocate(catchem_states%ChemState(im), stat=errflg)
            if (check_allocation_error('catchem_states%ChemState', errflg, errmsg)) return

            ! Initialize each ChemState using CC API
            do i = 1, im
                call CC_Allocate_ChemState(catchem_states%ChemState(i), config%Species_File, MetState(i), errflg)
                if (errflg /= CC_SUCCESS) then
                    errmsg = 'Error in CC_Allocate_ChemState'
                    return
                endif
            end do
        endif

        if (.not. allocated(catchem_states%DiagState)) then
            allocate(catchem_states%DiagState(im), stat=errflg)
            if (check_allocation_error('catchem_states%DiagState', errflg, errmsg)) return

            ! Initialize each DiagState using CC API
            do i = 1, im
                call CC_Allocate_DiagState(catchem_states%DiagState(i), config, MetState(i), errflg)
                if (errflg /= CC_SUCCESS) then
                    errmsg = 'Error in CC_Allocate_DiagState'
                    return
                endif
            end do
        endif

        if (.not. allocated(catchem_states%EmisState)) then
            allocate(catchem_states%EmisState(im), stat=errflg)
            if (check_allocation_error('catchem_states%EmisState', errflg, errmsg)) return

            ! Initialize each EmisState using CC API
            do i = 1, im
                call CC_Allocate_EmisState(catchem_states%EmisState(i), MetState(i), errflg)
                if (errflg /= CC_SUCCESS) then
                    errmsg = 'Error in CC_Allocate_EmisState'
                    return
                endif
            end do
        endif


    end subroutine allocate_catchem_container

    !> Transforms CCPP meteorological arrays into CATChem meteorological and chemistry states
    !!
    !! \param[in] im         Number of horizontal grid points
    !! \param[in] kme        Number of vertical levels
    !! \param[in] temp       3D temperature field (K)
    !! \param[in] spechum    3D specific humidity field (kg/kg)
    !! \param[in] pfull      3D full level pressure (Pa)
    !! \param[in] phalf      3D half level pressure (Pa)
    !! \param[in] u          3D zonal wind (m/s)
    !! \param[in] v          3D meridional wind (m/s)
    !! \param[in] delp       3D pressure thickness (Pa)
    !! \param[in] zh         3D geopotential height (m)
    !! \param[in] kh         3D vertical diffusivity (m2/s)
    !! \param[in] prsl       3D layer mean pressure (Pa)
    !! \param[in] prslk      3D Exner function
    !! \param[in] u10m       10m u wind (m/s)
    !! \param[in] v10m       10m v wind (m/s)
    !! \param[in] tskin      Surface skin temperature (K)
    !! \param[in] ps         Surface pressure (Pa)
    !! \param[in] precip     Precipitation rate (kg/m2/s)
    !! \param[in] slmsk      Land-sea mask (1=land,0=sea,2=ice)
    !! \param[in] snowh      Snow depth (m)
    !! \param[in] vegtype    Vegetation type
    !! \param[in] soiltyp    Soil type
    !! \param[in] hf         Sensible heat flux (W/m2)
    !! \param[in] ust        Friction velocity (m/s)
    !! \param[in] zpbl       PBL height (m)
    !! \param[in] coszen     Cosine of solar zenith angle
    !! \param[in] albedo     Surface albedo
    !! \param[in] emis       Surface emissivity
    !! \param[in] ustar      Friction velocity (m/s)
    !! \param[in] shflx      Surface sensible heat flux (W/m2)
    !! \param[in] lhflx      Surface latent heat flux (W/m2)
    !! \param[in] snowc      Snow cover fraction (0-1)
    !! \param[in] vegfrac    Green vegetation fraction (0-1)
    !! \param[in] swdn       Surface downward SW radiation (W/m2)
    !! \param[in] swup       Surface upward SW radiation (W/m2)
    !! \param[in] lwdn       Surface downward LW radiation (W/m2)
    !! \param[in] lwup       Surface upward LW radiation (W/m2)
    !! \param[in] swdnc      Clear-sky surface downward SW radiation (W/m2)
    !! \param[in] swupc      Clear-sky surface upward SW radiation (W/m2)
    !! \param[in] lwdnc      Clear-sky surface downward LW radiation (W/m2)
    !! \param[in] lwupc      Clear-sky surface upward LW radiation (W/m2)
    !! \param[inout] MetState   CATChem meteorology state
    !! \param[inout] ChemState  CATChem chemistry state
    !! \param[out] errmsg     Error message
    !! \param[out] errflg     Error flag
    !!
    !! \ingroup catchem_ccpp_group
    !!!>
    subroutine transform_ccpp_to_catchem(im, kme, kte, nsoil, nlndcat, nsoilcat, rlat, rlon, lat, lon, &      ! Grid Information
                                        ktau, dt, jdate, tile_num, garea, &  ! Grid Information
                                        aero_rad_freq_opt, aero_feedback_opt, plmrise_freq_opt, & ! Model Options
                                        lwi, dluse, &  ! Model Options
                                        temp, spechum, pfull,pfull_wet, phalf, rh, &  ! Meteorological Variables
                                        u, v, delp, zh, kh, prsl, prslk, &  ! Meteorological Variables
                                        u10m, v10m, tskin, ps, ts, precip, &  ! Meteorological Variables
                                        cldf, airden, delp_dry, &
                                        pfl_lsan, pfl_isan, &  ! precipitation variables
                                        slmsk, snowh, vegtype, soiltyp, soilmoist, &  ! Surface Variables
                                        hf, zpbl, coszen, emis, &  ! Surface Variables
                                        ustar, shflx, lhflx, &  ! Near-Surface Meteorology
                                        snowc, vegfrac, lai, frlanduse, frsoil, pores, resid, &  ! Surface Variables
                                        z0, landfrac, oceanfrac, lakefrac, landicefrac, seaicefrac, &
                                        swdn, swup, lwdn, lwup, &  ! Radiation Fluxes
                                        nirbmdi, nirdfdi, visbmdi, visdfdi, &  ! Radiation Fluxes
                                        swdnc, swupc, lwdnc, lwupc, &  ! Radiation Fluxes
                                        sfc_alb_nir_dir, sfc_alb_nir_dif, sfc_alb_uvvis_dir, sfc_alb_uvvis_dif,&  ! surface albedo
                                        MetState, ChemState, EmisState, DiagState, &  ! CATChem States
                                        errmsg, errflg)  ! Error Handling

      use CATChem, only: MetStateType, ChemStateType, DiagStateType, EmisStateType

      implicit none
      !! Transform CCPP meteorological arrays to CATChem states
      integer,  intent(in)    :: im              !> number of horizontal points
      integer,  intent(in)    :: kme             !> number of vertical interfaces
      integer,  intent(in)    :: kte             !> number of vertical levels
      !integer, intent(in)     :: nVert           !> number of vertical levels
      !integer, intent(in)     :: nVertInterface  !> number of vertical interfaces
      integer, intent(in)     :: nsoil           !> number of soil levels
      integer, intent(in)     :: nlndcat         !> number of land categories
      integer, intent(in)     :: nsoilcat        !> number of soil categories
      real(kind=phys), intent(in) :: rlat(:)     !> latitude (radian)
      real(kind=phys), intent(in) :: rlon(:)     !> longitude (radian)
      real(kind=phys), intent(in) :: lat(:)      !> latitude (degrees)
      real(kind=phys), intent(in) :: lon(:)      !> longitude (degrees)
      integer, intent(in)     :: ktau            !> number of timestep
      real(kind=phys), intent(in) :: dt          !> physics timestep (s)
      integer, intent(in)     :: jdate           !> current forecast date and time
      integer, intent(in)     :: tile_num        !> index of cubed sphere tile
      real(kind=phys), intent(in) :: garea(:)    !> grid cell area
      integer, intent(in) :: aero_rad_freq_opt   !> catchem aer radiation frequency
      integer, intent(in) :: aero_feedback_opt   !> catchem aerosol radiation feedback option
      integer, intent(in) :: plmrise_freq_opt    !> catchem plmrise frequency option
      integer, intent(in) :: lwi(:)               !> sea/land/ice=0/1/2 index
      integer, intent(in) :: dluse(:)            !> land use index
      real(kind=phys), intent(in) :: pores(30)   !> maximum soil moisture content
      real(kind=phys), intent(in) :: resid(30)   !> minimum soil moisture content

      ! 3D/Layer Variables (dim(:,:))
      real(kind=phys), intent(in) :: temp(:,:)       !> temperature (K)
      real(kind=phys), intent(in) :: spechum(:,:)    !> specific humidity (kg/kg)
      real(kind=phys), intent(in) :: pfull(:,:)      !> full level pressure of dry air (Pa)
      real(kind=phys), intent(in) :: pfull_wet(:,:)  !> full level pressure of moist air (Pa)
      real(kind=phys), intent(in) :: phalf(:,:)      !> half level pressure of dry air(Pa)
      real(kind=phys), intent(in) :: u(:,:)          !> zonal wind (m/s)
      real(kind=phys), intent(in) :: v(:,:)          !> meridional wind (m/s)
      real(kind=phys), intent(in) :: delp(:,:)       !> pressure thickness (Pa)
      real(kind=phys), intent(in) :: delp_dry(:,:)   !> dry air pressure thickness (Pa)
      real(kind=phys), intent(in) :: zh(:,:)         !> geopotential height (m)
      real(kind=phys), intent(in) :: kh(:,:)         !> vertical diffusivity (m2/s)
      real(kind=phys), intent(in) :: prsl(:,:)       !> layer mean pressure (Pa)
      real(kind=phys), intent(in) :: prslk(:,:)      !> Exner function
      real(kind=phys), intent(in) :: rh(:,:)         !> relative humidity [fraction]
      real(kind=phys), intent(in) :: pfl_lsan(:,:)   !>  liquid flux from large scale precipitation (kg/m2/s)
      real(kind=phys), intent(in) :: pfl_isan(:,:)   !>  ice flux from large scale precipitation (kg/m2/s)
        real(kind=phys), intent(in) :: airden(:,:)   !> dry air density [kg/m3]

      ! Surface Variables (dim(:))
      real(kind=phys), intent(in) :: ps(:)           !> surface pressure (Pa)
      real(kind=phys), intent(in) :: tskin(:)        !> skin temperature (K)
      real(kind=phys), intent(in) :: ts(:)           !> surface temperature (K)
      real(kind=phys), intent(in) :: slmsk(:)        !> land-sea mask (1=land,0=sea,2=ice)
      real(kind=phys), intent(in) :: snowh(:)        !> snow depth (m)
      real(kind=phys), intent(in) :: vegtype(:)      !> vegetation type
      real(kind=phys), intent(in) :: lai(:)          !> leaf area index (m2/m2)
      real(kind=phys), intent(in) :: z0(:)           !> surface roughness length (m)
      real(kind=phys), intent(in) :: soiltyp(:)      !> soil type
      real(kind=phys), intent(in) :: soilmoist(:)    !> soil moisture (m3/m3)
      real(kind=phys), intent(in) :: snowc(:)        !> snow cover fraction (0-1)
      real(kind=phys), intent(in) :: vegfrac(:)      !> green vegetation fraction (0-1)
      real(kind=phys), intent(in) :: frlanduse(:,:)  !> fraction of each land type
      real(kind=phys), intent(in) :: frsoil(:,:)     !> fraction of each soil type
      real(kind=phys), intent(in) :: landfrac(:)     !> fraction of land
      real(kind=phys), intent(in) :: oceanfrac(:)    !> fraction of ocean
      real(kind=phys), intent(in) :: lakefrac(:)     !> fraction of lake
      real(kind=phys), intent(in) :: landicefrac(:)  !> fraction of land ice
      real(kind=phys), intent(in) :: seaicefrac(:)   !> fraction of sea ice

      ! Near-Surface Meteorology (dim(:))
      real(kind=phys), intent(in) :: u10m(:)         !> 10m u wind (m/s)
      real(kind=phys), intent(in) :: v10m(:)         !> 10m v wind (m/s)
      real(kind=phys), intent(in) :: ustar(:)        !> friction velocity (m/s)
      real(kind=phys), intent(in) :: zpbl(:)         !> PBL height (m)

      ! Surface Fluxes (dim(:))
      real(kind=phys), intent(in) :: hf(:)           !> sensible heat flux (W/m2)
      real(kind=phys), intent(in) :: shflx(:)        !> surface sensible heat flux (W/m2)
      real(kind=phys), intent(in) :: lhflx(:)        !> surface latent heat flux (W/m2)
      real(kind=phys), intent(in) :: precip(:)       !> precipitation rate (kg/m2/s)

      ! Radiation Properties (dim(:))
      real(kind=phys), intent(in) :: coszen(:)            !> cosine of solar zenith angle
      real(kind=phys), intent(in) :: cldf(:)              !> fraction of grid box area in which updrafts occur
      real(kind=phys), intent(in) :: sfc_alb_nir_dir(:)   !> surface near-infrared direct albedo
      real(kind=phys), intent(in) :: sfc_alb_nir_dif(:)   !> surface near-infrared diffuse albedo
      real(kind=phys), intent(in) :: sfc_alb_uvvis_dir(:) !> surface visible + uv direct albedo
      real(kind=phys), intent(in) :: sfc_alb_uvvis_dif(:) !> surface visible + uv diffuse albedo
      real(kind=phys), intent(in) :: emis(:)              !> surface emissivity

      ! Radiation Fluxes (dim(:))
      real(kind=phys), intent(in) :: swdn(:)         !> downward shortwave radiation at surface (W/m2)
      real(kind=phys), intent(in) :: swup(:)         !> upward shortwave radiation at surface (W/m2)
      real(kind=phys), intent(in) :: lwdn(:)         !> downward longwave radiation at surface (W/m2)
      real(kind=phys), intent(in) :: lwup(:)         !> upward longwave radiation at surface (W/m2)
      real(kind=phys), intent(in) :: swdnc(:)        !> clear-sky downward shortwave radiation (W/m2)
      real(kind=phys), intent(in) :: swupc(:)        !> clear-sky upward shortwave radiation (W/m2)
      real(kind=phys), intent(in) :: lwdnc(:)        !> clear-sky downward longwave radiation (W/m2)2)
      real(kind=phys), intent(in) :: lwupc(:)        !> clear-sky upward longwave radiation (W/m2)
      real(kind=phys), intent(in) :: nirbmdi(:)      !> surface near-infrared beam shortwave radiation (W/m2)
      real(kind=phys), intent(in) :: nirdfdi(:)      !> surface near-infrared diffuse shortwave radiation (W/m2)
      real(kind=phys), intent(in) :: visbmdi(:)      !> surface visible + uv beam shortwave radiation (W/m2)
      real(kind=phys), intent(in) :: visdfdi(:)      !> surface visible + uv diffuse shortwave radiation (W/m2)


      ! Emissions
      !real(kind=phys), intent(in) :: emi_in(:)         !> emissions
      ! CATChem States
      type(MetStateType),  intent(inout) :: MetState(:)    !> CATChem meteorology state
      type(ChemStateType), intent(inout) :: ChemState(:)   !> CATChem chemistry state
      type(EmisStateType), intent(inout) :: EmisState(:)   !> CATChem emission state
      type(DiagStateType), intent(inout) :: DiagState(:)   !> CATChem diagnostic state

      ! Error handling
      character(len=*), intent(out) :: errmsg    !> error message
      integer,          intent(out) :: errflg    !> error flag

      ! Local variables
      integer :: i, k, k_rev
      real :: FRLAND_NOSNO_NOICE, FRWATER, FRICE, FRSNO
      real :: frac_temp
      ! Initialize error handling
      errmsg = ''
      errflg = 0

      ! Check dimensions
      if (size(MetState) /= im) then
        errmsg = 'MetState dimension mismatch'
        errflg = 1
        return
      endif

      ! Transform data for each horizontal point
      horiz: do i = 1, im
        ! Verify vertical dimension
        if (MetState(i)%nLEVS /= kte) then
            errmsg = 'Vertical dimension mismatch'
            errflg = 1
            return
        endif

        !Grid state variables
        MetState(i)%LAT = lat(i)
        MetState(i)%LON = lon(i)
        MetState(i)%AREA_M2 = garea(i)
        MetState(i)%TSTEP = dt
        !TODO: assume the 8 dimension of jdate is something like: (/ 2025, 5, 9, 14, 30, 0, 0, 0 /)
        MetState(i)%YMD = jdate(1) * 10000 + jdate(2) * 100 + jdate(3) !YYYYMMDD
        MetState(i)%HMS = jdate(4) * 10000 + jdate(5) * 100 + jdate(6) !HHMMSS

        ! 3D/Layer center Variables; reverse vertical index for CATChem
        !TODO: make sure layer index starts from 1 not 0
        k_rev = kte
        vert: do k = 1, kte
            MetState(i)%T(k_rev)        = temp(i,k)   !tk3d
            MetState(i)%SPHU(k_rev)     = spechum(i,k) /1000.0 !q3d (kg/kg in CC --> g/kg in MetState)
            MetState(i)%PMID_DRY(k_rev) = phalf(i,k) !prl3d_dry
            MetState(i)%U(k_rev)        = u(i,k) !us3d TODO: why use new_state in META???
            MetState(i)%V(k_rev)        = v(i,k)  !vs3d
            MetState(i)%DELP(k_rev)     = delp(i,k)
            MetState(i)%DELP_DRY(k_rev) = delp_dry(i,k)
            MetState(i)%ZMID(k_rev)     = zh(i, k) !geohlcl3d (geopotential height w.r.t local surface at layer center)
            !MetState(i)%kh(k_rev)       = kh(i,k) !TODO: MetState does not have kh defined yet.
            !MetState(i)%prsl(k_rev)    = prsl(i,k) !assume the same as phalf
            !MetState(i)%prslk(k_rev)   = prslk(i,k) !TODO: MetState does not have prslk defined yet.
            MetState(i)%AIRDEN(k_rev)  = airden(i, k)    !< Dry air density [kg/m3]
            MetState(i)%RH(k_rev)      = rh(i,k)
            MetState(i)%BXHEIGHT(k_rev)  = delp (i,k) / 9.80
            k_rev = k_rev - 1
        end do vert

        ! 3D/Layer interface Variables; reverse vertical index for CATChem
        k_rev = kme
        vert: do k = 1, kme
            MetState(i)%PEDGE_DRY(k_rev) = pfull(i,k) !pr3d_dry
            MetState(i)%PEDGE(k_rev)     = pfull_wet(i,k) !pr3d TODO: PEDGE needs to be defined in MetState
            if ( k == kme ) then
                MetState(i)%PFLLSAN(k_rev)   = 0.0 !TODO: vertical_layer_dimension not found vertical_interface_dimension
                MetState(i)%PFILSAN(k_rev)   = 0.0       ! so giving zero for the first layer for now
            else
                MetState(i)%PFLLSAN(k_rev)   = pfl_lsan(i,k)
                MetState(i)%PFILSAN(k_rev)   = pfl_isan(i,k)
            end if
            k_rev = k_rev - 1
        end do vert

        ! Surface Variables
        MetState(i)%PS       = ps(i)  !psfc
        MetState(i)%TS       = ts(i)  !ts
        MetState(i)%TSKIN    = tskin(i) !tskin
        MetState(i)%LWI      = lwi(i) ! TODO: slmsk same as lwi and slmsk can be deleted?
        MetState(i)%SNODP    = snowh(i) !snowdepth
        MetState(i)%DLUSE    = vegtype(i) !vtype
        MetState(i)%DSOILTYPE = soiltyp(i) !stype
        MetState(i)%FRVEG    = vegfrac(i) !gvf
        MetState(i)%LAI      = lai(i) !lai
        MetState(i)%FRLANDUSE(:) = frlanduse(i,:) !< Fraction of each land type
        MetState(i)%FRSOIL(:) = frsoil(i,:) !< Fraction of each soil type
        MetState(i)%FRLAND = landfrac(i)       !< Fraction of land [1]
        !https://dtcenter.ucar.edu/GMTB/v6.0.0/sci_doc/_c_c_p_psuite_nml_desp.html
        if (nlndcat == 20) then !modified NoahMP 20 categories; only choice for now (ivegsrc = 1)
            MetState%LUCNAME = 'NOAH'
        else
            errmsg = 'number of land category is not 20 and not supported'
            errflg = 1
            return
        endif
        MetState(i)%FRLAI(:) = frlanduse(i,:) * lai(i)
        MetState(i)%FRLAI(15:17) = 0.0 !manually give index 15(snow and ice), 16(barren), 17(water) zeros
        MetState(i)%FROCEAN = oceanfrac(i)         !< Fraction of ocean [1]
        MetState(i)%FRSEAICE = seaicefrac(i)  !fice
        MetState(i)%FRLANDIC = landicefrac(i)
        MetState(i)%FRLAKE = lakefrac(i)          !< Fraction of lake [1]
        MetState(i)%FRSNO = snowc(i) !snowfrac
        ! Water without sea ice
        FRWATER = MetState(i)%FRLAKE + MetState(i)%FROCEAN - MetState(i)%FRSEAICE
        ! Land and sea ice
        FRICE = MetState(i)%FRLANDIC + MetState(i)%FRSEAICE
        ! Land snow
        FRSNO = MetState(i)%FRSNO
        ! Set IsLand, IsWater, IsIce, IsSnow based on max fractional area (adoped from GEOS-Chem)
        MetState(i)%IsLand  = (FRLAND_NOSNO_NOICE > MAX(FRWATER, FRICE, FRSNO))
        MetState(i)%IsWater = (FRWATER > MAX(FRLAND_NOSNO_NOICE, FRICE, FRSNO))
        MetState(i)%IsIce   = (FRICE > MAX(FRLAND_NOSNO_NOICE, FRWATER, FRSNO))
        MetState(i)%IsSnow  = (FRSNO > MAX(FRLAND_NOSNO_NOICE, FRWATER, FRICE))
        !soil moisture
        MetState(i)%SOILM(:) = soilmoist(i, :)
        MetState(i)%GWETROOT = 0.0
        MetState(i)%GWETTOP  = 0.0
        if ( nsoilcat == 19) then ! 19 STATSGO soil type (isot = 1); only choice for now
            do k =1, nsoilcat
                if ( soiltyp(i) == k .and. (pores(k) - resid(k)) > 0.0 ) then
                    !here we asume nsoil is the lowest layer and index=1 is the top layer
                    MetState(i)%GWETROOT = (soilmoist(i,nsoil) - resid(k)) / (pores(k) - resid(k))
                    MetState(i)%GWETTOP = (soilmoist(i,1) - resid(k)) / (pores(k) - resid(k))
                    MetState(i)%GWETROOT = max(0.0, min(1.0, MetState(i)%GWETROOT))
                    MetState(i)%GWETTOP = max(0.0, min(1.0, MetState(i)%GWETTOP))
                end if
            end do
        else
            errmsg = 'number of soil category is not 19 and not supported'
            errflg = 1
            return
        end if

        !TODO: need to check the order of the five variables in dust_in
        MetState(i)%SSM = dust_in(i, 1) ! Sediment Supply Map [1]
        MetState(i)%RDRAG = dust_in(i, 2) !Drag Partition [1]
        MetState(i)%USTAR_THRESHOLD = dust_in(i, 3) !< Threshold friction velocity [m/s]
        MetState(i)%CLAYFRAC = dust_in(i, 4)        !< Fraction of clay [1]
        MetState(i)%SANDFRAC = dust_in(i, 5)       !< Fraction of sand [1]

        ! Near-Surface Meteorology
        MetState(i)%U10M     = u10m(i)
        MetState(i)%V10M     = v10m(i)
        MetState(i)%PBLH     = zpbl(i) !pblh
        MetState(i)%USTAR    = ustar(i)
        MetState(i)%CLDFRC   = cldf(i)

        ! Surface Fluxes
        MetState(i)%Z0 = z0(i) !znt
        MetState(i)%Z0H = z0(i) !< roughness height, for heat (thermal roughness) [m] TODO: give the same value as z0 for now
        MetState(i)%HFLUX    = shflx(i) !hf2d
        MetState(i)%EFLUX    = lhflx(i) !lf2d
        MetState(i)%PRECLSC  = precip(i) !rain_cplchm
        !calculate the Monin-Obukhov length following a simplified method as used in GOCART
        !in the ESMF/GOCART_GridComp/O3_GridComp/O3_GridCompMod.F90 file
        ! check MAPL/Python/MAPL/constants.py for MAPL_CP and MAPL_GRAV values
        if (abs(shflx(i)) > 1.00e-32) then
            MetState(i)%OBK  = - airden(i,kte) * 1004.68307 * temp(i,kte) * ustar(i)**3. / (0.40 * 9.80665 * shflx(i))
        else
            MetState(i)%OBK = 1.00e+05
        end if

        ! Radiation Properties
        MetState(i)%SUNCOS   = coszen(i) !xcosz
        frac_temp = visbmdi / (visbmdi + visdfdi)
        MetState(i)%ALBD_VIS   = frac_temp * sfc_alb_uvvis_dir(i) + (1.0 - frac_temp) * sfc_alb_uvvis_dif(i) !albedo_vis TODO: should be vis and uv
        frac_temp = nirbmdi / (nirbmdi + nirdfdi)
        MetState(i)%ALBD_NIR   = frac_temp * sfc_alb_nir_dir(i) + (1.0 - frac_temp) * sfc_alb_nir_dif(i) !albedo_nir
        !MetState(i)%emis     = emis(i) !emissivity have not defined in MetState yet

        ! Radiation Fluxes
        MetState(i)%PARDR = 0.47 * ( nirbmdi + visbmdi ) !TODO: empirical factor from papers
        MetState(i)%PARDF = 0.57 * ( nirdfdi + visdfdi )
        MetState(i)%SWGDN     = swdn(i) !dswsfc
        !MetState(i)%swup     = swup(i) !not definded in MetState
        !MetState(i)%lwdn     = lwdn(i)
        !MetState(i)%lwup     = lwup(i)
        !MetState(i)%swdnc    = swdnc(i)
        !MetState(i)%swupc    = swupc(i)
        !MetState(i)%lwdnc    = lwdnc(i)
        !MetState(i)%lwupc    = lwupc(i)

      end do horiz

    end subroutine transform_ccpp_to_catchem

    !> Checks for allocation errors and sets error message
    !!
    !! \param[in] state_name Name of the state being allocated
    !! \param[in] errflg     Error flag from allocation
    !! \param[out] errmsg    Output error message if allocation failed
    !! \return has_error     True if allocation error occurred
    !! \ingroup catchem_ccpp_group
    !!!>
    function check_allocation_error(state_name, errflg, errmsg) result(has_error)
        character(len=*), intent(in) :: state_name
        integer, intent(in) :: errflg
        character(len=*), intent(out) :: errmsg
        logical :: has_error

        has_error = (errflg /= 0)
        if (has_error) then
            errmsg = 'Error allocating '//trim(state_name)//' - catchem_wrapper_init'
        end if
    end function check_allocation_error

end module catchem_wrapper_utils
