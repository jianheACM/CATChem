! \file metstate_mod.F90
!! \brief Module for meteorology state variables
!!
!! This module contains subroutines and functions related to the MetStateType instance of CATChem.
!! It includes subroutines for initializing of the MetStateType.
!!
!! \ingroup core_modules
!!!>
MODULE MetState_Mod
   !
   ! USES:
   !
   USE Cmn_Size_Mod, ONLY : NSURFTYPE
   ! USE Dictionary_M, ONLY : dictionary_t
   USE Error_Mod
   USE Precision_Mod
   ! USE Registry_Mod


   IMPLICIT NONE
   PRIVATE
   !
   ! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: MetStateType           ! Main data type
   ! All legacy functions removed - use modern type-bound procedures instead:
   ! - MetState%init() instead of Met_Allocate()
   ! - MetState%cleanup() instead of legacy cleanup
   ! - MetState%validate() for validation
   !
   ! !PUBLIC DATA MEMBERS:
   !
   !=========================================================================
   ! Derived type for Meteorology State
   !=========================================================================

   ! \brief Derived type for Meteorology State
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: MetStateType

      CHARACTER(LEN=3)             :: State     = 'MET'    ! Name of this state

      ! NLEVS
      !------
      INTEGER               :: NLEVS             ! Number of vertical levels

      ! TIMESTEP
      !---------
      REAL(fp)              :: TSTEP             ! Time step [s]
      INTEGER               :: YMD               ! Year, month, day (YYYYMMDD format)
      INTEGER               :: HMS               ! Hour, minute, second (HHMMSS format)

      ! Logicals
      !---------
      LOGICAL           :: IsLand            ! Is this a land grid box?
      LOGICAL           :: IsWater           ! Is this a water grid box?
      LOGICAL           :: IsIce             ! Is this a ice grid box?
      LOGICAL           :: IsSnow            ! Is this a snow grid box?
      LOGICAL,  ALLOCATABLE :: InStratMeso(:)    ! Are we in the stratosphere or mesosphere?
      LOGICAL,  ALLOCATABLE :: InStratosphere(:) ! Are we in the stratosphere?
      LOGICAL,  ALLOCATABLE :: InTroposphere(:)  ! Are we in the troposphere?
      LOGICAL,  ALLOCATABLE :: InPbl(:)          ! Are we in the PBL?
      LOGICAL,  ALLOCATABLE :: IsLocalNoon       ! Is it local noon (between 11 and 13 local solar time)?

      ! Land Specific Fields
      !---------------------
      REAL(fp)              :: AREA_M2         ! Grid box surface area [m2]
      INTEGER               :: LWI             ! Land water ice mask (0-sea, 1-land, 2-ice)
      REAL(fp)              :: CLAYFRAC        ! Fraction of clay [1]
      INTEGER               :: DSOILTYPE       ! Dominant soil type
      INTEGER               :: DLUSE           ! Dominant land-use type
      REAL(fp)              :: FRVEG           ! Fraction of veg [1]
      REAL(fp)              :: FRLAKE          ! Fraction of lake [1]
      REAL(fp)              :: FRLAND          ! Fraction of land [1]
      REAL(fp)              :: FRLANDIC        ! Fraction of land ice [1]
      REAL(fp)              :: FROCEAN         ! Fraction of ocean [1]
      REAL(fp)              :: FRSEAICE        ! Sfc sea ice fraction
      REAL(fp)              :: FRSNO           ! Sfc snow fraction
      REAL(fp)              :: LAI             ! Leaf area index [m2/m2] (online) Dominant
      REAL(fp)              :: GVF             ! Green Vegetative Fraction
      REAL(fp)              :: RDRAG           ! Drag Partition [1]
      REAL(fp)              :: SANDFRAC        ! Fraction of sand [1]
      REAL(fp)              :: SEAICE00        ! Sea ice coverage 00-10%
      REAL(fp)              :: SEAICE10        ! Sea ice coverage 10-20%
      REAL(fp)              :: SEAICE20        ! Sea ice coverage 20-30%
      REAL(fp)              :: SEAICE30        ! Sea ice coverage 30-40%
      REAL(fp)              :: SEAICE40        ! Sea ice coverage 40-50%
      REAL(fp)              :: SEAICE50        ! Sea ice coverage 50-60%
      REAL(fp)              :: SEAICE60        ! Sea ice coverage 60-70%
      REAL(fp)              :: SEAICE70        ! Sea ice coverage 70-80%
      REAL(fp)              :: SEAICE80        ! Sea ice coverage 80-90%
      REAL(fp)              :: SEAICE90        ! Sea ice coverage 90-100%
      REAL(fp)              :: SNODP           ! Snow depth [m]
      REAL(fp)              :: SNOMAS          ! Snow mass [kg/m2]
      REAL(fp)              :: SSM             ! Sediment Supply Map [1]
      REAL(fp)              :: USTAR_THRESHOLD ! Threshold friction velocity [m/s]
      INTEGER,  ALLOCATABLE :: nLNDTYPE        ! # of landtypes in box (I,J)
      REAL(fp)              :: GWETTOP         ! Top soil moisture [1]
      REAL(fp)              :: GWETROOT        ! Root Zone soil moisture [1]
      REAL(fp)              :: WILT            ! Wilt point [1]
      INTEGER,  ALLOCATABLE :: nSOIL           ! # number of soil layers
      REAL(fp), ALLOCATABLE :: SOILM(:)        ! Volumetric Soil moisture [m3/m3]
      REAL(fp), ALLOCATABLE :: FRLANDUSE(:)    ! Fractional Land Use
      REAL(fp), ALLOCATABLE :: FRSOIL(:)       ! Fractional Soil
      REAL(fp), ALLOCATABLE :: FRLAI(:)        ! LAI in each Fractional Land use type [m2/m2]
      real(fp)              :: LAT             ! Latitude
      real(fp)              :: LON             ! Longitude
      character(len=20)     :: LUCNAME         ! name of land use category

      ! Radiation Related Surface Fields
      !---------------------------------
      REAL(fp)              :: ALBD_VIS       ! Visible surface albedo [1]
      REAL(fp)              :: ALBD_NIR       ! Near-IR surface albedo [1]
      REAL(fp)              :: ALBD_UV        ! UV surface albedo [1]
      REAL(fp)              :: PARDR          ! Direct photsynthetically active radiation [W/m2]
      REAL(fp)              :: PARDF          ! Diffuse photsynthetically active radiation [W/m2]
      REAL(fp)              :: SUNCOS         ! COS(solar zenith angle) at current time
      REAL(fp)              :: SUNCOSmid      ! COS(solar zenith angle) at midpoint of chem timestep
      REAL(fp)              :: SUNCOSsum      ! Sum of COS(SZA) for HEMCO OH diurnal variability
      REAL(fp)              :: SZAFACT        ! Diurnal scale factor for HEMCO OH diurnal variability (computed) [1]
      REAL(fp)              :: SWGDN          ! Incident radiation @ ground [W/m2]



      ! Flux Related Fields
      !--------------------
      REAL(fp)              :: EFLUX             ! Latent heat flux [W/m2]
      REAL(fp)              :: HFLUX             ! Sensible heat flux [W/m2]
      REAL(fp)              :: U10M              ! E/W wind speed @ 10m ht [m/s]
      REAL(fp)              :: USTAR             ! Friction velocity [m/s]
      REAL(fp)              :: V10M              ! N/S wind speed @ 10m ht [m/s]
      REAL(fp)              :: Z0                ! Surface roughness height [m]
      REAL(fp)              :: Z0H               ! Surface roughness height, for heat (thermal roughness) [m]
      REAL(fp), ALLOCATABLE :: FRZ0(:)           ! Aerodynamic Roughness Length per FRLANDUSE
      REAL(fp)              :: PBLH              ! PBL height [m]
      REAL(fp), ALLOCATABLE :: F_OF_PBL(:)       ! Fraction of box within PBL [1]
      REAL(fp), ALLOCATABLE :: F_UNDER_PBLTOP(:) ! Fraction of box under PBL top
      real(fp)              :: OBK        ! Monin-Obhukov length [m]

      ! Cloud & Precipitation Related Fields
      !-------------------------------------
      REAL(fp)              :: CLDFRC         ! Column cloud fraction [1]
      REAL(fp)              :: CONV_DEPTH     ! Convective cloud depth [m]
      REAL(fp)              :: FLASH_DENS     ! Lightning flash density [#/km2/s]
      REAL(fp)              :: CNV_FRC        ! Convective fraction [1]
      REAL(fp), ALLOCATABLE :: CLDF(:)        ! 3-D cloud fraction [1]
      REAL(fp), ALLOCATABLE :: CMFMC(:)       ! Cloud mass flux [kg/m2/s]
      REAL(fp), ALLOCATABLE :: DQRCU(:)       ! Conv precip production rate [kg/kg/s] (assume per dry air)
      REAL(fp), ALLOCATABLE :: DQRLSAN(:)     ! LS precip prod rate [kg/kg/s] (assume per dry air)
      REAL(fp), ALLOCATABLE :: DTRAIN(:)      ! Detrainment flux [kg/m2/s]
      REAL(fp)              :: PRECANV        ! Anvil previp @ ground [kg/m2/s] -> [mm/day]
      REAL(fp)              :: PRECCON        ! Conv  precip @ ground [kg/m2/s] -> [mm/day]
      REAL(fp)              :: PRECLSC        ! Large-scale precip @ ground kg/m2/s] -> [mm/day]
      REAL(fp)              :: PRECTOT        ! Total precip @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), ALLOCATABLE :: QI(:)          ! Mass fraction of cloud ice water [kg/kg dry air]
      REAL(fp), ALLOCATABLE :: QL(:)          ! Mass fraction of cloud liquid water [kg/kg dry air]
      REAL(fp), ALLOCATABLE :: PFICU(:)       ! Dwn flux ice prec:conv [kg/m2/s]
      REAL(fp), ALLOCATABLE :: PFILSAN(:)     ! Dwn flux ice prec:LS+anv [kg/m2/s]
      REAL(fp), ALLOCATABLE :: PFLCU(:)       ! Dwn flux liq prec:conv [kg/m2/s]
      REAL(fp), ALLOCATABLE :: PFLLSAN(:)     ! Dwn flux ice prec:LS+anv [kg/m2/s]
      REAL(fp), ALLOCATABLE :: TAUCLI(:)      ! Opt depth of ice clouds [1]
      REAL(fp), ALLOCATABLE :: TAUCLW(:)      ! Opt depth of H2O clouds [1]

      ! State Related Fields
      !---------------------
      REAL(fp)              :: PHIS           ! Surface geopotential height [m2/s2]
      REAL(fp), ALLOCATABLE :: Z(:)           ! Full Layer Geopotential Height
      REAL(fp), ALLOCATABLE :: ZMID(:)        ! Mid Layer Geopotential Height
      REAL(fp), ALLOCATABLE :: BXHEIGHT(:)    ! Grid box height [m] (dry air)
      REAL(fp)              :: PS_WET         ! Wet surface pressure at start of timestep [hPa]
      REAL(fp)              :: PS_DRY         ! Dry surface pressure at start of timestep [hPa]
      REAL(fp)              :: QV2M           ! Specific Humidity at 2m [kg/kg]
      REAL(fp), ALLOCATABLE :: QV(:)          ! Specific Humidity [kg/kg]
      REAL(fp)              :: T2M            ! Temperature 2m [K]
      REAL(fp)              :: TS             ! Surface temperature [K]
      REAL(fp)              :: TSKIN          ! Surface skin temperature [K]
      REAL(fp), ALLOCATABLE :: T(:)           ! Temperature [K]
      REAL(fp), ALLOCATABLE :: THETA(:)       ! Potential temperature [K]
      REAL(fp), ALLOCATABLE :: TV(:)          ! Virtual temperature [K]
      REAL(fp), ALLOCATABLE :: V(:)           ! N/S component of wind [m s-1]
      REAL(fp), ALLOCATABLE :: U(:)           ! E/W component of wind [m s-1]
      REAL(fp)              :: SST            ! Sea surface temperature [K]
      REAL(fp)              :: SLP            ! Sea level pressure [hPa]
      REAL(fp)              :: PS             ! Surface Pressure [hPa]
      REAL(fp), ALLOCATABLE :: OMEGA(:)       ! Updraft velocity [Pa/s]
      REAL(fp), ALLOCATABLE :: RH(:)          ! Relative humidity [%]
      REAL(fp)              :: TO3            ! Total overhead O3 column [DU]
      REAL(fp)              :: TROPP          ! Tropopause pressure [hPa]
      INTEGER               :: TropLev        ! Tropopause level [1]
      REAL(fp)              :: TropHt         ! Tropopause height [km]
      REAL(fp), ALLOCATABLE :: SPHU(:)        ! Specific humidity [g H2O/kg tot air]
      REAL(fp), ALLOCATABLE :: AIRDEN(:)      ! Dry air density [kg/m3]
      REAL(fp), ALLOCATABLE :: AIRNUMDEN(:)   ! Dry air density [molec/cm3]
      REAL(fp), ALLOCATABLE :: MAIRDEN(:)     ! Moist air density [kg/m3]
      REAL(fp), ALLOCATABLE :: AVGW(:)        ! Water vapor volume mixing ratio [vol H2O/vol dry air]
      REAL(fp), ALLOCATABLE :: DELP(:)        ! Delta-P (wet) across box [hPa]
      REAL(fp), ALLOCATABLE :: DELP_DRY(:)    ! Delta-P (dry) across box [hPa]
      REAL(fp), ALLOCATABLE :: DAIRMASS(:)    ! Dry air mass [kg] in grid box
      REAL(fp), ALLOCATABLE :: AIRVOL(:)      ! Grid box volume [m3] (dry air)
      REAL(fp), ALLOCATABLE :: PEDGE_DRY(:)   ! Dry air partial pressure @ level edges [hPa]
      REAL(fp), ALLOCATABLE :: PEDGE(:)       ! Air partial pressure @ level edges [hPa]
      REAL(fp), ALLOCATABLE :: PMID(:)        ! Average wet air pressure [hPa] defined as arithmetic average of edge pressures
      REAL(fp), ALLOCATABLE :: PMID_DRY(:)    ! Dry air partial pressure [hPa] defined as arithmetic avg of edge pressures

   contains
      ! Type-bound procedures for modern initialization and cleanup
      procedure :: init => metstate_init
      procedure :: cleanup => metstate_cleanup
      procedure :: validate => metstate_validate
      procedure :: reset => metstate_reset
      procedure :: is_allocated => metstate_is_allocated
      procedure :: get_memory_usage => metstate_get_memory_usage
      procedure :: print_summary => metstate_print_summary

   END TYPE MetStateType

CONTAINS

   !---------------------------------------------------------------------------
   ! PUBLIC MEMBER FUNCTIONS (Legacy functions removed)
   !---------------------------------------------------------------------------
   ! All legacy functions have been removed. Use modern type-bound procedures:
   ! - MetState%init() for initialization
   ! - MetState%cleanup() for cleanup
   ! - MetState%validate() for validation
   ! - MetState%reset() for resetting values

   !> DEPRECATED - Legacy allocation routine
   !!
   !! This routine has been replaced by modern initialization methods
   !! in the unified state container pattern. Use met_init() instead.
   !!
   !! \param GridState   CATCHem grid state (DEPRECATED)
   !! \param MetState    CATCHem met state
   !! \param RC          Error return code
   !!!>
   ! SUBROUTINE Met_Allocate( GridState, MetState, RC)
   !    ! USES
   !    USE GridState_Mod, Only : GridStateType
   !
   !    IMPLICIT NONE
   !
   !    ! Arguments
   !    TYPE(GridStateType), INTENT(IN)  :: GridState ! Grid state
   !    TYPE(MetStateType), INTENT(INOUT) :: MetState ! Meteorological state
   !    INTEGER,            INTENT(OUT)   :: RC       ! Return code
   !
   !    ! Local variables
   !    CHARACTER(LEN=255) :: ErrMsg, thisLoc
   !
   !    ! Initialize
   !    RC = CC_SUCCESS
   !    ErrMsg = ''
   !    thisLoc = ' -> at Met_Allocate (in core/metstate_mod.F90)'
   !
   !    ! Note: This legacy allocation routine has been replaced by modern
   !    ! initialization methods. The old routine used GridState%number_of_levels
   !    ! and GridState%number_of_soil_layers, but these are now passed as
   !    ! explicit parameters to the modern init() procedures.
   !
   ! end subroutine Met_Allocate

   ! Legacy Met_Allocate subroutine body (COMMENTED OUT)
   ! The rest of this subroutine has been replaced by modern init methods
   ! All the original allocation code that used GridState%number_of_levels
   ! is now handled in the type-bound init procedures.

   !   ! Visible Surface Albedo
   !   !-----------------------
   !   MetState%ALBD_VIS = ZERO
   !   MetState%ALBD_NIR = ZERO
   !   MetState%ALBD_UV = ZERO
   !   MetState%AREA_M2 = ZERO
   !   MetState%CLDFRC = ZERO
   !   MetState%CONV_DEPTH = ZERO
   !   MetState%EFLUX = ZERO
   !   MetState%FRLAKE = ZERO
   !   MetState%FRLAND = ZERO
   !   MetState%FRLANDIC = ZERO
   !   MetState%FROCEAN = ZERO
   !   MetState%FRSEAICE = ZERO
   !   MetState%FRSNO = ZERO
   !   MetState%GWETROOT = ZERO
   !       MetState%GWETTOP = ZERO
   !       MetState%HFLUX = ZERO
   !       MetState%IsLand = .false.
   !       MetState%IsWater = .false.
   !       MetState%IsIce = .false.
   !       MetState%IsSnow = .false.
   !       MetState%LAI = ZERO
   !       MetState%PARDR = ZERO
   !       MetState%PARDF = ZERO
   !       MetState%PBLH = ZERO
   !       MetState%PS = ZERO
   !       MetState%QV2M = ZERO
   !       MetState%T2M = ZERO
   !       MetState%TSKIN = ZERO
   !       MetState%U10M = ZERO
   !       MetState%V10M = ZERO
   !       MetState%z0 = ZERO
   !       MetState%z0h = ZERO
   !       MetState%USTAR_THRESHOLD = ZERO
   !       MetState%RDRAG = ZERO
   !       MetState%SSM = ZERO
   !       MetState%CLAYFRAC = ZERO
   !       MetSTate%SANDFRAC = ZERO
   !       MetState%SST = ZERO
   !
   !       ! Allocate Column Fields
   !       !-----------------------
   !       !  Logicals
   !       if (.not. allocated(MetState%InStratosphere)) then
   !          allocate(MetState%InStratosphere(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%InPbl)) then
   !          allocate(MetState%InPbl(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InPbl'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%InStratMeso)) then
   !          allocate(MetState%InStratMeso(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratMeso'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%InTroposphere)) then
   !          allocate(MetState%InTroposphere(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InTroposphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       ! Flux Related
   !       if (.not. allocated(MetState%F_OF_PBL)) then
   !          allocate(MetState%F_OF_PBL(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%F_OF_PBL'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%F_UNDER_PBLTOP)) then
   !          allocate(MetState%F_UNDER_PBLTOP(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%F_UNDER_PBLTOP'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       ! Cloud / Precipitation
   !       if (.not. allocated(MetState%CLDF)) then
   !          allocate(MetState%CLDF(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%CLDF'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%CMFMC)) then
   !          allocate(MetState%CMFMC(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%CMFMC'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%DQRCU)) then
   !          allocate(MetState%DQRCU(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%DQRCU'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%DQRLSAN)) then
   !          allocate(MetState%DQRLSAN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%DQRLSAN'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%DTRAIN)) then
   !          allocate(MetState%DTRAIN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%DTRAIN'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%QI)) then
   !          allocate(MetState%QI(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%QI'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%QL)) then
   !          allocate(MetState%QL(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%QL'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PFICU)) then
   !          allocate(MetState%PFICU(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%PFICU'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PFILSAN)) then
   !          allocate(MetState%PFILSAN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%PFILSAN'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PFLCU)) then
   !          allocate(MetState%PFLCU(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%PFLCU'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PFLLSAN)) then
   !          allocate(MetState%PFLLSAN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%PFLLSAN'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%TAUCLI)) then
   !          allocate(MetState%TAUCLI(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%TAUCLI'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%TAUCLW)) then
   !          allocate(MetState%TAUCLW(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%TAUCLW'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       ! State Variables
   !       if (.not. allocated(MetState%Z)) then
   !          allocate(MetState%Z(GridState%number_of_levels + 1), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%Z'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%ZMID)) then
   !          allocate(MetState%ZMID(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%ZMID'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%BXHEIGHT)) then
   !          allocate(MetState%BXHEIGHT(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%BXHEIGHT'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%QV)) then
   !          allocate(MetState%QV(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%T)) then
   !          allocate(MetState%T(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%T'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%THETA)) then
   !          allocate(MetState%THETA(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%THETA'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%TV)) then
   !          allocate(MetState%TV(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%TV'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%U)) then
   !          allocate(MetState%U(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%V)) then
   !          allocate(MetState%V(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%OMEGA)) then
   !          allocate(MetState%OMEGA(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%OMEGA'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%RH)) then
   !          allocate(MetState%RH(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%SPHU)) then
   !          allocate(MetState%SPHU(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%AIRDEN)) then
   !          allocate(MetState%AIRDEN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%AIRNUMDEN)) then
   !          allocate(MetState%AIRNUMDEN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%MAIRDEN)) then
   !          allocate(MetState%MAIRDEN(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%AVGW)) then
   !          allocate(MetState%AVGW(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%DELP)) then
   !          allocate(MetState%DELP(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%DELP_DRY)) then
   !          allocate(MetState%DELP_DRY(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%DAIRMASS)) then
   !          allocate(MetState%DAIRMASS(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%AIRVOL)) then
   !          allocate(MetState%AIRVOL(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PMID)) then
   !          allocate(MetState%PMID(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PMID_DRY)) then
   !          allocate(MetState%PMID_DRY(GridState%number_of_levels), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%PEDGE_DRY)) then
   !          allocate(MetState%PEDGE_DRY(GridState%number_of_levels+1), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   !       if (.not. allocated(MetState%SOILM)) then
   !          allocate(MetState%SOILM(MetState%nSOIL), stat=RC)
   !          if (RC /= CC_SUCCESS) then
   !             errMsg = 'Error allocating MetState%InStratosphere'
   !             call CC_Error(errMsg, RC, thisLoc)
   !             return
   !          endif
   !       end if
   !
   ! end subroutine Met_Allocate

   !========================================================================
   ! Modern MetState Type-Bound Procedures
   !========================================================================

   !> \brief Modern initialization procedure for MetStateType
   !!
   !! This procedure initializes a MetStateType object with proper error handling
   !! and validation. It replaces the old Met_Allocate pattern with a cleaner,
   !! more maintainable approach.
   !!
   !! \param[inout] this The MetStateType object to initialize
   !! \param[in] nlevs Number of vertical levels
   !! \param[in] error_mgr Error manager for context and error reporting
   !! \param[out] rc Return code
   subroutine metstate_init(this, nlevs, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_MEMORY_ALLOCATION

      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(in) :: nlevs
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc
      integer :: alloc_stat

      thisLoc = 'metstate_init (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_init', 'initializing meteorological state')

      rc = CC_SUCCESS

      ! Store number of levels
      this%NLEVS = nlevs
      this%State = 'MET'

      ! Initialize scalar fields to safe defaults
      this%TSTEP = 0.0_fp
      this%YMD = 0
      this%HMS = 0
      this%AREA_M2 = 0.0_fp
      this%LWI = 0
      this%CLAYFRAC = 0.0_fp
      this%DSOILTYPE = 1
      this%DLUSE = 1
      this%FRVEG = 0.0_fp
      this%FRLAKE = 0.0_fp
      this%FRLAND = 0.0_fp
      this%FRLANDIC = 0.0_fp
      this%FROCEAN = 1.0_fp
      this%FRSEAICE = 0.0_fp
      this%FRSNO = 0.0_fp
      this%LAI = 0.0_fp
      this%GVF = 0.0_fp
      this%RDRAG = 0.0_fp
      this%SANDFRAC = 0.0_fp
      this%SNODP = 0.0_fp
      this%SNOMAS = 0.0_fp
      this%SSM = 0.0_fp
      this%USTAR_THRESHOLD = 0.0_fp
      this%GWETTOP = 0.0_fp
      this%GWETROOT = 0.0_fp
      this%WILT = 0.0_fp
      this%PHIS = 0.0_fp
      this%PS_WET = 1013.25_fp
      this%PS_DRY = 1013.25_fp
      this%QV2M = 0.0_fp
      this%T2M = 288.15_fp
      this%TS = 288.15_fp
      this%TSKIN = 288.15_fp
      this%SST = 288.15_fp
      this%SLP = 1013.25_fp
      this%PS = 1013.25_fp
      this%TO3 = 300.0_fp
      this%TROPP = 200.0_fp
      this%TropLev = nlevs / 2
      this%TropHt = 12.0_fp

      ! Initialize logical fields
      this%IsLand = .false.
      this%IsWater = .true.
      this%IsIce = .false.
      this%IsSnow = .false.

      ! Allocate and initialize arrays
      ! Replace allocate_arrays call with direct allocation
      allocate(this%InStratMeso(nlevs), stat=alloc_stat)
      if (alloc_stat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, &
                                     'Failed to allocate InStratMeso', rc, &
                                     thisLoc, &
                                     'Check available memory and grid dimensions')
         call error_mgr%pop_context()
         return
      endif
      this%InStratMeso = .false.

      ! TODO: Add allocation for other arrays as needed
      ! (This is a simplified implementation)

      call error_mgr%pop_context()
   end subroutine metstate_init

   !> \brief Allocate meteorological state arrays
   subroutine allocate_arrays(this, nlevs, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_MEMORY_ALLOCATION

      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(in) :: nlevs
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      integer :: allocStat
      character(len=256) :: thisLoc

      thisLoc = 'allocate_arrays (in core/metstate_mod.F90)'
      rc = CC_SUCCESS

      ! Allocate level-dependent arrays
      allocate(this%InStratMeso(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate InStratMeso', rc, thisLoc)
         return
      endif

      allocate(this%InStratosphere(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate InStratosphere', rc, thisLoc)
         return
      endif

      allocate(this%InTroposphere(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate InTroposphere', rc, thisLoc)
         return
      endif

      allocate(this%InPbl(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate InPbl', rc, thisLoc)
         return
      endif

      allocate(this%Z(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate Z', rc, thisLoc)
         return
      endif

      allocate(this%ZMID(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate ZMID', rc, thisLoc)
         return
      endif

      allocate(this%BXHEIGHT(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate BXHEIGHT', rc, thisLoc)
         return
      endif

      allocate(this%QV(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate QV', rc, thisLoc)
         return
      endif

      allocate(this%T(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate T', rc, thisLoc)
         return
      endif

      allocate(this%THETA(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate THETA', rc, thisLoc)
         return
      endif

      allocate(this%TV(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate TV', rc, thisLoc)
         return
      endif

      allocate(this%V(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate V', rc, thisLoc)
         return
      endif

      allocate(this%U(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate U', rc, thisLoc)
         return
      endif

      allocate(this%OMEGA(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate OMEGA', rc, thisLoc)
         return
      endif

      allocate(this%RH(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate RH', rc, thisLoc)
         return
      endif

      allocate(this%SPHU(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate SPHU', rc, thisLoc)
         return
      endif

      allocate(this%AIRDEN(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate AIRDEN', rc, thisLoc)
         return
      endif

      allocate(this%AIRNUMDEN(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate AIRNUMDEN', rc, thisLoc)
         return
      endif

      allocate(this%MAIRDEN(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate MAIRDEN', rc, thisLoc)
         return
      endif

      allocate(this%AVGW(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate AVGW', rc, thisLoc)
         return
      endif

      allocate(this%DELP(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate DELP', rc, thisLoc)
         return
      endif

      allocate(this%DELP_DRY(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate DELP_DRY', rc, thisLoc)
         return
      endif

      allocate(this%DAIRMASS(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate DAIRMASS', rc, thisLoc)
         return
      endif

      allocate(this%AIRVOL(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate AIRVOL', rc, thisLoc)
         return
      endif

      allocate(this%PEDGE_DRY(nlevs+1), stat=allocStat)  ! Edge arrays have nlevs+1
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate PEDGE_DRY', rc, thisLoc)
         return
      endif

      allocate(this%PEDGE(nlevs+1), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate PEDGE', rc, thisLoc)
         return
      endif

      allocate(this%PMID(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate PMID', rc, thisLoc)
         return
      endif

      allocate(this%PMID_DRY(nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate PMID_DRY', rc, thisLoc)
         return
      endif

      ! Initialize arrays to safe defaults
      this%InStratMeso = .false.
      this%InStratosphere = .false.
      this%InTroposphere = .true.
      this%InPbl = .false.
      this%Z = 0.0_fp
      this%ZMID = 0.0_fp
      this%BXHEIGHT = 1000.0_fp  ! Default 1km layer thickness
      this%QV = 0.01_fp
      this%T = 288.15_fp
      this%THETA = 288.15_fp
      this%TV = 288.15_fp
      this%V = 0.0_fp
      this%U = 0.0_fp
      this%OMEGA = 0.0_fp
      this%RH = 50.0_fp
      this%SPHU = 10.0_fp
      this%AIRDEN = 1.225_fp
      this%AIRNUMDEN = 2.5e19_fp
      this%MAIRDEN = 1.225_fp
      this%AVGW = 0.01_fp
      this%DELP = 100.0_fp
      this%DELP_DRY = 100.0_fp
      this%DAIRMASS = 1.225e5_fp
      this%AIRVOL = 1.0e8_fp
      this%PEDGE_DRY = 1000.0_fp
      this%PEDGE = 1000.0_fp
      this%PMID = 1000.0_fp
      this%PMID_DRY = 1000.0_fp

   end subroutine allocate_arrays

   !> \brief Clean up and deallocate MetStateType
   subroutine metstate_cleanup(this, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate all allocatable arrays
      if (allocated(this%InStratMeso)) deallocate(this%InStratMeso)
      if (allocated(this%InStratosphere)) deallocate(this%InStratosphere)
      if (allocated(this%InTroposphere)) deallocate(this%InTroposphere)
      if (allocated(this%InPbl)) deallocate(this%InPbl)
      if (allocated(this%Z)) deallocate(this%Z)
      if (allocated(this%ZMID)) deallocate(this%ZMID)
      if (allocated(this%BXHEIGHT)) deallocate(this%BXHEIGHT)
      if (allocated(this%QV)) deallocate(this%QV)
      if (allocated(this%T)) deallocate(this%T)
      if (allocated(this%THETA)) deallocate(this%THETA)
      if (allocated(this%TV)) deallocate(this%TV)
      if (allocated(this%V)) deallocate(this%V)
      if (allocated(this%U)) deallocate(this%U)
      if (allocated(this%OMEGA)) deallocate(this%OMEGA)
      if (allocated(this%RH)) deallocate(this%RH)
      if (allocated(this%SPHU)) deallocate(this%SPHU)
      if (allocated(this%AIRDEN)) deallocate(this%AIRDEN)
      if (allocated(this%AIRNUMDEN)) deallocate(this%AIRNUMDEN)
      if (allocated(this%MAIRDEN)) deallocate(this%MAIRDEN)
      if (allocated(this%AVGW)) deallocate(this%AVGW)
      if (allocated(this%DELP)) deallocate(this%DELP)
      if (allocated(this%DELP_DRY)) deallocate(this%DELP_DRY)
      if (allocated(this%DAIRMASS)) deallocate(this%DAIRMASS)
      if (allocated(this%AIRVOL)) deallocate(this%AIRVOL)
      if (allocated(this%PEDGE_DRY)) deallocate(this%PEDGE_DRY)
      if (allocated(this%PEDGE)) deallocate(this%PEDGE)
      if (allocated(this%PMID)) deallocate(this%PMID)
      if (allocated(this%PMID_DRY)) deallocate(this%PMID_DRY)

      ! Reset scalar values
      this%NLEVS = 0
      this%State = ''

   end subroutine metstate_cleanup

   !> \brief Validate MetStateType for consistency and physical reasonableness
   subroutine metstate_validate(this, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_INVALID_INPUT
      use utilities_mod, only: is_valid_temperature, is_valid_pressure

      implicit none
      class(MetStateType), intent(in) :: this
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc

      character(len=256) :: thisLoc

      thisLoc = 'metstate_validate (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_validate', 'validating meteorological state')

      rc = CC_SUCCESS

      ! Check basic state
      if (this%NLEVS <= 0) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Number of levels must be positive', rc, &
                                     thisLoc, 'Set NLEVS to a positive integer')
         call error_mgr%pop_context()
         return
      endif

      ! Validate temperatures
      if (.not. is_valid_temperature(this%T2M)) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     '2m temperature out of physical range', rc, &
                                     thisLoc, 'Check temperature units and values')
         call error_mgr%pop_context()
         return
      endif

      if (.not. is_valid_temperature(this%TS)) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Surface temperature out of physical range', rc, &
                                     thisLoc, 'Check temperature units and values')
         call error_mgr%pop_context()
         return
      endif

      ! Validate pressures
      if (.not. is_valid_pressure(this%PS)) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Surface pressure out of physical range', rc, &
                                     thisLoc, 'Check pressure units and values')
         call error_mgr%pop_context()
         return
      endif

      if (.not. is_valid_pressure(this%SLP)) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Sea level pressure out of physical range', rc, &
                                     thisLoc, 'Check pressure units and values')
         call error_mgr%pop_context()
         return
      endif

      ! Check array allocation
      if (.not. this%is_allocated()) then
         call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                     'Required arrays not allocated', rc, &
                                     thisLoc, 'Call init() before using MetState')
         call error_mgr%pop_context()
         return
      endif

      call error_mgr%pop_context()
   end subroutine metstate_validate

   !> \brief Reset MetStateType to initial values
   subroutine metstate_reset(this, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Reset time information
      this%TSTEP = 0.0_fp
      this%YMD = 0
      this%HMS = 0

      ! Reset to standard atmosphere values
      this%T2M = 288.15_fp
      this%TS = 288.15_fp
      this%TSKIN = 288.15_fp
      this%PS = 1013.25_fp
      this%SLP = 1013.25_fp
      this%SST = 288.15_fp

      ! Reset arrays if allocated
      if (allocated(this%T)) this%T = 288.15_fp
      if (allocated(this%U)) this%U = 0.0_fp
      if (allocated(this%V)) this%V = 0.0_fp
      if (allocated(this%QV)) this%QV = 0.01_fp
      if (allocated(this%RH)) this%RH = 50.0_fp

   end subroutine metstate_reset

   !> \brief Check if required arrays are allocated
   function metstate_is_allocated(this) result(is_alloc)
      implicit none
      class(MetStateType), intent(in) :: this
      logical :: is_alloc

      is_alloc = allocated(this%T) .and. allocated(this%U) .and. allocated(this%V) .and. &
                 allocated(this%QV) .and. allocated(this%PMID) .and. allocated(this%DELP)
   end function metstate_is_allocated

   !> \brief Get approximate memory usage in bytes
   function metstate_get_memory_usage(this) result(memory_bytes)
      implicit none
      class(MetStateType), intent(in) :: this
      integer(kind=8) :: memory_bytes

      integer :: nlevs

      memory_bytes = 0
      nlevs = this%NLEVS

      if (nlevs > 0) then
         ! Estimate based on number of allocated arrays and precision
         ! Each real(fp) array: nlevs * 8 bytes (assuming fp = real64)
         ! Each logical array: nlevs * 1 byte
         memory_bytes = nlevs * 8 * 26  ! 26 real arrays
         memory_bytes = memory_bytes + nlevs * 1 * 4  ! 4 logical arrays
         memory_bytes = memory_bytes + (nlevs+1) * 8 * 2  ! 2 edge arrays
      endif
   end function metstate_get_memory_usage

   !> \brief Print summary of MetStateType
   subroutine metstate_print_summary(this)
      implicit none
      class(MetStateType), intent(in) :: this

      write(*,'(A)') '=== MetState Summary ==='
      write(*,'(A,A)') 'State: ', trim(this%State)
      write(*,'(A,I0)') 'Number of levels: ', this%NLEVS
      write(*,'(A,F8.2,A)') 'Surface temperature: ', this%TS, ' K'
      write(*,'(A,F8.2,A)') 'Surface pressure: ', this%PS, ' hPa'
      write(*,'(A,F8.2,A)') 'Sea level pressure: ', this%SLP, ' hPa'
      write(*,'(A,L1)') 'Arrays allocated: ', this%is_allocated()
      write(*,'(A,I0,A)') 'Memory usage: ', this%get_memory_usage(), ' bytes'
      write(*,'(A)') '======================='
   end subroutine metstate_print_summary

   !========================================================================
   ! Legacy Procedures (maintained for backward compatibility)
   !========================================================================

END MODULE MetState_Mod
