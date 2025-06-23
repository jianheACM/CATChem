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
   USE Error_Mod
   USE Precision_Mod
   USE GridGeometry_Mod, only: GridGeometryType
   USE TimeState_Mod, only: TimeStateType



   IMPLICIT NONE
   PRIVATE

   PUBLIC :: MetStateType           ! Main data type

   !=========================================================================
   ! Derived type for Meteorology State
   !=========================================================================

   ! \brief Derived type for Meteorology State
   !!
   !! Contains all meteorological state variables for CATChem including
   !! land, radiation, flux, cloud, and state-related fields. Use type-bound
   !! procedures for initialization, cleanup, validation, and memory usage.
   !!
   !! \ingroup core_modules
   !!!>
   TYPE, PUBLIC :: MetStateType
      CHARACTER(LEN=3)             :: State     = 'MET'    !< Name of this state
      INTEGER                      :: NLEVS             !< Number of vertical levels
      ! Grid flags (2D: nx, ny)
      LOGICAL, ALLOCATABLE         :: IsLand(:,:)       !< Is this a land grid box?
      LOGICAL, ALLOCATABLE         :: IsWater(:,:)      !< Is this a water grid box?
      LOGICAL, ALLOCATABLE         :: IsIce(:,:)        !< Is this an ice grid box?
      LOGICAL, ALLOCATABLE         :: IsSnow(:,:)       !< Is this a snow grid box?
      ! Vertical flags and arrays (3D: nx, ny, nz)
      LOGICAL,  ALLOCATABLE        :: InStratMeso(:,:,:)    !< Are we in the stratosphere or mesosphere?
      LOGICAL,  ALLOCATABLE        :: InStratosphere(:,:,:) !< Are we in the stratosphere?
      LOGICAL,  ALLOCATABLE        :: InTroposphere(:,:,:)  !< Are we in the troposphere?
      LOGICAL,  ALLOCATABLE        :: InPbl(:,:,:)          !< Are we in the PBL?
      LOGICAL,  ALLOCATABLE        :: IsLocalNoon(:,:)      !< Is it local noon (between 11 and 13 local solar time)?
      ! Surface properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: AREA_M2(:,:)      !< Grid box surface area [m2]
      INTEGER,  ALLOCATABLE        :: LWI(:,:)          !< Land water ice mask (0-sea, 1-land, 2-ice)
      REAL(fp), ALLOCATABLE        :: CLAYFRAC(:,:)     !< Fraction of clay [1]
      INTEGER,  ALLOCATABLE        :: DSOILTYPE(:,:)    !< Dominant soil type
      INTEGER,  ALLOCATABLE        :: DLUSE(:,:)        !< Dominant land-use type
      REAL(fp), ALLOCATABLE        :: FRVEG(:,:)        !< Fraction of veg [1]
      REAL(fp), ALLOCATABLE        :: FRLAKE(:,:)       !< Fraction of lake [1]
      REAL(fp), ALLOCATABLE        :: FRLAND(:,:)       !< Fraction of land [1]
      REAL(fp), ALLOCATABLE        :: FRLANDIC(:,:)     !< Fraction of land ice [1]
      REAL(fp), ALLOCATABLE        :: FROCEAN(:,:)      !< Fraction of ocean [1]
      REAL(fp), ALLOCATABLE        :: FRSEAICE(:,:)     !< Sfc sea ice fraction
      REAL(fp), ALLOCATABLE        :: FRSNO(:,:)        !< Sfc snow fraction
      REAL(fp), ALLOCATABLE        :: LAI(:,:)          !< Leaf area index [m2/m2] (online) Dominant
      REAL(fp), ALLOCATABLE        :: GVF(:,:)          !< Green Vegetative Fraction
      REAL(fp), ALLOCATABLE        :: RDRAG(:,:)        !< Drag Partition [1]
      REAL(fp), ALLOCATABLE        :: SANDFRAC(:,:)     !< Fraction of sand [1]
      ! Surface and ice properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: SEAICE00(:,:)     !< Sea ice coverage 00-10%
      REAL(fp), ALLOCATABLE        :: SEAICE10(:,:)     !< Sea ice coverage 10-20%
      REAL(fp), ALLOCATABLE        :: SEAICE20(:,:)     !< Sea ice coverage 20-30%
      REAL(fp), ALLOCATABLE        :: SEAICE30(:,:)     !< Sea ice coverage 30-40%
      REAL(fp), ALLOCATABLE        :: SEAICE40(:,:)     !< Sea ice coverage 40-50%
      REAL(fp), ALLOCATABLE        :: SEAICE50(:,:)     !< Sea ice coverage 50-60%
      REAL(fp), ALLOCATABLE        :: SEAICE60(:,:)     !< Sea ice coverage 60-70%
      REAL(fp), ALLOCATABLE        :: SEAICE70(:,:)     !< Sea ice coverage 70-80%
      REAL(fp), ALLOCATABLE        :: SEAICE80(:,:)     !< Sea ice coverage 80-90%
      REAL(fp), ALLOCATABLE        :: SEAICE90(:,:)     !< Sea ice coverage 90-100%
      REAL(fp), ALLOCATABLE        :: SNODP(:,:)        !< Snow depth [m]
      REAL(fp), ALLOCATABLE        :: SNOMAS(:,:)       !< Snow mass [kg/m2]
      REAL(fp), ALLOCATABLE        :: SSM(:,:)          !< Sediment Supply Map [1]
      REAL(fp), ALLOCATABLE        :: USTAR_THRESHOLD(:,:) !< Threshold friction velocity [m/s]
      ! Soil and land use arrays (2D for counts, 3D for fractions)
      INTEGER,  ALLOCATABLE        :: nLNDTYPE(:,:)     !< # of landtypes in box (I,J)
      REAL(fp), ALLOCATABLE        :: GWETTOP(:,:)      !< Top soil moisture [1]
      REAL(fp), ALLOCATABLE        :: GWETROOT(:,:)     !< Root Zone soil moisture [1]
      REAL(fp), ALLOCATABLE        :: WILT(:,:)         !< Wilt point [1]
      INTEGER,  ALLOCATABLE        :: nSOIL(:,:)        !< # number of soil layers
      REAL(fp), ALLOCATABLE        :: SOILM(:,:,:)      !< Volumetric Soil moisture [m3/m3] (nx,ny,nsoil)
      REAL(fp), ALLOCATABLE        :: FRLANDUSE(:,:,:)  !< Fractional Land Use (nx,ny,nlanduse)
      REAL(fp), ALLOCATABLE        :: FRSOIL(:,:,:)     !< Fractional Soil (nx,ny,nsoil)
      REAL(fp), ALLOCATABLE        :: FRLAI(:,:,:)      !< LAI in each Fractional Land use type [m2/m2] (nx,ny,nlanduse)
      ! Location arrays (1D for single point, 2D for grid)
      real(fp), ALLOCATABLE        :: LAT(:,:)         !< Latitude
      real(fp), ALLOCATABLE        :: LON(:,:)         !< Longitude
      character(len=20), ALLOCATABLE :: LUCNAME(:,:)   !< name of land use category
      ! Surface meteorological properties (2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: ALBD_VIS(:,:)     !< Visible surface albedo [1]
      REAL(fp), ALLOCATABLE        :: ALBD_NIR(:,:)     !< Near-IR surface albedo [1]
      REAL(fp), ALLOCATABLE        :: ALBD_UV(:,:)      !< UV surface albedo [1]
      REAL(fp), ALLOCATABLE        :: PARDR(:,:)        !< Direct photsynthetically active radiation [W/m2]
      REAL(fp), ALLOCATABLE        :: PARDF(:,:)        !< Diffuse photsynthetically active radiation [W/m2]
      REAL(fp), ALLOCATABLE        :: SUNCOS(:,:)       !< COS(solar zenith angle) at current time
      REAL(fp), ALLOCATABLE        :: SUNCOSmid(:,:)    !< COS(solar zenith angle) at midpoint of chem timestep
      REAL(fp), ALLOCATABLE        :: SUNCOSsum(:,:)    !< Sum of COS(SZA) for HEMCO OH diurnal variability
      REAL(fp), ALLOCATABLE        :: SZAFACT(:,:)      !< Diurnal scale factor for HEMCO OH diurnal variability (computed) [1]
      REAL(fp), ALLOCATABLE        :: SWGDN(:,:)        !< Incident radiation @ ground [W/m2]
      REAL(fp), ALLOCATABLE        :: EFLUX(:,:)        !< Latent heat flux [W/m2]
      REAL(fp), ALLOCATABLE        :: HFLUX(:,:)        !< Sensible heat flux [W/m2]
      REAL(fp), ALLOCATABLE        :: U10M(:,:)         !< E/W wind speed @ 10m ht [m/s]
      REAL(fp), ALLOCATABLE        :: USTAR(:,:)        !< Friction velocity [m/s]
      REAL(fp), ALLOCATABLE        :: V10M(:,:)         !< N/S wind speed @ 10m ht [m/s]
      REAL(fp), ALLOCATABLE        :: Z0(:,:)           !< Surface roughness height [m]
      REAL(fp), ALLOCATABLE        :: Z0H(:,:)          !< Surface roughness height, for heat (thermal roughness) [m]
      REAL(fp), ALLOCATABLE        :: FRZ0(:,:,:)       !< Aerodynamic Roughness Length per FRLANDUSE (nx,ny,nlanduse)
      REAL(fp), ALLOCATABLE        :: PBLH(:,:)         !< PBL height [m]
      ! 3D volumetric fields (3D: nx, ny, nz)
      REAL(fp), ALLOCATABLE        :: F_OF_PBL(:,:,:)       !< Fraction of box within PBL [1]
      REAL(fp), ALLOCATABLE        :: F_UNDER_PBLTOP(:,:,:) !< Fraction of box under PBL top
      real(fp), ALLOCATABLE        :: OBK(:,:)          !< Monin-Obhukov length [m]
      ! Cloud and precipitation properties (2D for surface, 3D for volumetric)
      REAL(fp), ALLOCATABLE        :: CLDFRC(:,:)       !< Column cloud fraction [1]
      REAL(fp), ALLOCATABLE        :: CONV_DEPTH(:,:)   !< Convective cloud depth [m]
      REAL(fp), ALLOCATABLE        :: FLASH_DENS(:,:)   !< Lightning flash density [#/km2/s]
      REAL(fp), ALLOCATABLE        :: CNV_FRC(:,:)      !< Convective fraction [1]
      REAL(fp), ALLOCATABLE        :: CLDF(:,:,:)       !< 3-D cloud fraction [1]
      REAL(fp), ALLOCATABLE        :: CMFMC(:,:,:)      !< Cloud mass flux [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: DQRCU(:,:,:)      !< Conv precip production rate [kg/kg/s] (assume per dry air)
      REAL(fp), ALLOCATABLE        :: DQRLSAN(:,:,:)    !< LS precip prod rate [kg/kg/s] (assume per dry air)
      REAL(fp), ALLOCATABLE        :: DTRAIN(:,:,:)     !< Detrainment flux [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PRECANV(:,:)      !< Anvil previp @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), ALLOCATABLE        :: PRECCON(:,:)      !< Conv  precip @ ground [kg/m2/s] -> [mm/day]
      REAL(fp), ALLOCATABLE        :: PRECLSC(:,:)      !< Large-scale precip @ ground kg/m2/s] -> [mm/day]
      REAL(fp), ALLOCATABLE        :: PRECTOT(:,:)      !< Total precip @ ground [kg/m2/s] -> [mm/day]
      ! 3D cloud and precipitation arrays
      REAL(fp), ALLOCATABLE        :: QI(:,:,:)         !< Mass fraction of cloud ice water [kg/kg dry air]
      REAL(fp), ALLOCATABLE        :: QL(:,:,:)         !< Mass fraction of cloud liquid water [kg/kg dry air]
      REAL(fp), ALLOCATABLE        :: PFICU(:,:,:)      !< Dwn flux ice prec:conv [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PFILSAN(:,:,:)    !< Dwn flux ice prec:LS+anv [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PFLCU(:,:,:)      !< Dwn flux liq prec:conv [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: PFLLSAN(:,:,:)    !< Dwn flux ice prec:LS+anv [kg/m2/s]
      REAL(fp), ALLOCATABLE        :: TAUCLI(:,:,:)     !< Opt depth of ice clouds [1]
      REAL(fp), ALLOCATABLE        :: TAUCLW(:,:,:)     !< Opt depth of H2O clouds [1]
      ! Surface scalars (now 2D: nx, ny)
      REAL(fp), ALLOCATABLE        :: PHIS(:,:)         !< Surface geopotential height [m2/s2]
      REAL(fp), ALLOCATABLE        :: PS_WET(:,:)       !< Wet surface pressure at start of timestep [hPa]
      REAL(fp), ALLOCATABLE        :: PS_DRY(:,:)       !< Dry surface pressure at start of timestep [hPa]
      REAL(fp), ALLOCATABLE        :: QV2M(:,:)         !< Specific Humidity at 2m [kg/kg]
      REAL(fp), ALLOCATABLE        :: T2M(:,:)          !< Temperature 2m [K]
      REAL(fp), ALLOCATABLE        :: TS(:,:)           !< Surface temperature [K]
      REAL(fp), ALLOCATABLE        :: TSKIN(:,:)        !< Surface skin temperature [K]
      REAL(fp), ALLOCATABLE        :: SST(:,:)          !< Sea surface temperature [K]
      REAL(fp), ALLOCATABLE        :: SLP(:,:)          !< Sea level pressure [hPa]
      REAL(fp), ALLOCATABLE        :: PS(:,:)           !< Surface Pressure [hPa]
      REAL(fp), ALLOCATABLE        :: TO3(:,:)          !< Total overhead O3 column [DU]
      REAL(fp), ALLOCATABLE        :: TROPP(:,:)        !< Tropopause pressure [hPa]
      INTEGER,  ALLOCATABLE        :: TropLev(:,:)      !< Tropopause level [1]
      REAL(fp), ALLOCATABLE        :: TropHt(:,:)       !< Tropopause height [km]
      ! 3D atmospheric variables (3D: nx, ny, nz)
      REAL(fp), ALLOCATABLE        :: Z(:,:,:)          !< Full Layer Geopotential Height
      REAL(fp), ALLOCATABLE        :: ZMID(:,:,:)       !< Mid Layer Geopotential Height
      REAL(fp), ALLOCATABLE        :: BXHEIGHT(:,:,:)   !< Grid box height [m] (dry air)
      REAL(fp), ALLOCATABLE        :: QV(:,:,:)         !< Specific Humidity [kg/kg]
      REAL(fp), ALLOCATABLE        :: T(:,:,:)          !< Temperature [K]
      REAL(fp), ALLOCATABLE        :: THETA(:,:,:)      !< Potential temperature [K]
      REAL(fp), ALLOCATABLE        :: TV(:,:,:)         !< Virtual temperature [K]
      REAL(fp), ALLOCATABLE        :: V(:,:,:)          !< N/S component of wind [m s-1]
      REAL(fp), ALLOCATABLE        :: U(:,:,:)          !< E/W component of wind [m s-1]
      REAL(fp), ALLOCATABLE        :: OMEGA(:,:,:)      !< Updraft velocity [Pa/s]
      REAL(fp), ALLOCATABLE        :: RH(:,:,:)         !< Relative humidity [%]
      REAL(fp), ALLOCATABLE        :: SPHU(:,:,:)       !< Specific humidity [g H2O/kg tot air]
      REAL(fp), ALLOCATABLE        :: AIRDEN(:,:,:)     !< Dry air density [kg/m3]
      REAL(fp), ALLOCATABLE        :: AIRNUMDEN(:,:,:)  !< Dry air density [molec/cm3]
      REAL(fp), ALLOCATABLE        :: MAIRDEN(:,:,:)    !< Moist air density [kg/m3]
      REAL(fp), ALLOCATABLE        :: AVGW(:,:,:)       !< Water vapor volume mixing ratio [vol H2O/vol dry air]
      REAL(fp), ALLOCATABLE        :: DELP(:,:,:)       !< Delta-P (wet) across box [hPa]
      REAL(fp), ALLOCATABLE        :: DELP_DRY(:,:,:)   !< Delta-P (dry) across box [hPa]
      REAL(fp), ALLOCATABLE        :: DAIRMASS(:,:,:)   !< Dry air mass [kg] in grid box
      REAL(fp), ALLOCATABLE        :: AIRVOL(:,:,:)     !< Grid box volume [m3] (dry air)
      REAL(fp), ALLOCATABLE        :: PEDGE_DRY(:,:,:)  !< Dry air partial pressure @ level edges [hPa] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: PEDGE(:,:,:)      !< Air partial pressure @ level edges [hPa] (nx,ny,nz+1)
      REAL(fp), ALLOCATABLE        :: PMID(:,:,:)       !< Average wet air pressure [hPa] defined as arithmetic average of edge pressures
      REAL(fp), ALLOCATABLE        :: PMID_DRY(:,:,:)   !< Dry air partial pressure [hPa] defined as arithmetic avg of edge pressures
      type(GridGeometryType) :: geometry !< Grid geometry configuration
   contains
      !> \copydoc metstate_init
      procedure :: init => metstate_init
      !> \copydoc metstate_cleanup
      procedure :: cleanup => metstate_cleanup
      !> \copydoc metstate_validate
      procedure :: validate => metstate_validate
      !> \copydoc metstate_reset
      procedure :: reset => metstate_reset
      !> \copydoc metstate_is_allocated
      procedure :: is_allocated => metstate_is_allocated
      !> \copydoc metstate_get_memory_usage
      procedure :: get_memory_usage => metstate_get_memory_usage
      !> \copydoc metstate_print_summary
      procedure :: print_summary => metstate_print_summary
      !> \copydoc metstate_get_dimensions
      procedure :: get_dimensions => metstate_get_dimensions
      !> \copydoc metstate_get_field_ptr
      procedure :: get_field_ptr => metstate_get_field_ptr
      !> \copydoc metstate_get_grid_geometry
      procedure :: get_grid_geometry => metstate_get_grid_geometry
      !> \copydoc metstate_get_column_ptr
      procedure :: get_column_ptr => metstate_get_column_ptr
      !> \copydoc allocate_arrays
      procedure, private :: allocate_arrays => allocate_metstate_arrays
   end type MetStateType

CONTAINS

   !> \brief Initialize a MetStateType object
   !!
   !! Initializes the meteorological state object, sets default values, and allocates required arrays.
   !!
   !! \param[inout] this      MetStateType object to initialize
   !! \param[in]    nlevs     Number of vertical levels
   !! \param[inout] error_mgr Error manager for context and error reporting
   !! \param[out]   rc        Return code (CC_SUCCESS or error code)
   subroutine metstate_init(this, nlevs, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_MEMORY_ALLOCATION

      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(in) :: nlevs
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      character(len=256) :: thisLoc

      thisLoc = 'metstate_init (in core/metstate_mod.F90)'
      call error_mgr%push_context('metstate_init', 'initializing meteorological state')

      rc = CC_SUCCESS

      this%NLEVS = nlevs
      this%State = 'MET'

      ! Call helper procedure to allocate arrays
      call this%allocate_arrays(nlevs, error_mgr, rc)

      call error_mgr%pop_context()
   end subroutine metstate_init

   !> \brief Allocate all arrays for MetStateType
   !!
   !! Helper procedure to allocate and initialize all arrays in the meteorological state.
   !!
   !! \param[inout] this      MetStateType object
   !! \param[in]    nlevs     Number of vertical levels
   !! \param[inout] error_mgr Error manager for context and error reporting
   !! \param[out]   rc        Return code
   subroutine allocate_metstate_arrays(this, nlevs, error_mgr, rc)
      use error_mod, only: ErrorManagerType, CC_SUCCESS, CC_FAILURE, ERROR_MEMORY_ALLOCATION

      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(in) :: nlevs
      type(ErrorManagerType), pointer, intent(inout) :: error_mgr
      integer, intent(out) :: rc
      integer :: nx, ny, allocStat
      character(len=256) :: thisLoc

      thisLoc = 'allocate_metstate_arrays (in core/metstate_mod.F90)'
      rc = CC_SUCCESS

      ! Get grid dimensions using public methods
      nx = this%geometry%get_nx()
      ny = this%geometry%get_ny()

      ! Allocate grid-dependent arrays
      allocate(this%IsLand(nx,ny), this%IsWater(nx,ny), this%IsIce(nx,ny), this%IsSnow(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate land/water/ice/snow flags', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif
      this%IsLand = .false.
      this%IsWater = .true.
      this%IsIce = .false.
      this%IsSnow = .false.

      allocate(this%AREA_M2(nx,ny), this%LWI(nx,ny), this%CLAYFRAC(nx,ny), this%DSOILTYPE(nx,ny), this%DLUSE(nx,ny), &
               this%FRVEG(nx,ny), this%FRLAKE(nx,ny), this%FRLAND(nx,ny), this%FRLANDIC(nx,ny), this%FROCEAN(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate surface properties 1', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      allocate(this%FRSEAICE(nx,ny), this%FRSNO(nx,ny), this%LAI(nx,ny), &
               this%GVF(nx,ny), this%RDRAG(nx,ny), &
               this%SANDFRAC(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate surface properties 2', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate more surface arrays
      allocate(this%SNODP(nx,ny), this%SNOMAS(nx,ny), this%SSM(nx,ny), this%USTAR_THRESHOLD(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate surface arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate sea ice coverage arrays
      allocate(this%SEAICE00(nx,ny), this%SEAICE10(nx,ny), this%SEAICE20(nx,ny), this%SEAICE30(nx,ny), &
               this%SEAICE40(nx,ny), this%SEAICE50(nx,ny), this%SEAICE60(nx,ny), this%SEAICE70(nx,ny), &
               this%SEAICE80(nx,ny), this%SEAICE90(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate sea ice arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate soil and land use arrays
      allocate(this%nLNDTYPE(nx,ny), this%GWETTOP(nx,ny), this%GWETROOT(nx,ny), this%WILT(nx,ny), &
               this%nSOIL(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate soil arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate location arrays
      allocate(this%LAT(nx,ny), this%LON(nx,ny), this%LUCNAME(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate location arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate surface radiation arrays
      allocate(this%ALBD_VIS(nx,ny), this%ALBD_NIR(nx,ny), this%ALBD_UV(nx,ny), this%PARDR(nx,ny), &
               this%PARDF(nx,ny), this%SUNCOS(nx,ny), this%SUNCOSmid(nx,ny), this%SUNCOSsum(nx,ny), &
               this%SZAFACT(nx,ny), this%SWGDN(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate radiation arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate surface flux and boundary layer arrays
      allocate(this%EFLUX(nx,ny), this%HFLUX(nx,ny), this%U10M(nx,ny), this%USTAR(nx,ny), &
               this%V10M(nx,ny), this%Z0(nx,ny), this%Z0H(nx,ny), this%PBLH(nx,ny), this%OBK(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate flux arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate cloud and precipitation 2D arrays
      allocate(this%CLDFRC(nx,ny), this%CONV_DEPTH(nx,ny), this%FLASH_DENS(nx,ny), this%CNV_FRC(nx,ny), &
               this%PRECANV(nx,ny), this%PRECCON(nx,ny), this%PRECLSC(nx,ny), this%PRECTOT(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate cloud 2D arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate 2D meteorological fields
      allocate(this%PHIS(nx,ny), this%PS_WET(nx,ny), this%PS_DRY(nx,ny), this%QV2M(nx,ny), this%T2M(nx,ny), &
               this%TS(nx,ny), this%TSKIN(nx,ny), this%SST(nx,ny), this%SLP(nx,ny), this%PS(nx,ny), &
               this%TO3(nx,ny), this%TROPP(nx,ny), this%TropLev(nx,ny), this%TropHt(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate 2D met fields', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate 3D vertical flag arrays
      allocate(this%InStratMeso(nx,ny,nlevs), this%InStratosphere(nx,ny,nlevs), &
               this%InTroposphere(nx,ny,nlevs), this%InPbl(nx,ny,nlevs), this%IsLocalNoon(nx,ny), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate 3D flag arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate 3D atmospheric arrays
      allocate(this%Z(nx,ny,nlevs), this%ZMID(nx,ny,nlevs), this%BXHEIGHT(nx,ny,nlevs), &
               this%T(nx,ny,nlevs), this%THETA(nx,ny,nlevs), this%TV(nx,ny,nlevs), &
               this%U(nx,ny,nlevs), this%V(nx,ny,nlevs), this%OMEGA(nx,ny,nlevs), &
               this%QV(nx,ny,nlevs), this%RH(nx,ny,nlevs), this%SPHU(nx,ny,nlevs), &
               this%AIRDEN(nx,ny,nlevs), this%AIRNUMDEN(nx,ny,nlevs), this%MAIRDEN(nx,ny,nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate 3D atmospheric arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate more 3D atmospheric arrays
      allocate(this%AVGW(nx,ny,nlevs), this%DELP(nx,ny,nlevs), this%DELP_DRY(nx,ny,nlevs), &
               this%DAIRMASS(nx,ny,nlevs), this%AIRVOL(nx,ny,nlevs), this%PMID(nx,ny,nlevs), &
               this%PMID_DRY(nx,ny,nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate 3D mass arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate edge pressure arrays (nlevs+1)
      allocate(this%PEDGE_DRY(nx,ny,nlevs+1), this%PEDGE(nx,ny,nlevs+1), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate edge pressure arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate 3D boundary layer arrays
      allocate(this%F_OF_PBL(nx,ny,nlevs), this%F_UNDER_PBLTOP(nx,ny,nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate PBL arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate 3D cloud arrays
      allocate(this%CLDF(nx,ny,nlevs), this%CMFMC(nx,ny,nlevs), this%DQRCU(nx,ny,nlevs), &
               this%DQRLSAN(nx,ny,nlevs), this%DTRAIN(nx,ny,nlevs), this%QI(nx,ny,nlevs), &
               this%QL(nx,ny,nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate 3D cloud arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Allocate 3D precipitation flux arrays
      allocate(this%PFICU(nx,ny,nlevs), this%PFILSAN(nx,ny,nlevs), this%PFLCU(nx,ny,nlevs), &
               this%PFLLSAN(nx,ny,nlevs), this%TAUCLI(nx,ny,nlevs), this%TAUCLW(nx,ny,nlevs), stat=allocStat)
      if (allocStat /= 0) then
         call error_mgr%report_error(ERROR_MEMORY_ALLOCATION, 'Failed to allocate 3D precip arrays', rc, thisLoc)
         call error_mgr%pop_context(); return
      endif

      ! Initialize to safe defaults
      this%T = 288.15_fp
      this%U = 0.0_fp
      this%V = 0.0_fp
      this%QV = 0.001_fp
      this%RH = 50.0_fp
      this%AIRDEN = 1.2_fp
      this%BXHEIGHT = 100.0_fp
      this%PS = 1013.25_fp
      this%SLP = 1013.25_fp
      this%T2M = 288.15_fp
      this%TS = 288.15_fp

   end subroutine allocate_metstate_arrays

   !> \brief Deallocate and clean up all arrays in MetStateType
   !!
   !! Deallocates all allocatable arrays and resets scalar values in the meteorological state.
   !!
   !! \param[inout] this MetStateType object
   !! \param[out]   rc   Return code (CC_SUCCESS)
   subroutine metstate_cleanup(this, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Deallocate all allocatable arrays
      ! Grid flags
      if (allocated(this%IsLand)) deallocate(this%IsLand)
      if (allocated(this%IsWater)) deallocate(this%IsWater)
      if (allocated(this%IsIce)) deallocate(this%IsIce)
      if (allocated(this%IsSnow)) deallocate(this%IsSnow)
      if (allocated(this%IsLocalNoon)) deallocate(this%IsLocalNoon)

      ! Vertical flags
      if (allocated(this%InStratMeso)) deallocate(this%InStratMeso)
      if (allocated(this%InStratosphere)) deallocate(this%InStratosphere)
      if (allocated(this%InTroposphere)) deallocate(this%InTroposphere)
      if (allocated(this%InPbl)) deallocate(this%InPbl)

      ! Surface properties
      if (allocated(this%AREA_M2)) deallocate(this%AREA_M2)
      if (allocated(this%LWI)) deallocate(this%LWI)
      if (allocated(this%CLAYFRAC)) deallocate(this%CLAYFRAC)
      if (allocated(this%DSOILTYPE)) deallocate(this%DSOILTYPE)
      if (allocated(this%DLUSE)) deallocate(this%DLUSE)
      if (allocated(this%FRVEG)) deallocate(this%FRVEG)
      if (allocated(this%FRLAKE)) deallocate(this%FRLAKE)
      if (allocated(this%FRLAND)) deallocate(this%FRLAND)
      if (allocated(this%FRLANDIC)) deallocate(this%FRLANDIC)
      if (allocated(this%FROCEAN)) deallocate(this%FROCEAN)
      if (allocated(this%FRSEAICE)) deallocate(this%FRSEAICE)
      if (allocated(this%FRSNO)) deallocate(this%FRSNO)
      if (allocated(this%LAI)) deallocate(this%LAI)
      if (allocated(this%GVF)) deallocate(this%GVF)
      if (allocated(this%RDRAG)) deallocate(this%RDRAG)
      if (allocated(this%SANDFRAC)) deallocate(this%SANDFRAC)

      ! Sea ice coverage
      if (allocated(this%SEAICE00)) deallocate(this%SEAICE00)
      if (allocated(this%SEAICE10)) deallocate(this%SEAICE10)
      if (allocated(this%SEAICE20)) deallocate(this%SEAICE20)
      if (allocated(this%SEAICE30)) deallocate(this%SEAICE30)
      if (allocated(this%SEAICE40)) deallocate(this%SEAICE40)
      if (allocated(this%SEAICE50)) deallocate(this%SEAICE50)
      if (allocated(this%SEAICE60)) deallocate(this%SEAICE60)
      if (allocated(this%SEAICE70)) deallocate(this%SEAICE70)
      if (allocated(this%SEAICE80)) deallocate(this%SEAICE80)
      if (allocated(this%SEAICE90)) deallocate(this%SEAICE90)
      if (allocated(this%SNODP)) deallocate(this%SNODP)
      if (allocated(this%SNOMAS)) deallocate(this%SNOMAS)
      if (allocated(this%SSM)) deallocate(this%SSM)
      if (allocated(this%USTAR_THRESHOLD)) deallocate(this%USTAR_THRESHOLD)

      ! Soil and land use
      if (allocated(this%nLNDTYPE)) deallocate(this%nLNDTYPE)
      if (allocated(this%GWETTOP)) deallocate(this%GWETTOP)
      if (allocated(this%GWETROOT)) deallocate(this%GWETROOT)
      if (allocated(this%WILT)) deallocate(this%WILT)
      if (allocated(this%nSOIL)) deallocate(this%nSOIL)
      if (allocated(this%SOILM)) deallocate(this%SOILM)
      if (allocated(this%FRLANDUSE)) deallocate(this%FRLANDUSE)
      if (allocated(this%FRSOIL)) deallocate(this%FRSOIL)
      if (allocated(this%FRLAI)) deallocate(this%FRLAI)
      if (allocated(this%FRZ0)) deallocate(this%FRZ0)

      ! Location arrays
      if (allocated(this%LAT)) deallocate(this%LAT)
      if (allocated(this%LON)) deallocate(this%LON)
      if (allocated(this%LUCNAME)) deallocate(this%LUCNAME)

      ! Surface radiation
      if (allocated(this%ALBD_VIS)) deallocate(this%ALBD_VIS)
      if (allocated(this%ALBD_NIR)) deallocate(this%ALBD_NIR)
      if (allocated(this%ALBD_UV)) deallocate(this%ALBD_UV)
      if (allocated(this%PARDR)) deallocate(this%PARDR)
      if (allocated(this%PARDF)) deallocate(this%PARDF)
      if (allocated(this%SUNCOS)) deallocate(this%SUNCOS)
      if (allocated(this%SUNCOSmid)) deallocate(this%SUNCOSmid)
      if (allocated(this%SUNCOSsum)) deallocate(this%SUNCOSsum)
      if (allocated(this%SZAFACT)) deallocate(this%SZAFACT)
      if (allocated(this%SWGDN)) deallocate(this%SWGDN)

      ! Surface fluxes
      if (allocated(this%EFLUX)) deallocate(this%EFLUX)
      if (allocated(this%HFLUX)) deallocate(this%HFLUX)
      if (allocated(this%U10M)) deallocate(this%U10M)
      if (allocated(this%USTAR)) deallocate(this%USTAR)
      if (allocated(this%V10M)) deallocate(this%V10M)
      if (allocated(this%Z0)) deallocate(this%Z0)
      if (allocated(this%Z0H)) deallocate(this%Z0H)
      if (allocated(this%PBLH)) deallocate(this%PBLH)
      if (allocated(this%OBK)) deallocate(this%OBK)

      ! Cloud and precipitation 2D
      if (allocated(this%CLDFRC)) deallocate(this%CLDFRC)
      if (allocated(this%CONV_DEPTH)) deallocate(this%CONV_DEPTH)
      if (allocated(this%FLASH_DENS)) deallocate(this%FLASH_DENS)
      if (allocated(this%CNV_FRC)) deallocate(this%CNV_FRC)
      if (allocated(this%PRECANV)) deallocate(this%PRECANV)
      if (allocated(this%PRECCON)) deallocate(this%PRECCON)
      if (allocated(this%PRECLSC)) deallocate(this%PRECLSC)
      if (allocated(this%PRECTOT)) deallocate(this%PRECTOT)

      ! 2D met fields
      if (allocated(this%PHIS)) deallocate(this%PHIS)
      if (allocated(this%PS_WET)) deallocate(this%PS_WET)
      if (allocated(this%PS_DRY)) deallocate(this%PS_DRY)
      if (allocated(this%QV2M)) deallocate(this%QV2M)
      if (allocated(this%T2M)) deallocate(this%T2M)
      if (allocated(this%TS)) deallocate(this%TS)
      if (allocated(this%TSKIN)) deallocate(this%TSKIN)
      if (allocated(this%SST)) deallocate(this%SST)
      if (allocated(this%SLP)) deallocate(this%SLP)
      if (allocated(this%PS)) deallocate(this%PS)
      if (allocated(this%TO3)) deallocate(this%TO3)
      if (allocated(this%TROPP)) deallocate(this%TROPP)
      if (allocated(this%TropLev)) deallocate(this%TropLev)
      if (allocated(this%TropHt)) deallocate(this%TropHt)

      ! 3D atmospheric variables
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

      ! PBL arrays
      if (allocated(this%F_OF_PBL)) deallocate(this%F_OF_PBL)
      if (allocated(this%F_UNDER_PBLTOP)) deallocate(this%F_UNDER_PBLTOP)

      ! 3D cloud arrays
      if (allocated(this%CLDF)) deallocate(this%CLDF)
      if (allocated(this%CMFMC)) deallocate(this%CMFMC)
      if (allocated(this%DQRCU)) deallocate(this%DQRCU)
      if (allocated(this%DQRLSAN)) deallocate(this%DQRLSAN)
      if (allocated(this%DTRAIN)) deallocate(this%DTRAIN)
      if (allocated(this%QI)) deallocate(this%QI)
      if (allocated(this%QL)) deallocate(this%QL)
      if (allocated(this%PFICU)) deallocate(this%PFICU)
      if (allocated(this%PFILSAN)) deallocate(this%PFILSAN)
      if (allocated(this%PFLCU)) deallocate(this%PFLCU)
      if (allocated(this%PFLLSAN)) deallocate(this%PFLLSAN)
      if (allocated(this%TAUCLI)) deallocate(this%TAUCLI)
      if (allocated(this%TAUCLW)) deallocate(this%TAUCLW)

      ! Reset scalar values
      this%NLEVS = 0
      this%State = ''

   end subroutine metstate_cleanup

   !> \brief Validate the MetStateType object for consistency and physical reasonableness
   !!
   !! Checks that the number of levels is positive, temperatures and pressures are within physical ranges,
   !! and that required arrays are allocated.
   !!
   !! \param[in]    this      MetStateType object
   !! \param[inout] error_mgr Error manager for context and error reporting
   !! \param[out]   rc        Return code (CC_SUCCESS or error code)
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

      ! Validate temperatures (use maxval/minval for array validation)
      if (allocated(this%T2M)) then
         if (maxval(this%T2M) > 400.0_fp .or. minval(this%T2M) < 100.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                        '2m temperature out of physical range', rc, &
                                        thisLoc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (allocated(this%TS)) then
         if (maxval(this%TS) > 400.0_fp .or. minval(this%TS) < 100.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                        'Surface temperature out of physical range', rc, &
                                        thisLoc, 'Check temperature units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      ! Validate pressures (use maxval/minval for array validation)
      if (allocated(this%PS)) then
         if (maxval(this%PS) > 1200.0_fp .or. minval(this%PS) < 10.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                        'Surface pressure out of physical range', rc, &
                                        thisLoc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
      endif

      if (allocated(this%SLP)) then
         if (maxval(this%SLP) > 1200.0_fp .or. minval(this%SLP) < 500.0_fp) then
            call error_mgr%report_error(ERROR_INVALID_INPUT, &
                                        'Sea level pressure out of physical range', rc, &
                                        thisLoc, 'Check pressure units and values')
            call error_mgr%pop_context()
            return
         endif
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
   !!
   !! Resets time, surface fields, and arrays to standard atmosphere values.
   !!
   !! \param[inout] this MetStateType object
   !! \param[out]   rc   Return code (CC_SUCCESS)
   subroutine metstate_reset(this, rc)
      implicit none
      class(MetStateType), intent(inout) :: this
      integer, intent(out) :: rc

      rc = CC_SUCCESS

      ! Reset to standard atmosphere values
      if (allocated(this%T2M)) this%T2M = 288.15_fp
      if (allocated(this%TS)) this%TS = 288.15_fp
      if (allocated(this%TSKIN)) this%TSKIN = 288.15_fp
      if (allocated(this%PS)) this%PS = 1013.25_fp
      if (allocated(this%SLP)) this%SLP = 1013.25_fp
      if (allocated(this%SST)) this%SST = 288.15_fp

      ! Reset arrays if allocated
      if (allocated(this%T)) this%T = 288.15_fp
      if (allocated(this%U)) this%U = 0.0_fp
      if (allocated(this%V)) this%V = 0.0_fp
      if (allocated(this%QV)) this%QV = 0.01_fp
      if (allocated(this%RH)) this%RH = 50.0_fp

   end subroutine metstate_reset

   !> \brief Check if required arrays are allocated in MetStateType
   !!
   !! Returns .true. if all required arrays are allocated, .false. otherwise.
   !!
   !! \param[in] this MetStateType object
   !! \return    Logical flag indicating allocation status
   function metstate_is_allocated(this) result(is_alloc)
      implicit none
      class(MetStateType), intent(in) :: this
      logical :: is_alloc

      is_alloc = allocated(this%T) .and. allocated(this%U) .and. allocated(this%V) .and. &
                 allocated(this%QV) .and. allocated(this%PMID) .and. allocated(this%DELP)
   end function metstate_is_allocated

   !> \brief Get approximate memory usage of MetStateType in bytes
   !!
   !! Estimates the memory usage of all arrays in the meteorological state object.
   !!
   !! \param[in] this MetStateType object
   !! \return    Estimated memory usage in bytes (integer(kind=8))
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

   !> \brief Print a summary of the MetStateType object to standard output
   !!
   !! Prints key fields, allocation status, and memory usage for diagnostics.
   !!
   !! \param[in] this MetStateType object
   subroutine metstate_print_summary(this)
      implicit none
      class(MetStateType), intent(in) :: this

      write(*,'(A)') '=== MetState Summary ==='
      write(*,'(A,A)') 'State: ', trim(this%State)
      write(*,'(A,I0)') 'Number of levels: ', this%NLEVS
      if (allocated(this%TS)) then
         write(*,'(A,F8.2,A)') 'Surface temperature: ', this%TS(1,1), ' K'
      endif
      if (allocated(this%PS)) then
         write(*,'(A,F8.2,A)') 'Surface pressure: ', this%PS(1,1), ' hPa'
      endif
      if (allocated(this%SLP)) then
         write(*,'(A,F8.2,A)') 'Sea level pressure: ', this%SLP(1,1), ' hPa'
      endif
      write(*,'(A,L1)') 'Arrays allocated: ', this%is_allocated()
      write(*,'(A,I0,A)') 'Memory usage: ', this%get_memory_usage(), ' bytes'
      write(*,'(A)') '======================='
   end subroutine metstate_print_summary

   !========================================================================
   !! Get grid dimensions for column interface support
   !========================================================================
   subroutine metstate_get_dimensions(this, nx, ny, nlev)
      class(MetStateType), intent(in) :: this
      integer, intent(out) :: nx, ny, nlev

      ! For current column-based architecture
      nx = 1
      ny = 1
      nlev = this%NLEVS

   end subroutine metstate_get_dimensions

   !========================================================================
   !! Get pointer to a meteorological field array by name
   !!
   !! Returns a pointer to the requested 1D field array (e.g., temperature, pressure) for the column model.
   !! Returns null() if the field is not found or not allocated.
   !!
   !! \param[in]  this       MetStateType object
   !! \param[in]  field_name Name of the field (e.g., 'T', 'temperature')
   !! \return     Pointer to the requested field array, or null()
   function metstate_get_field_ptr(this, field_name) result(field_ptr)
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      real(fp), pointer :: field_ptr(:)

      ! Return pointer to 1D column data
      field_ptr => null()

      select case (trim(field_name))
      case ('temperature', 'T')
         if (allocated(this%T)) field_ptr => this%T(1,1,:)
      case ('pressure', 'P', 'PMID')
         if (allocated(this%PMID)) field_ptr => this%PMID(1,1,:)
      case ('air_density', 'AIRDEN')
         if (allocated(this%AIRDEN)) field_ptr => this%AIRDEN(1,1,:)
      case ('u_wind', 'U')
         if (allocated(this%U)) field_ptr => this%U(1,1,:)
      case ('v_wind', 'V')
         if (allocated(this%V)) field_ptr => this%V(1,1,:)
      case ('specific_humidity', 'QV')
         if (allocated(this%QV)) field_ptr => this%QV(1,1,:)
      case ('relative_humidity', 'RH')
         if (allocated(this%RH)) field_ptr => this%RH(1,1,:)
      case ('box_height', 'BXHEIGHT')
         if (allocated(this%BXHEIGHT)) field_ptr => this%BXHEIGHT(1,1,:)
      case default
         ! Field not found
         field_ptr => null()
      end select

   end function metstate_get_field_ptr

   !========================================================================
   !! Get pointer to a vertical column for a given field at (i,j)
   !!
   !! Returns a pointer to the vertical profile for a given field at grid location (i,j).
   !! For column models, returns the full 1D array.
   !!
   !! \param[in]  this       MetStateType object
   !! \param[in]  field_name Name of the field (e.g., 'T', 'temperature')
   !! \param[in]  i          Grid column index (optional, default 1)
   !! \param[in]  j          Grid row index (optional, default 1)
   !! \return     Pointer to vertical profile (1D)
   function metstate_get_column_ptr(this, field_name, i, j) result(column_ptr)
      implicit none
      class(MetStateType), intent(in), target :: this
      character(len=*), intent(in) :: field_name
      integer, intent(in), optional :: i, j
      real(fp), pointer :: column_ptr(:)
      integer :: nx, ny, nlev, col_i, col_j
      column_ptr => null()
      call this%get_dimensions(nx, ny, nlev)
      col_i = 1; col_j = 1
      if (present(i)) col_i = max(1, min(i, nx))
      if (present(j)) col_j = max(1, min(j, ny))
      select case (trim(field_name))
      case ('temperature', 'T')
         if (allocated(this%T)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%T(1,1,:)
            else
               column_ptr => this%T(col_i, col_j, :)
            endif
         endif
      case ('pressure', 'P', 'PMID')
         if (allocated(this%PMID)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%PMID(1,1,:)
            else
               column_ptr => this%PMID(col_i, col_j, :)
            endif
         endif
      case ('air_density', 'AIRDEN')
         if (allocated(this%AIRDEN)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%AIRDEN(1,1,:)
            else
               column_ptr => this%AIRDEN(col_i, col_j, :)
            endif
         endif
      case ('u_wind', 'U')
         if (allocated(this%U)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%U(1,1,:)
            else
               column_ptr => this%U(col_i, col_j, :)
            endif
         endif
      case ('v_wind', 'V')
         if (allocated(this%V)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%V(1,1,:)
            else
               column_ptr => this%V(col_i, col_j, :)
            endif
         endif
      case ('specific_humidity', 'QV')
         if (allocated(this%QV)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%QV(1,1,:)
            else
               column_ptr => this%QV(col_i, col_j, :)
            endif
         endif
      case ('relative_humidity', 'RH')
         if (allocated(this%RH)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%RH(1,1,:)
            else
               column_ptr => this%RH(col_i, col_j, :)
            endif
         endif
      case ('box_height', 'BXHEIGHT')
         if (allocated(this%BXHEIGHT)) then
            if (nx == 1 .and. ny == 1) then
               column_ptr => this%BXHEIGHT(1,1,:)
            else
               column_ptr => this%BXHEIGHT(col_i, col_j, :)
            endif
         endif
      case default
         column_ptr => null()
      end select
   end function metstate_get_column_ptr

   !========================================================================
   !! Get the grid geometry object
   !========================================================================
   function metstate_get_grid_geometry(this) result(geometry)
      implicit none
      class(MetStateType), intent(in) :: this
      type(GridGeometryType) :: geometry

      geometry = this%geometry

   end function metstate_get_grid_geometry

!========================================================================
! Legacy Procedures (maintained for backward compatibility)
!========================================================================

END MODULE MetState_Mod
