!> \file init_mod.F90
!! \brief Initialization module for the program.
!!
!! This module contains subroutines and functions related to the initialization of the program.
!! It includes subroutines for initializing the grid, the time step, and the solution.
!!
!! \ingroup core_modules
!!!>
module init_mod

   implicit none

   PUBLIC :: Init_Met
   PUBLIC :: Init_Diag
   PUBLIC :: Init_Emis
   PUBLIC :: Init_Chem
   PUBLIC :: Init_Process

contains

   !> \brief Initialize the met state
   !!
   !! This subroutine allocates the met state.
   !!
   !! \param GridState The grid state containing information about the grid.
   !! \param MetState The met state to be initialized.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Met(GridState, MetState, RC)
      !
      use GridState_Mod, Only : GridStateType
      use MetState_Mod
      USE Error_Mod
      implicit none

      ! Arguments
      TYPE(GridStateType), INTENT(IN)  :: GridState
      TYPE(MetStateType), INTENT(INOUT) :: MetState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Init_Met (in core/state_mod.F90)'

      call Met_Allocate(GridState, MetState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating met state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Met

   !> \brief Initialize the emission state
   !!
   !! This subroutine allocates the emission state.
   !!
   !! \param GridState The grid state containing information about the grid.
   !! \param EmisState The emission state to be initialized.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Emis(GridState, EmisState, RC)
      !
      use GridState_Mod, Only : GridStateType
      use EmisState_Mod
      USE Error_Mod
      implicit none

      ! Arguments
      TYPE(GridStateType), INTENT(IN)  :: GridState
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Init_Emis (in core/init_mod.F90)'

      call Emis_Allocate(GridState, EmisState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating emission state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Emis

   !> \brief Initialize the chem state
   !!
   !! This subroutine allocates the chem state.
   !!
   !! \param GridState The grid state containing information about the grid.
   !! \param ChemState The chem state to be initialized.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Chem(GridState, ChemState, RC)
      !
      use GridState_Mod, Only : GridStateType
      use ChemState_Mod
      USE Error_Mod
      implicit none

      ! Arguments
      TYPE(GridStateType), INTENT(IN)  :: GridState
      TYPE(ChemStateType), INTENT(INOUT) :: ChemState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Init_Chem (in core/init_mod.F90)'

      call Chem_Allocate(GridState, ChemState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating Chem state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Chem

   !> \brief Initialize the diag state
   !!
   !! This subroutine allocates the diag state.
   !!
   !! \param Config_Opt The config.
   !! \param GridState The grid state containing information about the grid.
   !! \param DiagState The diag state to be initialized.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Init_Diag(Config, DiagState, ChemState, RC)
      use DiagState_Mod
      use Config_Opt_Mod, Only : ConfigType
      ! use GridState_Mod, Only : GridStateType
      use ChemState_Mod, Only : ChemStateType
      use Error_Mod

      implicit none

      ! Arguments
      TYPE(ConfigType),    INTENT(IN)    :: Config
      ! TYPE(GridStateType), INTENT(IN)    :: GridState
      TYPE(DiagStateType), INTENT(INOUT) :: DiagState
      TYPE(ChemStateType), INTENT(INOUT)    :: ChemState
      INTEGER,         INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = CC_SUCCESS
      ErrMsg = ''
      thisLoc = ' -> at Init_Diag (in core/init_mod.F90)'

      call Diag_Allocate(Config, DiagState, ChemState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error allocating diag state'
         call CC_Error(errMsg, RC , thisLoc)
      endif

   end subroutine Init_Diag

   subroutine Init_Process (config, ChemState, EmisState, DustState, SeaSaltState, DryDepState, RC)
      use Config_Opt_Mod, only: ConfigType
      use Error_Mod
      use ChemState_Mod,  only: ChemStateType    !< Chemical State
      use EmisState_Mod,  only: EmisStateType    !< Emission State
      use CCPr_Dust_Common_Mod, only: DustStateType !< Dust State
      use CCPr_SeaSalt_Common_Mod, only: SeaSaltStateType !< SeaSalt State
      use CCPr_DryDep_mod, only: DryDepStateType !< DryDep State
      use CCPr_Dust_mod, only: CCPr_Dust_Init
      use CCPr_SeaSalt_mod, only: CCPr_SeaSalt_Init
      use CCPr_DryDep_mod, only: CCPr_DryDep_Init
      implicit none

      type(ConfigType), intent(in) :: config
      type(ChemStateType), intent(inout) :: ChemState
      type(EmisStateType), intent(inout) :: EmisState
      type(DustStateType), intent(inout) :: DustState
      type(SeaSaltStateType), intent(inout) :: SeaSaltState
      type(DryDepStateType), intent(inout) :: DryDepState
      integer, intent(inout) :: RC

      ! Error handling
      !---------------
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc
      ThisLoc = ' -> at init_process (in core/init_mod.F90)'

      ! Initialize DustState
      call CCPR_Dust_Init(config, DustState, ChemState, EmisState, RC)
      if (RC /= CC_SUCCESS) then
         ErrMsg = 'Error in CCPR_Dust_Init'
         call CC_Error(ErrMsg, RC, ThisLoc )
         return
      end if

      ! Initialize SealSaltState
      call CCPr_SeaSalt_Init(config, SeaSaltState, ChemState, EmisState, RC)
      if (RC /= CC_SUCCESS) then
         ErrMsg = 'Error in CCPr_SeaSalt_Init'
         call CC_Error(ErrMsg, RC, ThisLoc )
         return
      end if

      ! Initialize DryDepState
      call CCPr_DryDep_Init(config, DryDepState, ChemState, RC)
      if (RC /= CC_SUCCESS) then
         ErrMsg = 'Error in CCPr_DryDep_Init'
         call CC_Error(ErrMsg, RC, ThisLoc )
         return
      end if

   end subroutine Init_Process

end module init_mod
