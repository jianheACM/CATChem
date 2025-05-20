!> \file run_mod.F90
!! \brief Run module for the program.
!!
!! This module contains subroutines and functions related to the run of the program.
!! It includes subroutines for run emissions and other processes.
!!
!! \ingroup core_modules
!!!>
module run_mod
   USE MetState_Mod, ONLY : MetStateType
   USE ChemState_Mod, ONLY : ChemStateType
   USE DiagState_Mod, ONLY : DiagStateType
   USE EmisState_Mod, ONLY : EmisStateType
   USE CCPr_Dust_Common_Mod, ONLY : DustStateType
   USE CCPr_SeaSalt_Common_Mod, ONLY : SeaSaltStateType
   USE CCPR_DryDep_Mod, ONLY : DryDepStateType
   USE Error_Mod

   implicit none

   PRIVATE

   PUBLIC :: Run_Process
   PUBLIC :: Finalize_Process

contains


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
   subroutine Run_Process(MetState, DiagState, ChemState, EmisState, DustState, SeaSaltState, DryDepState, RC)
      !
      use CCPR_DryDep_Mod, Only : CCPR_DryDep_Run

      implicit none

      ! Arguments
      type(MetStateType), INTENT(IN) :: MetState
      type(DiagStateType), INTENT(IN) :: DiagState
      type(ChemStateType), INTENT(IN) :: ChemState
      TYPE(EmisStateType), INTENT(INOUT) :: EmisState
      type(DustStateType), INTENT(INOUT) :: DustState
      type(SeaSaltStateType), INTENT(INOUT) :: SeaSaltState
      type(DryDepStateType), INTENT(INOUT) :: DryDepState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Run_Process (in core/run_mod.F90)'

      !run all the emission processes
      call Run_Emis(MetState, DiagState, ChemState, EmisState, DustState, SeaSaltState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error running emission processes.'
         call CC_Error(errMsg, RC , thisLoc)
      endif

      !run dry deposition
      call CCPr_DryDep_Run( MetState, DiagState, DryDepState, ChemState, RC )
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error running dry deposition.'
         call CC_Error(errMsg, RC , thisLoc)
      endif

      !run wet deposition

   end subroutine Run_Process


   !> \brief Run the emission processes
   !!
   !! This subroutine runs the emission processes.
   !!
   !! \param MetState The meteorological state.
   !! \param DiagState The diagnostic state.
   !! \param ChemState The chemical state.
   !! \param DustState The dust state.
   !! \param SeaSaltState The sea salt state.
   !! \param RC The return code.
   !!
   !! \ingroup core_modules
   !!!>
   subroutine Run_Emis(MetState, DiagState, ChemState, EmisState, DustState, SeaSaltState, RC)
      use CCPr_Dust_mod, ONLY : CCPR_Dust_Run
      use CCPr_SeaSalt_mod, ONLY : CCPR_SeaSalt_Run
      use EmisState_Mod, only: Apply_Emis_to_Chem
      implicit none
      ! Arguments
      type(ChemStateType), intent(inout) :: ChemState
      type(DiagStateType), intent(inout) :: DiagState
      type(MetStateType), intent(inout)  :: MetState
      type(EmisStateType), intent(inout) :: EmisState
      type(DustStateType), intent(inout) :: DustState
      type(SeaSaltStateType), intent(inout) :: SeaSaltState
      integer, intent(inout) :: RC

      !local variables
      !integer :: i, k

      ! Error handling
      CHARACTER(LEN=255)    :: ErrMsg
      CHARACTER(LEN=255)    :: ThisLoc
      ThisLoc = ' -> at Run_Emis (in core/run_mo.F90)'


      !TODO: make sure each emission process has filled the EmisState flux since cc_apply_emis_to_chem is based on the EmisState
      ! Call dust emission
      !TODO: dust per bin emission only works on the default 5 biins of the Ginux scheme despite what defined in the emission namelist
      call CCPR_Dust_Run(MetState, DiagState, DustState, EmisState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error running dust emissions'
         call CC_Error(errMsg, RC , thisLoc)
         return
      end if

      ! Call sea salt emission
      call CCPR_SeaSalt_Run(MetState, SeaSaltState, EmisState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error running seal salt emissions'
         call CC_Error(errMsg, RC , thisLoc)
         return
      end if

      !apply emissions to chemistry state
      call Apply_Emis_to_Chem(EmisState, MetState, ChemState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error applying emissions to chemistry state'
         call CC_Error(errMsg, RC , thisLoc)
         return
      end if


   end subroutine Run_Emis


   subroutine Finalize_Process(DustState, SeaSaltState, DryDepState, RC)
      !
      use CCPR_Dust_Mod, Only : CCPr_Dust_Finalize
      use CCPR_SeaSalt_Mod, Only : CCPr_SeaSalt_Finalize
      use CCPR_DryDep_Mod, Only : CCPr_DryDep_Finalize

      implicit none

      ! Arguments
      type(DustStateType), INTENT(INOUT) :: DustState
      type(SeaSaltStateType), INTENT(INOUT) :: SeaSaltState
      type(DryDepStateType), INTENT(INOUT) :: DryDepState
      INTEGER,        INTENT(OUT) :: RC

      ! Local variables
      CHARACTER(LEN=255) :: ErrMsg, thisLoc

      ! Initialize
      RC = 0
      ErrMsg = ''
      thisLoc = ' -> at Finalize_Process (in core/run_mod.F90)'

      !finalize all the emission processes
      call CCPr_Dust_Finalize(DustState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error finalizing dust emissions.'
         call CC_Error(errMsg, RC , thisLoc)
      endif
      call CCPr_SeaSalt_Finalize(SeaSaltState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error finalizing sea salt emissions.'
         call CC_Error(errMsg, RC , thisLoc)
      endif

      !finalize dry deposition
      call CCPr_DryDep_Finalize(DryDepState, RC)
      if (RC /= CC_SUCCESS) then
         errMsg = 'Error finalizing dry deposition.'
         call CC_Error(errMsg, RC , thisLoc)
      endif

      !finalize wet deposition

   end subroutine Finalize_Process

end module run_mod
