! \file catchem_nuopc_cap.F90
! \brief NUOPC cap for CATChem atmospheric chemistry model
!>
! \defgroup catchem_nuopc_group CATChem NUOPC Interface
! \ingroup catchem
!>
! \details
! This module provides the NUOPC (National Unified Operational Prediction
! Capability) cap for the CATChem (Configurable ATmospheric Chemistry) model.
! The cap enables CATChem to run within the NUOPC/ESMF framework as a
! component in coupled Earth system models.
!>
! The cap implements the standard NUOPC phases:
! - Initialize: Set up the component, advertise fields, and realize the grid
! - Run: Execute chemistry calculations and data exchange
! - Finalize: Clean up resources
!>
! \author Barry Baker
! \date 11/2024

module catchem_nuopc_cap

  use ESMF
  use NUOPC
  use NUOPC_Model, &
    modelSS        => SetServices, &
    model_label_Advance => label_Advance, &
    model_label_CheckImport => label_CheckImport, &
    model_label_SetRunClock => label_SetRunClock

  use CATChem
  use catchem_types, only: catchem_container_type
  use catchem_nuopc_interface
  use catchem_nuopc_utils

  implicit none

  private

  public :: SetServices

  ! CATChem component data
  type(catchem_container_type), save :: catchem_states
  type(ConfigType), save :: config
  type(DustStateType), save :: dustState
  type(SeaSaltStateType), save :: seaSaltState
  type(DryDepStateType), save :: dryDepState

  ! Component configuration
  character(len=256), save :: config_file = 'CATChem_config.yml'
  logical, save :: do_chemistry = .true.

contains

  ! Set services for the CATChem NUOPC cap
  !!
  !! This is the main entry point called by NUOPC to set up the component.
  !! It registers the initialize, run, and finalize phases.
  !!
  !! \param model  NUOPC model component
  !! \param   rc     ESMF return code
  !!
  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! Set the model services
    call NUOPC_CompDerive(model, modelSS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Set component metadata
    call ESMF_AttributeSet(model, name="model_name", value="CATChem", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_AttributeSet(model, name="model_version", value="1.0", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Set entry points for standard NUOPC phases
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_RUN, &
      phaseLabelList=(/"RunPhase1"/), userRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_FINALIZE, &
      phaseLabelList=(/"FinalizePhase1"/), userRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

  end subroutine SetServices

  ! Initialize Phase 1 - Advertise fields
  !!
  !! In this phase, the component advertises what fields it can import
  !! (receive from other components) and export (provide to other components).
  !!
  !! \param model  NUOPC model component
  !! \param   rc     ESMF return code
  !!
  subroutine InitializeP1(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    type(ESMF_State) :: importState, exportState
    character(len=*), parameter :: routine = 'InitializeP1'

    rc = ESMF_SUCCESS

    ! Get import and export states
    call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Advertise fields
    call catchem_advertise_fields(importState, exportState, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Log successful completion
    call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeP1

  ! Initialize Phase 2 - Realize fields and initialize model
  !!
  !! In this phase, the component creates the actual ESMF fields for the
  !! advertised imports and exports, and initializes the CATChem model.
  !!
  !! \param model  NUOPC model component
  !! \param   rc     ESMF return code
  !!
  subroutine InitializeP2(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    type(ESMF_State) :: importState, exportState
    type(ESMF_Grid) :: grid
    type(ESMF_Config) :: config_esmf
    character(len=*), parameter :: routine = 'InitializeP2'
    character(len=512) :: errmsg
    integer :: localPet, petCount
    integer :: im, jm  ! Grid dimensions

    rc = ESMF_SUCCESS

    ! Get component information
    call ESMF_GridCompGet(model, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get import and export states
    call NUOPC_ModelGet(model, importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get or create grid (this should be provided by the driver)
    call ESMF_StateGet(importState, "grid", grid, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      ! Create a simple grid if not provided
      grid = ESMF_GridCreateNoPeriDim(minIndex=(/1,1/), maxIndex=(/144,91/), &
        regDecomp=(/petCount,1/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=__FILE__)) return
    end if

    ! Get grid dimensions
    call ESMF_GridGet(grid, tile=1, staggerloc=ESMF_STAGGERLOC_CENTER, &
      computationalCount=(/im, jm/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Realize fields
    call catchem_realize_fields(importState, exportState, grid, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Read configuration
    config_esmf = ESMF_ConfigCreate(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_ConfigLoadFile(config_esmf, filename=config_file, rc=rc)
    if (rc /= ESMF_SUCCESS) then
      call ESMF_LogWrite("CATChem: Could not load config file, using defaults", &
        ESMF_LOGMSG_WARNING, rc=rc)
      rc = ESMF_SUCCESS
    end if

    ! Initialize CATChem using the interface
    call catchem_nuopc_init(config, catchem_states, dustState, seaSaltState, &
                           dryDepState, im*jm, config_file, grid, rc, errmsg)
    if (rc /= CC_SUCCESS) then
      call ESMF_LogWrite("CATChem: Failed to initialize - " // trim(errmsg), &
        ESMF_LOGMSG_ERROR, rc=rc)
      rc = ESMF_FAILURE
      return
    end if

    ! Clean up
    call ESMF_ConfigDestroy(config_esmf, rc=rc)

    ! Log successful completion
    call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)

  end subroutine InitializeP2

  ! Model advance routine - Run the chemistry model
  !!
  !! This routine is called at each time step to advance the chemistry model.
  !! It imports meteorological data, runs chemistry calculations, and exports
  !! the results.
  !!
  !! \param model  NUOPC model component
  !! \param   rc     ESMF return code
  !!
  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    type(ESMF_State) :: importState, exportState
    type(ESMF_Clock) :: clock
    type(ESMF_Time) :: currTime
    type(ESMF_TimeInterval) :: timeStep
    character(len=*), parameter :: routine = 'ModelAdvance'
    character(len=512) :: errmsg
    integer :: localPet, im, jm, kme
    real(ESMF_KIND_R8) :: dt_seconds

    rc = ESMF_SUCCESS

    ! Get component information
    call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get states and clock
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Get current time and time step
    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    call ESMF_TimeIntervalGet(timeStep, s_r8=dt_seconds, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (localPet == 0) then
      call ESMF_LogWrite("CATChem: Running chemistry for dt = " // &
        trim(adjustl(real_to_string(dt_seconds))) // " seconds", &
        ESMF_LOGMSG_INFO, rc=rc)
    end if

    ! Get grid dimensions (placeholder - should come from grid)
    im = catchem_states%im
    jm = 1  ! Assume 1D for now
    kme = 1  ! Will be updated based on actual grid

    ! Import meteorological data from other components
    call transform_nuopc_to_catchem(importState, catchem_states, im, kme, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    if (do_chemistry) then
      ! Run chemistry processes using new interface with current time
      call catchem_nuopc_run(config, catchem_states, dustState, seaSaltState, &
                            dryDepState, dt_seconds, currTime, rc, errmsg)
      if (rc /= CC_SUCCESS) then
        call ESMF_LogWrite("CATChem: Failed to run chemistry - " // trim(errmsg), &
          ESMF_LOGMSG_ERROR, rc=rc)
        rc = ESMF_FAILURE
        return
      end if
    end if

    ! Export results to other components
    call transform_catchem_to_nuopc(exportState, catchem_states, im, kme, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Log successful completion
    if (localPet == 0) then
      call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)
    end if

  end subroutine ModelAdvance

  ! Run chemistry processes
  !!
  !! \param  dt  Time step in seconds
  !! \param rc  Return code
  !!
  subroutine run_chemistry_processes(dt, rc)
    real(ESMF_KIND_R8), intent(in) :: dt
    integer, intent(out) :: rc

    integer :: i

    rc = CC_SUCCESS

    do i = 1, catchem_states%im
      ! Run dust process
      call cc_dust_run(dustState, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return

      ! Run seasalt process
      call cc_seasalt_run(seaSaltState, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return

      ! Run dry deposition process
      call cc_drydep_run(dryDepState, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return

      ! Run other chemistry processes as needed
      call cc_run_process(config, catchem_states%GridState, &
        catchem_states%MetState(i), catchem_states%ChemState(i), &
        catchem_states%EmisState(i), catchem_states%DiagState(i), &
        dt, rc)
      if (rc /= CC_SUCCESS) return
    end do

  end subroutine run_chemistry_processes

  ! Finalize the model component
  !!
  !! Clean up resources and finalize CATChem processes.
  !!
  !! \param model  NUOPC model component
  !! \param   rc     ESMF return code
  !!
  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    character(len=*), parameter :: routine = 'ModelFinalize'
    character(len=512) :: errmsg
    integer :: localPet

    rc = ESMF_SUCCESS

    ! Get component information
    call ESMF_GridCompGet(model, localPet=localPet, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return

    ! Finalize CATChem using the interface
    call catchem_nuopc_finalize(config, catchem_states, dustState, seaSaltState, &
                               dryDepState, rc, errmsg)
    if (rc /= CC_SUCCESS) then
      call ESMF_LogWrite("CATChem: Warning - " // trim(errmsg), &
        ESMF_LOGMSG_WARNING, rc=rc)
    end if

    ! Log successful completion
    if (localPet == 0) then
      call ESMF_LogWrite("CATChem: Completed "//routine, ESMF_LOGMSG_INFO, rc=rc)
    end if

    rc = ESMF_SUCCESS

  end subroutine ModelFinalize

  ! Convert real to string for logging
  !!
  !! \param val  Real value to convert
  !! \return         String representation
  !!
  function real_to_string(val) result(str)
    real(ESMF_KIND_R8), intent(in) :: val
    character(len=32) :: str

    write(str, '(f0.2)') val
    str = adjustl(str)
  end function real_to_string

end module catchem_nuopc_cap
