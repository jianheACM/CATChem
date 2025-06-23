!> \file ExternalEmissionCommon_Mod.F90
!! \brief Common utilities for ExternalEmission process
!! \ingroup externalemission_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! This module provides common utility functions for the ExternalEmission process.
!! It contains pure computational routines that are independent of state management.
!!
module ExternalEmissionCommon_Mod
   use precision_mod
   use error_mod

   implicit none
   private

   ! Public constants
   public :: EXTERNALEMISSION_VERSION
   public :: EXTERNALEMISSION_NAME

   ! Public utility functions

   ! Module constants
   real(fp), parameter :: EXTERNALEMISSION_VERSION = 1.0_fp
   character(len=*), parameter :: EXTERNALEMISSION_NAME = 'ExternalEmission'


contains


end module ExternalEmissionCommon_Mod