!> \file DryDepCommon_Mod.F90
!! \brief Common utilities for DryDep process
!! \ingroup drydep_modules
!!
!! \author CATChem Development Team
!! \date 2025
!! \version 1.0
!!
!! This module provides common utility functions for the DryDep process.
!! It contains pure computational routines that are independent of state management.
!!
module DryDepCommon_Mod
   use precision_mod
   use error_mod

   implicit none
   private

   ! Public constants
   public :: DRYDEP_VERSION
   public :: DRYDEP_NAME

   ! Public utility functions

   ! Module constants
   real(fp), parameter :: DRYDEP_VERSION = 1.0_fp
   character(len=*), parameter :: DRYDEP_NAME = 'DryDep'


contains


end module DryDepCommon_Mod