!Jian.He@noaa.gov, 08/2023
!Photolysis data for RACM

module phot_data_mod

use catchem_constants ,        only : kind_chem

IMPLICIT NONE

  integer, parameter ::  p_ph_o31d  = 1
  integer, parameter ::  p_ph_o33p  = 2
  integer, parameter ::  p_ph_no2   = 3
  integer, parameter ::  p_ph_no3o2 = 4
  integer, parameter ::  p_ph_no3o  = 5
  integer, parameter ::  p_ph_hno2  = 6
  integer, parameter ::  p_ph_hno3  = 7
  integer, parameter ::  p_ph_hno4  = 8
  integer, parameter ::  p_ph_h2o2  = 9
  integer, parameter ::  p_ph_ch2or = 10
  integer, parameter ::  p_ph_ch2om2 = 11
  integer, parameter ::  p_ph_ch3cho  = 12
  integer, parameter ::  p_ph_ch3coch3 = 13
  integer, parameter ::  p_ph_ch3coc2h5 = 14
  integer, parameter ::  p_ph_hcocho  = 15
  integer, parameter ::  p_ph_ch3cocho = 16
  integer, parameter ::  p_ph_hcochest = 17
  integer, parameter ::  p_ph_ch3o2h  = 18
  integer, parameter ::  p_ph_ch3coo2h = 19
  integer, parameter ::  p_ph_ch3ono2 = 20
  integer, parameter ::  p_ph_hcochob = 21
  integer, parameter ::  p_ph_macr = 22
  integer, parameter ::  p_ph_n2o5 = 23
  integer, parameter ::  p_ph_o2   = 24
  integer, parameter ::  p_ph_pan  = 25
  integer, parameter ::  p_ph_acet = 26
  integer, parameter ::  p_ph_mglo = 27
  integer, parameter ::  p_ph_hno4_2 = 28

end module phot_data_mod
