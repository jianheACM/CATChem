set(
  _srcs
  Precision_Mod.F90
  Constants.F90
  logging_mod.F90
  cmn_size_mod.F90
  TimeState_Mod.F90
  GridGeometry_Mod.F90
  GridManager_Mod.F90
  DiagnosticRegistry_Mod.F90
  EmissionConfigValidator_Mod.F90
  Species_Mod.F90
  ChemState_Mod.F90
  Met_Utilities_Mod.F90
  MetState_Mod.F90
  ExtEmisData_Mod.F90
  ConfigManager_Mod.F90
  EmisState_Mod.F90
  Error_Mod.F90
  Init_Mod.F90
  State_Interface.F90
  State_Mod.F90
  Utilities_Mod.F90
  ColumnInterface_Mod.F90
)

add_library(CATChem_core ${_srcs})
set_target_properties(
  CATChem_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_compile_options(CATChem_core PRIVATE ${_catchem_compiler_options})

# Link fyaml library for YAML configuration support
target_link_libraries(CATChem_core PRIVATE fyaml)
