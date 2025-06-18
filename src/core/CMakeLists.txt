set(
  _srcs
  charpak_mod.F90
  chemstate_mod.F90
  cmn_size_mod.F90
  # complete_integration_demo.F90  # Temporarily disabled
  config/ConfigManager_Mod.F90
  # config_mod.F90  # Temporarily disabled - replaced by ConfigManager
  # config_opt_mod.F90  # Temporarily disabled - replaced by ConfigManager
  constants.F90
  # core_integration_demo.F90  # Temporarily disabled
  diagstate_mod.F90
  emisstate_mod.F90  # Re-enabled after GridState cleanup
  error_mod.F90
  # error_handling_example.F90  # Temporarily disabled
  # gridstate_mod.F90  # Removed - replaced by modern state initialization
  init_mod.F90
  metstate_mod.F90
  precision_mod.F90
  # qfyaml_mod.F90  # Temporarily disabled due to EmisState dependency
  species_mod.F90
  # state_interface.F90  # Temporarily disabled
  state_mod.F90
  utilities_mod.F90
)

add_library(CATChem_core ${_srcs})
set_target_properties(
  CATChem_core
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_compile_options(CATChem_core PRIVATE ${_catchem_compiler_options})

# Link fyaml library for YAML configuration support
target_link_libraries(CATChem_core PRIVATE fyaml)
