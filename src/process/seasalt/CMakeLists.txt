# Modern SeaSalt Process CMakeLists.txt
# Following PROCESS_ARCHITECTURE_GUIDE.md structure (matching dust format)

# Define seasalt process sources
set(SEASALT_PROCESS_SOURCES
    SeaSaltCommon_Mod.F90
    SeaSaltProcess_Mod.F90
)

# Define seasalt scheme sources
set(SEASALT_SCHEME_SOURCES
    schemes/SeaSaltScheme_GEOS12_Mod.F90
    schemes/SeaSaltScheme_Gong03_Mod.F90
    schemes/SeaSaltScheme_Gong97_Mod.F90
)

# Combine all sources
set(SEASALT_ALL_SOURCES
    ${SEASALT_PROCESS_SOURCES}
    ${SEASALT_SCHEME_SOURCES}
)

# Create the seasalt process library
set(_lib CATChem_process_seasalt)
add_library(${_lib} ${SEASALT_ALL_SOURCES})

# Link with required libraries
target_link_libraries(${_lib} PUBLIC CATChem_core CATChem_process_interface)

# Set module directory
set_target_properties(
  ${_lib}
  PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)
target_link_libraries(${_lib} PUBLIC CATChem_process_seasalt_modern_schemes)
