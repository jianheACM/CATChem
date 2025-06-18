# Modern Dust Process CMakeLists.txt
# Following PROCESS_ARCHITECTURE_GUIDE.md structure

# Define dust process sources
set(DUST_PROCESS_SOURCES
    DustCommon_Mod.F90
    DustProcess_Mod.F90
    DustProcessCreator_Mod.F90
)

# Define dust scheme sources
set(DUST_SCHEME_SOURCES
    schemes/DustScheme_Fengsha_Mod.F90
    schemes/DustScheme_Ginoux_Mod.F90
)

# Combine all sources
set(DUST_ALL_SOURCES
    ${DUST_PROCESS_SOURCES}
    ${DUST_SCHEME_SOURCES}
)

# Create the dust process library
set(_lib CATChem_process_dust)
add_library(${_lib} ${DUST_ALL_SOURCES})

# Link with required libraries
target_link_libraries(${_lib} PUBLIC CATChem_core process_core)

# Set module directory
set_target_properties(${_lib} PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include
)

# Include directories
target_include_directories(${_lib} PRIVATE ${CMAKE_BINARY_DIR}/include)
