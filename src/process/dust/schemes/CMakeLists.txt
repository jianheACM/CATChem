# Dust Schemes CMakeLists.txt
# Following PROCESS_ARCHITECTURE_GUIDE.md structure

# Define scheme sources
set(DUST_SCHEME_SOURCES
    DustScheme_Fengsha_Mod.F90
    DustScheme_Ginoux_Mod.F90
)

# Add to parent scope
set(DUST_SCHEME_SOURCES ${DUST_SCHEME_SOURCES} PARENT_SCOPE)
