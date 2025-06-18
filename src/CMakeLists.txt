add_subdirectory(core)

# Always build external dependencies (fyaml is required for core)
add_subdirectory(external)

# Make process modules optional (some depend on external libraries)
option(BUILD_PROCESSES "Build process modules (may require external dependencies)" OFF)
if(BUILD_PROCESSES)
  add_subdirectory(process)
endif()

# API depends on process modules
option(BUILD_API "Build API (requires process modules)" OFF)
if(BUILD_API AND BUILD_PROCESSES)
  add_subdirectory(api)
endif()
