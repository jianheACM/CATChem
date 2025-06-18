# FYAMLConfig.cmake.in - Configuration file for FYAML
#
# This file is configured by CMake to create FYAMLConfig.cmake
# which can be used by other projects to find and use FYAML.


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was FYAMLConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

# Set FYAML version
set(FYAML_VERSION 0.2.0)

# Check if components are requested
set(_FYAML_SUPPORTED_COMPONENTS fyaml)
foreach(_comp ${FYAML_FIND_COMPONENTS})
    if(NOT _comp IN_LIST _FYAML_SUPPORTED_COMPONENTS)
        set(FYAML_FOUND False)
        set(FYAML_NOT_FOUND_MESSAGE "Unsupported component: ${_comp}")
    endif()
endforeach()

# Include the targets file
include("${CMAKE_CURRENT_LIST_DIR}/FYAMLTargets.cmake")

# Set variables for compatibility
set(FYAML_LIBRARIES FYAML::fyaml)
set(FYAML_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")

# Check that the targets actually exist
if(NOT TARGET FYAML::fyaml)
    set(FYAML_FOUND False)
    set(FYAML_NOT_FOUND_MESSAGE "FYAML targets not found")
    return()
endif()

# Print found message
if(NOT FYAML_FIND_QUIETLY)
    message(STATUS "Found FYAML: ${FYAML_VERSION}")
endif()

# Set FOUND to TRUE
set(FYAML_FOUND True)
