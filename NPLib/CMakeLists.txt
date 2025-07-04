cmake_minimum_required (VERSION 3.10)
include(CheckCXXCompilerFlag)
project(NPLib CXX)
set(CMAKE_BUILD_TYPE Release)

# Setting the policy to match Cmake version
cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})

# include the nptool standard CMake preamble
include("${CMAKE_SOURCE_DIR}/ressources/CMake/NPTool_CMake_Preamble.cmake")

# Major change in the Core/Physics (affecting the detector/analysis/simulation)
set(NPLIB_VERSION_MAJOR 3)
# Minor change in the Core/Physics (not affecting any other part)
set(NPLIB_VERSION_MINOR 0)
# Change any of the detector in NPA
set(NPLIB_VERSION_DETA 0)

#activate Multithreading (on by default)
if(NOT DEFINED NPMULTITHREADING)
  set(NPMULTITHREADING 1)
  else()
  set(NPMULTITHREADING ${NPMULTITHREADING})
endif()

if(NPMULTITHREADING)
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNPMULTITHREADING=1")
 message("Building application with MultiThreading Support")
 else()
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNPMULTITHREADING=0")
 message("Building application with no MutilThreading Support")
endif()

configure_file(scripts/build_dict.sh.in scripts/build_dict.sh @ONLY)


configure_file(Core/NPLibVersion.h.in Core/NPLibVersion.h @ONLY)

string(COMPARE EQUAL "${NPTOOL_DETLIST}" "" rdet)
if(rdet)
    message("Building all detectors")
  else()
    message("Building the following detectors ${NPTOOL_DETLIST}")
endif()

set(CMAKE_BINARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/include )

# A macro that return the list of directory to compile
MACRO(subdirlist result curdir)
  FILE(GLOB children RELATIVE ${curdir} ${curdir}/Detectors/*)
  SET(dirlist "")
  FOREACH(child ${children})
    IF(IS_DIRECTORY ${curdir}/${child})
      string(COMPARE EQUAL "${child}" "CMakeFiles" r0)
      string(COMPARE EQUAL "${child}" ".DS_Store" r1)
      set(det "")
      string(REPLACE "Detectors/" "" det ${child})

      IF(NOT r0 AND NOT r1)
        string(COMPARE EQUAL "${NPTOOL_DETLIST}" "" r2)
        IF(r2)
          LIST(APPEND dirlist ${child})
        ELSEIF(${NPTOOL_DETLIST} MATCHES ".*${det}.*" )
          LIST(APPEND dirlist ${child})
        ENDIF()
      ENDIF()
    ENDIF()
  ENDFOREACH()
  SET(${result} ${dirlist})
ENDMACRO()

# Add include directories to all target
include_directories("Core/")
include_directories("Physics/")
include_directories("TrackReconstruction/")
include_directories("Online/")
include_directories("Neutron/")

# Call the macro
subdirlist(SUB_DIRECTORY ${CMAKE_BINARY_DIR})
set(SUB_DIRECTORY ${SUB_DIRECTORY} Core Physics TrackReconstruction Calibration Online Utility Neutron)
# Add each sub folder to the project
set(TARGET_LIST "")
foreach(subdir ${SUB_DIRECTORY})
# add the subdirectory
  add_subdirectory(${subdir})
  set(target "")
  string(REPLACE "Detectors/" "" target ${subdir})
  string(COMPARE EQUAL "${target}" "Utility" r0)
  string(COMPARE EQUAL "${target}" "Maya" r1)

  if(NOT r0 AND NOT r1)
    LIST(APPEND TARGET_LIST NP${target})
  endif()
  if(r1 AND DEFINED ENV{ROOTSYS})
	LIST(APPEND TARGET_LIST NP${target})
  endif()

endforeach()

export(TARGETS ${TARGET_LIST} FILE NPLibTargets.cmake)
export(PACKAGE NPLib)
file(RELATIVE_PATH REL_INCLUDE_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}" "${CMAKE_INCLUDE_OUTPUT_DIRECTORY}")
set(CONF_INCLUDE_DIRS "./" "./")
configure_file(ressources/CMake/NPLibConfig.cmake.in "${PROJECT_BINARY_DIR}/NPLibConfig.cmake" @ONLY)
configure_file(scripts/nptool-cleaner.sh "${CMAKE_BINARY_OUTPUT_DIRECTORY}/nptool-cleaner" @ONLY)
configure_file(scripts/nptool-compilation.sh "${CMAKE_BINARY_OUTPUT_DIRECTORY}/npcompilation" @ONLY)
configure_file(scripts/nptool-wizard.sh "${CMAKE_BINARY_OUTPUT_DIRECTORY}/nptool-wizard" @ONLY)
