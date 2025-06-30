# This file find your active Cubix installation and use it for configuration

message("Looking for Cubix:")
set(NPTOOL_CUBIX_SEARCHPATH 
    $ENV{CUBIX_SYS}/bin 
    $ENV{PATH}
    # /usr/bin  # Debian install
    # /opt/local/bin # Macports install
    # /usr/local/bin
    )

find_program(NPTOOL_CUBIX_CONFIG NAMES cubix-config PATHS ${NPTOOL_CUBIX_SEARCHPATH} )
if(NPTOOL_CUBIX_CONFIG MATCHES "NOTFOUND")
  message(FATAL_ERROR "Cubix is not found, check your installation")
endif()

execute_process(OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND ${NPTOOL_CUBIX_CONFIG} "--version" OUTPUT_VARIABLE NPTOOL_CUBIX_VERSION)
message(" Cubix version ${NPTOOL_CUBIX_VERSION} found")
 
## Collect the different information ##

# List of Root dependencies
set(CUBIX_LIBRARIES CubixCore CubixTools)

# Lib directories
execute_process(OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND ${NPTOOL_CUBIX_CONFIG} "--libdir" OUTPUT_VARIABLE CUBIX_LIBRARY_DIR)
link_directories( ${CUBIX_LIBRARY_DIR})  

# Include directories
execute_process(OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND ${NPTOOL_CUBIX_CONFIG} "--incdir" OUTPUT_VARIABLE CUBIX_INCLUDE_DIR)
include_directories(SYSTEM ${CUBIX_INCLUDE_DIR})  
message(${CUBIX_INCLUDE_DIR})

# Get the compilator flag from root to assure consistancy
execute_process(OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND ${NPTOOL_CUBIX_CONFIG} "--cflags" OUTPUT_VARIABLE NPTOOL_CUBIX_CFLAGS)
message(${NPTOOL_CUBIX_CFLAGS})
if(CMAKE_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${NPTOOL_CUBIX_CFLAGS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS_RELEASE} ${NPTOOL_CUBIX_CFLAGS}")
endif()
