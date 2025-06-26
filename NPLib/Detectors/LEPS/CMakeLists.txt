add_custom_command(OUTPUT TLEPSPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TLEPSPhysics.h TLEPSPhysicsDict.cxx TLEPSPhysics.rootmap libNPLEPS.dylib DEPENDS TLEPSPhysics.h)
add_custom_command(OUTPUT TLEPSDataDict.cxx COMMAND ../../scripts/build_dict.sh TLEPSData.h TLEPSDataDict.cxx TLEPSData.rootmap libNPLEPS.dylib DEPENDS TLEPSData.h)
add_library(NPLEPS SHARED TLEPSSpectra.cxx TLEPSData.cxx TLEPSPhysics.cxx TLEPSDataDict.cxx TLEPSPhysicsDict.cxx )
target_link_libraries(NPLEPS ${ROOT_LIBRARIES} NPCore) 
install(FILES TLEPSData.h TLEPSPhysics.h TLEPSSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

