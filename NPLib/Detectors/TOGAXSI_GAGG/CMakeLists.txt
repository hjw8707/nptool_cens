add_custom_command(OUTPUT TTOGAXSI_GAGGPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TTOGAXSI_GAGGPhysics.h TTOGAXSI_GAGGPhysicsDict.cxx TTOGAXSI_GAGGPhysics.rootmap libNPTOGAXSI_GAGG.dylib DEPENDS TTOGAXSI_GAGGPhysics.h)
add_custom_command(OUTPUT TTOGAXSI_GAGGDataDict.cxx COMMAND ../../scripts/build_dict.sh TTOGAXSI_GAGGData.h TTOGAXSI_GAGGDataDict.cxx TTOGAXSI_GAGGData.rootmap libNPTOGAXSI_GAGG.dylib DEPENDS TTOGAXSI_GAGGData.h)
add_library(NPTOGAXSI_GAGG SHARED TTOGAXSI_GAGGSpectra.cxx TTOGAXSI_GAGGData.cxx TTOGAXSI_GAGGPhysics.cxx TTOGAXSI_GAGGDataDict.cxx TTOGAXSI_GAGGPhysicsDict.cxx )
target_link_libraries(NPTOGAXSI_GAGG ${ROOT_LIBRARIES} NPCore) 
install(FILES TTOGAXSI_GAGGData.h TTOGAXSI_GAGGPhysics.h TTOGAXSI_GAGGSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

