add_custom_command(OUTPUT TTOGAXSI_SIPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TTOGAXSI_SIPhysics.h TTOGAXSI_SIPhysicsDict.cxx TTOGAXSI_SIPhysics.rootmap libNPTOGAXSI_SI.dylib DEPENDS TTOGAXSI_SIPhysics.h)
add_custom_command(OUTPUT TTOGAXSI_SIDataDict.cxx COMMAND ../../scripts/build_dict.sh TTOGAXSI_SIData.h TTOGAXSI_SIDataDict.cxx TTOGAXSI_SIData.rootmap libNPTOGAXSI_SI.dylib DEPENDS TTOGAXSI_SIData.h)
add_library(NPTOGAXSI_SI SHARED TTOGAXSI_SISpectra.cxx TTOGAXSI_SIData.cxx TTOGAXSI_SIPhysics.cxx TTOGAXSI_SIDataDict.cxx TTOGAXSI_SIPhysicsDict.cxx )
target_link_libraries(NPTOGAXSI_SI ${ROOT_LIBRARIES} NPCore) 
install(FILES TTOGAXSI_SIData.h TTOGAXSI_SIPhysics.h TTOGAXSI_SISpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

