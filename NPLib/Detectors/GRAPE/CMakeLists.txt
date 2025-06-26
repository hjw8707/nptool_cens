add_custom_command(OUTPUT TGRAPEPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TGRAPEPhysics.h TGRAPEPhysicsDict.cxx TGRAPEPhysics.rootmap libNPGRAPE.dylib DEPENDS TGRAPEPhysics.h)
add_custom_command(OUTPUT TGRAPEDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TGRAPEData.h TGRAPEDataDict.cxx TGRAPEData.rootmap libNPGRAPE.dylib DEPENDS TGRAPEData.h)
add_library(NPGRAPE SHARED TGRAPEData.cxx TGRAPEPhysics.cxx TGRAPEDataDict.cxx TGRAPEPhysicsDict.cxx )
target_link_libraries(NPGRAPE ${ROOT_LIBRARIES} NPCore) 
install(FILES TGRAPEData.h TGRAPEPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

