add_custom_command(OUTPUT TGANILTacPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TGANILTacPhysics.h TGANILTacPhysicsDict.cxx TGANILTacPhysics.rootmap libNPGANIL.dylib DEPENDS TGANILTacPhysics.h)
add_custom_command(OUTPUT TGANILTriggerDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TGANILTriggerData.h TGANILTriggerDataDict.cxx TGANILTriggerData.rootmap libNPGANIL.dylib DEPENDS TGANILTriggerData.h)
add_custom_command(OUTPUT TGANILTacDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TGANILTacData.h TGANILTacDataDict.cxx TGANILTacData.rootmap libNPGANIL.dylib DEPENDS TGANILTacData.h)
add_library(NPGANIL SHARED TGANILTacPhysics.cxx TGANILTacData.cxx TGANILTriggerData.cxx TGANILTacDataDict.cxx TGANILTriggerDataDict.cxx TGANILTacPhysicsDict.cxx)
target_link_libraries(NPGANIL ${ROOT_LIBRARIES} NPCore) 
install(FILES TGANILTacData.h TGANILTriggerData.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

