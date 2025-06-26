add_custom_command(OUTPUT TASGARDPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TASGARDPhysics.h TASGARDPhysicsDict.cxx TASGARDPhysics.rootmap libNPASGARD.dylib DEPENDS TASGARDPhysics.h)
add_custom_command(OUTPUT TASGARDDataDict.cxx COMMAND ../../scripts/build_dict.sh TASGARDData.h TASGARDDataDict.cxx TASGARDData.rootmap libNPASGARD.dylib DEPENDS TASGARDData.h)
add_library(NPASGARD SHARED TASGARDData.cxx TASGARDPhysics.cxx TASGARDDataDict.cxx TASGARDPhysicsDict.cxx )
target_link_libraries(NPASGARD ${ROOT_LIBRARIES} NPCore) 
install(FILES TASGARDData.h TASGARDPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

