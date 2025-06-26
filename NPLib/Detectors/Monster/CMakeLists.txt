add_custom_command(OUTPUT TMonsterPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TMonsterPhysics.h TMonsterPhysicsDict.cxx TMonsterPhysics.rootmap libNPMonster.dylib DEPENDS TMonsterPhysics.h)
add_custom_command(OUTPUT TMonsterDataDict.cxx COMMAND ../../scripts/build_dict.sh TMonsterData.h TMonsterDataDict.cxx TMonsterData.rootmap libNPMonster.dylib DEPENDS TMonsterData.h)
add_library(NPMonster SHARED TMonsterSpectra.cxx TMonsterData.cxx TMonsterPhysics.cxx TMonsterDataDict.cxx TMonsterPhysicsDict.cxx )
target_link_libraries(NPMonster ${ROOT_LIBRARIES} NPCore) 
install(FILES TMonsterData.h TMonsterPhysics.h TMonsterSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

