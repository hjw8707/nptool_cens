add_custom_command(OUTPUT TGeLiPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TGeLiPhysics.h TGeLiPhysicsDict.cxx TGeLiPhysics.rootmap libNPGeLi.dylib DEPENDS TGeLiPhysics.h)
add_custom_command(OUTPUT TGeLiDataDict.cxx COMMAND ../../scripts/build_dict.sh TGeLiData.h TGeLiDataDict.cxx TGeLiData.rootmap libNPGeLi.dylib DEPENDS TGeLiData.h)
add_library(NPGeLi SHARED TGeLiSpectra.cxx TGeLiData.cxx TGeLiPhysics.cxx TGeLiDataDict.cxx TGeLiPhysicsDict.cxx )
target_link_libraries(NPGeLi ${ROOT_LIBRARIES} NPCore) 
install(FILES TGeLiData.h TGeLiPhysics.h TGeLiSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

