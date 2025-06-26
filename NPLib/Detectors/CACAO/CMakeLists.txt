add_custom_command(OUTPUT TCACAOPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TCACAOPhysics.h TCACAOPhysicsDict.cxx TCACAOPhysics.rootmap libNPCACAO.dylib DEPENDS TCACAOPhysics.h)
add_custom_command(OUTPUT TCACAODataDict.cxx COMMAND ../../scripts/build_dict.sh TCACAOData.h TCACAODataDict.cxx TCACAOData.rootmap libNPCACAO.dylib DEPENDS TCACAOData.h)
add_library(NPCACAO SHARED TCACAOSpectra.cxx TCACAOData.cxx TCACAOPhysics.cxx TCACAODataDict.cxx TCACAOPhysicsDict.cxx )
target_link_libraries(NPCACAO ${ROOT_LIBRARIES} NPCore)
install(FILES TCACAOData.h TCACAOPhysics.h TCACAOSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

