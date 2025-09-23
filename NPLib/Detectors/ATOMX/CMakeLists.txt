add_custom_command(OUTPUT TATOMXPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TATOMXPhysics.h TATOMXPhysicsDict.cxx TATOMXPhysics.rootmap libNPATOMX.dylib DEPENDS TATOMXPhysics.h)
add_custom_command(OUTPUT TATOMXDataDict.cxx COMMAND ../../scripts/build_dict.sh TATOMXData.h TATOMXDataDict.cxx TATOMXData.rootmap libNPATOMX.dylib DEPENDS TATOMXData.h)
add_library(NPATOMX SHARED TATOMXSpectra.cxx TATOMXData.cxx TATOMXPhysics.cxx TATOMXDataDict.cxx TATOMXPhysicsDict.cxx )
target_link_libraries(NPATOMX ${ROOT_LIBRARIES} NPCore) 
install(FILES TATOMXData.h TATOMXPhysics.h TATOMXSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

