add_custom_command(OUTPUT TSTARKPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TSTARKPhysics.h TSTARKPhysicsDict.cxx TSTARKPhysics.rootmap libNPSTARK.dylib DEPENDS TSTARKPhysics.h)
add_custom_command(OUTPUT TSTARKDataDict.cxx COMMAND ../../scripts/build_dict.sh TSTARKData.h TSTARKDataDict.cxx TSTARKData.rootmap libNPSTARK.dylib DEPENDS TSTARKData.h)
add_custom_command(OUTPUT TSTARKRawDict.cxx COMMAND ../../scripts/build_dict.sh TSTARKRaw.h TSTARKRawDict.cxx TSTARKRaw.rootmap libNPSTARK.dylib DEPENDS TSTARKRaw.h)
add_library(NPSTARK SHARED TSTARKData.cxx TSTARKDataDict.cxx TSTARKPhysics.cxx TSTARKPhysicsDict.cxx TSTARKRaw.cxx TSTARKRawDict.cxx )
target_link_libraries(NPSTARK ${ROOT_LIBRARIES} NPCore ${NPLib_LIBRARIES}) 
install(FILES TSTARKData.h TSTARKPhysics.h TSTARKRaw.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

