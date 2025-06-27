add_custom_command(OUTPUT TSTARKjrPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TSTARKjrPhysics.h TSTARKjrPhysicsDict.cxx TSTARKjrPhysics.rootmap libNPSTARKjr.dylib DEPENDS TSTARKjrPhysics.h)
add_custom_command(OUTPUT TSTARKjrDataDict.cxx COMMAND ../../scripts/build_dict.sh TSTARKjrData.h TSTARKjrDataDict.cxx TSTARKjrData.rootmap libNPSTARKjr.dylib DEPENDS TSTARKjrData.h)
add_library(NPSTARKjr SHARED TSTARKjrData.cxx TSTARKjrPhysics.cxx TSTARKjrDataDict.cxx TSTARKjrPhysicsDict.cxx)
target_link_libraries(NPSTARKjr ${ROOT_LIBRARIES} NPCore)
install(FILES TSTARKjrData.h TSTARKjrPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

