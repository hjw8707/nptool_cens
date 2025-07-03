add_custom_command(OUTPUT TPlungerPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TPlungerPhysics.h TPlungerPhysicsDict.cxx TPlungerPhysics.rootmap libNPPlunger.dylib DEPENDS TPlungerPhysics.h)
add_custom_command(OUTPUT TPlungerDataDict.cxx COMMAND ../../scripts/build_dict.sh TPlungerData.h TPlungerDataDict.cxx TPlungerData.rootmap libNPPlunger.dylib DEPENDS TPlungerData.h)
add_library(NPPlunger SHARED TPlungerData.cxx TPlungerPhysics.cxx TPlungerDataDict.cxx TPlungerPhysicsDict.cxx )
target_link_libraries(NPPlunger ${ROOT_LIBRARIES} NPCore)
install(FILES TPlungerData.h TPlungerPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

