add_custom_command(OUTPUT TEpicPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TEpicPhysics.h TEpicPhysicsDict.cxx TEpicPhysics.rootmap libNPEpic.dylib DEPENDS TEpicPhysics.h)
add_custom_command(OUTPUT TEpicDataDict.cxx COMMAND ../../scripts/build_dict.sh TEpicData.h TEpicDataDict.cxx TEpicData.rootmap libNPEpic.dylib DEPENDS TEpicData.h)
add_library(NPEpic SHARED TEpicSpectra.cxx TEpicData.cxx TEpicPhysics.cxx TEpicDataDict.cxx TEpicPhysicsDict.cxx )
target_link_libraries(NPEpic ${ROOT_LIBRARIES} NPCore) 
install(FILES TEpicData.h TEpicPhysics.h TEpicSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

