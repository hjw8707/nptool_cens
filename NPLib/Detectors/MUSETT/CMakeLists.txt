add_custom_command(OUTPUT TMUSETTPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TMUSETTPhysics.h TMUSETTPhysicsDict.cxx TMUSETTPhysics.rootmap libNPMUSETT.dylib DEPENDS TMUSETTPhysics.h)
add_custom_command(OUTPUT TMUSETTDataDict.cxx COMMAND ../../scripts/build_dict.sh TMUSETTData.h TMUSETTDataDict.cxx TMUSETTData.rootmap libNPMUSETT.dylib DEPENDS TMUSETTData.h)
add_library(NPMUSETT SHARED TMUSETTData.cxx TMUSETTPhysics.cxx TMUSETTDataDict.cxx TMUSETTPhysicsDict.cxx TMUSETTSpectra.h TMUSETTSpectra.cxx)
target_link_libraries(NPMUSETT ${ROOT_LIBRARIES} NPCore) 
install(FILES MUSETTMap.h TMUSETTSpectra.h TMUSETTData.h TMUSETTPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

