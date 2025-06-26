add_custom_command(OUTPUT TTinaPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TTinaPhysics.h TTinaPhysicsDict.cxx TTinaPhysics.rootmap libNPTina.dylib DEPENDS TTinaPhysics.h)
add_custom_command(OUTPUT TTinaDataDict.cxx COMMAND ../../scripts/build_dict.sh TTinaData.h TTinaDataDict.cxx TTinaData.rootmap libNPTina.dylib DEPENDS TTinaData.h)
add_library(NPTina SHARED TTinaData.cxx TTinaPhysics.cxx TTinaSpectra.cxx TTinaDataDict.cxx TTinaPhysicsDict.cxx )
target_link_libraries(NPTina ${ROOT_LIBRARIES} NPCore ${NPLib_LIBRARIES} -lNPReactionConditions) 
install(FILES TTinaData.h TTinaPhysics.h TTinaSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

