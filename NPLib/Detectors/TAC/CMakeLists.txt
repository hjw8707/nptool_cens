add_custom_command(OUTPUT TTACPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TTACPhysics.h TTACPhysicsDict.cxx TTACPhysics.rootmap libNPTAC.dylib DEPENDS TTACPhysics.h)
add_custom_command(OUTPUT TTACDataDict.cxx COMMAND ../../scripts/build_dict.sh TTACData.h TTACDataDict.cxx TTACData.rootmap libNPTAC.dylib DEPENDS TTACData.h)
add_custom_command(OUTPUT TTACPhysicsReaderDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TTACPhysicsReader.h TTACPhysicsReaderDict.cxx TTACPhysicsReader.rootmap libNPTAC.dylib DEPENDS TTACPhysicsReader.h)
add_library(NPTAC SHARED TTACSpectra.cxx TTACData.cxx TTACPhysics.cxx TTACPhysicsReader.cxx TTACDataDict.cxx TTACPhysicsDict.cxx TTACPhysicsReaderDict.cxx)
target_link_libraries(NPTAC ${ROOT_LIBRARIES} -lTreePlayer -lSpectrum NPCore NPPhysics) 
install(FILES TTACData.h TTACPhysics.h TTACSpectra.h TTACPhysicsReader.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

