add_custom_command(OUTPUT TDali2PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TDali2Physics.h TDali2PhysicsDict.cxx TDali2Physics.rootmap libNPDali2.dylib DEPENDS TDali2Physics.h)
add_custom_command(OUTPUT TDali2DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TDali2Data.h TDali2DataDict.cxx TDali2Data.rootmap libNPDali2.dylib DEPENDS TDali2Data.h)
add_library(NPDali2 SHARED TDali2Spectra.cxx TDali2Data.cxx TDali2Physics.cxx TDali2DataDict.cxx TDali2PhysicsDict.cxx )
target_link_libraries(NPDali2 ${ROOT_LIBRARIES} NPCore) 
install(FILES TDali2Data.h TDali2Physics.h TDali2Spectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

