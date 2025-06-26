add_custom_command(OUTPUT TSuperX3PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSuperX3Physics.h TSuperX3PhysicsDict.cxx TSuperX3Physics.rootmap libNPSuperX3.dylib DEPENDS TSuperX3Physics.h)
add_custom_command(OUTPUT TSuperX3DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TSuperX3Data.h TSuperX3DataDict.cxx TSuperX3Data.rootmap libNPSuperX3.dylib DEPENDS TSuperX3Data.h)
add_library(NPSuperX3 SHARED TSuperX3Spectra.cxx TSuperX3Data.cxx TSuperX3Physics.cxx TSuperX3DataDict.cxx TSuperX3PhysicsDict.cxx )
target_link_libraries(NPSuperX3 ${ROOT_LIBRARIES} NPCore) 
install(FILES TSuperX3Data.h TSuperX3Physics.h TSuperX3Spectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

