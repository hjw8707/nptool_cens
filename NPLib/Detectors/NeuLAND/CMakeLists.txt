add_custom_command(OUTPUT TNeuLANDPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TNeuLANDPhysics.h TNeuLANDPhysicsDict.cxx TNeuLANDPhysics.rootmap libNPNeuLAND.dylib DEPENDS TNeuLANDPhysics.h)
add_custom_command(OUTPUT TNeuLANDDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TNeuLANDData.h TNeuLANDDataDict.cxx TNeuLANDData.rootmap libNPNeuLAND.dylib DEPENDS TNeuLANDData.h)
add_library(NPNeuLAND SHARED TNeuLANDSpectra.cxx TNeuLANDData.cxx TNeuLANDPhysics.cxx TNeuLANDDataDict.cxx TNeuLANDPhysicsDict.cxx )
target_link_libraries(NPNeuLAND ${ROOT_LIBRARIES} NPCore) 
install(FILES TNeuLANDData.h TNeuLANDPhysics.h TNeuLANDSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

