add_custom_command(OUTPUT TNebulaPlusPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TNebulaPlusPhysics.h TNebulaPlusPhysicsDict.cxx TNebulaPlusPhysics.rootmap libNPNebulaPlus.dylib DEPENDS TNebulaPlusPhysics.h)
add_custom_command(OUTPUT TNebulaPlusDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TNebulaPlusData.h TNebulaPlusDataDict.cxx TNebulaPlusData.rootmap libNPNebulaPlus.dylib DEPENDS TNebulaPlusData.h)
add_library(NPNebulaPlus SHARED TNebulaPlusSpectra.cxx TNebulaPlusData.cxx TNebulaPlusPhysics.cxx TNebulaPlusDataDict.cxx TNebulaPlusPhysicsDict.cxx )
target_link_libraries(NPNebulaPlus ${ROOT_LIBRARIES} NPCore) 
install(FILES TNebulaPlusData.h TNebulaPlusPhysics.h TNebulaPlusSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

