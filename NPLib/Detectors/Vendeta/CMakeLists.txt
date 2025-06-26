add_custom_command(OUTPUT TVendetaPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TVendetaPhysics.h TVendetaPhysicsDict.cxx TVendetaPhysics.rootmap libNPVendeta.dylib DEPENDS TVendetaPhysics.h)
add_custom_command(OUTPUT TVendetaDataDict.cxx COMMAND ../../scripts/build_dict.sh TVendetaData.h TVendetaDataDict.cxx TVendetaData.rootmap libNPVendeta.dylib DEPENDS TVendetaData.h)
add_library(NPVendeta SHARED TVendetaSpectra.cxx TVendetaData.cxx TVendetaPhysics.cxx TVendetaDataDict.cxx TVendetaPhysicsDict.cxx )
target_link_libraries(NPVendeta ${ROOT_LIBRARIES} NPCore) 
install(FILES TVendetaData.h TVendetaPhysics.h TVendetaSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

