add_custom_command(OUTPUT TCsITlPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCsITlPhysics.h TCsITlPhysicsDict.cxx TCsITlPhysics.rootmap libNPCsITl.dylib DEPENDS TCsITlPhysics.h)
add_custom_command(OUTPUT TCsITlDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TCsITlData.h TCsITlDataDict.cxx TCsITlData.rootmap libNPCsITl.dylib DEPENDS TCsITlData.h)
add_library(NPCsITl SHARED TCsITlData.cxx TCsITlPhysics.cxx TCsITlDataDict.cxx TCsITlPhysicsDict.cxx )
target_link_libraries(NPCsITl ${ROOT_LIBRARIES} NPCore) 
install(FILES TCsITlData.h TCsITlPhysics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

