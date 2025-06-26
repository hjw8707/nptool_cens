add_custom_command(OUTPUT TQQQ3PhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TQQQ3Physics.h TQQQ3PhysicsDict.cxx TQQQ3Physics.rootmap libNPQQQ3.dylib DEPENDS TQQQ3Physics.h)
add_custom_command(OUTPUT TQQQ3DataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TQQQ3Data.h TQQQ3DataDict.cxx TQQQ3Data.rootmap libNPQQQ3.dylib DEPENDS TQQQ3Data.h)
add_library(NPQQQ3 SHARED TQQQ3Data.cxx TQQQ3Physics.cxx TQQQ3DataDict.cxx TQQQ3PhysicsDict.cxx )
target_link_libraries(NPQQQ3 ${ROOT_LIBRARIES} NPCore) 
install(FILES TQQQ3Data.h TQQQ3Physics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

