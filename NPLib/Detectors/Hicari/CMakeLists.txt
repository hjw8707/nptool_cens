add_custom_command(OUTPUT THicariPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh THicariPhysics.h THicariPhysicsDict.cxx THicariPhysics.rootmap libNPHicari.dylib DEPENDS THicariPhysics.h)
add_custom_command(OUTPUT THicariDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh THicariData.h THicariDataDict.cxx THicariData.rootmap libNPHicari.dylib DEPENDS THicariData.h)
add_library(NPHicari SHARED THicariData.cxx THicariPhysics.cxx THicariDataDict.cxx THicariPhysicsDict.cxx THicariSpectra.cxx)
target_link_libraries(NPHicari ${ROOT_LIBRARIES} NPCore) 
install(FILES THicariData.h THicariPhysics.h THicariSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

