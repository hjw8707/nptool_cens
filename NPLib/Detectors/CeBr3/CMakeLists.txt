add_custom_command(OUTPUT TCeBr3PhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TCeBr3Physics.h TCeBr3PhysicsDict.cxx TCeBr3Physics.rootmap libNPCeBr3.dylib DEPENDS TCeBr3Physics.h)
add_custom_command(OUTPUT TCeBr3DataDict.cxx COMMAND ../../scripts/build_dict.sh TCeBr3Data.h TCeBr3DataDict.cxx TCeBr3Data.rootmap libNPCeBr3.dylib DEPENDS TCeBr3Data.h)
add_library(NPCeBr3 SHARED TCeBr3Spectra.cxx TCeBr3Data.cxx TCeBr3Physics.cxx TCeBr3DataDict.cxx TCeBr3PhysicsDict.cxx )
target_link_libraries(NPCeBr3 ${ROOT_LIBRARIES} NPCore) 
install(FILES TCeBr3Data.h TCeBr3Physics.h TCeBr3Spectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

