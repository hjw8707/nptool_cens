add_custom_command(OUTPUT TquadranMSQ25PhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TquadranMSQ25Physics.h TquadranMSQ25PhysicsDict.cxx TquadranMSQ25Physics.rootmap libNPquadranMSQ25.dylib DEPENDS TquadranMSQ25Physics.h)
add_custom_command(OUTPUT TquadranMSQ25DataDict.cxx COMMAND ../../scripts/build_dict.sh TquadranMSQ25Data.h TquadranMSQ25DataDict.cxx TquadranMSQ25Data.rootmap libNPquadranMSQ25.dylib DEPENDS TquadranMSQ25Data.h)
add_library(NPquadranMSQ25 SHARED TquadranMSQ25Spectra.cxx TquadranMSQ25Data.cxx TquadranMSQ25Physics.cxx TquadranMSQ25DataDict.cxx TquadranMSQ25PhysicsDict.cxx )
target_link_libraries(NPquadranMSQ25 ${ROOT_LIBRARIES} NPCore) 
install(FILES TquadranMSQ25Spectra.h TquadranMSQ25Data.h TquadranMSQ25Physics.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

