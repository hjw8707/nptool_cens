add_custom_command(OUTPUT TExogamPhysicsDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TExogamPhysics.h TExogamPhysicsDict.cxx TExogamPhysics.rootmap libNPExogam.dylib DEPENDS TExogamPhysics.h)
add_custom_command(OUTPUT TExogamDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TExogamData.h TExogamDataDict.cxx TExogamData.rootmap libNPExogam.dylib DEPENDS TExogamData.h)
add_custom_command(OUTPUT TExogamCalDataDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TExogamCalData.h TExogamCalDataDict.cxx TExogamCalData.rootmap libNPExogam.dylib DEPENDS TExogamCalData.h)
add_custom_command(OUTPUT TExogamPhysicsReaderDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TExogamPhysicsReader.h TExogamPhysicsReaderDict.cxx TExogamPhysicsReader.rootmap libNPExogam.dylib DEPENDS TExogamPhysicsReader.h)
add_library(NPExogam SHARED TExogamCalData.cxx TExogamData.cxx TExogamPhysics.cxx TExogamPhysicsReader.cxx TExogamGeo.cxx TExogamCalDataDict.cxx TExogamDataDict.cxx TExogamPhysicsDict.cxx TExogamPhysicsReaderDict.cxx TExogamSpectra.cxx)
if(CUBIX)
    target_link_libraries(NPExogam ${ROOT_LIBRARIES} ${CUBIX_LIBRARIES} -lTreePlayer -lSpectrum NPCore NPPhysics) 
    else()
    target_link_libraries(NPExogam ${ROOT_LIBRARIES} -lTreePlayer -lSpectrum NPCore NPPhysics) 
endif()
install(FILES Geometry_Clover_Exogam.h TExogamCalData.h TExogamData.h TExogamPhysics.h TExogamSpectra.h TExogamPhysicsReader.h TExogamGeo.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

