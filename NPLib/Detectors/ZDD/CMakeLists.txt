add_custom_command(OUTPUT TZDDPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TZDDPhysics.h TZDDPhysicsDict.cxx TZDDPhysics.rootmap libNPZDD.dylib DEPENDS TZDDPhysics.h)
add_custom_command(OUTPUT TZDDDataDict.cxx COMMAND ../../scripts/build_dict.sh TZDDData.h TZDDDataDict.cxx TZDDData.rootmap libNPZDD.dylib DEPENDS TZDDData.h)
add_custom_command(OUTPUT TZDDPhysicsReaderDict.cxx COMMAND ${CMAKE_BINARY_DIR}/scripts/build_dict.sh TZDDPhysicsReader.h TZDDPhysicsReaderDict.cxx TZDDPhysicsReader.rootmap libNPZDD.dylib DEPENDS TZDDPhysicsReader.h)
add_custom_command(OUTPUT TZDDCalDataDict.cxx COMMAND ../../scripts/build_dict.sh TZDDCalData.h TZDDCalDataDict.cxx TZDDCalData.rootmap libNPZDD.dylib DEPENDS TZDDCalData.h)
add_library(NPZDD SHARED TZDDSpectra.cxx TZDDData.cxx TZDDCalData.cxx TZDDPhysics.cxx TZDDPhysicsReader.cxx TZDDCalDataDict.cxx TZDDDataDict.cxx TZDDPhysicsDict.cxx TZDDPhysicsReaderDict.cxx)
target_link_libraries(NPZDD ${ROOT_LIBRARIES} -lTreePlayer -lSpectrum NPCore NPPhysics) 
install(FILES TZDDData.h TZDDCalData.h TZDDPhysics.h TZDDSpectra.h TZDDPhysicsReader.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

