add_custom_command(OUTPUT TActarPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TActarPhysics.h TActarPhysicsDict.cxx TActarPhysics.rootmap libNPActar.dylib DEPENDS TActarPhysics.h)
add_custom_command(OUTPUT TActarDataDict.cxx COMMAND ../../scripts/build_dict.sh TActarData.h TActarDataDict.cxx TActarData.rootmap libNPActar.dylib DEPENDS TActarData.h)
add_custom_command(OUTPUT MEventReducedDict.cxx COMMAND ../../scripts/build_dict.sh MEventReduced.h MEventReducedDict.cxx MEventReduced.rootmap libNPActar.dylib MTreeStructureLinkDef.h DEPENDS MEventReduced.h)
add_custom_command(OUTPUT TActarPhysicsReaderDict.cxx COMMAND ../../scripts/build_dict.sh TActarPhysicsReader.h TActarPhysicsReaderDict.cxx TActarPhysicsReader.rootmap libNPActar.dylib DEPENDS TActarPhysicsReader.h)
add_custom_command(OUTPUT TActarScatteringDict.cxx COMMAND ../../scripts/build_dict.sh TActarScattering.h TActarScatteringDict.cxx TActarScattering.rootmap libNPActar.dylib DEPENDS TActarScattering.h)
add_custom_command(OUTPUT TActarBeamDict.cxx COMMAND ../../scripts/build_dict.sh TActarBeam.h TActarBeamDict.cxx TActarBeam.rootmap libNPActar.dylib DEPENDS TActarBeam.h)

add_library(NPActar SHARED TActarSpectra.cxx TActarData.cxx TActarPhysics.cxx TActarPhysicsReader.cxx TActarScatteringDict.cxx TActarScattering.cxx TActarBeam.cxx TActarBeamDict.cxx TActarDataDict.cxx TActarPhysicsDict.cxx TActarPhysicsReaderDict.cxx MEventReducedDict.cxx)
target_link_libraries(NPActar ${ROOT_LIBRARIES} NPCore NPTrackReconstruction) 
install(FILES TActarData.h TActarPhysics.h TActarSpectra.h MEventReduced.h TActarPhysicsReader.h TActarScattering.h TActarBeam.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

