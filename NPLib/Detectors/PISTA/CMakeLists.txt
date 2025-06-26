add_custom_command(OUTPUT TPISTAPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TPISTAPhysics.h TPISTAPhysicsDict.cxx TPISTAPhysics.rootmap libNPPISTA.dylib DEPENDS TPISTAPhysics.h)

add_custom_command(OUTPUT TPISTADataDict.cxx COMMAND ../../scripts/build_dict.sh TPISTAData.h TPISTADataDict.cxx TPISTAData.rootmap libNPPISTA.dylib DEPENDS TPISTAData.h)

add_custom_command(OUTPUT TFPMWDataDict.cxx COMMAND ../../scripts/build_dict.sh TFPMWData.h TFPMWDataDict.cxx TFPMWData.rootmap libNPPISTA.dylib DEPENDS TFPMWData.h)

add_custom_command(OUTPUT TFPMWPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TFPMWPhysics.h TFPMWPhysicsDict.cxx TFPMWPhysics.rootmap libNPPISTA.dylib DEPENDS TFPMWPhysics.h)

add_custom_command(OUTPUT TICDataDict.cxx COMMAND ../../scripts/build_dict.sh TICData.h TICDataDict.cxx TICData.rootmap libNPPISTA.dylib DEPENDS TICData.h)

add_custom_command(OUTPUT TICPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TICPhysics.h TICPhysicsDict.cxx TICPhysics.rootmap libNPPISTA.dylib DEPENDS TICPhysics.h)

add_custom_command(OUTPUT TVamosReconstructionDict.cxx COMMAND ../../scripts/build_dict.sh TVamosReconstruction.h TVamosReconstructionDict.cxx TVamosReconstruction.rootmap libNPPISTA.dylib DEPENDS TVamosReconstruction.h)

add_custom_command(OUTPUT TTimeDataDict.cxx COMMAND ../../scripts/build_dict.sh TTimeData.h TTimeDataDict.cxx TTimeData.rootmap libNPPISTA.dylib DEPENDS TTimeData.h)

add_library(NPPISTA SHARED TPISTASpectra.cxx TPISTAData.cxx TPISTAPhysics.cxx
    TPISTADataDict.cxx TPISTAPhysicsDict.cxx TFPMWData.cxx TFPMWDataDict.cxx
    TFPMWPhysics.cxx TFPMWPhysicsDict.cxx TICData.cxx TICDataDict.cxx
    TICPhysics.cxx TICPhysicsDict.cxx TVamosReconstruction.cxx TVamosReconstructionDict.cxx TTimeData.cxx TTimeDataDict.cxx)

target_link_libraries(NPPISTA ${ROOT_LIBRARIES} NPCore) 

install(FILES TPISTAData.h TPISTAPhysics.h TPISTASpectra.h TFPMWData.h TFPMWPhysics.h TICData.h TICPhysics.h
    TVamosReconstruction.h TTimeData.h TProfileEvaluator.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

