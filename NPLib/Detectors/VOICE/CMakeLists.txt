add_custom_command(OUTPUT TVOICEPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TVOICEPhysics.h TVOICEPhysicsDict.cxx TVOICEPhysics.rootmap libNPVOICE.dylib DEPENDS TVOICEPhysics.h)
add_custom_command(OUTPUT TVOICEDataDict.cxx COMMAND ../../scripts/build_dict.sh TVOICEData.h TVOICEDataDict.cxx TVOICEData.rootmap libNPVOICE.dylib DEPENDS TVOICEData.h)
add_library(NPVOICE SHARED TVOICESpectra.cxx TVOICEData.cxx TVOICEPhysics.cxx TVOICEDataDict.cxx TVOICEPhysicsDict.cxx )
target_link_libraries(NPVOICE ${ROOT_LIBRARIES} NPCore) 
install(FILES TVOICEData.h TVOICEPhysics.h TVOICESpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

