add_custom_command(OUTPUT TPlastic_BEDOPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TPlastic_BEDOPhysics.h TPlastic_BEDOPhysicsDict.cxx TPlastic_BEDOPhysics.rootmap libNPPlastic_BEDO.dylib DEPENDS TPlastic_BEDOPhysics.h)
add_custom_command(OUTPUT TPlastic_BEDODataDict.cxx COMMAND ../../scripts/build_dict.sh TPlastic_BEDOData.h TPlastic_BEDODataDict.cxx TPlastic_BEDOData.rootmap libNPPlastic_BEDO.dylib DEPENDS TPlastic_BEDOData.h)
add_library(NPPlastic_BEDO SHARED TPlastic_BEDOSpectra.cxx TPlastic_BEDOData.cxx TPlastic_BEDOPhysics.cxx TPlastic_BEDODataDict.cxx TPlastic_BEDOPhysicsDict.cxx )
target_link_libraries(NPPlastic_BEDO ${ROOT_LIBRARIES} NPCore) 
install(FILES TPlastic_BEDOData.h TPlastic_BEDOPhysics.h TPlastic_BEDOSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

