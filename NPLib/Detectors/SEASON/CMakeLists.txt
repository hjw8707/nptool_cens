add_custom_command(OUTPUT TSEASONPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TSEASONPhysics.h TSEASONPhysicsDict.cxx TSEASONPhysics.rootmap libNPSEASON.dylib DEPENDS TSEASONPhysics.h)
add_custom_command(OUTPUT TSEASONDataDict.cxx COMMAND ../../scripts/build_dict.sh TSEASONData.h TSEASONDataDict.cxx TSEASONData.rootmap libNPSEASON.dylib DEPENDS TSEASONData.h)
add_library(NPSEASON SHARED TSEASONSpectra.cxx TSEASONData.cxx TSEASONPhysics.cxx TSEASONDataDict.cxx TSEASONPhysicsDict.cxx )
target_link_libraries(NPSEASON ${ROOT_LIBRARIES} NPCore) 
install(FILES TSEASONData.h TSEASONPhysics.h TSEASONSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

