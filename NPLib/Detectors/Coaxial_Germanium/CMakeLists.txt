add_custom_command(OUTPUT TCoaxial_GermaniumPhysicsDict.cxx COMMAND ../../scripts/build_dict.sh TCoaxial_GermaniumPhysics.h TCoaxial_GermaniumPhysicsDict.cxx TCoaxial_GermaniumPhysics.rootmap libNPCoaxial_Germanium.dylib DEPENDS TCoaxial_GermaniumPhysics.h)
add_custom_command(OUTPUT TCoaxial_GermaniumDataDict.cxx COMMAND ../../scripts/build_dict.sh TCoaxial_GermaniumData.h TCoaxial_GermaniumDataDict.cxx TCoaxial_GermaniumData.rootmap libNPCoaxial_Germanium.dylib DEPENDS TCoaxial_GermaniumData.h)
add_library(NPCoaxial_Germanium SHARED TCoaxial_GermaniumSpectra.cxx TCoaxial_GermaniumData.cxx TCoaxial_GermaniumPhysics.cxx TCoaxial_GermaniumDataDict.cxx TCoaxial_GermaniumPhysicsDict.cxx )
target_link_libraries(NPCoaxial_Germanium ${ROOT_LIBRARIES} NPCore) 
install(FILES TCoaxial_GermaniumData.h TCoaxial_GermaniumPhysics.h TCoaxial_GermaniumSpectra.h DESTINATION ${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

