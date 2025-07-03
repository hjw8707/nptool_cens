add_library(NPSPlunger SHARED  Plunger.cc)
target_link_libraries(NPSPlunger NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPPlunger)
