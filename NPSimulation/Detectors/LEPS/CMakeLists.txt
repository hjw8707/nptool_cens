add_library(NPSLEPS SHARED  LEPS.cc)
target_link_libraries(NPSLEPS NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPLEPS)
