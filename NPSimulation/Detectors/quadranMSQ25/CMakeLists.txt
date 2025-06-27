add_library(NPSquadranMSQ25 SHARED  quadranMSQ25.cc)
target_link_libraries(NPSquadranMSQ25 NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPquadranMSQ25)
