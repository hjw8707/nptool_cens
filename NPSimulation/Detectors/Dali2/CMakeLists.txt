add_library(NPSDali2 SHARED  Dali2.cc)
target_link_libraries(NPSDali2 NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPDali2)
