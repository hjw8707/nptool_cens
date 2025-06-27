add_library(NPSAnnularS4 SHARED  AnnularS4.cc)
target_link_libraries(NPSAnnularS4 NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPAnnularS4)
