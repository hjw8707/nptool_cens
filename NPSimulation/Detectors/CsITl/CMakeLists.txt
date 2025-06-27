add_library(NPSCsITl SHARED  CsITl.cc)
target_link_libraries(NPSCsITl NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPCsITl)
