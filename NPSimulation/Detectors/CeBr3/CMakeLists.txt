add_library(NPSCeBr3 SHARED  CeBr3.cc)
target_link_libraries(NPSCeBr3 NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPCeBr3)
