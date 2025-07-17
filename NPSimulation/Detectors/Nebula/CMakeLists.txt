add_library(NPSNebula SHARED  Nebula.cc)
target_link_libraries(NPSNebula NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPNebula)
