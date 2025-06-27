add_library(NPSNeuLAND SHARED  NeuLAND.cc)
target_link_libraries(NPSNeuLAND NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPNeuLAND)
