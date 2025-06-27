add_library(NPSPlastic_BEDO SHARED  Plastic_BEDO.cc)
target_link_libraries(NPSPlastic_BEDO NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPPlastic_BEDO)
