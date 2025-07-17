add_library(NPSFoil SHARED  Foil.cc)
target_link_libraries(NPSFoil NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPFoil)
