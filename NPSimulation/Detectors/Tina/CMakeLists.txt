add_library(NPSTina SHARED  Tina.cc)
target_link_libraries(NPSTina NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPTina)
