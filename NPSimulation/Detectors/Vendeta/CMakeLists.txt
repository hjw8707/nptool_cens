add_library(NPSVendeta SHARED  Vendeta.cc)
target_link_libraries(NPSVendeta NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPVendeta)
