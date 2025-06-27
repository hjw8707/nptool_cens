add_library(NPSVamos SHARED  Vamos.cc)
target_link_libraries(NPSVamos NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPVamos)
