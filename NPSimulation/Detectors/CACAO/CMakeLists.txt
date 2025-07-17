add_library(NPSCACAO SHARED  CACAO.cc)
target_link_libraries(NPSCACAO NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPCACAO)
