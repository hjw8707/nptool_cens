add_library(NPSEdinburghDSSD SHARED  EdinburghDSSD.cc)
target_link_libraries(NPSEdinburghDSSD NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPEdinburghDSSD)
