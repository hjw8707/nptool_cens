add_library(NPSEpic SHARED  Epic.cc)
target_link_libraries(NPSEpic NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPEpic)
