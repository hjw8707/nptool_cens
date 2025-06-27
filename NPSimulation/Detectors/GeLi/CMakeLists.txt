add_library(NPSGeLi SHARED  GeLi.cc)
target_link_libraries(NPSGeLi NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPGeLi)
