add_library(NPSSTARKjr SHARED  STARKjr.cc)
target_link_libraries(NPSSTARKjr NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPSTARKjr)
