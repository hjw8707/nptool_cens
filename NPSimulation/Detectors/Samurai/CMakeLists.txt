add_library(NPSSamurai SHARED  Samurai.cc SamuraiFieldPropagation.cc SamuraiFDC2.cc SamuraiFDC1.cc SamuraiFDC0.cc SamuraiBDC.cc)
target_link_libraries(NPSSamurai NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPSamurai)
