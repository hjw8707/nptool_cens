add_library(NPSHicari SHARED  Hicari.cc CConvexPolyhedron.cc)
target_link_libraries(NPSHicari NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPHicari)
