add_library(NPSZDD SHARED  ZDD.cc)
target_link_libraries(NPSZDD NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPZDD)
