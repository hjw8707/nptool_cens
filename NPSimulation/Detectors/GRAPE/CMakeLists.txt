add_library(NPSGRAPE SHARED  GRAPE.cc)
target_link_libraries(NPSGRAPE NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPGRAPE)
