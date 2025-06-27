add_library(NPSSTARK SHARED  STARK.cc)
target_link_libraries(NPSSTARK NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPSTARK)
