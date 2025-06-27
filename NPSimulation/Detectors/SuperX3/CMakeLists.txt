add_library(NPSSuperX3 SHARED  SuperX3.cc )
target_link_libraries(NPSSuperX3 NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPSuperX3)
