add_library(NPSMonster SHARED  Monster.cc)
target_link_libraries(NPSMonster NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPMonster)
