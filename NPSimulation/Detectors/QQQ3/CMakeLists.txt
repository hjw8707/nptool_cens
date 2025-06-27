add_library(NPSQQQ3 SHARED  QQQ3.cc)
target_link_libraries(NPSQQQ3 NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPQQQ3)
