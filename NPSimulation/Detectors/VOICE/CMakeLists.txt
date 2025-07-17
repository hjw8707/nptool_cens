add_library(NPSVOICE SHARED  VOICE.cc)
target_link_libraries(NPSVOICE NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPVOICE)
