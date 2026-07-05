add_library(NPSCSSU SHARED CSSU.cc)
target_include_directories(NPSCSSU PRIVATE ${CMAKE_SOURCE_DIR}/../NPLib/Detectors/CSSU)
target_link_libraries(NPSCSSU NPSCore NPSScorers ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPCSSU)
