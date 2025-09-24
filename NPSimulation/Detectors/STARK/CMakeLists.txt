add_library(NPSSTARK SHARED  STARK.cc)
target_link_libraries(NPSSTARK NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPSTARK)
target_include_directories(NPSSTARK PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
