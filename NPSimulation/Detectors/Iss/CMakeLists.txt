# Check for gdml support
if(Geant4_gdml_FOUND)  
add_library(NPSIss SHARED  IssArray.cc)
target_link_libraries(NPSIss NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} NPIss)
endif()
