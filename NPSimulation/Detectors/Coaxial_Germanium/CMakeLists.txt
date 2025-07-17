add_library(NPSCoaxial_Germanium SHARED  Coaxial_Germanium.cc)
target_link_libraries(NPSCoaxial_Germanium NPSCore ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} -lNPCoaxial_Germanium)
