add_library(NPSEventGenerator OBJECT EventGeneratorBeam.cc EventGeneratorMultipleParticle.cc EventGeneratorCosmic.cc EventGeneratorIsotropic.cc EventGeneratorGEFReader.cc EventGeneratorAlphaDecay.cc VEventGenerator.cc) 
#target_link_libraries(NPSEventGenerator ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${NPLib_LIBRARIES} -lNPInitialConditions -lNPInteractionCoordinates )
