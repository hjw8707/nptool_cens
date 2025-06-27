#ifndef EventGeneratorAlphaDecay_h
#define EventGeneratorAlphaDecay_h
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This event Generator is used to simulated AlphaDecay ion Source           *
 *  Very usefull to figure out Geometric Efficacity of experimental Set-Up   *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>
#include <cmath>
using namespace std;

using namespace CLHEP;

// G4 headers
#include "G4Event.hh"

// NPS headers
#include "VEventGenerator.hh"
#include "ParticleStack.hh"
#include "NPInputParser.h"

// ROOT headers
#include "TString.h"
#include "TF1.h"
#include "TH1.h"

class EventGeneratorAlphaDecay : public NPS::VEventGenerator{
public:     // Constructor and destructor
    EventGeneratorAlphaDecay() ;
    virtual ~EventGeneratorAlphaDecay();

public:     // Inherit from VEventGenerator Class
    void ReadConfiguration(NPL::InputParser) ;
    void GenerateEvent(G4Event*) ;
    void InitializeRootOutput()  ;

private:    // Source parameter from input file
    G4double              m_EnergyLow        ;  // Lower limit of energy range
    G4double              m_EnergyHigh       ;  // Upper limit of energy range
    G4double              m_HalfOpenAngleMin ;  // Min Half open angle of the source
    G4double              m_HalfOpenAngleMax ;  // Max Half open angle of the source
    G4double              m_x0               ;  // Vertex Position X
    G4double              m_y0               ;  // Vertex Position Y
    G4double              m_z0               ;  // Vertex Position Z
    G4double              m_SigmaX           ;  // if m_SourceProfile==Flat, m_SigmaX=radius
    G4double              m_SigmaY           ;  // if m_SourceProfile==Flat, m_SigmaY=0
    G4double              m_SigmaZ           ;
    TString               m_direction        ;
    TString               m_SourceProfile    ;  // either "Gauss" (by default) or "Flat" on disk with radius of m_SigmaX
    G4ParticleDefinition* m_particle         ;  // Kind of particle to shoot isotropically
    G4double              m_ExcitationEnergy ;  // Excitation energy of the emitted particle
	  G4double              m_ActivityBq       ;  // if m_DecayLaw = "on" gives here the Activity in Bq of the sample
    G4double              m_TimeWindow       ;  // DAQ time window during we see the alpha pile-up
    
    ParticleStack*        m_ParticleStack    ;
};
#endif
