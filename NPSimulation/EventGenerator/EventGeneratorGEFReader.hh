#ifndef EventGeneratorGEFReader_h
#define EventGeneratorGEFReader_h
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
 *  This event Generator is used to simulated Isotropic ion Source           *
 *  Very usefull to figure out Geometric Efficacity of experimental Set-Up   *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>
#include <cmath>
#include <fstream>
using namespace std;

using namespace CLHEP;

// G4 headers
#include "G4Event.hh"

// NPS headers
#include "VEventGenerator.hh"
#include "ParticleStack.hh"
#include "Particle.hh"

// NPL headers
#include "NPReaction.h"
#include "NPInputParser.h"
#include "NPParticle.h"
#include "TFissionConditions.h"

// ROOT headers
#include "TString.h"
#include "TF1.h"
#include "TH1.h"

class EventGeneratorGEFReader : public NPS::VEventGenerator{
  public:     // Constructor and destructor
    EventGeneratorGEFReader() ;
    virtual ~EventGeneratorGEFReader();

  public:     // Inherit from VEventGenerator Class
    void ReadConfiguration(NPL::InputParser);
    void GenerateEvent(G4Event*) ;
    void InitializeRootOutput() ;
    TVector3 ShootParticle(double,double, TString);
    void GetBoostFromTwoBodyReaction(double Ex);

  private:    // Source parameter from input file
    G4int event_ID;
    struct SourceParameters {
      SourceParameters()                          ;
      G4double                 m_x0               ;  // Vertex Position X
      G4double                 m_y0               ;  // Vertex Position Y
      G4double                 m_z0               ;  // Vertex Position Z
      G4double                 m_SigmaX           ;
      G4double                 m_SigmaY           ;
      G4double                 m_SigmaZ           ;
      G4double                 m_Boost            ;
      TString                  m_direction        ;
      vector<string>           m_particleName     ;
      string                   m_FissioningSystemName ;
      string                   m_BeamProfile      ;
      double                   m_GEFversion       ;
    };
    vector<SourceParameters> m_Parameters       ;
    ParticleStack*           m_ParticleStack    ;

    ifstream fInputDataFile;
    bool HasInputDataFile;
    NPL::Particle* m_FissioningSystem;
    NPL::Reaction* m_TwoBodyReaction;
    bool m_isTwoBody;

    int version_shift;
    int version_ff_shift;
    vector<double> LastLine;

    vector<string> AllowedParticles;  
  private:
    TFissionConditions* m_FissionConditions;

  public:
      void AttachFissionConditions();
};
#endif

// ============================================================================================================================================================================================================================================== //
// FORMAT OF THE lmd FILE OUTPUT AS A FUNCTION OF THE GEF version
// ============================================================================================================================================================================================================================================== //
// For each event, a line details CN, FFl, FFh and possible pre-fission particle emission
// - 2023.11 :      Zsad Asad        Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ... 
//   2023.12 :      Zsad Asad        Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ... 
//   2023.21 :      Zsad Asad        Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ... 
//   2023.22 :      Zsad Asad        Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ... 
// - 2023.31 : Mode Zsad Asad        Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ... 
//   2023.32 : Mode Zsad Asad        Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ...  
// - 2023.33 : Mode Zsad Asad Qvalue Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ...  
//   2024.11 : Mode Zsad Asad Qvalue Z1 Z2 A1sci A2sci A1post A2post I1pre I2pre I1gs I2gs Elab1pre cos(theta1) phi1 Elab2pre cos(theta2) phi2 Eexc1 Eexc2 n1 n2 TKEpre TKEpost E@fission particle-list E1 cos(theta1) phi1 E2 cos(theta2) phi2 ...
//
// ATTENTION : OUTPUT FORMAT CHANGES:
//   @version 2023.33 fission mode appears in the first column; 
//   @version 2023.31 Q value appears on the 4th column 
// ============================================================================================================================================================================================================================================== //
// While running GEF, if you ask for neutron and/or gamma an additional lines will be added before the next event line
// ============================================================================================================================================================================================================================================== //
