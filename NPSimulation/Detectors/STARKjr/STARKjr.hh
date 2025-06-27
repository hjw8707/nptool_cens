#ifndef STARKjr_h
#define STARKjr_h 1
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describes Tina simulation                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include <string>
#include <vector>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "NPSVDetector.hh"
#include "NPInputParser.h"
#include "TSTARKjrData.h"
#include "TSTARKjrRaw.h"

////////////////////////////////////////////////////////////////////////////////
// namespace for STARKjr
//  for some constants
////////////////////////////////////////////////////////////////////////////////
namespace STARKjrNS
{
  const G4double EnergyThreshold = 0.1*MeV;
  ////////////////////////////////////////////////////////////
  // Resolution
  ////////////////////////////////////////////////////////////
  const G4double X6_TRes = 0.213;
  const G4double X6_ERes = 0.015;
  const G4double BB10_TRes = 0.213;
  const G4double BB10_ERes = 0.015;
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Geometry
  ////////////////////////////////////////////////////////////
  //
  // X6
  ////////////////////////////////////////////////////////////
  const G4double X6_PCBX    = 45.20*mm; 
  const G4double X6_PCBY    = 93.10*mm;
  const G4double X6_PCBZ    =  2.40*mm;
  const G4double X6_PCBSub1X = 43.60*mm;
  const G4double X6_PCBSub1Y = 78.30*mm;
  const G4double X6_PCBSub1Z =  1.20*mm;
  const G4double X6_PCBSub1XOffset = 0.0*mm;
  const G4double X6_PCBSub1YOffset = 6.2*mm;
  const G4double X6_PCBSub1ZOffset = 0.6*mm;
  const G4double X6_PCBSub2X = 42.20*mm;
  const G4double X6_PCBSub2Y = 76.90*mm;
  const G4double X6_PCBSub2Z =  1.20*mm;
  const G4double X6_PCBSub2XOffset = 0.0*mm;
  const G4double X6_PCBSub2YOffset = 6.2*mm;
  const G4double X6_PCBSub2ZOffset = -0.6*mm;

  const G4double X6_SiX = 43.30*mm;
  const G4double X6_SiY = 78.00*mm;
  const G4double X6_SiZ =  1.00*mm; 
  const G4double X6_SiXOffset =  0.0*mm;
  const G4double X6_SiYOffset =  6.2*mm;
  const G4double X6_SiZOffset =  0.5*mm; 
  const G4double X6_SiActiveX = 40.30*mm; 
  const G4double X6_SiActiveY = 75.00*mm;
  const G4double X6_SiActiveZ =  1.00*mm; // 1000 um
  
  const G4int X6_NFrontStrips = 8;
  const G4int X6_NBackStrips  = 4;
  ////////////////////////////////////////////////////////////
  //
  // BB10
  ////////////////////////////////////////////////////////////
  const G4double BB10_PCBX    = 45.20*mm; 
  const G4double BB10_PCBY    = 93.10*mm;
  const G4double BB10_PCBZ    =  2.40*mm;
  const G4double BB10_PCBSub1X = 43.60*mm;
  const G4double BB10_PCBSub1Y = 78.30*mm;
  const G4double BB10_PCBSub1Z =  1.20*mm;
  const G4double BB10_PCBSub1XOffset = 0.0*mm;
  const G4double BB10_PCBSub1YOffset = 6.5*mm; // only different to X6 PCB
  const G4double BB10_PCBSub1ZOffset = 0.6*mm;
  const G4double BB10_PCBSub2X = 42.20*mm;
  const G4double BB10_PCBSub2Y = 76.90*mm;
  const G4double BB10_PCBSub2Z =  1.20*mm;
  const G4double BB10_PCBSub2XOffset = 0.0*mm;
  const G4double BB10_PCBSub2YOffset = 6.5*mm; // only different to X6 PCB
  const G4double BB10_PCBSub2ZOffset = -0.6*mm;

  const G4double BB10_SiX = 43.30*mm;
  const G4double BB10_SiY = 78.00*mm;
  const G4double BB10_SiZ =  0.14*mm; 
  const G4double BB10_SiXOffset =  0.0*mm;
  const G4double BB10_SiYOffset =  6.5*mm;
  const G4double BB10_SiZOffset =  0.07*mm; 
  const G4double BB10_SiActiveX = 39.45*mm; 
  const G4double BB10_SiActiveY = 74.15*mm;
  const G4double BB10_SiActiveZ =  0.14*mm; // 140 um
  
  const G4int BB10_NFrontStrips = 8;
  const G4int BB10_NBackStrips  = 1;
  ////////////////////////////////////////////////////////////
  //
  // QQQ5
  const G4double QQQ5_PCBOutR = 86*mm;
  const G4double QQQ5_PCBInR  = 15*mm;
  const G4double QQQ5_PCBPhi0 =  0 * deg; // Starting point
  const G4double QQQ5_PCBPhi1 = 90 * deg; // ANGLE
  const G4double QQQ5_PCBT    = 3.4*mm;
  const G4double QQQ5_PCBCutX = 3.4*2*mm; // 3.4 mm gap from the arc center?
  const G4double QQQ5_PCBCutY = 86*mm;
  const G4double QQQ5_PCBCutZ =  4*mm;
  const G4double QQQ5_PCBCutXOffset =  0*mm;
  const G4double QQQ5_PCBCutYOffset = 43*mm;
  const G4double QQQ5_PCBCutZOffset =  0*mm;

  // QQQ Wafer
  const G4double QQQ5_SiOutR = 84.0*mm;
  const G4double QQQ5_SiInR  = 23.2*mm;
  const G4double QQQ5_SiT    = 1*mm;
  const G4double QQQ5_SiPhi0 =  0 * deg;
  const G4double QQQ5_SiPhi1 = 90 * deg;
  const G4double QQQ5_SiActiveOutR = 81.95*mm;
  const G4double QQQ5_SiActiveInR  = 25.25*mm;
  const G4double QQQ5_SiCut1X = (3.4+0.68)*2*mm;
  const G4double QQQ5_SiCut1Y = QQQ5_SiOutR;
  const G4double QQQ5_SiCut1Z =   2*mm;
  const G4double QQQ5_SiCut1XOffset =  0*mm;
  const G4double QQQ5_SiCut1YOffset = QQQ5_SiCut1Y/2;
  const G4double QQQ5_SiCut1ZOffset =  0*mm;
  const G4double QQQ5_SiCut2X = QQQ5_SiOutR;
  const G4double QQQ5_SiCut2Y = 0.92*2*mm;
  const G4double QQQ5_SiCut2Z =   2*mm;
  const G4double QQQ5_SiCut2XOffset = QQQ5_SiCut2X/2;
  const G4double QQQ5_SiCut2YOffset =  0*mm;
  const G4double QQQ5_SiCut2ZOffset =  0*mm;

    
  const G4int    QQQ5_NRStrip = 32 ;
  const G4int    QQQ5_NAStrip = 4 ;

  ////////////////////////////////////////////////////////////
  //
  // Connector
  const G4double Conn_X = 40.0*mm;
  const G4double Conn_Y =  5.0*mm;
  const G4double Conn_Z =  5.0*mm;
  ////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////



class STARKjr: public NPS::VDetector {
    
public:
  STARKjr();
  virtual ~STARKjr();
  
  void AddDetector(string Type, G4ThreeVector POS, int Flip, double Beta);
  G4AssemblyVolume* BuildX6Detector();
  G4AssemblyVolume* BuildBB10Detector();
  G4AssemblyVolume* BuildQQQ5Detector();
  
  // Inherited from NPS::VDetector class /////////////
public:
  // Read stream at Configfile to pick-up parameters of detector
  // called in DetectorConstruction::ReadDetectorConfiguration
  void ReadConfiguration(NPL::InputParser);
  
  // Construct detector and initialise sensitive part
  // called after DetectorConstruction::AddDetector
  void ConstructDetector(G4LogicalVolume* world);
  
  // Add detector branch to the EventTree
  // called after DetectorConstruction::AddDetector
  void InitializeRootOutput();
  
  // Read sensitive part and fill the Root tree
  // called in EventAction::EndOfEventAvtion
  void ReadSensitive(const G4Event* event);
  
  // Initialise all scorers used by the detector
  void InitializeScorers();
  G4MultiFunctionalDetector* m_X6Det;
  G4MultiFunctionalDetector* m_BB10Det;
  G4MultiFunctionalDetector* m_QQQ5Det;
  ////////////////////////////////////////////////////
  
private:
  G4int HCID_X6;
  G4int HCID_BB10;
  G4int HCID_QQQ5;
  
  G4AssemblyVolume* m_X6;
  G4AssemblyVolume* m_BB10;
  G4AssemblyVolume* m_QQQ5;

  // Event class to store data
  TSTARKjrData* m_Event;
  TSTARKjrRaw* m_Raw;

  // Type & Geometry
  vector<string>        m_Type;
  vector<G4ThreeVector> m_Pos;
  vector<int>           m_Flip;
  vector<G4double>      m_Beta;

  // Visualisation
  G4VisAttributes *m_VisX6, *m_VisX6PCB;
  G4VisAttributes *m_VisBB10, *m_VisBB10PCB;
  G4VisAttributes *m_VisQQQ5, *m_VisQQQ5PCB;
  G4VisAttributes *m_VisConn; 
    
public:
    // Dynamic loading of the library
    static NPS::VDetector* Construct();  
};
#endif
