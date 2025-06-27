#ifndef CACAO_h
#define CACAO_h 1
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jwhwang@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 4월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  CACAO simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TCACAOData.h"
#include "NPInputParser.h"

class CACAO : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
public:
  CACAO() ;
  virtual ~CACAO() ;

  ////////////////////////////////////////////////////
  /////// Specific Function of this Class ///////////
  ////////////////////////////////////////////////////
public:
  void AddDetector(G4ThreeVector Pos, G4RotationMatrix Rot,
		   G4ThreeVector Dim,
		   G4double ShieldThickness);  
  
  G4LogicalVolume* BuildDetector(G4int i);

  void DefineMaterials();
  void ConstructChamber(G4LogicalVolume* world); // should be called in 'ConstructDetector'

  ////////////////////////////////////////////////////
  //////  Inherite from NPS::VDetector class /////////
  ////////////////////////////////////////////////////
public:
  // Read stream at Configfile to pick-up parameters of detector (Position,...)
  // Called in DetecorConstruction::ReadDetextorConfiguration Method
  void ReadConfiguration(NPL::InputParser) ;

  // Construct detector and inialise sensitive part.
  // Called After DetecorConstruction::AddDetector Method
  void ConstructDetector(G4LogicalVolume* world) ;

  // Add Detector branch to the EventTree.
  // Called After DetecorConstruction::AddDetector Method
  void InitializeRootOutput() ;

  // Read sensitive part and fill the Root tree.
  // Called at in the EventAction::EndOfEventAvtion
  void ReadSensitive(const G4Event* event) ;

public:   // Scorer
  //   Initialize all Scorer used by the MUST2Array
  void InitializeScorers() ;

  //   Associated Scorer
  G4MultiFunctionalDetector* m_CACAOScorer ;
  ////////////////////////////////////////////////////
  ///////////Event class to store Data////////////////
  ////////////////////////////////////////////////////
private:
  TCACAOData* m_Event ;

  ////////////////////////////////////////////////////
  ///////////////Private intern Data//////////////////
  ////////////////////////////////////////////////////
  G4Material* m_matScint;
  G4Material* m_matShield;
  G4Material* m_matPCB;
  G4Material* m_matChamber;

private: // Geometry
  // Detector Coordinate 
  vector<G4ThreeVector> m_Pos;     // Detector Position
  vector<G4RotationMatrix> m_Rot;  // Detector Rotation
  vector<G4ThreeVector> m_Dim;     // Detector Dimension

  vector<double> m_ShieldThickness; 

  ////////////////////////////////////////////////////////////////////////////////
  // CACAO chamber
  G4bool m_ChamberFound;
  G4double m_ChamberTRmin, m_ChamberTRmax, m_ChamberTZ0, m_ChamberTZ1;
  G4double m_ChamberCRmin, m_ChamberCRmax, m_ChamberTZ2, m_ChamberTZ3;
  G4double m_ChamberLRmin, m_ChamberLRmax, m_ChamberLH;
  G4double m_ChamberBRmin, m_ChamberBRmax, m_ChamberBZ1, m_ChamberBZ2;
  G4double m_ChamberBZ3;
  ////////////////////////////////////////////////////////////////////////////////

  // Visualisation Attribute
  G4VisAttributes* m_VisScint;
  G4VisAttributes* m_VisPCB;

  // Needed for dynamic loading of the library
public:
  static NPS::VDetector* Construct();
};
#endif
