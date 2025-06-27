#ifndef Dali2_h
#define Dali2_h 1
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: elidiano.tronchin@studenti.unipd.it                        *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Dali2 simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
#include <cmath>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TDali2Data.h"
#include "NPInputParser.h"

class Dali2 : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
public:
  Dali2() ;
  virtual ~Dali2() ;
  
  ////////////////////////////////////////////////////
  /////// Specific Function of this Class ///////////
  ////////////////////////////////////////////////////
public:
  // Cartesian
  void AddDetector(G4ThreeVector POS, int Type);
  void AddDetector(G4ThreeVector POS, G4ThreeVector ANG, int Type);
  // Spherical
  //  void AddDetector(double R,double Theta,double Phi, int Type);  
  //Cylindrical
  void AddDetector(double R,double Alpha,double Zeta, int Type);
  
  //  G4LogicalVolume* BuildSquareDetector();
  //  G4LogicalVolume* BuildCylindricalDetector();

  G4LogicalVolume* BuildDali2Detector(int type);

private:
  
  G4LogicalVolume* m_Dali2SGDetector;
  G4LogicalVolume* m_Dali2SCDetector;
  G4LogicalVolume* m_Dali1Detector;

  G4Material* NaI_Tl;
  G4Material* Air;
  G4Material* Al;
  G4Material* MgO;
  G4Material* Vacuum;
  G4Material* muMetal;
  G4Material* BoroSili_Glass;
  
  ////////////////////////////////////////////////////
  //////  Inherite from NPS::VDetector class /////////
  ////////////////////////////////////////////////////
public:
  // Read stream at Configfile to pick-up parameters of detector (Position,...)
  // Called in DetecorConstruction::ReadDetextorConfiguration Method
  void ReadConfiguration(NPL::InputParser) ;
  
  // Definition of materials
  // Called in ConstructDetector()
  void DefinitionMaterials();
  
  // Construct detector and inialise sensitive part.
  // Called After DetecorConstruction::AddDetector Method
  void ConstructDetector(G4LogicalVolume* world) ;
  
  // Add Detector branch to the EventTree.
  // Called After DetecorConstruction::AddDetector Method
  void InitializeRootOutput() ;
  
  // Read sensitive part and fill the Root tree.
  // Called at in the EventAction::EndOfEventAvtion
  void ReadSensitive(const G4Event* ) ;
  
public:   // Scorer
  //   Initialize all Scorer used by the MUST2Array
  void InitializeScorers() ;
  
  //   Associated Scorer
  G4MultiFunctionalDetector* m_DaliScorer ;
  ////////////////////////////////////////////////////
  ///////////Event class to store Data////////////////
  ////////////////////////////////////////////////////
private:
  TDali2Data* m_Event ;
  
  ////////////////////////////////////////////////////
  ///////////////Private intern Data//////////////////
  ////////////////////////////////////////////////////
private: // Geometry
  // Detector Coordinate 
  vector<double>  m_x;
  vector<double>  m_y; 
  vector<double>  m_z; 
  
  vector<double> m_psi;
  vector<double> m_theta;
  vector<double> m_phi;
  vector<int>  m_Type; 
  
  // Visualisation Attribute
  G4VisAttributes* m_VisSquare;
  
  // Needed for dynamic loading of the library
public:
  static NPS::VDetector* Construct();
};
#endif
