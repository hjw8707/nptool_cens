#ifndef GRAPE_h
#define GRAPE_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe the GRAPE Silicon detector                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>
#include <vector>

// G4 header defining G4 types
#include "globals.hh"

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4VSolid.hh"

// NPSimulation header
#include "NPSVDetector.hh"

// NPLib
#include "TGRAPEData.h"
#include "NPInputParser.h"
using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class GRAPE : public NPS::VDetector
{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
public:
  GRAPE() ;
   ~GRAPE() ;
  
  ////////////////////////////////////////////////////
  //////// Specific Function of this Class ///////////
  ////////////////////////////////////////////////////
public:
  // Add clover at the standard position of the array
  // Take as argument the standard clover Id.
  void AddGRAPEStandard(vector<int> GRAPEId);
  
  // Add clover at a free position in space with coordinate
  // in spherical coordinate
  // Beta are the three angles of rotation in the GRAPE frame 
  void AddGRAPEFreePosition(int GRAPEId,
                             double R,
                             double Theta,
                             double Phi,
                             double BetaX,
                             double BetaY,
                             double BetaZ);
  
  // Return a clover in the configuration given by option (not use a the moment)
  void ConstructGRAPE();
  
  // Return a modeling of the Crystal
  //  G4LogicalVolume* ConstructCrystal();
  
  // Return a modeling of the Capsule
  //  G4LogicalVolume* ConstructCapsule();
 
  // Return a modeling of the Dewar
  //  G4LogicalVolume* ConstructDewar();

  // Return a G4VSolid modeling the BGO
  //  G4LogicalVolume* ConstructBGO();
  
  ////////////////////////////////////////////////////
  /////////  Inherite from NPS::VDetector class ///////////
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
  
  
  ////////////////////////////////////////////////////
  ///////////Event class to store Data////////////////
  ////////////////////////////////////////////////////
private:
  TGRAPEData*    m_GRAPEData ;
  
  ////////////////////////////////////////////////////
  ///////////////// Scorer Related ///////////////////
  ////////////////////////////////////////////////////
  
private:
  // Geometry related
  G4LogicalVolume* m_LogicGRAPE;
  
  
  //   Initialize all Scorer
  void InitializeScorers() ;
  
  //   Scorer Associate to the Silicon
  G4MultiFunctionalDetector*   m_GRAPEScorer ;
  //G4MultiFunctionalDetector*   m_BGOScorer ;

private:
  //    Initialize material used in detector definition
  void InitializeMaterial();
  
  //   List of material
  G4Material* m_MaterialVacuum  ;
  G4Material* m_MaterialGe; 
  G4Material* m_MaterialAl;
  G4Material* m_MaterialCu;
  G4Material* m_MaterialN2;
  G4Material* m_MaterialC;
  ////////////////////////////////////////////////////
  ///////////////Private intern Data//////////////////
  ////////////////////////////////////////////////////
private:
  // GRAPE Position
  vector<int>    m_GRAPEId;
  vector<double> m_R;
  vector<double> m_Theta;
  vector<double> m_Phi;
  vector<double> m_BetaX;
  vector<double> m_BetaY;
  vector<double> m_BetaZ;
    
  // Frame: true if the frame is to be done
  //  bool m_LeftFrame;
  //  bool m_RightFrame;
  
private:/// Visualisation Attribute:
  G4VisAttributes* BlueVisAtt;
  G4VisAttributes* GreenVisAtt;
  G4VisAttributes* RedVisAtt;
  G4VisAttributes* WhiteVisAtt;
  G4VisAttributes* TrGreyVisAtt;
public:
    static NPS::VDetector* Construct();
};
#endif
