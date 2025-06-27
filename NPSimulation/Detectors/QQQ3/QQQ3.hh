#ifndef QQQ3_h
#define QQQ3_h 1
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
 *  This class describe the QQQ3 Silicon detector                           *
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
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"

// NPSimulation header
#include "NPSVDetector.hh"

// NPLib
#include "NPInputParser.h"
#include "TQQQ3Data.h"
using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SHARC {
  // Energy and time Resolution
  const G4double ResoTime = 0;
  const G4double ResoEnergy = 0.035 * MeV; // = 82.25 keV of Resolution   //   Unit is MeV/2.35
  const G4double EnergyThreshold = 0.1 * MeV;
  // Geometry
  // QQQ //
  // QQQ PCB
  const G4double QQQ_PCB_Outer_Radius = 115.0 * mm;
  const G4double QQQ_PCB_Inner_Radius = 48.0 * mm;
  const G4double QQQ_PCB_Thickness = 2.4 * mm;

  // QQQ Wafer
  const G4double QQQ_Wafer_Outer_Radius = 99.0 * mm;
  const G4double QQQ_Wafer_Inner_Radius = 50.1 * mm;
  const G4double QQQ_Wafer_Starting_Phi = 3 * deg;
  const G4double QQQ_Wafer_Stopping_Phi = 87 * deg;
  const G4int QQQ_Wafer_NumberOf_RadialStrip = 16;
  const G4int QQQ_Wafer_NumberOf_AnnularStrip = 16;

} // namespace SHARC

using namespace SHARC;
class QQQ3 : public NPS::VDetector {
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
 public:
  QQQ3();
  ~QQQ3();

  ////////////////////////////////////////////////////
  //////// Specific Function of this Class ///////////
  ////////////////////////////////////////////////////
 public:
  // To add a Quadrant detector
  void AddDetector(G4ThreeVector Pos, G4double Thickness);
  void AddDetector(double R, double Theta, double Phi, double beta_u, double beta_v, double beta_w,
                   G4ThreeVector Offset, double Thickness);

  ////////////////////////////////////////////////////
  /////////  Inherite from NPS::VDetector class ///////////
  ////////////////////////////////////////////////////
 public:
  // Read stream at Configfile to pick-up parameters of detector (Position,...)
  // Called in DetecorConstruction::ReadDetextorConfiguration Method
  void ReadConfiguration(NPL::InputParser);

  // Construct detector and inialise sensitive part.
  // Called After DetecorConstruction::AddDetector Method
  void ConstructDetector(G4LogicalVolume* world);

  // Add Detector branch to the EventTree.
  // Called After DetecorConstruction::AddDetector Method
  void InitializeRootOutput();

  // Read sensitive part and fill the Root tree.
  // Called at in the EventAction::EndOfEventAvtion
  void ReadSensitive(const G4Event* event);

  ////////////////////////////////////////////////////
  ///////////Event class to store Data////////////////
  ////////////////////////////////////////////////////
 private:
  TQQQ3Data* m_Event;

  ////////////////////////////////////////////////////
  ///////////////// Scorer Related ///////////////////
  ////////////////////////////////////////////////////

 private:
  //   Initialize all Scorer
  void InitializeScorers();

  //   Scorer Associate to the Silicon
  G4MultiFunctionalDetector* m_QQQScorer;

 private:
  //    Initialize material used in detector definition
  void InitializeMaterial();

  //   List of material
  G4Material* m_MaterialSilicon;
  G4Material* m_MaterialVacuum;
  G4Material* m_MaterialPCB;

  ////////////////////////////////////////////////////
  ///////////////Private intern Data//////////////////
  ////////////////////////////////////////////////////
 private:
  // True if the detector is a Box, false if a quadrant
  vector<bool> m_Type;

  // Used for Quadrant detectors
  vector<G4ThreeVector> m_Pos, m_Offset; // R , Phi , Z
  vector<G4double> m_ThicknessQQQ;
  vector<double> m_R, m_Theta, m_Phi, m_beta_u, m_beta_v, m_beta_w;

  // Used for Box detectors
  vector<G4double> m_Z;
  vector<vector<G4double>> m_ThicknessPAD;

 private: /// Visualisation Attribute:
          // Dark Grey
  G4VisAttributes* SiliconVisAtt;
  // Green
  G4VisAttributes* PCBVisAtt;
  // Gold Yellow
  G4VisAttributes* PADVisAtt;
  // Light Grey
  G4VisAttributes* FrameVisAtt;

 public:
  static NPS::VDetector* Construct();
};
#endif
