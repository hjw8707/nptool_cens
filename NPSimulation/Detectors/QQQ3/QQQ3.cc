/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe the QQQ3 Silicon array                              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <limits>
#include <sstream>
// G4 Geometry object
#include "G4Box.hh"
#include "G4Tubs.hh"

// G4 sensitive
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4MaterialTable.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4UnionSolid.hh"
// NPS
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "QQQ3.hh"
// NPL
#include "NPOptionManager.h"

#include "RootOutput.h"
using namespace SHARC;

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// QQQ3 Specific Method
QQQ3::QQQ3() {
  InitializeMaterial();
  m_Event = new TQQQ3Data();
  // Dark Grey
  SiliconVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  // Green
  PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2));
  // Light Grey
  FrameVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));

  m_QQQScorer = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
QQQ3::~QQQ3() {
  // delete m_QQQScorer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QQQ3::AddDetector(G4ThreeVector Pos, G4double Thickness) {
  m_Type.push_back(true);
  m_Pos.push_back(Pos);
  m_ThicknessQQQ.push_back(Thickness);
  m_R.push_back(-1000);
  m_Theta.push_back(-1000);
  m_Phi.push_back(-1000);
  m_beta_u.push_back(-1000);
  m_beta_v.push_back(-1000);
  m_beta_w.push_back(-1000);
  m_Offset.push_back(G4ThreeVector(-1000, -1000, -1000));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QQQ3::AddDetector(double R, double Theta, double Phi, double beta_u, double beta_v, double beta_w,
                       G4ThreeVector Offset, double Thickness) {
  m_Type.push_back(false);
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_beta_u.push_back(beta_u);
  m_beta_v.push_back(beta_v);
  m_beta_w.push_back(beta_w);
  m_Offset.push_back(Offset);
  m_Pos.push_back(G4ThreeVector(-1000, -1000, -1000));
  m_ThicknessQQQ.push_back(Thickness);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void QQQ3::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("QQQ3");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> cyli = {"Z", "R", "Phi", "ThicknessDetector"};
  vector<string> sphe = {"R", "THETA", "PHI", "BETA", "ThicknessDetector"};

  for (unsigned int i = 0; i < blocks.size(); i++) {

    if (blocks[i]->HasTokenList(cyli)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  QQQ3 " << i + 1 << endl;
      double Z = blocks[i]->GetDouble("Z", "mm");
      double R = blocks[i]->GetDouble("R", "mm");
      double Phi = blocks[i]->GetDouble("Phi", "deg");
      double Thickness = blocks[i]->GetDouble("ThicknessDetector", "micrometer");
      AddDetector(G4ThreeVector(R, Phi, Z), Thickness);
    }
    else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  QQQ " << i + 1 << endl;
      double R = blocks[i]->GetDouble("R", "mm");
      double Theta = blocks[i]->GetDouble("THETA", "deg");
      double Phi = blocks[i]->GetDouble("PHI", "deg");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA", "deg");
      G4ThreeVector Offset = NPS::ConvertVector(blocks[i]->GetTVector3("Offset", "mm"));
      double Thickness = blocks[i]->GetDouble("ThicknessDetector", "micrometer");
      AddDetector(R, Theta, Phi, beta[0], beta[1], beta[2], Offset, Thickness);
    }
    else {
      cout << "Warning: check your input file formatting " << endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QQQ3::ConstructDetector(G4LogicalVolume* world) {
  // create the QQQ
  for (unsigned int i = 0; i < m_Type.size(); i++) {
    int DetNbr = i + 1;
    // Make the a single detector geometry
    G4Tubs* QQQDetector =
        new G4Tubs("QQQDetector", QQQ_PCB_Inner_Radius, QQQ_PCB_Outer_Radius, QQQ_PCB_Thickness * 0.5, 0., M_PI / 2.);

    G4Tubs* PCBFull =
        new G4Tubs("PCBFull", QQQ_PCB_Inner_Radius, QQQ_PCB_Outer_Radius, QQQ_PCB_Thickness * 0.5, 0., M_PI * 0.5);

    G4Tubs* WaferShape = new G4Tubs("WaferShape", QQQ_Wafer_Inner_Radius, QQQ_Wafer_Outer_Radius,
                                    QQQ_PCB_Thickness * 0.5 + 0.1 * mm, QQQ_Wafer_Starting_Phi, QQQ_Wafer_Stopping_Phi);

    G4Tubs* Wafer = new G4Tubs("Wafer", QQQ_Wafer_Inner_Radius, QQQ_Wafer_Outer_Radius, m_ThicknessQQQ[i] * 0.5,
                               QQQ_Wafer_Starting_Phi, QQQ_Wafer_Stopping_Phi);

    G4SubtractionSolid* PCB =
        new G4SubtractionSolid("PCB", PCBFull, WaferShape, new G4RotationMatrix, G4ThreeVector(0, 0, 0));

    // Master Volume
    G4LogicalVolume* logicQQQDetector = new G4LogicalVolume(QQQDetector, m_MaterialVacuum, "logicQQQDetector", 0, 0, 0);
    logicQQQDetector->SetVisAttributes(G4VisAttributes::GetInvisible());
    // Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(PCBVisAtt);

    // Sub Volume Wafer
    G4LogicalVolume* logicWafer = new G4LogicalVolume(Wafer, m_MaterialSilicon, "logicWafer", 0, 0, 0);
    logicWafer->SetVisAttributes(SiliconVisAtt);

    logicWafer->SetSensitiveDetector(m_QQQScorer);

    // Place the sub volume in the master volume
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "QQQ_PCB", logicQQQDetector,
                      false, DetNbr);

    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicWafer, "QQQ_Wafer", logicQQQDetector,
                      false, DetNbr);

    if (m_Type[i]) {
      // Place the masters volume in the world
      new G4PVPlacement(new G4RotationMatrix(0, 0, m_Pos[i].y()), G4ThreeVector(0, 0, m_Pos[i].z()), logicQQQDetector,
                        "QQQ", world, false, DetNbr);
    }
    else {
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      // Phi is angle between X axis and projection in (X,Y) plan
      // m_Theta[0] is angle between  position vector and z axis
      G4double wX = m_R[i] * sin(m_Theta[i] / rad) * cos(m_Phi[i] / rad);
      G4double wY = m_R[i] * sin(m_Theta[i] / rad) * sin(m_Phi[i] / rad);
      G4double wZ = m_R[i] * cos(m_Theta[i] / rad);
      auto Det_w = G4ThreeVector(wX, wY, wZ);
      // vector corresponding to the center of the module
      G4ThreeVector CT = Det_w;

      // vector parallel to one axis of silicon plane
      G4double ii = cos(m_Theta[i] / rad) * cos(m_Phi[i] / rad);
      G4double jj = cos(m_Theta[i] / rad) * sin(m_Phi[i] / rad);
      G4double kk = -sin(m_Theta[i] / rad);
      G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

      Det_w = Det_w.unit();
      auto Det_u = Det_w.cross(Y);
      auto Det_v = Det_w.cross(Det_u);
      Det_v = Det_v.unit();
      Det_u = Det_u.unit();
      // Passage Matrix from Lab Referential to Telescope Referential
      auto Det_rot = new G4RotationMatrix(Det_u, Det_v, Det_w);
      // Telescope is rotate of Beta angle around Det_v axis.
      Det_rot->rotate(m_beta_u[i], Det_u);
      Det_rot->rotate(m_beta_v[i], Det_v);
      Det_rot->rotate(m_beta_w[i], Det_w);
      // translation to place Telescope
      auto Det_pos = Det_w + CT + m_Offset[i];
      // Place the masters volume in the world
      new G4PVPlacement(Det_rot, Det_pos, logicQQQDetector, "QQQ", world, false, DetNbr);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void QQQ3::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("QQQ3")) {
    pTree->Branch("QQQ3", "TQQQ3Data", &m_Event);
  }
  pTree->SetBranchAddress("QQQ3", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void QQQ3::ReadSensitive(const G4Event*) {
  m_Event->Clear();
  ///////////
  // QQQ
  DSSDScorers::PS_Annular* QQQScorer = (DSSDScorers::PS_Annular*)m_QQQScorer->GetPrimitive(0);

  // Loop on the QQQ map
  unsigned int sizeRing = QQQScorer->GetRingMult();
  for (unsigned int i = 0; i < sizeRing; i++) {
    double Energy = QQQScorer->GetEnergyRing(i);
    if (Energy > EnergyThreshold) {
      double Time = QQQScorer->GetTimeRing(i);
      int DetNbr = QQQScorer->GetDetectorRing(i);
      int StripRing = QQQScorer->GetStripRing(i);
      m_Event->SetFront(DetNbr, QQQ_Wafer_NumberOf_AnnularStrip - StripRing + 1, RandGauss::shoot(Energy, ResoEnergy),
                        RandGauss::shoot(Time, ResoTime), RandGauss::shoot(Time, ResoTime));
    }
  }

  unsigned int sizeSector = QQQScorer->GetSectorMult();
  for (unsigned int i = 0; i < sizeSector; i++) {

    double Energy = QQQScorer->GetEnergySector(i);

    if (Energy > EnergyThreshold) {
      double Time = QQQScorer->GetTimeSector(i);
      int DetNbr = QQQScorer->GetDetectorSector(i);
      int StripSector = QQQScorer->GetStripSector(i);

      m_Event->SetBack(DetNbr, StripSector, RandGauss::shoot(Energy, ResoEnergy), RandGauss::shoot(Time, ResoTime),
                       RandGauss::shoot(Time, ResoTime));
    }
  }
  // clear map for next event
  QQQScorer->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void QQQ3::InitializeScorers() {
  //   Silicon Associate Scorer
  bool already_exist = false;
  m_QQQScorer = CheckScorer("QQQ3_Scorer", already_exist);
  // if the scorer were created previously nothing else need to be made
  if (already_exist)
    return;

  G4VPrimitiveScorer* QQQScorer = new DSSDScorers::PS_Annular(
      "QQQ3", 0, QQQ_Wafer_Inner_Radius, QQQ_Wafer_Outer_Radius, QQQ_Wafer_Starting_Phi, QQQ_Wafer_Stopping_Phi,
      QQQ_Wafer_NumberOf_AnnularStrip, QQQ_Wafer_NumberOf_RadialStrip, 1);

  G4VPrimitiveScorer* InterScorerQQQ = new InteractionScorers::PS_Interactions("QQQ3", ms_InterCoord, 0);

  // and register it to the multifunctionnal detector
  m_QQQScorer->RegisterPrimitive(QQQScorer);
  m_QQQScorer->RegisterPrimitive(InterScorerQQQ);

  G4SDManager::GetSDMpointer()->ListTree();
  //   Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_QQQScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
/////////////////Material Definition ///////////////////////////
////////////////////////////////////////////////////////////////
void QQQ3::InitializeMaterial() {
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* QQQ3::Construct() { return (NPS::VDetector*)new QQQ3(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_sharc {
 public:
  proxy_nps_sharc() {
    NPS::DetectorFactory::getInstance()->AddToken("QQQ3", "QQQ3");
    NPS::DetectorFactory::getInstance()->AddDetector("QQQ3", QQQ3::Construct);
  }
};

proxy_nps_sharc p_nps_sharc;
}
