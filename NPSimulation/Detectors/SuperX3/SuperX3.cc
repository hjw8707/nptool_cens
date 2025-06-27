/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 12/01/11                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription: Define the SuperX3 detector from Micron *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <sstream>
#include <string>

// G4 Geometry headers
#include "G4Box.hh"
#include "G4Tubs.hh"

// G4 various headers
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// NPTool headers
#include "DSSDScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"

#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "SuperX3.hh"
#include "TSuperX3Data.h"

// CLHEP
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace SuperX3SQUARE;

SuperX3::SuperX3() : m_Event(new TSuperX3Data) { InitializeMaterials(); }

SuperX3::~SuperX3() { delete m_Event; }

void SuperX3::AddDetector(G4ThreeVector X1_Y1, G4ThreeVector X16_Y1, G4ThreeVector X1_Y16, G4ThreeVector X16_Y16) {
  m_DefinitionType.push_back(true);

  m_X1_Y1.push_back(X1_Y1);
  m_X16_Y1.push_back(X16_Y1);
  m_X1_Y16.push_back(X1_Y16);
  m_X16_Y16.push_back(X16_Y16);

  m_R.push_back(0);
  m_Theta.push_back(0);
  m_Phi.push_back(0);
  m_beta_u.push_back(0);
  m_beta_v.push_back(0);
  m_beta_w.push_back(0);
}

void SuperX3::AddDetector(G4double R, G4double Theta, G4double Phi, G4double beta_u, G4double beta_v, G4double beta_w) {
  G4ThreeVector empty = G4ThreeVector(0, 0, 0);

  m_DefinitionType.push_back(false);

  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_beta_u.push_back(beta_u);
  m_beta_v.push_back(beta_v);
  m_beta_w.push_back(beta_w);

  m_X1_Y1.push_back(empty);
  m_X16_Y1.push_back(empty);
  m_X1_Y16.push_back(empty);
  m_X16_Y16.push_back(empty);
}

void SuperX3::VolumeMaker(G4int DetecNumber, G4ThreeVector position, G4RotationMatrix* rotation,
                          G4LogicalVolume* world) {
  G4double NbrTelescopes = DetecNumber;
  G4String DetectorNumber;
  ostringstream Number;
  Number << NbrTelescopes;
  DetectorNumber = Number.str();

  ////////////////////////////////////////////////////////////////
  ////////////// Starting Volume Definition //////////////////////
  ////////////////////////////////////////////////////////////////
  G4String Name = "SuperX3Square" + DetectorNumber;

  // Definition of the volume containing the sensitive detector
  G4Box* solidSuperX3 = new G4Box(Name, 0.5 * FaceFront, 0.5 * SiliconFaceLength * mm, 0.5 * Length);
  G4LogicalVolume* logicSuperX3 = new G4LogicalVolume(solidSuperX3, m_MaterialVacuum, Name, 0, 0, 0);

  new G4PVPlacement(G4Transform3D(*rotation, position), logicSuperX3, Name, world, false, DetecNumber);

  logicSuperX3->SetVisAttributes(G4VisAttributes::GetInvisible());
  if (m_non_sensitive_part_visiualisation)
    logicSuperX3->SetVisAttributes(G4VisAttributes(G4Colour(0.90, 0.90, 0.90)));

  // Aluminium dead layers
  G4ThreeVector positionAluStripFront = G4ThreeVector(0, 0, AluStripBack_PosZ);
  G4ThreeVector positionAluStripBack = G4ThreeVector(0, 0, AluStripFront_PosZ);

  G4Box* solidAluStrip = new G4Box("AluBox", 0.5 * SiliconFaceWidth, 0.5 * SiliconFaceLength, 0.5 * AluStripThickness);
  //   G4LogicalVolume* logicAluStrip = new G4LogicalVolume(solidAluStrip,
  //   m_MaterialAluminium, "logicAluStrip", 0, 0, 0);
  G4LogicalVolume* logicAluStrip = new G4LogicalVolume(solidAluStrip, m_MaterialVacuum, "logicAluStrip", 0, 0, 0);

  new G4PVPlacement(0, positionAluStripFront, logicAluStrip, Name + "_AluStripFront", logicSuperX3, false, DetecNumber);
  new G4PVPlacement(0, positionAluStripBack, logicAluStrip, Name + "_AluStripBack", logicSuperX3, false, DetecNumber);

  logicAluStrip->SetVisAttributes(G4VisAttributes(G4Colour(1.90, 0.0, 0.0)));

  // Silicon detector itself
  G4ThreeVector positionSilicon = G4ThreeVector(0, 0, Silicon_PosZ);

  G4Box* solidSilicon =
      new G4Box("solidSilicon", 0.5 * SiliconFaceWidth, 0.5 * SiliconFaceLength * mm, 0.5 * SiliconThickness);
  G4LogicalVolume* logicSilicon = new G4LogicalVolume(solidSilicon, m_MaterialSilicon, "logicSilicon", 0, 0, 0);

  new G4PVPlacement(0, positionSilicon, logicSilicon, Name + "_Silicon", logicSuperX3, false, DetecNumber);

  // Set Silicon strip sensible
  logicSilicon->SetSensitiveDetector(m_Scorer);

  /// Visualisation of Silicon Strip
  G4VisAttributes* SiliconVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  logicSilicon->SetVisAttributes(SiliconVisAtt);
}

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void SuperX3::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SuperX3");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;
  for (unsigned int i = 0; i < blocks.size(); i++) {
    // Cartesian Case
    vector<string> cart = {"X1_Y1", "X1_Y16", "X16_Y1", "X16_Y16", "VIS"};
    // Spherical Case
    vector<string> sphe = {"R", "THETA", "PHI", "BETA", "VIS"};

    if (blocks[i]->HasTokenList(cart)) {
      cout << endl << "////  SuperX3 " << i + 1 << endl;
      G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("X16_Y1", "mm"));
      G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y16", "mm"));
      G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("X16_Y16", "mm"));
      if (blocks[i]->GetInt("VIS"))
        m_non_sensitive_part_visiualisation = true;
      AddDetector(A, B, C, D);
    }

    else if (blocks[i]->HasTokenList(sphe)) {
      cout << endl << "////  SuperX3 " << i + 1 << endl;
      double Theta = blocks[i]->GetDouble("THETA", "deg");
      double Phi = blocks[i]->GetDouble("PHI", "deg");
      double R = blocks[i]->GetDouble("R", "mm");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA", "deg");
      if (blocks[i]->GetInt("VIS"))
        m_non_sensitive_part_visiualisation = true;

      AddDetector(R, Theta, Phi, beta[0], beta[1], beta[2]);
    }

    else {
      cout << "ERROR: Missing token for SuperX3 blocks, check your input file" << endl;
      exit(1);
    }
  }
}

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void SuperX3::ConstructDetector(G4LogicalVolume* world) {
  G4RotationMatrix* SuperX3rot = NULL;
  G4ThreeVector SuperX3pos = G4ThreeVector(0, 0, 0);
  G4ThreeVector SuperX3u = G4ThreeVector(0, 0, 0);
  G4ThreeVector SuperX3v = G4ThreeVector(0, 0, 0);
  G4ThreeVector SuperX3w = G4ThreeVector(0, 0, 0);
  G4ThreeVector SuperX3Center = G4ThreeVector(0, 0, 0);

  G4int NumberOfDetector = m_DefinitionType.size();
  for (G4int i = 0; i < NumberOfDetector; i++) {
    // By Point
    if (m_DefinitionType[i]) {
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      SuperX3u = m_X16_Y1[i] - m_X1_Y1[i];
      SuperX3u = SuperX3u.unit();

      SuperX3v = m_X1_Y16[i] - m_X1_Y1[i];
      SuperX3v = SuperX3v.unit();

      SuperX3w = SuperX3u.cross(SuperX3v);
      SuperX3w = SuperX3w.unit();

      SuperX3Center = (m_X1_Y1[i] + m_X1_Y16[i] + m_X16_Y1[i] + m_X16_Y16[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      SuperX3rot = new G4RotationMatrix(SuperX3u, SuperX3v, SuperX3w);
      // translation to place Telescope
      SuperX3pos = SuperX3w * Length * 0.5 + SuperX3Center;
    }

    // By Angle
    else {
      G4double Theta = m_Theta[i];
      G4double Phi = m_Phi[i];

      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      // Phi is angle between X axis and projection in (X,Y) plan
      // Theta is angle between  position vector and z axis
      G4double wX = m_R[i] * sin(Theta) * cos(Phi);
      G4double wY = m_R[i] * sin(Theta) * sin(Phi);
      G4double wZ = m_R[i] * cos(Theta);
      SuperX3w = G4ThreeVector(wX, wY, wZ);

      // vector corresponding to the center of the module
      SuperX3Center = SuperX3w;

      // vector parallel to one axis of silicon plane
      G4double ii = cos(Theta) * cos(Phi);
      G4double jj = cos(Theta) * sin(Phi);
      G4double kk = -sin(Theta);
      G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

      SuperX3w = SuperX3w.unit();
      SuperX3u = SuperX3w.cross(Y);
      SuperX3v = SuperX3w.cross(SuperX3u);
      SuperX3v = SuperX3v.unit();
      SuperX3u = SuperX3u.unit();

      // Passage Matrix from Lab Referential to Telescope Referential
      // MUST2
      SuperX3rot = new G4RotationMatrix(SuperX3u, SuperX3v, SuperX3w);
      // Telescope is rotate of Beta angle around SuperX3v axis.
      SuperX3rot->rotate(m_beta_u[i], SuperX3u);
      SuperX3rot->rotate(m_beta_v[i], SuperX3v);
      SuperX3rot->rotate(m_beta_w[i], SuperX3w);
      // translation to place Telescope
      SuperX3pos = SuperX3w * (SiliconThickness + AluStripThickness) * 0.5 + SuperX3Center;
    }

    VolumeMaker(i + 1, SuperX3pos, SuperX3rot, world);
  }

  delete SuperX3rot;
}

// Connect the GaspardTrackingData class to the output TTree
// of the simulation
void SuperX3::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("SuperX3")) {
    pTree->Branch("SuperX3", "TSuperX3Data", &m_Event);
  }
  pTree->SetBranchAddress("SuperX3", &m_Event);
}

// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void SuperX3::ReadSensitive(const G4Event*) {
  // Clear ROOT objects
  m_Event->Clear();

  auto resistive = (DSSDScorers::PS_Resistive*)m_Scorer->GetPrimitive(0);
  auto backstrip = (DSSDScorers::PS_Rectangle*)m_Scorer->GetPrimitive(1);
  auto sizeUp = resistive->GetUpMult();
  for (unsigned int i = 0; i < sizeUp; i++) {
    double energy = resistive->GetEnergyUp(i);
    double time = resistive->GetTimeUp(i);
    int det = resistive->GetDetectorUp(i);
    int strip = resistive->GetStripUp(i);
    if (energy > EnergyThresholdResitive) {
      m_Event->SetUpE(det, strip, energy);
      m_Event->SetUpT(det, strip, time);
    }
    else {
      m_Event->SetUpE(det, strip, 0);
      m_Event->SetUpT(det, strip, 0);
    }
  }
  auto sizeDown = resistive->GetDownMult();
  for (unsigned int i = 0; i < sizeDown; i++) {
    double energy = resistive->GetEnergyDown(i);
    double time = resistive->GetTimeDown(i);
    int det = resistive->GetDetectorDown(i);
    int strip = resistive->GetStripDown(i);
    if (energy > EnergyThresholdResitive) {
      m_Event->SetDownE(det, strip, energy);
      m_Event->SetDownT(det, strip, time);
    }

    else {
      m_Event->SetDownE(det, strip, 0);
      m_Event->SetDownT(det, strip, 0);
    }
  }
  auto sizeBack = backstrip->GetLengthMult();
  for (unsigned int i = 0; i < sizeBack; i++) {
    double energy = backstrip->GetEnergyLength(i);
    double time = backstrip->GetTimeLength(i);
    int det = backstrip->GetDetectorLength(i);
    int strip = backstrip->GetStripLength(i);
    if (energy > EnergyThresholdBack) {
      m_Event->SetBackE(det, strip, energy);
      m_Event->SetBackT(det, strip, time);
    }
    else {
      m_Event->SetBackE(det, strip, 0);
      m_Event->SetBackT(det, strip, 0);
    }
  }
  resistive->clear();
  backstrip->clear();
}

void SuperX3::InitializeMaterials() {
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  m_MaterialAluminium = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialIron = MaterialManager::getInstance()->GetMaterialFromLibrary("Fe");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
}

void SuperX3::InitializeScorers() {
  bool already_exist = false;
  // Associate Scorer
  m_Scorer = CheckScorer("ScorerSuperX3", already_exist);
  if (already_exist)
    return;

  //..... resistive starts..
  G4VPrimitiveScorer* resistivestrip =
      new DSSDScorers::PS_Resistive("resistivestrip", "y", 0, SiliconFaceLength, SiliconFaceWidth, NbStrips, 0);
  G4VPrimitiveScorer* backstrip =
      new DSSDScorers::PS_Rectangle("backstrip", 0, SiliconFaceLength, SiliconFaceWidth, 4, 1, 0);

  G4VPrimitiveScorer* interaction = new InteractionScorers::PS_Interactions("Interaction", ms_InterCoord, 0);
  //... resistive ends......
  // and register it to the multifunctionnal detector
  //.... resistive starts...
  m_Scorer->RegisterPrimitive(resistivestrip);
  m_Scorer->RegisterPrimitive(backstrip);
  m_Scorer->RegisterPrimitive(interaction);
  //.....resistive ends...

  //  Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_Scorer);
}
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SuperX3::Construct() { return (NPS::VDetector*)new SuperX3(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_w1 {
 public:
  proxy_nps_w1() {
    NPS::DetectorFactory::getInstance()->AddToken("SuperX3", "SuperX3");
    NPS::DetectorFactory::getInstance()->AddDetector("SuperX3", SuperX3::Construct);
  }
};

proxy_nps_w1 p_nps_w1;
}
