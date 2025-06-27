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

// C++ headers
#include <cmath>
#include <limits>
#include <sstream>
// G4 Geometry object
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4AssemblyVolume.hh"
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"

// NPTool header
#include "CACAO.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace CACAO_NS {
  // Geometry
  const double PCBThickness = 1.6 * mm;

  // Energy and time Resolution
  const double EnergyThreshold = 0.1 * MeV;
  const double ResoTime = 4.5 * ns;
  double ResoEnergy = 0.027455; // dE = Reso*Sqrt(E) where E in MeV
  const double Radius = 50 * mm;
  const double Width = 100 * mm;
  const double Thickness = 300 * mm;
  const string ScintMaterial = "CsI";
  const string ShieldMaterial = "Mylar";
  const string ChamberMaterial = "Al";
} // namespace CACAO_NS

using namespace CACAO_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// CACAO Specific Method
CACAO::CACAO() {
  m_Event = new TCACAOData();
  m_CACAOScorer = 0;

  m_matScint = 0;
  m_matShield = 0;

  m_ChamberFound = false;
  m_ChamberTRmin = m_ChamberTRmax = m_ChamberTZ0 = m_ChamberTZ1 = 0;
  m_ChamberCRmin = m_ChamberCRmax = m_ChamberTZ2 = m_ChamberTZ3 = 0;
  m_ChamberLRmin = m_ChamberLRmax = m_ChamberLH = 0;
  m_ChamberBRmin = m_ChamberBRmax = m_ChamberBZ1 = m_ChamberBZ2 = 0;
  m_ChamberBZ3 = 0;

  // RGB Color + Transparency
  m_VisScint = new G4VisAttributes(G4Colour(243 / 255., 198 / 255., 165 / 255., 1));
  m_VisPCB = new G4VisAttributes(G4Colour(71 / 255., 255 / 255., 87 / 255., 1));
}

CACAO::~CACAO() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CACAO::DefineMaterials() {

  G4Element* Tl = new G4Element("Thallium", "Tl", 81., 204.383 * g / mole);

  m_matScint = new G4Material("CsITl", 4.51 * g / cm3, 2);
  m_matScint->AddMaterial(MaterialManager::getInstance()->GetMaterialFromLibrary(ScintMaterial), 99.6 * perCent);
  m_matScint->AddElement(Tl, 0.4 * perCent);

  m_matShield = MaterialManager::getInstance()->GetMaterialFromLibrary(ShieldMaterial);
  m_matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  m_matChamber = MaterialManager::getInstance()->GetMaterialFromLibrary(ChamberMaterial);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CACAO::AddDetector(G4ThreeVector Pos, G4RotationMatrix Rot, G4ThreeVector Dim, G4double ShieldThickness) {
  m_Pos.push_back(Pos);
  m_Rot.push_back(Rot);
  m_Dim.push_back(Dim);
  m_ShieldThickness.push_back(ShieldThickness);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* CACAO::BuildDetector(G4int i) {

  G4Box* solidScint = new G4Box("Scint_Box", m_Dim[i].x() * 0.5, m_Dim[i].y() * 0.5, m_Dim[i].z() * 0.5);

  G4LogicalVolume* logicScint = new G4LogicalVolume(solidScint, m_matScint, "logicScint", 0, 0, 0);
  logicScint->SetVisAttributes(m_VisScint);
  logicScint->SetSensitiveDetector(m_CACAOScorer);

  G4VSolid* solidPCB = new G4Box("PCB_Box", m_Dim[i].x() * 0.5, m_Dim[i].y() * 0.5, PCBThickness * 0.5);
  G4LogicalVolume* logicPCB = new G4LogicalVolume(solidPCB, m_matPCB, "logicPCB", 0, 0, 0);
  logicPCB->SetVisAttributes(m_VisPCB);
  new G4PVPlacement(0, G4ThreeVector(0, 0, -m_Dim[i].z() * 0.5 - PCBThickness * 0.5), logicPCB, "solidPCB", logicScint,
                    false, 0);

  if (m_ShieldThickness[i] > 0) {
    G4Box* solidShieldAll =
        new G4Box("Shield_Box", m_Dim[i].x() * 0.5 + m_ShieldThickness[i], m_Dim[i].y() * 0.5 + m_ShieldThickness[i],
                  m_Dim[i].z() * 0.5 + m_ShieldThickness[i]);
    G4VSolid* solidShield = new G4SubtractionSolid("Shield", solidShieldAll, solidScint);
    G4LogicalVolume* logicShield = new G4LogicalVolume(solidShield, m_matShield, "logicShield", 0, 0, 0);
    new G4PVPlacement(0, G4ThreeVector(), logicShield, "solidShield", logicScint, false, 0);
  }

  return logicScint;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void CACAO::ReadConfiguration(NPL::InputParser parser) {
  ///////////////////////////////////////////////////////////////////////////////////////
  // Detector
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("CACAO");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> reso = {"Reso"};

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R", "Theta", "Phi"};
  vector<string> cyli = {"Rho", "Phi", "Z"};

  vector<string> cuboid = {"DIM", "ShieldThickness"};

  for (unsigned int i = 0; i < blocks.size(); i++) {

    ////////////////////////////////////////////////////////////
    // Resolution
    if (blocks[i]->HasTokenList(reso)) {
      CACAO_NS::ResoEnergy = blocks[i]->GetDouble("Reso", "void");
      continue;
    }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // checking the position items of the block
    G4ThreeVector Pos;
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO " << i + 1 << endl;
      Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS", "mm"));
    }
    else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO " << i + 1 << endl;
      Pos.setRThetaPhi(blocks[i]->GetDouble("R", "mm"), blocks[i]->GetDouble("Theta", "deg"),
                       blocks[i]->GetDouble("Phi", "deg"));
    }
    else if (blocks[i]->HasTokenList(cyli)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO " << i + 1 << endl;
      Pos.setRhoPhiZ(blocks[i]->GetDouble("Rho", "mm"), blocks[i]->GetDouble("Phi", "deg"),
                     blocks[i]->GetDouble("Z", "mm"));
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }

    G4RotationMatrix Rot;
    if (blocks[i]->HasToken("ANG")) {
      G4ThreeVector Ang;
      Ang = NPS::ConvertVector(blocks[i]->GetTVector3("ANG", "deg"));
      Rot.set(Ang.x(), Ang.y(), Ang.z());
    } // (phi, theta, psi)
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // checking the shape items of the block
    G4ThreeVector Dim;
    if (blocks[i]->HasTokenList(cuboid)) {
      Dim = NPS::ConvertVector(blocks[i]->GetTVector3("DIM", "mm"));
      double ShieldThickness = blocks[i]->GetDouble("ShieldThickness", "mm");
      AddDetector(Pos, Rot, Dim, ShieldThickness);
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
    ////////////////////////////////////////////////////////////
  }
  ///////////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////////////
  // Chamber
  blocks.clear();
  blocks = parser.GetAllBlocksWithToken("CACAOChamber");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " CACAO Chamber found " << endl;
  vector<string> token = {"TRmin", "TRmax", "TZ0", "TZ1", "CRmin", "CRmax", "TZ2", "TZ3"};
  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(token)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO Chamber " << i + 1 << endl;
      m_ChamberFound = true;
      m_ChamberTRmin = blocks[i]->GetDouble("TRmin", "mm");
      m_ChamberTRmax = blocks[i]->GetDouble("TRmax", "mm");
      m_ChamberTZ0 = blocks[i]->GetDouble("TZ0", "mm");
      m_ChamberTZ1 = blocks[i]->GetDouble("TZ1", "mm");
      m_ChamberCRmin = blocks[i]->GetDouble("CRmin", "mm");
      m_ChamberCRmax = blocks[i]->GetDouble("CRmax", "mm");
      m_ChamberTZ2 = blocks[i]->GetDouble("TZ2", "mm");
      m_ChamberTZ3 = blocks[i]->GetDouble("TZ3", "mm");
      m_ChamberLRmin = m_ChamberTRmin;
      m_ChamberLRmax = m_ChamberTRmax;
      m_ChamberLH = 200; // mm
      m_ChamberLRmin = blocks[i]->GetDouble("LRmin", "mm");
      m_ChamberLRmax = blocks[i]->GetDouble("LRmax", "mm");
      m_ChamberLH = blocks[i]->GetDouble("LH", "mm");
      m_ChamberBRmin = blocks[i]->GetDouble("BRmin", "mm");
      m_ChamberBRmax = blocks[i]->GetDouble("BRmax", "mm");
      m_ChamberBZ1 = blocks[i]->GetDouble("BZ1", "mm");
      m_ChamberBZ2 = blocks[i]->GetDouble("BZ2", "mm");      
      //m_ChamberBZ1 = blocks[i]->GetDouble("BZ2", "mm");
    }
    else {
      cout << "Warning: check your input file formatting " << endl;
    }
  }
  ///////////////////////////////////////////////////////////////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void CACAO::ConstructDetector(G4LogicalVolume* world) {
  DefineMaterials();
  if (m_ChamberFound)
    ConstructChamber(world);
  for (unsigned short i = 0; i < m_Pos.size(); i++) {
    new G4PVPlacement(G4Transform3D(m_Rot[i], m_Pos[i]), BuildDetector(i), "CACAO", world, false, i + 1);
  }
}

void CACAO::ConstructChamber(G4LogicalVolume* world) {

  G4Tubs* tubsChamber0All =
      new G4Tubs("tubsChamber0All", m_ChamberTRmin, m_ChamberTRmax, (m_ChamberTZ1 - m_ChamberTZ0) / 2, 0, 2 * pi);
  G4Tubs* tubsChamber0Sub = new G4Tubs("tubsChamber0Sub", 0, m_ChamberLRmax, m_ChamberTRmax, 0, 2 * pi);
  G4RotationMatrix* rot = new G4RotationMatrix;
  rot->rotateX(pi / 2.);
  G4SubtractionSolid* tubsChamber0Part =
      new G4SubtractionSolid("tubsChamber0Part", tubsChamber0All, tubsChamber0Sub, rot,
                             G4ThreeVector(0, m_ChamberTRmax, -(m_ChamberTZ1 + m_ChamberTZ0) / 2));
  G4Tubs* lh2Chamber = new G4Tubs("lh2Chamber", m_ChamberLRmin, m_ChamberLRmax, m_ChamberLH / 2, 0, 2 * pi);
  G4UnionSolid* tubsChamber0 = new G4UnionSolid(
      "tubsChamber0", tubsChamber0Part, lh2Chamber, rot,
      G4ThreeVector(0, m_ChamberLH / 2 + sqrt(m_ChamberTRmin * m_ChamberTRmin - m_ChamberLRmin * m_ChamberLRmin),
                    -(m_ChamberTZ1 + m_ChamberTZ0) / 2.));
  G4LogicalVolume* lTubsChamber0 = new G4LogicalVolume(tubsChamber0, m_matChamber, "lTubsChamber0", 0, 0, 0);

  G4Cons* consChamber = new G4Cons("consChamber", m_ChamberTRmin, m_ChamberTRmax, m_ChamberCRmin, m_ChamberCRmax,
                                   (m_ChamberTZ2 - m_ChamberTZ1) / 2, 0, 2 * pi);
  G4LogicalVolume* lConsChamber = new G4LogicalVolume(consChamber, m_matChamber, "lConsChamber", 0, 0, 0);
  G4Tubs* tubsChamber1 =
      new G4Tubs("tubsChamber1", m_ChamberCRmin, m_ChamberCRmax, (m_ChamberTZ3 - m_ChamberTZ2) / 2, 0, 2 * pi);
  G4LogicalVolume* lTubsChamber1 = new G4LogicalVolume(tubsChamber1, m_matChamber, "lTubsChamber1", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., (m_ChamberTZ1 + m_ChamberTZ0) / 2), lTubsChamber0, "pTubsChamber0", world,
                    false, 0);
  new G4PVPlacement(0, G4ThreeVector(0., 0., (m_ChamberTZ2 + m_ChamberTZ1) / 2), lConsChamber, "pConsChamber", world,
                    false, 0);
  new G4PVPlacement(0, G4ThreeVector(0., 0., (m_ChamberTZ3 + m_ChamberTZ2) / 2), lTubsChamber1, "pTubsChamber1", world,
                    false, 0);

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Backward chamber
  G4Cons* consBChamber = new G4Cons("consBChamber", m_ChamberBRmin, m_ChamberBRmax, m_ChamberTRmin, m_ChamberTRmax,
                                   (m_ChamberTZ0 - m_ChamberBZ1) / 2, 0, 2 * pi);
  G4LogicalVolume* lConsBChamber = new G4LogicalVolume(consBChamber, m_matChamber, "lConsBChamber", 0, 0, 0);
  G4Tubs* tubsBChamber1 =
      new G4Tubs("tubsBChamber1", m_ChamberBRmin, m_ChamberBRmax, (m_ChamberBZ1 - m_ChamberBZ2) / 2, 0, 2 * pi);
  G4LogicalVolume* lTubsBChamber1 = new G4LogicalVolume(tubsBChamber1, m_matChamber, "lTubsBChamber1", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0., 0., (m_ChamberTZ0 + m_ChamberBZ1) / 2), lConsBChamber, "pConsBChamber", world,
                    false, 0);
  new G4PVPlacement(0, G4ThreeVector(0., 0., (m_ChamberBZ1 + m_ChamberBZ2) / 2), lTubsBChamber1, "pTubsBChamber1", world,
                    false, 0);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(0., 1., 1.));
  lTubsChamber0->SetVisAttributes(ChamberVisAtt);
  lConsChamber->SetVisAttributes(ChamberVisAtt);
  lTubsChamber1->SetVisAttributes(ChamberVisAtt);
  lConsBChamber->SetVisAttributes(ChamberVisAtt);
  lTubsBChamber1->SetVisAttributes(ChamberVisAtt);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void CACAO::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("CACAO")) {
    pTree->Branch("CACAO", "TCACAOData", &m_Event);
  }
  pTree->SetBranchAddress("CACAO", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void CACAO::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer = (CalorimeterScorers::PS_Calorimeter*)m_CACAOScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  for (unsigned int i = 0; i < size; i++) {
    vector<unsigned int> level = Scorer->GetLevel(i);
    // double Energy = RandGauss::shoot(Scorer->GetEnergy(i),CACAO_NS::ResoEnergy);
    // double Energy = RandGauss::shoot(Scorer->GetEnergy(i), 100. * CACAO_NS::ResoEnergy * (Scorer->GetEnergy(i)));
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i), ResoEnergy * std::sqrt(Scorer->GetEnergy(i)));

    if (Energy > CACAO_NS::EnergyThreshold) {
      double Time = RandGauss::shoot(Scorer->GetTime(i), CACAO_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->Set(DetectorNbr, Energy, Time);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void CACAO::InitializeScorers() {
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  m_CACAOScorer = CheckScorer("CACAOScorer", already_exist);

  if (already_exist)
    return;

  // Otherwise the scorer is initialised
  vector<int> level;
  level.push_back(0);
  G4VPrimitiveScorer* Calorimeter = new CalorimeterScorers::PS_Calorimeter("Calorimeter", level, 0);
  G4VPrimitiveScorer* Interaction = new InteractionScorers::PS_Interactions("Interaction", ms_InterCoord, 0);
  // and register it to the multifunctionnal detector
  m_CACAOScorer->RegisterPrimitive(Calorimeter);
  m_CACAOScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CACAOScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* CACAO::Construct() { return (NPS::VDetector*)new CACAO(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_CACAO {
 public:
  proxy_nps_CACAO() {
    NPS::DetectorFactory::getInstance()->AddToken("CACAO", "CACAO");
    NPS::DetectorFactory::getInstance()->AddDetector("CACAO", CACAO::Construct);
  }
};

proxy_nps_CACAO p_nps_CACAO;
}
