/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    : 16/10/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: Define the AnnularS4 detector from Micron                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *  + 11/10/09: Change scorer philosophy, one scorer for the detector number *
 *              added (N. de Sereville)                                      *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <sstream>
#include <string>

// Geant4
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Material.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4PVDivision.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"

// NPTool headers
#include "AnnularS4.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "RootOutput.h"
#include "TAnnularS4Data.h"
// CLHEP
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace ANNULARS4;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AnnularS4::AnnularS4() {
  m_Event = new TAnnularS4Data();
  m_LogicalDetector = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
AnnularS4::~AnnularS4() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void AnnularS4::AddModule(G4double PosZ) { m_PosZ.push_back(PosZ); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* AnnularS4::ConstructVolume() {

  if (!m_LogicalDetector) {
    G4Material* Silicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    ////////////////////////////////////////////////////////////////
    ////////////// Starting Volume Definition //////////////////////
    ////////////////////////////////////////////////////////////////
    // Name of the module
    G4String Name = "AnnularS4";

    // Building the PCB
    // The PCB is a simple extruded volume from 8reference point
    vector<G4TwoVector> polygon;
    for (unsigned int i = 0; i < 4; i++) {
      G4TwoVector Point(PCBPointsX[i], PCBPointsY[i]);
      polygon.push_back(Point);
    }

    // Mast volume containing all the detector
    G4ExtrudedSolid* solidAnnularS4 =
        new G4ExtrudedSolid(Name, polygon, PCBThickness * 0.5, G4TwoVector(0, 0), 1, G4TwoVector(0, 0), 1);

    // Definition of the volume containing the sensitive detector
    m_LogicalDetector = new G4LogicalVolume(solidAnnularS4, Vacuum, Name, 0, 0, 0);
    m_LogicalDetector->SetVisAttributes(G4VisAttributes::GetInvisible());

    // PCB Base
    G4ExtrudedSolid* solidPCBBase =
        new G4ExtrudedSolid("PCBBase", polygon, PCBThickness * 0.5, G4TwoVector(0, 0), 1, G4TwoVector(0, 0), 1);

    // Wafer Shape to be substracted to the PCB
    G4Tubs* solidWaferShape = new G4Tubs("WaferShape", 0, WaferOutterRadius, PCBThickness, 0 * deg, 360 * deg);
    // G4Tubs* solidWaferShapeBase = new G4Tubs("WaferShape", 0, WaferOutterRadius, PCBThickness, 0 * deg, 360 * deg);

    // A no rotation matrix is always handy ;)
    G4RotationMatrix* norotation = new G4RotationMatrix();
    // Rotation of the box that make the Si cut
    // G4RotationMatrix* cutrotation = new G4RotationMatrix(G4ThreeVector(0, 0, 1), 45 * deg);
    // G4ThreeVector cutposition1(80 * mm + WaferRCut, 0, 0);
    // cutposition1.setPhi(45 * deg);
    // G4Transform3D transform1(*cutrotation, cutposition1);

    // G4Box* solidCutout = new G4Box("cuttout", 80 * mm, 80 * mm, 80 * mm);

    // G4SubtractionSolid* solidWaferShape1 =
        // new G4SubtractionSolid("WaferShape1", solidWaferShapeBase, solidCutout, transform1);

    // G4ThreeVector cutposition2(-80 * mm - WaferRCut, 0, 0);
    // cutposition2.setPhi(-135 * deg);
    // G4Transform3D transform2(*cutrotation, cutposition2);
    // G4SubtractionSolid* solidWaferShape =
        // new G4SubtractionSolid("WaferShape", solidWaferShape1, solidCutout, transform2);

    // PCB final
    G4SubtractionSolid* solidPCB = new G4SubtractionSolid("AnnularS4_PCB1", solidPCBBase, solidWaferShape);

    G4LogicalVolume* logicPCB = new G4LogicalVolume(solidPCB, Vacuum, "AnnularS4_PCB", 0, 0, 0);

    new G4PVPlacement(G4Transform3D(*norotation, G4ThreeVector()), logicPCB, "AnnularS4_PCB", m_LogicalDetector, false,
                      0);

    G4VisAttributes* PCBVisAtt = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2));
    logicPCB->SetVisAttributes(PCBVisAtt);

    // Wafer itself
    // G4Tubs* solidWaferBase =
    G4Tubs* solidWafer =
        new G4Tubs("Wafer", WaferInnerRadius, WaferOutterRadius, 0.5 * WaferThickness, 0 * deg, 360 * deg);

    // G4SubtractionSolid* solidWafer1 = new G4SubtractionSolid("Wafer1", solidWaferBase, solidCutout, transform1);

    // G4SubtractionSolid* solidWafer = new G4SubtractionSolid("Wafer", solidWafer1, solidCutout, transform2);

    G4LogicalVolume* logicWafer = new G4LogicalVolume(solidWafer, Silicon, "AnnularS4_Wafer", 0, 0, 0);
    new G4PVPlacement(G4Transform3D(*norotation, G4ThreeVector()), logicWafer, "AnnularS4_Wafer", m_LogicalDetector,
                      false, 0);

    G4VisAttributes* SiVisAtt = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
    logicWafer->SetVisAttributes(SiVisAtt);

    // Active Wafer
    // G4Tubs* solidActiveWaferBase = new G4Tubs("ActiveWafer", ActiveWaferInnerRadius, ActiveWaferOutterRadius,
    G4Tubs* solidActiveWafer = new G4Tubs("ActiveWafer", ActiveWaferInnerRadius, ActiveWaferOutterRadius,
                                              0.5 * WaferThickness, 0 * deg, 360 * deg);

    // G4ThreeVector activecutposition1(80 * mm + ActiveWaferRCut, 0, 0);
    // activecutposition1.setPhi(45 * deg);
    // G4Transform3D activetransform1(*cutrotation, activecutposition1);

    // G4SubtractionSolid* solidActiveWafer1 =
        // new G4SubtractionSolid("ActiveWafer1", solidActiveWaferBase, solidCutout, activetransform1);

    // G4ThreeVector activecutposition2(-80 * mm - ActiveWaferRCut, 0, 0);
    // activecutposition2.setPhi(-135 * deg);
    // G4Transform3D activetransform2(*cutrotation, activecutposition2);

    // G4SubtractionSolid* solidActiveWafer =
        // new G4SubtractionSolid("ActiveWafer", solidActiveWafer1, solidCutout, activetransform2);

    G4LogicalVolume* logicActiveWafer =
        new G4LogicalVolume(solidActiveWafer, Silicon, "AnnularS4_ActiveWafer", 0, 0, 0);
    new G4PVPlacement(G4Transform3D(*norotation, G4ThreeVector()), logicActiveWafer, "AnnularS4_ActiveWafer",
                      logicWafer, false, 0);

    logicActiveWafer->SetVisAttributes(SiVisAtt);

    // Set Silicon strip sensible
    logicActiveWafer->SetSensitiveDetector(m_Scorer);
  }
  return m_LogicalDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void AnnularS4::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("AnnularS4");
  cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> token = {"Z"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(token)) {
      double Z = blocks[i]->GetDouble("Z", "mm");
      AddModule(Z);
    }

    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void AnnularS4::ConstructDetector(G4LogicalVolume* world) {
  G4RotationMatrix* rotation = NULL;
  G4ThreeVector position = G4ThreeVector(0, 0, 0);

  G4int NumberOfModule = m_PosZ.size();

  for (G4int i = 0; i < NumberOfModule; i++) {
    // translation to position the module
    position = G4ThreeVector(0, 0, m_PosZ[i]);

    // Passage Matrix from Lab Referential to Module Referential
    // Identity matrix by default
    rotation = new G4RotationMatrix();
    if (position.z() < 0)
      rotation->rotateX(180 * deg);

    new G4PVPlacement(G4Transform3D(*rotation, position), ConstructVolume(), "AnnularS4", world, false, i + 1);
  }

  delete rotation;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Connect the GaspardTrackingData class to the output TTree
// of the simulation
void AnnularS4::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("AnnularS4")) {
    pTree->Branch("AnnularS4", "TAnnularS4Data", &m_Event);
  }
  pTree->SetBranchAddress("AnnularS4", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void AnnularS4::ReadSensitive(const G4Event*) {
  // Clear ROOT objects
  m_Event->Clear();

  DSSDScorers::PS_Annular* Scorer = (DSSDScorers::PS_Annular*)m_Scorer->GetPrimitive(0);

  // Loop on Silicon Sector Hit
  unsigned int sizeSector = Scorer->GetSectorMult();
  for (unsigned int i = 0; i < sizeSector; i++) {
    double Energy = Scorer->GetEnergyRing(i);

    if (Energy > EnergyThreshold) {
      double Time = Scorer->GetTimeRing(i);
      unsigned int DetNbr = Scorer->GetDetectorRing(i);
      ;
      unsigned int StripPhi = Scorer->GetStripSector(i);

      m_Event->SetS4_E_DetectorNbr(DetNbr);
      m_Event->SetS4_E_StripNbr(StripPhi);
      m_Event->SetS4_E_Energy(RandGauss::shoot(Energy, ResoEnergy));

      m_Event->SetS4_T_DetectorNbr(DetNbr);
      m_Event->SetS4_T_StripNbr(StripPhi);
      m_Event->SetS4_T_Time(RandGauss::shoot(Time, ResoTime));
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Initilize the Scorer use to read out the sensitive volume
void AnnularS4::InitializeScorers() {
  bool already_exist = false;
  // Associate Scorer
  m_Scorer = CheckScorer("AnnularS4_Scorer", already_exist);
  if (already_exist)
    return;

  G4VPrimitiveScorer* AnnularScorer =
      new DSSDScorers::PS_Annular("AnnularS4_Scorer", 2, ActiveWaferInnerRadius, ActiveWaferOutterRadius,
                                  -8 * 22.5 * deg, // MUST2 campaign 2009, See: Phd Sandra Giron
                                  +8 * 22.5 * deg, NbrRingStrips, NbrSectorStrips, NbQuadrant);

  m_Scorer->RegisterPrimitive(AnnularScorer);
  G4VPrimitiveScorer* InteractionScorer = new InteractionScorers::PS_Interactions("InteractionAnnularS4", ms_InterCoord, 2);
  m_Scorer->RegisterPrimitive(InteractionScorer);
  //  Add All Scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_Scorer);
}
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* AnnularS4::Construct() { return (NPS::VDetector*)new AnnularS4(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_annulars4 {
 public:
  proxy_nps_annulars4() {
    NPS::DetectorFactory::getInstance()->AddToken("AnnularS4", "AnnularS4");
    NPS::DetectorFactory::getInstance()->AddDetector("AnnularS4", AnnularS4::Construct);
  }
};

proxy_nps_annulars4 p_nps_annulars4;
}

