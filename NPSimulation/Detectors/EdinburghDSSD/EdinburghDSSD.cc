/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : mars 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  EdinburghDSSD simulation                             *
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
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Transform3D.hh"
#include "G4TwoVector.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"

// NPTool header
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "EdinburghDSSD.hh"
#include "EdinburghDSSDMap.h"
#include "NPCore.h"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/Random.h"

using namespace std;
using namespace CLHEP;
using namespace EDIN_NS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


using namespace EDIN_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// EdinburghDSSD Specific Method
EdinburghDSSD::EdinburghDSSD(){
  m_Event = new TEdinburghDSSDData();
  m_SquareScorer = 0;
  m_SquareDetector = 0;
  m_SquareDeadLayer =0;
  m_subWorld = 0;
  m_substract = 0;
  m_Foil = 0;
  m_FoilHolder = 0;
  m_FirstSledge = 0;
  m_SecondSledge = 0;
  m_Filler = 0;
  m_bar = 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 0, 1, 0.8));

  m_VisDeadLayer = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));
}

EdinburghDSSD::~EdinburghDSSD(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EdinburghDSSD::AddDetector(int DetectorNumber, G4ThreeVector PX1_Y1, G4ThreeVector PX1_Y16,G4ThreeVector PX16_Y1, G4ThreeVector PX16_Y16) {
  m_X1_Y1.push_back(PX1_Y1);         // Top Left Corner Position Vector
  m_X1_Y16.push_back(PX1_Y16);     // Bottom Left Corner Position Vector
  m_X16_Y1.push_back(PX16_Y1);     // Bottom Right Corner Position Vector
  m_X16_Y16.push_back(PX16_Y16); // Center Corner Position Vector
  m_DetectorNumber.push_back(DetectorNumber);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildSubWorld() {
  if(!m_subWorld)
  {
    G4Tubs* subWorldSolid = new G4Tubs("sub_world", 0., (BarHeight + 0.5) * 0.5, (BarLength + 0.5) * 0.5, 0., twopi);

    // Find or build the material for the SubWorld (e.g., vacuum)
    G4Material* subWorldMaterial = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");

    // Create the logical volume for the SubWorld
    m_subWorld = new G4LogicalVolume(subWorldSolid, subWorldMaterial, "sub_world");

    // Set visualization attributes for the SubWorld (optional)
    m_subWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  }
    return m_subWorld;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildSquareDetector() {
  if (!m_SquareDetector) {
    G4String Name = "EDIN";

    G4Box* solidSquare = new G4Box(Name, 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * SiliconThickness);

    G4LogicalVolume* logicSquare = new G4LogicalVolume(
    solidSquare, MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"), Name, 0, 0, 0);


    logicSquare->SetVisAttributes(m_VisSquare);

    G4ThreeVector positionFirstStage = G4ThreeVector(0, 0, 0);

    G4Box* solidFirstStage = new G4Box("solidFirstStage", 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * (SiliconThickness));
    G4LogicalVolume* logicFirstStage = new G4LogicalVolume(
        solidFirstStage, MaterialManager::getInstance()->GetMaterialFromLibrary("Si"), "logicFirstStage", 0, 0, 0);


    logicFirstStage->SetVisAttributes(m_VisSquare);
    new G4PVPlacement(0, positionFirstStage, logicFirstStage, Name + "_FirstStage", logicSquare, false, 0);

    m_SquareDetector = logicSquare;

    G4UserLimits* stepLimit = new G4UserLimits(10*um);
    logicFirstStage->SetUserLimits(stepLimit);
    // Set First Stage sensible
    logicFirstStage->SetSensitiveDetector(m_SquareScorer);


  }

  return m_SquareDetector;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildDeadLayer(){
  if (!m_SquareDeadLayer) {
    G4String Name = "EDIN";

    // Step 1: Create the Square Solid
    G4Box* solidSquare = new G4Box(Name, 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * DeadLayerThickness);
    G4LogicalVolume* logicSquare = new G4LogicalVolume(
      solidSquare, MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"), Name, 0, 0, 0);

    // Step 2: Create the Boron Dead Layer Solid
    G4Box* solidDeadLayer = new G4Box("solidDeadLayer", 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * DeadLayerThickness);

    // Step 3: Create Logical Volume for the Dead Layer
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* deadLayerMaterial = nistManager->FindOrBuildMaterial("G4_Si");

    G4LogicalVolume* logicDeadLayer = new G4LogicalVolume(solidDeadLayer, deadLayerMaterial, "logicDeadLayer");

    // Set the visualization attributes for the Dead Layer
    G4VisAttributes* DeadLayerVisAtt = new G4VisAttributes(G4Colour(0.1, 0.90, 0.90));
    //DeadLayerVisAtt->SetForceWireframe(true);
    logicDeadLayer->SetVisAttributes(DeadLayerVisAtt);

    // Step 4: Position the Dead Layer - Assuming it's placed on top of the first stage
    G4ThreeVector positionDeadLayer = G4ThreeVector(0, 0, 0.5 * SiliconThickness + 0.5 * DeadLayerThickness);
    new G4PVPlacement(0, positionDeadLayer, logicDeadLayer, (Name + "_DeadLayer").c_str(), logicSquare, false, 0);

    m_SquareDeadLayer = logicDeadLayer;

  }
  return m_SquareDeadLayer;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildFoil() {
  if(!m_Foil)
  {
    // Create the foil solid
    G4Tubs* solidFoil = new G4Tubs("Foil", 0., FoilRadius, FoilThickness * 0.5, 0., twopi);

    // Find or build the carbon material for the foil
    G4Material* C_material = G4NistManager::Instance()->FindOrBuildMaterial("G4_C");

    // Create the logical volume for the foil
    G4LogicalVolume* logicFoil = new G4LogicalVolume(solidFoil, C_material, "Foil");

    // Set visualization attributes for the foil
    logicFoil->SetVisAttributes(G4Colour(1.0, 0.0, 0.0, 1.0));  // Red color for the foil
    m_Foil = logicFoil;
  }
    return m_Foil;  // Return the logical volume
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildNozzle() {

    const G4int num_n = 8;
    G4double r_n[num_n] = {1.363, 1.932, 3.799, 6.839, 6.839, 16.986, 16.935, 1.363};
    G4double z_n[num_n] = {35., 35., 31.355, 4.262, 2.019, 1.969, 0.011, 0.011};

    G4VSolid* solidNozzle = new G4Polycone("Nozzle", 0., twopi, num_n, r_n, z_n);
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* HolderMaterial = nistManager->FindOrBuildMaterial("G4_Al");
    G4LogicalVolume* logicNozzle = new G4LogicalVolume(solidNozzle, HolderMaterial, "Nozzel");


    return logicNozzle;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildFoilHolder() {
  if(!m_FoilHolder)
  {
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* FoilHolderMaterial = nistManager->FindOrBuildMaterial("G4_Al");

  double halfHeight = FoilHolderThickness / 2;
  G4Tubs* Cylinder = new G4Tubs("Cylinder", FoilHolderInnerRadius, FoilHolderOuterRadius, halfHeight, 0.*deg, 360.*deg);
  G4LogicalVolume* logicFoilHolder = new G4LogicalVolume(Cylinder, FoilHolderMaterial, "FoilHolderLV");
  logicFoilHolder->SetVisAttributes(G4Colour(1.0, 0.5, 0.0));
  m_FoilHolder = logicFoilHolder;
  }
  return m_FoilHolder;  // Return the logical volume
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildFirstSledge() {
  if(!m_FirstSledge)
  {
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* SledgeMaterial = nistManager->FindOrBuildMaterial("G4_Al");

  double halfHeight = FirstSledgeThickness / 2;
  G4Tubs* Cylinder = new G4Tubs("Cylinder", FirstSledgeInnerRadius, FirstSledgeOuterRadius, halfHeight, 0.*deg, 360.*deg);
  G4LogicalVolume* logicFirstSledge = new G4LogicalVolume(Cylinder, SledgeMaterial, "FirtSledgeLV");
  m_FirstSledge = logicFirstSledge;
  }
  return m_FirstSledge;  // Return the logical volume
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildSecondSledge() {
  if(!m_SecondSledge)
  {
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* SledgeMaterial = nistManager->FindOrBuildMaterial("G4_Al");

  double halfHeight = SecondSledgeThickness / 2;
  double Rmin1 = SecondSledgeOuterRadius;
  double Rmax1 = SecondSledgeOuterRadius;
  double Rmin2 = SecondSledgeInnerRadius;
  double Rmax2 = SecondSledgeOuterRadius;
  G4Cons* cutCylinder = new G4Cons("CutCylinder", Rmin1,Rmax1,Rmin2,Rmax2, halfHeight, 0.*deg, 360.*deg);
  G4LogicalVolume* logicSecondSledge = new G4LogicalVolume(cutCylinder, SledgeMaterial, "SecondSledgeLV");
  m_SecondSledge = logicSecondSledge;
  }
  return m_SecondSledge;  // Return the logical volume
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildFiller() {
  if(!m_Filler)
  {
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* FillerMaterial = nistManager->FindOrBuildMaterial("G4_Al");

  double halfHeight = FillerThickness / 2;
  G4Tubs* Cylinder = new G4Tubs("Cylinder", FillerInnerRadius, FillerOuterRadius, halfHeight, 0.*deg, 360.*deg);
  G4LogicalVolume* logicFiller = new G4LogicalVolume(Cylinder, FillerMaterial, "FillerLV");

  m_Filler = logicFiller;
  }
  return m_Filler;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* EdinburghDSSD::BuildBarHoler() {
  if(!m_substract)
  {
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* HolderMaterial = nistManager->FindOrBuildMaterial("G4_Al");


    G4VSolid* sBarMain = new G4Box("bar_main", BarHeight * 0.5, BarLength * 0.5, BarWidth * 0.5);


    G4Tubs* hole = new G4Tubs("hole", 0, FirstSledgeOuterRadius, BarWidth, 0.*deg, 360.*deg);




    //G4VSolid* hole = new G4Polycone("hole", 0., twopi, numHole, holeRadius, holeZ);

    G4VSolid* subtract = new G4SubtractionSolid("Bar-Hole", sBarMain, hole, 0, G4ThreeVector(0., 0.,0.0));
    m_substract = new G4LogicalVolume(subtract, HolderMaterial, "Subtract");

    m_substract->SetVisAttributes(G4Colour(0.2, 0.2, 0.2));
  }
    return m_substract;  // Return the logical volume
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
G4LogicalVolume* EdinburghDSSD::BuildBar(){
  if(!m_bar)
  {
  G4NistManager* nistManager = G4NistManager::Instance();
  G4Material* HolderMaterial = nistManager->FindOrBuildMaterial("G4_Al");

  G4VSolid* solidBar = new G4Box("solidBar1",
           Bar1Height*0.5,Bar1Length*0.5,Bar1Width*0.5);
  m_bar = new G4LogicalVolume(solidBar,           //shape
                          HolderMaterial,              //material
                          "Bar1");                 //
  }
  return m_bar;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void EdinburghDSSD::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("EdinburghDSSD");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " EDIN found " << endl;

  // Cartesian Case
  vector<string> cart = {"X1_Y1", "X1_Y16", "X16_Y1", "X16_Y16"};
  // Spherical Case
  vector<string> sphe = {"R", "THETA", "PHI", "BETA"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  EDIN Detector " << i + 1 << endl;
        TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
        TVector3 B = blocks[i]->GetTVector3("X16_Y1", "mm");
        TVector3 C = blocks[i]->GetTVector3("X1_Y16", "mm");
        TVector3 D = blocks[i]->GetTVector3("X16_Y16", "mm");

        G4ThreeVector g4_A(A.X(), A.Y(), A.Z());
        G4ThreeVector g4_B(B.X(), B.Y(), B.Z());
        G4ThreeVector g4_C(C.X(), C.Y(), C.Z());
        G4ThreeVector g4_D(D.X(), D.Y(), D.Z());

        // Now, you can use g4_A, g4_B, g4_C, and g4_D as G4ThreeVector objects
        // For example, to call AddDetector
  AddDetector(i, g4_A, g4_B, g4_C, g4_D);
    }

    else {
      cout << "ERROR: Missing token for EdinburghDSSD blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void EdinburghDSSD::ConstructDetector(G4LogicalVolume* world) {

  for (unsigned short i = 0; i < m_DetectorNumber.size(); i++) {
      G4RotationMatrix* rot = NULL;
      G4ThreeVector pos = G4ThreeVector(0, 0, 0);
      G4ThreeVector u = G4ThreeVector(0, 0, 0);
      G4ThreeVector v = G4ThreeVector(0, 0, 0);
      G4ThreeVector w = G4ThreeVector(0, 0, 0);
      G4ThreeVector Center = G4ThreeVector(0, 0, 0);
      // (u,v,w) unitary vector associated to telescope referencial
      // (u,v) // to silicon plan
      // w perpendicular to (u,v) plan and pointing ThirdStage
      u = (m_X1_Y16[i] + m_X16_Y16[i] - m_X1_Y1[i] - m_X16_Y1[i]);
      u = u.unit();

      v = m_X16_Y1[i] - m_X1_Y1[i];
      v = v.unit();

      w = u.cross(v);
      w = w.unit();

      Center = (m_X1_Y1[i] + m_X1_Y16[i] + m_X16_Y1[i] + m_X16_Y16[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      rot = new G4RotationMatrix(u, v, w);
      // translation to place Telescope
      pos = w * SiliconThickness * 0.5 + Center;
      new G4PVPlacement(G4Transform3D(*rot, pos), BuildSquareDetector(), "EDINSquare", world, false,m_DetectorNumber[i]);
      G4ThreeVector deadLayerPos = pos + w * (SiliconThickness * 0.5 + DeadLayerThickness * 0.5);
      new G4PVPlacement(G4Transform3D(*rot, deadLayerPos), BuildDeadLayer(), "EDINDeadLayer", world, false, m_DetectorNumber[i]);


    }


    G4ThreeVector global_shift = G4ThreeVector(-1.9,35.5,-0.99);

    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO
    //                    Construct Nozzle
    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO

    G4RotationMatrix* yRotNoz = new G4RotationMatrix();
    yRotNoz->rotateZ(180.0 * deg);
    yRotNoz->rotateX(90*deg);
    G4double shift_nozzle = 66.396 - 21.694;

    G4ThreeVector posNoz = G4ThreeVector(0, shift_nozzle, 0.);
    posNoz += global_shift;
    new G4PVPlacement(yRotNoz, posNoz, BuildNozzle(), "Nozzle", world, false, 0);

    G4ThreeVector positionTarget = G4ThreeVector(0, 0, 0);
    positionTarget +=global_shift;

    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO
    //                    Construct subWorld
    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO

    G4RotationMatrix* yRot = new G4RotationMatrix();
    yRot->rotateY(90.*deg);
    yRot->rotateZ(FoilAngle+90*deg);

    G4LogicalVolume* subWorldLV = BuildSubWorld();
    G4ThreeVector positionSubWorld = G4ThreeVector(0, 0, 0);
    positionSubWorld+=global_shift;
    new G4PVPlacement(yRot, positionSubWorld,subWorldLV,"SubWorld",world,false,0);

    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO
    //                    Construct Foil
    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO

    G4LogicalVolume* foilLV = BuildFoil();

    G4RotationMatrix* Rot_foil = new G4RotationMatrix();
    Rot_foil->rotateY(90. * deg);

    G4ThreeVector foilPosition = G4ThreeVector(0., 0, 0.);
    new G4PVPlacement(Rot_foil, foilPosition, foilLV, "Foil", BuildSubWorld(), false, 0);


    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO
    //                Construct Foil Holder
    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO

    G4RotationMatrix* rotHolder = new G4RotationMatrix();
    rotHolder->rotateY(90.*deg);
    double offsetHolder = FoilHolderThickness/2 + FoilThickness/2;
    double offsetFirstSledge = -FirstSledgeThickness/2 - FoilThickness/2;
    double offsetSecondSledge = -FirstSledgeThickness - FoilThickness/2 - SecondSledgeThickness/2;
    double offsetFiller = 0;
    new G4PVPlacement(rotHolder, G4ThreeVector(offsetHolder,0., 0.), BuildFoilHolder(), "Holder", BuildSubWorld(), false, 0);
    new G4PVPlacement(rotHolder, G4ThreeVector(offsetFirstSledge,0., 0.), BuildFirstSledge(), "FirstSledge", BuildSubWorld(), false, 0);
    new G4PVPlacement(rotHolder, G4ThreeVector(offsetSecondSledge ,0., 0.), BuildSecondSledge(), "SecondSledge", BuildSubWorld(), false, 0);
    new G4PVPlacement(rotHolder, G4ThreeVector(offsetFiller ,0., 0.), BuildFiller(), "Filler", BuildSubWorld(), false, 0);

    G4RotationMatrix* rotSub = new G4RotationMatrix();
    rotSub->rotateY(90.*deg);
    rotSub->rotateZ(90.*deg);


    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO
    //                Construct bar holder
    // OOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOOOOoooOOOooOOOO

    //G4ThreeVector(0.,0.,FirstSledgeThickness- FoilThickness/2),
    new G4PVPlacement(rotSub,                      //no rotation
             G4ThreeVector(-BarWidth/2+ FoilThickness/2,0.,0.),
             BuildBarHoler(),                //logical volume
             "BarHolder",                    //name
             BuildSubWorld(),                      //mother  volume
             false,                       //no boolean operation
               0);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void EdinburghDSSD::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("EdinburghDSSD")){
    pTree->Branch("EdinburghDSSD", "TEdinburghDSSDData", &m_Event) ;
  }
  pTree->SetBranchAddress("EdinburghDSSD", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EdinburghDSSD::ProcessStrip(std::map<unsigned int, std::pair<double, double>>& map, int key, double energy, double time) {
    auto& entry = map[key];
    if (entry.second == 0) {  // Check if this is the first time this strip is hit
        entry = std::make_pair(energy, time);
        //std::cout << "Initial hit processed for key " << key << ": energy = " << energy << ", time = " << time << std::endl;
    } else {
        entry.first += energy;  // Accumulate energy
        //std::cout << "Energy accumulated for key " << key << ": total energy = " << entry.first << std::endl;
    }
}

void EdinburghDSSD::UpdateMap(bool X_Y,std::map<unsigned int, std::pair<double, double>>& map,
                       std::map<unsigned int, bool>& mapInterstrip, int key, double energy,
                       double time, double t0, bool interstrip) {
    // X_Y = 0 for X, 1 for Y.
    int det = key / static_cast<int>(1e6);
    int strip = key % static_cast<int>(1e6);

    // Check if the key already exists in the map
    auto it = map.find(key);
    if (it == map.end()) {
        // If key does not exist, initialize it directly
        map[key] = std::make_pair(energy, time);
        mapInterstrip[key] = interstrip;
        //std::cout << "New entry added for (det,strip) = (" << det << "," << strip << ") with energy " << energy << " and time " << time << std::endl;
        return; // Exit the function after setting new data
    }

    // Key exists, process based on time difference
    double last_energy = it->second.first;
    double last_time = it->second.second;

    //std::cout << "In UpdateMap, (det,strip) = (" << det << "," << strip << ")";
    //std::cout << "∆T = " << time - last_time << std::endl;

    if (time - last_time < T_sum) {
        //std::cout << "Adding energy " << energy << std::endl;
        it->second.first += energy;  // Accumulate energy if within time sum
        if (interstrip) {
            mapInterstrip[key] = true;
        }
        if(it->second.first > 8)
        {
          std::cout << "UPDATE MAP : energy of event: " << it->second.first ;
          std::cout << ", (det,strip) = (" << det << "," << strip << ")";

        }
        //std::cout << "After adding energy " << map[key].first << std::endl;
    } else {
        // Set event data if time difference is greater than T_sum

        //std::cout << "Push event " << std::endl;
        //std::cout << "energy of event: " << last_energy << std::endl;
        if(X_Y==0)
        {
        m_Event->SetDSSDX_Tstamp(t0);
        m_Event->SetDSSDX_Interstrip(mapInterstrip[key]);
        m_Event->SetDSSDXE(false, det, strip, last_energy);
        m_Event->SetDSSDXT(false, det, strip, last_time);
        }
        else
        {
          m_Event->SetDSSDY_Tstamp(t0);
          m_Event->SetDSSDY_Interstrip(mapInterstrip[key]);
          m_Event->SetDSSDYE(false, det, strip, last_energy);
          m_Event->SetDSSDYT(false, det, strip, last_time);
        }
        // Update the map with the new energy and current time
        map[key] = std::make_pair(energy, time);
        mapInterstrip[key] = interstrip;
    }
}


// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void EdinburghDSSD::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  G4double t0 =  RandFlat::shoot()*1;
  //std::cout << "T0 = " << t0 << std::endl;
  ///////////
  // Square
  DSSDScorers::PS_Images* SiScorer = (DSSDScorers::PS_Images*)m_SquareScorer->GetPrimitive(0);

  // Loop on the Square map


  unsigned int sizeFront = SiScorer->GetFrontMult();
  unsigned int sizeBack = SiScorer->GetBackMult();

  std::map<unsigned int, std::pair<double, double>> mapFront;
  std::map<unsigned int, std::pair<double, double>> mapBack;

  std::map<unsigned int, std::pair<double, double>> mapTempo;

  std::map<unsigned int, bool> mapInterstripFront;
  std::map<unsigned int, bool> mapInterstripBack;
  std::set<std::pair<int, int>> interstripCouples;  // Stores interstrip pairs


  std::map<unsigned int, std::pair<double, double>>::iterator it;


  int strip1,strip2;
  double E_b, E_g;
  int refTrack;
  unsigned int j, a,r, g, b;

  double energy,time;
  int det, b_key, g_key;
  int key, partner_key,partner_det, partner_strip;

  double E_tot,rand, E1,E2;

  for (unsigned int i = 0; i < sizeFront; i++) {
      refTrack = SiScorer->GetTrackId(i);
      //std::cout << std::endl  << "i = " << i << std::endl;
      //std::cout << "Processing Track ID = " << refTrack << std::endl;

      j = i;

      while (j < sizeFront && SiScorer->GetTrackId(j) == refTrack) {
          energy = SiScorer->GetEnergyFront(j);
          det = SiScorer->GetDetectorFront(j);
          time = SiScorer->GetTimeFront(j);
          SiScorer->GetARGBFront(j, a, r, g, b);

          b_key = b + det * 1e6;
          g_key = g + det * 1e6;
          if (r == 0) { // no interstrip
              ProcessStrip(mapTempo,b_key, energy, time);
          } else { // interstrip
              if (g > 0 && g < 17 && b > 0 && b < 17) { // both in detector
                  interstripCouples.insert(std::make_pair(b_key, g_key));
                  ProcessStrip(mapTempo,b_key, energy, time);
                  ProcessStrip(mapTempo,g_key, 0, time); // To set Time
              }
              else if ((g<= 0 || g >= 17)&&(b>0 && b < 17))
              {
                ProcessStrip(mapTempo,b_key, energy, time);
              }
              else if ((b<= 0 || b >= 17) && (g>0 && g < 17))
              {
                ProcessStrip(mapTempo,g_key, energy, time);
              }
          }
          j++;
      }

      std::set<int> treatedKeys; // Tracks which keys have been processedx

      // Processing entries in mapTempo
      for (const auto& entry : mapTempo) {
          key = entry.first;
          if (treatedKeys.find(key) != treatedKeys.end()) {
              continue; // Skip this key if it has already been processed
          }

          auto it = std::find_if(interstripCouples.begin(), interstripCouples.end(),
              [key](const std::pair<int, int>& pair) {
                  return pair.first == key || pair.second == key;
              }); // function to find if pair in couples

          if (it != interstripCouples.end()) { // It is an interstrip couple

              partner_key = (it->first == key) ? it->second : it->first;
              E_tot = entry.second.first + mapTempo[partner_key].first;
              rand = G4UniformRand();
              E1 = rand * E_tot;
              E2 = E_tot - E1;
              //std::cout << "IS || Track ID = " << refTrack << ", ";
              //std::cout << "(det,strip) = (" << det << "," << strip << ")";
              //std::cout << " || PARTNER = (" << partner_det << "," << partner_strip << ")" ;
              //std::cout << "E = " << entry.second.first << ", E(part) =" <<  mapTempo[partner_key].first << std::endl;
              //std::cout << "T = " << entry.second.second << ", T(part) =" <<  mapTempo[partner_key].second << std::endl;
              if (E1 < E2) {
                  E1 *= -1;
              } else {
                  E2 *= -1;
              }

              UpdateMap(0,mapFront, mapInterstripFront, key, E1,mapTempo[key].second, t0, true);
              UpdateMap(0,mapFront, mapInterstripFront, partner_key, E2,mapTempo[partner_key].second, t0, true);

              treatedKeys.insert(key);
              treatedKeys.insert(partner_key);

              //std::cout << "Processed interstrip pair: det " << det << ", strip " << strip
              //          << " and det " << partner_det << ", strip " << partner_strip
              //          << " with energies: " << E1 << ", " << E2 << std::endl;
          } else {
              //int det = key / static_cast<int>(1e6);
              //int strip = key % static_cast<int>(1e6);
              //std::cout << "Not IS (det,strip) = " << det << "," << strip ;
              //std::cout << ", Energy = " << entry.second.first << std::endl;
              UpdateMap(0,mapFront, mapInterstripFront,key,
                entry.second.first, entry.second.second,t0, false);
              //std::cout << "Transferred non-interstrip det " << det << ", strip " << strip << " with energy " << entry.second.first << std::endl;
          }
      }

      mapTempo.clear(); // Clear temporary storage after processing
      interstripCouples.clear();
      //std::cout << "Cleared temporary storage for next track." << std::endl;
      //std::cout << "After while j = " << j  << std::endl;
      i = j-1;
  }


  double energyX, timeX;
  unsigned int detX, stripX;
  bool bool_interstrip;
  for (it = mapFront.begin(); it != mapFront.end(); it++) {
    energyX = RandGauss::shoot(it->second.first, SigmaEnergy);
    timeX = RandGauss::shoot(it->second.second, SigmaTime);
    bool_interstrip = mapInterstripFront[it->first];
    stripX = it->first - 1000000 * (it->first / 1000000);
    detX = it->first / 1000000;

    m_Event->SetDSSDX_Tstamp(t0);
    m_Event->SetDSSDX_Interstrip(bool_interstrip);
    m_Event->SetDSSDXE(false,detX, stripX,energyX);
    m_Event->SetDSSDXT(false,detX, stripX,timeX);
  }

  mapTempo.clear();

  for (unsigned int i = 0; i < sizeBack; i++) {
      refTrack = SiScorer->GetTrackId(i);
      //std::cout << std::endl  << "i = " << i << std::endl;
      //std::cout << "Processing Track ID = " << refTrack << std::endl;

      j = i;

      while (j < sizeBack && SiScorer->GetTrackId(j) == refTrack) {
          energy = SiScorer->GetEnergyBack(j);
          det = SiScorer->GetDetectorBack(j);
          time = SiScorer->GetTimeBack(j);
          SiScorer->GetARGBBack(j, a, r, g, b);

          b_key = b + det * 1e6;
          g_key = g + det * 1e6;

          if (r == 0) { // no interstrip
              ProcessStrip(mapTempo,b_key, energy, time);
          } else { // interstrip
              if (g > 0 && g < 17 && b > 0 && b < 17) { // both in detector
                  interstripCouples.insert(std::make_pair(b_key, g_key));
                  ProcessStrip(mapTempo,b_key, energy, time);
                  ProcessStrip(mapTempo,g_key, 0, time); // To set Time
              }
              else if ((g<= 0 || g >= 17)&&(b>0 && b < 17))
              {
                ProcessStrip(mapTempo,b_key, energy, time);
              }
              else if ((b<= 0 || b >= 17) && (g>0 && g < 17))
              {
                ProcessStrip(mapTempo,g_key, energy, time);
              }
          }
          j++;
      }

      std::set<int> treatedKeys; // Tracks which keys have been processedx

      // Processing entries in mapTempo
      for (const auto& entry : mapTempo) {
          int key = entry.first;
          if (treatedKeys.find(key) != treatedKeys.end()) {
              continue; // Skip this key if it has already been processed
          }

          auto it = std::find_if(interstripCouples.begin(), interstripCouples.end(),
              [key](const std::pair<int, int>& pair) {
                  return pair.first == key || pair.second == key;
              }); // function to find if pair in couples

          if (it != interstripCouples.end()) { // It is an interstrip couple

              partner_key = (it->first == key) ? it->second : it->first;
              E_tot = entry.second.first + mapTempo[partner_key].first;
              rand = G4UniformRand();
              E1 = rand * E_tot;
              E2 = E_tot - E1;

              UpdateMap(1,mapBack, mapInterstripBack, key, E1,mapTempo[key].second, t0, true);
              UpdateMap(1,mapBack, mapInterstripBack, partner_key, E2,mapTempo[partner_key].second, t0, true);

              treatedKeys.insert(key);
              treatedKeys.insert(partner_key);

              //std::cout << "Processed interstrip pair: det " << det << ", strip " << strip
              //          << " and det " << partner_det << ", strip " << partner_strip
              //          << " with energies: " << E1 << ", " << E2 << std::endl;
          } else {
              UpdateMap(1,mapBack, mapInterstripBack,key,
              entry.second.first, entry.second.second,t0, false);
              //std::cout << "Transferred non-interstrip det " << det << ", strip " << strip << " with energy " << entry.second.first << std::endl;
          }
      }

      mapTempo.clear(); // Clear temporary storage after processing
      interstripCouples.clear();
      //std::cout << "Cleared temporary storage for next track." << std::endl;
      //std::cout << "After while j = " << j  << std::endl;
      i = j-1;
  }

  double energyY, timeY;
  unsigned int stripY,detY;


  for (it = mapBack.begin(); it != mapBack.end(); it++) {
      energyY = RandGauss::shoot(it->second.first, SigmaEnergy);
      timeY = RandGauss::shoot(it->second.second, SigmaTime);
      bool_interstrip = mapInterstripBack[it->first];
      stripY = it->first - 1000000 * (it->first / 1000000);
      detY = it->first / 1000000;
      m_Event->SetDSSDY_Tstamp(t0);
      m_Event->SetDSSDY_Interstrip(bool_interstrip);
      m_Event->SetDSSDYE(false,detY, stripY,energyY);
      m_Event->SetDSSDYT(false,detY, stripY,timeY);
    }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void EdinburghDSSD::InitializeScorers() {
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  m_SquareScorer = CheckScorer("SquareScorer", already_exist);

  if (already_exist)
    return;

  string nptool = getenv("NPTOOL");
  G4VPrimitiveScorer* SiScorer = new DSSDScorers::PS_Images(
        "SquareScorer", nptool + "/NPLib/Detectors/EdinburghDSSD/ressources/maskFront.png",
        nptool + "/NPLib/Detectors/EdinburghDSSD/ressources/maskBack.png", 49.50 / 12800, 49.50 / 12800, 0, 0, 0xffff0000, 1, true);

  //G4VPrimitiveScorer* InterScorer = new InteractionScorers::PS_Interactions("SiScorer", ms_InterCoord, 0);

  m_SquareScorer->RegisterPrimitive(SiScorer);
  //m_SquareScorer->RegisterPrimitive(InterScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SquareScorer);


}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* EdinburghDSSD::Construct(){
  return  (NPS::VDetector*) new EdinburghDSSD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_EdinburghDSSD{
    public:
      proxy_nps_EdinburghDSSD(){
        NPS::DetectorFactory::getInstance()->AddToken("EdinburghDSSD","EdinburghDSSD");
        NPS::DetectorFactory::getInstance()->AddDetector("EdinburghDSSD",EdinburghDSSD::Construct);
      }
  };

  proxy_nps_EdinburghDSSD p_nps_EdinburghDSSD;
}
