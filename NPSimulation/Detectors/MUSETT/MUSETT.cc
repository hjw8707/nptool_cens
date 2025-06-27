/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr*
 *                                                                           *
 * Creation Date  : January 2024                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  MUSETT simulation for double alpha decay            *
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
#include "MUSETT.hh"
#include "MUSETTMap.h"
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
using namespace MUSETT_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// MUSETT Specific Method
MUSETT::MUSETT() {
  m_Event = new TMUSETTData();
  InitializeMaterial();
  m_SquareScorer = 0;
  m_SquareDetector = 0;
  m_SquareDeadLayer =0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 0, 1, 0.8));

  m_VisDeadLayer = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));
}

MUSETT::~MUSETT() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void MUSETT::AddDetector(int DetectorNumber, G4ThreeVector PX1_Y1, G4ThreeVector PX1_Y128,G4ThreeVector PX128_Y1, G4ThreeVector PX128_Y128) {
  m_X1_Y1.push_back(PX1_Y1);         // Top Left Corner Position Vector
  m_X1_Y128.push_back(PX1_Y128);     // Bottom Left Corner Position Vector
  m_X128_Y1.push_back(PX128_Y1);     // Bottom Right Corner Position Vector
  m_X128_Y128.push_back(PX128_Y128); // Center Corner Position Vector
  m_DetectorNumber.push_back(DetectorNumber);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* MUSETT::BuildSquareDetector() {
  if (!m_SquareDetector) {
    G4String Name = "MUSETT";


    G4Box* solidSquare = new G4Box(Name, 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * SiliconThickness);

    G4LogicalVolume* logicSquare = new G4LogicalVolume(
    solidSquare, m_MaterialVacuum, Name, 0, 0, 0);


    logicSquare->SetVisAttributes(m_VisSquare);

    G4ThreeVector positionFirstStage = G4ThreeVector(0, 0, 0);

    G4Box* solidFirstStage = new G4Box("solidFirstStage", 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * SiliconThickness);
    G4LogicalVolume* logicFirstStage = new G4LogicalVolume(
        solidFirstStage, m_MaterialSilicon, "logicFirstStage", 0, 0, 0);


    logicFirstStage->SetVisAttributes(m_VisSquare);
    new G4PVPlacement(0, positionFirstStage, logicFirstStage, Name + "_FirstStage", logicSquare, false, 0);

    m_SquareDetector = logicSquare;
    // Set First Stage sensible
    logicFirstStage->SetSensitiveDetector(m_SquareScorer);

    G4UserLimits* stepLimit = new G4UserLimits(1*um);

    logicFirstStage->SetUserLimits(stepLimit);

  }

  return m_SquareDetector;
}

G4LogicalVolume* MUSETT::BuildDeadLayer(){
  if (!m_SquareDeadLayer) {
    G4String Name = "MUSETT";

    // Step 1: Create the Square Solid
    G4Box* solidSquare = new G4Box(Name, 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * DeadLayerThickness);
    G4LogicalVolume* logicSquare = new G4LogicalVolume(
      solidSquare, m_MaterialVacuum, Name, 0, 0, 0);

    // Step 2: Create the Boron Dead Layer Solid
    G4Box* solidDeadLayer = new G4Box("solidDeadLayer", 0.5 * SquareLength, 0.5 * SquareLength, 0.5 * DeadLayerThickness);


    G4LogicalVolume* logicDeadLayer = new G4LogicalVolume(solidDeadLayer, m_MaterialBoron, "logicDeadLayer");

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


G4LogicalVolume* MUSETT::BuildTarget() {
  // Define target properties
  G4String Name = "FoilTarget";


  // Create the solid foil target (cylinder)
 solidTarget = new G4Tubs(Name, 0.0, targetRadius, targetHalfThickness, 0.0, 360.0 * deg);


 // Create the logical volume for the foil target
 G4LogicalVolume* logicTarget = new G4LogicalVolume(solidTarget, m_MaterialCarbon, Name);

 // Set the visualization attributes to make it red
 G4VisAttributes* targetVisAtt = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0)); // Red color
 logicTarget->SetVisAttributes(targetVisAtt);

 return logicTarget;
}

G4LogicalVolume* MUSETT::BuildTargetHolder() {
  // Define TargetHolder properties
  G4String Name = "TargetHolder";


  // Position of the lower part

  // Create the solid TargetHolder (box)
  G4Box* solidTargetHolder = new G4Box(Name,0.5 * holderLength, 0.5 * holderWidth , 0.5 * holderHeight);

G4ThreeVector holePosition(0, 0, -holderHeight/2+lowerPartZPosition);
//G4ThreeVector holePosition(0, 0, 0);

// Create a tubular hole centered at the origin with a specified radius
G4Tubs* solidHole = new G4Tubs("Hole", 0, targetRadius, holderWidth, 0, 360 * degree);

// Create a rotation matrix to rotate the hole by 90 degrees around the x-axis
G4RotationMatrix* rotation = new G4RotationMatrix();
rotation->rotateX(90 * degree);

// Apply the rotation to the solidHole
G4Transform3D transform(*rotation, holePosition);

// Create a subtraction solid to remove the hole from the TargetHolder
G4SubtractionSolid* solidTargetHolderWithHole = new G4SubtractionSolid("TargetHolderWithHole", solidTargetHolder, solidHole, transform);

//G4SubtractionSolid* solidTargetHolderWithHole = new G4SubtractionSolid("TargetHolderWithHole", solidTargetHolder, solidHole);

// Create a logical volume for the TargetHolder
G4LogicalVolume* logicTargetHolderWithHole = new G4LogicalVolume(solidTargetHolderWithHole, m_MaterialAluminium, "TargetHolderWithHole");




  // Create the logical volume for the TargetHolder
  //G4LogicalVolume* logicTargetHolder = new G4LogicalVolume(solidTargetHolder, holderMaterial, Name);


  // Set the visualization attributes to make it light gray
  G4VisAttributes* holderVisAtt = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8)); // Light gray color
  logicTargetHolderWithHole->SetVisAttributes(holderVisAtt);

  // Position the TargetHolder // Replace "yourMotherLogicalVolume" with the actual parent logical volume

  return logicTargetHolderWithHole;
}

G4LogicalVolume* MUSETT::BuildAlLattice() {

    // Define dimensions of Xstick and Ystick

    // Create Xstick and Ystick solids
    G4Box* solidXstick = new G4Box("Xstick", 0.5 * XstickWidth, 0.5 * XstickLength, 0.5 * XstickHeight);
    G4Box* solidYstick = new G4Box("Ystick", 0.5 * YstickWidth, 0.5 * YstickLength, 0.5 * YstickHeight);
    G4Box* solidBondPad = new G4Box("BondPad", 0.5 * BondPadSize, 0.5 * BondPadSize, 0.5 * BondPadHeight);

    // Create logical volume for Xstick and Ystick
    G4LogicalVolume* logicXstick = new G4LogicalVolume(solidXstick, m_MaterialAluminium, "Xstick");
    G4LogicalVolume* logicYstick = new G4LogicalVolume(solidYstick, m_MaterialAluminium, "Ystick");
    G4LogicalVolume* logicBondPad = new G4LogicalVolume(solidBondPad, m_MaterialAluminium, "BondPad");

    // Set visualization attributes
    //logicXstick->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 0.0, 0.0))); // Red color for Xsticks
    //logicYstick->SetVisAttributes(new G4VisAttributes(G4Colour(0.0, 0.0, 1.0))); // Blue color for Ysticks

    // Define lattice parameters


    //int numXsticks = static_cast<int> (SquareLength / (XstickWidth + spacingX)) ;

    //G4int numXsticks = 6;

    //G4int numYsticks = 2;
    //G4cout << '@@@@@@@'<< G4endl;
    //G4cout << "number X sticks = " << numXsticks << G4endl;
    //G4cout << "number Y sticks = " << numYsticks << G4endl;

    // Calculate total width of Xsticks and total length of Ysticks
    const G4double totalXstickWidth = SquareLength ;
    const G4double totalYstickLength =  SquareLength ;
    // Create lattice volume
    G4Box* solidLattice = new G4Box("Lattice", totalXstickWidth/2, totalYstickLength/2, XstickHeight/2);
    G4LogicalVolume* logicLattice = new G4LogicalVolume(solidLattice, m_MaterialVacuum, "Lattice");


    // Placement of Ysticks in the lattice
    const G4double x0 = -totalXstickWidth/2 + YstickWidth/2;
    const G4double x_up_down = XstickWidth + YstickWidth;
    const G4double x_down_down = x_up_down + InterStrip + YstickWidth;
    const G4double y0_y = 0;
    const G4double y0_x = -totalYstickLength/2;
    const G4double step_y = 2890/3*um + XstickLength;
    G4int cnt_y = 0;
    G4int cnt_x = 0;
    G4int cnt_pad = 0;
    for (G4int i = 0; i < numYsticks ; i++) {
    //for (G4int i = 0; i < 2 ; i++) {
            G4double xDown = x0 + i*x_down_down ;
            G4double Xup = xDown +x_up_down;
            new G4PVPlacement(0, G4ThreeVector(xDown, y0_y, 0), logicYstick, "Ystick", logicLattice, false, cnt_y++);
            new G4PVPlacement(0, G4ThreeVector(Xup, y0_y, 0), logicYstick, "Ystick", logicLattice, false, cnt_y++);
            G4double xPos = (xDown + Xup)/2;

            // First and last sticks are different from the others
            G4double yPos = y0_x +XstickLength/2;//;+ 550*um + XstickLength/2;
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos, 0), logicXstick, "Xstick", logicLattice, false, cnt_x++);
            yPos += 550*um + XstickLength/2;
            G4double new_y0 = yPos+ XstickLength/2;
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos, 0), logicXstick, "Xstick", logicLattice, false, cnt_x++);
            G4double y =0;

            for (G4int j = 1; j < numXsticks+1; j++)
            {
              G4double yPos = new_y0 + j * step_y;
              new G4PVPlacement(0, G4ThreeVector(xPos, yPos, 0), logicXstick, "Xstick", logicLattice, false, cnt_x++);
            }

            yPos = y0_x + totalYstickLength - XstickLength/2;
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos, 0), logicXstick, "Xstick", logicLattice, false, cnt_x++);
            G4double YposPadLeft = y0_x + XstickLength +BondPadSize/2;
            G4double YposPadRight = y0_x + totalYstickLength - XstickLength  -BondPadSize/2;
            new G4PVPlacement(0, G4ThreeVector(xPos, YposPadLeft, 0), logicBondPad, "BondPad", logicLattice, false, cnt_pad++);
            new G4PVPlacement(0, G4ThreeVector(xPos, YposPadRight, 0), logicBondPad, "BondPad", logicLattice, false, cnt_pad++);

    }


    G4VisAttributes* latticeVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.1)); // White color with 50% transparency
    latticeVisAtt->SetForceSolid(false); // Allow transparency
    logicLattice->SetVisAttributes(latticeVisAtt);
    return logicLattice;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void MUSETT::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("MUSETT");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " MUSETT found " << endl;

  // Cartesian Case
  vector<string> cart = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128"};
  // Spherical Case
  vector<string> sphe = {"R", "THETA", "PHI", "BETA"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  MUSETT Detector " << i + 1 << endl;
        TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
        TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
        TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
        TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");

        G4ThreeVector g4_A(A.X(), A.Y(), A.Z());
        G4ThreeVector g4_B(B.X(), B.Y(), B.Z());
        G4ThreeVector g4_C(C.X(), C.Y(), C.Z());
        G4ThreeVector g4_D(D.X(), D.Y(), D.Z());

        // Now, you can use g4_A, g4_B, g4_C, and g4_D as G4ThreeVector objects
        // For example, to call AddDetector
  AddDetector(i, g4_A, g4_B, g4_C, g4_D);
    }

    else {
      cout << "ERROR: Missing token for MUSETT blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After Detecor
// onstruction::AddDetector Method
void MUSETT::ConstructDetector(G4LogicalVolume* world) {
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
      u = (m_X1_Y128[i] + m_X128_Y128[i] - m_X1_Y1[i] - m_X128_Y1[i]);
      u = u.unit();

      v = m_X128_Y1[i] - m_X1_Y1[i];
      v = v.unit();

      w = u.cross(v);
      w = w.unit();

      Center = (m_X1_Y1[i] + m_X1_Y128[i] + m_X128_Y1[i] + m_X128_Y128[i]) / 4;

      // Passage Matrix from Lab Referential to Telescope Referential
      rot = new G4RotationMatrix(u, v, w);
      // translation to place Telescope
      pos = w * SiliconThickness * 0.5 + Center;


      new G4PVPlacement(G4Transform3D(*rot, pos), BuildSquareDetector(), "MUSETTSquare", world, false,                  m_DetectorNumber[i]);





      G4ThreeVector deadLayerPos = pos + w * (SiliconThickness * 0.5 + DeadLayerThickness * 0.5);
      new G4PVPlacement(G4Transform3D(*rot, deadLayerPos), BuildDeadLayer(), "MUSETTDeadLayer", world, false, m_DetectorNumber[i]);


      // Define the rotation matrix for rotating the lattice by 90 degrees around the z-axis
      G4RotationMatrix* rotation90Z = new G4RotationMatrix();
      //rotation90Z->rotateZ(90 * degree);

      // Combine the rotations
      G4RotationMatrix* combinedRotation = new G4RotationMatrix(*rotation90Z);
      combinedRotation->transform(*rot);

      // Apply the combined rotation
      G4ThreeVector LatticePos = deadLayerPos + w * (DeadLayerThickness * 0.5 + YstickHeight*0.5);
      new G4PVPlacement(G4Transform3D(*combinedRotation, LatticePos), BuildAlLattice(), "MUSETT_lattice", world, false, m_DetectorNumber[i]);

    }


    double rotationAngle = 0.00000000 * deg;
    G4RotationMatrix* TargetAngle = new G4RotationMatrix();
    TargetAngle->rotateZ(rotationAngle);

    double offsetTargetX = -3.00000000 * mm;
    double offsetTargetY = 6.00000000 * mm;
    //offsetTargetY += targetHalfThickness;
    G4ThreeVector positionTarget = G4ThreeVector(offsetTargetX, offsetTargetY, 0);


    G4RotationMatrix* initialTarget = new G4RotationMatrix();
    initialTarget->rotateX(90.0 * deg);

    G4RotationMatrix* TargetTotalRotation = new G4RotationMatrix(*TargetAngle);
    TargetTotalRotation->transform(*initialTarget);

    new G4PVPlacement(TargetTotalRotation, positionTarget, BuildTarget(), "MUSETT_Target", world, false, 0);

    double holder_shift = 0.50000000 * mm;
    double offsetTargetHolderX = holder_shift * TMath::Sin(rotationAngle);
    double offsetTargetHolderY = holder_shift * TMath::Cos(rotationAngle);
    G4ThreeVector positionTargetHolder = G4ThreeVector(offsetTargetX-offsetTargetHolderX,offsetTargetY-offsetTargetHolderY, -lowerPartZPosition + 0.5 * holderHeight);

    new G4PVPlacement(TargetAngle, positionTargetHolder, BuildTargetHolder(), "Musett_TargetHolder", world, false, 0);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void MUSETT::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("MUSETT")) {
    pTree->Branch("MUSETT", "TMUSETTData", &m_Event);
  }
  pTree->SetBranchAddress("MUSETT", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
void MUSETT::ProcessStrip(std::map<unsigned int, std::pair<double, double>>& map, int key, double energy, double time) {
    auto& entry = map[key];
    if (entry.second == 0) {  // Check if this is the first time this strip is hit
        entry = std::make_pair(energy, time);
        //std::cout << "Initial hit processed for key " << key << ": energy = " << energy << ", time = " << time << std::endl;
    } else {
        entry.first += energy;  // Accumulate energy
        //std::cout << "Energy accumulated for key " << key << ": total energy = " << entry.first << std::endl;
    }
}

void MUSETT::UpdateMap(bool X_Y,std::map<unsigned int, std::pair<double, double>>& map,
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
          std::cout << std::endl;

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



// Called at in the EventAction::EndOfEventAvtion
void MUSETT::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  G4double t0 =  RandFlat::shoot()*1;
  //std::cout << "✅ New event" << std::endl;
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
              if (g > 0 && g < 129 && b > 0 && b < 129) { // both in detector
                  interstripCouples.insert(std::make_pair(b_key, g_key));
                  ProcessStrip(mapTempo,b_key, energy, time);
                  ProcessStrip(mapTempo,g_key, 0, time); // To set Time
              }
              else if ((g<= 0 || g >= 129)&&(b>0 && b < 129))
              {
                ProcessStrip(mapTempo,b_key, energy, time);
              }
              else if ((b<= 0 || b >= 129) && (g>0 && g < 129))
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
              if (g > 0 && g < 129 && b > 0 && b < 129) { // both in detector
                  interstripCouples.insert(std::make_pair(b_key, g_key));
                  ProcessStrip(mapTempo,b_key, energy, time);
                  ProcessStrip(mapTempo,g_key, 0, time); // To set Time
              }
              else if ((g<= 0 || g >= 129)&&(b>0 && b < 129))
              {
                ProcessStrip(mapTempo,b_key, energy, time);
              }
              else if ((b<= 0 || b >= 129) && (g>0 && g < 129))
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
void MUSETT::InitializeScorers() {
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  m_SquareScorer = CheckScorer("MUSETT", already_exist);

  if (already_exist)
    return;

  string nptool = getenv("NPTOOL");
  G4VPrimitiveScorer* SiScorer = new DSSDScorers::PS_Images(
        "SquareScorer", nptool + "/NPLib/Detectors/MUSETT/ressources/maskFront.png",
        nptool + "/NPLib/Detectors/MUSETT/ressources/maskBack.png", 97.22 / 12800, 97.22 / 12800, 0, 0, 0xffff0000, 1,true);

  G4VPrimitiveScorer* InterScorer = new InteractionScorers::PS_Interactions("SiScorer", ms_InterCoord, 0);

  m_SquareScorer->RegisterPrimitive(SiScorer);
  m_SquareScorer->RegisterPrimitive(InterScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SquareScorer);



}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void MUSETT::InitializeMaterial(){
  m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_Si");
  m_MaterialAluminium = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_Al");
  m_MaterialBoron = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_B");
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_MaterialCarbon = MaterialManager::getInstance()->GetMaterialFromLibrary("G4_C");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* MUSETT::Construct() { return (NPS::VDetector*)new MUSETT(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_MUSETT {
 public:
  proxy_nps_MUSETT() {
    NPS::DetectorFactory::getInstance()->AddToken("MUSETT", "MUSETT");
    NPS::DetectorFactory::getInstance()->AddDetector("MUSETT", MUSETT::Construct);
  }
};

proxy_nps_MUSETT p_nps_MUSETT;
}
