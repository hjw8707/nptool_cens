/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emmanuel Rey-herme                                       *
 * contact address: marine.vandebrouck@cea.fr                                *
 *                                                                           *
 * Creation Date  : septembre 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  SEASON simulation                             *
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
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"
#include "Randomize.hh"

// NPTool header
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"
#include "SEASON.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

#include "TMath.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SEASON_NS {
  // Energy and time Resolution
  const double ResoTime = 8.5 * ns;
  const double ResoAlpha = 6.37 * keV;    // Equivalent to 15 keV FWHM
  const double ResoElectron = 2.97 * keV; // Equivalent to 7 keV FWHM;

  //
  const double Thickness = 1 * mm;             // Thickness of DSSDs
  const double FoilWheelCenterDist = 16. * cm; // distance from foil center to wheel center
  const double WheelRadius = 11.5 * cm;
  const double WheelThickness = 2 * mm;
  const double FinThickness = 0.5 * mm;
  const double FinLenght = 6. * cm;
  const double FinWidth = 3. * cm;
  const int FoilNbr = 11;
  const double DistDSSD1_AlWindow = 9 * mm;

  const double StripWidth = 2 * mm;
  const double StripPitch = 1.925 * mm;
  const double StripLength = 63.96 * mm;
  const double DetSize = 67.975 * mm;
  const int NumberOfStripsX = 32;
  const int NumberOfStripsY = 32;

  const double DeadLayerThickness = 50 * nm;
  const double GridThickness = 300 * nm;
  const double GridWidth = 30 * um;

  G4Material* FoilMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("C");
  G4Material* DeadLayerMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* Grid_Material = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
} // namespace SEASON_NS

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// SEASON Specific Method
SEASON::SEASON() {
  m_Event = new TSEASONData();
  m_SEASONScorer = 0;

  m_DistDSSD1 = 2 * mm;
  m_DistTunnel = 2 * mm;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0.5, 1));
  m_VisWheel = new G4VisAttributes(G4Colour(1, 1, 0, 1));
  m_VisFoil = new G4VisAttributes(G4Colour(1, 0, 0, 1));
  m_VisGe = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5));
  m_VisWindow = new G4VisAttributes(G4VisAttributes::GetInvisible());
  m_VisGrid = new G4VisAttributes(G4VisAttributes::GetInvisible());
  m_VisChamber = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5));
}

SEASON::~SEASON() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SEASON::AddDetector(G4ThreeVector X1_Y1, G4ThreeVector X1_YMax, G4ThreeVector XMax_Y1, G4ThreeVector XMax_YMax) {
  m_X1_Y1.push_back(X1_Y1);
  m_X1_YMax.push_back(X1_YMax);
  m_XMax_Y1.push_back(XMax_Y1);
  m_XMax_YMax.push_back(XMax_YMax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SEASON::AddDetector(double DistDSSD1, double DistTunnel, bool Deloc, bool EnableWheel, bool RotateWheel,
                         bool EnableChamber) {

  m_EnableWheel = EnableWheel;
  m_RotateWheel = RotateWheel;
  m_EnableChamber = EnableChamber;

  m_DistDSSD1 = DistDSSD1;
  m_DistTunnel = DistTunnel;

  double pos = SEASON_NS::DetSize * 0.5;
  double interX = 0.2 * mm;
  double interY = 0.2 * mm;

  // Add DSSD1
  AddDetector(G4ThreeVector(-pos, -pos, DistDSSD1), G4ThreeVector(-pos, pos, DistDSSD1),
              G4ThreeVector(pos, -pos, DistDSSD1), G4ThreeVector(pos, pos, DistDSSD1));

  // Add Tunnel detectors
  AddDetector(G4ThreeVector(pos + interX, -pos, -DistTunnel), G4ThreeVector(pos + interX, pos, -DistTunnel),
              G4ThreeVector(pos + interX, -pos, -DistTunnel - SEASON_NS::DetSize),
              G4ThreeVector(pos + interX, pos, -DistTunnel - SEASON_NS::DetSize));
  AddDetector(G4ThreeVector(pos, pos + interY, -DistTunnel), G4ThreeVector(-pos, pos + interY, -DistTunnel),
              G4ThreeVector(pos, pos + interY, -DistTunnel - SEASON_NS::DetSize),
              G4ThreeVector(-pos, pos + interY, -DistTunnel - SEASON_NS::DetSize));
  AddDetector(G4ThreeVector(-pos - interX, pos, -DistTunnel), G4ThreeVector(-pos - interX, -pos, -DistTunnel),
              G4ThreeVector(-pos - interX, pos, -DistTunnel - SEASON_NS::DetSize),
              G4ThreeVector(-pos, -pos - interX, -DistTunnel - SEASON_NS::DetSize));
  AddDetector(G4ThreeVector(-pos, -pos - interY, -DistTunnel), G4ThreeVector(pos, -pos - interY, -DistTunnel),
              G4ThreeVector(-pos, -pos - interY, -DistTunnel - SEASON_NS::DetSize),
              G4ThreeVector(pos, -pos - interY, -DistTunnel - SEASON_NS::DetSize));

  // Add Delocalized station detectors if activated
  if (Deloc) {
    G4double angle = -4. * 2. * TMath::Pi() / SEASON_NS::FoilNbr - TMath::Pi() * 0.5;
    G4ThreeVector* DelocFoilPos = new G4ThreeVector(TMath::Cos(angle) * SEASON_NS::FoilWheelCenterDist,
                                                    (TMath::Sin(angle) + 1) * SEASON_NS::FoilWheelCenterDist, 0 * cm);

    double DelocFoilPosX = DelocFoilPos->getX();
    double DelocFoilPosY = DelocFoilPos->getY();

    AddDetector(G4ThreeVector(DelocFoilPosX - pos, DelocFoilPosY - pos, DistDSSD1),
                G4ThreeVector(DelocFoilPosX - pos, DelocFoilPosY + pos, DistDSSD1),
                G4ThreeVector(DelocFoilPosX + pos, DelocFoilPosY - pos, DistDSSD1),
                G4ThreeVector(DelocFoilPosX + pos, DelocFoilPosY + pos, DistDSSD1));
    AddDetector(G4ThreeVector(DelocFoilPosX - pos, DelocFoilPosY - pos, -DistTunnel),
                G4ThreeVector(DelocFoilPosX - pos, DelocFoilPosY + pos, -DistTunnel),
                G4ThreeVector(DelocFoilPosX + pos, DelocFoilPosY - pos, -DistTunnel),
                G4ThreeVector(DelocFoilPosX + pos, DelocFoilPosY + pos, -DistTunnel));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SEASON::BuildSquareDetector() {

  G4Box* box = new G4Box("SEASON_Box", SEASON_NS::DetSize * 0.5, SEASON_NS::DetSize * 0.5, SEASON_NS::Thickness * 0.5);
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4LogicalVolume* m_SquareDetector = new G4LogicalVolume(box, DetectorMaterial, "logic_SEASON_Box", 0, 0, 0);
  m_SquareDetector->SetVisAttributes(m_VisSquare);
  m_SquareDetector->SetSensitiveDetector(m_SEASONScorer);
  return m_SquareDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SEASON::BuildWindow() {

  G4Box* solidWindow =
      new G4Box("WindowBox", SEASON_NS::DetSize * 0.5, SEASON_NS::DetSize * 0.5, SEASON_NS::DeadLayerThickness * 0.5);
  G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow, SEASON_NS::DeadLayerMaterial, "logicWindow", 0, 0, 0);
  // logicWindow->SetVisAttributes(new G4VisAttributes(G4Colour(0, 1, 1, 0.5)));
  logicWindow->SetVisAttributes(m_VisWindow);
  return logicWindow;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SEASON::BuildChamber() {
  G4Box* plainChamber = new G4Box("PlainChamber", 300 * mm, 372.75 * mm, 4 * mm);

  //  make holes
  G4Tubs* hole =
      new G4Tubs("Hole", 0, 141.5 * mm + 1 * micrometer, 4 * mm + 100 * micrometer, 0 * degree, 360 * degree);

  G4SubtractionSolid* DrilledOnceChamber = new G4SubtractionSolid("DrilledOnceChamber", plainChamber, hole, NULL,
                                                                  G4ThreeVector(70.9 * mm, -159.05 * mm, 1.5 * mm));
  G4SubtractionSolid* DrilledChamber = new G4SubtractionSolid("DrilledChamber", DrilledOnceChamber, hole, NULL,
                                                              G4ThreeVector(-50 * mm, 105.75 * mm, 1.5 * mm));

  G4Material* ChamberMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4LogicalVolume* m_Chamber = new G4LogicalVolume(DrilledChamber, ChamberMaterial, "logic_SEASON_Chamber");
  m_Chamber->SetVisAttributes(m_VisChamber);

  return m_Chamber;
}

G4LogicalVolume* SEASON::BuildGrid() {
  G4Box* plainStrip = new G4Box("plainStrip", SEASON_NS::StripPitch * 0.5, SEASON_NS::StripLength * 0.5,
                                SEASON_NS::GridThickness * 0.5);
  G4Box* hole =
      new G4Box("hole", SEASON_NS::StripPitch * 0.5 - SEASON_NS::GridWidth,
                SEASON_NS::StripLength * 0.5 - SEASON_NS::GridWidth, SEASON_NS::GridThickness * 0.5 + 10 * um);
  G4SubtractionSolid* solidStrip = new G4SubtractionSolid("SEASON_Strip", plainStrip, hole, 0, G4ThreeVector(0, 0, 0));

  G4LogicalVolume* logicStrip = new G4LogicalVolume(solidStrip, SEASON_NS::Grid_Material, "logicStrip");
  // logicStrip->SetVisAttributes(new G4VisAttributes(G4Colour(1, 0, 0, 1)));
  logicStrip->SetVisAttributes(m_VisGrid);
  return logicStrip;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  Wheel with the mounting implantation foils
G4LogicalVolume* SEASON::BuildWheel() {
  G4double WheelRadius = SEASON_NS::WheelRadius;
  G4double WheelThickness = SEASON_NS::WheelThickness;
  G4double FinLenght = SEASON_NS::FinLenght;
  G4double FinWidth = SEASON_NS::FinWidth;
  G4double FinThickness = SEASON_NS::FinThickness;
  G4double FoilNbr = SEASON_NS::FoilNbr;
  // Make wheel
  G4Tubs* plainWheel = new G4Tubs("PlainWheel", 0, WheelRadius, WheelThickness * 0.5, 0 * degree, 360 * degree);

  // Make fins
  G4Box* plainFin = new G4Box("PlainFin", FinWidth * 0.5, FinLenght * 0.5 + 0.5 * cm, FinThickness * 0.5);

  //  make holes
  G4Tubs* hole = new G4Tubs("Hole", 0, m_FoilRadius + 1 * micrometer, FinThickness * 0.5 + 100 * micrometer, 0 * degree,
                            360 * degree);

  G4UnionSolid* wheel =
      new G4UnionSolid("SEASON_Wheel", plainWheel, plainFin, 0,
                       G4ThreeVector(0, WheelRadius + FinLenght * 0.5 - 0.5 * cm, -FinThickness * 0.5));
  G4SubtractionSolid* wheel2 = new G4SubtractionSolid(
      "SEASON_Wheel", wheel, hole, NULL, G4ThreeVector(0, SEASON_NS::FoilWheelCenterDist, -FinThickness * 0.5));

  for (int i = 1; i < FoilNbr; i++) {
    G4RotationMatrix zRot;
    zRot.rotateZ(i * 2 * TMath::Pi() / FoilNbr * rad);
    G4UnionSolid* wheelbis = new G4UnionSolid(
        "SEASON_Wheel", wheel2, plainFin, &zRot,
        G4ThreeVector(TMath::Sin(i * 2 * TMath::Pi() / FoilNbr) * (WheelRadius + FinLenght * 0.5 - 0.5 * cm),
                      TMath::Cos(i * 2 * TMath::Pi() / FoilNbr) * (WheelRadius + FinLenght * 0.5 - 0.5 * cm),
                      -FinThickness * 0.5));
    wheel = wheelbis;
    G4SubtractionSolid* wheelbis2 = new G4SubtractionSolid(
        "SEASON_Wheel", wheel, hole, NULL,
        G4ThreeVector(TMath::Sin(i * 2 * TMath::Pi() / FoilNbr) * SEASON_NS::FoilWheelCenterDist,
                      TMath::Cos(i * 2 * TMath::Pi() / FoilNbr) * SEASON_NS::FoilWheelCenterDist, -FinThickness * 0.5));
    wheel2 = wheelbis2;
  }

  G4Material* WheelMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4LogicalVolume* m_Wheel = new G4LogicalVolume(wheel2, WheelMaterial, "logic_SEASON_Wheel");
  m_Wheel->SetVisAttributes(m_VisWheel);

  return m_Wheel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  Wheel with the mounting implantation foils
G4LogicalVolume* SEASON::BuildRotatedWheel() {
  G4double WheelRadius = SEASON_NS::WheelRadius;
  G4double WheelThickness = SEASON_NS::WheelThickness;
  G4double FinLenght = SEASON_NS::FinLenght;
  G4double FinWidth = SEASON_NS::FinWidth;
  G4double FinThickness = SEASON_NS::FinThickness;
  G4double FoilNbr = SEASON_NS::FoilNbr;
  // Make wheel
  G4Tubs* plainWheel = new G4Tubs("PlainWheel", 0, WheelRadius, WheelThickness * 0.5, 0 * degree, 360 * degree);

  // Make fins
  G4Box* plainFin = new G4Box("PlainFin", FinWidth * 0.5, FinLenght * 0.5 + 0.5 * cm, FinThickness * 0.5);

  //  make holes
  G4Tubs* hole = new G4Tubs("Hole", 0, m_FoilRadius + 1 * micrometer, FinThickness * 0.5 + 100 * micrometer, 0 * degree,
                            360 * degree);

  G4RotationMatrix zRot0;
  zRot0.rotateZ(TMath::Pi() / FoilNbr * rad);

  G4UnionSolid* wheel = new G4UnionSolid(
      "SEASON_Wheel", plainWheel, plainFin, &zRot0,
      G4ThreeVector(TMath::Sin(TMath::Pi() / FoilNbr) * (WheelRadius + FinLenght * 0.5 - 0.5 * cm),
                    TMath::Cos(TMath::Pi() / FoilNbr) * (WheelRadius + FinLenght * 0.5 - 0.5 * cm), -FinThickness));

  G4SubtractionSolid* wheel2 = new G4SubtractionSolid(
      "SEASON_Wheel", wheel, hole, NULL,
      G4ThreeVector(TMath::Sin(TMath::Pi() / FoilNbr) * SEASON_NS::FoilWheelCenterDist,
                    TMath::Cos(TMath::Pi() / FoilNbr) * SEASON_NS::FoilWheelCenterDist, -FinThickness));

  for (int i = 1; i < FoilNbr; i++) {
    G4RotationMatrix zRot;
    zRot.rotateZ((i + 1 * 0.5) * 2 * TMath::Pi() / FoilNbr * rad);
    G4UnionSolid* wheelbis = new G4UnionSolid(
        "SEASON_Wheel", wheel2, plainFin, &zRot,
        G4ThreeVector(
            TMath::Sin((i + 1 * 0.5) * 2 * TMath::Pi() / FoilNbr) * (WheelRadius + FinLenght * 0.5 - 0.5 * cm),
            TMath::Cos((i + 1 * 0.5) * 2 * TMath::Pi() / FoilNbr) * (WheelRadius + FinLenght * 0.5 - 0.5 * cm),
            -FinThickness));
    wheel = wheelbis;
    G4SubtractionSolid* wheelbis2 = new G4SubtractionSolid(
        "SEASON_Wheel", wheel, hole, NULL,
        G4ThreeVector(TMath::Sin((i + 1 * 0.5) * 2 * TMath::Pi() / FoilNbr) * SEASON_NS::FoilWheelCenterDist,
                      TMath::Cos((i + 1 * 0.5) * 2 * TMath::Pi() / FoilNbr) * SEASON_NS::FoilWheelCenterDist,
                      -FinThickness));
    wheel2 = wheelbis2;
  }

  G4Material* WheelMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4LogicalVolume* m_Wheel = new G4LogicalVolume(wheel2, WheelMaterial, "logic_SEASON_Wheel");
  m_Wheel->SetVisAttributes(m_VisWheel);

  return m_Wheel;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//  The implantation foil (Carbon)
G4LogicalVolume* SEASON::BuildFoil() {
  G4Tubs* foil = new G4Tubs("ImplantationFoil", 0, m_FoilRadius, m_FoilThickness * 0.5, 0 * degree, 360 * degree);
  G4LogicalVolume* m_Foil = new G4LogicalVolume(foil, SEASON_NS::FoilMaterial, "logic_SEASON_Foil");

  m_Foil->SetVisAttributes(m_VisFoil); // red G4Colour(0.5, 0.5, 0.5, 1)); //grey

  return m_Foil;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void SEASON::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SEASON");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> cart = {"X1_Y1", "X1_YMax", "XMax_Y1", "XMax_YMax"};
  vector<string> cartfirst = {"X1_Y1", "X1_YMax", "XMax_Y1", "XMax_YMax", "EnergyThreshold"};
  vector<string> full = {"DSSD1Dist",  "TunnelDist",      "UseDelocStation", "UseWheel",
                         "UseChamber", "EnergyThreshold", "FoilRadius",      "FoilThickness"};
  double FoilThickness;
  for (unsigned int i = 0; i < blocks.size(); i++) {
    // cout << "!!! !!! !!! OK1 " << endl;
    if (blocks[i]->HasTokenList(cartfirst)) {
      // cout << "!!! !!! !!! Cart ok" << endl;
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SEASON " << i + 1 << endl;

      G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("X1_YMax", "mm"));
      G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("XMax_Y1", "mm"));
      G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("XMax_YMax", "mm"));
      m_EnergyThreshold = blocks[i]->GetDouble("EnergyThreshold", "keV");
      m_FoilThickness = (FoilThickness / 0.226) * nm;
      AddDetector(A, B, C, D);
    }
    else if (blocks[i]->HasTokenList(cart)) {
      // cout << "!!! !!! !!! Cart ok" << endl;
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SEASON " << i + 1 << endl;

      G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("X1_Y1", "mm"));
      G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("X1_YMax", "mm"));
      G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("XMax_Y1", "mm"));
      G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("XMax_YMax", "mm"));
      AddDetector(A, B, C, D);
    }
    else if (blocks[i]->HasTokenList(full)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SEASON " << i + 1 << endl;

      double DistDSSD1 = blocks[i]->GetDouble("DSSD1Dist", "mm");
      double DistTunnel = blocks[i]->GetDouble("TunnelDist", "mm");

      m_EnergyThreshold = blocks[i]->GetDouble("EnergyThreshold", "keV");
      m_FoilRadius = blocks[i]->GetDouble("FoilRadius", "mm");
      FoilThickness = blocks[i]->GetDouble("FoilThickness", "void");
      m_FoilThickness = (FoilThickness / 0.226) * nm;

      string Deloc = blocks[i]->GetString("UseDelocStation");
      bool EnableDeloc;
      if (Deloc == "True")
        EnableDeloc = true;
      else if (Deloc == "False")
        EnableDeloc = false;
      else {
        cout << "ERROR: UseDelocStation must be either True or False, please check input file formatting" << endl;
        exit(1);
      }

      string Wheel = blocks[i]->GetString("UseWheel");
      bool EnableWheel = false;
      bool RotateWheel = false;
      if (Wheel == "Rotated")
        RotateWheel = true;
      else if (Wheel == "True")
        EnableWheel = true;
      else if (Wheel == "False")
        EnableWheel = false;
      else {
        cout << "ERROR: UseWheel must be either True, False or Rotated, please check input file formatting" << endl;
        exit(1);
      }

      string Chamber = blocks[i]->GetString("UseChamber");
      bool EnableChamber;
      if (Chamber == "True")
        EnableChamber = true;
      else if (Chamber == "False")
        EnableChamber = false;
      else {
        cout << "ERROR: UseChamber must be either True or False, please check input file formatting" << endl;
        exit(1);
      }

      AddDetector(DistDSSD1, DistTunnel, EnableDeloc, EnableWheel, RotateWheel, EnableChamber);
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  cout << " Foil Thickness : " << FoilThickness << " -> " << m_FoilThickness / nm
       << " nm with a carbon density rho = 2.26 g/cm3" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void SEASON::ConstructDetector(G4LogicalVolume* world) {

  G4bool checkOverlaps = true;

  //
  // Construct DSSDs

  for (unsigned short i = 0; i < m_X1_Y1.size(); i++) {

    G4RotationMatrix* Rot = NULL;
    G4ThreeVector Det_pos = G4ThreeVector(0, 0, 0);
    G4ThreeVector u = G4ThreeVector(0, 0, 0);
    G4ThreeVector v = G4ThreeVector(0, 0, 0);
    G4ThreeVector w = G4ThreeVector(0, 0, 0);
    G4ThreeVector Center = G4ThreeVector(0, 0, 0);

    G4ThreeVector Front_Window_pos = G4ThreeVector(0, 0, 0);
    G4ThreeVector Grid_pos = G4ThreeVector(0, 0, 0);
    // (u,v,w) unitary vector associated to telescope referencial
    // (u,v) // to silicon plan
    // w perpendicular to (u,v) plan and pointing CsI
    u = m_XMax_Y1[i] - m_X1_Y1[i];
    u = u.unit();
    v = m_X1_YMax[i] - m_X1_Y1[i];
    v = v.unit();
    w = u.cross(v);
    w = w.unit();

    Center = (m_X1_Y1[i] + m_X1_YMax[i] + m_XMax_Y1[i] + m_XMax_YMax[i]) / 4;

    // Passage Matrix from Lab Referential to Telescope Referential
    Rot = new G4RotationMatrix(u, v, w);
    if (w == G4ThreeVector(0, 0, -1) && Center.getZ() > 0)
      w = -w;
    Det_pos = w * SEASON_NS::Thickness * 0.5 + Center;

    Front_Window_pos = -w * SEASON_NS::DeadLayerThickness * 0.5 + Center;
    Grid_pos = -w * (SEASON_NS::DeadLayerThickness + SEASON_NS::GridThickness * 0.5) + Center;

    if (SEASON_NS::DeadLayerThickness > 0) {
      new G4PVPlacement(G4Transform3D(*Rot, Front_Window_pos), BuildWindow(), "SEASON", world, false,
                        m_X1_Y1.size() + i + 1, FALSE);
    }
    else if (i == 0)
      cout << "No entrance window " << endl;

    if (SEASON_NS::GridThickness > 0 and SEASON_NS::GridWidth > 0) {
      for (int strip = 0; strip < SEASON_NS::NumberOfStripsX; strip++) {
        G4ThreeVector Strip_pos = u * ((-strip + (SEASON_NS::NumberOfStripsX - 1) / 2) * SEASON_NS::StripWidth);
        new G4PVPlacement(G4Transform3D(*Rot, Grid_pos + Strip_pos), BuildGrid(), "SEASON", world, false,
                          2 * m_X1_Y1.size() + i + 1, FALSE);
      }
    }
    else if (i == 0)
      cout << "No aluminium grid" << endl;

    new G4PVPlacement(G4Transform3D(*Rot, Det_pos), BuildSquareDetector(), "SEASON", world, false, i + 1,
                      checkOverlaps);
  }

  if (m_RotateWheel) {
    //
    //  Build Wheel
    G4RotationMatrix zRot;
    zRot.rotateZ(TMath::Pi() * rad);

    new G4PVPlacement(G4Transform3D(zRot, G4ThreeVector(0, SEASON_NS::FoilWheelCenterDist, 0)), BuildRotatedWheel(),
                      "Wheel", world, false, 0, checkOverlaps);

    // Build foils in the holes
    for (int i = 0; i < SEASON_NS::FoilNbr; i++) {
      new G4PVPlacement(
          G4Transform3D(G4RotationMatrix(0, 0, 0),
                        G4ThreeVector(TMath::Sin((i + 1 * 0.5) * 2 * TMath::Pi() / SEASON_NS::FoilNbr + TMath::Pi()) *
                                          SEASON_NS::FoilWheelCenterDist,
                                      TMath::Cos((i + 1 * 0.5) * 2 * TMath::Pi() / SEASON_NS::FoilNbr + TMath::Pi()) *
                                              SEASON_NS::FoilWheelCenterDist +
                                          SEASON_NS::FoilWheelCenterDist,
                                      m_FoilThickness * 0.5)),
          BuildFoil(), "Foil", world, false, 0, checkOverlaps);
    }
  }
  else if (m_EnableWheel) {
    //
    //  Build Wheel
    G4RotationMatrix zRot;
    zRot.rotateZ(TMath::Pi() * rad);

    new G4PVPlacement(G4Transform3D(zRot, G4ThreeVector(0, SEASON_NS::FoilWheelCenterDist, 0)), BuildWheel(), "Wheel",
                      world, false, 0, checkOverlaps);

    // Build foils in the holes
    for (int i = 0; i < SEASON_NS::FoilNbr; i++) {
      new G4PVPlacement(G4Transform3D(G4RotationMatrix(0, 0, 0),
                                      G4ThreeVector(TMath::Sin(i * 2 * TMath::Pi() / SEASON_NS::FoilNbr + TMath::Pi()) *
                                                        SEASON_NS::FoilWheelCenterDist,
                                                    TMath::Cos(i * 2 * TMath::Pi() / SEASON_NS::FoilNbr + TMath::Pi()) *
                                                            SEASON_NS::FoilWheelCenterDist +
                                                        SEASON_NS::FoilWheelCenterDist,
                                                    m_FoilThickness * 0.5)),
                        BuildFoil(), "Foil", world, false, 0, checkOverlaps);
    }
  }

  if (m_EnableChamber) {
    //
    //  Build Wheel
    G4RotationMatrix zRot;
    new G4PVPlacement(
        G4Transform3D(zRot, G4ThreeVector(-70.9 * mm, 159.05 * mm, m_DistDSSD1 + SEASON_NS::DistDSSD1_AlWindow)),
        BuildChamber(), "Chamber", world, false, 0, checkOverlaps);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SEASON::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("SEASON")) {
    pTree->Branch("SEASON", "TSEASONData", &m_Event);
  }
  pTree->SetBranchAddress("SEASON", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void SEASON::ReadSensitive(const G4Event* event) {
  m_Event->Clear();

  ///////////
  // Strip scorer
  DSSDScorers::PS_Images* SiScorer = (DSSDScorers::PS_Images*)m_SEASONScorer->GetPrimitive(0);
  map<unsigned int, pair<double, double>> mapFront;
  map<unsigned int, unsigned short int> mapRFront;
  map<unsigned int, pair<double, double>>::iterator it;

  // X-side reading
  for (unsigned int i = 0; i < SiScorer->GetFrontMult(); i++) {
    double EnergyX = SiScorer->GetEnergyFront(i);
    unsigned int DetNbrX = SiScorer->GetDetectorFront(i);
    double TimeX = SiScorer->GetTimeFront(i);

    unsigned int a, r, g, b;
    SiScorer->GetARGBFront(i, a, r, g, b);

    if (r == 0) {
      mapFront[b + DetNbrX * 1e6].first += EnergyX;
      mapFront[b + DetNbrX * 1e6].second += TimeX;
      mapRFront[b + DetNbrX * 1e6] += r + i;
    }
    else {
      double rand = G4UniformRand();
      if (rand > 0.5) {
        double energy1 = EnergyX * rand;
        double energy2 = EnergyX * (rand - 1);
        mapFront[b + DetNbrX * 1e6].first += energy1;
        mapFront[b + DetNbrX * 1e6].second = TimeX;
        mapRFront[b + DetNbrX * 1e6] += r + i;
        mapFront[g + DetNbrX * 1e6].first += energy2;
        mapFront[g + DetNbrX * 1e6].second = TimeX;
        mapRFront[g + DetNbrX * 1e6] += r + i;
      }
      else {
        double energy1 = -EnergyX * rand;
        double energy2 = EnergyX * (1 - rand);
        mapFront[b + DetNbrX * 1e6].first += energy1;
        mapFront[b + DetNbrX * 1e6].second = TimeX;
        mapRFront[b + DetNbrX * 1e6] += r + i;
        mapFront[g + DetNbrX * 1e6].first += energy2;
        mapFront[g + DetNbrX * 1e6].second = TimeX;
        mapRFront[g + DetNbrX * 1e6] += r + i;
      }
    }
  }

  for (it = mapFront.begin(); it != mapFront.end(); it++) {
    double EnergyX = it->second.first;
    if (abs(EnergyX) > m_EnergyThreshold) {
      if (EnergyX > 2.)
        EnergyX = RandGauss::shoot(EnergyX, SEASON_NS::ResoAlpha);
      else
        EnergyX = RandGauss::shoot(EnergyX, SEASON_NS::ResoElectron);
      double TimeX = RandGauss::shoot(it->second.second, SEASON_NS::ResoTime);
      unsigned int strip = it->first - 1000000 * (it->first / 1000000);
      unsigned int det = it->first / 1000000;
      m_Event->SetXEnergy(det, strip, EnergyX);
      m_Event->SetXTime(det, strip, TimeX);
      m_Event->SetXParticleID(mapRFront[it->first]);
    }
  }

  map<unsigned int, pair<double, double>> mapBack;
  map<unsigned int, unsigned short int> mapRBack;
  // Y-side reading
  for (unsigned int i = 0; i < SiScorer->GetBackMult(); i++) {
    double EnergyY = SiScorer->GetEnergyBack(i);
    unsigned int DetNbrY = SiScorer->GetDetectorBack(i);
    double TimeY = SiScorer->GetTimeBack(i);

    unsigned int a, r, g, b;
    SiScorer->GetARGBBack(i, a, r, g, b);

    if (r == 0) {
      mapBack[b + DetNbrY * 1e6].first += EnergyY;
      mapBack[b + DetNbrY * 1e6].second += TimeY;
      mapRBack[b + DetNbrY * 1e6] += r + i;
    }
    else {
      double rand = G4UniformRand();
      double energy1 = rand * EnergyY;
      double energy2 = (1 - rand) * EnergyY;
      mapBack[b + DetNbrY * 1e6].first += energy1;
      mapBack[b + DetNbrY * 1e6].second += TimeY;
      mapRBack[b + DetNbrY * 1e6] += r + i;
      mapBack[g + DetNbrY * 1e6].first += energy2;
      mapBack[g + DetNbrY * 1e6].second += TimeY;
      mapRBack[g + DetNbrY * 1e6] += r + i;
    }
  }

  for (it = mapBack.begin(); it != mapBack.end(); it++) {
    double EnergyY = it->second.first;
    if (EnergyY > m_EnergyThreshold) {
      if (EnergyY > 2.)
        EnergyY = RandGauss::shoot(EnergyY, SEASON_NS::ResoAlpha);
      else
        EnergyY = RandGauss::shoot(EnergyY, SEASON_NS::ResoElectron);
      double TimeY = RandGauss::shoot(it->second.second, SEASON_NS::ResoTime);
      unsigned int strip = it->first - 1000000 * (it->first / 1000000);
      unsigned int det = it->first / 1000000;
      m_Event->SetYEnergy(det, strip, EnergyY);
      m_Event->SetYTime(det, strip, TimeY);
      m_Event->SetYParticleID(mapRBack[it->first]);
    }
  }
  SiScorer->clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void SEASON::InitializeScorers() {
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;

  if (already_exist)
    return;

  m_SEASONScorer = CheckScorer("SEASONScorer", already_exist);

  cout << "Size of scorers :\n";
  cout << SEASON_NS::DetSize << " " << SEASON_NS::DetSize;
  cout << " " << SEASON_NS::NumberOfStripsX << " " << SEASON_NS::NumberOfStripsY << "\n";

  string nptool = getenv("NPTOOL");
  G4VPrimitiveScorer* Scorer = new DSSDScorers::PS_Images(
      "Scorer", nptool + "/NPLib/Detectors/SEASON/ressources/maskFront.png",
      nptool + "/NPLib/Detectors/SEASON/ressources/maskBack.png", 0.005, 0.005, 0, 0, 0xffff0000, 0);

  m_SEASONScorer->RegisterPrimitive(Scorer);

  // Interaction scorer
  m_SEASONScorer->RegisterPrimitive(
      new InteractionScorers::PS_Interactions("InteractionScorerSEASON", ms_InterCoord, 0));
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SEASONScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SEASON::Construct() { return (NPS::VDetector*)new SEASON(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_SEASON {
 public:
  proxy_nps_SEASON() {
    NPS::DetectorFactory::getInstance()->AddToken("SEASON", "SEASON");
    NPS::DetectorFactory::getInstance()->AddDetector("SEASON", SEASON::Construct);
  }
};

proxy_nps_SEASON p_nps_SEASON;
}
