/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Tina simulation                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>

// Geant4
#include "G4Box.hh"
#include "G4Colour.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4NistManager.hh"
#include "G4ProductionCuts.hh"
#include "G4SDManager.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"

// NPTool
#include "BeamReaction.hh"
#include "Decay.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"
#include "STARK.hh"
#include "STARKScorers.hh"

// CLHEP
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace STARKNS;

namespace STARKNS {
  G4double EnergyThreshold = 0.1 * MeV;
  ////////////////////////////////////////////////////////////
  // Resolution
  ////////////////////////////////////////////////////////////
  G4double X6_TRes = 0.213;
  G4double X6_ERes = 0.5; // [%] in sigma
  G4double BB10_TRes = 0.213;
  G4double BB10_ERes = 0.5; // [%] in sigma
  G4double QQQ5_TRes = 0.213;
  G4double QQQ5_ERes = 0.5; // [%] in sigma
  G4double CsI_TRes = 0.213;
  G4double CsI_ERes = 0.5; // [%] in sigma
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Geometry
  ////////////////////////////////////////////////////////////
  //
  // X6
  ////////////////////////////////////////////////////////////
  const G4double X6_PCBX = 45.20 * mm;
  const G4double X6_PCBY = 93.10 * mm;
  const G4double X6_PCBZ = 2.40 * mm;
  const G4double X6_PCBSub1X = 43.60 * mm;
  const G4double X6_PCBSub1Y = 78.30 * mm;
  const G4double X6_PCBSub1Z = 1.20 * mm;
  const G4double X6_PCBSub1XOffset = 0.0 * mm;
  const G4double X6_PCBSub1YOffset = 6.2 * mm;
  const G4double X6_PCBSub1ZOffset = 0.6 * mm;
  const G4double X6_PCBSub2X = 42.20 * mm;
  const G4double X6_PCBSub2Y = 76.90 * mm;
  const G4double X6_PCBSub2Z = 1.20 * mm;
  const G4double X6_PCBSub2XOffset = 0.0 * mm;
  const G4double X6_PCBSub2YOffset = 6.2 * mm;
  const G4double X6_PCBSub2ZOffset = -0.6 * mm;

  const G4double X6_SiX = 43.30 * mm;
  const G4double X6_SiY = 78.00 * mm;
  const G4double X6_SiZ = 1.00 * mm;
  const G4double X6_SiXOffset = 0.0 * mm;
  const G4double X6_SiYOffset = 6.2 * mm;
  const G4double X6_SiZOffset = 0.5 * mm;
  const G4double X6_SiActiveX = 40.30 * mm;
  const G4double X6_SiActiveY = 75.00 * mm;
  const G4double X6_SiActiveZ = 1.00 * mm; // 1000 um

  const G4int X6_NFrontStrips = 8;
  const G4int X6_NBackStrips = 4;
  ////////////////////////////////////////////////////////////
  //
  // BB10
  ////////////////////////////////////////////////////////////
  const G4double BB10_PCBX = 45.20 * mm;
  const G4double BB10_PCBY = 93.10 * mm;
  const G4double BB10_PCBZ = 2.40 * mm;
  const G4double BB10_PCBSub1X = 43.60 * mm;
  const G4double BB10_PCBSub1Y = 78.30 * mm;
  const G4double BB10_PCBSub1Z = 1.20 * mm;
  const G4double BB10_PCBSub1XOffset = 0.0 * mm;
  const G4double BB10_PCBSub1YOffset = 6.5 * mm; // only different to X6 PCB
  const G4double BB10_PCBSub1ZOffset = 0.6 * mm;
  const G4double BB10_PCBSub2X = 42.20 * mm;
  const G4double BB10_PCBSub2Y = 76.90 * mm;
  const G4double BB10_PCBSub2Z = 1.20 * mm;
  const G4double BB10_PCBSub2XOffset = 0.0 * mm;
  const G4double BB10_PCBSub2YOffset = 6.5 * mm; // only different to X6 PCB
  const G4double BB10_PCBSub2ZOffset = -0.6 * mm;

  const G4double BB10_SiX = 43.30 * mm;
  const G4double BB10_SiY = 78.00 * mm;
  const G4double BB10_SiZ = 0.14 * mm;
  const G4double BB10_SiXOffset = 0.0 * mm;
  const G4double BB10_SiYOffset = 6.5 * mm;
  const G4double BB10_SiZOffset = 0.07 * mm;
  const G4double BB10_SiActiveX = 39.45 * mm;
  const G4double BB10_SiActiveY = 74.15 * mm;
  const G4double BB10_SiActiveZ = 0.14 * mm; // 140 um

  const G4int BB10_NFrontStrips = 8;
  const G4int BB10_NBackStrips = 1;
  ////////////////////////////////////////////////////////////
  //
  // QQQ5
  ////////////////////////////////////////////////////////////
  const G4double QQQ5_PCBOutR = 86 * mm;
  const G4double QQQ5_PCBInR = 15 * mm;
  const G4double QQQ5_PCBPhi0 = 0 * deg;  // Starting point
  const G4double QQQ5_PCBPhi1 = 90 * deg; // ANGLE
  const G4double QQQ5_PCBT = 3.4 * mm;
  const G4double QQQ5_PCBCutX = 3.4 * 2 * mm; // 3.4 mm gap from the arc center?
  const G4double QQQ5_PCBCutY = 86 * mm;
  const G4double QQQ5_PCBCutZ = 4 * mm;
  const G4double QQQ5_PCBCutXOffset = 0 * mm;
  const G4double QQQ5_PCBCutYOffset = 43 * mm;
  const G4double QQQ5_PCBCutZOffset = 0 * mm;

  ////////////////////////////////////////////////////////////
  //
  // QQQ Wafer
  ////////////////////////////////////////////////////////////
  const G4double QQQ5_SiOutR = 84.0 * mm;
  const G4double QQQ5_SiInR = 23.2 * mm;
  const G4double QQQ5_SiT = 1 * mm;
  const G4double QQQ5_SiPhi0 = 0 * deg;
  const G4double QQQ5_SiPhi1 = 90 * deg;
  const G4double QQQ5_SiActiveOutR = 81.95 * mm;
  const G4double QQQ5_SiActiveInR = 25.25 * mm;
  const G4double QQQ5_SiCut1X = (3.4 + 0.68) * 2 * mm;
  const G4double QQQ5_SiCut1Y = QQQ5_SiOutR;
  const G4double QQQ5_SiCut1Z = 2 * mm;
  const G4double QQQ5_SiCut1XOffset = 0 * mm;
  const G4double QQQ5_SiCut1YOffset = QQQ5_SiCut1Y / 2;
  const G4double QQQ5_SiCut1ZOffset = 0 * mm;
  const G4double QQQ5_SiCut2X = QQQ5_SiOutR;
  const G4double QQQ5_SiCut2Y = 0.92 * 2 * mm;
  const G4double QQQ5_SiCut2Z = 2 * mm;
  const G4double QQQ5_SiCut2XOffset = QQQ5_SiCut2X / 2;
  const G4double QQQ5_SiCut2YOffset = 0 * mm;
  const G4double QQQ5_SiCut2ZOffset = 0 * mm;

  const G4int QQQ5_NRStrip = 32;
  const G4int QQQ5_NAStrip = 4;

  ////////////////////////////////////////////////////////////
  //
  // Connector
  ////////////////////////////////////////////////////////////
  const G4double Conn_X = 40.0 * mm;
  const G4double Conn_Y = 5.0 * mm;
  const G4double Conn_Z = 5.0 * mm;

  ////////////////////////////////////////////////////////////
  //
  // CsI
  ////////////////////////////////////////////////////////////
  const G4double CsI_X = 40.80 * mm;
  const G4double CsI_Y = 40.80 * mm;
  const G4double CsI_Z = 25.50 * mm;
  const G4double CsI_X6_XOffset = 0;
  const G4double CsI_X6_YOffset1 = 0.5 * X6_SiY - 0.5 * CsI_Y;
  const G4double CsI_X6_YOffset2 = CsI_X6_YOffset1 - CsI_Y;
  const G4double CsI_X6_ZOffset = 5.0 * mm;
  const G4double CsI_BB10_XOffset = 0;
  const G4double CsI_BB10_YOffset1 = 0.5 * BB10_SiY - 0.5 * CsI_Y;
  const G4double CsI_BB10_YOffset2 = CsI_BB10_YOffset1 - CsI_Y;
  const G4double CsI_BB10_ZOffset = 5.0 * mm;

  ////////////////////////////////////////////////////////////
  //
  // ANASEN CsI for QQQ3
  ////////////////////////////////////////////////////////////
  const G4double CsI_QQQ5_XOffset = 0;
  const G4double CsI_QQQ5_YOffset = 7 * mm;
  const G4double CsI_QQQ5_ZOffset = 5.0 * mm;
  const G4double m_ANASENQQQ3CsIhypotenuse = 56.8 * mm;
  const G4double m_ANASENQQQ3CsIThickness = 26.0 * mm;
  const G4double m_ANASENQQQ3CsIWidthBot = 19.2 * mm;
  const G4double m_ANASENQQQ3CsIWidthTop = 41.3 * mm;
  const G4double m_ANASENQQQ3CsITotHyp =
      (m_ANASENQQQ3CsIhypotenuse * m_ANASENQQQ3CsIWidthTop) / (m_ANASENQQQ3CsIWidthTop - m_ANASENQQQ3CsIWidthBot);
  const G4double m_ANASENQQQ3CsITotHeight =
      sqrt(m_ANASENQQQ3CsITotHyp * m_ANASENQQQ3CsITotHyp - m_ANASENQQQ3CsIWidthTop * m_ANASENQQQ3CsIWidthTop);
  const G4double m_ANASENQQQ3CsIHeight =
      m_ANASENQQQ3CsITotHeight * (m_ANASENQQQ3CsIWidthTop - m_ANASENQQQ3CsIWidthBot) / m_ANASENQQQ3CsIWidthTop;
  // m_ANASENQQQ3CsITotHyp     106.147
  // m_ANASENQQQ3CsITotHeight  97.7825
  // m_ANASENQQQ3CsIHeight     52.3243

  ////////////////////////////////////////////////////////////

} // namespace STARKNS

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
STARK::STARK() {
  HCID_X6 = HCID_BB10 = HCID_QQQ5 = HCID_CsI = -1;

  m_X6 = m_BB10 = m_QQQ5 = m_Target = nullptr;
  for (auto i = 0; i < 4; ++i) {
    m_X6_wCsI[i] = nullptr;
    m_BB10_wCsI[i] = nullptr;
    m_QQQ5_wCsI[i] = nullptr;
  }
  m_X6Det = m_BB10Det = m_QQQ5Det = nullptr;
  m_logicSquareCsI1 = nullptr;
  m_logicSquareCsI2 = nullptr;
  m_logicANASENQQQ3CsI = nullptr;

  m_VisX6 = new G4VisAttributes(G4Colour(0., 0.5, 0.5));
  m_VisX6PCB = new G4VisAttributes(G4Colour(0.8, 0.5, 0.5));
  m_VisBB10 = new G4VisAttributes(G4Colour(0., 0.5, 0.7));
  m_VisBB10PCB = new G4VisAttributes(G4Colour(0.8, 0.5, 0.7));
  m_VisQQQ5 = new G4VisAttributes(G4Colour(0., 0.5, 0.3));
  m_VisQQQ5PCB = new G4VisAttributes(G4Colour(0.8, 0.5, 0.3));
  m_VisConn = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
  m_VisTarget = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.1));
  m_VisCsI = new G4VisAttributes(G4Colour(0.5, 0.5, 0.1, 0.5));
  m_VisCsI->SetForceWireframe(true);

  m_ReactionRegion = nullptr;
  m_useTarget = false;
  m_TargetMaterial = "";
  m_Pressure = 760.0;     // Torr
  m_Temperature = 293.15; // K
  m_Radius = 100.0;       // mm
  m_Z = 100.0;            // mm

  m_Event = new TSTARKData();
  m_Raw = new TSTARKRaw();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
STARK::~STARK() {}

void STARK::AddDetector(string Type, G4ThreeVector Pos, bool Flip, bool Rev, double Beta, int csi, int Group,
                        string mvName) {
  m_Type.push_back(Type);
  m_Pos.push_back(Pos);
  m_Rot.push_back(G4ThreeVector());
  m_AutoRotateFacingBeamAxis.push_back(false);
  m_Flip.push_back(Flip);
  m_Rev.push_back(Rev);
  m_Beta.push_back(Beta);
  m_Group.push_back(Group);
  m_CsI.push_back(csi);
  m_MVName.push_back(mvName);
}

void STARK::AddDetector(string Type, G4ThreeVector Pos, G4ThreeVector Rot, int csi, int Group, string mvName) {
  m_Type.push_back(Type);
  m_Pos.push_back(Pos);
  m_Rot.push_back(Rot);
  m_AutoRotateFacingBeamAxis.push_back(true);
  m_Flip.push_back(false);
  m_Rev.push_back(false);
  m_Beta.push_back(0);
  m_Group.push_back(Group);
  m_CsI.push_back(csi);
  m_MVName.push_back(mvName);
}

void STARK::BuildSquareCsI() {
  if (m_logicSquareCsI1)
    return;
  G4Box* solidCsI = new G4Box("solidCsI", CsI_X / 2., CsI_Y / 2., CsI_Z / 2.);
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* matCsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  m_logicSquareCsI1 = new G4LogicalVolume(solidCsI, matCsI, "logicCsI1", 0, 0, 0);
  m_logicSquareCsI2 = new G4LogicalVolume(solidCsI, matCsI, "logicCsI2", 0, 0, 0);
  m_logicSquareCsI1->SetVisAttributes(m_VisCsI);
  m_logicSquareCsI2->SetVisAttributes(m_VisCsI);
  m_logicSquareCsI1->SetSensitiveDetector(m_CsIDet);
  m_logicSquareCsI2->SetSensitiveDetector(m_CsIDet);
}

void STARK::BuildANASENQQQ3CsI() {
  if (m_logicANASENQQQ3CsI)
    return;
  G4Trap* solidANASENQQQ3CsI =
      new G4Trap("DetectorBody", 0.5 * m_ANASENQQQ3CsIHeight, 0, 0, 0.5 * m_ANASENQQQ3CsIThickness,
                 0.5 * m_ANASENQQQ3CsIWidthTop, 0.5 * m_ANASENQQQ3CsIWidthTop, 0, 0.5 * m_ANASENQQQ3CsIThickness,
                 0.5 * m_ANASENQQQ3CsIWidthBot, 0.5 * m_ANASENQQQ3CsIWidthBot, 0);
  G4NistManager* nist = G4NistManager::Instance();
  G4Material* matCsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  m_logicANASENQQQ3CsI = new G4LogicalVolume(solidANASENQQQ3CsI, matCsI, "logicCsI3", 0, 0, 0);
  m_logicANASENQQQ3CsI->SetVisAttributes(m_VisCsI);
  m_logicANASENQQQ3CsI->SetSensitiveDetector(m_CsIDet);
}

//////////////////////////////////////////////////////////////////////
//
// BuildX6Detector: making a G4AssemblyVolume for X6 (only once)
//
G4AssemblyVolume* STARK::BuildX6Detector(int buildCsI) {
  if (buildCsI == 0) {
    if (m_X6)
      return m_X6;
  }
  else {
    if (m_X6_wCsI[buildCsI])
      return m_X6_wCsI[buildCsI];
  }

  ////////////////////////////////////////////////////////////
  // material definition
  G4Material* matSi = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  ////////////////////////////////////////////////////////////

  G4Box* solidX6PCBAll = new G4Box("solidX6PCBAll", X6_PCBX / 2., X6_PCBY / 2., X6_PCBZ / 2.);
  G4Box* solidX6PCBSub1 = new G4Box("solidX6PCBSub1", X6_PCBSub1X / 2., X6_PCBSub1Y / 2.,
                                    X6_PCBSub1Z / 2. + 0.01 * mm); // +0.01 mm for the perfect subtraction of solid
  G4Box* solidX6PCBSub2 = new G4Box("solidX6PCBSub2", X6_PCBSub2X / 2., X6_PCBSub2Y / 2.,
                                    X6_PCBSub2Z / 2. + 0.01 * mm); // +0.01 mm for the perfect subtraction of solid

  G4VSolid* solidX6PCBTemp =
      new G4SubtractionSolid("solidX6PCBTemp", solidX6PCBAll, solidX6PCBSub1, 0,
                             G4ThreeVector(X6_PCBSub1XOffset, X6_PCBSub1YOffset, X6_PCBSub1ZOffset));
  G4VSolid* solidX6PCB = new G4SubtractionSolid("solidX6PCB", solidX6PCBTemp, solidX6PCBSub2, 0,
                                                G4ThreeVector(X6_PCBSub2XOffset, X6_PCBSub2YOffset, X6_PCBSub2ZOffset));

  G4LogicalVolume* logicX6PCB = new G4LogicalVolume(solidX6PCB, matPCB, "logicX6PCB", 0, 0, 0);
  logicX6PCB->SetVisAttributes(m_VisX6PCB);

  G4Box* solidX6Conn = new G4Box("solidX6Conn", Conn_X / 2., Conn_Y / 2., Conn_Z / 2.);
  G4LogicalVolume* logicX6Conn = new G4LogicalVolume(solidX6Conn, matPCB, "logicX6Conn", 0, 0, 0);
  logicX6Conn->SetVisAttributes(m_VisConn);

  G4Box* solidX6Si = new G4Box("solidX6Si", X6_SiX / 2., X6_SiY / 2., X6_SiZ / 2.);
  G4LogicalVolume* logicX6Si = new G4LogicalVolume(solidX6Si, matSi, "logicX6Si", 0, 0, 0);
  logicX6Si->SetVisAttributes(m_VisX6);
  logicX6Si->SetSensitiveDetector(m_X6Det);

  G4AssemblyVolume* assembly;
  assembly = new G4AssemblyVolume();
  G4ThreeVector Pos;
  assembly->AddPlacedVolume(logicX6Si, Pos, 0); // reference = center of the X6 Si wafer
  Pos = G4ThreeVector(-X6_SiXOffset, -X6_SiYOffset, -X6_SiZOffset);
  assembly->AddPlacedVolume(logicX6PCB, Pos, 0);
  Pos = G4ThreeVector(-X6_SiXOffset, -X6_SiYOffset - X6_PCBY / 2. + Conn_Y / 2., X6_PCBZ / 2. + Conn_Z / 2.);
  assembly->AddPlacedVolume(logicX6Conn, Pos, 0);

  if (buildCsI != 0) {
    BuildSquareCsI();
    G4double offZ = 0.5 * CsI_Z + CsI_X6_ZOffset;
    if (buildCsI == 2)
      offZ = -offZ;
    G4ThreeVector Pos1(CsI_X6_XOffset, CsI_X6_YOffset1, offZ);
    G4ThreeVector Pos2(CsI_X6_XOffset, CsI_X6_YOffset2, offZ);
    assembly->AddPlacedVolume(m_logicSquareCsI1, Pos1, 0);
    assembly->AddPlacedVolume(m_logicSquareCsI2, Pos2, 0);
  }

  if (buildCsI == 0)
    m_X6 = assembly;
  else
    m_X6_wCsI[buildCsI] = assembly;

  return assembly;
}

//////////////////////////////////////////////////////////////////////
//
// BuildBB10Detector: making a G4AssemblyVolume for BB10 (only once)
//
G4AssemblyVolume* STARK::BuildBB10Detector(int buildCsI) {
  if (buildCsI == 0) {
    if (m_BB10)
      return m_BB10;
  }
  else {
    if (m_BB10_wCsI[buildCsI])
      return m_BB10_wCsI[buildCsI];
  }

  ////////////////////////////////////////////////////////////
  // material definition
  G4Material* matSi = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  ////////////////////////////////////////////////////////////

  G4Box* solidBB10PCBAll = new G4Box("solidBB10PCBAll", BB10_PCBX / 2., BB10_PCBY / 2., BB10_PCBZ / 2.);
  G4Box* solidBB10PCBSub1 = new G4Box("solidBB10PCBSub1", BB10_PCBSub1X / 2., BB10_PCBSub1Y / 2.,
                                      BB10_PCBSub1Z / 2. + 0.01 * mm); // +0.01 mm for the perfect subtraction of solid
  G4Box* solidBB10PCBSub2 = new G4Box("solidBB10PCBSub2", BB10_PCBSub2X / 2., BB10_PCBSub2Y / 2.,
                                      BB10_PCBSub2Z / 2. + 0.01 * mm); // +0.01 mm for the perfect subtraction of solid

  G4VSolid* solidBB10PCBTemp =
      new G4SubtractionSolid("solidBB10PCBTemp", solidBB10PCBAll, solidBB10PCBSub1, 0,
                             G4ThreeVector(BB10_PCBSub1XOffset, BB10_PCBSub1YOffset, BB10_PCBSub1ZOffset));
  G4VSolid* solidBB10PCB =
      new G4SubtractionSolid("solidBB10PCB", solidBB10PCBTemp, solidBB10PCBSub2, 0,
                             G4ThreeVector(BB10_PCBSub2XOffset, BB10_PCBSub2YOffset, BB10_PCBSub2ZOffset));

  G4LogicalVolume* logicBB10PCB = new G4LogicalVolume(solidBB10PCB, matPCB, "logicBB10PCB", 0, 0, 0);
  logicBB10PCB->SetVisAttributes(m_VisBB10PCB);

  G4Box* solidBB10Conn = new G4Box("solidBB10Conn", Conn_X / 2., Conn_Y / 2., Conn_Z / 2.);
  G4LogicalVolume* logicBB10Conn = new G4LogicalVolume(solidBB10Conn, matPCB, "logicBB10Conn", 0, 0, 0);
  logicBB10Conn->SetVisAttributes(m_VisConn);

  G4Box* solidBB10Si = new G4Box("solidBB10Si", BB10_SiX / 2., BB10_SiY / 2., BB10_SiZ / 2.);
  G4LogicalVolume* logicBB10Si = new G4LogicalVolume(solidBB10Si, matSi, "logicBB10Si", 0, 0, 0);
  logicBB10Si->SetVisAttributes(m_VisBB10);
  logicBB10Si->SetSensitiveDetector(m_BB10Det);

  G4AssemblyVolume* assembly = new G4AssemblyVolume();
  G4ThreeVector Pos;
  assembly->AddPlacedVolume(logicBB10Si, Pos, 0); // reference = center of the BB10 Si wafer
  Pos = G4ThreeVector(-BB10_SiXOffset, -BB10_SiYOffset, -BB10_SiZOffset);
  assembly->AddPlacedVolume(logicBB10PCB, Pos, 0);
  Pos = G4ThreeVector(-BB10_SiXOffset, -BB10_SiYOffset - BB10_PCBY / 2. + Conn_Y / 2., BB10_PCBZ / 2. + Conn_Z / 2.);
  assembly->AddPlacedVolume(logicBB10Conn, Pos, 0);

  if (buildCsI != 0) {
    BuildSquareCsI();
    G4double offZ = 0.5 * CsI_Z + CsI_BB10_ZOffset;
    if (buildCsI == 2)
      offZ = -offZ;
    G4ThreeVector Pos1(CsI_BB10_XOffset, CsI_BB10_YOffset1, offZ);
    G4ThreeVector Pos2(CsI_BB10_XOffset, CsI_BB10_YOffset2, offZ);
    assembly->AddPlacedVolume(m_logicSquareCsI1, Pos1, 0);
    assembly->AddPlacedVolume(m_logicSquareCsI2, Pos2, 0);
  }

  if (buildCsI == 0)
    m_BB10 = assembly;
  else
    m_BB10_wCsI[buildCsI] = assembly;

  return assembly;
}

//////////////////////////////////////////////////////////////////////
//
// BuildQQQ5Detector: making a G4AssemblyVolume for QQQ5 (only once)
//
G4AssemblyVolume* STARK::BuildQQQ5Detector(int buildCsI) {
  if (m_QQQ5)
    return m_QQQ5;

  ////////////////////////////////////////////////////////////
  // material definition
  G4Material* matSi = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  ////////////////////////////////////////////////////////////

  // Make the a single detector geometry
  G4Tubs* solidQQQ5PCBAll =
      new G4Tubs("solidQQQ5PCBAll", QQQ5_PCBInR, QQQ5_PCBOutR, QQQ5_PCBT * 0.5, QQQ5_PCBPhi0, QQQ5_PCBPhi1);
  G4Box* solidQQQ5PCBCut = new G4Box("solidQQQ5PCBCut", QQQ5_PCBCutX / 2., QQQ5_PCBCutY / 2., QQQ5_PCBCutZ / 2.);
  G4VSolid* solidQQQ5PCBSub1 =
      new G4SubtractionSolid("solidQQQ5PCBSub1", solidQQQ5PCBAll, solidQQQ5PCBCut, 0,
                             G4ThreeVector(QQQ5_PCBCutXOffset, QQQ5_PCBCutYOffset, QQQ5_PCBCutZOffset));

  G4Tubs* solidQQQ5SiAllForPCBCut =
      new G4Tubs("solidQQQ5SiAllForPCBCut", QQQ5_SiInR, QQQ5_SiOutR, QQQ5_PCBCutZ / 2., QQQ5_SiPhi0, QQQ5_SiPhi1);
  G4Box* solidQQQ5SiCut1ForPCBCut =
      new G4Box("solidQQQ5SiCut1ForPCBCut", QQQ5_SiCut1X / 2., QQQ5_SiCut1Y / 2., QQQ5_PCBCutZ / 2. + 1 * mm);
  G4Box* solidQQQ5SiCut2ForPCBCut =
      new G4Box("solidQQQ5SiCut2ForPCBCut", QQQ5_SiCut2X / 2., QQQ5_SiCut2Y / 2., QQQ5_PCBCutZ / 2. + 1 * mm);
  G4VSolid* solidQQQ5SiSub1ForPCBCut =
      new G4SubtractionSolid("solidQQQ5SiSub1ForPCBCut", solidQQQ5SiAllForPCBCut, solidQQQ5SiCut1ForPCBCut, 0,
                             G4ThreeVector(QQQ5_SiCut1XOffset, QQQ5_SiCut1YOffset, QQQ5_SiCut1ZOffset));
  G4VSolid* solidQQQ5SiForPCBCut =
      new G4SubtractionSolid("solidQQQ5SiForPCBCut", solidQQQ5SiSub1ForPCBCut, solidQQQ5SiCut2ForPCBCut, 0,
                             G4ThreeVector(QQQ5_SiCut2XOffset, QQQ5_SiCut2YOffset, QQQ5_SiCut2ZOffset));
  G4VSolid* solidQQQ5PCB =
      new G4SubtractionSolid("solidQQQ5PCB", solidQQQ5PCBSub1, solidQQQ5SiForPCBCut, 0, G4ThreeVector(0, 0, 0));

  G4Tubs* solidQQQ5SiAll =
      new G4Tubs("solidQQQ5SiAll", QQQ5_SiInR, QQQ5_SiOutR, QQQ5_SiT * 0.5, QQQ5_SiPhi0, QQQ5_SiPhi1);
  G4Box* solidQQQ5SiCut1 = new G4Box("solidQQQ5SiCut1", QQQ5_SiCut1X / 2., QQQ5_SiCut1Y / 2., QQQ5_SiCut1Z / 2.);
  G4Box* solidQQQ5SiCut2 = new G4Box("solidQQQ5SiCut2", QQQ5_SiCut2X / 2., QQQ5_SiCut2Y / 2., QQQ5_SiCut2Z / 2.);
  G4VSolid* solidQQQ5SiSub1 =
      new G4SubtractionSolid("solidQQQ5SiSub1", solidQQQ5SiAll, solidQQQ5SiCut1, 0,
                             G4ThreeVector(QQQ5_SiCut1XOffset, QQQ5_SiCut1YOffset, QQQ5_SiCut1ZOffset));
  G4VSolid* solidQQQ5Si =
      new G4SubtractionSolid("solidQQQ5Si", solidQQQ5SiSub1, solidQQQ5SiCut2, 0,
                             G4ThreeVector(QQQ5_SiCut2XOffset, QQQ5_SiCut2YOffset, QQQ5_SiCut2ZOffset));

  G4LogicalVolume* logicQQQ5PCB = new G4LogicalVolume(solidQQQ5PCB, matPCB, "logicQQQ5PCB", 0, 0, 0);
  logicQQQ5PCB->SetVisAttributes(m_VisQQQ5PCB);

  G4Box* solidQQQ5Conn = new G4Box("solidQQQ5Conn", Conn_X / 2., Conn_Y / 2., Conn_Z / 2.);
  G4LogicalVolume* logicQQQ5Conn = new G4LogicalVolume(solidQQQ5Conn, matPCB, "logicQQQ5Conn", 0, 0, 0);
  logicQQQ5Conn->SetVisAttributes(m_VisConn);

  G4LogicalVolume* logicQQQ5Si = new G4LogicalVolume(solidQQQ5Si, matSi, "logicQQQ5", 0, 0, 0);
  logicQQQ5Si->SetVisAttributes(m_VisQQQ5);
  logicQQQ5Si->SetSensitiveDetector(m_QQQ5Det);

  G4AssemblyVolume* assembly = new G4AssemblyVolume();
  G4ThreeVector Pos;
  assembly->AddPlacedVolume(logicQQQ5Si, Pos, 0);
  assembly->AddPlacedVolume(logicQQQ5PCB, Pos, 0);
  Pos = G4ThreeVector(Conn_Y / 2., (QQQ5_SiOutR + QQQ5_SiInR) / 2., QQQ5_PCBT / 2. + Conn_Z / 2.);
  G4RotationMatrix* rot = new G4RotationMatrix;
  rot->rotateZ(90 * deg);
  assembly->AddPlacedVolume(logicQQQ5Conn, Pos, rot);

  if (buildCsI != 0) {
    BuildANASENQQQ3CsI();
    G4double m_ANASENQQQ3CsI_Z = 50;
    G4double offZ = 0.5 * m_ANASENQQQ3CsI_Z + CsI_QQQ5_ZOffset;
    if (buildCsI == 2)
      offZ = -offZ;
    for (auto iCsI : {0, 1, 2, 3}) {
      G4double angle1 = -360 / 32. * deg;
      G4double angle2 = -360 / 16. * deg;
      G4double angleZ = angle1 + iCsI * angle2;
      G4double angleO = 90 * deg + angle1 + iCsI * angle2;
      G4RotationMatrix* Rot = new G4RotationMatrix(0, 0, 0);
      Rot->rotateX(90 * deg);
      // Rot->rotateZ(180*deg);
      Rot->rotateZ(angleZ);
      G4double offR = CsI_QQQ5_YOffset + m_ANASENQQQ3CsITotHeight - m_ANASENQQQ3CsIHeight / 2.;
      G4ThreeVector Pos1(CsI_QQQ5_XOffset, 0, offZ);
      G4ThreeVector direction(std::cos(angleO), std::sin(angleO), 0);
      Pos1 = Pos1 + offR * direction;
      assembly->AddPlacedVolume(m_logicANASENQQQ3CsI, Pos1, Rot);
    }
  }

  if (buildCsI == 0)
    m_QQQ5 = assembly;
  else
    m_QQQ5_wCsI[buildCsI] = assembly;

  return assembly;
}

G4AssemblyVolume* STARK::BuildTarget() {
  if (m_Target)
    return m_Target;

  G4Material* matTarget =
      MaterialManager::getInstance()->GetGasFromLibrary(m_TargetMaterial, m_Pressure, m_Temperature);
  G4Tubs* solidTarget = new G4Tubs("solidTarget", 0, m_Radius / 2., m_Z / 2., 0, 360 * deg);
  m_logicTarget = new G4LogicalVolume(solidTarget, matTarget, "logicTarget", 0, 0, 0);
  m_logicTarget->SetVisAttributes(m_VisTarget);
  G4ThreeVector Pos;

  m_Target = new G4AssemblyVolume();
  m_Target->AddPlacedVolume(m_logicTarget, Pos, 0);
  return m_Target;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ReadConfiguration: reading a geometry file for each Si detector
void STARK::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("STARK");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> reso = {"Type", "Reso"}; // for put resolutions as an input parameter
  vector<string> targ = {"TargetMaterial", "Pressure", "Temperature", "Radius", "Z"};
  vector<string> cart = {"Type", "POS"};
  vector<string> sphe = {"Type", "R", "Theta", "Phi"};
  vector<string> cyld = {"Type", "Rho", "Phi", "Z"};
  // vector<string> car2 = {"Type", "POS", "RotateXYZ"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    ////////////////////////////////////////////////////////////
    // Resolution
    if (blocks[i]->HasTokenList(reso)) {
      string Type = blocks[i]->GetString("Type");
      if (Type.compare("X6") == 0)
        X6_ERes = blocks[i]->GetDouble("Reso", "void");
      else if (Type.compare("BB10") == 0)
        BB10_ERes = blocks[i]->GetDouble("Reso", "void");
      else if (Type.compare("QQQ5") == 0)
        QQQ5_ERes = blocks[i]->GetDouble("Reso", "void");
      else if (Type.compare("CsI") == 0)
        CsI_ERes = blocks[i]->GetDouble("Reso", "void");
      continue; // only block for resolution
    }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Target
    if (blocks[i]->HasTokenList(targ)) {
      m_useTarget = true;
      m_TargetMaterial = blocks[i]->GetString("TargetMaterial");
      m_Pressure = blocks[i]->GetDouble("Pressure", "Torr");
      m_Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
      m_Radius = blocks[i]->GetDouble("Radius", "mm");
      m_Z = blocks[i]->GetDouble("Z", "mm");
      cout << "////  TargetMaterial " << m_TargetMaterial << endl;
      cout << "////  Pressure " << m_Pressure << endl;
      cout << "////  Temperature " << m_Temperature << endl;
      cout << "////  Radius " << m_Radius << endl;
      cout << "////  Z " << m_Z << endl;
    }
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    // Cartesian coordinate or Spherical coordinate or Cylindrical coordinate
    else if (blocks[i]->HasTokenList(cart) || blocks[i]->HasTokenList(sphe) || blocks[i]->HasTokenList(cyld)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARK " << i + 1 << endl;
      string Type = blocks[i]->GetString("Type");
      G4ThreeVector Pos = G4ThreeVector();
      if (blocks[i]->HasTokenList(cart))
        Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS", "mm"));
      else if (blocks[i]->HasTokenList(sphe))
        Pos.setRThetaPhi(blocks[i]->GetDouble("R", "mm"), blocks[i]->GetDouble("Theta", "deg"),
                         blocks[i]->GetDouble("Phi", "deg"));
      else if (blocks[i]->HasTokenList(cyld))
        Pos.setRhoPhiZ(blocks[i]->GetDouble("Rho", "mm"), blocks[i]->GetDouble("Phi", "deg"),
                       blocks[i]->GetDouble("Z", "mm"));
      G4ThreeVector Rot = G4ThreeVector();
      bool AutoRotateFacingBeamAxis = true;
      bool Flip = false;
      bool Rev = false;
      double Beta = 0;
      ////////////////////////////////////////////////////////////
      // User defined rotation
      if (blocks[i]->HasToken("RotateXYZ")) {
        Rot = NPS::ConvertVector(blocks[i]->GetTVector3("RotateXYZ", "deg"));
        AutoRotateFacingBeamAxis = false;
      }
      else { // Auto rotate facing beam axis
        AutoRotateFacingBeamAxis = true;
        if (blocks[i]->HasToken("Flip"))
          Flip = blocks[i]->GetInt("Flip");
        if (blocks[i]->HasToken("Rev"))
          Rev = blocks[i]->GetInt("Rev");
        if (blocks[i]->HasToken("Beta"))
          Beta = blocks[i]->GetDouble("Beta", "deg");
      }
      ////////////////////////////////////////////////////////////
      int Group = 0;
      if (blocks[i]->HasToken("Group"))
        Group = blocks[i]->GetInt("Group");
      bool csi = false;
      if (blocks[i]->HasToken("CsI"))
        csi = blocks[i]->GetBool("CsI");
      string mvName;
      if (blocks[i]->HasToken("MotherVolume"))
        mvName = blocks[i]->GetString("MotherVolume");
      cout << endl << "////  Type " << Type << endl;
      cout << "////  Pos " << Pos << endl;
      cout << "////  Rot " << Rot << endl;
      cout << "////  AutoRotateFacingBeamAxis " << (AutoRotateFacingBeamAxis ? "True" : "False") << endl;
      cout << "////  Flip " << (Flip ? "True" : "False") << endl;
      cout << "////  Rev " << (Rev ? "True" : "False") << endl;
      cout << "////  Beta " << Beta << endl;
      cout << "////  CsI " << (csi ? "True" : "False") << endl;
      if (Group != 0) {
        cout << "////  Group " << Group << endl;
      }
      if (mvName.empty() == false) {
        cout << "////  MotherVolume " << mvName << endl;
      }
      if (AutoRotateFacingBeamAxis) {
        AddDetector(Type, Pos, Flip, Rev, Beta, csi, Group, mvName);
      }
      else {
        AddDetector(Type, Pos, Rot, csi, Group, mvName);
      }
    }
    ////////////////////////////////////////////////////////////
    else {
      cout << "Error: check your input file formatting" << endl;
      exit(1);
    }
  }
  std::cout << "read complete" << std::endl;
}

void STARK::SetReactionRegion(G4LogicalVolume* world) {
  if (m_useTarget) {
    if (!m_ReactionRegion) {
      G4ProductionCuts* productionCuts = new G4ProductionCuts();
      productionCuts->SetProductionCut(1000 * mm, "e-");
      m_ReactionRegion = new G4Region("NPSimulationProcess");
      m_ReactionRegion->SetProductionCuts(productionCuts);
      m_ReactionRegion->AddRootLogicalVolume(m_logicTarget);
      m_ReactionRegion->SetUserLimits(new G4UserLimits(0.5 * mm));
    }
    new NPS::BeamReaction("BeamReaction", m_ReactionRegion);
    new NPS::Decay("Decay", m_ReactionRegion);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::ConstructDetector(G4LogicalVolume* world) {
  std::cout << "start constuct detector" << std::endl;
  for (unsigned short i = 0; i < m_Pos.size(); i++) {
    G4double phi = m_Pos[i].getPhi();
    G4RotationMatrix* Rot = new G4RotationMatrix(0, 0, 0);
    if (m_AutoRotateFacingBeamAxis[i]) {
      Rot->rotateX(m_Rot[i].x());
      Rot->rotateY(m_Rot[i].y());
      Rot->rotateZ(m_Rot[i].z());
    }
    else {
      if (m_Rev[i])
        Rot->rotateZ(180 * deg);
      if (m_Flip[i])
        Rot->rotateY(180 * deg);
      Rot->rotateX(90 * deg);
      Rot->rotateZ(90 * deg + phi);
    }

    G4AssemblyVolume* det;
    if (m_Type[i] == "X6")
      det = BuildX6Detector(m_CsI[i]);
    else if (m_Type[i] == "BB10")
      det = BuildBB10Detector(m_CsI[i]);
    else if (m_Type[i] == "QQQ5") {
      det = BuildQQQ5Detector(m_CsI[i]);
      Rot = new G4RotationMatrix;
      Rot->rotateZ(m_Beta[i]);
      if (m_Flip[i])
        Rot->rotateY(180 * deg);
    }
    else {
      std::cerr << "no type " << m_Type[i] << " exists.\n";
      continue;
    }

    if (m_MVName[i].empty() == false) {
      for (const auto& lv : *G4LogicalVolumeStore::GetInstance()) {
        if (std::string(lv->GetName()) == m_MVName[i]) {
          world = lv;
        }
      }
    }
    det->MakeImprint(world, m_Pos[i], Rot, i + 1, true);

    // iterator is equal to fPVStore.begin()
    std::vector<G4VPhysicalVolume*>::iterator it = det->GetVolumesIterator();
    unsigned int NbrImprints = det->GetImprintsCount();
    unsigned int NbrTotalPV = det->TotalImprintedVolumes();
    unsigned int NbrComponents = NbrTotalPV / NbrImprints;
    // set copy numbers of components of assembly volume to the current detector number
    int countComponents = 0;
    for (it += (NbrImprints - 1) * NbrComponents; it <= det->GetVolumesIterator() + NbrTotalPV - 1; it++)
      (*it)->SetCopyNo(i + 1 + (countComponents++) * 100);
  }

  if (m_useTarget) {
    G4AssemblyVolume* target = BuildTarget();
    G4ThreeVector Pos;
    target->MakeImprint(world, Pos, nullptr, 0, true);
    SetReactionRegion(world);
  }
  std::cout << "construct complete" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::InitializeRootOutput() {
  // add detector branch to the EventTree
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("STARK"))
    pTree->Branch("STARK", "TSTARKData", &m_Event);
  if (!pTree->FindBranch("Raw"))
    pTree->Branch("Raw", "TSTARKRaw", &m_Raw);
  pTree->SetBranchAddress("STARK", &m_Event);
  pTree->SetBranchAddress("Raw", &m_Raw);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::ReadSensitive(const G4Event* event) {
  // read sensitive part and fill the Root tree
  m_Event->Clear();
  m_Raw->Clear();

  auto HCE = event->GetHCofThisEvent();
  if (!HCE)
    return;

  if (HCID_X6 == -1)
    HCID_X6 = G4SDManager::GetSDMpointer()->GetCollectionID("X6Det/X6Scorer");
  if (HCID_BB10 == -1)
    HCID_BB10 = G4SDManager::GetSDMpointer()->GetCollectionID("BB10Det/BB10Scorer");
  if (HCID_QQQ5 == -1)
    HCID_QQQ5 = G4SDManager::GetSDMpointer()->GetCollectionID("QQQ5Det/QQQ5Scorer");
  if (HCID_CsI == -1)
    HCID_CsI = G4SDManager::GetSDMpointer()->GetCollectionID("CsIDet/CsIScorer");

  /////////////////////////////////////////////////////////////////////////////////
  // loop for the event map
  /////////////////////////////////////
  //  for X6
  auto evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_X6));
  map<G4int, G4double**>::iterator it;
  for (it = evtMap->GetMap()->begin(); it != evtMap->GetMap()->end(); it++) {
    // energy smearing
    G4double enSmear2 = RandGauss::shoot((*(it->second))[2], (*(it->second))[2] * X6_ERes / 100.);
    G4double enSmear0 = RandGauss::shoot((*(it->second))[0], (*(it->second))[0] * X6_ERes / 100.);
    G4double enSmear1 = RandGauss::shoot((*(it->second))[1], (*(it->second))[1] * X6_ERes / 100.);
    m_Event->Set(0,
                 (*(it->second))[4],  // detector number
                 (*(it->second))[5],  // front strip number
                 (*(it->second))[6],  // back strip number
                 enSmear2,            // frontside energy
                 enSmear2,            // backside energy
                 enSmear0,            // upstream energy
                 enSmear1,            // downstream energy
                 (*(it->second))[3],  // global time
                 (*(it->second))[7],  // hit position X
                 (*(it->second))[8],  // hit position Y
                 (*(it->second))[9]); // hit position Z
    m_Raw->SetRaw(m_Event);
  }
  /////////////////////////////////////
  // for BB10
  evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_BB10));
  for (it = evtMap->GetMap()->begin(); it != evtMap->GetMap()->end(); it++) {
    // energy smearing
    G4double enSmear0 = RandGauss::shoot((*(it->second))[0], (*(it->second))[0] * BB10_ERes / 100.);
    G4double enSmear1 = RandGauss::shoot((*(it->second))[1], (*(it->second))[1] * BB10_ERes / 100.);
    m_Event->Set(1,
                 (*(it->second))[3],  // detector number
                 (*(it->second))[4],  // front strip number
                 1,                   // back strip number (only 1)
                 enSmear0,            // frontside energy
                 enSmear1,            // backside energy
                 0,                   // upstream energy
                 0,                   // downstream energy
                 (*(it->second))[2],  // global time
                 (*(it->second))[5],  // hit position X
                 (*(it->second))[6],  // hit position Y
                 (*(it->second))[7]); // hit position Z
  }
  /////////////////////////////////////
  // for QQQ5
  evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_QQQ5));
  for (it = evtMap->GetMap()->begin(); it != evtMap->GetMap()->end(); it++) {
    // energy smearing
    G4double enSmear0 = RandGauss::shoot((*(it->second))[0], (*(it->second))[0] * QQQ5_ERes / 100.);
    G4double enSmear1 = RandGauss::shoot((*(it->second))[1], (*(it->second))[1] * QQQ5_ERes / 100.);
    m_Event->Set(2,
                 (*(it->second))[3],  // detector number
                 (*(it->second))[4],  // front strip number
                 (*(it->second))[5],  // back strip number
                 enSmear0,            // frontside energy
                 enSmear1,            // backside energy
                 0,                   // upstream energy
                 0,                   // downstream energy
                 (*(it->second))[2],  // global time
                 (*(it->second))[6],  // hit position X
                 (*(it->second))[7],  // hit position Y
                 (*(it->second))[8]); // hit position Z
  }
  /////////////////////////////////////
  // for CsI
  evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_CsI));
  for (it = evtMap->GetMap()->begin(); it != evtMap->GetMap()->end(); it++) {
    // energy smearing
    G4double enSmear0 = RandGauss::shoot((*(it->second))[1], (*(it->second))[1] * CsI_ERes / 100.);
    m_Event->Set(3,
                 (*(it->second))[0], // detector number
                 0,                  // front strip number
                 0,                  // back strip number
                 enSmear0,           // frontside energy
                 0,                  // backside energy
                 0,                  // upstream energy
                 0,                  // downstream energy
                 (*(it->second))[2], // global time
                 0,                  // hit position X
                 0,                  // hit position Y
                 0);                 // hit position Z
  }
  /////////////////////////////////////////////////////////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::InitializeScorers() {
  // Check the detectors initialized
  bool already_exist_X6 = false;
  bool already_exist_BB10 = false;
  bool already_exist_QQQ5 = false;
  bool already_exist_CsI = false;
  m_X6Det = CheckScorer("X6Det", already_exist_X6);       // MultiFunctionalDetector
  m_BB10Det = CheckScorer("BB10Det", already_exist_BB10); // MultiFunctionalDetector
  m_QQQ5Det = CheckScorer("QQQ5Det", already_exist_QQQ5); // MultiFunctionalDetector
  m_CsIDet = CheckScorer("CsIDet", already_exist_CsI);    // MultiFunctionalDetector

  // if not, create them
  if (!already_exist_X6) {
    G4VPrimitiveScorer* X6Scorer =
        new STARKSCORERS::PS_STARK_X6("X6Scorer", 0, X6_SiActiveX, X6_SiActiveY, X6_NFrontStrips, X6_NBackStrips, 0);
    m_X6Det->RegisterPrimitive(X6Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_X6Det);
  }

  if (!already_exist_BB10) {
    G4VPrimitiveScorer* BB10Scorer =
        new STARKSCORERS::PS_STARK_BB10("BB10Scorer", 0, BB10_SiActiveX, BB10_SiActiveY, BB10_NFrontStrips, 0);
    m_BB10Det->RegisterPrimitive(BB10Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_BB10Det);
  }

  if (!already_exist_QQQ5) {
    G4VPrimitiveScorer* QQQ5Scorer = new STARKSCORERS::PS_STARK_QQQ5(
        "QQQ5Scorer", 0, QQQ5_SiActiveInR, QQQ5_SiActiveOutR, QQQ5_SiPhi0, QQQ5_SiPhi1, QQQ5_NAStrip, QQQ5_NRStrip, 0);
    m_QQQ5Det->RegisterPrimitive(QQQ5Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_QQQ5Det);
  }

  if (!already_exist_CsI) {
    G4VPrimitiveScorer* CsIScorer = new STARKSCORERS::PS_STARK_CsI("CsIScorer", 0, 0);
    m_CsIDet->RegisterPrimitive(CsIScorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_CsIDet);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// construct method to be passed to the DetectorFactory
NPS::VDetector* STARK::Construct() { return (NPS::VDetector*)new STARK(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// register the construct method to the factory
extern "C" {
class proxy_nps_STARK {
 public:
  proxy_nps_STARK() {
    NPS::DetectorFactory::getInstance()->AddToken("STARK", "STARK");
    NPS::DetectorFactory::getInstance()->AddDetector("STARK", STARK::Construct);
  }
};
proxy_nps_STARK p_nps_STARK;
}
