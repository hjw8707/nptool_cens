/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  TOGAXSI_SI simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Trap.hh"

#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4TwoVector.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// NPTool header
#include "TOGAXSI_SI.hh"
//#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"

#include "DSSDScorers.hh"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace TOGAXSI_SI_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 1*keV;
  const double ResoEnergy = 0.015*MeV ;

  ////////////////////
  // Inner Detector //
  ////////////////////
  // WaferX Parameter
  double InnerX_Wafer_Length = 78.4*mm ; 
  double InnerX_Wafer_Width = 51*mm ;
  double InnerX_Wafer_Thickness = 0.1*mm ;
  double InnerX_Wafer_LongitudinalStrips = 255;
//  double InnerX_Wafer_LongitudinalStrips = 64;
  double InnerX_Wafer_TransverseStrips = 1;
//  double InnerX_Wafer_LongitudinalStrips = 128;
//  double InnerX_Wafer_TransverseStrips = 128;

  // PCB parameter
  double InnerX_PCB_Thickness=3*mm;

  double InnerX_PCB_Length = 84.9*mm;
  double InnerX_PCB_Width = 55.83*mm;

  // WaferZ Parameter
  double InnerZ_Wafer_Length = 78.4*mm ; 
  double InnerZ_Wafer_Width = 51*mm ;
  double InnerZ_Wafer_Thickness = 0.1*mm ;
  double InnerZ_Wafer_LongitudinalStrips = 1;
  double InnerZ_Wafer_TransverseStrips = 392;
//  double InnerZ_Wafer_TransverseStrips = 98;
//  double InnerZ_Wafer_LongitudinalStrips = 128;
//  double InnerZ_Wafer_TransverseStrips = 128;

  // PCB parameter
  double InnerZ_PCB_Thickness=3*mm;

  double InnerZ_PCB_Length = 84.9*mm;
  double InnerZ_PCB_Width = 55.83*mm;

  ////////////////////
  // Outer Detector //
  ////////////////////
  // WaferX Parameter
  double OuterX_Wafer_Length = 78.4*mm ; 
  double OuterX_Wafer_Width = 51*mm ;
//  double OuterX_Wafer_Length = 51*mm ; 
//  double OuterX_Wafer_Width = 78.4*mm ;
  double OuterX_Wafer_Thickness = 0.1*mm ;
  double OuterX_Wafer_LongitudinalStrips = 392;
//  double OuterX_Wafer_LongitudinalStrips = 98;
  double OuterX_Wafer_TransverseStrips = 1;
//  double OuterX_Wafer_LongitudinalStrips = 1;
//  double OuterX_Wafer_TransverseStrips = 392;

  // PCB parameter
  double OuterX_PCB_Thickness=3*mm;

  double OuterX_PCB_Length = 107.16*mm;
  double OuterX_PCB_Width = 84.9*mm;
  double OuterX_PCB_gap = 3.5*mm;


  // WaferZ Parameter
  double OuterZ_Wafer_Length = 78.4*mm ; 
  double OuterZ_Wafer_Width = 51*mm ;
//  double OuterZ_Wafer_Length = 51*mm ; 
//  double OuterZ_Wafer_Width = 78.4*mm ;
  double OuterZ_Wafer_Thickness = 0.1*mm ;
  double OuterZ_Wafer_LongitudinalStrips = 1;
  double OuterZ_Wafer_TransverseStrips = 255;
//  double OuterZ_Wafer_TransverseStrips = 64;
//  double OuterZ_Wafer_LongitudinalStrips = 255;
//  double OuterZ_Wafer_TransverseStrips = 1;

  // PCB parameter
  double OuterZ_PCB_Thickness=3*mm;

  double OuterZ_PCB_Length = 108.6*mm;
  double OuterZ_PCB_Width = 83.3*mm;
  double OuterZ_PCB_gap = 3.5*mm;


  /////////////////////
  // ClusterDetector //
  /////////////////////
  double ClusterInner_Wafer_Base = 56.*mm; 
  double ClusterInner_Wafer_Top = 7.*mm;
  double ClusterInner_Wafer_Height = 80.*mm;
  double ClusterInner_Wafer_Thickness = 0.1*mm;
 
  double ClusterInner_ActiveWafer_Base = 88.*mm; 
  double ClusterInner_ActiveWafer_Top = 12.*mm;
  double ClusterInner_ActiveWafer_Height = 66.*mm;
  double ClusterInner_ActiveWafer_Thickness = 0.1*mm;

  double ClusterInner_Wafer_LongitudinalStrips = 128;   
  double ClusterInner_Wafer_TransverseStrips = 128;   

  //Cluster Demonstrator X
  double ClusterX_Wafer_Length = 78.4*mm ; 
  double ClusterX_Wafer_Width = 51*mm ;
  double ClusterX_Wafer_Thickness = 0.1*mm ;
  double ClusterX_Wafer_LongitudinalStrips = 392;
  double ClusterX_Wafer_TransverseStrips = 1;

  double ClusterX_PCB_Thickness=3*mm;
  double ClusterX_PCB_Length = 107.16*mm;
  double ClusterX_PCB_Width = 84.9*mm;
  double ClusterX_PCB_gap = 3.5*mm;

  //Cluster Demonstrator Y
  double ClusterY_Wafer_Length = 78.4*mm ; 
  double ClusterY_Wafer_Width = 51*mm ;
  double ClusterY_Wafer_Thickness = 0.1*mm ;
  double ClusterY_Wafer_LongitudinalStrips = 1;
  double ClusterY_Wafer_TransverseStrips = 255;

  double ClusterY_PCB_Thickness=3*mm;
  double ClusterY_PCB_Length = 108.6*mm;
  double ClusterY_PCB_Width = 83.3*mm;
  double ClusterY_PCB_gap = 3.5*mm;

}

using namespace TOGAXSI_SI_NS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// TOGAXSI_SI Specific Method
TOGAXSI_SI::TOGAXSI_SI(){
  InitializeMaterial();
  m_ReactionRegion = NULL;
  m_Event = new TTOGAXSI_SIData() ;
//  m_ReactionRegion = NULL;
  m_InnerXScorer = 0;
  m_InnerXDetector = 0;
  m_InnerZScorer = 0;
  m_InnerZDetector = 0;

  m_OuterXScorer1 = 0;
  m_OuterXScorer2 = 0;
  m_OuterZScorer1 = 0;
  m_OuterZScorer2 = 0;
  m_OuterXDetector = 0; 
  m_OuterZDetector = 0; 

  m_ClusterInnerScorer = 0;
  m_ClusterInnerDetector = 0;

  m_ClusterX1Scorer1 = 0;
  m_ClusterX1Scorer2 = 0;
  m_ClusterX1Detector = 0;
  m_ClusterY1Scorer1 = 0;
  m_ClusterY1Scorer2 = 0;
  m_ClusterY1Detector = 0;
  m_ClusterX2Scorer1 = 0;
  m_ClusterX2Scorer2 = 0;
  m_ClusterX2Detector = 0;
  m_ClusterY2Scorer1 = 0;
  m_ClusterY2Scorer2 = 0;
  m_ClusterY2Detector = 0;

  m_Target = 0;
  m_TargetCell = 0;
  m_TargetCellScorer = 0;

  // Dark Grey
  m_VisSi = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));  
  // Green
  m_VisPCB = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2));

  //Transparent blue
  m_VisTarget = new G4VisAttributes(G4Colour(0.15,0.85,0.85,0.1));
  m_VisTargetCell = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.5));

}

TOGAXSI_SI::~TOGAXSI_SI(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddInnerXDetector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_InnerX_R.push_back(R);
  m_InnerX_Z.push_back(Z);
  m_InnerX_Phi.push_back(Phi);
  m_InnerX_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddInnerZDetector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_InnerZ_R.push_back(R);
  m_InnerZ_Z.push_back(Z);
  m_InnerZ_Phi.push_back(Phi);
  m_InnerZ_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddOuterXDetector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_OuterX_R.push_back(R);
  m_OuterX_Z.push_back(Z);
  m_OuterX_Phi.push_back(Phi);
  m_OuterX_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddOuterZDetector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_OuterZ_R.push_back(R);
  m_OuterZ_Z.push_back(Z);
  m_OuterZ_Phi.push_back(Phi);
  m_OuterZ_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddClusterInnerDetector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_ClusterInner_R.push_back(R);
  m_ClusterInner_Z.push_back(Z);
  m_ClusterInner_Phi.push_back(Phi);
  m_ClusterInner_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddClusterX1Detector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_ClusterX1_R.push_back(R);
  m_ClusterX1_Z.push_back(Z);
  m_ClusterX1_Phi.push_back(Phi);
  m_ClusterX1_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddClusterY1Detector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_ClusterY1_R.push_back(R);
  m_ClusterY1_Z.push_back(Z);
  m_ClusterY1_Phi.push_back(Phi);
  m_ClusterY1_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddClusterX2Detector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_ClusterX2_R.push_back(R);
  m_ClusterX2_Z.push_back(Z);
  m_ClusterX2_Phi.push_back(Phi);
  m_ClusterX2_Ref.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddClusterY2Detector(double R, double Z, double Phi, G4ThreeVector Ref) {
  m_ClusterY2_R.push_back(R);
  m_ClusterY2_Z.push_back(Z);
  m_ClusterY2_Phi.push_back(Phi);
  m_ClusterY2_Ref.push_back(Ref);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_SI::AddTarget(double R, double L, string MaterialName, string CellMaterialName, double CellThickness, G4ThreeVector Pos) {
  
  m_Target_R.push_back(R);
  m_Target_L.push_back(L);
  m_Target_MaterialName.push_back(MaterialName);
  m_Target_CellMaterialName.push_back(CellMaterialName);
  m_Target_CellThickness.push_back(CellThickness);
  m_Target_Pos.push_back(Pos);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildInnerXDetector(){
  if(!m_InnerXDetector){
    
    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * InnerX_PCB_Width, 0.5 * InnerX_PCB_Thickness, 0.5 * InnerX_PCB_Length);
  
  
    // Master volume inner Detector
    m_InnerXDetector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_InnerXDetector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", InnerX_Wafer_Width * 0.5, InnerX_PCB_Thickness * 0.5 + 0.1 * mm, InnerX_Wafer_Length * 0.5);


   G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCBFull, HoleShape);

   //Sub Volume PCB
   G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
   logicPCB->SetVisAttributes(m_VisPCB);
   
   new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_InnerX_PCB", m_InnerXDetector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * InnerX_Wafer_Width, 0.5 * InnerX_Wafer_Thickness, 0.5 * InnerX_Wafer_Length);

    G4LogicalVolume* logicWafer = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer->SetVisAttributes(m_VisSi);
//    logicWafer->SetSensitiveDetector(m_InnerScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicWafer, "TOGAXSI_SI_InnerX_Wafer", m_InnerXDetector, false, 1);



   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("InnerActiveXWaferShape",0.5 * InnerX_Wafer_Width, 0.5 * InnerX_Wafer_Thickness, 0.5 * InnerX_Wafer_Length);

    G4LogicalVolume* logicActiveWafer = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer", 0, 0, 0);
    logicActiveWafer->SetVisAttributes(m_VisSi);
    logicActiveWafer->SetSensitiveDetector(m_InnerXScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer, "TOGAXSI_SI_Active_InnerX_Wafer", logicWafer, false, 1);


  }
  return m_InnerXDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildInnerZDetector(){
  if(!m_InnerZDetector){
    
    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * InnerZ_PCB_Width, 0.5 * InnerZ_PCB_Thickness, 0.5 * InnerZ_PCB_Length);
  
  
    // Master volume inner Detector
    m_InnerZDetector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_InnerZDetector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", InnerZ_Wafer_Width * 0.5, InnerZ_PCB_Thickness * 0.5 + 0.1 * mm, InnerZ_Wafer_Length * 0.5);


   G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCBFull, HoleShape);

   //Sub Volume PCB
   G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
   logicPCB->SetVisAttributes(m_VisPCB);
   
   new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_InnerZ_PCB", m_InnerZDetector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * InnerZ_Wafer_Width, 0.5 * InnerZ_Wafer_Thickness, 0.5 * InnerZ_Wafer_Length);

    G4LogicalVolume* logicWafer = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer->SetVisAttributes(m_VisSi);
//    logicWafer->SetSensitiveDetector(m_InnerScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicWafer, "TOGAXSI_SI_InnerZ_Wafer", m_InnerZDetector, false, 1);



   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("InnerActiveXWaferShape",0.5 * InnerZ_Wafer_Width, 0.5 * InnerZ_Wafer_Thickness, 0.5 * InnerZ_Wafer_Length);

    G4LogicalVolume* logicActiveWafer = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer", 0, 0, 0);
    logicActiveWafer->SetVisAttributes(m_VisSi);
    logicActiveWafer->SetSensitiveDetector(m_InnerZScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer, "TOGAXSI_SI_Active_InnerZ_Wafer", logicWafer, false, 1);


  }
  return m_InnerZDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildOuterXDetector(){
  if(!m_OuterXDetector){


    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * OuterX_PCB_Width, 0.5 * OuterX_PCB_Thickness, 0.5 * OuterX_PCB_Length);
  
  
    // Master volume inner Detector
    m_OuterXDetector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_OuterXDetector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", OuterX_Wafer_Length * 0.5, OuterX_PCB_Thickness*0.5 + 0.1 * mm, OuterX_Wafer_Width * 0.5);

    G4ThreeVector HoleShift = G4ThreeVector(0,0,0.5 * OuterX_Wafer_Width + 0.5 * OuterX_PCB_gap);

    G4SubtractionSolid* PCB_temp = new G4SubtractionSolid("PCB", PCBFull, HoleShape, new G4RotationMatrix,HoleShift);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_temp, HoleShape,new G4RotationMatrix,-HoleShift);


    //Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(m_VisPCB);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_OuterX_PCB", m_OuterXDetector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * OuterX_Wafer_Length, 0.5 * OuterX_Wafer_Thickness, 0.5 * OuterX_Wafer_Width);

    G4LogicalVolume* logicWafer1 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer1->SetVisAttributes(m_VisSi);
//    logicWafer1->SetSensitiveDetector(m_OuterScorer1);
    G4LogicalVolume* logicWafer2 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer2->SetVisAttributes(m_VisSi);
//    logicWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),-HoleShift, logicWafer1, "TOGAXSI_SI_OuterX_Wafer 1", m_OuterXDetector, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),HoleShift, logicWafer2, "TOGAXSI_SI_OuterX_Wafer 2", m_OuterXDetector, false, 1);


   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("OuterXActiveWaferShape",0.5 * OuterX_Wafer_Length, 0.5 * OuterX_Wafer_Thickness, 0.5 * OuterX_Wafer_Width);

    G4LogicalVolume* logicActiveWafer1 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(m_VisSi);
    logicActiveWafer1->SetSensitiveDetector(m_OuterXScorer1);
    G4LogicalVolume* logicActiveWafer2 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(m_VisSi);
    logicActiveWafer2->SetSensitiveDetector(m_OuterXScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer1, "TOGAXSI_SI_Inner_Wafer", logicWafer1, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer2, "TOGAXSI_SI_Inner_Wafer", logicWafer2, false, 1);


  }
  return m_OuterXDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildOuterZDetector(){
  if(!m_OuterZDetector){


    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * OuterZ_PCB_Width, 0.5 * OuterZ_PCB_Thickness, 0.5 * OuterZ_PCB_Length);
  
  
    // Master volume inner Detector
    m_OuterZDetector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_OuterZDetector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", OuterZ_Wafer_Length * 0.5, OuterZ_PCB_Thickness * 0.5 + 0.1 * mm, OuterZ_Wafer_Width * 0.5);

    G4ThreeVector HoleShift = G4ThreeVector(0,0,0.5 * OuterZ_Wafer_Width + 0.5 * OuterZ_PCB_gap);

    G4SubtractionSolid* PCB_temp = new G4SubtractionSolid("PCB", PCBFull, HoleShape, new G4RotationMatrix,HoleShift);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_temp, HoleShape,new G4RotationMatrix,-HoleShift);


    //Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(m_VisPCB);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_OuterZ_PCB", m_OuterZDetector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * OuterZ_Wafer_Length, 0.5 * OuterZ_Wafer_Thickness, 0.5 * OuterZ_Wafer_Width);

    G4LogicalVolume* logicWafer1 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer1->SetVisAttributes(m_VisSi);
//    logicWafer1->SetSensitiveDetector(m_OuterScorer1);
    G4LogicalVolume* logicWafer2 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer2->SetVisAttributes(m_VisSi);
//    logicWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),-HoleShift, logicWafer1, "TOGAXSI_SI_OuterZ_Wafer 1", m_OuterZDetector, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),HoleShift, logicWafer2, "TOGAXSI_SI_OuterZ_Wafer 2", m_OuterZDetector, false, 1);


   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("OuterZActiveWaferShape",0.5 * OuterZ_Wafer_Length, 0.5 * OuterZ_Wafer_Thickness, 0.5 * OuterZ_Wafer_Width);

    G4LogicalVolume* logicActiveWafer1 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(m_VisSi);
    logicActiveWafer1->SetSensitiveDetector(m_OuterZScorer1);
    G4LogicalVolume* logicActiveWafer2 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(m_VisSi);
    logicActiveWafer2->SetSensitiveDetector(m_OuterZScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer1, "TOGAXSI_SI_Inner_Wafer", logicWafer1, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer2, "TOGAXSI_SI_Inner_Wafer", logicWafer2, false, 1);


  }
  return m_OuterZDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildClusterInnerDetector(){
  if(!m_ClusterInnerDetector){

    vector<G4TwoVector> coord;

    coord.push_back(G4TwoVector(0.5 * ClusterInner_Wafer_Base,0.5 * ClusterInner_Wafer_Height)); 
    coord.push_back(G4TwoVector(-0.5 * ClusterInner_Wafer_Base, 0.5 * ClusterInner_Wafer_Height)); 
    coord.push_back(G4TwoVector(-0.5 * ClusterInner_Wafer_Top, -0.5 * ClusterInner_Wafer_Height)); 
    coord.push_back(G4TwoVector(0.5 * ClusterInner_Wafer_Top, -0.5 * ClusterInner_Wafer_Height)); 

    cout << "_____________Coordinates______________________" << endl;
    cout << coord[0][0] << endl;
    cout << coord[1][0] << endl;
    cout << coord[2][0] << endl;
    cout << coord[3][0] << endl;
    cout << coord[0][1] << endl;
    cout << coord[1][1] << endl;
    cout << coord[2][1] << endl;
    cout << coord[3][1] << endl;

    // for (int i = 0; i++; i<4) cout << coord[i][0] << endl;
    //for (int i = 0; i++; i<4) cout << coord[i][1] << endl;

    G4ExtrudedSolid* WaferShape = new G4ExtrudedSolid("TOGAXSI_SI_ClusterInnerWaferShape", coord,  0.5 * ClusterInner_Wafer_Thickness, G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);
   
    // Master volume ClusterInner Detector
    m_ClusterInnerDetector = new G4LogicalVolume(WaferShape,m_MaterialVacuum,"logicClusterInnerDetector",0,0,0);
    m_ClusterInnerDetector->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Sub Volume Wafer
    G4LogicalVolume* logicWaferShape = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafershape", 0, 0, 0);
    logicWaferShape->SetVisAttributes(m_VisSi);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicWaferShape, "TOGAXSI_ClusterInner", m_ClusterInnerDetector, false, 0);

   // Active Wafer

    vector<G4TwoVector> coord_active;

    coord_active.push_back(G4TwoVector(0.5 * ClusterInner_Wafer_Base,0.5 * ClusterInner_Wafer_Height)); 
    coord_active.push_back(G4TwoVector(-0.5 * ClusterInner_Wafer_Base, 0.5 * ClusterInner_Wafer_Height)); 
    coord_active.push_back(G4TwoVector(-0.5 * ClusterInner_Wafer_Top, -0.5 * ClusterInner_Wafer_Height)); 
    coord_active.push_back(G4TwoVector(0.5 * ClusterInner_Wafer_Top, -0.5 * ClusterInner_Wafer_Height)); 

    G4ExtrudedSolid* ActiveWaferShape = new G4ExtrudedSolid("ClusterInnerActiveWaferShape", coord_active, 0.5 * ClusterInner_ActiveWafer_Thickness, G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    G4LogicalVolume* logicActiveWafer = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer->SetVisAttributes(m_VisSi);
//    logicActiveWafer->SetSensitiveDetector(m_ClusterScorer);
  
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer, "TOGAXSI_SI_Inner_Wafer", logicWaferShape, false, 1);


  }
  return m_ClusterInnerDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildClusterX1Detector(){
  if(!m_ClusterX1Detector){


    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * ClusterX_PCB_Width, 0.5 * ClusterX_PCB_Length, 0.5 * ClusterX_PCB_Thickness);
  
  
    // Master volume inner Detector
    m_ClusterX1Detector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_ClusterX1Detector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", ClusterX_Wafer_Length * 0.5, ClusterX_Wafer_Width * 0.5, ClusterX_PCB_Thickness * 0.5 + 0.1 * mm);

    G4ThreeVector HoleShift = G4ThreeVector(0,0.5 * ClusterX_Wafer_Width + 0.5 * ClusterX_PCB_gap,0);

    G4SubtractionSolid* PCB_temp = new G4SubtractionSolid("PCB", PCBFull, HoleShape, new G4RotationMatrix,HoleShift);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_temp, HoleShape,new G4RotationMatrix,-HoleShift);


    //Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(m_VisPCB);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_ClusterX_PCB", m_ClusterX1Detector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * ClusterX_Wafer_Length, 0.5 * ClusterX_Wafer_Width, 0.5 * ClusterX_Wafer_Thickness);

    G4LogicalVolume* logicWafer1 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer1->SetVisAttributes(m_VisSi);
//    logicWafer1->SetSensitiveDetector(m_OuterScorer1);
    G4LogicalVolume* logicWafer2 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer2->SetVisAttributes(m_VisSi);
//    logicWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),-HoleShift, logicWafer1, "TOGAXSI_SI_ClusterX_Wafer 1", m_ClusterX1Detector, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),HoleShift, logicWafer2, "TOGAXSI_SI_ClusterX_Wafer 2", m_ClusterX1Detector, false, 1);


   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("ClusterXActiveWaferShape",0.5 * ClusterX_Wafer_Length, 0.5 * ClusterX_Wafer_Width, 0.5 * ClusterX_Wafer_Thickness);

    G4LogicalVolume* logicActiveWafer1 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(m_VisSi);
    logicActiveWafer1->SetSensitiveDetector(m_ClusterX1Scorer1);
    G4LogicalVolume* logicActiveWafer2 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(m_VisSi);
    logicActiveWafer2->SetSensitiveDetector(m_ClusterX1Scorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer1, "TOGAXSI_SI_ClusterX_Wafer", logicWafer1, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer2, "TOGAXSI_SI_ClusterX_Wafer", logicWafer2, false, 1);


  }
  return m_ClusterX1Detector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildClusterY1Detector(){
  if(!m_ClusterY1Detector){


    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * ClusterY_PCB_Width, 0.5 * ClusterY_PCB_Length, 0.5 * ClusterY_PCB_Thickness);
  
  
    // Master volume inner Detector
    m_ClusterY1Detector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_ClusterY1Detector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", ClusterY_Wafer_Length * 0.5, ClusterY_Wafer_Width * 0.5, ClusterY_PCB_Thickness * 0.5 + 0.1 * mm);

    G4ThreeVector HoleShift = G4ThreeVector(0,0.5 * ClusterY_Wafer_Width + 0.5 * ClusterY_PCB_gap,0);

    G4SubtractionSolid* PCB_temp = new G4SubtractionSolid("PCB", PCBFull, HoleShape, new G4RotationMatrix,HoleShift);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_temp, HoleShape,new G4RotationMatrix,-HoleShift);


    //Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(m_VisPCB);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_ClusterY_PCB", m_ClusterY1Detector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * ClusterY_Wafer_Length, 0.5 * ClusterY_Wafer_Width, 0.5 * ClusterY_Wafer_Thickness);

    G4LogicalVolume* logicWafer1 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer1->SetVisAttributes(m_VisSi);
//    logicWafer1->SetSensitiveDetector(m_OuterScorer1);
    G4LogicalVolume* logicWafer2 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer2->SetVisAttributes(m_VisSi);
//    logicWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),-HoleShift, logicWafer1, "TOGAXSI_SI_ClusterY_Wafer 1", m_ClusterY1Detector, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),HoleShift, logicWafer2, "TOGAXSI_SI_ClusterY_Wafer 2", m_ClusterY1Detector, false, 1);


   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("ClusterYActiveWaferShape",0.5 * ClusterY_Wafer_Length, 0.5 * ClusterY_Wafer_Width, 0.5 * ClusterY_Wafer_Thickness);

    G4LogicalVolume* logicActiveWafer1 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(m_VisSi);
    logicActiveWafer1->SetSensitiveDetector(m_ClusterY1Scorer1);
    G4LogicalVolume* logicActiveWafer2 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(m_VisSi);
    logicActiveWafer2->SetSensitiveDetector(m_ClusterY1Scorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer1, "TOGAXSI_SI_ClusterY_Wafer", logicWafer1, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer2, "TOGAXSI_SI_ClusterY_Wafer", logicWafer2, false, 1);


  }
  return m_ClusterY1Detector;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildClusterX2Detector(){
  if(!m_ClusterX2Detector){


    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * ClusterX_PCB_Width, 0.5 * ClusterX_PCB_Length, 0.5 * ClusterX_PCB_Thickness);
  
  
    // Master volume inner Detector
    m_ClusterX2Detector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_ClusterX2Detector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", ClusterX_Wafer_Length * 0.5, ClusterX_Wafer_Width * 0.5, ClusterX_PCB_Thickness * 0.5 + 0.1 * mm);

    G4ThreeVector HoleShift = G4ThreeVector(0,0.5 * ClusterX_Wafer_Width + 0.5 * ClusterX_PCB_gap,0);

    G4SubtractionSolid* PCB_temp = new G4SubtractionSolid("PCB", PCBFull, HoleShape, new G4RotationMatrix,HoleShift);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_temp, HoleShape,new G4RotationMatrix,-HoleShift);


    //Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(m_VisPCB);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_ClusterX_PCB", m_ClusterX2Detector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * ClusterX_Wafer_Length, 0.5 * ClusterX_Wafer_Width, 0.5 * ClusterX_Wafer_Thickness);

    G4LogicalVolume* logicWafer1 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer1->SetVisAttributes(m_VisSi);
//    logicWafer1->SetSensitiveDetector(m_OuterScorer1);
    G4LogicalVolume* logicWafer2 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer2->SetVisAttributes(m_VisSi);
//    logicWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),-HoleShift, logicWafer1, "TOGAXSI_SI_ClusterX_Wafer 1", m_ClusterX2Detector, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),HoleShift, logicWafer2, "TOGAXSI_SI_ClusterX_Wafer 2", m_ClusterX2Detector, false, 1);


   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("ClusterXActiveWaferShape",0.5 * ClusterX_Wafer_Length, 0.5 * ClusterX_Wafer_Width, 0.5 * ClusterX_Wafer_Thickness);

    G4LogicalVolume* logicActiveWafer1 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(m_VisSi);
    logicActiveWafer1->SetSensitiveDetector(m_ClusterX2Scorer1);
    G4LogicalVolume* logicActiveWafer2 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(m_VisSi);
    logicActiveWafer2->SetSensitiveDetector(m_ClusterX2Scorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer1, "TOGAXSI_SI_ClusterX_Wafer", logicWafer1, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer2, "TOGAXSI_SI_ClusterX_Wafer", logicWafer2, false, 1);


  }
  return m_ClusterX2Detector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildClusterY2Detector(){
  if(!m_ClusterY2Detector){


    G4Box* PCBFull = new G4Box("TOGAXSI_SI_PCBFull", 0.5 * ClusterY_PCB_Width, 0.5 * ClusterY_PCB_Length, 0.5 * ClusterY_PCB_Thickness);
  
  
    // Master volume inner Detector
    m_ClusterY2Detector = new G4LogicalVolume(PCBFull,m_MaterialVacuum,"logicBoxDetector",0,0,0);
    m_ClusterY2Detector->SetVisAttributes(G4VisAttributes::GetInvisible());


    // PCB frame
    G4Box* HoleShape = new G4Box("HoleShape", ClusterY_Wafer_Length * 0.5, ClusterY_Wafer_Width * 0.5, ClusterY_PCB_Thickness * 0.5 + 0.1 * mm);

    G4ThreeVector HoleShift = G4ThreeVector(0,0.5 * ClusterY_Wafer_Width + 0.5 * ClusterY_PCB_gap,0);

    G4SubtractionSolid* PCB_temp = new G4SubtractionSolid("PCB", PCBFull, HoleShape, new G4RotationMatrix,HoleShift);
    G4SubtractionSolid* PCB = new G4SubtractionSolid("PCB", PCB_temp, HoleShape,new G4RotationMatrix,-HoleShift);


    //Sub Volume PCB
    G4LogicalVolume* logicPCB = new G4LogicalVolume(PCB, m_MaterialPCB, "logicPCB", 0, 0, 0);
    logicPCB->SetVisAttributes(m_VisPCB);
   
    new G4PVPlacement(new G4RotationMatrix(0, 0, 0), G4ThreeVector(0, 0, 0), logicPCB, "TOGAXSI_ClusterY_PCB", m_ClusterY2Detector, false, 0);


   /////////////////////////////////////////////////
   // Si Wafer 
   // Subvolume Wafer
    G4Box* WaferShape = new G4Box("WaferShape",0.5 * ClusterY_Wafer_Length, 0.5 * ClusterY_Wafer_Width, 0.5 * ClusterY_Wafer_Thickness);

    G4LogicalVolume* logicWafer1 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer1->SetVisAttributes(m_VisSi);
//    logicWafer1->SetSensitiveDetector(m_OuterScorer1);
    G4LogicalVolume* logicWafer2 = new G4LogicalVolume(WaferShape, m_MaterialSi, "logicWafer", 0, 0, 0);
    logicWafer2->SetVisAttributes(m_VisSi);
//    logicWafer2->SetSensitiveDetector(m_OuterScorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),-HoleShift, logicWafer1, "TOGAXSI_SI_ClusterY_Wafer 1", m_ClusterY2Detector, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),HoleShift, logicWafer2, "TOGAXSI_SI_ClusterY_Wafer 2", m_ClusterY2Detector, false, 1);


   // Active Wafer
    G4Box* ActiveWaferShape = new G4Box("ClusterYActiveWaferShape",0.5 * ClusterY_Wafer_Length, 0.5 * ClusterY_Wafer_Width, 0.5 * ClusterY_Wafer_Thickness);

    G4LogicalVolume* logicActiveWafer1 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer1", 0, 0, 0);
    logicActiveWafer1->SetVisAttributes(m_VisSi);
    logicActiveWafer1->SetSensitiveDetector(m_ClusterY2Scorer1);
    G4LogicalVolume* logicActiveWafer2 = new G4LogicalVolume(ActiveWaferShape, m_MaterialSi, "logicActiveWafer2", 0, 0, 0);
    logicActiveWafer2->SetVisAttributes(m_VisSi);
    logicActiveWafer2->SetSensitiveDetector(m_ClusterY2Scorer2);

    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer1, "TOGAXSI_SI_ClusterY_Wafer", logicWafer1, false, 1);
    new G4PVPlacement(new G4RotationMatrix(0,0,0),G4ThreeVector(0,0,0), logicActiveWafer2, "TOGAXSI_SI_ClusterY_Wafer", logicWafer2, false, 1);


  }
  return m_ClusterY2Detector;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildTarget(int i) {

  if( i > 0) {
    cout << "ERROR: Multiple TOGAXSI target block defined in detector file" << endl;
  }

  G4Material* TargetMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_Target_MaterialName[i]);
  G4Tubs* solidTarget = new G4Tubs("Target", 0., m_Target_R[i], m_Target_L[i] / 2., 0, 360.);
  m_Target = new G4LogicalVolume(solidTarget, TargetMaterial, "Target");
  m_Target->SetVisAttributes(m_VisTarget);

  return m_Target;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_SI::BuildTargetCell(int i) {
  if (i>0) {
    cout << "ERROR: Multiple TOGAXSI target blocks defined in detector file" << endl;
  }

  G4Material* CellMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_Target_CellMaterialName[i]);
  G4Tubs* CellFull = new G4Tubs("CellFull", 0., m_Target_R[i] + m_Target_CellThickness[i], m_Target_L[i] / 2 + m_Target_CellThickness[i], 0, 360.);

  G4Tubs* CellInterior = new G4Tubs("CellInterior", 0., m_Target_R[i], m_Target_L[i] / 2., 0, 360.);

  G4SubtractionSolid* solidCell = new G4SubtractionSolid("solidCell", CellFull, CellInterior, new G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0));

  m_TargetCell = new G4LogicalVolume(solidCell, CellMaterial, "TargetCell");
  m_TargetCell->SetVisAttributes(m_VisTargetCell);
  m_TargetCell->SetSensitiveDetector(m_TargetCellScorer);


  return m_TargetCell;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void TOGAXSI_SI::ReadConfiguration(NPL::InputParser parser){

  //InnerX Si Tracker
  vector<NPL::InputBlock*> blocks_innerX = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","InnerX");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_innerX.size() << " innerX detectors found " << endl; 

  vector<string> wafer_innerX = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_innerX.size() ; i++){
    if(blocks_innerX[i]->HasTokenList(wafer_innerX)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_innerX[i]->GetDouble("Radius","mm");
      double Z = blocks_innerX[i]->GetDouble("Z","mm");
      double Phi = blocks_innerX[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_innerX[i]->GetTVector3("Ref","mm"));
      AddInnerXDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  //InnerZ Si Tracker
  vector<NPL::InputBlock*> blocks_innerZ = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","InnerZ");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_innerZ.size() << " innerZ detectors found " << endl; 

  vector<string> wafer_innerZ = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_innerZ.size() ; i++){
    if(blocks_innerZ[i]->HasTokenList(wafer_innerZ)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_innerZ[i]->GetDouble("Radius","mm");
      double Z = blocks_innerZ[i]->GetDouble("Z","mm");
      double Phi = blocks_innerZ[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_innerZ[i]->GetTVector3("Ref","mm"));
      AddInnerZDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  //OuterX Si Tracker
  vector<NPL::InputBlock*> blocks_outerX = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","OuterX");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_outerX.size() << " outerX detectors found " << endl; 

  vector<string> wafer_outerX = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_outerX.size() ; i++){
    if(blocks_outerX[i]->HasTokenList(wafer_outerX)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_outerX[i]->GetDouble("Radius","mm");
      double Z = blocks_outerX[i]->GetDouble("Z","mm");
      double Phi = blocks_outerX[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_outerX[i]->GetTVector3("Ref","mm"));
      AddOuterXDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  //OuterZ Si Tracker
  vector<NPL::InputBlock*> blocks_outerZ = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","OuterZ");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_outerZ.size() << " outerZ detectors found " << endl; 

  vector<string> wafer_outerZ = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_outerZ.size() ; i++){
    if(blocks_outerZ[i]->HasTokenList(wafer_outerZ)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_outerZ[i]->GetDouble("Radius","mm");
      double Z = blocks_outerZ[i]->GetDouble("Z","mm");
      double Phi = blocks_outerZ[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_outerZ[i]->GetTVector3("Ref","mm"));
      AddOuterZDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }


  //ClusterInner Si Tracker
  vector<NPL::InputBlock*> blocks_clusterInner = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterInner");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterInner.size() << " cluster detectors found " << endl; 

  vector<string> wafer_clusterInner = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_clusterInner.size() ; i++){
    if(blocks_clusterInner[i]->HasTokenList(wafer_clusterInner)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_clusterInner[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterInner[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterInner[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_clusterInner[i]->GetTVector3("Ref","mm"));
      AddClusterInnerDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }


  //ClusterX Si Tracker
  vector<NPL::InputBlock*> blocks_clusterX1 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterX1");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterX1.size() << " ClusterX1 detectors found " << endl; 

  vector<string> wafer_cluster = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_clusterX1.size() ; i++){
    if(blocks_clusterX1[i]->HasTokenList(wafer_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_clusterX1[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterX1[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterX1[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_clusterX1[i]->GetTVector3("Ref","mm"));
      AddClusterX1Detector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  //ClusterY Si Tracker
  vector<NPL::InputBlock*> blocks_clusterY1 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterY1");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterY1.size() << " ClusterY1 detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_clusterY1.size() ; i++){
    if(blocks_clusterY1[i]->HasTokenList(wafer_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_clusterY1[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterY1[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterY1[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_clusterY1[i]->GetTVector3("Ref","mm"));
      AddClusterY1Detector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }


  //ClusterX2 Si Tracker
  vector<NPL::InputBlock*> blocks_clusterX2 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterX2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterX2.size() << " ClusterX2 detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_clusterX2.size() ; i++){
    if(blocks_clusterX2[i]->HasTokenList(wafer_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_clusterX2[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterX2[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterX2[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_clusterX2[i]->GetTVector3("Ref","mm"));
      AddClusterX2Detector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  //ClusterY2 Si Tracker
  vector<NPL::InputBlock*> blocks_clusterY2 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterY2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterY2.size() << " ClusterY2 detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_clusterY2.size() ; i++){
    if(blocks_clusterY2[i]->HasTokenList(wafer_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_clusterY2[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterY2[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterY2[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_clusterY2[i]->GetTVector3("Ref","mm"));
      AddClusterY2Detector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }



  //Target
  vector<NPL::InputBlock*> blocks_target = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","Target");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_target.size() << " target found " << endl; 

  vector<string> targettoken = {"Radius","Length","TargetMaterial","CellMaterial","CellThickness","Pos"};

  for(unsigned int i = 0 ; i < blocks_target.size() ; i++){
    if(blocks_target[i]->HasTokenList(targettoken)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI " << i+1 <<  endl;
      double R = blocks_target[i]->GetDouble("Radius","mm");
      double L = blocks_target[i]->GetDouble("Length","mm");
      string TargetMaterialName = blocks_target[i]->GetString("TargetMaterial");
      string CellMaterialName = blocks_target[i]->GetString("CellMaterial");
      double CellThickness = blocks_target[i]->GetDouble("CellThickness","mm");
      G4ThreeVector Pos = NPS::ConvertVector(blocks_target[i]->GetTVector3("Pos","mm"));  
      AddTarget(R,L,TargetMaterialName,CellMaterialName,CellThickness,Pos);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void TOGAXSI_SI::ConstructDetector(G4LogicalVolume* world){

  //InnerX Si tracker
  cout << "______________Size: " << m_InnerX_R.size() << endl;
  for (unsigned short i = 0; i < m_InnerX_R.size(); i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_InnerX_R[i], m_InnerX_Z[i]);
   
    Det_pos.rotate(-m_InnerX_Phi[i], G4ThreeVector(0,0,1));
//    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
   G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_InnerX_Phi[i]);
//    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    cout << m_InnerX_Phi[i] << endl; 
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_InnerX_Ref[i]), BuildInnerXDetector(),"TOGAXSI_SI_InnerX",world, false, i + 1);
 
  }	  

  //InnerZ Si tracker
  cout << "______________Size: " << m_InnerZ_R.size() << endl;
  for (unsigned short i = 0; i < m_InnerZ_R.size(); i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_InnerZ_R[i], m_InnerZ_Z[i]);
   
    Det_pos.rotate(-m_InnerZ_Phi[i], G4ThreeVector(0,0,1));
//    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
   G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_InnerZ_Phi[i]);
//    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    cout << m_InnerZ_Phi[i] << endl; 
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_InnerZ_Ref[i]), BuildInnerZDetector(),"TOGAXSI_SI_InnerZ",world, false, i + 1);
 
  }	  

  //OuterX Si tracker
  for (unsigned short i = 0; i < m_OuterX_R.size(); i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_OuterX_R[i], m_OuterX_Z[i]);
   
    Det_pos.rotate(-m_OuterX_Phi[i], G4ThreeVector(0,0,1));
//    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
   G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_OuterX_Phi[i]);
//    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    cout << m_OuterX_Phi[i] << endl; 
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_OuterX_Ref[i]), BuildOuterXDetector(),"TOGAXSI_SI_OuterX",world, false, i + 1);
 
  }	  

  //OuterZ Si tracker
  for (unsigned short i = 0; i < m_OuterZ_R.size(); i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_OuterZ_R[i], m_OuterZ_Z[i]);
   
    Det_pos.rotate(-m_OuterZ_Phi[i], G4ThreeVector(0,0,1));
//    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
   G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_OuterZ_Phi[i]);
//    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    cout << m_OuterZ_Phi[i] << endl; 
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_OuterZ_Ref[i]), BuildOuterZDetector(),"TOGAXSI_SI_OuterZ",world, false, i + 1);
 
  }	  


  //Cluster Si tracker
  for (unsigned short i = 0; i < m_ClusterInner_R.size(); i++) {

    G4ThreeVector Det_pos = G4ThreeVector(0,m_ClusterInner_R[i], m_ClusterInner_Z[i]);
   
    Det_pos.rotate(-m_ClusterInner_Phi[i], G4ThreeVector(0,0,1));
//    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
   G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_ClusterInner_Phi[i]);
//    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    cout << m_ClusterInner_Phi[i] << endl; 
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_ClusterInner_Ref[i]), BuildClusterInnerDetector(),"TOGAXSI_SI_ClusterInner",world, false, i + 1);
 
  }	  

  //ClusterX1 Si tracker
  for (unsigned short i = 0; i < m_ClusterX1_R.size(); i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_ClusterX1_R[i],0, m_ClusterX1_Z[i]);
   
    Det_pos.rotate(-m_ClusterX1_Phi[i], G4ThreeVector(0,0,1));

    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_ClusterX1_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_ClusterX1_Ref[i]), BuildClusterX1Detector(),"TOGAXSI_SI_ClusterX1",world, false, i + 1); 
  }	  

  //ClusterY1 Si tracker
  for (unsigned short i = 0; i < m_ClusterY1_R.size(); i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_ClusterY1_R[i],0, m_ClusterY1_Z[i]);
   
    Det_pos.rotate(-m_ClusterY1_Phi[i], G4ThreeVector(0,0,1));

    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_ClusterY1_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_ClusterY1_Ref[i]), BuildClusterY1Detector(),"TOGAXSI_SI_ClusterY1",world, false, i + 1); 

  }	  

  
  //ClusterX1 Si tracker
  for (unsigned short i = 0; i < m_ClusterX2_R.size(); i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_ClusterX2_R[i],0, m_ClusterX2_Z[i]);
   
    Det_pos.rotate(-m_ClusterX2_Phi[i], G4ThreeVector(0,0,1));

    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_ClusterX2_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_ClusterX2_Ref[i]), BuildClusterX2Detector(),"TOGAXSI_SI_ClusterX2",world, false, i + 1); 
	   
  }	  

  //ClusterY2 Si tracker
  for (unsigned short i = 0; i < m_ClusterY2_R.size(); i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_ClusterY2_R[i],0, m_ClusterY2_Z[i]);
   
    Det_pos.rotate(-m_ClusterY2_Phi[i], G4ThreeVector(0,0,1));

    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_ClusterY2_Phi[i]);

    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_ClusterY2_Ref[i]), BuildClusterY2Detector(),"TOGAXSI_SI_ClusterY2",world, false, i + 1); 
 
  }	  


  // Target
  G4LogicalVolume* logicTarget[m_Target_R.size()];
  G4LogicalVolume* logicTargetCell[m_Target_R.size()];

  for (unsigned short i = 0; i < m_Target_R.size(); i++) {
    G4ThreeVector Tar_pos = m_Target_Pos[i];
    cout << "TargetPos" << m_Target_Pos[i].z() << endl;
    G4RotationMatrix* Rot = new G4RotationMatrix();
    logicTarget[i] = BuildTarget(i);
    new G4PVPlacement(Rot, Tar_pos, logicTarget[i], "TOGAXSI_SI_Target", world, false, i + 1);
    logicTargetCell[i] = BuildTargetCell(i);
    new G4PVPlacement(Rot, Tar_pos, logicTargetCell[i], "TOGAXSI_SI_TargetCell", world, false, i + 1);

    if (!m_ReactionRegion) {
      m_ReactionRegion = new G4Region("NPSimulationProcess");
    }

    m_ReactionRegion->AddRootLogicalVolume(m_Target);
    m_ReactionRegion->SetUserLimits(new G4UserLimits(1. * mm));

    G4FastSimulationManager* mng = m_ReactionRegion->GetFastSimulationManager();
    unsigned int size = m_ReactionModel.size();
    for (unsigned int o = 0; o < size; o++) {
      mng->RemoveFastSimulationModel(m_ReactionModel[o]);
    }
    m_ReactionModel.clear();
    
    G4VFastSimulationModel* fsm;
    fsm = new NPS::BeamReaction("BeamReaction", m_ReactionRegion);
    ((NPS::BeamReaction*)fsm)->SetStepSize(1. * mm);
    m_ReactionModel.push_back(fsm);

    fsm = new NPS::Decay("Decay", m_ReactionRegion);
    m_ReactionModel.push_back(fsm);

  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void TOGAXSI_SI::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("TOGAXSI_SI")){
    pTree->Branch("TOGAXSI_SI", "TTOGAXSI_SIData", &m_Event) ;
  }
  pTree->SetBranchAddress("TOGAXSI_SI", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void TOGAXSI_SI::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  //InnerX Si tracker scorer
  DSSDScorers::PS_Rectangle* InnerXScorer = (DSSDScorers::PS_Rectangle*)m_InnerXScorer->GetPrimitive(0);

  unsigned int size = InnerXScorer->GetLengthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(InnerXScorer->GetEnergyLength(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = InnerXScorer->GetDetectorLength(i);
      int StripLongitudinal = InnerXScorer->GetStripLength(i);
      m_Event->SetInnerXE(DetNbr, StripLongitudinal, Energy);
    }
  }
  InnerXScorer->clear();
  
  ///////////
  //InnerZ Si tracker scorer
  DSSDScorers::PS_Rectangle* InnerZScorer = (DSSDScorers::PS_Rectangle*)m_InnerZScorer->GetPrimitive(0);

  size = InnerZScorer->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(InnerZScorer->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = InnerZScorer->GetDetectorWidth(i);
      int StripTransverse = InnerZScorer->GetStripWidth(i);
      m_Event->SetInnerZE(DetNbr, StripTransverse, Energy);
    }
  }

  InnerZScorer->clear();
  

  ///////////
  //OuterX first Si tracker scorer
  DSSDScorers::PS_Rectangle* OuterXScorer1 = (DSSDScorers::PS_Rectangle*)m_OuterXScorer1->GetPrimitive(0);

  size = OuterXScorer1->GetLengthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(OuterXScorer1->GetEnergyLength(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = OuterXScorer1->GetDetectorLength(i);
      int StripLongitudinal = OuterXScorer1->GetStripLength(i);
      m_Event->SetOuterXE(DetNbr, StripLongitudinal, Energy);
    }
  }

  OuterXScorer1->clear();

  ///////////
  //OuterX second Si tracker scorer
  DSSDScorers::PS_Rectangle* OuterXScorer2 = (DSSDScorers::PS_Rectangle*)m_OuterXScorer2->GetPrimitive(0);

  size = OuterXScorer2->GetLengthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(OuterXScorer2->GetEnergyLength(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = OuterXScorer2->GetDetectorLength(i);
      int StripLongitudinal = OuterXScorer2->GetStripLength(i);
      m_Event->SetOuterXE(DetNbr, StripLongitudinal, Energy);
    }
  }

  OuterXScorer2->clear();


  ///////////
  //OuterZ first Si tracker scorer
  DSSDScorers::PS_Rectangle* OuterZScorer1 = (DSSDScorers::PS_Rectangle*)m_OuterZScorer1->GetPrimitive(0);

  size = OuterZScorer1->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(OuterZScorer1->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = OuterZScorer1->GetDetectorWidth(i);
      int StripTransverse = OuterZScorer1->GetStripWidth(i);
      m_Event->SetOuterZE(DetNbr, StripTransverse, Energy);
    }
  }

  OuterZScorer1->clear();

  ///////////
  //OuterZ second Si tracker scorer
  DSSDScorers::PS_Rectangle* OuterZScorer2 = (DSSDScorers::PS_Rectangle*)m_OuterZScorer2->GetPrimitive(0);

  size = OuterZScorer2->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(OuterZScorer2->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = OuterZScorer2->GetDetectorWidth(i);
      int StripTransverse = OuterZScorer2->GetStripWidth(i) + OuterZ_Wafer_TransverseStrips;
      m_Event->SetOuterZE(DetNbr, StripTransverse, Energy);
    }
  }

  OuterZScorer2->clear();


  ///////////
  //ClusterX1 first Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterX1Scorer1 = (DSSDScorers::PS_Rectangle*)m_ClusterX1Scorer1->GetPrimitive(0);

  size = ClusterX1Scorer1->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterX1Scorer1->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterX1Scorer1->GetDetectorWidth(i);
      int StripLongitudinal = ClusterX1Scorer1->GetStripLength(i);
      m_Event->SetClusterX1E(DetNbr, StripLongitudinal, Energy);
    }
  }

  ClusterX1Scorer1->clear();

  ///////////
  //ClusterX1 second Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterX1Scorer2 = (DSSDScorers::PS_Rectangle*)m_ClusterX1Scorer2->GetPrimitive(0);

  size = ClusterX1Scorer2->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterX1Scorer2->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterX1Scorer2->GetDetectorWidth(i);
      int StripLongitudinal = ClusterX1Scorer2->GetStripLength(i);
      m_Event->SetClusterX1E(DetNbr, StripLongitudinal, Energy);
    }
  }

  ClusterX1Scorer2->clear();


  ///////////
  //ClusterY1 first Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterY1Scorer1 = (DSSDScorers::PS_Rectangle*)m_ClusterY1Scorer1->GetPrimitive(0);

  size = ClusterY1Scorer1->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterY1Scorer1->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterY1Scorer1->GetDetectorWidth(i);
      int StripTransverse = ClusterY1Scorer1->GetStripWidth(i);
      m_Event->SetClusterY1E(DetNbr, StripTransverse, Energy);
    }
  }

  ClusterY1Scorer1->clear();

  ///////////
  //ClusterY1 second Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterY1Scorer2 = (DSSDScorers::PS_Rectangle*)m_ClusterY1Scorer2->GetPrimitive(0);

  size = ClusterY1Scorer2->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterY1Scorer2->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterY1Scorer2->GetDetectorWidth(i);
      int StripTransverse = ClusterY1Scorer2->GetStripWidth(i) + ClusterY_Wafer_TransverseStrips;
      m_Event->SetClusterY1E(DetNbr, StripTransverse, Energy);
    }
  }

  ClusterY1Scorer2->clear();



  ///////////
  //ClusterX2 first Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterX2Scorer1 = (DSSDScorers::PS_Rectangle*)m_ClusterX2Scorer1->GetPrimitive(0);

  size = ClusterX2Scorer1->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterX2Scorer1->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterX2Scorer1->GetDetectorWidth(i);
      int StripLongitudinal = ClusterX2Scorer1->GetStripLength(i);
      m_Event->SetClusterX2E(DetNbr, StripLongitudinal, Energy);
    }
  }

  ClusterX2Scorer1->clear();

  ///////////
  //ClusterX1 second Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterX2Scorer2 = (DSSDScorers::PS_Rectangle*)m_ClusterX2Scorer2->GetPrimitive(0);

  size = ClusterX2Scorer2->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterX2Scorer2->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterX2Scorer2->GetDetectorWidth(i);
      int StripLongitudinal = ClusterX2Scorer2->GetStripLength(i);
      m_Event->SetClusterX2E(DetNbr, StripLongitudinal, Energy);
    }
  }

  ClusterX2Scorer2->clear();


  ///////////
  //ClusterY1 first Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterY2Scorer1 = (DSSDScorers::PS_Rectangle*)m_ClusterY2Scorer1->GetPrimitive(0);

  size = ClusterY2Scorer1->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterY2Scorer1->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterY2Scorer1->GetDetectorWidth(i);
      int StripTransverse = ClusterY2Scorer1->GetStripWidth(i);
      m_Event->SetClusterY2E(DetNbr, StripTransverse, Energy);
    }
  }

  ClusterY2Scorer1->clear();

  ///////////
  //ClusterY2 second Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterY2Scorer2 = (DSSDScorers::PS_Rectangle*)m_ClusterY2Scorer2->GetPrimitive(0);

  size = ClusterY2Scorer2->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterY2Scorer2->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterY2Scorer2->GetDetectorWidth(i);
      int StripTransverse = ClusterY2Scorer2->GetStripWidth(i) + ClusterY_Wafer_TransverseStrips;
      m_Event->SetClusterY2E(DetNbr, StripTransverse, Energy);
    }
  }

  ClusterY2Scorer2->clear();


  ///////////
  //ClusterInner Si tracker scorer
  DSSDScorers::PS_Rectangle* ClusterInnerScorer = (DSSDScorers::PS_Rectangle*)m_ClusterInnerScorer->GetPrimitive(0);

  size = ClusterInnerScorer->GetWidthMult();
  for (unsigned int i = 0; i < size; i++) {
    double Energy = RandGauss::shoot(ClusterInnerScorer->GetEnergyWidth(i), ResoEnergy);
    if (Energy > EnergyThreshold) {
      int DetNbr = ClusterInnerScorer->GetDetectorWidth(i);
      int StripTransverse = ClusterInnerScorer->GetStripWidth(i);
      m_Event->SetClusterInnerE(DetNbr, StripTransverse, Energy);
    }
  }
  ClusterInnerScorer->clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void TOGAXSI_SI::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_InnerXScorer = CheckScorer("InnerXScorer", already_exist);
  m_InnerZScorer = CheckScorer("InnerZScorer", already_exist);

  m_OuterXScorer1 = CheckScorer("OuterXScorer1", already_exist);
  m_OuterXScorer2 = CheckScorer("OuterXScorer2", already_exist);

  m_OuterZScorer1 = CheckScorer("OuterZScorer1", already_exist);
  m_OuterZScorer2 = CheckScorer("OuterZScorer2", already_exist);

  m_ClusterInnerScorer = CheckScorer("ClusterInnerScorer", already_exist);

  m_ClusterX1Scorer1 = CheckScorer("ClusterX1Scorer1", already_exist);
  m_ClusterX1Scorer2 = CheckScorer("ClusterX1Scorer2", already_exist);
  m_ClusterY1Scorer1 = CheckScorer("ClusterY1Scorer1", already_exist);
  m_ClusterY1Scorer2 = CheckScorer("ClusterY1Scorer2", already_exist);
  m_ClusterX2Scorer1 = CheckScorer("ClusterX2Scorer1", already_exist);
  m_ClusterX2Scorer2 = CheckScorer("ClusterX2Scorer2", already_exist);
  m_ClusterY2Scorer1 = CheckScorer("ClusterY2Scorer1", already_exist);
  m_ClusterY2Scorer2 = CheckScorer("ClusterY2Scorer2", already_exist);

  m_TargetCellScorer = CheckScorer("TargetCellScorer", already_exist);

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  //Si tracker InnerX
  G4VPrimitiveScorer* InnerXScorer = new DSSDScorers::PS_Rectangle("InnerXScorer", 2, InnerX_Wafer_Width, InnerX_Wafer_Length, InnerX_Wafer_LongitudinalStrips, InnerX_Wafer_TransverseStrips, 0, "xz");

  G4VPrimitiveScorer* InteractionInnerX = new InteractionScorers::PS_Interactions("InteractionInnerX", ms_InterCoord, 0);

  //Si tracker InnerZ
  G4VPrimitiveScorer* InnerZScorer = new DSSDScorers::PS_Rectangle("InnerZScorer", 2, InnerZ_Wafer_Width, InnerZ_Wafer_Length, InnerZ_Wafer_LongitudinalStrips, InnerZ_Wafer_TransverseStrips, 0, "xz");

  G4VPrimitiveScorer* InteractionInnerZ = new InteractionScorers::PS_Interactions("InteractionInnerZ", ms_InterCoord, 0);


  //Si tracker OuterX
  G4VPrimitiveScorer* OuterXScorer1 = new DSSDScorers::PS_Rectangle("OuterXScorer1", 2, OuterX_Wafer_Length, OuterX_Wafer_Width, OuterX_Wafer_LongitudinalStrips, OuterX_Wafer_TransverseStrips, 0, "xz");
  G4VPrimitiveScorer* OuterXScorer2 = new DSSDScorers::PS_Rectangle("OuterXScorer2", 2, OuterX_Wafer_Length, OuterX_Wafer_Width, OuterX_Wafer_LongitudinalStrips, OuterX_Wafer_TransverseStrips, 0, "xz");

//  G4VPrimitiveScorer* OuterXScorer1 = new DSSDScorers::PS_Rectangle("OuterXScorer1", 2, OuterX_Wafer_Width, OuterX_Wafer_Length, OuterX_Wafer_LongitudinalStrips, OuterX_Wafer_TransverseStrips, 0, "xz");
//  G4VPrimitiveScorer* OuterXScorer2 = new DSSDScorers::PS_Rectangle("OuterXScorer2", 2, OuterX_Wafer_Width, OuterX_Wafer_Length, OuterX_Wafer_LongitudinalStrips, OuterX_Wafer_TransverseStrips, 0, "xz");

  G4VPrimitiveScorer* InteractionOuterX1 = new InteractionScorers::PS_Interactions("InteractionOuterX1", ms_InterCoord, 0);
  G4VPrimitiveScorer* InteractionOuterX2 = new InteractionScorers::PS_Interactions("InteractionOuterX2", ms_InterCoord, 0);

  //Si tracker OuterZ
  G4VPrimitiveScorer* OuterZScorer1 = new DSSDScorers::PS_Rectangle("OuterZScorer1", 2, OuterZ_Wafer_Length, OuterZ_Wafer_Width, OuterZ_Wafer_LongitudinalStrips, OuterZ_Wafer_TransverseStrips, 0, "xz");
  G4VPrimitiveScorer* OuterZScorer2 = new DSSDScorers::PS_Rectangle("OuterZScorer2", 2, OuterZ_Wafer_Length, OuterZ_Wafer_Width, OuterZ_Wafer_LongitudinalStrips, OuterZ_Wafer_TransverseStrips, 0, "xz");

  G4VPrimitiveScorer* InteractionOuterZ1 = new InteractionScorers::PS_Interactions("InteractionOuterZ1", ms_InterCoord, 0);
  G4VPrimitiveScorer* InteractionOuterZ2 = new InteractionScorers::PS_Interactions("InteractionOuterZ2", ms_InterCoord, 0);


  //Si tracker Cluster
  G4VPrimitiveScorer* ClusterInnerScorer = new DSSDScorers::PS_Rectangle("ClusterInnerScorer", 2, ClusterInner_Wafer_Height, ClusterInner_Wafer_Base, ClusterInner_Wafer_LongitudinalStrips, ClusterInner_Wafer_TransverseStrips, 0, "xz");

  G4VPrimitiveScorer* InteractionClusterInner = new InteractionScorers::PS_Interactions("InteractionClusterInner", ms_InterCoord, 0);


  //Si tracker Cluster1 demonstrator
  G4VPrimitiveScorer* ClusterX1Scorer1 = new DSSDScorers::PS_Rectangle("ClusterX1Scorer1", 2, ClusterX_Wafer_Length, ClusterX_Wafer_Width, ClusterX_Wafer_LongitudinalStrips, ClusterX_Wafer_TransverseStrips, 0, "xy");
  G4VPrimitiveScorer* ClusterX1Scorer2 = new DSSDScorers::PS_Rectangle("ClusterX1Scorer2", 2, ClusterX_Wafer_Length, ClusterX_Wafer_Width, ClusterX_Wafer_LongitudinalStrips, ClusterX_Wafer_TransverseStrips, 0, "xy");

  G4VPrimitiveScorer* InteractionClusterX1_1 = new InteractionScorers::PS_Interactions("InteractionClusterX1_1", ms_InterCoord, 0);
  G4VPrimitiveScorer* InteractionClusterX1_2 = new InteractionScorers::PS_Interactions("InteractionClusterX1_2", ms_InterCoord, 0);

  G4VPrimitiveScorer* ClusterY1Scorer1 = new DSSDScorers::PS_Rectangle("ClusterY1Scorer1", 2, ClusterY_Wafer_Length, ClusterY_Wafer_Width, ClusterY_Wafer_LongitudinalStrips, ClusterY_Wafer_TransverseStrips, 0, "xy");
  G4VPrimitiveScorer* ClusterY1Scorer2 = new DSSDScorers::PS_Rectangle("ClusterY1Scorer2", 2, ClusterY_Wafer_Length, ClusterY_Wafer_Width, ClusterY_Wafer_LongitudinalStrips, ClusterY_Wafer_TransverseStrips, 0, "xy");

  G4VPrimitiveScorer* InteractionClusterY1_1 = new InteractionScorers::PS_Interactions("InteractionClusterY1_1", ms_InterCoord, 0);
  G4VPrimitiveScorer* InteractionClusterY1_2 = new InteractionScorers::PS_Interactions("InteractionClusterY1_2", ms_InterCoord, 0);

  //Si tracker Cluster2 demonstrator
  G4VPrimitiveScorer* ClusterX2Scorer1 = new DSSDScorers::PS_Rectangle("ClusterX2Scorer1", 2, ClusterX_Wafer_Length, ClusterX_Wafer_Width, ClusterX_Wafer_LongitudinalStrips, ClusterX_Wafer_TransverseStrips, 0, "xy");
  G4VPrimitiveScorer* ClusterX2Scorer2 = new DSSDScorers::PS_Rectangle("ClusterX2Scorer2", 2, ClusterX_Wafer_Length, ClusterX_Wafer_Width, ClusterX_Wafer_LongitudinalStrips, ClusterX_Wafer_TransverseStrips, 0, "xy");

  G4VPrimitiveScorer* InteractionClusterX2_1 = new InteractionScorers::PS_Interactions("InteractionClusterX2_1", ms_InterCoord, 0);
  G4VPrimitiveScorer* InteractionClusterX2_2 = new InteractionScorers::PS_Interactions("InteractionClusterX2_2", ms_InterCoord, 0);

  G4VPrimitiveScorer* ClusterY2Scorer1 = new DSSDScorers::PS_Rectangle("ClusterY2Scorer1", 2, ClusterY_Wafer_Length, ClusterY_Wafer_Width, ClusterY_Wafer_LongitudinalStrips, ClusterY_Wafer_TransverseStrips, 0, "xy");
  G4VPrimitiveScorer* ClusterY2Scorer2 = new DSSDScorers::PS_Rectangle("ClusterY2Scorer2", 2, ClusterY_Wafer_Length, ClusterY_Wafer_Width, ClusterY_Wafer_LongitudinalStrips, ClusterY_Wafer_TransverseStrips, 0, "xy");

  G4VPrimitiveScorer* InteractionClusterY2_1 = new InteractionScorers::PS_Interactions("InteractionClusterY2_1", ms_InterCoord, 0);
  G4VPrimitiveScorer* InteractionClusterY2_2 = new InteractionScorers::PS_Interactions("InteractionClusterY2_2", ms_InterCoord, 0);

  //TargetCell
//  G4VPrimitiveScorer* TargetCellScorer = new PS_CalorimeterScorers::PS_Calorimeter("ClusterY2Scorer1", 2, ClusterY_Wafer_Length, ClusterY_Wafer_Width, ClusterY_Wafer_LongitudinalStrips, ClusterY_Wafer_TransverseStrips, 0, "xy");

  G4VPrimitiveScorer* InteractionTargetCell = new InteractionScorers::PS_Interactions("InteractionTargetCell", ms_InterCoord, 0);



  //and register it to the multifunctionnal detector
  m_InnerXScorer->RegisterPrimitive(InnerXScorer);
  m_InnerXScorer->RegisterPrimitive(InteractionInnerX);
  m_InnerZScorer->RegisterPrimitive(InnerZScorer);
  m_InnerZScorer->RegisterPrimitive(InteractionInnerZ);

  m_OuterXScorer1->RegisterPrimitive(OuterXScorer1);
  m_OuterXScorer1->RegisterPrimitive(InteractionOuterX1);
  m_OuterXScorer2->RegisterPrimitive(OuterXScorer2);
  m_OuterXScorer2->RegisterPrimitive(InteractionOuterX2);

  m_OuterZScorer1->RegisterPrimitive(OuterZScorer1);
  m_OuterZScorer1->RegisterPrimitive(InteractionOuterZ1);
  m_OuterZScorer2->RegisterPrimitive(OuterZScorer2);
  m_OuterZScorer2->RegisterPrimitive(InteractionOuterZ2);

  m_ClusterInnerScorer->RegisterPrimitive(ClusterInnerScorer);
  m_ClusterInnerScorer->RegisterPrimitive(InteractionClusterInner);

  m_ClusterX1Scorer1->RegisterPrimitive(ClusterX1Scorer1);
  m_ClusterX1Scorer1->RegisterPrimitive(InteractionClusterX1_1);
  m_ClusterX1Scorer2->RegisterPrimitive(ClusterX1Scorer2);
  m_ClusterX1Scorer2->RegisterPrimitive(InteractionClusterX1_2);
  m_ClusterY1Scorer1->RegisterPrimitive(ClusterY1Scorer1);
  m_ClusterY1Scorer1->RegisterPrimitive(InteractionClusterY1_1);
  m_ClusterY1Scorer2->RegisterPrimitive(ClusterY1Scorer2);
  m_ClusterY1Scorer2->RegisterPrimitive(InteractionClusterY1_2);

  m_ClusterX2Scorer1->RegisterPrimitive(ClusterX2Scorer1);
  m_ClusterX2Scorer1->RegisterPrimitive(InteractionClusterX2_1);
  m_ClusterX2Scorer2->RegisterPrimitive(ClusterX2Scorer2);
  m_ClusterX2Scorer2->RegisterPrimitive(InteractionClusterX2_2);
  m_ClusterY2Scorer1->RegisterPrimitive(ClusterY2Scorer1);
  m_ClusterY2Scorer1->RegisterPrimitive(InteractionClusterY2_1);
  m_ClusterY2Scorer2->RegisterPrimitive(ClusterY2Scorer2);
  m_ClusterY2Scorer2->RegisterPrimitive(InteractionClusterY2_2);

//  m_TargetCellScorer->RegisterPrimitive(TargetCellScorer);
  m_TargetCellScorer->RegisterPrimitive(InteractionTargetCell);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerXScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_InnerZScorer);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterXScorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterXScorer2);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterZScorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_OuterZScorer2);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterInnerScorer);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterX1Scorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterX1Scorer2);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterY1Scorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterY1Scorer2);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterX2Scorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterX2Scorer2);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterY2Scorer1);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ClusterY2Scorer2);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_TargetCellScorer);

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* TOGAXSI_SI::Construct(){
  return  (NPS::VDetector*) new TOGAXSI_SI();
}

void TOGAXSI_SI::InitializeMaterial() {
  m_MaterialSi = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");//SI
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");//PCB 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_TOGAXSI_SI{
    public:
      proxy_nps_TOGAXSI_SI(){
        NPS::DetectorFactory::getInstance()->AddToken("TOGAXSI_SI","TOGAXSI_SI");
        NPS::DetectorFactory::getInstance()->AddDetector("TOGAXSI_SI",TOGAXSI_SI::Construct);
      }
  };

  proxy_nps_TOGAXSI_SI p_nps_TOGAXSI_SI;
}
