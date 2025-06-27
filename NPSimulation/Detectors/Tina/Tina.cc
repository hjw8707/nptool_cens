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
#include <sstream>
#include <cmath>
#include <limits>

// Geant4 
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Randomize.hh"

// NPTool 
#include "Tina.hh"
#include "CalorimeterScorers.hh"
#include "SiliconScorers.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
#include "NPCore.h"

// CLHEP 
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace TINA;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Tina::Tina(){
  m_Event = new TTinaData() ;
  m_TTTPad = 0;
  m_YY1CsI = 0;
  m_TTTScorer = 0;
  m_PadScorer = 0;
  m_YY1Scorer = 0;
  m_CsIScorer = 0;

  // set RGB colour and opacity
  /*  m_VisTTT = new G4VisAttributes(G4Colour(0,1,0,0.5));
      m_VisPad = new G4VisAttributes(G4Colour(0,1,0,0.5)); 
      m_VisYY1 = new G4VisAttributes(G4Colour(0,1,0,0.5)); 
      m_VisCsI = new G4VisAttributes(G4Colour(0,1,0,0.5)); //*/
  m_VisTTT = new G4VisAttributes(G4Colour(0., 0.5, 0.5));
  m_VisPad = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
  m_VisYY1 = new G4VisAttributes(G4Colour(0., 0.5, 0.5));
  m_VisCsI = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Tina::~Tina(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::AddDetector(G4ThreeVector POS, double Alpha, string Shape){
  // convert the POS value to r theta phi as spherical coordinates are easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_T.push_back(0.);
  m_Alpha.push_back(Alpha);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::AddDetector(double R, double Theta, double Phi, double T, double Alpha, string Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_T.push_back(T);
  m_Alpha.push_back(Alpha);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/* G4LogicalVolume* Tina::BuildTTTDetector(){
  if(!m_TTT){
    G4Box* solidTTT = new G4Box("Tina_TTT",XwidthTTT/2.,YwidthTTT/2.,ZwidthTTT/2.);
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    m_TTT = new G4LogicalVolume(solidTTT,DetectorMaterial,"logic_Tina_TTT",0,0,0);
    m_TTT->SetVisAttributes(m_VisTTT);
    m_TTT->SetSensitiveDetector(m_TTTScorer);
  }
  return m_TTT;
} */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/* G4LogicalVolume* Tina::BuildPadSiDetector(){
  if(!m_PadSi){
    // G4Box* solidPadSi1 = new G4Box("Tina_PadSi1",Tina_NS::Width*26.,Tina_NS::Width*27.5,Tina_NS::Thickness*0.15);
    // G4Box* solidPadSi2 = new G4Box("Tina_PadSi2",Tina_NS::Width*26.,Tina_NS::Width*27.5,Tina_NS::Thickness*0.15);
    // G4Box* solidPadSi3 = new G4Box("Tina_PadSi3",Tina_NS::Width*26.,Tina_NS::Width*27.5,Tina_NS::Thickness*0.15);
    // G4Box* solidPadSi4 = new G4Box("Tina_PadSi4",Tina_NS::Width*26.,Tina_NS::Width*27.5,Tina_NS::Thickness*0.15);
    // G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    // replace PadSi at first by one, then by four thick CsI blocks enclosed in single box
    // G4Box* solidPadSi0 = new G4Box("Tina_PadSi0",Tina_NS::Width*62.,Tina_NS::Width*57.5,Tina_NS::Thickness*9.75);
    G4Box* solidPadSi1 = new G4Box("Tina_PadSi1",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadSi2 = new G4Box("Tina_PadSi2",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadSi3 = new G4Box("Tina_PadSi3",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadSi4 = new G4Box("Tina_PadSi4",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
    // G4LogicalVolume* logicPadSi0 = new G4LogicalVolume(solidPadSi0,DetectorMaterial,"logic_Tina_PadSi0",0,0,0);
    G4LogicalVolume* logicPadSi1 = new G4LogicalVolume(solidPadSi1,DetectorMaterial,"logic_Tina_PadSi1",0,0,0);
    G4LogicalVolume* logicPadSi2 = new G4LogicalVolume(solidPadSi2,DetectorMaterial,"logic_Tina_PadSi2",0,0,0);
    G4LogicalVolume* logicPadSi3 = new G4LogicalVolume(solidPadSi3,DetectorMaterial,"logic_Tina_PadSi3",0,0,0);
    G4LogicalVolume* logicPadSi4 = new G4LogicalVolume(solidPadSi4,DetectorMaterial,"logic_Tina_PadSi4",0,0,0);

    // G4Box* solidPadSi = new G4Box("Tina_PadSi",Tina_NS::Width*64.,Tina_NS::Width*59.5,Tina_NS::Thickness*2.75);
    G4Box* solidPadSi = new G4Box("Tina_PadSi",XwidthBox/2.,YwidthBox/2.,ZwidthBox/2.);
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    m_PadSi = new G4LogicalVolume(solidPadSi,Vacuum,"logic_Tina_PadSi",0,0,0);

    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    // x = XdistbetwPads/2 + XwidthPad/2 = 29.5
    // y = YdistbetwPads/2 + YwidthPad/2 = 28.5
    // half-thickness of frame is 2.75 and middle of detector sits 1+0.15 deep from the surface therefore z=-1.6
    // new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(29.5,28.5,-1.6)),logicPadSi1,"Tina_PadSi_PadSi1",m_PadSi,false,1,false);
    // new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(-29.5,28.5,-1.6)),logicPadSi2,"Tina_PadSi_PadSi2",m_PadSi,false,2,false);
    // new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(29.5,-28.5,-1.6)),logicPadSi3,"Tina_PadSi_PadSi3",m_PadSi,false,3,false);
    // new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(-29.5,-28.5,-1.6)),logicPadSi4,"Tina_PadSi_PadSi4",m_PadSi,false,4,false);
    // new G4PVPlacement(G4Transform3D(*Rot,G4ThreeVector(0.,0.,0.)),logicPadSi0,"Tina_PadSi_PadSi0",m_PadSi,false,1,false);
    new G4PVPlacement(
         G4Transform3D(*Rot,G4ThreeVector( (XwidthPad+XdistbetwPads)/2., (YwidthPad+YdistbetwPads)/2.,ZoffsetPad)),
         logicPadSi1,"Tina_PadSi_PadSi1",m_PadSi,false,1,false);
    new G4PVPlacement(
         G4Transform3D(*Rot,G4ThreeVector(-(XwidthPad+XdistbetwPads)/2., (YwidthPad+YdistbetwPads)/2.,ZoffsetPad)),
         logicPadSi2,"Tina_PadSi_PadSi2",m_PadSi,false,2,false);
    new G4PVPlacement(
         G4Transform3D(*Rot,G4ThreeVector( (XwidthPad+XdistbetwPads)/2.,-(YwidthPad+YdistbetwPads)/2.,ZoffsetPad)),
         logicPadSi3,"Tina_PadSi_PadSi3",m_PadSi,false,3,false);
    new G4PVPlacement(
         G4Transform3D(*Rot,G4ThreeVector(-(XwidthPad+XdistbetwPads)/2.,-(YwidthPad+YdistbetwPads)/2.,ZoffsetPad)),
         logicPadSi4,"Tina_PadSi_PadSi4",m_PadSi,false,4,false);

    // logicPadSi0->SetVisAttributes(m_VisPadSi);
    logicPadSi1->SetVisAttributes(m_VisPadSi);
    logicPadSi2->SetVisAttributes(m_VisPadSi);
    logicPadSi3->SetVisAttributes(m_VisPadSi);
    logicPadSi4->SetVisAttributes(m_VisPadSi);
    // logicPadSi0->SetSensitiveDetector(m_PadScorer);    
    logicPadSi1->SetSensitiveDetector(m_PadScorer);    
    logicPadSi2->SetSensitiveDetector(m_PadScorer);    
    logicPadSi3->SetSensitiveDetector(m_PadScorer);       
    logicPadSi4->SetSensitiveDetector(m_PadScorer);    

    m_PadSi->SetVisAttributes(G4Colour(1,1,0,0.25));         
  }
  return m_PadSi;
} */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Tina::BuildTTTPadDetector(){
  if(!m_TTTPad){

    // no frame defined for TTT since it is asymmetrical and hence not straightforward
    G4Box* solidTTT = new G4Box("Tina_TTT",XwidthTTT/2.,YwidthTTT/2.,ZwidthTTT/2.);
    G4Material* DetectorMaterial1 = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    G4LogicalVolume* logicTTT = new G4LogicalVolume(solidTTT,DetectorMaterial1,"logic_Tina_TTT",0,0,0);

    // PadSi replaced by four thick CsI blocks enclosed in single box
    G4Box* solidPadSi1 = new G4Box("Tina_PadSi1",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadSi2 = new G4Box("Tina_PadSi2",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadSi3 = new G4Box("Tina_PadSi3",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadSi4 = new G4Box("Tina_PadSi4",XwidthPad/2.,YwidthPad/2.,ZwidthPad/2.);
    G4Box* solidPadBox = new G4Box("Tina_PadFrame",XwidthBox/2.,YwidthBox/2.,ZwidthPad/2.-0.5);
    G4Material* DetectorMaterial2 = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
    //    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4LogicalVolume* logicPadSi1 = new G4LogicalVolume(solidPadSi1,DetectorMaterial2,"logic_Tina_PadSi1",0,0,0);
    G4LogicalVolume* logicPadSi2 = new G4LogicalVolume(solidPadSi2,DetectorMaterial2,"logic_Tina_PadSi2",0,0,0);
    G4LogicalVolume* logicPadSi3 = new G4LogicalVolume(solidPadSi3,DetectorMaterial2,"logic_Tina_PadSi3",0,0,0);
    G4LogicalVolume* logicPadSi4 = new G4LogicalVolume(solidPadSi4,DetectorMaterial2,"logic_Tina_PadSi4",0,0,0);
    G4LogicalVolume* logicPadBox = new G4LogicalVolume(solidPadBox,DetectorMaterial2,"logic_Tina_PadBox",0,0,0);
    //    G4LogicalVolume* logicPadBox = new G4LogicalVolume(solidPadBox,Vacuum,"logic_Tina_PadBox",0,0,0);

    m_TTTPad = new G4AssemblyVolume();

    // detector is constituted at centre of front face of TTT
    G4ThreeVector PosTTT(0,0,ZwidthTTT/2.);
    // there is a 2 mm offset in x between centre of TTT and centre of Pad
    // ZoffsetPad is distance between front face TTT and front face Pad
    G4ThreeVector PosPadSi1( (XwidthPad+XdistbetwPads)/2.-2., (YwidthPad+YdistbetwPads)/2.,ZwidthPad/2.+ZoffsetPad);
    G4ThreeVector PosPadSi2(-(XwidthPad+XdistbetwPads)/2.-2., (YwidthPad+YdistbetwPads)/2.,ZwidthPad/2.+ZoffsetPad);
    G4ThreeVector PosPadSi3(-(XwidthPad+XdistbetwPads)/2.-2.,-(YwidthPad+YdistbetwPads)/2.,ZwidthPad/2.+ZoffsetPad);
    G4ThreeVector PosPadSi4( (XwidthPad+XdistbetwPads)/2.-2.,-(YwidthPad+YdistbetwPads)/2.,ZwidthPad/2.+ZoffsetPad);
    // PadBox is pushed back by same amount as the pads themselves
    G4ThreeVector PosPadBox(-2.,0,ZwidthPad/2.+ZoffsetPad);
    //    G4ThreeVector PosPadBox(0,0,ZoffsetPad/2.);
    m_TTTPad->AddPlacedVolume(logicTTT,PosTTT,0);
    m_TTTPad->AddPlacedVolume(logicPadSi1,PosPadSi1,0);
    m_TTTPad->AddPlacedVolume(logicPadSi2,PosPadSi2,0);
    m_TTTPad->AddPlacedVolume(logicPadSi3,PosPadSi3,0);
    m_TTTPad->AddPlacedVolume(logicPadSi4,PosPadSi4,0);
    //    m_TTTPad->AddPlacedVolume(logicPadBox,PosPadBox,0);//*/

    logicTTT->SetVisAttributes(m_VisTTT);
    logicPadSi1->SetVisAttributes(m_VisPad);
    logicPadSi2->SetVisAttributes(m_VisPad);
    logicPadSi3->SetVisAttributes(m_VisPad);
    logicPadSi4->SetVisAttributes(m_VisPad);
    logicPadBox->SetVisAttributes(G4Colour(1,1,0,0.25));         
    logicTTT->SetSensitiveDetector(m_TTTScorer);
    logicPadSi1->SetSensitiveDetector(m_PadScorer);    
    logicPadSi2->SetSensitiveDetector(m_PadScorer);    
    logicPadSi3->SetSensitiveDetector(m_PadScorer);       
    logicPadSi4->SetSensitiveDetector(m_PadScorer);    
  }
  return m_TTTPad;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*G4AssemblyVolume* Tina::BuildYY1CsIDetector(){
  if(!m_YY1CsI){

    G4Box* solidType1A = new G4Box("Tina_Type1A",XwidthT1A/2.,YwidthT1A/2.,ZwidthT1A/2.);
    G4Trd* solidType1B = new G4Trd("Tina_Type1B",FwidthT1B/2.,BwidthT1B/2.,FwidthT1B/2.,BwidthT1B/2.,ZwidthT1B/2.);
    G4Box* solidType3  = new G4Box("Tina_Type3",XwidthT3/2.,YwidthT3/2.,ZwidthT3/2.);
    G4Material* DetectorMaterial1 = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
    G4LogicalVolume* logicType1A = new G4LogicalVolume(solidType1A,DetectorMaterial1,"logic_Tina_Type1A",0,0,0);
    G4LogicalVolume* logicType1B = new G4LogicalVolume(solidType1B,DetectorMaterial1,"logic_Tina_Type1B",0,0,0);
    G4LogicalVolume* logicType3  = new G4LogicalVolume(solidType3,DetectorMaterial1,"logic_Tina_Type3",0,0,0);//
    //  G4Trap* solidType1 = new G4Trap("Tina_Type1",Length/2, 0*deg, 0*deg, Height/2, BaseLarge/2, BaseSmall/2, 0*deg, Height/2, BaseLarge/2, BaseSmall/2, 0*deg);

    //  G4Material* DetectorMaterial1 = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
    //  G4LogicalVolume* logicType1 = new G4LogicalVolume(solidType1,DetectorMaterial1,"logic_Tina_Type1",0,0,0);
      
    // active area could also be defined as G4Trd with front=99.9x0.3 back=34.0x0.3 and z=85.0
    // G4Trd* solidYY1 = new G4Trd("Tina_YY1",Tina_NS::Width*49.95,Tina_NS::Width*17.,Tina_NS::Thickness*.15,Tina_NS::Thickness*.15,Tina_NS::Width*42.50);
    // dimensions of YY1 as G4Tubs are approximate and ignore that upper corners are cut away
    // YY1 detector is then constituted at phi=90°
    G4Tubs* solidYY1 = new G4Tubs("Tina_YY1",RinnerYY1,RouterYY1,ZwidthYY1/2.,phiYY1,alphaYY1);
    G4Tubs* solidFrame = new G4Tubs("Tina_Frame",RinnerFrame,RouterFrame,ZwidthFrame/2.,phiFrame,alphaFrame);
      
    G4Material* DetectorMaterial2 = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4LogicalVolume* logicYY1 = new G4LogicalVolume(solidYY1,DetectorMaterial2,"logic_Tina_YY1",0,0,0);
    G4LogicalVolume* logicFrame = new G4LogicalVolume(solidFrame,Vacuum,"logic_Tina_Frame",0,0,0);

    m_YY1CsI = new G4AssemblyVolume();

    // detector is constituted at centre of front face of Type1A
    G4ThreeVector PosType1A(0,0,ZwidthT1A/2.);
    G4ThreeVector PosType1B(0,0,ZwidthT1A+ZwidthT1B/2.);
    G4ThreeVector PosType3(0,-YwidthT1A/2.-YwidthT3/2.,ZwidthT3/2.);
    // shift YY1 upwards by 88.54-108.29=-19.75 in case it is constructed as G4Trd
    // G4ThreeVector PosYY1(0,-19.75,-6.85);   
    G4ThreeVector PosYY1(0,-YoffsetYY1,-ZwidthYY1/2.-ZoffsetYY1);   
    // rotate YY1 in case it is constructed as G4Trd
    // G4ThreeVector u = G4ThreeVector(1,0,0);
    // G4ThreeVector v = G4ThreeVector(0,0,1);   
    // G4ThreeVector w = G4ThreeVector(0,-1,0);
    // G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
    G4ThreeVector PosFrame(0,-YoffsetFrame,-ZwidthFrame/2.);//
      // detector is constituted at centre of front face of Type1A
    //  G4ThreeVector PosYY1(0,0,-ZwidthYY1/2.-ZoffsetYY1);
    //  G4ThreeVector PosType1(0,YoffsetYY1,ZwidthT1A/2.);
    //  G4ThreeVector PosType1A(0,YoffsetYY1,ZwidthT1A/2.);
    //  G4ThreeVector PosType1B(0,YoffsetYY1,ZwidthT1A+ZwidthT1B/2.);
    //  G4ThreeVector PosType3(0,YoffsetYY1-YwidthT1A/2.-YwidthT3/2.,ZwidthT3/2.);//
      // shift YY1 upwards by 88.54-108.29=-19.75 in case it is constructed as G4Trd
      // G4ThreeVector PosYY1(0,-19.75,-6.85);
      // rotate YY1 in case it is constructed as G4Trd
      // G4ThreeVector u = G4ThreeVector(1,0,0);
      // G4ThreeVector v = G4ThreeVector(0,0,1);
      // G4ThreeVector w = G4ThreeVector(0,-1,0);
      // G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
  //  G4ThreeVector PosFrame(0,0,-ZwidthFrame/2.);
  //  m_YY1CsI->AddPlacedVolume(logicType1,PosType1,0);
  //  m_YY1CsI->AddPlacedVolume(logicYY1,PosYY1,0);
  //  m_YY1CsI->AddPlacedVolume(logicFrame,PosFrame,0);

  //  logicType1->SetVisAttributes(m_VisCsI);
  //  logicYY1->SetVisAttributes(m_VisYY1);
  //  logicFrame->SetVisAttributes(G4Colour(1,1,0,0.25));
  //  logicType1->SetSensitiveDetector(m_CsIScorer);
  //  logicYY1->SetSensitiveDetector(m_YY1Scorer);
      m_YY1CsI->AddPlacedVolume(logicType1A,PosType1A,0);
    m_YY1CsI->AddPlacedVolume(logicType1B,PosType1B,0);
    m_YY1CsI->AddPlacedVolume(logicType3,PosType3,0);
    m_YY1CsI->AddPlacedVolume(logicYY1,PosYY1,0);
    m_YY1CsI->AddPlacedVolume(logicFrame,PosFrame,0);

    logicType1A->SetVisAttributes(m_VisCsI);
    logicType1B->SetVisAttributes(G4Colour(1,1,0,0.25));
    logicType3->SetVisAttributes(m_VisCsI);
    logicYY1->SetVisAttributes(m_VisYY1); 
    logicFrame->SetVisAttributes(G4Colour(1,1,0,0.25));
    logicType1A->SetSensitiveDetector(m_CsIScorer);
    logicType3->SetSensitiveDetector(m_CsIScorer);
    logicYY1->SetSensitiveDetector(m_YY1Scorer);//
  }
  return m_YY1CsI;
}//*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Tina::BuildYY1CsIDetector(){
  if(!m_YY1CsI){
    ////////////////////////////////////////////////////////////////
    ////////////// Starting Volume Definition //////////////////////
    ////////////////////////////////////////////////////////////////
    G4Material* DetectorMaterial1 = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
    G4Material* DetectorMaterial2 = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");
    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");

    // Definition of the volume containing the sensitive detector
    G4Trap* solidFrame = new G4Trap("solidFrame",
				    SolidWidth/2, 0*deg, 0*deg,
				    SolidHeight/2, SolidBaseLarge/2, SolidBaseSmall/2, 0*deg,
				    SolidHeight/2, SolidBaseLarge/2, SolidBaseSmall/2, 0*deg);//*/

    G4ThreeVector PosFrame(0,0,SolidWidth/2.);//

    G4LogicalVolume* logicFrame = new G4LogicalVolume(solidFrame, Vacuum, "logicFrame", 0, 0, 0);
        
    m_YY1CsI = new G4AssemblyVolume();
        
    G4VisAttributes* FrameVisAtt = new G4VisAttributes(G4Colour(0.90, 0.90, 0.90));
    FrameVisAtt->SetForceWireframe(true);
    logicFrame->SetVisAttributes(FrameVisAtt);
    //	m_YY1CsI->AddPlacedVolume(logicFrame,PosFrame,0);

    ////////////////////////////////////////////////////////////////
    /////////////////// Silicon detector Construction////////////////////
    ////////////////////////////////////////////////////////////////
    G4ThreeVector  PosYY1 = G4ThreeVector(0, (RinnerYY1+RouterYY1)/2., ZwidthYY1/2.);

    G4Tubs* solidYY1 = new G4Tubs("Tina_YY1",RinnerYY1,RouterYY1,ZwidthYY1/2.,phiYY1,alphaYY1);
    G4LogicalVolume* logicYY1 = new G4LogicalVolume(solidYY1, DetectorMaterial1, "logicYY1", 0, 0, 0);
        
    // Set First Stage sensible
    logicYY1->SetSensitiveDetector(m_YY1Scorer);
        
    ///Visualisation of FirstStage Strip
    logicYY1->SetVisAttributes(m_VisYY1);
    m_YY1CsI->AddPlacedVolume(logicYY1,PosYY1,0);

    ////////////////////////////////////////////////////////////////
    //////////////// CsI  Construction ////////////////////
    ////////////////////////////////////////////////////////////////
    double XShift=((CsIBaseLarge/2. - CsIBaseSmall/2.)/2. + CsIBaseSmall/2.)*2.+0.5;
    double Zrot=atan((CsIBaseLarge/2. - CsIBaseSmall/2.)*2./CsIHeight);
    G4ThreeVector  PosCsI0 = G4ThreeVector(-XShift,sin(Zrot)*CsIBaseSmall,ZwidthYY1 + CsIWidth/2. + ZoffsetYY1);
    G4ThreeVector  PosCsI1 = G4ThreeVector(0,0,ZwidthYY1 + CsIWidth/2. + ZoffsetYY1);
    G4ThreeVector  PosCsI2 = G4ThreeVector(XShift,sin(Zrot)*CsIBaseSmall,ZwidthYY1 + CsIWidth/2. + ZoffsetYY1);


    G4Trap* solidCsI0 = new G4Trap("solidCsI0",
				   CsIWidth/2, 0*deg, 0*deg,
				   CsIHeight/2, CsIBaseLarge/2, CsIBaseSmall/2, 0*deg,
				   CsIHeight/2, CsIBaseLarge/2, CsIBaseSmall/2, 0*deg);
    G4Trap* solidCsI1 = new G4Trap("solidCsI1",
				   CsIWidth/2, 0*deg, 0*deg,
				   CsIHeight/2, CsIBaseLarge/2, CsIBaseSmall/2, 0*deg,
				   CsIHeight/2, CsIBaseLarge/2, CsIBaseSmall/2, 0*deg);
    G4Trap* solidCsI2 = new G4Trap("solidCsI2",
				   CsIWidth/2, 0*deg, 0*deg,
				   CsIHeight/2, CsIBaseLarge/2, CsIBaseSmall/2, 0*deg,
				   CsIHeight/2, CsIBaseLarge/2, CsIBaseSmall/2, 0*deg);
        
        
    G4LogicalVolume* logicCsI0 = new G4LogicalVolume(solidCsI0,DetectorMaterial2,"logicCsI0",0,0,0);
    G4LogicalVolume* logicCsI1 = new G4LogicalVolume(solidCsI1,DetectorMaterial2,"logicCsI1",0,0,0);
    G4LogicalVolume* logicCsI2 = new G4LogicalVolume(solidCsI2,DetectorMaterial2,"logicCsI2",0,0,0);
        
    /// Placing the CsI on the Trapezoidal telescope
    m_YY1CsI->AddPlacedVolume(logicCsI1,PosCsI1,0);
        
    G4ThreeVector w = G4ThreeVector(0,0,-1);
    G4ThreeVector v = G4ThreeVector(-sin(Zrot),cos(Zrot),0);
    G4ThreeVector u = G4ThreeVector(cos(Zrot),sin(Zrot),0);
    m_YY1CsI->AddPlacedVolume(logicCsI2,PosCsI2,new G4RotationMatrix(u,v,w));
        
    v = G4ThreeVector(sin(Zrot),cos(Zrot),0);
    u = G4ThreeVector(cos(Zrot),-sin(Zrot),0);
    m_YY1CsI->AddPlacedVolume(logicCsI0,PosCsI0,new G4RotationMatrix(u,v,w));

    // Set Second Stage sensible
    logicCsI0->SetSensitiveDetector(m_CsIScorer);
    logicCsI1->SetSensitiveDetector(m_CsIScorer);
    logicCsI2->SetSensitiveDetector(m_CsIScorer);

    ///Visualisation of SecondStage Strip
    logicCsI0->SetVisAttributes(m_VisCsI);
    logicCsI1->SetVisAttributes(m_VisCsI);
    logicCsI2->SetVisAttributes(m_VisCsI);
  }
  return m_YY1CsI;
}
 //*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Tina");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","T","Shape"};
  // theta will be measured down from z, phi in xy plane
  // T is translation along z

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Tina " << i+1 <<  endl;    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      double Alpha = 0;
      if (blocks[i]->HasToken("Alpha")) blocks[i]->GetDouble("Alpha","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Alpha,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Tina " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double T = blocks[i]->GetDouble("T","mm");
      double Alpha = 0;
      if (blocks[i]->HasToken("Alpha")) blocks[i]->GetDouble("Alpha","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,T,Alpha,Shape);
    }
    else{
      cout << "Error: check your input file formatting" << endl;
      exit(1);
    }
  }
  std::cout << "read complete" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::ConstructDetector(G4LogicalVolume* world){
  std::cout << "start constuct detector" << std::endl;
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double Det_pos_x = m_R[i] * sin(m_Theta[i]) * cos(m_Phi[i]);
    G4double Det_pos_y = m_R[i] * sin(m_Theta[i]) * sin(m_Phi[i]);
    G4double Det_pos_z = m_R[i] * cos(m_Theta[i]);
    G4ThreeVector Det_pos = G4ThreeVector(Det_pos_x,Det_pos_y,Det_pos_z);

    /* if(m_Shape[i] == "TTT"){
    // move Det_pos from face to middle of active detector
    // G4ThreeVector Trans = G4ThreeVector(ZwidthTTT/2.*cos(m_Phi[i]),ZwidthTTT/2.*sin(m_Phi[i]),0);
    // for wall geometry use this 
    G4ThreeVector Trans = G4ThreeVector(0,0,ZwidthTTT/2.);
    Det_pos += Trans;
    G4ThreeVector u = G4ThreeVector(-sin(m_Phi[i]),cos(m_Phi[i]),0);
    G4ThreeVector v = G4ThreeVector(0,0,1);    
    G4ThreeVector w = G4ThreeVector(cos(m_Phi[i]),sin(m_Phi[i]),0);
    // G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
    // for wall geometry use this 
    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),BuildTTTDetector(),"Tina_TTT",world,false,i+1,false);
    }
    else if(m_Shape[i] == "PadSi"){
    // move Det_pos from face to middle of detector box 
    // G4ThreeVector Trans = G4ThreeVector(ZwidthBox/2.*cos(m_Phi[i]),ZwidthBox/2.*sin(m_Phi[i]),0);
    // for wall geometry use this 
    G4ThreeVector Trans = G4ThreeVector(0,0,ZwidthBox/2.);
    Det_pos += Trans;
    G4ThreeVector u = G4ThreeVector(-sin(m_Phi[i]),cos(m_Phi[i]),0);
    G4ThreeVector v = G4ThreeVector(0,0,1);    
    G4ThreeVector w = G4ThreeVector(cos(m_Phi[i]),sin(m_Phi[i]),0);
    // G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
    // for wall geometry use this 
    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),BuildPadSiDetector(),"Tina_PadSi",world,false,i+1,false);
    } */

    if(m_Shape[i] == "TTTPad"){
      // for box geometry use this
      G4ThreeVector u = G4ThreeVector(0,0,-1);
      G4ThreeVector v = G4ThreeVector(-sin(m_Phi[i]),cos(m_Phi[i]),0);
      G4ThreeVector w = G4ThreeVector(cos(m_Phi[i]),sin(m_Phi[i]),0);
      //      G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
      //       Rot->
      // the following is the same but perhaps more transparent
      /* G4RotationMatrix* Flip = new G4RotationMatrix;
	 Flip->rotateY(90.*deg);
	 G4RotationMatrix* Rot = new G4RotationMatrix;
	 Rot->rotateZ(m_Phi[i]);//*/
      // for cross geometry disable next line so there is no flip but only rotation
      // *Rot *= *Flip;
      // for wall geometry use this instead to disable also the rotation 
      G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
      G4ThreeVector Trans = G4ThreeVector(0,0,m_T[i]);
      Det_pos += Trans;
      BuildTTTPadDetector();

      m_TTTPad->MakeImprint(world,Det_pos,Rot,i+1,true);
      // iterator is equal to fPVStore.begin()
      std::vector< G4VPhysicalVolume * >::iterator it = m_TTTPad->GetVolumesIterator();
      unsigned int NbrImprints = m_TTTPad->GetImprintsCount();
      unsigned int NbrTotalPV = m_TTTPad->TotalImprintedVolumes(); 
      unsigned int NbrComponents = NbrTotalPV/NbrImprints;
      // set copy numbers of components of assembly volume to the current detector number
      for (it += (NbrImprints-1)*NbrComponents ; it <= m_TTTPad->GetVolumesIterator()+NbrTotalPV-1 ; it++){
        (*it)->SetCopyNo(i+1);
      }
    }
    else if(m_Shape[i] == "YY1CsI"){
      // central position b was taken ignoring 1 mm shift to position a        
      // build detector reference frame (u,v,w) with w pointing away from origin
      // Q is an intermediate vector that is shifted to theta+90° for same phi
      // warning: Q is the same as v so why is it defined again?
      G4double Qx = cos(m_Theta[i]) * cos(m_Phi[i]);
      G4double Qy = cos(m_Theta[i]) * sin(m_Phi[i]);
      G4double Qz = -sin(m_Theta[i]);
      G4ThreeVector Q = G4ThreeVector(Qx,Qy,Qz);
      G4ThreeVector u = Q.cross(Det_pos);
      G4ThreeVector v = Det_pos.cross(u);
      G4ThreeVector w = Det_pos.unit();
      u = u.unit();
      v = v.unit();
      G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
      // additional rotation to put YY1CsI in backwards direction 
      // warning: not clear why this got to be Z axis
      // G4RotationMatrix* RotBackwards = new G4RotationMatrix;
      // RotBackwards->rotateZ(180.*deg);
      // *Rot *= *RotBackwards;
      // once rotation is done, add translation T to Det_pos 
      G4ThreeVector Trans = G4ThreeVector(0,0,m_T[i]);
      Det_pos += Trans;
      BuildYY1CsIDetector();

      m_YY1CsI->MakeImprint(world,Det_pos,Rot,i+1,true);
      // iterator is equal to fPVStore.begin()
      std::vector< G4VPhysicalVolume * >::iterator it = m_YY1CsI->GetVolumesIterator();
      unsigned int NbrImprints = m_YY1CsI->GetImprintsCount();
      unsigned int NbrTotalPV = m_YY1CsI->TotalImprintedVolumes(); 
      unsigned int NbrComponents = NbrTotalPV/NbrImprints;
      // set copy numbers of components of assembly volume to the current detector number
      for (it += (NbrImprints-1)*NbrComponents ; it <= m_YY1CsI->GetVolumesIterator()+NbrTotalPV-1 ; it++){
        (*it)->SetCopyNo(i+1);
      }
    }
  }
  std::cout << "construct complete" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::InitializeRootOutput(){
  // add detector branch to the EventTree
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Tina")){
    pTree->Branch("Tina","TTinaData",&m_Event) ;
  }
  pTree->SetBranchAddress("Tina",&m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::ReadSensitive(const G4Event*){
  // read sensitive part and fill the Root tree
  m_Event->Clear();
 //   cout << "ReadSensitive" << endl;
  DSSDScorers::PS_Images* TTTScorer = (DSSDScorers::PS_Images*)m_TTTScorer->GetPrimitive(0);
  // in Must2 trig is the list of telescopes that got a Si trigger but probably not used anymore
  set<int> trig;
  unsigned int sizeFront = TTTScorer->GetFrontMult();
  unsigned int sizeBack  = TTTScorer->GetBackMult();
  //  m_Event->SetKE(TTTScorer->GetParticleKE());
  //  cout << "sizeFront=" << sizeFront << endl;
  //  cout << "sizeBack=" << sizeBack << endl;

  for (unsigned int i = 0; i < sizeBack; i++) {
    double energyX     = RandGauss::shoot(TTTScorer->GetEnergyBack(i), ResoTTT);
    int    detectorNbr = TTTScorer->GetDetectorBack(i);
    double timeX       = RandGauss::shoot(TTTScorer->GetTimeBack(i), ResoTime);
    if (energyX > 0.1 * keV) { // above threshold
      // pixel value at interaction point
      // sf warning: the png file numbers the strips blue=0...123 left to right so b+3 gives 3...126 with dead strips 0,1,127,128 
      unsigned int a, r, g, b;
      TTTScorer->GetARGBBack(i, a, r, g, b);
      //     b = b + 2;
      //   g = g + 2;
       if (r == 0) {
        trig.insert(detectorNbr);
	//   m_Event->SetTTTfrontEnergy(detectorNbr, b + 1, energyX);
	//  m_Event->SetTTTfrontTime(detectorNbr, b + 1, timeX);
        m_Event->SetTTTbackEnergy(detectorNbr, b+1, energyX);
	m_Event->SetTTTbackTime(detectorNbr, b+1, timeX);
      }
      else { 
        // interstrip X, keep maximum shared energy
        double rand = G4UniformRand();
        if (rand > 0.5) {
          energyX = rand * energyX;
          if (energyX > 0.1 * keV) {
            trig.insert(detectorNbr);
            m_Event->SetTTTbackEnergy(detectorNbr, b + 1, energyX);
            m_Event->SetTTTbackTime(detectorNbr, b + 1, timeX);
          }
        }
	else {
          energyX = (1 - rand) * energyX;
          if (energyX > 0.1 * keV) {
            trig.insert(detectorNbr);
            // sf warning: don't really understand where g+1 goes
            m_Event->SetTTTbackEnergy(detectorNbr, g + 1, energyX);
            m_Event->SetTTTbackTime(detectorNbr, g + 1, timeX);
          }
        }
      }
    }
  }

  for (unsigned int i = 0; i < sizeFront; i++) {
    double energyY     = RandGauss::shoot(TTTScorer->GetEnergyFront(i), ResoTTT);
    int    detectorNbr = TTTScorer->GetDetectorFront(i);
    double timeY       = RandGauss::shoot(TTTScorer->GetTimeFront(i), ResoTime);
    if (energyY > 0.1 * keV) { // above threshold
      // pixel value at interaction point
      // sf warning: the png file numbers the strips blue=0...123 top to bottom so b+3 gives 3...126 with dead strips 0,1,127,128 
      unsigned int a, r, g, b;
      TTTScorer->GetARGBFront(i, a, r, g, b);
      //     b = b + 2;
      //     g = g + 2;
      if (r == 0) {
        trig.insert(detectorNbr);
        m_Event->SetTTTfrontEnergy(detectorNbr, b + 1, energyY);
        m_Event->SetTTTfrontTime(detectorNbr, b + 1, timeY);
      } else { 
        // interstrip Y, keep both strips with shared energy
        double rand     = G4UniformRand();
        double energyY1 = rand * energyY;
        if (energyY1 > 0.1 * keV) {
          trig.insert(detectorNbr);
          m_Event->SetTTTfrontEnergy(detectorNbr, b + 1, energyY1);
          m_Event->SetTTTfrontTime(detectorNbr, b + 1, timeY);
        }
        if (energyY1 > 0.1 * keV) {
          trig.insert(detectorNbr);
          double energyY2 = (1 - rand) * energyY;
          // sf warning: don't really understand where g+1 goes
          m_Event->SetTTTfrontEnergy(detectorNbr, g + 1, energyY2);
          m_Event->SetTTTfrontTime(detectorNbr, g + 1, timeY);
        }
      }
    }
  }

  /* NPS::HitsMap<G4double*>* TTTHitMap;
  std::map<G4int, G4double**>::iterator TTT_itr;
  G4int TTTCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TinaTTTScorer/TTT");
  TTTHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(TTTCollectionID));

  for (TTT_itr = TTTHitMap->GetMap()->begin() ; TTT_itr != TTTHitMap->GetMap()->end() ; TTT_itr++){
    G4double* Info = *(TTT_itr->second);
    double Energy = Info[0];
    if(Energy>EnergyThreshold){
      double Time       = Info[1];
      int DetNbr        = (int) Info[7];
      int StripFront    = (int) Info[8];
      int StripBack     = (int) Info[9];    
      m_Event->SetTTTfrontEnergy(DetNbr,StripFront,RandGauss::shoot(Energy,ResoTTT));
      m_Event->SetTTTfrontTime(DetNbr,StripFront,RandGauss::shoot(Time,ResoTime));
      m_Event->SetTTTbackEnergy(DetNbr,StripBack,RandGauss::shoot(Energy,ResoTTT));
      m_Event->SetTTTbackTime(DetNbr,StripBack,RandGauss::shoot(Time,ResoTime));
      // ms_InterCoord->SetDetectedPositionX(Info[2]) ;
      // ms_InterCoord->SetDetectedPositionY(Info[3]) ;
      // ms_InterCoord->SetDetectedPositionZ(Info[4]) ;
      // ms_InterCoord->SetDetectedAngleTheta(Info[5]/deg) ;
      // ms_InterCoord->SetDetectedAnglePhi(Info[6]/deg) ;
    }
  }
  TTTHitMap->clear(); */

  CalorimeterScorers::PS_Calorimeter* PadScorer = (CalorimeterScorers::PS_Calorimeter*)m_PadScorer->GetPrimitive(0);
  unsigned int sizePad = PadScorer->GetMult();
  //  cout << "sizePad=" << sizePad << endl;
  for (unsigned int i = 0; i < sizePad; i++) {
      vector<unsigned int> level = PadScorer->GetLevel(i);
      //   cout << "level.size()=" << level.size() << endl;
      // for(int i=0;i<level.size();i++) cout << "level[" << i << "]=" << level[i] << " " ;
      // cout << endl;
      // detector number is contained in level [0], pad number in level[1]
      double EPad = RandGauss::shoot(PadScorer->GetEnergy(i), ResoPad);
      m_Event->SetPadEnergy(level[0], level[1], EPad); 
      double timePad = RandGauss::shoot(PadScorer->GetTime(i), ResoTime);
      m_Event->SetPadTime(level[0], level[1], timePad);
  }

  /* NPS::HitsMap<G4double*>* PadHitMap;
  std::map<G4int, G4double**>::iterator Pad_itr;
  G4int PadCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TinaPadScorer/Pad");
  PadHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(PadCollectionID));

  for (Pad_itr = PadHitMap->GetMap()->begin() ; Pad_itr != PadHitMap->GetMap()->end() ; Pad_itr++){
  // instead in Must2: Look for SiLi data in triggered telescope //////////////////////////////////////
  // std::set<int>::iterator itr;
  // for(itr=trig.begin();itr!=trig.end();itr++){
  //   for(SiLi_itr = SiLiHitMap->GetMap()->begin(); SiLi_itr!=SiLiHitMap->GetMap()->end() ; SiLi_itr++){
    G4double* Info = *(Pad_itr->second);
    // instead in Must2: matching telescope number /////////////////////////////////////////////////
    // if(Info[7]==*itr){
    double Energy = Info[0];
    if(Energy>EnergyThreshold){
      double Time       = Info[1];
      int DetNbr        = (int) Info[7];  
      m_Event->SetPadEnergy(DetNbr,RandGauss::shoot(Energy,ResoEnergy));
      m_Event->SetPadTime(DetNbr,RandGauss::shoot(Time,ResoTime));
      // m_Event->SetPadEnergy(Info[7],Info[8],E)); // Info8 is pad number
      // m_Event->SetPadTime(Info[7],Info[8],T));
      // ms_InterCoord->SetDetectedPositionX(Info[2]) ;
      // ms_InterCoord->SetDetectedPositionY(Info[3]) ;
      // ms_InterCoord->SetDetectedPositionZ(Info[4]) ;
      // ms_InterCoord->SetDetectedAngleTheta(Info[5]/deg) ;
      // ms_InterCoord->SetDetectedAnglePhi(Info[6]/deg) ;
    }
  }
  PadHitMap->clear(); */

  DSSDScorers::PS_Annular* YY1Scorer = (DSSDScorers::PS_Annular*)m_YY1Scorer->GetPrimitive(0);
  unsigned int sizeRing   = YY1Scorer->GetRingMult();
  unsigned int sizeSector = YY1Scorer->GetSectorMult();
  //  int MDetNbr=-1;
  //  if (sizeRing > 0 || sizeSector >> 0) 
  //    std::cout << "YY1: " << sizeRing << ", " << sizeSector << std::endl;
  for (unsigned int i = 0; i < sizeRing; i++) {
    double energyR     = RandGauss::shoot(YY1Scorer->GetEnergyRing(i), ResoYY1);
    int    detectorNbr = YY1Scorer->GetDetectorRing(i);
    int    ringNbr     = YY1Scorer->GetStripRing(i);
    double timeR       = RandGauss::shoot(YY1Scorer->GetTimeRing(i), ResoTime);
    //    std::cout << "  " << detectorNbr << ", " << ringNbr << std::endl;
    //  MDetNbr=detectorNbr;
    //if(detectorNbr==8) cout << "Before threshold, detectorNbr=" << detectorNbr << endl;
    if (energyR > 0.1 * keV) { // above threshold
       // if(detectorNbr==8) cout << "Above threshold, detectorNbr=" << detectorNbr << endl;
      if (ringNbr != 0) {
        m_Event->SetYY1ringEnergy(detectorNbr, ringNbr, energyR);
        m_Event->SetYY1ringTime(detectorNbr, ringNbr, timeR);
      }
      else { 
        // like TTT interstrip X, keep maximum shared energy
          // bm warning: obsolete, interstrip passes the first loop.
          // interstrip treatment should use the sizeRing variable.
          double rand = G4UniformRand();
        if (rand > 0.5) {
          energyR = rand * energyR;
          if (energyR > 0.1 * keV) {
            // sf warning: how to get ringNbr?
            // m_Event->SetYY1ringEnergy(detectorNbr, ringNbr, energyR);
            // m_Event->SetYY1ringTime(detectorNbr, ringNbr, timeR);
          }
        }
	else {
          energyR = (1 - rand) * energyR;
          if (energyR > 0.1 * keV) {
            // sf warning: how to get ringNbr?
            // m_Event->SetYY1ringEnergy(detectorNbr, ringNbr, energyR);
            // m_Event->SetYY1ringTime(detectorNbr, ringNbr, timeR);
          }
        }
      }
    }
  }

  // in fact YY1 is SSD so there is only 1 sector and hence no interstrip
  for (unsigned int i = 0; i < sizeSector; i++) {
    double energyS     = RandGauss::shoot(YY1Scorer->GetEnergySector(i), ResoYY1);
    int    detectorNbr = YY1Scorer->GetDetectorSector(i);
    int    sectorNbr   = YY1Scorer->GetStripSector(i);
    double timeS       = RandGauss::shoot(YY1Scorer->GetTimeSector(i), ResoTime);
   // if(detectorNbr==5) cout << "Before threshold, detectorNbr=" << detectorNbr << endl;
    if (energyS > 0.1 * keV) { // above threshold
     // if(detectorNbr==5) cout << "Above threshold, detectorNbr=" << detectorNbr << endl;
      if (sectorNbr != 0) {
        m_Event->SetYY1sectorEnergy(detectorNbr, sectorNbr, energyS);
        m_Event->SetYY1sectorTime(detectorNbr, sectorNbr, timeS);
      } else { 
        cout << "Warning: no interstrip treated for YY1 sectors" << endl;
      }
    }
  }

  /* NPS::HitsMap<G4double*>* YY1HitMap;
  std::map<G4int, G4double**>::iterator YY1_itr;
  G4int YY1CollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TinaYY1Scorer/YY1Scorer");
  YY1HitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(YY1CollectionID));
  
  // loop on the YY1 map
  for (YY1_itr = YY1HitMap->GetMap()->begin() ; YY1_itr != YY1HitMap->GetMap()->end() ; YY1_itr++){
    G4double* Info = *(YY1_itr->second);
    double Energy = Info[0];
    if(Energy>EnergyThreshold){
      double Time       = Info[1];
      int DetNbr        = (int) Info[7];
      int StripFront    = (int) Info[8];
      int StripBack     = (int) Info[9];  

      m_Event->SetYY1Energy(DetNbr,StripFront,RandGauss::shoot(Energy,ResoYY1));
      m_Event->SetYY1Time(DetNbr,StripFront,RandGauss::shoot(Time,ResoTime));
      // m_Event->SetBackEnergy(DetNbr,StripBack,RandGauss::shoot(Energy,ResoYY1));
      // m_Event->SetBackTime(DetNbr,StripBack,RandGauss::shoot(Time,ResoTime));
      // Interaction coordinates
      // ms_InterCoord->SetDetectedPositionX(Info[2]) ;
      // ms_InterCoord->SetDetectedPositionY(Info[3]) ;
      // ms_InterCoord->SetDetectedPositionZ(Info[4]) ;
      // ms_InterCoord->SetDetectedAngleTheta(Info[5]/deg) ;
      // ms_InterCoord->SetDetectedAnglePhi(Info[6]/deg) ;
    }
  }
  // clear map for next event
  YY1HitMap->clear(); */

  CalorimeterScorers::PS_Calorimeter* CsIScorer = (CalorimeterScorers::PS_Calorimeter*)m_CsIScorer->GetPrimitive(0);
  unsigned int sizeCsI = CsIScorer->GetMult();
  for (unsigned int i = 0; i < sizeCsI; i++) {
    vector<unsigned int> level = CsIScorer->GetLevel(i);
    double ECsI = RandGauss::shoot(CsIScorer->GetEnergy(i), ResoCsI);
   // if(level[0]==5) cout << "CsI level[0]=" << level[0] << "MDetNbr=" << MDetNbr << endl;
    m_Event->SetCsIEnergy(level[0], ECsI);
    //m_Event->SetCsIEnergy(level[0], level[1], NPL::EnergyToADC(ECsI, 0, 250, 8192, 16384)); // in case of crystal number
    double timeCsI = RandGauss::shoot(CsIScorer->GetTime(i), ResoTime);
    m_Event->SetCsITime(level[0], timeCsI);
    //m_Event->SetCsITime(level[0], level[1], NPL::EnergyToADC(timeCsI, 0, 1000, 16384, 8192));
  }

  /* NPS::HitsMap<G4double*>* CsIHitMap;
  std::map<G4int, G4double**>::iterator CsI_itr;
  G4int CsICollectionID = G4SDManager::GetSDMpointer()->GetCollectionID("TinaCsIScorer/CsI");
  CsIHitMap = (NPS::HitsMap<G4double*>*)(event->GetHCofThisEvent()->GetHC(CsICollectionID));

  // loop on the CsI map
  for (CsI_itr = CsIHitMap->GetMap()->begin() ; CsI_itr != CsIHitMap->GetMap()->end() ; CsI_itr++){
    G4double* Info = *(CsI_itr->second);
    double Energy = Info[0];
    if(Energy>EnergyThreshold){
      double Time = Info[1];
      int DetNbr  = (int) Info[7];

      m_Event->SetCsIEnergy(DetNbr,RandGauss::shoot(Energy,ResoEnergy));
      m_Event->SetCsITime(DetNbr,RandGauss::shoot(Time,ResoTime)); 
    }
  }
  // clear map for next event
  CsIHitMap->clear(); */
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Tina::InitializeScorers() {
  // this check is necessary in case the geometry is reloaded
  // if the multifunctional detector does not exist then it is created
  bool already_exist = false; 
  m_TTTScorer = CheckScorer("TinaTTTScorer",already_exist);
  m_PadScorer = CheckScorer("TinaPadScorer",already_exist);
  m_YY1Scorer = CheckScorer("TinaYY1Scorer",already_exist);
  m_CsIScorer = CheckScorer("TinaCsIScorer",already_exist);
  if(already_exist) 
    return ;

  // otherwise the scorer is initialised and registered to the multifunctional detector

  // Level=0, StripPlaneLength=100.5, StripPlaneWidth=100.5, 
  // NumberOfStripLength=128, NumberOfStripWidth=128, depth=0, axis=xy in local coordinates relative to centre of volume
  // G4VPrimitiveScorer* TTTScorer = new DSSDSCORERS::PS_Silicon_Rectangle("TTTScorer",0,100.5*mm,100.5*mm,128,128,0,"xy");
  string nptool = getenv("NPTOOL");
  // the png files define front strips as X and back strips as Y 
  // however in PixelToIndex in npl/Core/NPImage.cxx, index was running vertically, changed this to horizontally
  // therefore when mapping 3-dim (xyz) onto 2-dim (xy), x<>-x and y <>-y ie pixel (1,1) comes at bottom right
  // scalingFront=0.01, scalingBack=0.01, centerOffsetX=0, centerOffsetY=0, depth=0
  G4VPrimitiveScorer* TTTScorer = new DSSDScorers::PS_Images(
      "TTTScorer", nptool + "/NPLib/Detectors/Tina/resources/maskX.png",
      nptool + "/NPLib/Detectors/Tina/resources/maskY.png", 0.01*YwidthActiveTTT/Nstrips,0.01*XwidthActiveTTT/Nstrips,0,0,0xffff0000,0);
  G4VPrimitiveScorer* InterScorer = new InteractionScorers::PS_Interactions("TTTScorer", ms_InterCoord,0);
  // 1 means one level up ie the mother volume, 0 means current volume
  // CalorimeterScorers.cc GetCopyNumber(m_NestingLevel[i])) will first read GetCopyNumber(1) ie mother then GetCopyNumber(0) ie pad
//  vector<int> Pad_nesting = {1,0};
  vector<int> Pad_nesting = {0,1};
  G4VPrimitiveScorer* PadScorer = new CalorimeterScorers::PS_Calorimeter("PadScorer",Pad_nesting,0);

  G4VPrimitiveScorer* YY1Scorer = new DSSDScorers::PS_Annular(
      "YY1Scorer",0,RinnerYY1,RouterYY1,phiYY1,phiYY1+alphaYY1,Nfrontrings,Nbacksectors,1,0);

  vector<int> CsI_nesting = {0};
  G4VPrimitiveScorer* CsIScorer = new CalorimeterScorers::PS_Calorimeter("CsIScorer",CsI_nesting,0);

  m_TTTScorer->RegisterPrimitive(TTTScorer);
  m_TTTScorer->RegisterPrimitive(InterScorer);
  m_PadScorer->RegisterPrimitive(PadScorer);
  m_YY1Scorer->RegisterPrimitive(YY1Scorer);
  m_CsIScorer->RegisterPrimitive(CsIScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_TTTScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_PadScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_YY1Scorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CsIScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// construct method to be passed to the DetectorFactory
NPS::VDetector* Tina::Construct(){
  return  (NPS::VDetector*) new Tina();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// register the construct method to the factory                
extern"C" {
  class proxy_nps_Tina{
    public:
      proxy_nps_Tina(){
        NPS::DetectorFactory::getInstance()->AddToken("Tina","Tina");
        NPS::DetectorFactory::getInstance()->AddDetector("Tina",Tina::Construct);
      }
  };
  proxy_nps_Tina p_nps_Tina;
}
