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
 *  This class describe the GRAPE Germanium array                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
#include <fstream>
//G4 Geometry object
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4ThreeVector.hh"

// NPS
#include "GRAPE.hh"
#include "CalorimeterScorers.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"

// NPL
#include "NPOptionManager.h"
#include "RootOutput.h"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace GRAPENS {
  const double EnergyThreshold = 10*keV;
  //const double ResoTime = 4.5*ns ;  //not used
  const double ResoEnergy = 2.*keV ;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// GRAPE Specific Method
GRAPE::GRAPE(){
  InitializeMaterial();
  m_GRAPEData = new TGRAPEData();

  BlueVisAtt   = new G4VisAttributes(G4Colour(0, 0, 1)) ;
  GreenVisAtt  = new G4VisAttributes(G4Colour(0, 1, 0)) ;
  RedVisAtt    = new G4VisAttributes(G4Colour(1, 0, 0)) ;
  WhiteVisAtt  = new G4VisAttributes(G4Colour(1, 1, 1)) ;
  TrGreyVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.5)) ;

  m_LogicGRAPE = 0;

}

GRAPE::~GRAPE(){
  //  if (m_MaterialVacuum)  delete m_MaterialVacuum;
  //  if (m_MaterialGe) delete m_MaterialGe;
  //  if (m_MaterialAl) delete m_MaterialAl;
  //  if (m_MaterialCu) delete m_MaterialCu;
  //  if (m_MaterialN2) delete m_MaterialN2;
  //  if (m_MaterialC) delete m_MaterialC;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void GRAPE::ReadConfiguration(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GRAPE");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " free clovers found " << endl; 

  vector<string> token = {"GRAPEID","R","Theta","Phi","Beta"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      vector<double> beta = blocks[i]->GetVectorDouble("Beta","deg");
      int     id = blocks[i]->GetInt("GRAPEID");
      AddGRAPEFreePosition(id,R,Theta,Phi,beta[0],beta[1],beta[2]);
    }

    else{
      cout << "Warning: check your input file formatting " << endl;
    }
  }

  blocks.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Return a G4VSolid modeling the Crystal
//G4LogicalVolume* GRAPE::ConstructCrystal(){
//  return NULL;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Return a G4VSolid modeling the Capsule
//G4LogicalVolume* GRAPE::ConstructCapsule(){
//  return NULL;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//G4LogicalVolume* GRAPE::ConstructDewar(){
//  return NULL; }
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Return a G4VSolid modeling the BGO
//G4LogicalVolume* GRAPE::ConstructBGO(){
//  return NULL;}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Return a clover in the configuration given by option (not use a the moment)
void GRAPE::ConstructGRAPE(){	

  G4bool checkOverlaps = true;
  
  if (!m_LogicGRAPE) {
    G4VSolid *sHausingHexBase1; { // outer
      const G4int nz = 2;
      const G4int nss = 6;
      G4double z[nz]   = {-32.5*mm, 32.5*mm};
      G4double rin[nz] = {0.*mm, 0.*mm};
      G4double rout[nz] = {50.*mm, 50.*mm};
      sHausingHexBase1 = new G4Polyhedra("sHausingHexBase1", 30.*deg, 360.*deg,
					 nss, nz, z, rin, rout); }
    G4VSolid *sHausingTrd = // trd hausing
      new G4Trd("sHausingTrd",28.*mm,28.*mm,
		50.*2/sqrt(3.)*mm, 50.*2/sqrt(3.)*40./100.,30.*mm);
    // hausing hex 
    G4RotationMatrix rmY90;
    rmY90.rotateY(90.*deg);

    G4VSolid *sHausingHex =
      new G4UnionSolid("sHausingHex",sHausingHexBase1,sHausingTrd,
		       G4Transform3D(rmY90,G4ThreeVector(30.*mm, 0., 0.)));

    G4VSolid *sPreAmpHausing; {    // hausing pre-amp.
      const G4int nz = 15;
      G4double z[nz] = 
	{ 50.*mm, 60.*mm, 60.*mm, 62.*mm, 62.*mm,
	  100.*mm, 150.*mm, 150.*mm, 290.*mm, 290.*mm, 
	  334.*mm,334.*mm,340.*mm,340.*mm,607.*mm };
      G4double rin[nz] = 
	{ 0.*mm,0.*mm,0.*mm,0.*mm,0.*mm,
	  0.*mm,0.*mm,0.*mm,0.*mm,0.*mm,
	  0.*mm,0.*mm,0.*mm,0.*mm,0.*mm };
      G4double rout[nz] = 
	{ 20.*mm, 20.*mm, 23.09*mm, 23.09*mm, 20.*mm,
	  20.*mm, 50.*mm, 70.*mm, 70.*mm, 27.5*mm, 
	  27.5*mm, 30.*mm, 30.*mm, 111.*mm, 111.*mm };
      sPreAmpHausing = new G4Polycone("sPreAmpHausing",0.*deg, 360.*deg,
				      nz, z, rin, rout); }
    G4VSolid *sHausing =
      //new G4UnionSolid("sHausing",sHausingHex,sPreAmpHausing, // orig
      new G4UnionSolid("sHausing",sPreAmpHausing,sHausingHex, // orig
		       G4Transform3D(rmY90,G4ThreeVector()));
   
    ////////////////////////////////////////////////////////////
    // Vacuum 
    G4VSolid *sHeadVacuum; {// same as sHauseingHexBase2 
      const G4int nz = 2;
      const G4int nss = 6;
      G4double z[nz]   = {-30.5*mm, 30.5*mm};
      G4double rin[nz] = {0.*mm, 0.*mm};
      G4double rout[nz] = {48.*mm, 48.*mm};
      sHeadVacuum = new G4Polyhedra("sHeadVacuum", 30.*deg, 360.*deg,
				    nss, nz, z, rin, rout); }
      
    G4VSolid *sPreAmpVacuum; {
      const G4int nz = 10;
      G4double z[nz] =
	{ 157.*mm, 158.*mm, 158.*mm, 166.*mm, 166.*mm,
	  245.*mm, 245.*mm, 252.*mm, 252.*mm, 285.*mm };
      G4double rin[nz] = 
	{ 14.5*mm, 14.5*mm, 14.5*mm, 14.5*mm, 14.5*mm,
	  14.5*mm, 14.5*mm, 14.5*mm, 14.5*mm, 14.5*mm };
      G4double rout[nz] =
	{ 68.*mm, 68.*mm, 68.*mm, 68.*mm, 68.*mm, 
	  68.*mm, 68.*mm, 68.*mm, 65.*mm, 65.*mm } ;
      sPreAmpVacuum = 
	new G4Polycone("sPreAmpVacuum",0.*deg,360.*deg,nz,z,rin,rout);}
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Cooling parts
    G4VSolid *sDewarColor = new G4Tubs("sDewarColor",106.*mm,111.*mm,133.5*mm,
				       0.*deg,360.*deg);
    G4VSolid *sInnerDewar = new G4Tubs("sInnerDewar",0.,75.*mm,105.*mm,
				       0.*deg,360.*deg); // inner dewar
    G4VSolid *sLiquid = new G4Tubs("sLiquid",0.,70.*mm,100.*mm,0.*deg,360.*deg);
      
    G4VSolid *sColdFinger; { // cold finger (Cu)
      const G4int nz=2;
      G4double z[nz] = { 50.*mm, 368.5*mm };
      G4double rin[nz] = { 0., 0. };
      G4double rout[nz] = { 9.*mm, 9.*mm };
      sColdFinger = new G4Polycone("sColdFinger",0.*deg,360.*deg,
				   nz,z,rin,rout); }
    G4VSolid *sColdCase = new G4Tubs("sColdCase",35.*mm,40.*mm,20.*mm,
				     0.*deg,360.*deg); // cold case
    ////////////////////////////////////////////////////////////   

    ////////////////////////////////////////////////////////////
    // Ge Grystal
    G4VSolid *sGeCrystal = new G4Tubs("sGeCrystal",0.*mm,35.*mm,23.*mm,
				      0.*deg, 360.*deg);
    // Ge segment
    G4double rmaxSD = 32.*mm;
    G4VSolid *sGeSensitive = new G4Tubs("sGeSensitive",0.*mm,rmaxSD,20.*mm,
                                        0.*deg, 360.*deg);
    // cathode? 
    // todo: what is the material and thickness of the cathode?
    // is there other material between the crystals?
    G4VSolid *sGeCathode = new G4Tubs("sGeCathode",0.*mm, rmaxSD,0.1*mm,
				      0.*deg, 360.*deg);
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Logical Volumes
    G4LogicalVolume *lHausing =
      new G4LogicalVolume(sHausing,m_MaterialAl,"lHausing",0,0,0);
    G4VisAttributes *hausingVisAtt = new G4VisAttributes();
    hausingVisAtt->SetForceWireframe(true);
      
    G4LogicalVolume *lHeadVacuum = 
      new G4LogicalVolume(sHeadVacuum,m_MaterialVacuum,"lHeadVacuum",0,0,0);
    lHeadVacuum->SetVisAttributes(hausingVisAtt);
    G4LogicalVolume *lPreAmpVacuum =
      new G4LogicalVolume(sPreAmpVacuum,m_MaterialVacuum,"lPreAmpVacuum",0,0,0);
    G4LogicalVolume *lInnerDewar =
      new G4LogicalVolume(sInnerDewar,m_MaterialAl,"lInnerDewar",0,0,0);
    G4LogicalVolume *lDewarColor =
      new G4LogicalVolume(sDewarColor,m_MaterialAl,"lDewarColor",0,0,0);
    G4LogicalVolume *lLiquid =
      new G4LogicalVolume(sLiquid,m_MaterialVacuum,"lLiquid",0,0,0);
    G4LogicalVolume *lColdFinger =
      new G4LogicalVolume(sColdFinger,m_MaterialCu,"lColdFinger",0,0,0);
    G4LogicalVolume *lColdCase =
      new G4LogicalVolume(sColdCase,m_MaterialAl,"lColdCase",0,0,0);
    lColdCase->SetVisAttributes(hausingVisAtt);

    G4LogicalVolume *lGeCrystal =
      //new G4LogicalVolume(sGeCrystal,Ge,"lGeCrystal",0,0,0);
      new G4LogicalVolume(sGeCrystal,m_MaterialVacuum,"lGeCrystal",0,0,0);
    G4LogicalVolume *lGeCathode =
      new G4LogicalVolume(sGeCathode,m_MaterialC,"lGeCathode",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),lGeCathode, "pCathode",
		      lGeCrystal,false,0, checkOverlaps);
      
    //    G4LogicalVolume *lGeSensitive =
    //      new G4LogicalVolume(sGeSensitive,m_MaterialGe,"lGeSensitive",0,0,0);
    //    new G4PVPlacement(G4Transform3D(),lGeSensitive,"GeSensitive",lGeCrystal,true,0);
    //    lGeSensitive->SetSensitiveDetector(m_HPGeScorer);

    // segmentation
    G4double offs = 20.*mm + (rmaxSD-30.*mm)/2.;
    G4ThreeVector *segPos = new G4ThreeVector[18];
    for (G4int i = 0 ; i < 18 ; i++) 
      segPos[i] = G4ThreeVector(-((i%3)-1)*offs*mm,
				(1 - abs(int((i - 8.5)/3)))*offs*mm,
				((i/9)*2-1)*12.50*mm);
    G4double longLength = 10.*mm + (rmaxSD - 30.*mm)/2.;
    G4double defaultLength = 10.*mm;
    for (G4int i = 0 ; i < 18 ; i++) {
      G4VSolid *sSegBox;
      if      (i%9 == 4)
	sSegBox = new G4Box("sSegBox",defaultLength,defaultLength,defaultLength);
      else if ((i%9)%2 == 0)
	sSegBox = new G4Box("sSegBox",longLength,longLength,defaultLength);
      else if (i%9 == 1 || i %9 == 7)
	sSegBox = new G4Box("sSegBox",defaultLength,longLength,defaultLength);
      else 
	sSegBox = new G4Box("sSegBox",longLength,defaultLength,defaultLength);

      char sName[20],lName[20],pName[20];
      sprintf(sName,"sSeg%02d", i);
      sprintf(lName,"lSeg%02d", i);
      sprintf(pName,"pSeg%02d", i);
      G4VSolid *sSeg = new G4IntersectionSolid(sName,sGeSensitive,sSegBox,0,segPos[i]);
      G4LogicalVolume *lSeg =new G4LogicalVolume(sSeg,m_MaterialGe,lName,0,0,0);
      lSeg->SetSensitiveDetector(m_GRAPEScorer);
      new G4PVPlacement(0,G4ThreeVector(),lSeg,pName,lGeCrystal,true, i+1, checkOverlaps);}
    // visualize
    lDewarColor->SetVisAttributes(new G4VisAttributes(G4Colour(0.71,0.93,0.71)));
   
    // placements
    new G4PVPlacement(0,G4ThreeVector(),lLiquid,"Liquid",lInnerDewar,false,0, checkOverlaps);
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(),G4ThreeVector(0.,0.,473.5*mm)),
		      lInnerDewar,"InnerDewar",lHausing,false,0, checkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(),lGeCrystal,"GeCrystal",
		      lHeadVacuum,false,0,checkOverlaps); // orig
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(),G4ThreeVector()),
		      lColdCase,"ColdCase", lHeadVacuum,false,0, checkOverlaps); // orig
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(),G4ThreeVector()),
		      lColdFinger,"ColdFinger",lHausing,false,0, checkOverlaps);
    new G4PVPlacement(G4Transform3D(rmY90,G4ThreeVector()),
		      lHeadVacuum,"HeadVacuum", lHausing,false,0, checkOverlaps);
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(),G4ThreeVector()),
		      lPreAmpVacuum,"PreAmpVacuum",lHausing,false,0);
    new G4PVPlacement(G4Transform3D(G4RotationMatrix(),G4ThreeVector(0.,0.,473.5*mm)),
    		      lDewarColor,"DewarColor",lHausing,false,0, checkOverlaps);

    m_LogicGRAPE = lHausing; }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void GRAPE::ConstructDetector(G4LogicalVolume* world){
  ConstructGRAPE();

  //  std::ofstream fout("detpos.txt");
  G4RotationMatrix* DetectorRotation = new G4RotationMatrix(0,0,0);
  for (unsigned int i = 0 ;  i < m_GRAPEId.size(); i++) {
    
    // Constructing the Detector referential and the transition matrix
    G4ThreeVector U,V,W;
    G4double wX = sin(m_Theta[i]) * cos(m_Phi[i]) ;
    G4double wY = sin(m_Theta[i]) * sin(m_Phi[i]) ;
    G4double wZ = cos(m_Theta[i]);
    W = G4ThreeVector(wX, wY, wZ) ;

    // vector parallel to one axis of the entrance plane
    G4double vX = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double vY = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double vZ = -sin(m_Theta[i]);
    V = G4ThreeVector(vX, vY, vZ);

    W = W.unit();
    U = V.cross(W);
    U = U.unit();
    V = W.cross(U);
    V = V.unit();
    // Passage Matrix from Lab Referential to GRAPE Referential
    delete DetectorRotation;
    DetectorRotation = new G4RotationMatrix(U, V, W);

    DetectorRotation->rotate(m_BetaX[i], U);
    DetectorRotation->rotate(m_BetaY[i], V);
    DetectorRotation->rotate(m_BetaZ[i], W);
    G4ThreeVector DetectorPosition = m_R[i]*W;
  
    G4VPhysicalVolume *temp = new G4PVPlacement(G4Transform3D(*DetectorRotation, DetectorPosition),
						m_LogicGRAPE,"GRAPE",world,false,m_GRAPEId[i]);
    /*
    G4ThreeVector pos = temp->GetTranslation();
    //    std::cout << "( " << pos.x() <<", " << pos.y() <<", " << pos.z() << ")" << std::endl;

    fout << " Segment Position for Det " << (i+1) << std::endl;
    G4double offs = 20.*mm + (32.*mm-30.*mm)/2.;
    G4ThreeVector *segPos = new G4ThreeVector[18];
    for (G4int j = 0 ; j < 18 ; j++) {
      segPos[j] = G4ThreeVector(-((j%3)-1)*offs*mm,
				(1 - abs(int((j - 8.5)/3)))*offs*mm,
				((j/9)*2-1)*12.50*mm);
      G4ThreeVector rotPos = segPos[j].transform(*DetectorRotation);
      rotPos += DetectorPosition;
      fout << j << ": ( " << rotPos.x() <<", " << rotPos.y() <<", " << rotPos.z() << ")" << std::endl;
      }*/
    
  }
  //  fout.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add clover at the standard position of the array
// Take as argument the standard clover Id.
void GRAPE::AddGRAPEStandard(vector<int> GRAPEId){

  //////////////////////////////////////////////////
  // standard location
  double R_std[18] = {148, 138, 138, 138, 138, 148,
		      126, 101, 101, 101, 101, 126,
		      148, 138, 138, 138, 138, 148 };
  double Theta_std[18] = { 130, 125, 125, 125, 125, 130,
			    90,  90,  90,  90,  90,  90,
			    55,  55,  55,  55,  55,  55 };
  double Phi_std[18] = { 30, 90, 150, 210, 270, 330,
			 0, 60, 120, 180, 240, 300,
			 30, 90, 150, 210, 270, 330 };
  double betaX_std[18] = { 17, 20, 20, 20, 20, 17,
			   0, 0, 0, 0, 0, 0,
			   -17, -20, -20, -20, -20, -17 };
  double betaY_std[18] = { 0. };
  double betaZ_std[18] = { 0. };
  //////////////////////////////////////////////////
  
  for (unsigned int i = 0 ;  i < GRAPEId.size(); i++) {
    if (GRAPEId[i] < 1 || GRAPEId[i] > 18) continue;
    m_GRAPEId.push_back(GRAPEId[i]);
    m_R.push_back(R_std[GRAPEId[i] - 1]);
    m_Theta.push_back(Theta_std[GRAPEId[i] - 1]*deg);
    m_Phi.push_back(Phi_std[GRAPEId[i] - 1]*deg);
    m_BetaX.push_back(betaX_std[GRAPEId[i] - 1]);
    m_BetaY.push_back(betaY_std[GRAPEId[i] - 1]);
    m_BetaZ.push_back(betaZ_std[GRAPEId[i] - 1]);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add clover at a free position in space with coordinate
// in spherical coordinate
// Beta are the three angles of rotation in the GRAPE frame
void GRAPE::AddGRAPEFreePosition(int GRAPEId,double R,double Theta,double Phi,double BetaX,double BetaY,double BetaZ){

  m_GRAPEId.push_back(GRAPEId);
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_BetaX.push_back(BetaX);
  m_BetaY.push_back(BetaY);
  m_BetaZ.push_back(BetaZ);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void GRAPE::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("GRAPE")){
    pTree->Branch("GRAPE", "TGRAPEData", &m_GRAPEData) ;
  }
  pTree->SetBranchAddress("GRAPE", &m_GRAPEData) ;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void GRAPE::ReadSensitive(const G4Event*){ // event){
  m_GRAPEData->Clear();

  ////////////////////////////////////////////////////////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_GRAPEScorer->GetPrimitive(0);
  unsigned int size = Scorer->GetMult();
  //  if (size > 0) cout << size << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i), GRAPENS::ResoEnergy);
    if ( Energy > GRAPENS::EnergyThreshold ){
      double Time = Scorer->GetTime(i);
      int SegmentNbr = Scorer->GetLevel(i)[0];
      int GRAPENbr = Scorer->GetLevel(i)[1];
      m_GRAPEData->SetGeGRAPENbr(GRAPENbr);
      m_GRAPEData->SetGeCrystalNbr((SegmentNbr <= 9 ? 0 : 1));
      m_GRAPEData->SetGeSegmentNbr(SegmentNbr);
      m_GRAPEData->SetGeEnergy(Energy);
      m_GRAPEData->SetGeTimeCFD(Time);}}
  ////////////////////////////////////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void GRAPE::InitializeScorers(){
	//Look for previous definition of the scorer (geometry reload)
	//n.b. calls new G4MultiFunctionalDetector("GRAPE_CoreScorer");
  bool already_exist = false;
  m_GRAPEScorer = CheckScorer("GRAPEScorer",already_exist);

  // if the scorer were created previously nothing else need to be made
  if(already_exist) return;

  vector<int> level({0, 3}); 
  
  m_GRAPEScorer->RegisterPrimitive(new CalorimeterScorers::PS_Calorimeter("GRAPECry",level, 0)); 
  G4SDManager::GetSDMpointer()->AddNewDetector(m_GRAPEScorer);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
/////////////////Material Definition ///////////////////////////
////////////////////////////////////////////////////////////////
void GRAPE::InitializeMaterial(){
  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_MaterialGe= MaterialManager::getInstance()->GetMaterialFromLibrary("Ge");
  m_MaterialAl= MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_MaterialCu= MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
  m_MaterialC = MaterialManager::getInstance()->GetMaterialFromLibrary("C"); 
  m_MaterialN2= MaterialManager::getInstance()->GetMaterialFromLibrary("N2_liquid"); 
}


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* GRAPE::Construct(){
  return  (NPS::VDetector*) new GRAPE();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_grape{
    public:
      proxy_nps_grape(){
        NPS::DetectorFactory::getInstance()->AddToken("GRAPE","GRAPE");
        NPS::DetectorFactory::getInstance()->AddDetector("GRAPE",GRAPE::Construct);
      }
  };

  proxy_nps_grape p_nps_grape;
}
