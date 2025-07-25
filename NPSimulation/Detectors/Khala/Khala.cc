/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: M. Labiche  contact address: marc.labiche@stfc.ac.uk     *
 *                                                                           *
 * Creation Date  : December 2009                                            *
 * Last update    : December 2014                                            *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe the Khala scintillator array                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
using namespace std;

//Geant4
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4Colour.hh"

// NPS
#include "Khala.hh"
using namespace KHALA;

#include "CalorimeterScorers.hh"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
// NPL
#include "NPOptionManager.h"
#include "RootOutput.h"

//ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

// CLHEP header
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandFlat.h"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Khala Specific Method
Khala::Khala(){
  m_Event = new TKhalaData();

  // Blue
  m_LaBr3VisAtt = new G4VisAttributes(G4Colour(0, 0.5, 1));

  // Grey
  m_PMTVisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 0.6));

  // Grey wireframe
  m_DetectorCasingVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5,0.2));



  m_LogicalDetector = 0;
  m_LaBr3Scorer = 0 ;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Khala::~Khala(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Khala::AddDetector(G4ThreeVector Pos, G4RotationMatrix* r){
	m_Pos.push_back(Pos);
	m_Rot.push_back(r);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Khala::AddDetector(G4ThreeVector Pos1, G4ThreeVector Pos2, G4ThreeVector Pos3, G4ThreeVector Pos4){
  G4ThreeVector Pos=(Pos1+Pos2+Pos3+Pos4)/4.;
  G4ThreeVector u = Pos1-Pos2;
  G4ThreeVector v = Pos1-Pos4;
  u = u.unit(); v = v.unit();
  G4ThreeVector w = Pos.unit();
  Pos = Pos + w*Length*0.5;

  m_Pos.push_back(Pos);
  m_Rot.push_back(new G4RotationMatrix(u,v,w));
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Khala::AddDetector(G4ThreeVector Pos, double beta_u, double beta_v, double beta_w){
  double Theta = Pos.theta();
  double Phi = Pos.phi();

  // vector parallel to one axis of silicon plane
  G4double ii = cos(Theta / rad) * cos(Phi / rad);
  G4double jj = cos(Theta / rad) * sin(Phi / rad);
  G4double kk = -sin(Theta / rad);
  G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

  G4ThreeVector w = Pos.unit();
  G4ThreeVector u = w.cross(Y);
  G4ThreeVector v = w.cross(u);
  v = v.unit();
  u = u.unit();

  G4RotationMatrix* r = new G4RotationMatrix(u,v,w);
  r->rotate(beta_u,u);
  r->rotate(beta_v,v);
  r->rotate(beta_w,w);

  m_Pos.push_back(Pos);
  m_Rot.push_back(r);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Khala::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Khala");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 
  for(unsigned int i  = 0 ; i < blocks.size() ; i++){
    // Cartesian Case
    vector<string> cart = {"A","B","C","D"};

    // Spherical Case
    vector<string> sphe= {"R","THETA","PHI","BETA"};

    // Dynamite Case
    vector<string> dyna= {"X", "Y", "Z", "THETA", "PHI"};

    if(blocks[i]->HasTokenList(cart)){
      cout << endl << "////  Khala " << i+1 << endl;
      G4ThreeVector A = NPS::ConvertVector(blocks[i]->GetTVector3("A","mm"));
      G4ThreeVector B = NPS::ConvertVector(blocks[i]->GetTVector3("B","mm"));
      G4ThreeVector C = NPS::ConvertVector(blocks[i]->GetTVector3("C","mm"));
      G4ThreeVector D = NPS::ConvertVector(blocks[i]->GetTVector3("D","mm"));
      AddDetector(A,B,C,D) ;
    }

    else if(blocks[i]->HasTokenList(sphe)){
      cout << endl << "////  Khala " << i+1 << endl;
      double Theta = blocks[i]->GetDouble("THETA","deg");
      double Phi= blocks[i]->GetDouble("PHI","deg");
      double R = blocks[i]->GetDouble("R","mm");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA","deg");
      R = R +  0.5*Length;
      G4ThreeVector Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
      //G4ThreeVector Pos(R*cos(Theta)*cos(Phi),R*cos(Theta)*sin(Phi),R*sin(Theta));
      AddDetector(Pos,beta[0],beta[1],beta[2]);
    }else if(blocks[i]->HasTokenList(dyna)){
	double x = blocks[i]->GetDouble("X","mm");
	double y = blocks[i]->GetDouble("Y","mm");
	double z = blocks[i]->GetDouble("Z","mm");
	//z = z + 0.5*Length;
      	double Theta = blocks[i]->GetDouble("THETA","deg");
	double Phi= blocks[i]->GetDouble("PHI","deg");
	G4RotationMatrix* DetRot = new G4RotationMatrix();
	DetRot->rotateX(Theta);
	DetRot->rotateY(Phi);
	G4ThreeVector Pos(x,y,z);
	AddDetector(Pos,DetRot);
    }
    else{
      cout << "ERROR: Missing token for Khala blocks, check your input file" << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Khala::ConstructDetector(G4LogicalVolume* world){
  // Auxilliary elements
  double focus_po=0;
  G4Material* detector_mat = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4Material* Si_mat= MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Element* elN  = new G4Element("Nitrogen","N" ,  7.,14.01*g/mole);
  G4Material *N2= new G4Material("N2Gas",(1.251/1000)*g/cm3,1,kStateGas,300*kelvin,1.*atmosphere);
  N2->AddElement(elN,2);

  G4Element* C = new G4Element("Carbon", "C", 6., 12.011*g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", 1., 1.00794*g/mole);

  G4Material* plastic_mat = new G4Material("Plastic",  1.032*g/cm3, 2);
  plastic_mat->AddElement(C, 9);
  plastic_mat->AddElement(H, 10);

  G4Material* abs= new G4Material("abs",  1.05*g/cm3, 3);
  abs->AddElement(C, 15);
  abs->AddElement(H, 17);
  abs->AddElement(elN, 1);

  // 200-um thick Al chamber box
  G4Box* Chamber1=new G4Box("chamber1",100/2,100/2,180/2);
  G4Box* Chamber2=new G4Box("chamber2",100/2-0.2,100/2-0.2,180/2-0.2);
  G4SubtractionSolid* Chamber=new G4SubtractionSolid("Chamber",Chamber1,Chamber2);

  G4LogicalVolume* lv_chamber=
    new G4LogicalVolume(Chamber,detector_mat,"Chamber");
//  new G4PVPlacement(0,G4ThreeVector(focus_po,focus_po,focus_po),lv_chamber,"Chamber",world,false,0,true);

  G4VisAttributes* Att3= new G4VisAttributes(G4Colour(1.0,1.0,1.0,0.2));
  lv_chamber->SetVisAttributes(Att3);

  // space for LN2 gas volume
  G4Box* N2_box=new G4Box("N2_box",100/2-0.2,100/2-0.2,180/2-0.2);

  G4LogicalVolume* lv_N2= new G4LogicalVolume(N2_box,N2,"N2_gas");
  //  new G4PVPlacement(0,G4ThreeVector(focus_po,focus_po,focus_po),lv_N2,"N2_gas",world,false,0,true);


  G4VisAttributes* Att4= new G4VisAttributes(G4Colour(1.0,1.0,0.0,0.2));
  lv_N2->SetVisAttributes(Att4);


  // Placing KHALA detectors
  unsigned int mysize = m_Pos.size();
  for(unsigned int i = 0 ; i < mysize ; i++)
    new G4PVPlacement(G4Transform3D(*m_Rot[i], m_Pos[i]), ConstructDetector(),  "KhalaDetector", world, false, i); 

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Khala::ConstructDetector(){
  if(!m_LogicalDetector){

    G4Material* Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4Material* Alu = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* LaBr3 = MaterialManager::getInstance()->GetMaterialFromLibrary("LaBr3_Ce");

    // Mother Volume
    G4Tubs* solidKhalaDetector = 
      new G4Tubs("Khala",0, 0.5*FaceFront, 0.5*Length, 0.*deg, 360.*deg);
    m_LogicalDetector = 
      new G4LogicalVolume(solidKhalaDetector, Vacuum, "Khala", 0, 0, 0);

    m_LogicalDetector->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Detector construction
    // LaBr3
    G4ThreeVector  positionLaBr3 = G4ThreeVector(0, 0, LaBr3_PosZ);

    G4Tubs* solidLaBr3 = new G4Tubs("solidLaBr3", 0., 0.5*LaBr3Face, 0.5*LaBr3Thickness, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLaBr3 = new G4LogicalVolume(solidLaBr3, LaBr3, "logicLaBr3", 0, 0, 0);

    new G4PVPlacement(0, 
		      positionLaBr3, 
		      logicLaBr3, 
		      "Khala_LaBr3", 
		      m_LogicalDetector, 
		      false, 
		      0);


    // Set LaBr3 sensible
    logicLaBr3->SetSensitiveDetector(m_LaBr3Scorer);

    // Visualisation of LaBr3 Strip
    logicLaBr3->SetVisAttributes(m_LaBr3VisAtt);

    // Aluminium can around LaBr3
    // LaBr3 Can
    G4ThreeVector  positionLaBr3Can = G4ThreeVector(0, 0, LaBr3Can_PosZ);

    G4Tubs* solidLaBr3Can = new G4Tubs("solidLaBr3Can", 0.5*CanInnerDiameter, 0.5*CanOuterDiameter, 0.5*CanLength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLaBr3Can = new G4LogicalVolume(solidLaBr3Can, Alu, "logicLaBr3Can", 0, 0, 0);

    new G4PVPlacement(0, 
		      positionLaBr3Can, 
		      logicLaBr3Can, 
		      "Khala_LaBr3Can", 
		      m_LogicalDetector, 
		      false, 
		      0);

    // Visualisation of LaBr3Can
    //		logicLaBr3Can->SetVisAttributes(m_DetectorCasingVisAtt);
    logicLaBr3Can->SetVisAttributes(m_PMTVisAtt);

    // Aluminium window in front of LaBr3
    // LaBr3 Window
    G4ThreeVector  positionLaBr3Win = G4ThreeVector(0, 0, LaBr3Win_PosZ);

    G4Tubs* solidLaBr3Win = new G4Tubs("solidLaBr3Win", 0.5*WinInnerDiameter, 0.5*WinOuterDiameter, 0.5*WinLength, 0.*deg, 360.*deg);
    G4LogicalVolume* logicLaBr3Win = new G4LogicalVolume(solidLaBr3Win, Alu, "logicLaBr3Win", 0, 0, 0);

    new G4PVPlacement(0, 
		      positionLaBr3Win, 
		      logicLaBr3Win, 
		      "Khala_LaBr3Win", 
		      m_LogicalDetector, 
		      false, 
		      0);   

    // Visualisation of LaBr3Win
    //		logicLaBr3Win->SetVisAttributes(m_DetectorCasingVisAtt);
    logicLaBr3Win->SetVisAttributes(m_PMTVisAtt);


    // PMT
    G4ThreeVector  positionPMT = G4ThreeVector(0, 0, PMT_PosZ);

    G4Tubs* solidPMT = new G4Tubs("solidPMT", 0.5*PMTIn, 0.5*PMTOut, 0.5*PMTThickness, 0.*deg, 360.*deg);

    G4LogicalVolume* logicPMT = new G4LogicalVolume(solidPMT, Alu, "logicPMT", 0, 0, 0); 

    new G4PVPlacement(0, 
		      positionPMT, 
		      logicPMT, 
		      "Khala_PMT", 
		      m_LogicalDetector, 
		      false, 
		      0);

    // Visualisation of PMT Strip
    logicPMT->SetVisAttributes(m_PMTVisAtt);

    // PMT Cover
    G4ThreeVector  positionPMTCover = G4ThreeVector(0, 0, PMTCover_PosZ);

    G4Tubs* solidPMTCover = new G4Tubs("solidPMTCover", 0.5*PMTCoverInnerDiameter, 0.5*PMTCoverOuterDiameter, 0.5*PMTCoverLength, 0.*deg, 360.*deg);

    G4LogicalVolume* logicPMTCover = new G4LogicalVolume(solidPMTCover, Alu, "logicPMTCover", 0, 0, 0); 

    new G4PVPlacement(0, 
		      positionPMTCover, 
		      logicPMTCover, 
		      "Khala_PMTCover", 
		      m_LogicalDetector, 
		      false, 
		      0);

    // Visualisation of PMT Cover
    logicPMTCover->SetVisAttributes(m_PMTVisAtt);


    // PM Base
    G4ThreeVector  positionPM = G4ThreeVector(0, 0, PM_PosZ);

    G4Tubs* solidPM = new G4Tubs("solidPM", 0.5*PMIn, 0.5*PMOut, 0.5*PMThickness, 0.*deg, 360.*deg);

    G4LogicalVolume* logicPM = new G4LogicalVolume(solidPM, Alu, "logicPM", 0, 0, 0); 

    new G4PVPlacement(0, 
		      positionPM, 
		      logicPM, 
		      "Khala_PM", 
		      m_LogicalDetector, 
		      false, 
		      0);

    // Visualisation of PM Base
    logicPM->SetVisAttributes(m_PMTVisAtt);

    // PM Base Cover
    G4ThreeVector  positionPMCover = G4ThreeVector(0, 0, PMCover_PosZ);

    G4Tubs* solidPMCover = new G4Tubs("solidPMCover", 0.5*PMCoverInnerDiameter, 0.5*PMCoverOuterDiameter, 0.5*PMCoverLength, 0.*deg, 360.*deg);

    G4LogicalVolume* logicPMCover = new G4LogicalVolume(solidPMCover, Alu, "logicPMCover", 0, 0, 0); 

    new G4PVPlacement(0, 
		      positionPMCover, 
		      logicPMCover, 
		      "Khala_PMCover", 
		      m_LogicalDetector, 
		      false, 
		      0);

    // Visualisation of PM Base Cover
    logicPMCover->SetVisAttributes(m_PMTVisAtt);


  }

  return m_LogicalDetector;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Khala::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Khala")){
    pTree->Branch("Khala", "TKhalaData", &m_Event) ;
  } 
  pTree->SetBranchAddress("Khala", &m_Event) ;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Khala::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // LaBr3
  CalorimeterScorers::PS_Calorimeter* Scorer_labr3= (CalorimeterScorers::PS_Calorimeter*) m_LaBr3Scorer->GetPrimitive(0);
  unsigned int size_labr3 = Scorer_labr3->GetMult();

  // JP: simulate 138La decay beforehand
  // first get the decay time (large for beta decays!)
  double time_decay_sum = 0;
  int ntime_decay = 0;
  double time_for_bgd;
  for(unsigned int i = 0 ; i < size_labr3 ; i++){
    vector<unsigned int> level = Scorer_labr3->GetLevel(i); 
    time_decay_sum+=Scorer_labr3->GetTime(i);
    ntime_decay++;
  }
  
  if (ntime_decay > 0) // get the average time for decay occuring (may contain long-lived beta-decay time)
    time_for_bgd = time_decay_sum/ntime_decay;
  else 
    time_for_bgd = 0;
  if (time_for_bgd < DECAY_WINDOW)  // not larger than the decay window, likely isomeric decay - start from 0
    time_for_bgd = 0;
 

  // add 138La decay gammas
  for(unsigned int i = 0 ; i < NDET_KHALA ; i++){    
    if (RandFlat::shoot() < DECAY_PROB){ // we have a 138La decay
      double Energy = hist_138La->GetRandom();
      double Time = time_for_bgd+RandFlat::shoot()*DECAY_WINDOW; // add to the (beta) decay time
      m_Event->SetKhalaLaBr3E(i, Energy);
      m_Event->SetKhalaLaBr3T(i, Time);
      //      G4cout<<i<<" "<<Energy<<" "<<time_for_bgd<<" "<<Time/ns<<G4endl;
    }
  }
  // JP: ... then record the actual gamma events
  for(unsigned int i = 0 ; i < size_labr3 ; i++){
    vector<unsigned int> level = Scorer_labr3->GetLevel(i); 

    double E = Scorer_labr3->GetEnergy(i)*1000; 
    double sigma_energy = pow(E,-0.51)*E/2.355;
    double Energy = RandGauss::shoot(E,sigma_energy);

    if(Energy>EnergyThreshold){
      double sigma_time = (Jitter_A/sqrt((Energy + Jitter_B))+Jitter_C)/2.355;
      double Time = RandGauss::shoot(Scorer_labr3->GetTime(i), sigma_time)-DIST_KHALA/c_light; // correct for the time taken by gamma to travel to the detector
      int DetectorNbr =Scorer_labr3-> GetLevel(i)[0];
      // G4cout<<Time<<G4endl;
      // Revised geometry, Aug. 1, 2022
      // 9 detectors to be placed at 40 and 140 degree-rings with alternate placements, 17 cm from the center
      // if (DetectorNbr >= 36) // relative timing correction 
      // 	Time-=3*cm/c_light; // 3 cm further away for the 12 detectors at theta = 40 degrees
      
      m_Event->SetKhalaLaBr3E(DetectorNbr,Energy);
      m_Event->SetKhalaLaBr3T(DetectorNbr,Time);
    }
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Khala::InitializeScorers(){
  vector<G4int> NestingLevel;
  NestingLevel.push_back(1);

  bool already_exist_labr3 = false;
  
  m_LaBr3Scorer = CheckScorer("Khala_LaBr3Scorer",already_exist_labr3);

  // if the scorers were created previously nothing else need to be made
  if(already_exist_labr3) return;
  
   //   LaBr3 associated Scorer
  G4VPrimitiveScorer* LaBr3Scorer =
    new  CalorimeterScorers::PS_Calorimeter("KhalaLaBr3",NestingLevel);
  //and register it to the multifunctional detector
  m_LaBr3Scorer->RegisterPrimitive(LaBr3Scorer);

  
  //   Add the scorer to the Global Scorer Manager
  G4SDManager::GetSDMpointer()->AddNewDetector(m_LaBr3Scorer);

  // Read in 138La background spectrum file and histogram
  file_138La = new TFile("la138_bgd.root");
  hist_138La = (TH1D*)file_138La->Get("hbgd_138La");
  
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Khala::Construct(){
  return  (NPS::VDetector*) new Khala();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_khala{
  public:
    proxy_nps_khala(){
      NPS::DetectorFactory::getInstance()->AddToken("Khala","Khala");
      NPS::DetectorFactory::getInstance()->AddDetector("Khala",Khala::Construct);
    }
  };

  proxy_nps_khala p_nps_khala;
}
