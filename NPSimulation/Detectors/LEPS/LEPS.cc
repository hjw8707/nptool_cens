/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : March 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  LEPS simulation                             *
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
#include "G4Polycone.hh"

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
#include "LEPS.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace LEPS_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.01*MeV;
  const double ResoTime = 0*ns ;
  // const double ResoEnergy = 0.34*keV ; // 0.8 keV FWHM, quite optimist
  const double ResoEnergy = 1e-9*keV ; 
  const double Radius = 50*mm ; 
  const double Width = 100*mm ;
  const double Thickness = 300*mm ;
  const string Material_Shell = "Al";
  const double radii_internal[16] = {50/2*mm-1*mm, 50/2*mm-1*mm, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  const double radii_external[16] = {50/2*mm, 50/2.*mm, 52/2.*mm, 52/2.*mm, 33.7/2.*mm, 33.7/2.*mm, 95/2.*mm, 95/2.*mm, 33.7/2.*mm, 33.7/2.*mm, 95/2.*mm, 95/2.*mm, 33/2.*mm, 33/2.*mm, 222/2.*mm, 222/2.*mm};
  const double length_external[16] = {0, 85*mm, 85*mm, 157*mm, 157*mm, 207*mm, 207*mm, 208*mm, 208*mm, 298*mm, 298*mm, 393*mm, 393*mm, 453*mm, 453*mm, 720*mm};
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// LEPS Specific Method
LEPS::LEPS(){
  m_Event = new TLEPSData() ;
  m_LEPSScorer = 0;
}

LEPS::~LEPS(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LEPS::AddDetector(G4ThreeVector POS, string  Shape){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void LEPS::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}


G4LogicalVolume* LEPS::BuildDetector(G4int DetNumber, G4ThreeVector Det_pos, G4RotationMatrix* Det_rot, G4LogicalVolume* world){
  // Materials
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(LEPS_NS::Material_Shell);
  G4Material* m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* m_MaterialCarbon = MaterialManager::getInstance()->GetMaterialFromLibrary("C");
  G4Material* m_MaterialGermanium = MaterialManager::getInstance()->GetMaterialFromLibrary("Germanium");

  G4VisAttributes* light_GreyAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.7));
  G4VisAttributes* RedAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.6));
  G4VisAttributes* GreenAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.6));

  // Global volume
  // Origin axis at the front of detector
  G4Tubs* solidLEPS = new G4Tubs("solidLEPS", 0, 25, 720*0.5*mm, 0, 360*deg);
  G4LogicalVolume* logicLEPS = new G4LogicalVolume(solidLEPS, m_MaterialVacuum, "logicLEPS", 0, 0);
  new G4PVPlacement(G4Transform3D(*Det_rot, Det_pos), logicLEPS, "LEPS", world, false, DetNumber);
  logicLEPS->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Enveloppe
  G4Polycone* LEPS_Cyl = new G4Polycone("LEPS_Cyl", 0, 360*deg, 16, LEPS_NS::length_external, LEPS_NS::radii_internal, LEPS_NS::radii_external);
  G4LogicalVolume* vol_LEPS = new G4LogicalVolume(LEPS_Cyl, DetectorMaterial, "logic_LEPS_Cyl", 0, 0, 0);
  G4ThreeVector LEPS_cyl_Pos = G4ThreeVector(0, 0, 0);
  new G4PVPlacement(0, LEPS_cyl_Pos, vol_LEPS, "LEPS_cyl", logicLEPS, false, DetNumber);
  vol_LEPS->SetVisAttributes(light_GreyAtt);

  // Carbon window
  G4Tubs* LEPS_CWindow = new G4Tubs("LEPS_CWindow", 0, 24, 0.6*0.5*mm, 0, 360*deg);
  G4LogicalVolume* vol_CWindow = new G4LogicalVolume(LEPS_CWindow, m_MaterialCarbon, "logic_LEPS_Window", 0, 0, 0);
  G4ThreeVector CWindow_Pos = G4ThreeVector(0, 0, 0);
  // G4ThreeVector CWindow_Pos = G4ThreeVector(0, 0, 720*0.5*mm);
  new G4PVPlacement(0, CWindow_Pos, vol_CWindow, "LEPS_CWindow", logicLEPS, false, DetNumber);
  vol_CWindow->SetVisAttributes(GreenAtt);
  
  // Germanium crystal
  G4Tubs* LEPS_crys = new G4Tubs("LEPS_crys",0 , 13*mm, 10.2*0.5*mm, 0, 360*deg);
  G4LogicalVolume* vol_crys = new G4LogicalVolume(LEPS_crys, m_MaterialGermanium, "logic_LEPS_crys", 0, 0, 0);
  G4ThreeVector crys_Pos = G4ThreeVector(0, 0, 0.6*1.5*mm+11*0.5*mm);
  new G4PVPlacement(0, crys_Pos, vol_crys, "LEPS_crys", logicLEPS, false, DetNumber);
  vol_crys->SetVisAttributes(RedAtt);

  vol_crys->SetSensitiveDetector(m_LEPSScorer);
  
  return vol_LEPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void LEPS::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("LEPS");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  LEPS " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  LEPS " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,Shape);
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
void LEPS::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    // Det_pos+=Det_pos.unit()*LEPS_NS::Thickness*0.5;
    // Building Detector reference frame
    G4double ii = cos(m_Theta[i]) * cos(m_Phi[i]);
    G4double jj = cos(m_Theta[i]) * sin(m_Phi[i]);
    G4double kk = -sin(m_Theta[i]);
    G4ThreeVector Y(ii,jj,kk);
    G4ThreeVector w = Det_pos.unit();
    G4ThreeVector u = w.cross(Y);
    G4ThreeVector v = w.cross(u);
    v = v.unit();
    u = u.unit();

    G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);

    BuildDetector(i+1, Det_pos, Rot, world);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void LEPS::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("LEPS")){
    pTree->Branch("LEPS", "TLEPSData", &m_Event) ;
  }
  pTree->SetBranchAddress("LEPS", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void LEPS::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_LEPSScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),LEPS_NS::ResoEnergy);
    if(Energy>LEPS_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),LEPS_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void LEPS::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_LEPSScorer = CheckScorer("LEPSScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_LEPSScorer->RegisterPrimitive(Calorimeter);
  m_LEPSScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_LEPSScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* LEPS::Construct(){
  return  (NPS::VDetector*) new LEPS();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_LEPS{
    public:
      proxy_nps_LEPS(){
        NPS::DetectorFactory::getInstance()->AddToken("LEPS","LEPS");
        NPS::DetectorFactory::getInstance()->AddDetector("LEPS",LEPS::Construct);
      }
  };

  proxy_nps_LEPS p_nps_LEPS;
}
