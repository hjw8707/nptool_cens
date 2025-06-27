/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com        *
 *                                                                           *
 * Creation Date  : 11월 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  CsITl simulation                             *
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
#include "CsITl.hh"
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
namespace CsITl_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.*MeV;
  const double ResoTime = 4.5*ns ;
  double ResoEnergy = 0.1; // [%] in sigma
  const double Radius = 50*mm ; 
  const double Width = 100*mm ;
  const double Thickness = 300*mm ;
  const string Material = "CsITl";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// CsITl Specific Method
CsITl::CsITl(){
  m_Event = new TCsITlData() ;
  m_CsITlScorer = 0;

  m_matCsITl = 0;
  m_matLead = 0;

  // RGB Color + Transparency
  m_VisCsITl = new G4VisAttributes(G4Colour(0.7, 0.9, 0.7, 1.0));   

}

CsITl::~CsITl(){
}


void CsITl::DefineMaterials() {

  G4Element *Tl = new G4Element("Thallium","Tl",81.,204.383*g/mole );
  
  m_matCsITl = new G4Material("CsITl",3.6667*g/cm3, 2);
  m_matCsITl->AddMaterial(MaterialManager::getInstance()->GetMaterialFromLibrary("CsI"), 99.6*perCent);
  m_matCsITl->AddElement(Tl, 0.4*perCent);
  
  m_matLead = MaterialManager::getInstance()->GetMaterialFromLibrary("Pb");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CsITl::AddDetector(G4ThreeVector Pos, G4RotationMatrix Rot,
			 G4double Horizontal, G4double Vertical, G4double CsIThickness,
			 G4double LeadThickness) {
  m_Pos.push_back(Pos);
  m_Rot.push_back(Rot);
  m_CsIHorizontal.push_back(Horizontal);
  m_CsIVertical.push_back(Vertical);
  m_CsIThickness.push_back(CsIThickness);
  m_LeadThickness.push_back(LeadThickness);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* CsITl::BuildDetector(G4int i){

  G4Box* solidCsITl = new G4Box("CsITl_Box",
				m_CsIHorizontal[i]*0.5,
				m_CsIVertical[i]*0.5,
				m_CsIThickness[i]*0.5);

  G4LogicalVolume* logicCsITl = new G4LogicalVolume(solidCsITl, m_matCsITl, "logicCsITl", 0, 0, 0);
  logicCsITl->SetVisAttributes(m_VisCsITl);
  logicCsITl->SetSensitiveDetector(m_CsITlScorer);                      
  
  return logicCsITl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void CsITl::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("CsITl");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> reso = {"Reso"};
  
  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};
  vector<string> cyli = {"Rho","Phi","Z"};

  vector<string> cuboid= {"Horizontal","Vertical","Thickness","LeadThickness"};


  for(unsigned int i = 0 ; i < blocks.size() ; i++){

    ////////////////////////////////////////////////////////////
    // Resolution
    if (blocks[i]->HasTokenList(reso)) {
      CsITl_NS::ResoEnergy = blocks[i]->GetDouble("Reso","void");
      continue; }
    ////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////
    // checking the position items of the block
    G4ThreeVector Pos;
    if (blocks[i]->HasTokenList(cart)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CsITl " << i+1 <<  endl;
      Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm")); }
    else if (blocks[i]->HasTokenList(sphe)) {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CsITl " << i+1 <<  endl;
      Pos.setRThetaPhi(blocks[i]->GetDouble("R","mm"),
		       blocks[i]->GetDouble("Theta","deg"),
		       blocks[i]->GetDouble("Phi","deg")); }
    else if (blocks[i]->HasTokenList(cyli)) {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CsITl " << i+1 <<  endl;
      Pos.setRhoPhiZ(blocks[i]->GetDouble("Rho","mm"),
		     blocks[i]->GetDouble("Phi","deg"),
		     blocks[i]->GetDouble("Z","mm")); }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);}

    G4RotationMatrix Rot;
    if (blocks[i]->HasToken("ANG")) {
      G4ThreeVector Ang;
      Ang = NPS::ConvertVector(blocks[i]->GetTVector3("ANG","deg"));
      Rot.set(Ang.x(), Ang.y(), Ang.z()); } // (phi, theta, psi)
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // checking the shape items of the block
    if(blocks[i]->HasTokenList(cuboid)){
      double Horizontal = blocks[i]->GetDouble("Horizontal","mm");
      double Vertical = blocks[i]->GetDouble("Vertical","mm");
      double Thickness = blocks[i]->GetDouble("Thickness","mm");
      double LeadThickness = blocks[i]->GetDouble("LeadThickness","mm");
      AddDetector(Pos, Rot,
		  Horizontal, Vertical, Thickness,
		  LeadThickness); }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
    ////////////////////////////////////////////////////////////
    
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void CsITl::ConstructDetector(G4LogicalVolume* world){

  DefineMaterials();
  for (unsigned short i = 0 ; i < m_Pos.size() ; i++) {
    new G4PVPlacement(G4Transform3D(m_Rot[i],m_Pos[i]),
		      BuildDetector(i),
		      "CsITl",world,false,i+1);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void CsITl::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("CsITl")){
    pTree->Branch("CsITl", "TCsITlData", &m_Event) ;
  }
  pTree->SetBranchAddress("CsITl", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void CsITl::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_CsITlScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    //    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),CsITl_NS::ResoEnergy);
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i), 100.*CsITl_NS::ResoEnergy*(Scorer->GetEnergy(i)));
    if(Energy > CsITl_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),CsITl_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->Set(DetectorNbr,Energy,Time);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void CsITl::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_CsITlScorer = CheckScorer("CsITlScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_CsITlScorer->RegisterPrimitive(Calorimeter);
  m_CsITlScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CsITlScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* CsITl::Construct(){
  return  (NPS::VDetector*) new CsITl();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_CsITl{
    public:
      proxy_nps_CsITl(){
        NPS::DetectorFactory::getInstance()->AddToken("CsITl","CsITl");
        NPS::DetectorFactory::getInstance()->AddDetector("CsITl",CsITl::Construct);
      }
  };

  proxy_nps_CsITl p_nps_CsITl;
}
