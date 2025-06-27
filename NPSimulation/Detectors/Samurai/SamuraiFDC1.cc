/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Omar Nasr  contact address: omar.h.nasr@outlook.com      *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    : septembre 2021                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Samurai simulation                                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
//G4 Geometry object
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
#include "G4RegionStore.hh"

// NPTool header
#include "SamuraiFDC1.hh"
//#include "CalorimeterScorers.hh"
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
namespace SamuraiFDC1_NS{
  // Samurai magnet construction paramethers

  //Main outer box
  const double FDC1_Width = 1000*mm; //(x)
  const double FDC1_Height = 696*mm;//(y)
  const double FDC1_Depth = 336*mm;//(z)
  const string FDC1_Material = "G4_Galactic"; 

  //Detector Number
  const short int FDC1_DetectorNumber = 1;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Samurai Specific Method
SamuraiFDC1::SamuraiFDC1(){

  //Visualization attributes
  m_VisFDC1 = new G4VisAttributes(G4Colour(1,1,0,0.5));
  //Logical volumes
  m_FDC1 = NULL;
  //Scorer
  m_FDC1Scorer = NULL;
  //Ideal Data event
  m_Event = new TSamuraiIdealData;

}

SamuraiFDC1::~SamuraiFDC1(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC1::AddDetector(G4ThreeVector Mag_Pos, G4ThreeVector Offset){

  m_Pos = Mag_Pos + Offset;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* SamuraiFDC1::BuildFDC1(){
  if(!m_FDC1){
    //Shape - G4Box
    G4Box* box = new G4Box("FDC1_Box",SamuraiFDC1_NS::FDC1_Width*0.5,
			   SamuraiFDC1_NS::FDC1_Height*0.5,SamuraiFDC1_NS::FDC1_Depth*0.5);
  
    //Material - vacuum
    G4Material* VacuumMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(SamuraiFDC1_NS::FDC1_Material);

    //Logical Volume
    m_FDC1 = new G4LogicalVolume(box, VacuumMaterial, "logic_SamuraiFDC1_box",0,0,0);
    m_FDC1->SetVisAttributes(m_VisFDC1);
    m_FDC1->SetSensitiveDetector(m_FDC1Scorer);
  }
  return m_FDC1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// (Called in DetecorConstruction::ReadDetectorConfiguration Method)
void SamuraiFDC1::ReadConfiguration(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Samurai");
  vector<NPL::InputBlock*> blocks2 = parser.GetAllBlocksWithToken("SAMURAIFDC1");

  if(blocks.size()==1 && blocks2.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "/////// Samurai FDC1 found ///////" << endl;
    }
    vector<string> cart = {"POS","ANGLE"};
    vector<string> sphe = {"R","Theta","Phi","ANGLE"};

    G4ThreeVector Mag_Pos;

    if(blocks[0]->HasTokenList(cart)){
      Mag_Pos = NPS::ConvertVector(blocks[0]->GetTVector3("POS", "cm"));
    }
    else if(blocks[0]->HasTokenList(sphe)){
      double R = blocks[0]->GetDouble("R","mm");
      double Theta = blocks[0]->GetDouble("Theta","deg");
      double Phi = blocks[0]->GetDouble("Phi","deg");
      Mag_Pos.setMag(R);
      Mag_Pos.setTheta(Theta);
      Mag_Pos.setPhi(Phi);
    }
    
    G4ThreeVector Offset = NPS::ConvertVector(blocks2[0]->GetTVector3("Offset", "mm"));

    //string xml = blocks2[0]->GetString("XML");
    //bool invert_x = blocks2[0]->GetBool("InvertX");
    //bool invert_y = blocks2[0]->GetBool("InvertY");
    //bool invert_z = blocks2[0]->GetBool("InvertD");
  
    AddDetector(Mag_Pos, Offset);

  }
  else{
    cout << "ERROR: there should be only one Samurai magnet, check your input file" << endl;
    exit(1);
  }
    
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// (Called After DetectorConstruction::AddDetector Method)
void SamuraiFDC1::ConstructDetector(G4LogicalVolume* world){

  new G4PVPlacement(0, m_Pos,
          BuildFDC1(), "SamuraiFDC1", world, false, 0);

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SamuraiFDC1::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("IdealDataFDC1")){
    pTree->Branch("IdealDataFDC1", "TSamuraiIdealData", &m_Event) ;
  }
  pTree->SetBranchAddress("IdealDataFDC1", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Read sensitive part and fill the Root tree.
// (Called at in the EventAction::EndOfEventAvtion)

void SamuraiFDC1::ReadSensitive(const G4Event* event){
  
  m_Event->Clear();
  //Interaction Scorer
  InteractionScorers::PS_Interactions* Scorer= (InteractionScorers::PS_Interactions*) m_FDC1Scorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    //vector<unsigned int> level = Scorer->GetLevel(i); 
    //double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SamuraiFDC1_NS::ResoEnergy);
    //double Energy = Scorer->GetEnergy(i);
    
    double energy = Scorer->GetEnergy(i);
    double brho = Scorer->GetBrho(i);
    double posx = Scorer->GetPositionX(i);
    double posy = Scorer->GetPositionY(i);
    double posz = Scorer->GetPositionZ(i);
    double mom_mag = brho*Scorer->GetCharge(i);
    double theta = Scorer->GetTheta(i);
    double phi = Scorer->GetPhi(i);
    m_Event->SetData(SamuraiFDC1_NS::FDC1_DetectorNumber, energy, posx, posy, posz, mom_mag, theta, phi, brho);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  

void SamuraiFDC1::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_FDC1Scorer = CheckScorer("FDC1Scorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  //Interaction 
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ; 
  m_FDC1Scorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_FDC1Scorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////

NPS::VDetector* SamuraiFDC1::Construct(){
  return  (NPS::VDetector*) new SamuraiFDC1();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_samuraiFDC1{
    public:
      proxy_nps_samuraiFDC1(){
        NPS::DetectorFactory::getInstance()->AddToken("SAMURAIFDC1","SAMURAIFDC1");
        NPS::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC1",SamuraiFDC1::Construct);
      }
  };
  
  proxy_nps_samuraiFDC1 p_nps_samuraiFDC1;
}


