/******************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project          *
 *                                                                            *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                 *
 * For the list of contributors see $NPTOOL/Licence/Contributors              *
 ******************************************************************************/

/******************************************************************************
 * Original Author: Pierre Morfouace contact address: pierre.morfouace2@cea.fr*
 *                                                                            *
 * Creation Date  : November 2020                                             *
 * Last update    :                                                           *
 *----------------------------------------------------------------------------*
 * Decription:                                                                *
 *  This class describe a simple SofTofW setup for simulation                   *
 *                                                                            *
 *----------------------------------------------------------------------------*
 * Comment:                                                                   *
 *                                                                            *
 ******************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"

// NPTool header
#include "SofTofW.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

//CAD Mesh
#include "CADMesh.hh"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SofTofW_NS{
  // Energy and time Resolution
  const double EnergyThreshold       = 0.1*MeV;
  const double ResoTime              = 0.007*ns;
  const double ResoEnergy            = 1.0*MeV;
  
  const double tof_plastic_height    = 660*mm;
  const double tof_plastic_width     = 32*mm;
  const double tof_plastic_thickness = 0.5*mm;
  const string Material              = "BC400";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// SofTofW Specific Method
SofTofW::SofTofW(){
  m_Event = new TSofTofWData() ;
  m_TofScorer = 0;
  m_PlasticTof = 0;
  m_VacuumPipe= 0;
  m_TofWall = 0;

  m_VacuumPipeX= 0;
  m_VacuumPipeY= 0;
  m_VacuumPipeZ= 0;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0.53, 0.81, 0.98, 1));   
  m_VisKapton = new G4VisAttributes(G4Colour(1, 0.4, 0., 0.6));   

}

SofTofW::~SofTofW(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SofTofW::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SofTofW::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* SofTofW::BuildTOFDetector(){
  m_TofWall = new G4AssemblyVolume();

  if(!m_PlasticTof){
    G4Box* box = new G4Box("SofTofW_Box",SofTofW_NS::tof_plastic_height*0.5,
        SofTofW_NS::tof_plastic_width*0.5,SofTofW_NS::tof_plastic_thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(SofTofW_NS::Material);
    m_PlasticTof = new G4LogicalVolume(box,DetectorMaterial,"logic_SofTofW_Box",0,0,0);
    m_PlasticTof->SetVisAttributes(m_VisSquare);
    m_PlasticTof->SetSensitiveDetector(m_TofScorer);

    G4RotationMatrix* Rv = new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0);
    Tv.setY(0);
    Tv.setZ(0);
    for(unsigned int i=0; i<28; i++){
      int k = -14+i;
      Tv.setY(k*(SofTofW_NS::tof_plastic_width+0.1));
      m_TofWall->AddPlacedVolume(m_PlasticTof, Tv, Rv);
    }

  }
  return m_TofWall;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* SofTofW::BuildVacuumPipe(){
  if(!m_VacuumPipe){
    m_VacuumPipe = new G4AssemblyVolume;

    G4Tubs* tube = new G4Tubs("tube",8.5*cm,9.*cm,160./2*cm,0,360*deg);
    G4Tubs* Mylartube = new G4Tubs("Mylartube",0*cm,9.*cm,75./2*um,0,360*deg);
    
    G4Material* tube_mat = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    G4Material* mylar_mat = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");

    G4LogicalVolume* tube_vol = new G4LogicalVolume(tube,tube_mat,"logic_tube",0,0,0);
    G4LogicalVolume* LogicMylar = new G4LogicalVolume(Mylartube,mylar_mat,"logic_mylar",0,0,0);

    G4VisAttributes* VisTube = new G4VisAttributes(G4Colour(0., 0.7, 0.7));   
    tube_vol->SetVisAttributes(VisTube);

    LogicMylar->SetVisAttributes(m_VisKapton);

    G4ThreeVector Pos = G4ThreeVector(0,0,0);
    G4RotationMatrix* Rot = new G4RotationMatrix();
    m_VacuumPipe->AddPlacedVolume(tube_vol,Pos,Rot);

    // Entrance window
    Pos.setZ(-160./2*cm);
    m_VacuumPipe->AddPlacedVolume(LogicMylar,Pos,Rot);
    // Exit window
    Pos.setZ(160./2*cm);
    m_VacuumPipe->AddPlacedVolume(LogicMylar,Pos,Rot);
  }
  
  return m_VacuumPipe;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void SofTofW::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SofTofW");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Build_VacuumPipe"};
  vector<string> sphe = {"R","Theta","Phi","Build_VacuumPipe"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTofW " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      m_Build_VacuumPipe = blocks[i]->GetInt("Build_VacuumPipe");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SofTofW " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      m_Build_VacuumPipe = blocks[i]->GetInt("Build_VacuumPipe");
      m_VacuumPipeX = blocks[i]->GetDouble("VacuumPipeX","m");
      m_VacuumPipeY = blocks[i]->GetDouble("VacuumPipeY","m");
      m_VacuumPipeZ = blocks[i]->GetDouble("VacuumPipeZ","m");

      AddDetector(R,Theta,Phi);
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
void SofTofW::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*SofTofW_NS::tof_plastic_thickness*0.5;
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
    BuildTOFDetector()->MakeImprint(world,Det_pos,Rot);
  }


 
  if(m_Build_VacuumPipe==1){
    G4ThreeVector Tube_Pos = G4ThreeVector(m_VacuumPipeX,m_VacuumPipeY,m_VacuumPipeZ);
    BuildVacuumPipe()->MakeImprint(world,Tube_Pos,0);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SofTofW::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("SofTofW")){
    pTree->Branch("SofTofW", "TSofTofWData", &m_Event) ;
  }
  pTree->SetBranchAddress("SofTofW", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void SofTofW::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // TOF scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_TofScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),SofTofW_NS::ResoEnergy);
    if(Energy>SofTofW_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),SofTofW_NS::ResoTime);
      //int DetectorNbr = level[0];
      int PlasticNbr = level[1]-1;
      //m_Event->SetDetectorNbr(DetectorNbr);
      m_Event->SetPlasticNbr(PlasticNbr);
      m_Event->SetEnergy(Energy);
      m_Event->SetCoarseTime(Time); 
    }
  }
  Scorer->clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void SofTofW::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_TofScorer = CheckScorer("TofScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; 
  level.push_back(1);
  level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_TofScorer->RegisterPrimitive(Calorimeter);
  m_TofScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_TofScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SofTofW::Construct(){
  return  (NPS::VDetector*) new SofTofW();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_SofTofW{
    public:
      proxy_nps_SofTofW(){
        NPS::DetectorFactory::getInstance()->AddToken("SofTofW","SofTofW");
        NPS::DetectorFactory::getInstance()->AddDetector("SofTofW",SofTofW::Construct);
      }
  };

  proxy_nps_SofTofW p_nps_SofTofW;
}
