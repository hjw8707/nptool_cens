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
 *  This class describe  TOGAXSI_GAGG simulation                             *
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
#include "TOGAXSI_GAGG.hh"
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
namespace TOGAXSI_GAGG_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 1.0*MeV ;

  //Size of GAGG crystal
  const double Crystal_Length = 3.5*cm ; 
  const double Crystal_Width = 3.5*cm ;
  const double Crystal_Height = 12*cm ;
}

using namespace TOGAXSI_GAGG_NS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// TOGAXSI_GAGG Specific Method
TOGAXSI_GAGG::TOGAXSI_GAGG(){
  InitializeMaterial();
  m_ReactionRegion = NULL;
  m_Event = new TTOGAXSI_GAGGData() ;

  m_RecoilArray = 0;
  m_ClusterArray = 0;

  m_RecoilArrayScorer = 0;
  m_ClusterArrayScorer = 0;

  // RGB Color + Transparency
  // Yellow
  m_VisGAGG = new G4VisAttributes(G4Colour(1, 1, 0, 0.5));   

  m_VisFrame = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));   

  //Transparent blue
  m_VisTarget = new G4VisAttributes(G4Colour(0.15,0.85,0.85,0.1));
}

TOGAXSI_GAGG::~TOGAXSI_GAGG(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_GAGG::AddGAGGRecoilArray(G4ThreeVector Pos, double Phi, G4ThreeVector Ref){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_Pos_RecoilArray.push_back(Pos);
  m_Phi_RecoilArray.push_back(Phi);
  m_Ref_RecoilArray.push_back(Ref);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TOGAXSI_GAGG::AddGAGGClusterArray(G4ThreeVector Pos, double Phi, G4ThreeVector Ref){
  m_Pos_ClusterArray.push_back(Pos);
  m_Phi_ClusterArray.push_back(Phi);
  m_Ref_ClusterArray.push_back(Ref);
}

void TOGAXSI_GAGG::AddTarget(double R, double L, string MaterialName, string CellMaterialName, double CellThickness, G4ThreeVector Pos) {

  m_Target_R.push_back(R);
  m_Target_L.push_back(L);
  m_Target_MaterialName.push_back(MaterialName);
  m_Target_CellMaterialName.push_back(CellMaterialName);
  m_Target_CellThickness.push_back(CellThickness);
  m_Target_Pos.push_back(Pos);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_GAGG::BuildGAGGRecoilArray(){
  if(!m_RecoilArray) {

    G4Box* RecoilArrayFull = new G4Box("RecoilArrayFull", Crystal_Length*0.5, Crystal_Width*0.5, Crystal_Height*0.5);

    //Master volume Recoil Array
    m_RecoilArray = new G4LogicalVolume(RecoilArrayFull, m_MaterialVacuum,"logicRecoilArray",0,0,0);
    m_RecoilArray->SetVisAttributes(G4VisAttributes::GetInvisible());

    //Crystals 
    G4Box* Crystal = new G4Box("GAGG_Crystal", Crystal_Length*0.5, Crystal_Width*0.5, Crystal_Height*0.5);

    G4LogicalVolume* logicCrystal = new G4LogicalVolume(Crystal,m_MaterialGAGG,"logic_Crystal",0,0,0);
    logicCrystal->SetVisAttributes(m_VisGAGG);
    logicCrystal->SetSensitiveDetector(m_RecoilArrayScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0), logicCrystal, "GAGG_RecoilArray_Crystal", m_RecoilArray, false, 1);

    //Active area needed?
/*
    G4Box* ActiveCrystal = new G4Box("RecoilArrayActiveCrystal", 0.5 * Crystal_Length, 0.5 * Crystal_Width, 0.5 * Crystal_Height);
    G4LogicalVolume* logicActiveCrystal = new G4LogicalVolume(ActiveCrystal, m_MaterialGAGG, "logicActiveCrystal", 0, 0, 0);
    logicActiveCrystal->SetVisAttributes(m_VisCrystal);
    ActiveCrystal->SetSensitiveDetector(m_RecoilArrayScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0), logicActiveCrystal, "RecoilArray_Active_Crystal", logicCrystal, false, 1);
*/
    }

  return m_RecoilArray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_GAGG::BuildGAGGClusterArray(){
  if(!m_ClusterArray) {

    G4Box* ClusterArrayFull = new G4Box("ClusterArrayFull", Crystal_Length*0.5, Crystal_Width*0.5, Crystal_Height*0.5);

    //Master volume Recoil Array
    m_ClusterArray = new G4LogicalVolume(ClusterArrayFull, m_MaterialVacuum,"logicClusterArray",0,0,0);
    m_ClusterArray->SetVisAttributes(G4VisAttributes::GetInvisible());

    //Crystals 
    G4Box* Crystal = new G4Box("GAGG_Crystal", Crystal_Length*0.5, Crystal_Width*0.5, Crystal_Height*0.5);

    G4LogicalVolume* logicCrystal = new G4LogicalVolume(Crystal,m_MaterialGAGG,"logic_Crystal",0,0,0);
    logicCrystal->SetVisAttributes(m_VisGAGG);
    logicCrystal->SetSensitiveDetector(m_ClusterArrayScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0), logicCrystal, "GAGG_ClusterArray_Crystal", m_ClusterArray, false, 1);

    //Active area needed?
/*
    G4Box* ActiveCrystal = new G4Box("RecoilArrayActiveCrystal", 0.5 * Crystal_Length, 0.5 * Crystal_Width, 0.5 * Crystal_Height);
    G4LogicalVolume* logicActiveCrystal = new G4LogicalVolume(ActiveCrystal, m_MaterialGAGG, "logicActiveCrystal", 0, 0, 0);
    logicActiveCrystal->SetVisAttributes(m_VisCrystal);
    logicActiveCrystal->SetSensitiveDetector(m_RecoilArrayScorer);

    new G4PVPlacement(new G4RotationMatrix(0,0,0), G4ThreeVector(0,0,0), logicActiveCrystal, "RecoilArray_Active_Crystal", logicCrystal, false, 1);
*/
    }

  return m_ClusterArray;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* TOGAXSI_GAGG::BuildTarget(int i) {

  if(i>0) {
    cout << "ERROR: Multiple TOGAXSI target block defined in detector file" << endl;
  }

  G4Material* TargetMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_Target_MaterialName[i]);
  G4Tubs* solidTarget = new G4Tubs("Target",0.,m_Target_R[i], m_Target_L[i] / 2., 0, 360.);
  m_Target = new G4LogicalVolume(solidTarget, TargetMaterial, "Target");
  m_Target->SetVisAttributes(m_VisTarget);

  return m_Target;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void TOGAXSI_GAGG::ReadConfiguration(NPL::InputParser parser){

  //GAGG Recoil
  vector<NPL::InputBlock*> blocks_recoil = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_GAGG", "RecoilArray");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_recoil.size() << " detectors found " << endl; 

  vector<string> gagg_recoilarray = {"Pos","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_recoil.size() ; i++){
    if(blocks_recoil[i]->HasTokenList(gagg_recoilarray)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_GAGG " << i+1 <<  endl;
            
      G4ThreeVector Pos = NPS::ConvertVector(blocks_recoil[i]->GetTVector3("Pos","mm"));
      double Phi = blocks_recoil[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_recoil[i]->GetTVector3("Ref","mm"));
      AddGAGGRecoilArray(Pos,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  //GAGG Cluster
  vector<NPL::InputBlock*> blocks_cluster = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_GAGG", "ClusterArray");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_cluster.size() << " detectors found " << endl; 

  vector<string> gagg_clusterarray = {"Pos","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_cluster.size() ; i++){
    if(blocks_cluster[i]->HasTokenList(gagg_clusterarray)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_GAGG " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks_cluster[i]->GetTVector3("Pos","mm"));
      double Phi = blocks_cluster[i]->GetDouble("Phi","deg");
      G4ThreeVector Ref = NPS::ConvertVector(blocks_cluster[i]->GetTVector3("Ref","mm"));
      AddGAGGClusterArray(Pos,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }



  //Target
  vector<NPL::InputBlock*> blocks_target = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_GAGG","Target");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_target.size() << " target found " << endl; 

  vector<string> targettoken = {"Radius","Length","TargetMaterial","CellMaterial","CellThickness","Pos"};

  for(unsigned int i = 0 ; i < blocks_target.size() ; i++){
    if(blocks_target[i]->HasTokenList(targettoken)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_GAGG " << i+1 <<  endl;

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
void TOGAXSI_GAGG::ConstructDetector(G4LogicalVolume* world){
  //RecoilArray	
  for (unsigned short i = 0; i < m_Phi_RecoilArray.size(); i++) {

    G4ThreeVector Det_pos = m_Pos_RecoilArray[i];

    Det_pos.rotate(-m_Phi_RecoilArray[i], G4ThreeVector(0,0,1));
    //    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_Phi_RecoilArray[i]);
    //    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_Ref_RecoilArray[i]), BuildGAGGRecoilArray(),"TOGAXSI_GAGG",world, false, i + 1);

  }  

  cout << "________Test______" << endl;
  //ClusterArray	
  for (unsigned short i = 0; i < m_Phi_ClusterArray.size(); i++) {

    G4ThreeVector Det_pos = m_Pos_ClusterArray[i];

    Det_pos.rotate(-m_Phi_ClusterArray[i], G4ThreeVector(0,0,1));
    //    Det_pos.rotate( 0, G4ThreeVector(0,0,1));
    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, m_Phi_ClusterArray[i]);
    //    G4RotationMatrix* Rot = new G4RotationMatrix(0*deg,0*deg, 60*deg);
    new G4PVPlacement(G4Transform3D(*Rot, Det_pos + m_Ref_ClusterArray[i]), BuildGAGGClusterArray(),"TOGAXSI_GAGG",world, false, i + 1);

  }  

  cout << "________Test______" << endl;
  // Target
  G4LogicalVolume* logicTarget[m_Target_R.size()];
  //  G4LogicalVolume* logicTargetCell[m_Target_R.size()];
 
  for (unsigned short i = 0; i < m_Target_R.size(); i++) {
    G4ThreeVector Tar_pos = m_Target_Pos[i];
    cout << "TargetPos" << m_Target_Pos[i].z() << endl;
    G4RotationMatrix* Rot = new G4RotationMatrix();
    logicTarget[i] = BuildTarget(i);
    new G4PVPlacement(Rot, Tar_pos, logicTarget[i], "TOGAXSI_SI_Target", world, false, i + 1);
    //    logicTargetCell[i] = BuildTargetCell(i);
    //    new G4PVPlacement(Rot, Tar_pos, logicTargetCell[i], "Strasse_TargetCell", world, false, i + 1);
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
void TOGAXSI_GAGG::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("TOGAXSI_GAGG")){
    pTree->Branch("TOGAXSI_GAGG", "TTOGAXSI_GAGGData", &m_Event) ;
  }
  pTree->SetBranchAddress("TOGAXSI_GAGG", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void TOGAXSI_GAGG::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer for RecoilArray 
  CalorimeterScorers::PS_Calorimeter* RecoilArrayScorer= (CalorimeterScorers::PS_Calorimeter*) m_RecoilArrayScorer->GetPrimitive(0);

  unsigned int size = RecoilArrayScorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = RecoilArrayScorer->GetLevel(i); 
    double Energy = RandGauss::shoot(RecoilArrayScorer->GetEnergy(i),ResoEnergy);
    if(Energy>EnergyThreshold){
      double Time = RandGauss::shoot(RecoilArrayScorer->GetTime(i),ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetRecoilEnergy(DetectorNbr,Energy);
      m_Event->SetRecoilTime(DetectorNbr,Time); 
    }
  }
  RecoilArrayScorer->clear(); 

  // Calorimeter scorer for RecoilArray 
  CalorimeterScorers::PS_Calorimeter* ClusterArrayScorer= (CalorimeterScorers::PS_Calorimeter*) m_ClusterArrayScorer->GetPrimitive(0);

  size = ClusterArrayScorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = ClusterArrayScorer->GetLevel(i); 
    double Energy = RandGauss::shoot(ClusterArrayScorer->GetEnergy(i),ResoEnergy);
    if(Energy>EnergyThreshold){
      double Time = RandGauss::shoot(ClusterArrayScorer->GetTime(i),ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetClusterEnergy(DetectorNbr,Energy);
      m_Event->SetClusterTime(DetectorNbr,Time); 
    }
  }
  ClusterArrayScorer->clear(); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void TOGAXSI_GAGG::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_RecoilArrayScorer = CheckScorer("RecoilArrayScorer",already_exist) ;
  m_ClusterArrayScorer = CheckScorer("ClusterArrayScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  // RecoilArray
  vector<int> level; level.push_back(1);
  G4VPrimitiveScorer* Calorimeter_RecoilArray= new CalorimeterScorers::PS_Calorimeter("Calorimeter_RecoilArray",level, 0);
  G4VPrimitiveScorer* Interaction_RecoilArray= new InteractionScorers::PS_Interactions("Interaction_RecoilArray",ms_InterCoord, 0);
  G4VPrimitiveScorer* Calorimeter_ClusterArray= new CalorimeterScorers::PS_Calorimeter("Calorimeter_ClusterArray",level, 0);
  G4VPrimitiveScorer* Interaction_ClusterArray= new InteractionScorers::PS_Interactions("Interaction_ClusterArray",ms_InterCoord, 0);

  //and register it to the multifunctionnal detector
  m_RecoilArrayScorer->RegisterPrimitive(Calorimeter_RecoilArray);
  m_RecoilArrayScorer->RegisterPrimitive(Interaction_RecoilArray);
  m_ClusterArrayScorer->RegisterPrimitive(Calorimeter_ClusterArray);
  m_ClusterArrayScorer->RegisterPrimitive(Interaction_ClusterArray);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_RecoilArrayScorer);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_RecoilArrayScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* TOGAXSI_GAGG::Construct(){
  return  (NPS::VDetector*) new TOGAXSI_GAGG();
}

void TOGAXSI_GAGG::InitializeMaterial() {

  G4Element* elementGd = new G4Element("Gadolinium", "Gd", 64., 157.25*g/mole);
  G4Element* elementAl = new G4Element("Aluminum", "Gd", 13., 26.982*g/mole);
  G4Element* elementGa = new G4Element("Gallium", "Ga", 31., 69.723*g/mole);
  G4Element* elementO = new G4Element("Oxygen", "O", 8., 15.999*g/mole);
  m_MaterialGAGG = new G4Material("Material GAGG", 6.63*g/cm3, 4);
  m_MaterialGAGG->AddElement(elementGd,3);
  m_MaterialGAGG->AddElement(elementAl,2);
  m_MaterialGAGG->AddElement(elementGa,3);
  m_MaterialGAGG->AddElement(elementO,12);


  m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
//  m_MaterialFrame = MaterialManager::getInstance()->GetMaterialFromLibrary("");
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_TOGAXSI_GAGG{
    public:
      proxy_nps_TOGAXSI_GAGG(){
        NPS::DetectorFactory::getInstance()->AddToken("TOGAXSI_GAGG","TOGAXSI_GAGG");
        NPS::DetectorFactory::getInstance()->AddDetector("TOGAXSI_GAGG",TOGAXSI_GAGG::Construct);
      }
  };

  proxy_nps_TOGAXSI_GAGG p_nps_TOGAXSI_GAGG;
}
