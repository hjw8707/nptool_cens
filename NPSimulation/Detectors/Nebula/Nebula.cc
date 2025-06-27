/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Nebula simulation                                   *
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
#include "G4MaterialPropertiesTable.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// NPTool header
#include "Nebula.hh"
#include "PlasticBar.hh"
#include "InteractionScorers.hh"
#include "ProcessScorers.hh"
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
namespace Nebula_NS{
  // Energy and time Resolution
  const double LightThreshold  = 0.1*MeV;
  const double ResoTime        = 0.75/2.355*ns; //0.75
  const double ResoEnergy      = 0.1/2.355*MeV; 
  const double ResoLight       = 0.1/2.355*MeV; 
  const double ResoPosition    = 1.0*um; //1.0
  const double ModuleWidth     = 120*mm ;
  const double ModuleLength    = 120*mm ;
  const double ModuleHeight    = 1800*mm ;
  const double InterModule     = 1*mm ;
  const double VetoWidth       = 320*mm ;
  const double VetoLength      = 10*mm ;
  const double VetoHeight      = 1900*mm ;
  const double InterVeto       = 1*mm ;
  const int    VetoPerWall     = 12;
  const int    VetoPerExpand   = 6;
  const double WallToVeto      = 10*cm;
  const double MaterialIndex   = 1.58;
  const double Attenuation     = 6680*mm; 

  const string Material = "BC400";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Nebula Specific Method
Nebula::Nebula(){
  m_Event = new TNebulaData() ;
  m_ModuleScorer = 0;
  m_VetoScorer = 0;
  m_Module = 0;
  m_Veto = 0;

  // RGB Color + Transparency
  m_VisModule = new G4VisAttributes(G4Colour(0.263, 0.682, 0.639, 1));   
  //m_VisModule = new G4VisAttributes(G4Colour(0.145, 0.384, 0.596, 1));   
  m_VisVeto   = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.8));   
  m_VisPMT    = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 1));   
  m_VisFrame  = new G4VisAttributes(G4Colour(0, 0.3, 1, 0.5));   

}

Nebula::~Nebula(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Nebula::AddWall(G4ThreeVector Pos, int NbrModule, bool Veto, bool Frame){
  // Convert the Pos value to R theta Phi as Spherical coordinate is easier in G4 
  m_Pos.push_back(Pos);
  m_NbrModule.push_back(NbrModule);
  m_HasVeto.push_back(Veto);
  m_HasFrame.push_back(Frame);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Nebula::BuildModule(){
  if(!m_Module){
    G4Box* box = new G4Box("Nebula_Module",Nebula_NS::ModuleWidth*0.5,
        Nebula_NS::ModuleHeight*0.5,Nebula_NS::ModuleLength*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Nebula_NS::Material);
    m_Module = new G4LogicalVolume(box,DetectorMaterial,"logic_Nebula_Module",0,0,0);
    m_Module->SetVisAttributes(m_VisModule);
    m_Module->SetSensitiveDetector(m_ModuleScorer);
  }
  return m_Module;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Nebula::BuildVeto(){
  if(!m_Veto){
    G4Box* box = new G4Box("Nebula_Veto",Nebula_NS::VetoWidth*0.5,
        Nebula_NS::VetoHeight*0.5,Nebula_NS::VetoLength*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Nebula_NS::Material);
    m_Veto = new G4LogicalVolume(box,DetectorMaterial,"logic_Nebula_Veto",0,0,0);
    m_Veto->SetVisAttributes(m_VisVeto);
    m_Veto->SetSensitiveDetector(m_VetoScorer);
  }
  return m_Veto;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Nebula::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEBULA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  // define an entire wall
  vector<string> wall = {"Pos","NumberOfModule","Veto","Frame"};

  // use an experiment xml file to position bars individually
  vector<string> xml= {"XML","Offset","InvertX","InvertY"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(wall)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Nebula " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("Pos","mm"));
      int NbrModule = blocks[i]->GetInt("NumberOfModule");
      bool Veto = blocks[i]->GetInt("Veto");
      bool Frame= blocks[i]->GetInt("Frame");
      AddWall(Pos,NbrModule,Veto,Frame);
    }
    else if(blocks[i]->HasTokenList(xml)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Nebula XML file " << i+1 <<  endl;
      std::string xml_file = blocks[i]->GetString("XML"); 
      G4ThreeVector Offset = NPS::ConvertVector(blocks[i]->GetTVector3("Offset","mm"));
      bool InvertX = blocks[i]->GetInt("InvertX"); 
      bool InvertY = blocks[i]->GetInt("InvertY"); 
      ReadXML(xml_file,Offset,InvertX,InvertY);
    }
      else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  std::for_each(m_NbrModule.begin(), m_NbrModule.end(), [&] (int n) {
     m_TotalModule += n;
  });  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Nebula::ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY){ 
  NPL::XmlParser xml;
  xml.LoadFile(xml_file);
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("NEBULA");  
  int NumberOfBars=0;
  for(unsigned int i = 0 ; i < b.size() ; i++){
    NumberOfBars++;
    unsigned int id = b[i]->AsInt("ID");
    // position
    auto PositionX = b[i]->AsDouble("PosX"); 
    auto PositionY = b[i]->AsDouble("PosY"); 
    auto PositionZ = b[i]->AsDouble("PosZ"); 
    // SubLayer 0 is use for Veto
    auto SubLayer  = b[i]->AsInt("SubLayer");
    // Name "NoUseX" is used to silence bars
    auto nousestr  = b[i]->AsString("NAME");
    // Remove unused bar
    if(nousestr.find("NoUse")==std::string::npos && PositionX!=PositionZ){
      if(InvertX)
        PositionX*=-1;
      if(InvertY)
        PositionY*=-1;
      m_PositionBar[id]= G4ThreeVector(PositionX,PositionY,PositionZ)+offset;
      if(SubLayer)
        m_IsVetoBar[id]= false;
      else
        m_IsVetoBar[id]=true;
    }
  } 
  cout << " -> " << NumberOfBars << " bars found" << endl;
} 

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Nebula::ConstructDetector(G4LogicalVolume* world){
  // Start with XML case
  G4RotationMatrix* Rot = new G4RotationMatrix();
  for(auto pos : m_PositionBar){
    if(!m_IsVetoBar[pos.first]){
      new G4PVPlacement(G4Transform3D(*Rot,pos.second),
            BuildModule(),
            "NebulaModule",world,false,pos.first);
    }
    else{
      new G4PVPlacement(G4Transform3D(*Rot,pos.second),
            BuildVeto(),
            "NebulaModule",world,false,pos.first);
      }
    }

  unsigned int nbrM = 1 ;
  unsigned int nbrV = 1 ;
  for (unsigned short i = 0 ; i < m_Pos.size() ; i++) {
    for (unsigned short m = 0 ; m < m_NbrModule[i] ; m++) {
      double offset = (Nebula_NS::ModuleWidth+Nebula_NS::InterModule)*(-m_NbrModule[i]*0.5+m)+Nebula_NS::ModuleWidth*0.5;
      G4ThreeVector Offset(offset,0,0);
      new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
          BuildModule(),
          "NebulaModule",world,false,nbrM++);
    }

    if(m_HasVeto[i]){
      if(m_NbrModule[i] > 15){
        for (unsigned short m = 0 ; m < Nebula_NS::VetoPerWall ; m++) {
          double offset = (Nebula_NS::VetoWidth+Nebula_NS::InterVeto)*(-Nebula_NS::VetoPerWall*0.5+m)+Nebula_NS::VetoWidth*0.5;
          G4ThreeVector Offset(offset,0,-Nebula_NS::WallToVeto);
          new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
            BuildVeto(),
            "NebulaVeto",world,false,nbrV++);
        }
      }
      else{
        for (unsigned short m = 0 ; m < Nebula_NS::VetoPerExpand ; m++) {
          double offset = (Nebula_NS::VetoWidth+Nebula_NS::InterVeto)*(-Nebula_NS::VetoPerExpand*0.5+m)+Nebula_NS::VetoWidth*0.5;
          G4ThreeVector Offset(offset,0,-Nebula_NS::WallToVeto);
          new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
            BuildVeto(),
            "NebulaVeto",world,false,nbrV++);
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetectorConstruction::AddDetector Method
void Nebula::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Nebula")){
    pTree->Branch("Nebula", "TNebulaData",&m_Event) ;
  }
  pTree->SetBranchAddress("Nebula", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Nebula::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // PlasticBar scorer
  PlasticBar::PS_PlasticBar* PlasticScorer_Module = (PlasticBar::PS_PlasticBar*) m_ModuleScorer->GetPrimitive(0);
  PlasticBar::PS_PlasticBar* PlasticScorer_Veto = (PlasticBar::PS_PlasticBar*) m_VetoScorer->GetPrimitive(0);
  // Should we put a ProcessScorer here to get the info if the particle is first neutron and give it to NebulaData ?
  
  double Time_up, Time_down;
  double Energy_tmp, Light_tmp;

  //////////// TRIAL TO GET THE OPTICAL INDEX FROM MATERIAL PROPERTIES /////////////
  //Trying to get Optical Index from Material directly
  //const G4Material* aMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Nebula_NS::Material);
  //G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  //if(!aMaterialPropertiesTable->PropertyExists("RINDEX")){
  //  MaterialIndex = !aMaterialPropertiesTable->GetConstProperty("RINDEX"); 
  //}
  //else{
  //  MaterialIndex = 0; 
  //}
  //cout << MaterialManager::getInstance()->GetMaterialFromLibrary(Nebula_NS::Material)->GetMaterialPropertiesTable()->GetMaterialPropertyNames()[0] << endl;
  //////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////// MODULE SCORER //////////////////////////////////
  unsigned int ModuleHits_size = PlasticScorer_Module->GetMult(); 
  for(unsigned int i = 0 ; i < ModuleHits_size ; i++){
    vector<unsigned int> level = PlasticScorer_Module->GetLevel(i); 
    Energy_tmp = PlasticScorer_Module->GetEnergy(i);
    Light_tmp = PlasticScorer_Module->GetLight(i);
    Energy = RandGauss::shoot(Energy_tmp, Energy_tmp*Nebula_NS::ResoEnergy);
    Light = RandGauss::shoot(Light_tmp, Light_tmp*Nebula_NS::ResoLight);
    
    if(Light>Nebula_NS::LightThreshold){
      int DetectorNbr = level[0];
      double Position = RandGauss::shoot(PlasticScorer_Module->GetPositionY(i),Nebula_NS::ResoPosition);
      
      m_Event->SetChargeUp(DetectorNbr,Light*exp(-(Nebula_NS::ModuleHeight/2-Position)/Nebula_NS::Attenuation));
      m_Event->SetChargeDown(DetectorNbr,Light*exp(-(Nebula_NS::ModuleHeight/2+Position)/Nebula_NS::Attenuation));
      
      // Take Time and Y Position and compute Tup and Tdown
      double Time = RandGauss::shoot(PlasticScorer_Module->GetTime(i),Nebula_NS::ResoTime);

      Time_up = (Nebula_NS::ModuleHeight/2-Position)/(c_light/Nebula_NS::MaterialIndex) + Time;
      m_Event->SetTimeUp(DetectorNbr,Time_up);
      
      Time_down = (Nebula_NS::ModuleHeight/2+Position)/(c_light/Nebula_NS::MaterialIndex) + Time;
      m_Event->SetTimeDown(DetectorNbr,Time_down);
    }
  }
  //cout << endl;

  ///////////////////////////////// VETO SCORER //////////////////////////////////
  unsigned int VetoHits_size = PlasticScorer_Veto->GetMult();
  for(unsigned int i = 0 ; i < VetoHits_size ; i++){
    vector<unsigned int> level = PlasticScorer_Veto->GetLevel(i); 
    Energy_tmp = PlasticScorer_Veto->GetEnergy(i);
    Light_tmp = PlasticScorer_Veto->GetLight(i);
    Energy = RandGauss::shoot(Energy_tmp, Energy_tmp*Nebula_NS::ResoEnergy);
    Light = RandGauss::shoot(Light_tmp, Light_tmp*Nebula_NS::ResoLight);
 
    if(Light>Nebula_NS::LightThreshold){
      double Time = RandGauss::shoot(PlasticScorer_Veto->GetTime(i),Nebula_NS::ResoTime);
      //cout << "Time is " << Time << endl;
      double Position = RandGauss::shoot(PlasticScorer_Veto->GetPositionY(i),Nebula_NS::ResoPosition);
      //cout << "Position is " << Position << endl;
      int DetectorNbr = level[0] + m_TotalModule;
      // cout << "Veto ID: " << DetectorNbr << endl;
      
      m_Event->SetChargeUp(DetectorNbr,Light*exp(-(Nebula_NS::VetoHeight/2-Position)/Nebula_NS::Attenuation));
      m_Event->SetChargeDown(DetectorNbr,Light*exp(-(Nebula_NS::VetoHeight/2+Position)/Nebula_NS::Attenuation));
      
      Time_up = (Nebula_NS::VetoHeight/2-Position)/(c_light/Nebula_NS::MaterialIndex) + Time;
      //cout << "Time_up is " << Time_up << endl;
      m_Event->SetTimeUp(DetectorNbr,Time_up);
      
      Time_down = (Nebula_NS::VetoHeight/2+Position)/(c_light/Nebula_NS::MaterialIndex) + Time;
      //cout << "Time_down is " << Time_down << endl;
      m_Event->SetTimeDown(DetectorNbr,Time_down);
      
    }
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Nebula::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_ModuleScorer = CheckScorer("NebulaModuleScorer",already_exist) ;
  m_VetoScorer = CheckScorer("NebulaVetoScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialise
  // Module 
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* ModulePlasticBar= new PlasticBar::PS_PlasticBar("ModulePlasticBar",level, 0);
  G4VPrimitiveScorer* ModuleInteraction= new InteractionScorers::PS_Interactions("ModuleInteraction",ms_InterCoord, 0);
  G4VPrimitiveScorer* ModuleProcess= new ProcessScorers::PS_Process("ModuleProcess", 0);
  //and register it to the multifunctionnal detector
  m_ModuleScorer->RegisterPrimitive(ModulePlasticBar);
  m_ModuleScorer->RegisterPrimitive(ModuleInteraction);
  m_ModuleScorer->RegisterPrimitive(ModuleProcess);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ModuleScorer) ;

  // Veto 
  G4VPrimitiveScorer* VetoPlasticBar= new PlasticBar::PS_PlasticBar("VetoPlasticBar",level, 0);
  G4VPrimitiveScorer* VetoInteraction= new InteractionScorers::PS_Interactions("VetoInteraction",ms_InterCoord, 0);
  G4VPrimitiveScorer* VetoProcess= new ProcessScorers::PS_Process("ModuleProcess", 0);
  //and register it to the multifunctionnal detector
  m_VetoScorer->RegisterPrimitive(VetoPlasticBar);
  m_VetoScorer->RegisterPrimitive(VetoInteraction);
  m_VetoScorer->RegisterPrimitive(VetoProcess);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_VetoScorer) ;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Nebula::Construct(){
  return  (NPS::VDetector*) new Nebula();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Nebula{
    public:
      proxy_nps_Nebula(){
        NPS::DetectorFactory::getInstance()->AddToken("NEBULA","NEBULA");
        NPS::DetectorFactory::getInstance()->AddDetector("NEBULA",Nebula::Construct);
      }
  };

  proxy_nps_Nebula p_nps_Nebula;
}
