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
 *  This class describe  NeuLAND simulation                                   *
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
#include "NeuLAND.hh"
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
namespace NeuLAND_NS{
  // Energy and time Resolution
  const double LightThreshold  = 0.1*MeV;
  const double ResoTime        = 0.75/2.355*ns; //0.75
  const double ResoEnergy      = 0.1/2.355*MeV; 
  const double ResoLight       = 0.1/2.355*MeV; 
  const double ResoPosition    = 1.0*um; //1.0
  const double ModuleWidth     = 50*mm ;
  const double ModuleLength    = 50*mm ;
  const double ModuleHeight    = 2500*mm ;
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
// NeuLAND Specific Method
NeuLAND::NeuLAND(){
  m_Event = new TNeuLANDData() ;
  m_ModuleScorer = 0;
  m_VetoScorer = 0;
  m_Module = 0;
  m_Veto = 0;


  // RGB Color + Transparency
  m_VisModule = new G4VisAttributes(G4Colour(0.263, 0.682, 0.639, 0.3));   
  //m_VisModule = new G4VisAttributes(G4Colour(0.145, 0.384, 0.596, 1));   
  m_VisVeto   = new G4VisAttributes(G4Colour(0.4, 0.4, 0.4, 0.8));   
  m_VisPMT    = new G4VisAttributes(G4Colour(0.1, 0.1, 0.1, 1));   
  m_VisFrame  = new G4VisAttributes(G4Colour(0, 0.3, 1, 0.5));   

}

NeuLAND::~NeuLAND(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void NeuLAND::AddWall(G4ThreeVector Pos, int NbrModule, bool Veto, bool Frame){
  // Convert the Pos value to R theta Phi as Spherical coordinate is easier in G4 
  m_Pos.push_back(Pos);
  m_NbrModule.push_back(NbrModule);
  m_HasVeto.push_back(Veto);
  m_HasFrame.push_back(Frame);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* NeuLAND::BuildModule(){
  if(!m_Module){
    G4Box* box = new G4Box("NeuLAND_Module",NeuLAND_NS::ModuleWidth*0.5,
        NeuLAND_NS::ModuleHeight*0.5,NeuLAND_NS::ModuleLength*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(NeuLAND_NS::Material);
    m_Module = new G4LogicalVolume(box,DetectorMaterial,"logic_NeuLAND_Module",0,0,0);
    m_Module->SetVisAttributes(m_VisModule);
    m_Module->SetSensitiveDetector(m_ModuleScorer);
  }
  return m_Module;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* NeuLAND::BuildVeto(){
  if(!m_Veto){
    G4Box* box = new G4Box("NeuLAND_Veto",NeuLAND_NS::VetoWidth*0.5,
        NeuLAND_NS::VetoHeight*0.5,NeuLAND_NS::VetoLength*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(NeuLAND_NS::Material);
    
    m_Veto = new G4LogicalVolume(box,DetectorMaterial,"logic_NeuLAND_Veto",0,0,0);
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
// Called in DetectorConstruction::ReadDetectorConfiguration Method
void NeuLAND::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEULAND");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  // define an entire wall
  vector<string> wall = {"Pos","NumberOfModule","Veto","Frame"};

  // use an experiment xml file to position bars individually
  vector<string> xml= {"XML","Offset","InvertX","InvertY"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(wall)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  NeuLAND " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("Pos","mm"));
      int NbrModule = blocks[i]->GetInt("NumberOfModule");
      bool Veto = blocks[i]->GetInt("Veto");
      bool Frame= blocks[i]->GetInt("Frame");
      AddWall(Pos,NbrModule,Veto,Frame);
    }
    else if(blocks[i]->HasTokenList(xml)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  NeuLAND XML file " << i+1 <<  endl;
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
void NeuLAND::ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY){ 
  NPL::XmlParser xml;
  xml.LoadFile(xml_file);
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("NEULAND");  
  int NumberOfBars=0;
  for(unsigned int i = 0 ; i < b.size() ; i++){
    NumberOfBars++;
    unsigned int id = b[i]->AsInt("ID");
    
    // position
    auto PositionX = b[i]->AsDouble("xpos");
    //cout << "//// " << PositionX << " position en x " << endl;
    auto PositionY = b[i]->AsDouble("ypos");
    //cout << "//// " << PositionY << " position en y " << endl;
    auto PositionZ = b[i]->AsDouble("zpos");
    //cout << "//// " << PositionZ << " position en z " << endl;
    //cout <<"============================" << endl;
    
    // SubLayer 0 is used for Veto
    auto SubLayer  = b[i]->AsInt("SubLayer");
    // Name "NoUseX" is used to silence bars
    auto nousestr  = b[i]->AsString("NAME");

    auto DirectionBar = b[i]->AsString("direction");
    //cout << "direction des barres -> " << DirectionBar << endl;
    // cout << "id = " << id << endl;
    
    // Remove unused bar
    if(nousestr.find("NoUse")==std::string::npos && PositionX!=PositionZ){
      if(InvertX)
        PositionX*=-1;
      if(InvertY)
        PositionY*=-1;
      m_PositionBar[id]= G4ThreeVector(PositionX,PositionY,PositionZ)+offset;

      // Direction bar
      if(DirectionBar.find("V")==std::string::npos)
	{m_IsHorizontal[id] = true;
	  // cout << DirectionBar << endl; 
	}
      else m_IsHorizontal[id] = false;

      // cout << " m_IsHorizontal[id] = " << m_IsHorizontal[id] << endl;
      
      if(SubLayer)
        m_IsVetoBar[id]= false;
      else
        m_IsVetoBar[id]= true;
    }
    
  } //end of for 
  cout << " -> " << NumberOfBars << " bars found" << endl;

  
} //end of ReadXML

// Construct detector and inialise sensitive part.
// Called After DetectorConstruction::AddDetector Method
void NeuLAND::ConstructDetector(G4LogicalVolume* world){
  
  // Start with XML case
  G4RotationMatrix* Rot = new G4RotationMatrix();
    
  for(auto pos : m_PositionBar){

    //cout << "m_IsHorizontal " <<  m_IsHorizontal[pos.first] << endl;
    
    if(m_IsHorizontal[pos.first])
      {
         Rot->rotateZ(90.*deg);
      }

    
    
    if(!m_IsVetoBar[pos.first]){
      new G4PVPlacement(G4Transform3D(*Rot,pos.second),
            BuildModule(),
            "NeuLANDModule",world,false,pos.first);
    }
    else{
      new G4PVPlacement(G4Transform3D(*Rot,pos.second),
            BuildVeto(),
	    "NeuLANDModule",world,false,pos.first); //"NeuLANDVeto"
      }

    Rot->set(0,0,0);
    }


  //Not XML
  unsigned int nbrM = 1 ;
  unsigned int nbrV = 1 ;
  
  for (unsigned short i = 0 ; i < m_Pos.size() ; i++) {
    for (unsigned short m = 0 ; m < m_NbrModule[i] ; m++) {
      double offset = (NeuLAND_NS::ModuleWidth+NeuLAND_NS::InterModule)*(-m_NbrModule[i]*0.5+m)+NeuLAND_NS::ModuleWidth*0.5;
      G4ThreeVector Offset(offset,0,0);
      new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
          BuildModule(),
          "NeuLANDModule",world,false,nbrM++);
    }

    if(m_HasVeto[i]){
      if(m_NbrModule[i] > 15){
        for (unsigned short m = 0 ; m < NeuLAND_NS::VetoPerWall ; m++) {
          double offset = (NeuLAND_NS::VetoWidth+NeuLAND_NS::InterVeto)*(-NeuLAND_NS::VetoPerWall*0.5+m)+NeuLAND_NS::VetoWidth*0.5;
          G4ThreeVector Offset(offset,0,-NeuLAND_NS::WallToVeto);
          new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
            BuildVeto(),
            "NeuLANDVeto",world,false,nbrV++);
        }
      }
      else{
        for (unsigned short m = 0 ; m < NeuLAND_NS::VetoPerExpand ; m++) {
          double offset = (NeuLAND_NS::VetoWidth+NeuLAND_NS::InterVeto)*(-NeuLAND_NS::VetoPerExpand*0.5+m)+NeuLAND_NS::VetoWidth*0.5;
          G4ThreeVector Offset(offset,0,-NeuLAND_NS::WallToVeto);
          new G4PVPlacement(G4Transform3D(*Rot,m_Pos[i]+Offset),
            BuildVeto(),
            "NeuLANDVeto",world,false,nbrV++);
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetectorConstruction::AddDetector Method
void NeuLAND::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("NeuLAND")){
    pTree->Branch("NeuLAND", "TNeuLANDData",&m_Event) ;
  }
  pTree->SetBranchAddress("NeuLAND", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void NeuLAND::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // PlasticBar scorer
  PlasticBar::PS_PlasticBar* PlasticScorer_Module = (PlasticBar::PS_PlasticBar*) m_ModuleScorer->GetPrimitive(0);
  PlasticBar::PS_PlasticBar* PlasticScorer_Veto = (PlasticBar::PS_PlasticBar*) m_VetoScorer->GetPrimitive(0);
  // Should we put a ProcessScorer here to get the info if the particle is first neutron and give it to NeuLANDData ?
  
  double Time_up, Time_down;
  double Energy_tmp, Light_tmp;

  //////////// TRIAL TO GET THE OPTICAL INDEX FROM MATERIAL PROPERTIES /////////////
  //Trying to get Optical Index from Material directly
  //const G4Material* aMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(NeuLAND_NS::Material);
  //G4MaterialPropertiesTable* aMaterialPropertiesTable = aMaterial->GetMaterialPropertiesTable();
  //if(!aMaterialPropertiesTable->PropertyExists("RINDEX")){
  //  MaterialIndex = !aMaterialPropertiesTable->GetConstProperty("RINDEX"); 
  //}
  //else{
  //  MaterialIndex = 0; 
  //}
  //cout << MaterialManager::getInstance()->GetMaterialFromLibrary(NeuLAND_NS::Material)->GetMaterialPropertiesTable()->GetMaterialPropertyNames()[0] << endl;
  //////////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////// MODULE SCORER //////////////////////////////////
  unsigned int ModuleHits_size = PlasticScorer_Module->GetMult(); 
  for(unsigned int i = 0 ; i < ModuleHits_size ; i++){
    vector<unsigned int> level = PlasticScorer_Module->GetLevel(i); 
    Energy_tmp = PlasticScorer_Module->GetEnergy(i);
    Light_tmp = PlasticScorer_Module->GetLight(i);
    Energy = RandGauss::shoot(Energy_tmp, Energy_tmp*NeuLAND_NS::ResoEnergy);
    Light = RandGauss::shoot(Light_tmp, Light_tmp*NeuLAND_NS::ResoLight);

    if(Light>NeuLAND_NS::LightThreshold){
      int DetectorNbr = level[0];
      double PositionY = RandGauss::shoot(PlasticScorer_Module->GetPositionY(i),NeuLAND_NS::ResoPosition);
      double PositionX = RandGauss::shoot(PlasticScorer_Module->GetPositionX(i),NeuLAND_NS::ResoPosition);

      
      double Position;

      if((1<=DetectorNbr and DetectorNbr <=50)  or (101<=DetectorNbr and DetectorNbr<=150) or (201<=DetectorNbr and DetectorNbr<=250) or (301<=DetectorNbr and DetectorNbr<=350) )
	{
	  Position = PositionX;
	}

      else
	{
	  Position = PositionY;
	}

      //cout << "--------------" << endl;
/*      m_Event->SetChargeUp(DetectorNbr,Light*exp(-(NeuLAND_NS::ModuleHeight/2-Position)/NeuLAND_NS::Attenuation));
      m_Event->SetChargeDown(DetectorNbr,Light*exp(-(NeuLAND_NS::ModuleHeight/2+Position)/NeuLAND_NS::Attenuation));
     
      // Take Time and Position and compute Tup and Tdown
      double Time = RandGauss::shoot(PlasticScorer_Module->GetTime(i),NeuLAND_NS::ResoTime);

      Time_up = (NeuLAND_NS::ModuleHeight/2-Position)/(c_light/NeuLAND_NS::MaterialIndex) + Time;
      m_Event->SetTimeUp(DetectorNbr,Time_up);
      
      Time_down = (NeuLAND_NS::ModuleHeight/2+Position)/(c_light/NeuLAND_NS::MaterialIndex) + Time;
      m_Event->SetTimeDown(DetectorNbr,Time_down);
      */
    }
  }

  ///////////////////////////////// VETO SCORER //////////////////////////////////
  unsigned int VetoHits_size = PlasticScorer_Veto->GetMult(); 
  for(unsigned int i = 0 ; i < VetoHits_size ; i++){
    vector<unsigned int> level = PlasticScorer_Veto->GetLevel(i); 
    Energy_tmp = PlasticScorer_Veto->GetEnergy(i);
    Light_tmp = PlasticScorer_Veto->GetLight(i);
    Energy = RandGauss::shoot(Energy_tmp, Energy_tmp*NeuLAND_NS::ResoEnergy);
    Light = RandGauss::shoot(Light_tmp, Light_tmp*NeuLAND_NS::ResoLight);

    if(Light>NeuLAND_NS::LightThreshold){
      double Time = RandGauss::shoot(PlasticScorer_Veto->GetTime(i),NeuLAND_NS::ResoTime);
      //cout << "Time is " << Time << endl;
      double Position = RandGauss::shoot(PlasticScorer_Veto->GetPositionY(i),NeuLAND_NS::ResoPosition);
      //cout << "Position is " << Position << endl;
      int DetectorNbr = level[0] + m_TotalModule;
      //cout << "Veto ID: " << DetectorNbr << endl;
      /*
      m_Event->SetChargeUp(DetectorNbr,Light*exp(-(NeuLAND_NS::VetoHeight/2-Position)/NeuLAND_NS::Attenuation));
      m_Event->SetChargeDown(DetectorNbr,Light*exp(-(NeuLAND_NS::VetoHeight/2+Position)/NeuLAND_NS::Attenuation));
      
      Time_up = (NeuLAND_NS::VetoHeight/2-Position)/(c_light/NeuLAND_NS::MaterialIndex) + Time;
      //cout << "Time_up is " << Time_up << endl;
      m_Event->SetTimeUp(DetectorNbr,Time_up);
      
      Time_down = (NeuLAND_NS::VetoHeight/2+Position)/(c_light/NeuLAND_NS::MaterialIndex) + Time;
      //cout << "Time_down is " << Time_down << endl;
      m_Event->SetTimeDown(DetectorNbr,Time_down);
    */
    }
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void NeuLAND::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_ModuleScorer = CheckScorer("NeuLANDModuleScorer",already_exist) ;
  m_VetoScorer = CheckScorer("NeuLANDVetoScorer",already_exist) ;

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
NPS::VDetector* NeuLAND::Construct(){
  return  (NPS::VDetector*) new NeuLAND();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_NeuLAND{
    public:
      proxy_nps_NeuLAND(){
        NPS::DetectorFactory::getInstance()->AddToken("NEULAND","NEULAND");
        NPS::DetectorFactory::getInstance()->AddDetector("NEULAND",NeuLAND::Construct);
      }
  };

  proxy_nps_NeuLAND p_nps_NeuLAND;
}
