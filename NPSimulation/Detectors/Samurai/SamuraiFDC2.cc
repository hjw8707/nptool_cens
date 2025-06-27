/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey ANNE  contact address: anne@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    : august 2024                                              *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  SamuraiFDC2 simulation                              *
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

#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

// NPTool header
#include "SamuraiFDC2.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

#include "WireScorers.hh"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SamuraiFDC2_NS{
  // Samurai magnet construction paramethers

  //Main outer box
  const double FDC2_Width = 2296*mm; //(x)
  const double FDC2_Height = 836*mm;//(y)
  const double FDC2_Depth = 650*mm;//(z)
  const string FDC2_Material_Void = "G4_Galactic";

  const string FDC2_Gas = "CH4_60_He_40";
  const double FDC2_Temperature = 298.15; // K
  const double FDC2_Pressure = 1.0; //atm

  //Detector Number
  const short int FDC2_DetectorNumber = 2;

  //Wires
  const double Wire_length = 3000*mm;
  const double Wire_diameter = 19.99*mm;
  const string Wire_Gas = "CH4_60_He_40";
  const double Wire_Temperature = 298.15; // K 
  const double Wire_Pressure = 1.0; //atm

  const int Number_Of_Layer = 14; // from 0 to 13 in the xml file
  const int Number_Of_Wire_By_Layer = 112; // from 0 to 111 in the xml file
  const double Drift_Speed = 1.33e-4; //------> FIX ME!!! not important (yet)

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Samurai Specific Method
SamuraiFDC2::SamuraiFDC2(){

  //Visualization attributes
  m_VisFDC2 = new G4VisAttributes(G4Colour(1,0,1,0.3));
  //Logical volumes
  m_FDC2 = NULL;

  //Wire
  m_VisWire = new G4VisAttributes(G4Colour(1.0,0,0.0,0.7));
  //Logical volumes
  m_Wire = NULL;
  //Scorer
  m_WireScorerFDC2 = NULL;

  //Data event
  m_Event = new TSamuraiFDC2Data;
  
}

SamuraiFDC2::~SamuraiFDC2(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC2::AddDetector(G4ThreeVector Mag_Pos, double Mag_Angle, G4ThreeVector Offset, double Off_Angle){

  m_Angle = Mag_Angle + (90.*deg - Off_Angle);

  Offset.rotateY(-(m_Angle));
  m_Pos = Mag_Pos + Offset;
  
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* SamuraiFDC2::BuildFDC2(){
  if(!m_FDC2){
    //Shape - G4Box
    G4Box* box = new G4Box("FDC2_Box",SamuraiFDC2_NS::FDC2_Width*0.5,
			   SamuraiFDC2_NS::FDC2_Height*0.5,SamuraiFDC2_NS::FDC2_Depth*0.5);
    
    //Material
    G4Material* Gas = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiFDC2_NS::FDC2_Gas,SamuraiFDC2_NS::FDC2_Pressure, SamuraiFDC2_NS::FDC2_Temperature);

    //Logical Volume
    m_FDC2 = new G4LogicalVolume(box, Gas, "logic_SamuraiFDC2_box",0,0,0);
    m_FDC2->SetVisAttributes(m_VisFDC2);
  }
  return m_FDC2;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiFDC2::BuildWire(){
if(!m_Wire){

  //Shape - G4Tubs
  G4Tubs* solidWire = new G4Tubs("FDC2_Wire", 0.0, SamuraiFDC2_NS::Wire_diameter*0.5,SamuraiFDC2_NS::Wire_length*0.5, 0., 360.);

  
  //material
  G4Material* GasWire = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiFDC2_NS::Wire_Gas, SamuraiFDC2_NS::Wire_Pressure,SamuraiFDC2_NS::Wire_Temperature);
  
  
  //Logical Volume
  m_Wire = new G4LogicalVolume(solidWire, GasWire, "logic_SamuraiWire_tub",0,0,0);
  m_Wire->SetVisAttributes(m_VisWire);
  m_Wire->SetSensitiveDetector(m_WireScorerFDC2);
  
  }
  return m_Wire;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// (Called in DetecorConstruction::ReadDetectorConfiguration Method)
void SamuraiFDC2::ReadConfiguration(NPL::InputParser parser){

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Samurai");
  vector<NPL::InputBlock*> blocks2 = parser.GetAllBlocksWithToken("SAMURAIFDC2");

  G4ThreeVector Mag_Pos;
  double Mag_Angle;
  G4ThreeVector Offset;
  double Off_Angle;
  
  if(blocks.size()==1 && blocks2.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "/////// Samurai FDC2 found with Samurai Magnet ///////" << endl;
    }
    vector<string> cart = {"POS","ANGLE"};
    vector<string> sphe = {"R","Theta","Phi","ANGLE"};
    vector<string> xml = {"XML","Offset", "InvertX","InvertY","InvertD"};
    
    if(blocks[0]->HasTokenList(cart)){
      Mag_Pos = NPS::ConvertVector(blocks[0]->GetTVector3("POS", "cm"));
      Mag_Angle = blocks[0]->GetDouble("ANGLE","deg");
    }
    
    if(blocks[0]->HasTokenList(sphe)){
      double R = blocks[0]->GetDouble("R","mm");
      double Theta = blocks[0]->GetDouble("Theta","deg");
      double Phi = blocks[0]->GetDouble("Phi","deg");
      Mag_Pos.setMag(R);
      Mag_Pos.setTheta(Theta);
      Mag_Pos.setPhi(Phi);
      Mag_Angle = blocks[0]->GetDouble("ANGLE","deg");
    }

    if(blocks2[0]->HasTokenList(xml)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FDC2 XML file " <<  endl;
      string xml_file = blocks2[0]->GetString("XML");
      Offset = NPS::ConvertVector(blocks2[0]->GetTVector3("Offset", "mm"));
      Off_Angle = blocks2[0]->GetDouble("OffAngle","deg");
      bool invert_x = blocks2[0]->GetInt("InvertX");
      bool invert_y = blocks2[0]->GetInt("InvertY");
      bool invert_z = blocks2[0]->GetInt("InvertD");
      ReadXML(xml_file,Offset,invert_x,invert_y);
    }    

    AddDetector(Mag_Pos, Mag_Angle, Offset, Off_Angle);

  }// SamuraiFDC2 and SamuraiMagnet

  else if(blocks.size()==0 && blocks2.size()==1)
    {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
	{
	  cout << "/////// Samurai FDC2 ///////" << endl;
	}

      vector<string> xml = {"XML","Offset", "InvertX","InvertY","InvertD"};
      if(blocks2[0]->HasTokenList(xml)){//--------------> nouveau
	if(NPOptionManager::getInstance()->GetVerboseLevel())
	  cout << endl << "////  FDC2 XML file " <<  endl;//--------------> nouveau
	string xml_file = blocks2[0]->GetString("XML");//--------------> nouveau
	Offset = NPS::ConvertVector(blocks2[0]->GetTVector3("Offset", "mm"));//-----------> nouveau
	Off_Angle = blocks2[0]->GetDouble("OffAngle","deg");//-----------> nouveau
	bool invert_x = blocks2[0]->GetInt("InvertX");//--------------> nouveau
	bool invert_y = blocks2[0]->GetInt("InvertY");//--------------> nouveau
	bool invert_z = blocks2[0]->GetInt("InvertD");//--------------> nouveau
	ReadXML(xml_file,Offset,invert_x,invert_y);//--------------> nouveau
      }    

       m_Pos = Mag_Pos + Offset;

    }//SamuraiFDC2

  
  else{
    cout << "ERROR: there should be only one Samurai magnet, check your input file" << endl;
    exit(1);
  }
    
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC2::ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY){ 
  NPL::XmlParser xml;
  xml.LoadFile(xml_file);

  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("SAMURAIFDC2");
  int NumberOfWires=0;
  for(unsigned int i = 0 ; i < b.size() ; i++){
    NumberOfWires++;

    // Wire ID in general
    auto WireId  = b[i]->AsInt("ID");

    // Layer Number
    auto LayerNbr = b[i]->AsInt("layer");
    
    // Wire direction
    auto DirectionWire = b[i]->AsString("anodedir");
    //Wire id in the layer
    auto WireNbr  = b[i]->AsInt("wireid");

    //Wire x or y position
    auto PositionX = b[i]->AsDouble("wirepos");
    auto PositionZ = b[i]->AsDouble("wirez");

    // position wire and direction wire
    // X , U or V
    m_PositionWire[WireId] = G4ThreeVector(PositionX,0.0,PositionZ);
    m_IsHorizontal[WireId] = false;
    
    m_WireAngle[WireId] = 0.0*deg;
    
    if(DirectionWire.find("X")!=std::string::npos){ //X
      m_WireAngle[WireId] = 0.0*deg;
    }
    
    else if(DirectionWire.find("U")!=std::string::npos){ //U
      m_WireAngle[WireId] = -30.0*deg;
    }
    
    else if (DirectionWire.find("V")!=std::string::npos){ //V
      m_WireAngle[WireId] = 30.0*deg;
    }
    
    m_LayerNbr[WireId] = LayerNbr;
    m_WireNbr[WireId] = WireNbr;
    
  } //end of for loop
  cout << " -> " << NumberOfWires << " wires found" << endl;

}//end of ReadXML


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// (Called After DetectorConstruction::AddDetector Method)
void SamuraiFDC2::ConstructDetector(G4LogicalVolume* world){

  // Box
  G4RotationMatrix* Rot = new G4RotationMatrix();
  Rot->rotateY(m_Angle);//rotation in volume's reference frame
  new G4PVPlacement(Rot, m_Pos,
          BuildFDC2(), "SamuraiFDC2", world, false, 0);

  // Wires
  G4RotationMatrix* RotWire = new G4RotationMatrix();// wire's rotation at the end
  RotWire->set(0,0,0);

  G4RotationMatrix* BoxRot = new G4RotationMatrix();
  BoxRot->set(0,0,0);

  
  //material
  G4Material* GasWire = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiFDC2_NS::Wire_Gas, SamuraiFDC2_NS::Wire_Pressure,SamuraiFDC2_NS::Wire_Temperature);

  
  for(auto pos : m_PositionWire){ //loop on wires
    
      RotWire->rotateX(-90.*deg);    
      RotWire->rotateZ(m_WireAngle[pos.first]);

      G4ThreeVector position = pos.second;
      double X = position[0];
      double Y = position[1];
      double Z = position[2];

      double newX, newY;
      newX = X*cos(-m_WireAngle[pos.first])+Y*sin(-m_WireAngle[pos.first]);
      newY = -X*sin(-m_WireAngle[pos.first]) + Y*cos(-m_WireAngle[pos.first]);

      m_PositionWire[pos.first] = G4ThreeVector(newX, newY, Z);
      G4ThreeVector positionWire(newX, newY, Z);

      
      //////////////////////////////////////////////////////
      G4Box* box = new G4Box("FDC2_Box",SamuraiFDC2_NS::FDC2_Width*0.5,
			     SamuraiFDC2_NS::FDC2_Depth*0.5,SamuraiFDC2_NS::FDC2_Height*0.5);
      G4Tubs* solidWire = new G4Tubs("FDC2_Wire", 0.0, SamuraiFDC2_NS::Wire_diameter*0.5,SamuraiFDC2_NS::Wire_length*0.5, 0., 360.);

      BoxRot->rotateY(m_WireAngle[pos.first]);
      
      G4ThreeVector BoxPosition(-X,0, 0);
      G4Transform3D TransformBox(*BoxRot, BoxPosition);
      
      G4VSolid* Intersection = new G4IntersectionSolid("FDC2_Wire", solidWire, box, TransformBox);
      m_Wire = new G4LogicalVolume(Intersection, GasWire, "logic_SamuraiWire_tub",0,0,0);
      m_Wire->SetVisAttributes(m_VisWire);
      m_Wire->SetSensitiveDetector(m_WireScorerFDC2);

      // intersection
      new G4PVPlacement(G4Transform3D(*RotWire, positionWire),
			m_Wire,
			"WireModule",m_FDC2,false,pos.first);

      RotWire->set(0,0,0);
      BoxRot->set(0,0,0);
   
  }// end for loop on wires
  
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SamuraiFDC2::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("SamuraiFDC2")){
    pTree->Branch("SamuraiFDC2", "TSamuraiFDC2Data", &m_Event) ;
  }
  pTree->SetBranchAddress("SamuraiFDC2", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Read sensitive part and fill the Root tree.
// (Called at in the EventAction::EndOfEventAvtion)

void SamuraiFDC2::ReadSensitive(const G4Event* event){
  
  m_Event->Clear();

  WireScorers::PS_Wire* WireScorerFDC2= (WireScorers::PS_Wire*) m_WireScorerFDC2->GetPrimitive(0);

  unsigned int sizeWire =WireScorerFDC2->GetMult();
  for(unsigned int i = 0 ; i < sizeWire ; i++)
    {
      int layer = WireScorerFDC2->GetLayerNumber(i);
      int wireNbr = WireScorerFDC2->GetWireNumber(i);
      double time = WireScorerFDC2->GetTime(i);
      int edge = WireScorerFDC2->GetEdge(i);
      double DriftLength = WireScorerFDC2->GetDriftLength(i);    
      
      m_Event->SetData(SamuraiFDC2_NS::FDC2_DetectorNumber, layer,
		       wireNbr, DriftLength, edge);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  

void SamuraiFDC2::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_WireScorerFDC2 = CheckScorer("WireScorerFDC2",already_exist);

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  //Wire
  G4VPrimitiveScorer* InteractionWire= new WireScorers::PS_Wire("WS_FDC2",level, SamuraiFDC2_NS::Number_Of_Layer,SamuraiFDC2_NS::Number_Of_Wire_By_Layer,SamuraiFDC2_NS::Drift_Speed, m_WireAngle) ;
  m_WireScorerFDC2->RegisterPrimitive(InteractionWire);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_WireScorerFDC2) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////

NPS::VDetector* SamuraiFDC2::Construct(){
  return  (NPS::VDetector*) new SamuraiFDC2();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_samuraiFDC2{
    public:
      proxy_nps_samuraiFDC2(){
        NPS::DetectorFactory::getInstance()->AddToken("SAMURAIFDC2","SAMURAIFDC2");
        NPS::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC2",SamuraiFDC2::Construct);
      }
  };
  
  proxy_nps_samuraiFDC2 p_nps_samuraiFDC2;
}


