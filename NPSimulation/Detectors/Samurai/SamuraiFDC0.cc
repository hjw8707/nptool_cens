/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Anne  contact address: anne@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : january 2024                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  SamuraiFDC0 simulation                              *
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
#include "G4Tubs.hh"

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
#include "SamuraiFDC0.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

#include "WireScorers.hh"

#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"

using namespace std;
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace SamuraiFDC0_NS{

  //Main outer box
  // const double FDC0_Width = 160*mm; //(x)
  // const double FDC0_Height = 160*mm;//(y)
  // const double FDC0_Depth = 160*mm;//(z)

  const double FDC0_Width = 170*mm; //(x)
  const double FDC0_Height = 170*mm;//(y)
  const double FDC0_Depth = 170*mm;//(z)

  const string FDC0_Gas = "CH4_60_He_40";
  const double FDC0_Temperature = 298.15; // K
  const double FDC0_Pressure = 1.0; //atm

  //Detector Number
  const short int FDC0_DetectorNumber = 0;

  //Wires
  const double Wire_Length = 150*mm;
  const double Wire_Diameter = 5*mm;
  const string Wire_Gas = "CH4_60_He_40";
  const double Wire_Temperature = 298.15; // K 
  const double Wire_Pressure = 1.0; //atm 

  const int Number_Of_Layer = 8; // from 0 to 7 in the xml file
  const int Number_Of_Wire_By_Layer = 32; // from 0 to 31 in the xml file
  const double Drift_Speed = 1.33e-4; //------> FIX ME!!! not important (yet)

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Samurai Specific Method
SamuraiFDC0::SamuraiFDC0(){

  //Visualization attributes

  //Main box
  m_VisFDC0 = new G4VisAttributes(G4Colour(0.0,1.0,0,0.1));
  //Logical volumes
  m_FDC0 = NULL;

  //Wire
  m_VisWire = new G4VisAttributes(G4Colour(1.0,0,0.0,1.0));
  //Logical volumes
  m_Wire = NULL;
  //Scorer
  m_WireScorerFDC0 = NULL;

  //Data event
  m_Event = new TSamuraiFDC0Data; 
  
  //Hexagon
  m_Hexagon = NULL;
}

SamuraiFDC0::~SamuraiFDC0(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiFDC0::AddDetector(G4ThreeVector Mag_Pos, G4ThreeVector Offset){

  m_Pos = Mag_Pos + Offset;

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiFDC0::BuildFDC0(){
  if(!m_FDC0){
    //Shape - G4Box
    G4Box* box = new G4Box("FDC0_Box",SamuraiFDC0_NS::FDC0_Width*0.5,
			   SamuraiFDC0_NS::FDC0_Height*0.5,SamuraiFDC0_NS::FDC0_Depth*0.5);

    //Material - Gas
    G4Material* Gas = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiFDC0_NS::FDC0_Gas, SamuraiFDC0_NS::FDC0_Pressure,SamuraiFDC0_NS::FDC0_Temperature);

    //Logical Volume
    m_FDC0 = new G4LogicalVolume(box, Gas , "logic_SamuraiFDC0_box",0,0,0);
    m_FDC0->SetVisAttributes(m_VisFDC0);
  }
  return m_FDC0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiFDC0::BuildWire(){
if(!m_Wire){
  //Shape - G4Tubs
  G4Tubs* solidWire = new G4Tubs("FDC0_Wire", 0.0, SamuraiFDC0_NS::Wire_Diameter*0.5,SamuraiFDC0_NS::Wire_Length*0.5, 0., 360.);

  //Material - Gas
  G4Material* GasWire = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiFDC0_NS::Wire_Gas, SamuraiFDC0_NS::Wire_Pressure,SamuraiFDC0_NS::Wire_Temperature);
  //Logical Volume
  m_Wire = new G4LogicalVolume(solidWire, GasWire, "logic_SamuraiWire_tub",0,0,0);
  m_Wire->SetVisAttributes(m_VisWire);
  m_Wire->SetSensitiveDetector(m_WireScorerFDC0);
  
  }
  return m_Wire;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiFDC0::BuildHexagon(){
if(!m_Hexagon){
  //Shape 
 
  const G4int nsect = 6;
  std::vector<G4TwoVector> polygon(nsect);
  G4double ang = twopi/nsect;

  G4double dz = SamuraiFDC0_NS::Wire_Length*0.5;
  G4double rmax = SamuraiFDC0_NS::Wire_Diameter*0.5;
  
  for (G4int i = 0; i < nsect; ++i)
    {
      G4double phi = i*ang;
      G4double cosphi = std::cos(phi);
      G4double sinphi = std::sin(phi);
      polygon[i].set(rmax*cosphi, rmax*sinphi);
    }
  
  G4TwoVector offsetA(0,0), offsetB(0,0);
  G4double scaleA = 1, scaleB = 1;
  G4VSolid* solidHexagon = new G4ExtrudedSolid("Extruded", polygon, dz, offsetA, scaleA, offsetB, scaleB);
  

  
  //Material - Gas
  G4Material* GasWire = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiFDC0_NS::Wire_Gas, SamuraiFDC0_NS::Wire_Pressure,SamuraiFDC0_NS::Wire_Temperature);
  
  //Logical Volume
  m_Hexagon = new G4LogicalVolume(solidHexagon, GasWire, "logic_SamuraiWire_tub",0,0,0);
  m_Hexagon->SetVisAttributes(m_VisWire);
  m_Hexagon->SetSensitiveDetector(m_WireScorerFDC0);
  
  }
  return m_Hexagon;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// (Called in DetecorConstruction::ReadDetectorConfiguration Method)
void SamuraiFDC0::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks2 = parser.GetAllBlocksWithToken("SAMURAIFDC0");

  if(blocks2.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "/////// Samurai FDC0 found ///////" << endl;
    }

    vector<string> xml = {"XML","Offset", "InvertX","InvertY","InvertD"};

    G4ThreeVector Mag_Pos;
    G4ThreeVector Offset;
     
    if(blocks2[0]->HasTokenList(xml)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FDC0 XML file " <<  endl;
      string xml_file = blocks2[0]->GetString("XML");
      Offset = NPS::ConvertVector(blocks2[0]->GetTVector3("Offset", "mm"));
      bool invert_x = blocks2[0]->GetInt("InvertX");
      bool invert_y = blocks2[0]->GetInt("InvertY");
      bool invert_z = blocks2[0]->GetInt("InvertD");
      ReadXML(xml_file,Offset,invert_x,invert_y);
    }
    
    AddDetector(Mag_Pos, Offset);

  }
  else{
    cout << "ERROR for FDC0: check your input file " << endl;
    exit(1);
  }
     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//////////////////// nouveau ///////////////////
void SamuraiFDC0::ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY){ 
  NPL::XmlParser xml;
  xml.LoadFile(xml_file);

  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("SAMURAIFDC0");
  int NumberOfWires=0;
  for(unsigned int i = 0 ; i < b.size() ; i++){
    NumberOfWires++;

    // Wire ID in general
    auto ID  = b[i]->AsInt("ID");

    // Layer Number
    auto LayerNbr = b[i]->AsInt("layer");
    
    // Wire direction
    auto DirectionWire = b[i]->AsString("anodedir");
    //Wire id in the layer
    auto WireNbr  = b[i]->AsInt("wireid");

    //Wire x or y position
    auto PositionXY = b[i]->AsDouble("wirepos");
    auto PositionZ = b[i]->AsDouble("wirez");

    // position wire and direction wire
    if(DirectionWire.find("Y")==std::string::npos){ // if not horizontal
      m_PositionWire[ID] = G4ThreeVector(PositionXY,0.0,PositionZ);
      m_IsHorizontal[ID] = false;
    }
    else { // if horizontal
      m_PositionWire[ID] = G4ThreeVector(0.0,PositionXY,PositionZ);
      m_IsHorizontal[ID] = true;
    }

    m_LayerNbr[ID] = LayerNbr;
    m_WireNbr[ID] = WireNbr;
        
  } //end of for loop
  cout << " -> " << NumberOfWires << " wires found for FDC0" << endl;

}//end of ReadXML

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Construct detector and inialise sensitive part.
// (Called After DetectorConstruction::AddDetector Method)
void SamuraiFDC0::ConstructDetector(G4LogicalVolume* world){

  G4RotationMatrix* Rot = new G4RotationMatrix();
  
  new G4PVPlacement(0, m_Pos,
          BuildFDC0(), "SamuraiFDC0", world, false, 0);
  
  for(auto pos : m_PositionWire){ //loop on FDC0 wires
    
    if(m_IsHorizontal[pos.first]){
	Rot->rotateY(90.*deg);

	
      }
    else Rot->rotateX(90.*deg);
	
    new G4PVPlacement(G4Transform3D(*Rot,pos.second),BuildWire(),
		      "WireModule",m_FDC0,false,pos.first);
    Rot->set(0,0,0);
  }// end for loop on wires
  
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SamuraiFDC0::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  
  if(!pTree->FindBranch("SamuraiFDC0")){
    pTree->Branch("SamuraiFDC0", "TSamuraiFDC0Data", &m_Event) ;
  }
  pTree->SetBranchAddress("SamuraiFDC0", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// (Called at in the EventAction::EndOfEventAvtion)
void SamuraiFDC0::ReadSensitive(const G4Event* event){
  
  m_Event->Clear();

  WireScorers::PS_Wire* WireScorerFDC0= (WireScorers::PS_Wire*) m_WireScorerFDC0->GetPrimitive(0);

  unsigned int sizeWire =WireScorerFDC0->GetMult();
  for(unsigned int i = 0 ; i < sizeWire ; i++)
    {
      int layer = WireScorerFDC0->GetLayerNumber(i);
      int wireNbr = WireScorerFDC0->GetWireNumber(i);
      double time = WireScorerFDC0->GetTime(i);
      int edge = WireScorerFDC0->GetEdge(i);
      double DriftLength = WireScorerFDC0->GetDriftLength(i);
      
      m_Event->SetData(SamuraiFDC0_NS::FDC0_DetectorNumber, layer,
		       wireNbr, DriftLength, edge);// driftlength instead of time
    }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  
void SamuraiFDC0::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_WireScorerFDC0 = CheckScorer("WireScorerFDC0",already_exist);
  if(already_exist) 
    return;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  //Wire
  G4VPrimitiveScorer* InteractionWire= new WireScorers::PS_Wire("WS_FDC0",level, SamuraiFDC0_NS::Number_Of_Layer,SamuraiFDC0_NS::Number_Of_Wire_By_Layer,SamuraiFDC0_NS::Drift_Speed, m_WireAngle) ;
  m_WireScorerFDC0->RegisterPrimitive(InteractionWire);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_WireScorerFDC0) ;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* SamuraiFDC0::Construct(){
  return  (NPS::VDetector*) new SamuraiFDC0();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_samuraiFDC0{
    public:
      proxy_nps_samuraiFDC0(){
        NPS::DetectorFactory::getInstance()->AddToken("SAMURAIFDC0","SAMURAIFDC0");
        NPS::DetectorFactory::getInstance()->AddDetector("SAMURAIFDC0",SamuraiFDC0::Construct);
      }
  };
  
  proxy_nps_samuraiFDC0 p_nps_samuraiFDC0;
}


