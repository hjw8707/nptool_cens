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
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Samurai BDCs simulation                             *
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
#include "SamuraiBDC.hh"
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
namespace SamuraiBDC_NS{
  // Samurai magnet construction paramethers

  //Main outer box
  const double BDC_Width = 94*mm; //(x)
  const double BDC_Height = 94*mm;//(y)
  const double BDC_Depth = 90*mm;//(z)
  const string BDC_Gas = "CH4_60_He_40";

  const double BDC_Temperature = 298.15; // K
  const double BDC_Pressure = 1.0; //atm

  //Detector Number
  const short int BDC_DetectorNumber = 1;

  //Wires
  const double Wire_Length = 80*mm;
  const double Wire_Diameter = 4.8*mm;
  const string Wire_Gas = "CH4_60_He_40";
  const double Wire_Temperature = 298.15; // K 
  const double Wire_Pressure = 1.0; //atm
  const int Number_Of_Layer = 8; // from 0 to 7 in the xml file
  const int Number_Of_Wire_By_Layer = 16; // from 0 to 15 in the xml file
  const double Drift_Speed = 1.33e-4;//------> FIX ME!!! not important (yet)

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Samurai Specific Method
SamuraiBDC::SamuraiBDC(){

  //Visualization attributes
  m_VisBDC = new G4VisAttributes(G4Colour(1,1,0,0.5));
  //Logical volumes
  m_BDC1 = NULL;
  m_BDC2 = NULL;
 
 //Wire
  m_VisWire = new G4VisAttributes(G4Colour(1.0,0,0.0,1.0));
  //Logical volumes
  m_Wire = NULL;
  //Scorer
  m_WireScorerBDC = NULL;

  //Data event
  m_Event = new TSamuraiBDCData;

  //Hexagon
  m_Hexagon = NULL;
  
}

SamuraiBDC::~SamuraiBDC(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SamuraiBDC::AddDetector(G4ThreeVector Mag_Pos, G4ThreeVector Offset, unsigned int det){

    m_position[det] = Mag_Pos + Offset;
  
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* SamuraiBDC::BuildBDC1(){
  if(!m_BDC1){
    //Shape - G4Box
    G4Box* box = new G4Box("BDC_Box",SamuraiBDC_NS::BDC_Width*0.5,
			   SamuraiBDC_NS::BDC_Height*0.5,SamuraiBDC_NS::BDC_Depth*0.5);

    //Material - Gas
    G4Material* Gas = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiBDC_NS::BDC_Gas, SamuraiBDC_NS::BDC_Pressure,SamuraiBDC_NS::BDC_Temperature);

    //Logical Volume
    m_BDC1 = new G4LogicalVolume(box, Gas, "logic_SamuraiBDC_box",0,0,0);
    m_BDC1->SetVisAttributes(m_VisBDC);
  }
  return m_BDC1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* SamuraiBDC::BuildBDC2(){
  if(!m_BDC2){
    //Shape - G4Box
    G4Box* box = new G4Box("BDC_Box",SamuraiBDC_NS::BDC_Width*0.5,
			   SamuraiBDC_NS::BDC_Height*0.5,SamuraiBDC_NS::BDC_Depth*0.5);

    //Material - Gas
    G4Material* Gas = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiBDC_NS::BDC_Gas, SamuraiBDC_NS::BDC_Pressure,SamuraiBDC_NS::BDC_Temperature);

    //Logical Volume
    m_BDC2 = new G4LogicalVolume(box, Gas, "logic_SamuraiBDC_box",0,0,0);
    m_BDC2->SetVisAttributes(m_VisBDC);
  }
  return m_BDC2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiBDC::BuildWire(){
if(!m_Wire){
  //Shape - G4Tubs
  G4Tubs* solidWire = new G4Tubs("BDC_Wire", 0.0, SamuraiBDC_NS::Wire_Diameter*0.5,SamuraiBDC_NS::Wire_Length*0.5, 0., 360.);

  //Material - Gas
  G4Material* GasWire = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiBDC_NS::Wire_Gas, SamuraiBDC_NS::Wire_Pressure,SamuraiBDC_NS::Wire_Temperature);
  
  //Logical Volume
  m_Wire = new G4LogicalVolume(solidWire, GasWire, "logic_SamuraiWire_tub",0,0,0);
  m_Wire->SetVisAttributes(m_VisWire);
  m_Wire->SetSensitiveDetector(m_WireScorerBDC);
  
  }
  return m_Wire;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* SamuraiBDC::BuildHexagon(){
if(!m_Hexagon){
  //Shape 
 
  const G4int nsect = 6;
  std::vector<G4TwoVector> polygon(nsect);
  G4double ang = twopi/nsect;

  G4double dz = SamuraiBDC_NS::Wire_Length*0.5;
  G4double rmax = SamuraiBDC_NS::Wire_Diameter*0.5;
  
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
  G4Material* GasWire = MaterialManager::getInstance()->GetGasFromLibrary(SamuraiBDC_NS::Wire_Gas, SamuraiBDC_NS::Wire_Pressure,SamuraiBDC_NS::Wire_Temperature);
  
  //Logical Volume
  m_Hexagon = new G4LogicalVolume(solidHexagon, GasWire, "logic_SamuraiWire_tub",0,0,0);
  m_Hexagon->SetVisAttributes(m_VisWire);
  m_Hexagon->SetSensitiveDetector(m_WireScorerBDC);
  
  }
  return m_Hexagon;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// (Called in DetecorConstruction::ReadDetectorConfiguration Method)
void SamuraiBDC::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SAMURAIBDC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " BDC detector(s) found " << endl;
  
  vector<string> xml = {"XML","Offset", "InvertX","InvertY","InvertD"};

  G4ThreeVector Mag_Pos;
  G4ThreeVector Offset;
  
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(xml))
      {
	cout << endl << "////  Samurai BDC (" << i+1 << ")" << endl;
	unsigned int det = std::atoi(blocks[i]->GetMainValue().c_str());//1 or 2
	string xml_file = blocks[i]->GetString("XML");
	Offset = NPS::ConvertVector(blocks[i]->GetTVector3("Offset", "mm"));
	bool invert_x = blocks[i]->GetInt("InvertX"); 
	bool invert_y = blocks[i]->GetInt("InvertY"); 
	bool invert_d = blocks[i]->GetInt("InvertD"); 
	m_offset[det] = Offset;
	m_invertX[det] = invert_x;
	m_invertY[det] = invert_y;
	m_invertD[det] = invert_d;
	ReadXML(xml_file,Offset,invert_y,invert_y,det);
	AddDetector(Mag_Pos, Offset, det);
      }

    else{
      cout << "ERROR: check your XML input file for BDC detector" << endl;
      exit(1);
    }//end exception  

  }
     
}//end read configuration


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//////////////////// nouveau ///////////////////
void SamuraiBDC::ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY,unsigned int det){ 
  NPL::XmlParser xml;
  xml.LoadFile(xml_file);

  std::string name = "SAMURAIBDC"+NPL::itoa(det);
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName(name);

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

    if(det==1)
      {
	// position wire and direction wire for BDC1
	if(DirectionWire.find("Y")==std::string::npos){ // if not horizontal
	  m_PositionWire1[ID] = G4ThreeVector(PositionXY,0.0,PositionZ);
	  m_IsHorizontal1[ID] = false;
	}
	else {
	  m_PositionWire1[ID] = G4ThreeVector(0.0,PositionXY,PositionZ);
	  m_IsHorizontal1[ID] = true;
	}
	
	m_LayerNbr1[ID] = LayerNbr;
	m_WireNbr1[ID] = WireNbr;
      }//end det==1

    else if(det==2)
      {
	// position wire and direction wire for BDC2
	if(DirectionWire.find("Y")==std::string::npos){ // if not horizontal
	  m_PositionWire2[ID] = G4ThreeVector(PositionXY,0.0,PositionZ);
	  m_IsHorizontal2[ID] = false;
	}
	else {
	  m_PositionWire2[ID] = G4ThreeVector(0.0,PositionXY,PositionZ);
	  m_IsHorizontal2[ID] = true;
	}
	
	m_LayerNbr2[ID] = LayerNbr;
	m_WireNbr2[ID] = WireNbr;

      }//end det==2
    
  } //end of for loop
  cout << " -> " << NumberOfWires << " wires found for BDC" <<det << endl;

}//end of ReadXML


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// (Called After DetectorConstruction::AddDetector Method)
void SamuraiBDC::ConstructDetector(G4LogicalVolume* world){

  G4RotationMatrix* Rot = new G4RotationMatrix();
    
  new G4PVPlacement(0, m_position[1],
		    BuildBDC1(), "SAMURAIBDC1" , world, false, 1);

  new G4PVPlacement(0, m_position[2],
		    BuildBDC2(), "SAMURAIBDC2" , world, false, 2);
  
  for(auto pos1 : m_PositionWire1){    
    if(m_IsHorizontal1[pos1.first])
      {
	Rot->rotateY(90.*deg);
      }
    else {
      Rot->rotateX(90.*deg);
      }

    new G4PVPlacement(G4Transform3D(*Rot,pos1.second),
            BuildWire(),
		      "WireModule", m_BDC1,false,pos1.first);
    Rot->set(0,0,0);
  }// end for loop on wires BDC1

  for(auto pos2 : m_PositionWire2){  
    if(m_IsHorizontal2[pos2.first])
      {
	Rot->rotateY(90.*deg);
      }
    else Rot->rotateX(90.*deg);
    	
    new G4PVPlacement(G4Transform3D(*Rot,pos2.second),
            BuildWire(),
            "WireModule", m_BDC2 ,false,pos2.first);
    Rot->set(0,0,0);
  }// end for loop on wires BDC2
  
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void SamuraiBDC::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("SamuraiBDC")){
    pTree->Branch("SamuraiBDC", "TSamuraiBDCData", &m_Event) ;
  }
  pTree->SetBranchAddress("SamuraiBDC", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Read sensitive part and fill the Root tree.
// (Called at in the EventAction::EndOfEventAvtion)

void SamuraiBDC::ReadSensitive(const G4Event* event){
  
  m_Event->Clear();
  //Interaction Scorer
  WireScorers::PS_Wire* WireScorerBDC= (WireScorers::PS_Wire*) m_WireScorerBDC->GetPrimitive(0);

  unsigned int sizeWire =WireScorerBDC->GetMult();
  for(unsigned int i = 0 ; i < sizeWire ; i++)
    {
      int layer = WireScorerBDC->GetLayerNumber(i);
      int wireNbr = WireScorerBDC->GetWireNumber(i);
      double time = WireScorerBDC->GetTime(i);
      int edge = WireScorerBDC->GetEdge(i);
      double DriftLength = WireScorerBDC->GetDriftLength(i);
      int DetectorNumber = WireScorerBDC->GetDetectorNumber(i);    
      m_Event->SetData(DetectorNumber, layer, wireNbr, DriftLength, edge);
    }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  
void SamuraiBDC::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_WireScorerBDC = CheckScorer("WireScorerBDC",already_exist) ;
 
  if(already_exist) 
    return;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  
  //Wire
  G4VPrimitiveScorer* InteractionWireBDC= new WireScorers::PS_Wire("WS_BDC",level, SamuraiBDC_NS::Number_Of_Layer,SamuraiBDC_NS::Number_Of_Wire_By_Layer,SamuraiBDC_NS::Drift_Speed, m_WireAngle1) ;
  m_WireScorerBDC->RegisterPrimitive(InteractionWireBDC);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_WireScorerBDC) ;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////

NPS::VDetector* SamuraiBDC::Construct(){
  return  (NPS::VDetector*) new SamuraiBDC();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_samuraiBDC{
    public:
      proxy_nps_samuraiBDC(){
        NPS::DetectorFactory::getInstance()->AddToken("SAMURAIBDC","SAMURAIBDC");
        NPS::DetectorFactory::getInstance()->AddDetector("SAMURAIBDC",SamuraiBDC::Construct);
      }
  };
  
  proxy_nps_samuraiBDC p_nps_samuraiBDC;
}


