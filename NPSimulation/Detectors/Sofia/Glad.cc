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
 *  This class describe a simple Glad setup for simulation                   *
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
#include "G4UserLimits.hh"

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
#include "G4FastSimulationManager.hh"
#include "G4Region.hh"

// NPTool header
#include "Glad.hh"
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
namespace Glad_NS{
  const double GLAD_height = 2000*mm; // y
  const double GLAD_width  = 4000*mm; // x
  const double GLAD_Depth  = 6000*mm; //z

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Glad Specific Method
Glad::Glad(){
  m_Magnet= NULL;
  m_GLAD_STL= NULL;
  m_PropagationRegion= NULL;
  m_GladScorer= NULL;
  m_GLAD_TiltAngle= -14*deg;
  m_Current= 2185;
  m_StepSize= 50*mm;

  m_GLAD_X = 0;
  m_GLAD_Y = 0;
  m_GLAD_Z = -1113.5*mm;

  m_MW3_R = 4440*mm;

  // RGB Color + Transparency
  m_VisSquare = new G4VisAttributes(G4Colour(0.53, 0.81, 0.98, 1));   
  m_VisGLAD = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.9));   
  m_VisField = new G4VisAttributes(G4Colour(0.1, 0.6, 0.9, 0.5));   
  m_VisKapton = new G4VisAttributes(G4Colour(1, 0.4, 0., 0.6));   

}

Glad::~Glad(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Glad::AddMagnet(G4ThreeVector POS, string fieldmap){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_X = POS.x();
  m_Y = POS.y();
  m_Z = POS.z();
  m_R = POS.mag();
  m_Theta = POS.theta();
  m_Phi = POS.phi();

  m_FieldMapFile = fieldmap;

  m_MW3_POS = TVector3(-m_MW3_R*sin(20*deg),0,m_MW3_R*cos(20*deg));
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Glad::BuildGLADFromSTL()
{
  if(!m_GLAD_STL){
    string basepath = getenv("NPTOOL");
    string path = basepath + "/NPSimulation/Detectors/Sofia/stl/GLAD_only.stl";

    auto mesh = CADMesh::TessellatedMesh::FromSTL((char*) path.c_str());
    mesh->SetScale(mm);

    G4Material* GLAD_Material = MaterialManager::getInstance()->GetMaterialFromLibrary("Inox");

    auto cad_solid = mesh->GetSolid();
    m_GLAD_STL = new G4LogicalVolume(cad_solid,GLAD_Material,"GLAD_Magnet",0,0,0);

    m_GLAD_STL->SetVisAttributes(m_VisGLAD);
  }

  return m_GLAD_STL;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Glad::BuildMagnet()
{
  if(!m_Magnet){
    // Shape - Box
    G4Box* box = new G4Box("glad_Box",Glad_NS::GLAD_width*0.5,Glad_NS::GLAD_height*0.5,Glad_NS::GLAD_Depth*0.5);

    // Material
    G4Material* vac = MaterialManager::getInstance()->GetMaterialFromLibrary("Vaccuum");

    // Logical Volume
    m_Magnet = new G4LogicalVolume(box,vac,"logic_GLAD_box",0,0,0);
    m_Magnet->SetVisAttributes(m_VisField);
  }
  return m_Magnet;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Glad::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Glad");

  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// GLAD Magnet found ////" << endl; 

  vector<string> cart = {"POS","TILT_ANGLE","FIELDMAP","CURRENT","STEPSIZE"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Glad information: " <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      m_GLAD_TiltAngle = blocks[i]->GetDouble("TILT_ANGLE","deg");
      string fieldmap = blocks[i]->GetString("FIELDMAP");
      m_Current = blocks[i]->GetDouble("CURRENT","A");
      m_StepSize = blocks[i]->GetDouble("STEPSIZE","mm");
      m_MW3_R = blocks[i]->GetDouble("MW3_DISTANCE","mm");
      m_MW3_CentralTheta = blocks[i]->GetDouble("MW3_CENTRAL_THETA","deg");
      
      AddMagnet(Pos,fieldmap);
    }
    else{
      cout << "ERROR: check your GLAD input file formatting " << endl;
      exit(1);
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Glad::ConstructDetector(G4LogicalVolume* world){

  //G4double wX = m_R * sin(m_Theta) * cos(m_Phi);
  //G4double wY = m_R * sin(m_Theta) * sin(m_Phi);
  //G4double wZ = m_R * cos(m_Theta);
  G4double wX = m_X;// - 1500*sin(m_GLAD_TiltAngle);
  G4double wY = m_Y;
  G4double wZ = m_Z;// - 1500*(1-cos(m_GLAD_TiltAngle));
  G4ThreeVector Mag_pos = G4ThreeVector(wX,wY,wZ);

  G4ThreeVector u (cos(m_GLAD_TiltAngle), 0, -sin(m_GLAD_TiltAngle));
  G4ThreeVector v (0,1,0);
  G4ThreeVector w (sin(m_GLAD_TiltAngle), 0, cos(m_GLAD_TiltAngle));
  //G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
  G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);

  new G4PVPlacement(Rot, Mag_pos,
      BuildMagnet(), "Glad", world, false, 0);

  SetPropagationRegion();

  /*G4RotationMatrix* RotGLAD = new G4RotationMatrix();
    RotGLAD->rotateX(90*deg);
  //RotGLAD->rotateZ(-(90-m_GLAD_TiltAngle)*deg);
  RotGLAD->rotateZ(-90*deg);
  RotGLAD->rotateZ(m_GLAD_TiltAngle);
  new G4PVPlacement(RotGLAD, GLAD_pos,
  BuildGLADFromSTL(),
  "GLAD",
  world, false, 0);
   */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Glad::SetPropagationRegion()
{
  if(!m_PropagationRegion){
    m_PropagationRegion= new G4Region("NPGladFieldPropagation");
    //m_PropagationRegion -> AddRootLogicalVolume(BuildPropvol());
    m_PropagationRegion -> AddRootLogicalVolume(BuildMagnet());
    m_PropagationRegion->SetUserLimits(new G4UserLimits(m_StepSize));
  }

  G4FastSimulationManager* mng = m_PropagationRegion->GetFastSimulationManager();
  //To make sure no other models are present in the region
  unsigned int size = m_PropagationModel.size();
  for(unsigned int i = 0 ; i < size ; i++)
    mng->RemoveFastSimulationModel(m_PropagationModel[i]);
  m_PropagationModel.clear();

  G4VFastSimulationModel* fsm;
  fsm = new NPS::GladFieldPropagation("GladFieldPropagation", m_PropagationRegion);
  ((NPS::GladFieldPropagation*) fsm)->SetGladEntrance(m_GLAD_X, m_GLAD_Y, m_GLAD_Z);
  ((NPS::GladFieldPropagation*) fsm)->SetGladTiltAngle(m_GLAD_TiltAngle);
  ((NPS::GladFieldPropagation*) fsm)->SetCurrent(m_Current);
  ((NPS::GladFieldPropagation*) fsm)->SetStepSize(m_StepSize);
  ((NPS::GladFieldPropagation*) fsm)->SetFieldMap(m_FieldMapFile);
  ((NPS::GladFieldPropagation*) fsm)->Set_MWPC3_Position(m_MW3_POS.X(),m_MW3_POS.Y(),m_MW3_POS.Z());
  ((NPS::GladFieldPropagation*) fsm)->Set_MWPC3_CentralTheta(m_MW3_CentralTheta);
  /*if(m_Method == NPS::RungeKutta){
    double r_max = sqrt( 
    Glad_NS::GLAD_Width * Glad_NS::GLAD_Width /4. + 
    Glad_NS::GLAD_Depth * Glad_NS::GLAD_Depth /4. );//sqrt(x^2 + z^2)
    ((NPS::GladFieldPropagation*) fsm)->SetRmax(r_max);
    }*/
  m_PropagationModel.push_back(fsm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Glad::InitializeRootOutput(){
  //  RootOutput *pAnalysis = RootOutput::getInstance();
  //  TTree *pTree = pAnalysis->GetTree();
  //  if(!pTree->FindBranch("IdealGladData")){
  //    pTree->Branch("IdealGladData", "TGladIdealData", &m_Event) ;
  //  }
  //  pTree->SetBranchAddress("IdealGladData", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Glad::ReadSensitive(const G4Event* ){
  //m_Event->Clear();

  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_GladScorer->GetPrimitive(0);

  /*unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    double energy = Scorer->GetEnergy(i);
  }*/
  Scorer->clear();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Glad::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_GladScorer = CheckScorer("GladScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; 
  level.push_back(0);
  //G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  //and register it to the multifunctionnal detector
  //m_GladScorer->RegisterPrimitive(Interaction);
  m_GladScorer->RegisterPrimitive(Calorimeter);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_GladScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Glad::Construct(){
  return  (NPS::VDetector*) new Glad();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_glad{
    public:
      proxy_nps_glad(){
        NPS::DetectorFactory::getInstance()->AddToken("Glad","Glad");
        NPS::DetectorFactory::getInstance()->AddDetector("Glad",Glad::Construct);
      }
  };

  proxy_nps_glad p_nps_glad;
}
