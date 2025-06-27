/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto  contact address: pilottoelia@gmail.com     *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    : septembre 2021                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Samurai simulation                                 *
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
#include "G4Trap.hh"
#include "G4UnionSolid.hh"
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
#include "G4RegionStore.hh"
#include "G4FastSimulationManager.hh"
#include "G4UserLimits.hh"

// NPTool header
#include "Samurai.hh"
#include "SamuraiFieldPropagation.hh"
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
namespace Samurai_NS{
  // Samurai magnet construction paramethers

  //Main outer box
  const double Magnet_Width = 6700*mm; //(x)
  const double Magnet_Height = 4640*mm;//(y)
  const double Magnet_Depth = 4950*mm;//(z) //If depth is larger, it overlaps on FDC2
  const double Magnet_Depth_eps = 0.2*mm; //(delta-z)
  const string Magnet_Material = "G4_Galactic"; //G4_Galactic

  //Gap between top and bottom yoke slabs
  const double Magnet_Gap = 880*mm;//(y)

  //Iron yoke - top and bottom slabs
  const double Yoke_Width = 6700*mm;//(x)
  const double Yoke_Height = ( Magnet_Height - Magnet_Gap ) / 2.;//(y)
  const double Yoke_Depth = 3500*mm;//(z)
  const string Yoke_Material = "G4_Fe";

  //Iron yoke - left and right pillars
  //Gap inside yoke pillars
  const double RYHoriz_Gap = 400*mm;//(z)
  //Box
  const double RYBox_Width = 750.*mm;//(x)
  const double RYBox_Height = Magnet_Gap;//(y)
  const double RYBox_Depth = (Yoke_Depth - RYHoriz_Gap) / 2.; //1550*mm(z)
  //Trap
  const double RYTrap_Height = RYBox_Width; //(x)
  const double RYTrap_Depth = RYBox_Height; //(y)
  const double RYTrap_Width = RYBox_Depth; //(lower_z)
  const double RYTrap_Width_LTX = RYBox_Depth*0.4; //(upper_z)
  const string RYoke_Material = "G4_Fe";

  //Propagation volume
  const double pv_main_Width = Magnet_Width;
  const double pv_main_Height = Magnet_Gap;
  const double pv_main_Depth = Magnet_Depth;
  const string Propvol_Material = "G4_Galactic";

  //Ideal Data saving location
  const short int Magnet_DetectorNumber = -1;

} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Samurai Specific Method
Samurai::Samurai(){

  //Visualization attributes

  m_VisMagnet = new G4VisAttributes(G4Colour(0,0,1,0.2));
  m_VisYokes = new G4VisAttributes(G4Colour(0,1,0,0.2));
  m_VisRYokes = new G4VisAttributes(G4Colour(1,0,0,0.2));
  m_VisPropvol = new G4VisAttributes(G4Colour(0,1,1,0.2));

  // m_VisMagnet = new G4VisAttributes(G4Colour(0,0,1,0.0));
  // m_VisYokes = new G4VisAttributes(G4Colour(0,1,0,0.0));
  // m_VisRYokes = new G4VisAttributes(G4Colour(1,0,0,0.0));
  // m_VisPropvol = new G4VisAttributes(G4Colour(0,1,1,0.0));
  
  //Logical volumes
  m_Magnet = NULL;
  m_Yoke = NULL;
  m_RYoke = NULL;
  m_Propvol = NULL;
  m_RYokeSolid = NULL;

  //Propagation region
  m_PropagationRegion = NULL;

  //Scorer
  m_SamuraiScorer = NULL;

  //Data
  m_Event = new TSamuraiIdealData;
}

Samurai::~Samurai(){
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Samurai::AddMagnet(G4ThreeVector POS, double Angle, int method, string fieldmap,
			double n_steps){

  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R = POS.mag();
  m_Theta = POS.theta();
  m_Phi = POS.phi();
  m_Angle = Angle;

  switch (method) {
    case 0:
      m_Method = NPS::RungeKutta;
      break;
    case 1:
      m_Method = NPS::EliaOmar;//FIXME
      break;
    default:
      cout << endl;
      cout << "//////WARNING///////" << endl;
      cout << "In Samurai.cc:\n";
      cout << "Propagation method " << method << " not defined\n";
      cout << endl;
    }

  m_FieldMapFile = fieldmap;
  m_StepSize = 1.*meter / n_steps;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Samurai::AddMagnet(double R, double Theta, double Phi, double Angle, int method,
			string fieldmap, double n_steps){

  m_R = R;
  m_Theta = Theta;
  m_Phi = Phi;
  m_Angle = Angle;

  switch (method) {
    case 0:
      m_Method = NPS::RungeKutta;
      break;
    case 1:
      m_Method = NPS::EliaOmar;//FIXME
      break;
    default:
      cout << endl;
      cout << "//////WARNING///////" << endl;
      cout << "In Samurai.cc:\n";
      cout << "Propagation method " << method << " not defined\n";
      cout << endl;
    }

  m_FieldMapFile = fieldmap;
  m_StepSize = 1.* meter / n_steps;
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Samurai::BuildMagnet(){
  if(!m_Magnet){
    //Shape - G4Box
    G4Box* box = new G4Box("Samurai_Box",Samurai_NS::Magnet_Width*0.5,
			   Samurai_NS::Magnet_Height*0.5,Samurai_NS::Magnet_Depth*0.5 + Samurai_NS::Magnet_Depth_eps*0.5);
  
    //Material - vacuum
    G4Material* VacuumMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(Samurai_NS::Magnet_Material);

    //Logical Volume
    m_Magnet = new G4LogicalVolume(box, VacuumMaterial, "logic_Samurai_box",0,0,0);
    m_Magnet->SetVisAttributes(m_VisMagnet);
    m_Magnet->SetSensitiveDetector(m_SamuraiScorer);
    cout << "ATTACHING SCORER\n";
    
  }
  return m_Magnet;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Samurai::BuildYoke(){
  if(!m_Yoke){
    //Shape - G4Box
    G4Box* yoke = new G4Box("Samurai_yoke_1",Samurai_NS::Yoke_Width*0.5,
			    Samurai_NS::Yoke_Height*0.5,Samurai_NS::Yoke_Depth*0.5);

    //Material
    G4Material* YokeMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(Samurai_NS::Yoke_Material);

    //Logical Volume
    m_Yoke = new G4LogicalVolume(yoke, YokeMaterial, "logic_Samurai_yoke",0,0,0);
    m_Yoke->SetVisAttributes(m_VisYokes);
  
  }
  return m_Yoke;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VSolid* Samurai::BuildRYokeSolid(){
  if(!m_RYokeSolid){
    //Shape - G4Box
    G4Box* ryokebox = new G4Box("Samurai_ryoke_box",Samurai_NS::RYBox_Width*0.5,
			    Samurai_NS::RYBox_Height*0.5,Samurai_NS::RYBox_Depth*0.5);
    //Shape -G4Trap
    G4Trap* ryoketrap = new G4Trap("Samurai_ryoke_trap", Samurai_NS::RYTrap_Depth, Samurai_NS::RYTrap_Height,
                                    Samurai_NS::RYTrap_Width, Samurai_NS::RYTrap_Width_LTX);
                      
    //Union
    G4ThreeVector rytrans( 
      0.5*(Samurai_NS::RYBox_Width + Samurai_NS::RYTrap_Height),
      0.,
      0.25*(Samurai_NS::RYTrap_Width - Samurai_NS::RYTrap_Width_LTX));
    G4RotationMatrix* traprot = new G4RotationMatrix();
    traprot->rotateX(-90. *deg);
    traprot->rotateZ(90.*deg);

    double displacement_z = Samurai_NS::RYBox_Depth + Samurai_NS::RYHoriz_Gap;
    G4RotationMatrix* rotation1 = new G4RotationMatrix();
    rotation1->rotateY(180.*deg);
    rotation1->rotateZ(180.*deg);
    G4RotationMatrix* rotation2 = new G4RotationMatrix();
    rotation2->rotateY(180.*deg);
    G4UnionSolid* ryoke1 = new G4UnionSolid("Samurai_ryoke_1", ryokebox, ryoketrap, traprot, rytrans);
    G4UnionSolid* ryoke2 = new G4UnionSolid("Samurai_ryoke_2", ryoke1, ryoke1, rotation1, G4ThreeVector(0,0,displacement_z));
    m_RYokeSolid = new G4UnionSolid("solid_Samurai_ryoke", ryoke2, ryoke2, rotation2, G4ThreeVector(Samurai_NS::Yoke_Width - (Samurai_NS::RYBox_Width),0.,displacement_z));

  }
  return m_RYokeSolid;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Samurai::BuildRYoke(){
  if(!m_RYoke){

    //Material
    G4Material* RYokeMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(Samurai_NS::RYoke_Material);

    //Logical Volume
    m_RYoke = new G4LogicalVolume(BuildRYokeSolid(), RYokeMaterial, "logic_Samurai_ryoke",0,0,0);
    m_RYoke->SetVisAttributes(m_VisRYokes);
  }
  return m_RYoke;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Samurai::BuildPropvol(){
  if(!m_Propvol){
    //Shape - G4Box (Main)
    G4Box* pvmain = new G4Box(
      "Samurai_propvol_box",
      Samurai_NS::pv_main_Width*0.5,
      Samurai_NS::pv_main_Height*0.5,
      Samurai_NS::pv_main_Depth*0.5   );
    //Subtraction1
    double displacement_x = (Samurai_NS::Yoke_Width - Samurai_NS::RYBox_Width) / 2.0;
    double displacement_z = (Samurai_NS::Yoke_Depth - Samurai_NS::RYBox_Depth) / 2.0;
    G4ThreeVector trans1 (-displacement_x,0., -displacement_z);
    G4SubtractionSolid* propvol = new G4SubtractionSolid("solid_Samurai_propvol", pvmain, BuildRYokeSolid(), 0, trans1);

    //Material
    G4Material* PropVolMaterial = MaterialManager::getInstance()
      ->GetMaterialFromLibrary(Samurai_NS::Propvol_Material);

    //Logical Volume
    m_Propvol = new G4LogicalVolume(propvol, PropVolMaterial, "logic_Samurai_propvol",0,0,0);
    m_Propvol->SetVisAttributes(m_VisPropvol);
  }
  return m_Propvol;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetectorConfiguration Method
void Samurai::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Samurai");

  if(blocks.size()==1){
    if(NPOptionManager::getInstance()->GetVerboseLevel()) {
      cout << "/////// Samurai magnet found ///////" << endl;
    }
    vector<string> cart = {"POS","ANGLE","METHOD","FIELDMAP","STEPS_PER_METER"};
    vector<string> sphe = {"R","Theta","Phi","ANGLE","METHOD","FIELDMAP","STEPS_PER_METER"};

    if(blocks[0]->HasTokenList(cart)){
      G4ThreeVector Pos = NPS::ConvertVector(blocks[0]->GetTVector3("POS", "cm"));
      double Angle = blocks[0]->GetDouble("ANGLE","deg");
      int method = blocks[0]->GetInt("METHOD");
      string fieldmap = blocks[0]->GetString("FIELDMAP");
      double n_steps = (double) blocks[0]->GetInt("STEPS_PER_METER");
      AddMagnet(Pos,Angle,method,fieldmap,n_steps);
    }
    else if(blocks[0]->HasTokenList(sphe)){
      double R = blocks[0]->GetDouble("R","mm");
      double Theta = blocks[0]->GetDouble("Theta","deg");
      double Phi = blocks[0]->GetDouble("Phi","deg");
      double Angle = blocks[0]->GetDouble("ANGLE","deg");
      int method = blocks[0]->GetInt("METHOD");
      string fieldmap = blocks[0]->GetString("FIELDMAP");
      double n_steps = (double) blocks[0]->GetInt("STEPS_PER_METER");
      AddMagnet(R,Theta,Phi,Angle,method,fieldmap,n_steps);
    }
  }
  else{
    cout << "ERROR: there should be only one Samurai magnet, check your input file" << endl;
    exit(1);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetectorConstruction::AddDetector Method
void Samurai::ConstructDetector(G4LogicalVolume* world){

  //Yokes placement
  double distance = ( Samurai_NS::Magnet_Height - Samurai_NS::Yoke_Height ) * 0.5;
  G4ThreeVector upper (0, distance,0);
  G4ThreeVector lower (0,-distance,0);
  
  new G4PVPlacement(0, upper, BuildYoke(), "Upper_Yoke", BuildMagnet(), false, 0);
  new G4PVPlacement(0, lower, BuildYoke(), "Lower_Yoke", BuildMagnet(), false, 0);
  
  //RYokes placement
  double displacement_x = (Samurai_NS::Yoke_Width - Samurai_NS::RYBox_Width) / 2.0;
  double displacement_z = (Samurai_NS::Yoke_Depth - Samurai_NS::RYBox_Depth) / 2.0;
  G4ThreeVector trans1 (-displacement_x,0., -displacement_z);

  new G4PVPlacement(0, trans1, BuildRYoke(), "R_Yoke1", BuildMagnet(), false, 0);

  //Propagation Volume placement
  new G4PVPlacement(0, G4ThreeVector(0,0,0), BuildPropvol(), "Propagation_Volume", BuildMagnet(), false, 0);

  //Magnet placement
  G4double wX = m_R * sin(m_Theta) * cos(m_Phi);
  G4double wY = m_R * sin(m_Theta) * sin(m_Phi);
  G4double wZ = m_R * cos(m_Theta);
  G4ThreeVector Mag_pos = G4ThreeVector(wX, wY, wZ);
  G4ThreeVector u ( cos(m_Angle) , 0, -sin(m_Angle) );
  G4ThreeVector v (0,1,0);
  G4ThreeVector w ( sin(m_Angle) , 0,  cos(m_Angle) );
  
  G4RotationMatrix* Rot = new G4RotationMatrix(u,v,w);
  
  new G4PVPlacement(Rot, Mag_pos,
          BuildMagnet(), "Samurai", world, false, 0);

  // Set region were magnetic field is active
  SetPropagationRegion();

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Set region were magnetic field is active
// Called in Samurai::ConstructDetector
void Samurai::SetPropagationRegion(){
  
  if(!m_PropagationRegion){
    m_PropagationRegion= new G4Region("NPSamuraiFieldPropagation");
    m_PropagationRegion -> AddRootLogicalVolume(BuildPropvol());
    m_PropagationRegion->SetUserLimits(new G4UserLimits(m_StepSize));    
  }
  
  G4FastSimulationManager* mng = m_PropagationRegion->GetFastSimulationManager();
  //To make sure no other models are present in the region
  unsigned int size = m_PropagationModel.size();
  for(unsigned int i = 0 ; i < size ; i++)
      mng->RemoveFastSimulationModel(m_PropagationModel[i]);
  m_PropagationModel.clear();
 
  G4VFastSimulationModel* fsm;
  fsm = new NPS::SamuraiFieldPropagation("SamuraiFieldPropagation", m_PropagationRegion);
  ((NPS::SamuraiFieldPropagation*) fsm)->SetStepSize(m_StepSize);
  ((NPS::SamuraiFieldPropagation*) fsm)->SetFieldMap(m_FieldMapFile);
  ((NPS::SamuraiFieldPropagation*) fsm)->SetMethod(m_Method);
  if(m_Method == NPS::RungeKutta){
    double r_max = sqrt( 
      Samurai_NS::Magnet_Width * Samurai_NS::Magnet_Width /4. + 
      Samurai_NS::Magnet_Depth * Samurai_NS::Magnet_Depth /4. );//sqrt(x^2 + z^2)
    ((NPS::SamuraiFieldPropagation*) fsm)->SetRmax(r_max);
  }
  m_PropagationModel.push_back(fsm);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Samurai::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("IdealDataMagnet")){
    pTree->Branch("IdealDataMagnet", "TSamuraiIdealData", &m_Event) ;
  }
  pTree->SetBranchAddress("IdealDataMagnet", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Samurai::ReadSensitive(const G4Event* ){

  //cout << "Read Sensitive\n";
  m_Event->Clear();
  //Interaction Scorer
  InteractionScorers::PS_Interactions* Scorer= (InteractionScorers::PS_Interactions*) m_SamuraiScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  //cout << "size " << size << endl;
  for(unsigned int i = 0 ; i < size ; i++){
    double energy = Scorer->GetEnergy(i);
    double brho = Scorer->GetBrho(i);
    double posx = Scorer->GetPositionX(i);
    double posy = Scorer->GetPositionY(i);
    double posz = Scorer->GetPositionZ(i);
    double mom_mag = brho*Scorer->GetCharge(i);
    double theta = Scorer->GetTheta(i);
    double phi = Scorer->GetPhi(i);
    m_Event->SetData(Samurai_NS::Magnet_DetectorNumber, energy, posx, posy, posz, mom_mag, theta, phi, brho);
    //cout << Samurai_NS::Magnet_DetectorNumber << endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   

void Samurai::InitializeScorers() {
  cout << "I'M IN InitializeScorers" << endl;
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_SamuraiScorer = CheckScorer("SamuraiScorer",already_exist) ;

  if (already_exist) cout << "//////////////// already exist ///////// " << endl;
  else cout << "//////////////// does not already exist ///////// " << endl;
  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_SamuraiScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_SamuraiScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Samurai::Construct(){
  return  (NPS::VDetector*) new Samurai();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_samurai{
    public:
      proxy_nps_samurai(){
        NPS::DetectorFactory::getInstance()->AddToken("Samurai","Samurai");
        NPS::DetectorFactory::getInstance()->AddDetector("Samurai",Samurai::Construct);
      }
  };
  
  proxy_nps_samurai p_nps_samurai;
}


