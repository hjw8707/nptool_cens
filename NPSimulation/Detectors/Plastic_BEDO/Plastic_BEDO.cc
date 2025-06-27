/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Léo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Plastic_BEDO simulation                             *
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

// NPTool header
#include "Plastic_BEDO.hh"
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
namespace Plastic_BEDO_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.*MeV;
  const double ResoTime = 0.*ns ;
  const double ResoEnergy = 0.*MeV ;
  const double Radius = 25.4*mm ; 
  const double Thickness = 3*mm ;
  const double Length = 98*mm ;
  const string Material = "BC400";
  const double Alu_Thickness = 50*um;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Plastic_BEDO Specific Method
Plastic_BEDO::Plastic_BEDO(){
  m_Event = new TPlastic_BEDOData() ;
  m_Plastic_BEDOScorer = 0;
  m_CylindricalDetector = 0;


  // RGB Color + Transparency
  m_VisPlastic = new G4VisAttributes(G4Colour(0, 1, 0));   

}

Plastic_BEDO::~Plastic_BEDO(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Plastic_BEDO::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Plastic_BEDO::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* Plastic_BEDO::BuildCylindricalDetector(G4int DetNumber, G4ThreeVector Det_pos, G4RotationMatrix* Det_rot, G4LogicalVolume* world){

  G4Material* m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* m_MaterialAluminium = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Plastic_BEDO_NS::Material);
  G4VisAttributes* light_GreyAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.7));
  G4ThreeVector Pos0 = G4ThreeVector(0, 0, 0);

  // Global volume
  G4Tubs* solidPlastic = new G4Tubs("solidPlastic", 0, Plastic_BEDO_NS::Radius,Plastic_BEDO_NS::Length*0.5,0,360*deg);
  G4LogicalVolume* logicPlastic = new G4LogicalVolume(solidPlastic, m_MaterialVacuum, "logicPlastic", 0, 0);
  new G4PVPlacement(G4Transform3D(*Det_rot, Det_pos), logicPlastic, "Plastic", world, false, DetNumber);
  logicPlastic->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4ThreeVector dzAlu(0., 0., Plastic_BEDO_NS::Alu_Thickness);
  G4Tubs* Al_tubout = new G4Tubs("Al_tubout", 0, Plastic_BEDO_NS::Radius-Plastic_BEDO_NS::Thickness,(Plastic_BEDO_NS::Length-Plastic_BEDO_NS::Thickness)*0.5, 0, 360*deg);
  G4Tubs* Al_tubin = new G4Tubs("Al_tubin", 0, Plastic_BEDO_NS::Radius-Plastic_BEDO_NS::Alu_Thickness-Plastic_BEDO_NS::Thickness, (Plastic_BEDO_NS::Length-Plastic_BEDO_NS::Thickness)*0.5, 0, 360*deg);
  G4SubtractionSolid* Al_tub = new G4SubtractionSolid("Al_tub", Al_tubout, Al_tubin, 0, dzAlu);
  G4LogicalVolume* Al_tub_log = new G4LogicalVolume(Al_tub, m_MaterialAluminium, "Al_tub_log", 0, 0);
  G4ThreeVector PosAlu = G4ThreeVector(0, 0, Plastic_BEDO_NS::Thickness*0.5);
  new G4PVPlacement(0, PosAlu, Al_tub_log, "Al_tub", logicPlastic, false, DetNumber);
  Al_tub_log->SetVisAttributes(light_GreyAtt);

  // Substract two tubs to have a semi-closed plastic
  G4ThreeVector zTrans(0., 0., Plastic_BEDO_NS::Thickness);
  G4Tubs* tubOut = new G4Tubs("Plastic_BEDO",0,Plastic_BEDO_NS::Radius,Plastic_BEDO_NS::Length*0.5,0,360*deg);
  G4Tubs* tubIn = new G4Tubs("Plastic_BEDO",0,Plastic_BEDO_NS::Radius-Plastic_BEDO_NS::Thickness,Plastic_BEDO_NS::Length*0.5,0,360*deg);
  G4SubtractionSolid* tub = new G4SubtractionSolid("tub", tubOut, tubIn, 0, zTrans);
  m_CylindricalDetector = new G4LogicalVolume(tub,DetectorMaterial,"logic_Plastic_BEDO_tub",0,0,0);
  new G4PVPlacement(0, Pos0, m_CylindricalDetector, "Plastic_BEDO", logicPlastic, false, DetNumber);
  m_CylindricalDetector->SetVisAttributes(m_VisPlastic);
  
  m_CylindricalDetector->SetSensitiveDetector(m_Plastic_BEDOScorer);

  return m_CylindricalDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Plastic_BEDO::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Plastic_BEDO");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Plastic_BEDO " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Plastic_BEDO " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
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
void Plastic_BEDO::ConstructDetector(G4LogicalVolume* world){
  
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*Plastic_BEDO_NS::Thickness*0.5;
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
    
    G4RotationMatrix* Rot = new G4RotationMatrix(u, v, w);
   
    BuildCylindricalDetector(i+1, Det_pos, Rot, world);
    
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Plastic_BEDO::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Plastic_BEDO")){
    pTree->Branch("Plastic_BEDO", "TPlastic_BEDOData", &m_Event) ;
  }
  pTree->SetBranchAddress("Plastic_BEDO", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Plastic_BEDO::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_Plastic_BEDOScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Plastic_BEDO_NS::ResoEnergy);
    if(Energy>Plastic_BEDO_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Plastic_BEDO_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Plastic_BEDO::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_Plastic_BEDOScorer = CheckScorer("Plastic_BEDOScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_Plastic_BEDOScorer->RegisterPrimitive(Calorimeter);
  m_Plastic_BEDOScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_Plastic_BEDOScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Plastic_BEDO::Construct(){
  return  (NPS::VDetector*) new Plastic_BEDO();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Plastic_BEDO{
    public:
      proxy_nps_Plastic_BEDO(){
        NPS::DetectorFactory::getInstance()->AddToken("Plastic_BEDO","Plastic_BEDO");
        NPS::DetectorFactory::getInstance()->AddDetector("Plastic_BEDO",Plastic_BEDO::Construct);
      }
  };

  proxy_nps_Plastic_BEDO p_nps_Plastic_BEDO;
}
