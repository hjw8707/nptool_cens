/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : December 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  CeBr3 simulation                             *
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
#include "G4Polycone.hh"

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
#include "CeBr3.hh"
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
namespace CeBr3_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.02*MeV;
  const double ResoTime = 19*ns ;
  const double ResoEnergy = 3.8*keV ;
  const double RadiusInternal[6] = {0.*mm, 0.*mm, 27.5*mm, 27.5*mm, 0.*mm, 0.*mm} ; 
  const double RadiusExternal[6] = {28.5*mm, 28.5*mm, 28.5*mm, 28.5*mm, 29.5*mm, 29.5*mm} ; 
  const double Lengths[6] = {0.*mm, 1*mm, 1*mm, 42.5*mm, 42.5*mm, 185*mm} ; 
  const double CrystalRadius = 25.5*mm ; 
  const double CrystalLength = 51*mm ; 
  const string Material = "CeBr3";
  const string WindowMaterial = "Al";
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// CeBr3 Specific Method
CeBr3::CeBr3(){
  m_Event = new TCeBr3Data() ;
  m_CeBr3Scorer = 0;
  m_CylindricalDetector = 0;


  // RGB Color + Transparency
  m_VisCylinder = new G4VisAttributes(G4Colour(1, 0.5, 0, 0.5));   

}

CeBr3::~CeBr3(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CeBr3::AddDetector(G4ThreeVector POS, string  Shape){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void CeBr3::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* CeBr3::BuildDetector(G4int DetNumber, G4ThreeVector Det_pos, G4RotationMatrix* Det_rot, G4LogicalVolume* world){
  if(!m_CylindricalDetector){
    G4VisAttributes* light_GreyAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.7));
    G4VisAttributes* RedAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.6));

    // Global volume
    G4Tubs* globalCeBr = new G4Tubs("CeBr3_Cyl", 0, 60, CeBr3_NS::Lengths[3]*0.5 , 0, 360*deg);
    G4LogicalVolume* global_logical = new G4LogicalVolume(globalCeBr, MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum"), "global_logical", 0, 0, 0);
    new G4PVPlacement(G4Transform3D(*Det_rot, Det_pos), global_logical, "CeBr3", world, false, DetNumber);
    global_logical->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Enveloppe
    G4Material* WMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(CeBr3_NS::WindowMaterial);
    G4Polycone* CeBrtub = new G4Polycone("CeBrtub", 0, 360*deg, 6, CeBr3_NS::Lengths, CeBr3_NS::RadiusInternal, CeBr3_NS::RadiusExternal);
    G4LogicalVolume* tub_log = new G4LogicalVolume(CeBrtub, WMaterial, "tub_log", 0, 0, 0);
    G4ThreeVector tub_pos = G4ThreeVector(0, 0, 0);
    new G4PVPlacement(0, tub_pos, tub_log, "CeBrtub", global_logical, false, DetNumber);
    tub_log->SetVisAttributes(light_GreyAtt);

    // Crystal
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(CeBr3_NS::Material);
    G4Tubs* CeBr_crys = new G4Tubs("CeBr_crys", 0, CeBr3_NS::CrystalRadius, CeBr3_NS::CrystalLength*0.5*mm, 0, 360*deg);
    m_CylindricalDetector = new G4LogicalVolume(CeBr_crys, DetectorMaterial, "logic_CeBr3_crys", 0, 0, 0);
    G4ThreeVector crys_Pos = G4ThreeVector(0, 0, CeBr3_NS::CrystalLength*0.5);
    new G4PVPlacement(0, crys_Pos, m_CylindricalDetector, "CeBr_crys", global_logical, false, DetNumber);
    m_CylindricalDetector->SetVisAttributes(RedAtt);

    m_CylindricalDetector->SetSensitiveDetector(m_CeBr3Scorer);
  }

  return m_CylindricalDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void CeBr3::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("CeBr3");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CeBr3 " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CeBr3 " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,Shape);
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
void CeBr3::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    // Det_pos+=Det_pos.unit()*CeBr3_NS::Thickness*0.5;
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
    
    BuildDetector(i+1, Det_pos, Rot, world);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void CeBr3::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("CeBr3")){
    pTree->Branch("CeBr3", "TCeBr3Data", &m_Event) ;
  }
  pTree->SetBranchAddress("CeBr3", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void CeBr3::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_CeBr3Scorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),CeBr3_NS::ResoEnergy);
    if(Energy>CeBr3_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),CeBr3_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void CeBr3::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_CeBr3Scorer = CheckScorer("CeBr3Scorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_CeBr3Scorer->RegisterPrimitive(Calorimeter);
  m_CeBr3Scorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_CeBr3Scorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* CeBr3::Construct(){
  return  (NPS::VDetector*) new CeBr3();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_CeBr3{
    public:
      proxy_nps_CeBr3(){
        NPS::DetectorFactory::getInstance()->AddToken("CeBr3","CeBr3");
        NPS::DetectorFactory::getInstance()->AddDetector("CeBr3",CeBr3::Construct);
      }
  };

  proxy_nps_CeBr3 p_nps_CeBr3;
}
