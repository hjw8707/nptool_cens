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
 *  This class describe  Coaxial_Germanium simulation                             *
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
#include "Coaxial_Germanium.hh"
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
namespace Coaxial_Germanium_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.01*MeV;
  const double ResoTime = 0*ns ;
  const double ResoEnergy = 1e-9*keV ; 
  const string Material_Shell = "Al";
  const double radii_internal[7] =  {0       , 65/2.*mm-1*mm, 90/2*mm-1*mm, 90/2*mm-1*mm, 90/2*mm-1*mm, 222/2.*mm-1*mm, 0.            };
  const double radii_external[7] =  {65/2.*mm, 65/2.*mm     , 90/2.*mm    , 90/2.*mm    , 222/2.*mm     , 222/2.*mm     , 222/2.*mm     };
  const double length_external[7] = {0       , 0            , 100*mm      , 450*mm      , 450*mm        , 700*mm        , 700*mm        };
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Coaxial_Germanium Specific Method
Coaxial_Germanium::Coaxial_Germanium(){
  m_Event = new TCoaxial_GermaniumData() ;
  m_Coaxial_GermaniumScorer = 0;
}

Coaxial_Germanium::~Coaxial_Germanium(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Coaxial_Germanium::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Coaxial_Germanium::AddDetector(double  R, double  Theta, double  Phi){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
}


G4LogicalVolume* Coaxial_Germanium::BuildDetector(G4int DetNumber, G4ThreeVector Det_pos, G4RotationMatrix* Det_rot, G4LogicalVolume* world){
  // Materials
  G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(Coaxial_Germanium_NS::Material_Shell);
  G4Material* m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* m_MaterialCarbon = MaterialManager::getInstance()->GetMaterialFromLibrary("C");
  G4Material* m_MaterialGermanium = MaterialManager::getInstance()->GetMaterialFromLibrary("Germanium");

  G4VisAttributes* light_GreyAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.7));
  G4VisAttributes* RedAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.6));
  G4VisAttributes* GreenAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.6));

  // Global volume
  // Origin axis at the front of detector
  G4Tubs* solidCoaxial_Germanium = new G4Tubs("solidCoaxial_Germanium", 0, 25, 720*0.5*mm, 0, 360*deg);
  G4LogicalVolume* logicCoaxial_Germanium = new G4LogicalVolume(solidCoaxial_Germanium, m_MaterialVacuum, "logicCoaxial_Germanium", 0, 0);
  new G4PVPlacement(G4Transform3D(*Det_rot, Det_pos), logicCoaxial_Germanium, "Coaxial_Germanium", world, false, DetNumber);
  logicCoaxial_Germanium->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Enveloppe
  G4Polycone* Coaxial_Germanium_Cyl = new G4Polycone("Coaxial_Germanium_Cyl", 0, 360*deg, 7, Coaxial_Germanium_NS::length_external, Coaxial_Germanium_NS::radii_internal, Coaxial_Germanium_NS::radii_external);
  G4LogicalVolume* vol_Coaxial_Germanium = new G4LogicalVolume(Coaxial_Germanium_Cyl, DetectorMaterial, "logic_Coaxial_Germanium_Cyl", 0, 0, 0);
  G4ThreeVector Coaxial_Germanium_cyl_Pos = G4ThreeVector(0, 0, 0);
  new G4PVPlacement(0, Coaxial_Germanium_cyl_Pos, vol_Coaxial_Germanium, "Coaxial_Germanium_cyl", logicCoaxial_Germanium, false, DetNumber);
  vol_Coaxial_Germanium->SetVisAttributes(light_GreyAtt);

  // Germanium crystal
  G4Tubs* Coaxial_Germanium_crys = new G4Tubs("Coaxial_Germanium_crys",0 , 25.*mm, 35.*mm, 0, 360*deg);
  G4LogicalVolume* vol_crys = new G4LogicalVolume(Coaxial_Germanium_crys, m_MaterialGermanium, "logic_Coaxial_Germanium_crys", 0, 0, 0);
  G4ThreeVector crys_Pos = G4ThreeVector(0, 0, 50.*mm);
  new G4PVPlacement(0, crys_Pos, vol_crys, "Coaxial_Germanium_crys", logicCoaxial_Germanium, false, DetNumber);
  vol_crys->SetVisAttributes(RedAtt);

  vol_crys->SetSensitiveDetector(m_Coaxial_GermaniumScorer);
  
  return vol_Coaxial_Germanium;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Coaxial_Germanium::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Coaxial_Germanium");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Coaxial_Germanium " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Coaxial_Germanium " << i+1 <<  endl;
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
void Coaxial_Germanium::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    // Det_pos+=Det_pos.unit()*Coaxial_Germanium_NS::Thickness*0.5;
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
void Coaxial_Germanium::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Coaxial_Germanium")){
    pTree->Branch("Coaxial_Germanium", "TCoaxial_GermaniumData", &m_Event) ;
  }
  pTree->SetBranchAddress("Coaxial_Germanium", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Coaxial_Germanium::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_Coaxial_GermaniumScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Coaxial_Germanium_NS::ResoEnergy);
    if(Energy>Coaxial_Germanium_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Coaxial_Germanium_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Coaxial_Germanium::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_Coaxial_GermaniumScorer = CheckScorer("Coaxial_GermaniumScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_Coaxial_GermaniumScorer->RegisterPrimitive(Calorimeter);
  m_Coaxial_GermaniumScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_Coaxial_GermaniumScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Coaxial_Germanium::Construct(){
  return  (NPS::VDetector*) new Coaxial_Germanium();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Coaxial_Germanium{
    public:
      proxy_nps_Coaxial_Germanium(){
        NPS::DetectorFactory::getInstance()->AddToken("Coaxial_Germanium","Coaxial_Germanium");
        NPS::DetectorFactory::getInstance()->AddDetector("Coaxial_Germanium",Coaxial_Germanium::Construct);
      }
  };

  proxy_nps_Coaxial_Germanium p_nps_Coaxial_Germanium;
}
