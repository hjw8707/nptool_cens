/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Monster simulation                             *
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
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// NPTool header
#include "Monster.hh"
#include "ProcessScorers.hh"
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

namespace Monster_NS
{
  // Energy and time Resolution
  const double EnergyThreshold = 0.*MeV;
  const double Thickness = 5.*cm ;
  const double Radius = 10.*cm ;
}

// Monster Specific Method
Monster::Monster()
{
  m_Event = new TMonsterData() ;
  m_MonsterScorer = 0;
  m_MonsterDetector = 0;

  // RGB Color + Transparency  
  m_VisEJ309   = new G4VisAttributes(G4Colour(0.2, 0.85, 0.85, 1));      

  // Material definition
  m_EJ309   = MaterialManager::getInstance()->GetMaterialFromLibrary("EJ309");
}

Monster::~Monster()
{
}

void Monster::AddDetector(G4ThreeVector POS, string  Shape){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Monster::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}

G4LogicalVolume* Monster::BuildMonsterDetector()
{
  if(!m_MonsterDetector)
  {
    G4Tubs* monster_det = new G4Tubs("Monster_scin", 0, Monster_NS::Radius, Monster_NS::Thickness*0.5, 0., 360*deg);
    m_MonsterDetector = new G4LogicalVolume(monster_det,m_EJ309,"logic_Monster_scin",0,0,0);
    m_MonsterDetector->SetVisAttributes(m_VisEJ309);
    m_MonsterDetector->SetSensitiveDetector(m_MonsterScorer);
  }
  return m_MonsterDetector;
}

// Virtual Method of NPS::VDetector class
// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetectorConstruction::ReadDetectorConfiguration Method
void Monster::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Monster");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Monster " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Monster " << i+1 <<  endl;
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

// Construct detector and inialise sensitive part.
// Called After DetectorConstruction::AddDetector Method
void Monster::ConstructDetector(G4LogicalVolume* world)
{
  for (unsigned short i = 0 ; i < m_R.size() ; i++)
  {
    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*Monster_NS::Thickness*0.5;
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

    if(m_Shape[i] == "Cylindrical")
    {
      new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildMonsterDetector(),
          "Monster",world,false,i+1);
    }
  }
}

// Add Detector branch to the EventTree.
// Called After DetectorConstruction::AddDetector Method
void Monster::InitializeRootOutput()
{
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Monster"))
  {
    pTree->Branch("Monster", "TMonsterData", &m_Event) ;
  }
  pTree->SetBranchAddress("Monster", &m_Event) ;
}

// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAction
void Monster::ReadSensitive(const G4Event* )
{
  m_Event->Clear();

  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_MonsterScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++)
  {
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy =Scorer->GetEnergy(i);
    if( Energy > Monster_NS::EnergyThreshold)
    {
      double Time = Scorer->GetTime(i);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

void Monster::InitializeScorers()
{ 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_MonsterScorer = CheckScorer("MonsterScorer",already_exist) ;

  if(already_exist) return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_MonsterScorer->RegisterPrimitive(Calorimeter);
  m_MonsterScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_MonsterScorer) ;
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Monster::Construct()
{
  return  (NPS::VDetector*) new Monster();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" 
{
  class proxy_nps_Monster
  {
    public:
      proxy_nps_Monster()
      {
        NPS::DetectorFactory::getInstance()->AddToken("Monster","Monster");
        NPS::DetectorFactory::getInstance()->AddDetector("Monster",Monster::Construct);
      }
  };

  proxy_nps_Monster p_nps_Monster;
}
