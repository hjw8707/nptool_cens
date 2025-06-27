/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  quadranMSQ25 simulation                             *
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
#include "quadranMSQ25.hh"
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
namespace quadranMSQ25_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0.*MeV;
  const double ResoTime = 0.*ns ;
  const double ResoEnergy = 0.*MeV ;
  const G4double DetectorSize      = 71.12*mm           ;
  const G4double FrameThickness    = 3.2*mm          ; 
  const G4double SiliconSize       = 50.8*mm           ; //unchecked
  const G4int  NumberOfStripH       = 2              ;
  const G4int  NumberOfStripL       = 2              ;
}
using namespace quadranMSQ25_NS;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// quadranMSQ25 Specific Method
quadranMSQ25::quadranMSQ25(){
  m_Event = new TquadranMSQ25Data() ;
  m_quadranMSQ25Scorer = 0;


}

quadranMSQ25::~quadranMSQ25(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void quadranMSQ25::AddDetector(G4ThreeVector POS, double Thickness){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Thickness.push_back(Thickness);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void quadranMSQ25::AddDetector(double  R, double  Theta, double  Phi, double Thickness){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Thickness.push_back(Thickness);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4LogicalVolume* quadranMSQ25::BuildDetector(G4int DetNumber, G4ThreeVector Det_pos, G4RotationMatrix* Det_rot, G4LogicalVolume* world, double thick){
  // Materials
  G4Material* m_MaterialVacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  G4Material* m_MaterialSilicon = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* m_MaterialPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");

  G4VisAttributes* VisAtt1 = new G4VisAttributes(G4Colour(0.2, 0.5, 0.2));

  G4double NbrTelescopes = DetNumber;
  G4String DetectorNumber;
  std::ostringstream Number;
  Number << NbrTelescopes;
  DetectorNumber = Number.str();
  G4String Name = "quadranMSQ25" + DetectorNumber;
  ////////////////////////////////////////////////////////////////
  /////////General Geometry Parameter Definition /////////////////
  ////////////////////////////////////////////////////////////////
  /////// Starting Volume Definition ///////

  G4Box* solidquadranMSQ25 = new G4Box(Name + "Solid", 0.5 * DetectorSize, 0.5 * DetectorSize, 0.5 * FrameThickness);
  G4LogicalVolume* logicquadranMSQ25 = new G4LogicalVolume(solidquadranMSQ25, m_MaterialVacuum, Name + "_logic", 0, 0);
  new G4PVPlacement(G4Transform3D(*Det_rot, Det_pos), logicquadranMSQ25, Name, world, false, DetNumber);

  logicquadranMSQ25->SetVisAttributes(G4VisAttributes::GetInvisible());
  // Frame is made of 4 thick box (2 Horizontal and 2 Vertical)
  G4Box* solidFrameHorizontal =
      new G4Box(Name + "_Frame", 0.5 * SiliconSize, 0.5 * (DetectorSize - SiliconSize) / 2, 0.5 * FrameThickness * mm);
  G4Box* solidFrameVertical =
      new G4Box(Name + "_Frame", 0.5 * (DetectorSize - SiliconSize) / 2, 0.5 * DetectorSize, 0.5 * FrameThickness * mm);

  G4LogicalVolume* logicFrameHorizontal = new G4LogicalVolume(solidFrameHorizontal, m_MaterialPCB, Name, 0, 0);
  logicFrameHorizontal->SetVisAttributes(VisAtt1);

  G4LogicalVolume* logicFrameVertical = new G4LogicalVolume(solidFrameVertical, m_MaterialPCB, Name, 0, 0);
  logicFrameVertical->SetVisAttributes(VisAtt1);

  G4ThreeVector FrameTopPosition = G4ThreeVector(0, 0.5 * SiliconSize + 0.5 * (DetectorSize - SiliconSize) / 2, 0);
  G4ThreeVector FrameBottomPosition = G4ThreeVector(0, -0.5 * SiliconSize - 0.5 * (DetectorSize - SiliconSize) / 2, 0);
  G4ThreeVector FrameLeftPosition = G4ThreeVector(0.5 * SiliconSize + 0.5 * (DetectorSize - SiliconSize) / 2, 0, 0);
  G4ThreeVector FrameRightPosition = G4ThreeVector(-0.5 * SiliconSize - 0.5 * (DetectorSize - SiliconSize) / 2, 0, 0);

  new G4PVPlacement(0, FrameTopPosition, logicFrameHorizontal, Name + "_FrameT", logicquadranMSQ25, false, DetNumber);
  new G4PVPlacement(0, FrameBottomPosition, logicFrameHorizontal, Name + "_FrameH", logicquadranMSQ25, false, DetNumber);
  new G4PVPlacement(0, FrameLeftPosition, logicFrameVertical, Name + "_FrameL", logicquadranMSQ25, false, DetNumber);
  new G4PVPlacement(0, FrameRightPosition, logicFrameVertical, Name + "_FrameR", logicquadranMSQ25, false, DetNumber);

  G4ThreeVector posSi = G4ThreeVector(0, 0, 0);
  G4Box* solidSi = new G4Box("quadranMSQ25", 0.5 * SiliconSize, 0.5 * SiliconSize, 0.5 * thick);
  G4LogicalVolume* logicSi = new G4LogicalVolume(solidSi, m_MaterialSilicon, "logicSi", 0, 0, 0);

  new G4PVPlacement(0, posSi, logicSi, Name + "_Si", logicquadranMSQ25, true, DetNumber);

  logicSi->SetSensitiveDetector(m_quadranMSQ25Scorer);

  return logicSi;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void quadranMSQ25::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("quadranMSQ25");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Thickness"};
  vector<string> sphe = {"R","Theta","Phi","Thickness"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  quadranMSQ25 " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      double Thickness = blocks[i]->GetDouble("Thickness","um");
      AddDetector(Pos,Thickness);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  quadranMSQ25 " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double Thickness = blocks[i]->GetDouble("Thickness","um");
      AddDetector(R,Theta,Phi,Thickness);
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
void quadranMSQ25::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*quadranMSQ25_NS::FrameThickness*0.5;
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
    
    BuildDetector(i+1, Det_pos, Rot, world, m_Thickness[i]);
   
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void quadranMSQ25::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("quadranMSQ25")){
    pTree->Branch("quadranMSQ25", "TquadranMSQ25Data", &m_Event) ;
  }
  pTree->SetBranchAddress("quadranMSQ25", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void quadranMSQ25::ReadSensitive(const G4Event* ){
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_quadranMSQ25Scorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult(); 
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),quadranMSQ25_NS::ResoEnergy);
    if(Energy>quadranMSQ25_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),quadranMSQ25_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void quadranMSQ25::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_quadranMSQ25Scorer = CheckScorer("quadranMSQ25Scorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_quadranMSQ25Scorer->RegisterPrimitive(Calorimeter);
  m_quadranMSQ25Scorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_quadranMSQ25Scorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* quadranMSQ25::Construct(){
  return  (NPS::VDetector*) new quadranMSQ25();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_quadranMSQ25{
    public:
      proxy_nps_quadranMSQ25(){
        NPS::DetectorFactory::getInstance()->AddToken("quadranMSQ25","quadranMSQ25");
        NPS::DetectorFactory::getInstance()->AddDetector("quadranMSQ25",quadranMSQ25::Construct);
      }
  };

  proxy_nps_quadranMSQ25 p_nps_quadranMSQ25;
}
