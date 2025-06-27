/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Tina simulation                                     * 
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>

// Geant4 
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Randomize.hh"

// NPTool 
#include "STARK.hh"
#include "STARKScorers.hh"
#include "DSSDScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
#include "NPCore.h"

// CLHEP 
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;
using namespace STARKNS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
STARK::STARK(){
  HCID_X6 = HCID_BB10 = HCID_QQQ5 = -1;
  
  m_X6 = m_BB10 = m_QQQ5 = NULL;
  m_X6Det = m_BB10Det = m_QQQ5Det = NULL;

  m_VisX6       = new G4VisAttributes(G4Colour(0., 0.5, 0.5));
  m_VisX6PCB    = new G4VisAttributes(G4Colour(0.8, 0.5, 0.5));
  m_VisBB10     = new G4VisAttributes(G4Colour(0., 0.5, 0.7));
  m_VisBB10PCB  = new G4VisAttributes(G4Colour(0.8, 0.5, 0.7));
  m_VisQQQ5     = new G4VisAttributes(G4Colour(0., 0.5, 0.3));
  m_VisQQQ5PCB  = new G4VisAttributes(G4Colour(0.8, 0.5, 0.3));
  m_VisConn     = new G4VisAttributes(G4Colour(0.8, 0.8, 0.8));

  m_Event = new TSTARKData();
  m_Raw = new TSTARKRaw(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
STARK::~STARK(){}

void STARK::AddDetector(string Type, G4ThreeVector Pos, int Flip, int Rev, double Beta) {
  m_Type.push_back(Type);
  m_Pos.push_back(Pos);
  m_Flip.push_back(Flip);
  m_Rev.push_back(Rev);
  m_Beta.push_back(Beta);}

//////////////////////////////////////////////////////////////////////
//
// BuildX6Detector: making a G4AssemblyVolume for X6 (only once)
//
G4AssemblyVolume* STARK::BuildX6Detector(){
  if (m_X6) return m_X6;

  ////////////////////////////////////////////////////////////
  // material definition
  G4Material* matSi  = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  ////////////////////////////////////////////////////////////

  G4Box* solidX6PCBAll  = new G4Box("solidX6PCBAll",  X6_PCBX/2.,X6_PCBY/2.,X6_PCBZ/2.);
  G4Box* solidX6PCBSub1 = new G4Box("solidX6PCBSub1", X6_PCBSub1X/2.,X6_PCBSub1Y/2.,X6_PCBSub1Z/2.+0.01*mm); // +0.01 mm for the perfect subtraction of solid
  G4Box* solidX6PCBSub2 = new G4Box("solidX6PCBSub2", X6_PCBSub2X/2.,X6_PCBSub2Y/2.,X6_PCBSub2Z/2.+0.01*mm); // +0.01 mm for the perfect subtraction of solid
  
  G4VSolid* solidX6PCBTemp = new G4SubtractionSolid("solidX6PCBTemp",
						    solidX6PCBAll, solidX6PCBSub1, 0,
						    G4ThreeVector(X6_PCBSub1XOffset,
								  X6_PCBSub1YOffset,
								  X6_PCBSub1ZOffset));
  G4VSolid* solidX6PCB = new G4SubtractionSolid("solidX6PCB",
						solidX6PCBTemp, solidX6PCBSub2, 0,
						G4ThreeVector(X6_PCBSub2XOffset,
							      X6_PCBSub2YOffset,
							      X6_PCBSub2ZOffset));
  
 
  G4LogicalVolume* logicX6PCB = new G4LogicalVolume(solidX6PCB, matPCB,"logicX6PCB",0,0,0);
  logicX6PCB->SetVisAttributes(m_VisX6PCB);

  G4Box* solidX6Conn = new G4Box("solidX6Conn", Conn_X/2.,Conn_Y/2.,Conn_Z/2.);
  G4LogicalVolume* logicX6Conn = new G4LogicalVolume(solidX6Conn, matPCB,"logicX6Conn",0,0,0);
  logicX6Conn->SetVisAttributes(m_VisConn);
  
  G4Box* solidX6Si = new G4Box("solidX6Si",
			       X6_SiX/2.,X6_SiY/2.,X6_SiZ/2.);
  G4LogicalVolume* logicX6Si = new G4LogicalVolume(solidX6Si, matSi,"logicX6Si",0,0,0);
  logicX6Si->SetVisAttributes(m_VisX6);
  logicX6Si->SetSensitiveDetector(m_X6Det); 

  m_X6 = new G4AssemblyVolume();
  G4ThreeVector Pos;
  m_X6->AddPlacedVolume(logicX6Si,Pos,0); // reference = center of the X6 Si wafer
  Pos = G4ThreeVector(-X6_SiXOffset,-X6_SiYOffset,-X6_SiZOffset);
  m_X6->AddPlacedVolume(logicX6PCB,Pos,0);
  Pos = G4ThreeVector(-X6_SiXOffset,
		      -X6_SiYOffset - X6_PCBY/2. + Conn_Y/2.,
		      X6_PCBZ/2. + Conn_Z/2.);
  m_X6->AddPlacedVolume(logicX6Conn,Pos,0);
			
  return m_X6; }

//////////////////////////////////////////////////////////////////////
//
// BuildBB10Detector: making a G4AssemblyVolume for BB10 (only once)
//
G4AssemblyVolume* STARK::BuildBB10Detector(){
  if (m_BB10) return m_BB10;

  ////////////////////////////////////////////////////////////
  // material definition
  G4Material* matSi = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  ////////////////////////////////////////////////////////////

    G4Box* solidBB10PCBAll  = new G4Box("solidBB10PCBAll",  BB10_PCBX/2.,BB10_PCBY/2.,BB10_PCBZ/2.);
  G4Box* solidBB10PCBSub1 = new G4Box("solidBB10PCBSub1", BB10_PCBSub1X/2.,BB10_PCBSub1Y/2.,BB10_PCBSub1Z/2.+0.01*mm); // +0.01 mm for the perfect subtraction of solid
  G4Box* solidBB10PCBSub2 = new G4Box("solidBB10PCBSub2", BB10_PCBSub2X/2.,BB10_PCBSub2Y/2.,BB10_PCBSub2Z/2.+0.01*mm); // +0.01 mm for the perfect subtraction of solid
  
  G4VSolid* solidBB10PCBTemp = new G4SubtractionSolid("solidBB10PCBTemp",
						    solidBB10PCBAll, solidBB10PCBSub1, 0,
						    G4ThreeVector(BB10_PCBSub1XOffset,
								  BB10_PCBSub1YOffset,
								  BB10_PCBSub1ZOffset));
  G4VSolid* solidBB10PCB = new G4SubtractionSolid("solidBB10PCB",
						solidBB10PCBTemp, solidBB10PCBSub2, 0,
						G4ThreeVector(BB10_PCBSub2XOffset,
							      BB10_PCBSub2YOffset,
							      BB10_PCBSub2ZOffset));
  
 
  G4LogicalVolume* logicBB10PCB = new G4LogicalVolume(solidBB10PCB, matPCB,"logicBB10PCB",0,0,0);
  logicBB10PCB->SetVisAttributes(m_VisBB10PCB);

  G4Box* solidBB10Conn = new G4Box("solidBB10Conn", Conn_X/2.,Conn_Y/2.,Conn_Z/2.);
  G4LogicalVolume* logicBB10Conn = new G4LogicalVolume(solidBB10Conn, matPCB,"logicBB10Conn",0,0,0);
  logicBB10Conn->SetVisAttributes(m_VisConn);
  
  G4Box* solidBB10Si = new G4Box("solidBB10Si",
			       BB10_SiX/2.,BB10_SiY/2.,BB10_SiZ/2.);
  G4LogicalVolume* logicBB10Si = new G4LogicalVolume(solidBB10Si, matSi,"logicBB10Si",0,0,0);
  logicBB10Si->SetVisAttributes(m_VisBB10);
  logicBB10Si->SetSensitiveDetector(m_BB10Det); 

  m_BB10 = new G4AssemblyVolume();
  G4ThreeVector Pos;
  m_BB10->AddPlacedVolume(logicBB10Si,Pos,0); // reference = center of the BB10 Si wafer
  Pos = G4ThreeVector(-BB10_SiXOffset,-BB10_SiYOffset,-BB10_SiZOffset);
  m_BB10->AddPlacedVolume(logicBB10PCB,Pos,0);
  Pos = G4ThreeVector(-BB10_SiXOffset,
		      -BB10_SiYOffset - BB10_PCBY/2. + Conn_Y/2.,
		      BB10_PCBZ/2. + Conn_Z/2.);
  m_BB10->AddPlacedVolume(logicBB10Conn,Pos,0);
			
  return m_BB10; }

//////////////////////////////////////////////////////////////////////
//
// BuildQQQ5Detector: making a G4AssemblyVolume for QQQ5 (only once)
//
G4AssemblyVolume* STARK::BuildQQQ5Detector(){
  if (m_QQQ5) return m_QQQ5;

  ////////////////////////////////////////////////////////////
  // material definition
  G4Material* matSi = MaterialManager::getInstance()->GetMaterialFromLibrary("Si");
  G4Material* matPCB = MaterialManager::getInstance()->GetMaterialFromLibrary("PCB");
  ////////////////////////////////////////////////////////////

  // Make the a single detector geometry
  G4Tubs* solidQQQ5PCBAll = new G4Tubs("solidQQQ5PCBAll", QQQ5_PCBInR, QQQ5_PCBOutR,
				       QQQ5_PCBT*0.5, QQQ5_PCBPhi0, QQQ5_PCBPhi1);
  G4Box* solidQQQ5PCBCut = new G4Box("solidQQQ5PCBCut",
				     QQQ5_PCBCutX/2., QQQ5_PCBCutY/2., QQQ5_PCBCutZ/2.);
  G4VSolid* solidQQQ5PCBSub1 =  new G4SubtractionSolid("solidQQQ5PCBSub1",
						       solidQQQ5PCBAll, solidQQQ5PCBCut,0,
						       G4ThreeVector(QQQ5_PCBCutXOffset,
								     QQQ5_PCBCutYOffset,
								     QQQ5_PCBCutZOffset));
  
  G4Tubs* solidQQQ5SiAllForPCBCut = new G4Tubs("solidQQQ5SiAllForPCBCut",
					       QQQ5_SiInR, QQQ5_SiOutR,
					       QQQ5_PCBCutZ/2., QQQ5_SiPhi0, QQQ5_SiPhi1);
  G4Box* solidQQQ5SiCut1ForPCBCut = new G4Box("solidQQQ5SiCut1ForPCBCut",
					      QQQ5_SiCut1X/2., QQQ5_SiCut1Y/2., QQQ5_PCBCutZ/2.+1*mm);
  G4Box* solidQQQ5SiCut2ForPCBCut = new G4Box("solidQQQ5SiCut2ForPCBCut",
					      QQQ5_SiCut2X/2., QQQ5_SiCut2Y/2., QQQ5_PCBCutZ/2.+1*mm);
  G4VSolid* solidQQQ5SiSub1ForPCBCut =  new G4SubtractionSolid("solidQQQ5SiSub1ForPCBCut",
							       solidQQQ5SiAllForPCBCut,
							       solidQQQ5SiCut1ForPCBCut,0,
							       G4ThreeVector(QQQ5_SiCut1XOffset,
									     QQQ5_SiCut1YOffset,
									     QQQ5_SiCut1ZOffset));
  G4VSolid* solidQQQ5SiForPCBCut = new G4SubtractionSolid("solidQQQ5SiForPCBCut",
							  solidQQQ5SiSub1ForPCBCut,
							  solidQQQ5SiCut2ForPCBCut,0,
							  G4ThreeVector(QQQ5_SiCut2XOffset,
									QQQ5_SiCut2YOffset,
									QQQ5_SiCut2ZOffset));  
  G4VSolid* solidQQQ5PCB =  new G4SubtractionSolid("solidQQQ5PCB",
						   solidQQQ5PCBSub1, solidQQQ5SiForPCBCut,0,
						   G4ThreeVector(0, 0, 0));
  
  
  
  
  G4Tubs* solidQQQ5SiAll = new G4Tubs("solidQQQ5SiAll", QQQ5_SiInR, QQQ5_SiOutR,
				      QQQ5_SiT*0.5, QQQ5_SiPhi0, QQQ5_SiPhi1);
  G4Box* solidQQQ5SiCut1 = new G4Box("solidQQQ5SiCut1",
				     QQQ5_SiCut1X/2., QQQ5_SiCut1Y/2., QQQ5_SiCut1Z/2.);
  G4Box* solidQQQ5SiCut2 = new G4Box("solidQQQ5SiCut2",
				     QQQ5_SiCut2X/2., QQQ5_SiCut2Y/2., QQQ5_SiCut2Z/2.);
  G4VSolid* solidQQQ5SiSub1 =  new G4SubtractionSolid("solidQQQ5SiSub1",
						      solidQQQ5SiAll, solidQQQ5SiCut1,0,
						      G4ThreeVector(QQQ5_SiCut1XOffset,
								    QQQ5_SiCut1YOffset,
								    QQQ5_SiCut1ZOffset));
  G4VSolid* solidQQQ5Si = new G4SubtractionSolid("solidQQQ5Si",
						 solidQQQ5SiSub1, solidQQQ5SiCut2,0,
						 G4ThreeVector(QQQ5_SiCut2XOffset,
							       QQQ5_SiCut2YOffset,
							       QQQ5_SiCut2ZOffset));  
  

  G4LogicalVolume* logicQQQ5PCB = new G4LogicalVolume(solidQQQ5PCB, matPCB,"logicQQQ5PCB",0,0,0);
  logicQQQ5PCB->SetVisAttributes(m_VisQQQ5PCB);

  G4Box* solidQQQ5Conn = new G4Box("solidQQQ5Conn", Conn_X/2.,Conn_Y/2.,Conn_Z/2.);
  G4LogicalVolume* logicQQQ5Conn = new G4LogicalVolume(solidQQQ5Conn, matPCB,"logicQQQ5Conn",0,0,0);
  logicQQQ5Conn->SetVisAttributes(m_VisConn);

  G4LogicalVolume* logicQQQ5Si = new G4LogicalVolume(solidQQQ5Si, matSi,"logicQQQ5",0,0,0);
  logicQQQ5Si->SetVisAttributes(m_VisQQQ5);
  logicQQQ5Si->SetSensitiveDetector(m_QQQ5Det); 

  m_QQQ5 = new G4AssemblyVolume();
  G4ThreeVector Pos;
  m_QQQ5->AddPlacedVolume(logicQQQ5Si,Pos,0);
  m_QQQ5->AddPlacedVolume(logicQQQ5PCB,Pos,0);
  Pos = G4ThreeVector(Conn_Y/2.,
		      (QQQ5_SiOutR+QQQ5_SiInR)/2.,
		      QQQ5_PCBT/2.+Conn_Z/2.);
  G4RotationMatrix *rot = new G4RotationMatrix;
  rot->rotateZ(90*deg);
  m_QQQ5->AddPlacedVolume(logicQQQ5Conn,Pos,rot);

  return m_QQQ5; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ReadConfiguration: reading a geometry file for each Si detector
void STARK::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("STARK");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> reso = {"Type", "Reso"}; // for put resolutions as an input parameter
  vector<string> cart = {"Type", "POS"};
  vector<string> sphe = {"Type", "R","Theta","Phi"};
  vector<string> cyld = {"Type", "Rho","Phi","Z"};

  for (unsigned int i = 0 ; i < blocks.size() ; i++){
    ////////////////////////////////////////////////////////////
    // Resolution
    if (blocks[i]->HasTokenList(reso)) {
      string Type = blocks[i]->GetString("Type");
      if      (Type.compare("X6")   == 0) X6_ERes   = blocks[i]->GetDouble("Reso","void");
      else if (Type.compare("BB10") == 0) BB10_ERes = blocks[i]->GetDouble("Reso","void");
      else if (Type.compare("QQQ5") == 0) QQQ5_ERes = blocks[i]->GetDouble("Reso","void");
      continue; // only block for resolution
    }
    ////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////
    // Cartesian coordinate
    if (blocks[i]->HasTokenList(cart)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARK " << i+1 <<  endl;
      string Type = blocks[i]->GetString("Type");
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      int Flip = 0;
      if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
      int Rev = 0;
      if (blocks[i]->HasToken("Rev")) Rev = blocks[i]->GetInt("Rev");      
      double Beta = 0;
      if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");
      cout << endl << "////  Type " << Type << endl;
      cout << "////  Pos " << Pos <<  endl;
      cout << "////  Flip " << Flip <<  endl;
      cout << "////  Rev " << Rev <<  endl;
      cout << "////  Beta " << Beta <<  endl;
      AddDetector(Type, Pos, Flip, Rev, Beta); }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Spherical coordinate
    else if (blocks[i]->HasTokenList(sphe)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARK " << i+1 <<  endl;
      string Type = blocks[i]->GetString("Type");
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      int Flip = 0;
      if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
      int Rev = 0;
      if (blocks[i]->HasToken("Rev")) Rev = blocks[i]->GetInt("Rev");            
      double Beta = 0;
      if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");
      
      G4ThreeVector Pos;
      Pos.setRThetaPhi(R, Theta, Phi);
      cout << endl << "////  Type " << Type << endl;
      cout << "////  Pos " << Pos <<  endl;
      cout << "////  Flip " << Flip <<  endl;
      cout << "////  Rev " << Rev <<  endl;
      cout << "////  Beta " << Beta <<  endl;
      AddDetector(Type, Pos, Flip, Rev, Beta); }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Cylindrical coordinate
    else if (blocks[i]->HasTokenList(cyld)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARK " << i+1 <<  endl;
      string Type = blocks[i]->GetString("Type");
      double Rho = blocks[i]->GetDouble("Rho","mm");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double Z   = blocks[i]->GetDouble("Z","mm");
      int Flip = 0;
      if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
      int Rev = 0;
      if (blocks[i]->HasToken("Rev")) Rev = blocks[i]->GetInt("Rev");      
      double Beta = 0;
      if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");

      G4ThreeVector Pos;
      Pos.setRhoPhiZ(Rho, Phi, Z);
      cout << endl << "////  Type " << Type << endl;
      cout << "////  Pos " << Pos <<  endl;
      cout << "////  Flip " << Flip <<  endl;
      cout << "////  Rev " << Rev <<  endl;      
      cout << "////  Beta " << Beta <<  endl;
      AddDetector(Type, Pos, Flip, Rev, Beta); }
    ////////////////////////////////////////////////////////////
    
    else {
      cout << "Error: check your input file formatting" << endl;
      exit(1); }
  }
  std::cout << "read complete" << std::endl;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::ConstructDetector(G4LogicalVolume* world){
  std::cout << "start constuct detector" << std::endl;
  for (unsigned short i = 0 ; i < m_Pos.size() ; i++) {
    G4double phi = m_Pos[i].getPhi();
    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    if (m_Rev[i]) Rot->rotateZ(180*deg);
    if (m_Flip[i]) Rot->rotateY(180*deg);
    Rot->rotateX(90*deg);
    Rot->rotateZ(90*deg + phi);

    G4AssemblyVolume* det;
    if      (m_Type[i] == "X6"  ) det = BuildX6Detector  ();
    else if (m_Type[i] == "BB10") det = BuildBB10Detector();
    else if (m_Type[i] == "QQQ5") {
      det = BuildQQQ5Detector();
      Rot = new G4RotationMatrix;
      Rot->rotateZ(m_Beta[i]);
      if (m_Flip[i]) Rot->rotateY(180*deg);    }
    else {
      std::cerr << "no type " << m_Type[i] << " exists.\n";
      continue; }
    
    det->MakeImprint(world,m_Pos[i],Rot,i+1,true);

    // iterator is equal to fPVStore.begin()
    std::vector< G4VPhysicalVolume* >::iterator it = det->GetVolumesIterator();
    unsigned int NbrImprints   = det->GetImprintsCount();
    unsigned int NbrTotalPV    = det->TotalImprintedVolumes(); 
    unsigned int NbrComponents = NbrTotalPV/NbrImprints;
    // set copy numbers of components of assembly volume to the current detector number
    for (it += (NbrImprints-1) * NbrComponents ;
	 it <= det->GetVolumesIterator() + NbrTotalPV-1 ;
	 it++)
      (*it)->SetCopyNo(i+1);
  }
  std::cout << "construct complete" << std::endl; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::InitializeRootOutput(){
  // add detector branch to the EventTree
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("STARK")) pTree->Branch("STARK","TSTARKData",&m_Event);
  if (!pTree->FindBranch("Raw"))   pTree->Branch("Raw","TSTARKRaw",&m_Raw);
  pTree->SetBranchAddress("STARK",&m_Event);
  pTree->SetBranchAddress("Raw",&m_Raw);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::ReadSensitive(const G4Event* event){
  // read sensitive part and fill the Root tree
  m_Event->Clear();
  m_Raw->Clear();

  auto HCE = event->GetHCofThisEvent();
  if (!HCE) return;

  if (HCID_X6 == -1)
    HCID_X6 = G4SDManager::GetSDMpointer()->GetCollectionID("X6Det/X6Scorer");
  if (HCID_BB10 == -1)
    HCID_BB10 = G4SDManager::GetSDMpointer()->GetCollectionID("BB10Det/BB10Scorer");  
  if (HCID_QQQ5 == -1)
    HCID_QQQ5 = G4SDManager::GetSDMpointer()->GetCollectionID("QQQ5Det/QQQ5Scorer");  

  
  /////////////////////////////////////////////////////////////////////////////////
  // loop for the event map
  /////////////////////////////////////
  //  for X6
  auto evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_X6));
  map<G4int, G4double**>::iterator it;
  for (it = evtMap->GetMap()->begin() ; it != evtMap->GetMap()->end() ; it++) {
    // energy smearing
    G4double enSmear2  = RandGauss::shoot((*(it->second))[2], (*(it->second))[2]*X6_ERes/100.);
    G4double enSmear0 = RandGauss::shoot((*(it->second))[0], (*(it->second))[0]*X6_ERes/100.);
    G4double enSmear1 = RandGauss::shoot((*(it->second))[1], (*(it->second))[1]*X6_ERes/100.);
    m_Event->Set(0,
		 (*(it->second))[4], // detector number
		 (*(it->second))[5], // front strip number
		 (*(it->second))[6], // back strip number
		 enSmear2, // frontside energy
		 enSmear2, // backside energy
		 enSmear0, // upstream energy
		 enSmear1, // downstream energy
		 (*(it->second))[3], // global time
		 (*(it->second))[7], // hit position X
		 (*(it->second))[8], // hit position Y
		 (*(it->second))[9]); // hit position Z
    m_Raw->SetRaw(m_Event); }
  /////////////////////////////////////
  // for BB10
  evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_BB10));
  for (it = evtMap->GetMap()->begin() ; it != evtMap->GetMap()->end() ; it++) {
    // energy smearing
    G4double enSmear0  = RandGauss::shoot((*(it->second))[0], (*(it->second))[0]*BB10_ERes/100.);
    G4double enSmear1  = RandGauss::shoot((*(it->second))[1], (*(it->second))[1]*BB10_ERes/100.);
    m_Event->Set(1,
		 (*(it->second))[3], // detector number
		 (*(it->second))[4], // front strip number
		 1, // back strip number (only 1)
		 enSmear0, // frontside energy
		 enSmear1, // backside energy
		 0, // upstream energy
		 0, // downstream energy
		 (*(it->second))[2], // global time
		 (*(it->second))[5], // hit position X
		 (*(it->second))[6], // hit position Y
		 (*(it->second))[7]); // hit position Z
  }
  /////////////////////////////////////
  // for QQQ5
  evtMap = static_cast<NPS::HitsMap<G4double*>*>(HCE->GetHC(HCID_QQQ5));
  for (it = evtMap->GetMap()->begin() ; it != evtMap->GetMap()->end() ; it++) {
    // energy smearing
    G4double enSmear0  = RandGauss::shoot((*(it->second))[0], (*(it->second))[0]*BB10_ERes/100.);
    G4double enSmear1  = RandGauss::shoot((*(it->second))[1], (*(it->second))[1]*BB10_ERes/100.);    
    m_Event->Set(2,
		 (*(it->second))[3], // detector number
		 (*(it->second))[4], // front strip number
		 (*(it->second))[5], // back strip number
		 enSmear0, // frontside energy
		 enSmear1, // backside energy
		 0, // upstream energy
		 0, // downstream energy
		 (*(it->second))[2], // global time
		 (*(it->second))[6], // hit position X
		 (*(it->second))[7], // hit position Y
		 (*(it->second))[8]); // hit position Z
  }
  /////////////////////////////////////////////////////////////////////////////////
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void STARK::InitializeScorers() {
  // Check the detectors initialized
  bool already_exist_X6 = false;
  bool already_exist_BB10 = false;
  bool already_exist_QQQ5 = false; 
  m_X6Det   = CheckScorer("X6Det",already_exist_X6); // MultiFunctionalDetector
  m_BB10Det = CheckScorer("BB10Det",already_exist_BB10); // MultiFunctionalDetector
  m_QQQ5Det = CheckScorer("QQQ5Det",already_exist_QQQ5); // MultiFunctionalDetector

  // if not, create them
  if (!already_exist_X6) {
    G4VPrimitiveScorer* X6Scorer = new STARKSCORERS::PS_STARK_X6
      ("X6Scorer", 0, X6_SiActiveX, X6_SiActiveY, X6_NFrontStrips,X6_NBackStrips,0);
    m_X6Det->RegisterPrimitive(X6Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_X6Det); }

  if (!already_exist_BB10) {
    G4VPrimitiveScorer* BB10Scorer = new STARKSCORERS::PS_STARK_BB10
      ("BB10Scorer", 0, BB10_SiActiveX, BB10_SiActiveY, BB10_NFrontStrips, 0);
    m_BB10Det->RegisterPrimitive(BB10Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_BB10Det); }

  if (!already_exist_QQQ5) {
    G4VPrimitiveScorer* QQQ5Scorer =  new STARKSCORERS::PS_STARK_QQQ5
      ("QQQ5Scorer", 0, QQQ5_SiActiveInR, QQQ5_SiActiveOutR, 
       QQQ5_SiPhi0, QQQ5_SiPhi1,
       QQQ5_NAStrip, QQQ5_NRStrip, 0);
    m_QQQ5Det->RegisterPrimitive(QQQ5Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_QQQ5Det); }}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// construct method to be passed to the DetectorFactory
NPS::VDetector* STARK::Construct(){
  return  (NPS::VDetector*) new STARK();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// register the construct method to the factory                
extern"C" {
  class proxy_nps_STARK{
    public:
      proxy_nps_STARK(){
        NPS::DetectorFactory::getInstance()->AddToken("STARK","STARK");
        NPS::DetectorFactory::getInstance()->AddDetector("STARK",STARK::Construct);
      }
  };
  proxy_nps_STARK p_nps_STARK;
}
