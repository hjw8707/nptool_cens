/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon contact address: audrey.chatillon@cea.fr*
 *                                                                           *
 * Creation Date  : décembre 2024                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describes the Epic fission chamber simulation                 *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: EPIC stands for Even Plutonium Isotopes fission Cross section    *
 *          measurement                                                      *
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
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"

// NPTool header
#include "Epic.hh"
#include "GaseousDetectorScorers.hh"
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
namespace Epic_NS{
  
  // =====================================
  // Thickness of the foils for the anodes
  const double Cu_Thickness = 6*um;
  const double Kapton_Thickness = 6*um;
  
  // =================
  // Energy Resolution
  const double ResoEnergyPerCent = 0.1 ;

  // ===============
  // Time Resolution
  const double Gauss_Sigma     = 10.*ns ;
  const int    Gauss_nSigma    = 3;
  const double Gauss_HalfRange = Gauss_Sigma * Gauss_nSigma ; 


}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Epic Specific Method
Epic::Epic(){
  m_Event = new TEpicData() ;
  m_EpicScorer = 0;
  m_EpicVolume = 0;

  // RGB Color + Transparency
  m_VisFCWall = new G4VisAttributes(G4Colour(0.1,0.5,0.7,1));
  m_VisAl = new G4VisAttributes(G4Colour(0.839,0.803,0.803,1));
  m_VisTi = new G4VisAttributes(G4Colour(0.776,0.662,0.662,0.5));
  m_VisGasSensitive = new G4VisAttributes(G4Colour(0.576,0.662,0.662,0.5));
  m_VisGas = new G4VisAttributes(G4Colour(0.576,0.662,0.662,0.1));
  m_VisCu = new G4VisAttributes(G4Colour(0.70, 0.40, 0. ,1));
  m_VisRogers4003C = new G4VisAttributes(G4Colour(0.60, 0.60, 0.2 ,1));
  m_VisSample235U = new G4VisAttributes(G4Colour(1.0 , 0.80, 0.0 ,1));
  m_VisSample238U = new G4VisAttributes(G4Colour(0.76, 0.88, 0.0 ,1));
  m_VisSamplePu   = new G4VisAttributes(G4Colour(0.54, 0.88, 0.0 ,1));
  m_VisSample252Cf = new G4VisAttributes(G4Colour(0.9, 0.00, 0.0 ,1));
}

Epic::~Epic(){
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Epic::AddDetector(G4ThreeVector POS){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetectorConstruction::ReadDetectorConfiguration Method
void Epic::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Epic");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","GasId","Pressure"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Epic " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      AddDetector(Pos);

      m_GasId            = blocks[i]->GetString("GasId");
      m_Pressure         = blocks[i]->GetDouble("Pressure","bar");
      m_IonisationEnergy = blocks[i]->GetDouble("IonisationEnergy","eV");
      m_DriftStepLength  = blocks[i]->GetDouble("DriftStepLength","um");
      m_DriftVelocity    = blocks[i]->GetDouble("DriftVelocity","mm_ns"); // FIXME WARNING: Unknown unit mm_ns 
      m_DriftWindowMax   = blocks[i]->GetDouble("DriftWindowMax","ns");
      m_GasMaterial = MaterialManager::getInstance()->GetGasFromLibrary(m_GasId, m_Pressure, 300*kelvin);
      m_SegmentedAnode   = blocks[i]->GetInt("SegmentedAnode");
      m_nA               = blocks[i]->GetInt("nAnodes");
      m_Distance_AK      = blocks[i]->GetDouble("Distance_AK","mm");
      m_InterDistance_KK = blocks[i]->GetDouble("InterDistance_KK","mm"); 
      m_Thickness_K      = blocks[i]->GetDouble("Thickness_K","um"); 
      m_KId              = blocks[i]->GetString("KId");
      m_KMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(m_KId);
      m_nSamplesPerA     = blocks[i]->GetVectorInt("nSamplesPerAnode");
      m_SampleMaterial   = blocks[i]->GetVectorString("SampleMaterialPerAnode");
      m_SampleThickness  = blocks[i]->GetVectorDouble("SampleThickness","um");
      if(m_SegmentedAnode) m_RadiusAnode = 20.0*mm;
      else                 m_RadiusAnode = 37.0*mm;
      m_DriftStepTime    = m_DriftStepLength / m_DriftVelocity ;
      m_Gauss_nbins      = (2*Epic_NS::Gauss_HalfRange) / m_DriftStepTime;
      m_Convo_nbins      = 2.*m_Gauss_nbins + ceil((m_DriftWindowMax + m_Distance_AK/m_DriftVelocity)/m_DriftStepTime);     
      //cout << "m_DriftStepLength = " << m_DriftStepLength << " mm" << endl;
      //cout << "m_DriftVelocity   = " << m_DriftVelocity << " mm/ns" << endl;
      //cout << "m_DriftStepTime = " << m_DriftStepTime << " ns" << endl;
      //cout << "m_DriftWindowMax = " << m_DriftWindowMax << " ns" << endl;
      //cout << "m_Gauss_nbins   = " << m_Gauss_nbins << endl;
      //cout << "m_Convo_nbins   = " << m_Convo_nbins << endl;
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
void Epic::ConstructDetector(G4LogicalVolume* world){
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    
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

    G4RotationMatrix* Rot = new G4RotationMatrix(0,0,0);
    BuildEpic()->MakeImprint(world, Det_pos, Rot, i);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* Epic::BuildEpic(){
  cout << "Enter in Epic::BuildEpic" << endl;


  m_EpicVolume = new G4AssemblyVolume();

  G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(0);

  int cpt_AddVolumes = 0;

  // --- Get Material
  G4Material* Al_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  G4Material* Cu_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
  G4Material* Ti_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Ti");
  G4Material* Rogers_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Rogers4003C");


  // --- posY = 0mm @center of the samples, 75 mm above the upper side of the PCB
  // --- PCB (6 layers of Cu)
  double PCB_width  = 180.0*mm;
  double PCB_length = 330.0*mm;
  double PCB_Rogers_height = 1.6*mm;      double posY_PCB_Rogers = -75.0*mm - 0.5*PCB_Rogers_height*mm;    
  double PCB_Cu_height = 6*0.035*mm;      double posY_PCB_Cu     = -75.0*mm - PCB_Rogers_height*mm - 0.5*PCB_Cu_height*mm;  
  //double posY_PCB   = -85.0*mm;    // @the 6-Cu layers

  G4Box* PCB_Rogers_solid = new G4Box("PCB_Rogers",0.5*PCB_width,0.5*PCB_Rogers_height,0.5*PCB_length); 
  G4Box* PCB_Cu_solid = new G4Box("PCB_Cu",0.5*PCB_width,0.5*PCB_Cu_height,0.5*PCB_length); 
  
  G4LogicalVolume* PCB_Rogers_vol = new G4LogicalVolume(PCB_Rogers_solid, Rogers_material,"PCB_Rogers_log",0,0,0);
  G4LogicalVolume* PCB_Cu_vol     = new G4LogicalVolume(PCB_Cu_solid, Cu_material,"PCB_Cu_log",0,0,0);
  
  PCB_Rogers_vol->SetVisAttributes(m_VisRogers4003C);
  Tv.setY(posY_PCB_Rogers); 
  m_EpicVolume->AddPlacedVolume(PCB_Rogers_vol, Tv, Rv);
  cpt_AddVolumes++; 
  PCB_Cu_vol->SetVisAttributes(m_VisCu);
  Tv.setY(posY_PCB_Cu);
  m_EpicVolume->AddPlacedVolume(PCB_Cu_vol, Tv, Rv);
  cpt_AddVolumes++; 

  // --- FLANGE IN ALUMINIUM
  double flange_full_width  = PCB_width;
  double flange_full_length = PCB_length;
  double flange_full_height =   5.*mm;       double posY_flange = posY_PCB_Cu - 0.5*PCB_Cu_height - 0.5*flange_full_height ; 
  double flange_open_width   = 150.0*mm; 
  double flange_open_height  =   5.1*mm; 
  double flange_open1_length =  74.5*mm;     double posZ_open1 = -89.5*mm;
  double flange_open2_length =  85.0*mm;     double posZ_open2 =   0.0*mm;
  double flange_open3_length =  42.0*mm;     double posZ_open3 =  73.5*mm;

  G4Box* flange_full  = new G4Box("flange_full" , 0.5*flange_full_width , 0.5*flange_full_height, 0.5*flange_full_length);
  G4Box* flange_open1 = new G4Box("flange_open1", 0.5*flange_open_width , 0.5*flange_open_height, 0.5*flange_open1_length);
  G4Box* flange_open2 = new G4Box("flange_open1", 0.5*flange_open_width , 0.5*flange_open_height, 0.5*flange_open2_length);
  G4Box* flange_open3 = new G4Box("flange_open1", 0.5*flange_open_width , 0.5*flange_open_height, 0.5*flange_open3_length);

  G4VSolid* flange_int1  = (G4VSolid*) new G4SubtractionSolid("flange_int1" ,flange_full,flange_open1,0,G4ThreeVector(0,0,posZ_open1));
  G4VSolid* flange_int2  = (G4VSolid*) new G4SubtractionSolid("flange_int2" ,flange_int1,flange_open2,0,G4ThreeVector(0,0,posZ_open2));
  G4VSolid* flange_final = (G4VSolid*) new G4SubtractionSolid("flange_final",flange_int2,flange_open3,0,G4ThreeVector(0,0,posZ_open3));

  G4LogicalVolume* flange_vol = new G4LogicalVolume(flange_final, Al_material, "flange_log", 0,0,0);
  flange_vol->SetVisAttributes(m_VisAl);
  Tv.setY(posY_flange);
  m_EpicVolume->AddPlacedVolume(flange_vol, Tv, Rv);
  cpt_AddVolumes++; 


  // --- FRAME IN ALUMINIUM
  double frame_foot_full_width  = flange_full_width;
  double frame_foot_full_length = flange_full_length;
  double frame_foot_full_height = 2.*mm;                
  double posY_frame_foot = posY_PCB_Rogers + 0.5*PCB_Rogers_height + 0.5*frame_foot_full_height ; 
  double frame_foot_open_width  = flange_full_width  - 2.*13.*mm;
  double frame_foot_open_length = flange_full_length - 2.*13.*mm;
  double frame_foot_open_height = 2.05*mm;   
  
  G4Box* frame_foot_full = new G4Box("frame_foot_full", 0.5*frame_foot_full_width, 0.5*frame_foot_full_height, 0.5*frame_foot_full_length);
  G4Box* frame_foot_open = new G4Box("frame_foot_open", 0.5*frame_foot_open_width, 0.5*frame_foot_open_height, 0.5*frame_foot_open_length);
 
  G4VSolid* frame_foot = (G4VSolid*) new G4SubtractionSolid("frame_int1",frame_foot_full,frame_foot_open,0,G4ThreeVector(0,0,0));
 
  G4LogicalVolume* frame_foot_final_vol = new G4LogicalVolume(frame_foot, Al_material, "frame_foot_log", 0,0,0);
  frame_foot_final_vol->SetVisAttributes(m_VisAl);
  Tv.setY(posY_frame_foot);
  m_EpicVolume->AddPlacedVolume(frame_foot_final_vol, Tv, Rv);
  cpt_AddVolumes++; 

  double frame_full_width  = frame_foot_open_width;
  double frame_full_length = frame_foot_open_length; 
  double frame_full_height = frame_foot_full_height + 155.5*mm; 
  double posY_frame = 0.5*frame_full_height - 75.*mm; 
  double frame_open_width  = frame_full_width  - 2.*2.*mm; 
  double frame_open_length = frame_full_length - 2.*2.*mm; 
  double frame_open_height = frame_full_height - 2.*mm; 
  double frame_open_side_width  = frame_foot_full_width;
  double frame_open_side_length = frame_full_length - 2.*15.*mm;
  double frame_open_side_height = frame_full_height - 2.*15.*mm;
  double frame_open_front_width  = frame_full_width - 2.*15.*mm;
  double frame_open_front_length = frame_foot_full_length ;
  double frame_open_front_height = frame_full_height - 2.*15.*mm;
  double frame_open_top_width  = frame_full_width  - 2.*15.*mm;
  double frame_open_top_length = frame_full_length - 2.*15.*mm;
  double frame_open_top_height = 2.1*mm;
 
  G4Box* frame_full  = new G4Box("frame_full" , 0.5*frame_full_width, 0.5*frame_full_height, 0.5*frame_full_length);
  G4Box* frame_open1 = new G4Box("frame_open1", 0.5*frame_open_width, 0.5*frame_open_height, 0.5*frame_open_length);
  G4Box* frame_open2 = new G4Box("frame_open2", 0.5*frame_open_side_width, 0.5*frame_open_side_height, 0.5*frame_open_side_length);
  G4Box* frame_open3 = new G4Box("frame_open3", 0.5*frame_open_front_width, 0.5*frame_open_front_height, 0.5*frame_open_front_length);
  G4Box* frame_open4 = new G4Box("frame_open4", 0.5*frame_open_top_width, 0.5*frame_open_top_height, 0.5*frame_open_top_length);
  
  G4VSolid* frame_int1 = (G4VSolid*) new G4SubtractionSolid("frame_int1",frame_full,frame_open1,0,G4ThreeVector(0,-1.01*mm,0));
  G4VSolid* frame_int2 = (G4VSolid*) new G4SubtractionSolid("frame_int2",frame_int1,frame_open2,0,G4ThreeVector(0,0,0));
  G4VSolid* frame_int3 = (G4VSolid*) new G4SubtractionSolid("frame_int3",frame_int2,frame_open3,0,G4ThreeVector(0,0,0));
  G4VSolid* frame_int4 = (G4VSolid*) new G4SubtractionSolid("frame_int4",frame_int3,frame_open4,0,G4ThreeVector(0,0.5*frame_full_height-1.*mm,0));
 
  G4LogicalVolume* frame_final_vol = new G4LogicalVolume(frame_int4, Al_material, "frame_log", 0,0,0);
  frame_final_vol->SetVisAttributes(m_VisAl);
  Tv.setY(posY_frame);
  m_EpicVolume->AddPlacedVolume(frame_final_vol, Tv, Rv);
  cpt_AddVolumes++; 

  //--- TITANE WINDOWS
  double window_front_width  = 142.0*mm;
  double window_front_length =   0.1*mm;
  double window_front_height = 150.0*mm;
  G4Box* window_front = new G4Box("window_front", 0.5*window_front_width, 0.5*window_front_height, 0.5*window_front_length);
  G4LogicalVolume* window_front_vol = new G4LogicalVolume(window_front, Ti_material,"window_front_log",0,0,0);
  window_front_vol->SetVisAttributes(m_VisTi);
  Tv.setY(posY_frame);
  Tv.setZ(-0.5*frame_open_length*mm + 0.5*window_front_length);
  m_EpicVolume->AddPlacedVolume(window_front_vol, Tv, Rv);
  cpt_AddVolumes++; 
  Tv.setZ(0.5*frame_open_length*mm - 0.5*window_front_length);
  m_EpicVolume->AddPlacedVolume(window_front_vol, Tv, Rv);
  cpt_AddVolumes++; 

  double window_side_width  =   0.1*mm;
  double window_side_length = 290.0*mm;
  double window_side_height = 150.0*mm;
  G4Box* window_side = new G4Box("window_side", 0.5*window_side_width, 0.5*window_side_height, 0.5*window_side_length);
  G4LogicalVolume* window_side_vol = new G4LogicalVolume(window_side, Ti_material,"window_side_log",0,0,0);
  window_side_vol->SetVisAttributes(m_VisTi);
  Tv.setX(0.5*frame_open_width*mm - 0.5*window_side_width);
  Tv.setZ(0);
  m_EpicVolume->AddPlacedVolume(window_side_vol, Tv, Rv);
  cpt_AddVolumes++; 
  Tv.setX(-0.5*frame_open_width*mm + 0.5*window_side_width);
  m_EpicVolume->AddPlacedVolume(window_side_vol, Tv, Rv);
  cpt_AddVolumes++; 

  double window_top_width  = 142.0*mm;
  double window_top_length = 290.0*mm;
  double window_top_height =   0.1*mm;
  G4Box* window_top = new G4Box("window_top", 0.5*window_top_width, 0.5*window_top_height, 0.5*window_top_length);
  G4LogicalVolume* window_top_vol = new G4LogicalVolume(window_top, Ti_material,"window_top_log",0,0,0);
  window_top_vol->SetVisAttributes(m_VisTi);
  Tv.setX(0);
  Tv.setY(posY_frame + 0.5*frame_open_height - 1.*mm - 0.5*window_top_height);
  Tv.setZ(0);
  m_EpicVolume->AddPlacedVolume(window_top_vol, Tv, Rv);
  cpt_AddVolumes++; 



  //--- CENTRAL PART OF THE FISSION CHAMBER
  //  - Cathodes with/without the actinide deposit 
  //  - Anodes with the gas volume upstream and downstream
  //  - gas volume between the pair of cathodes (back-to-back)
  Tv.setX(0); Tv.setY(0); Tv.setZ(0);

  // definition of the gas volume located @ pair of cathodes
  G4Tubs* gas_KK_solid = new G4Tubs("gas_KK",0,37.*mm,0.5*m_InterDistance_KK*mm,0,360*deg);
  G4LogicalVolume* gas_KK_vol = new G4LogicalVolume(gas_KK_solid, m_GasMaterial,"gas_KK_log",0,0,0);
  gas_KK_vol->SetVisAttributes(m_VisGas);

  // position along the beam axis of the first anode
  double thickness_A = 2*Epic_NS::Cu_Thickness*mm + Epic_NS::Kapton_Thickness*mm;
  double posZ_anode = (1-m_nA)*(m_Distance_AK*mm + m_Thickness_K*mm) - std::trunc(0.5*m_nA)*(m_InterDistance_KK*mm + thickness_A); 
  double posZ_first_cathode = posZ_anode - 0.5*thickness_A*mm - m_Distance_AK*mm - 0.5*m_Thickness_K*mm;

  // build the stack of cathodes / actinide samples / anodes to do build sensitive gas
  G4UserLimits * user_limits = new G4UserLimits(m_DriftStepLength, 33.*mm, 1.*us, 0, 0);
  for(int i=0; i<m_nA; i++){

    // cathode upstream the anode i with the samples
    BuildCathode(posZ_anode - 0.5*thickness_A*mm - m_Distance_AK*mm - 0.5*m_Thickness_K*mm);
    cpt_AddVolumes++; 
    BuildSample(posZ_anode - 0.5*thickness_A*mm - m_Distance_AK*mm + 0.5*m_SampleThickness.at(i)*mm,i);
    cpt_AddVolumes++; 

    // anode i
    BuildAnode(posZ_anode);
    cpt_AddVolumes+=3; 
    
    // cathode downstream anode i with or without the samples with the corresponding sensitive gas thickness
    if(m_nSamplesPerA.at(i)==1){ 
      // sensitive gas upstream anode i only
      // FIXME : attention gas should surround the sample
      G4Tubs* gas_KA = new G4Tubs("gas_KA",0,m_RadiusAnode,0.5*(m_Distance_AK*mm-m_SampleThickness.at(i)*mm),0,360*deg);
      G4LogicalVolume* gas_KA_vol = new G4LogicalVolume(gas_KA, m_GasMaterial,"gas_KA_log",0,0,0);
      gas_KA_vol->SetSensitiveDetector(m_EpicScorer);
      gas_KA_vol->SetVisAttributes(m_VisGasSensitive);
      gas_KA_vol->SetUserLimits(user_limits);
      Tv.setZ(posZ_anode - 0.5*thickness_A - 0.5*(m_Distance_AK*mm-m_SampleThickness.at(i)*mm));
      m_EpicVolume->AddPlacedVolume(gas_KA_vol, Tv, Rv);
      cpt_AddVolumes++; 
      m_mapping_A[cpt_AddVolumes] = i;
      
      // if segemented anode, sensitive gas upstream external anode i
      if (m_SegmentedAnode){
      	G4Tubs* gas_KAfull = new G4Tubs("gas_KAfull",0,37.*mm,0.5*m_Distance_AK*mm,0,360*deg);
      	G4Tubs* openKAint  = new G4Tubs("openKAint",0,m_RadiusAnode+1.*mm,0.5*m_Distance_AK*mm,0,360*deg);
      	G4VSolid* gas_KAext = (G4VSolid*) new G4SubtractionSolid("gas_KAext",gas_KAfull,openKAint,0,G4ThreeVector(0,0,0));
      	G4LogicalVolume* gas_KAext_vol = new G4LogicalVolume(gas_KAext, m_GasMaterial,"gas_KAext_log",0,0,0);
      	gas_KAext_vol->SetSensitiveDetector(m_EpicScorer);
      	gas_KAext_vol->SetVisAttributes(m_VisGasSensitive);
      	gas_KAext_vol->SetUserLimits(user_limits);
      	Tv.setZ(posZ_anode - 0.5*thickness_A - 0.5*m_Distance_AK*mm);
      	m_EpicVolume->AddPlacedVolume(gas_KAext_vol, Tv, Rv);
      	cpt_AddVolumes++; 
      	m_mapping_A[cpt_AddVolumes] = 100+i;
      } 
      
      // no sample on the seconde cathode downstream anode i
      G4Tubs* gas_AK = new G4Tubs("gas_AK",0,m_RadiusAnode,0.5*m_Distance_AK*mm,0,360*deg);
      G4LogicalVolume* gas_AK_vol = new G4LogicalVolume(gas_AK, m_GasMaterial,"gas_A_log",0,0,0);
      gas_AK_vol->SetVisAttributes(m_VisGas);
      Tv.setZ(posZ_anode + 0.5*thickness_A + 0.5*m_Distance_AK*mm);
      m_EpicVolume->AddPlacedVolume(gas_AK_vol, Tv, Rv);
      cpt_AddVolumes++; 
    }
    else if(m_nSamplesPerA.at(i)==2){
      // sensitive gas upstream and downstream segmented anode i
      // FIXME : attention gas should surround the sample
      G4Tubs* gas   = new G4Tubs("gas_around_A",0,m_RadiusAnode,0.5*thickness_A+m_Distance_AK*mm-m_SampleThickness.at(i)*mm,0,360*deg);
      G4Tubs* openA = new G4Tubs("openA",0,38.*mm,0.5*thickness_A,0,360*deg);
      // TO CHECK : 
      //G4Tubs* openS = new G4Tubs("openS",0,33.*mm,m_SampleThickness.at(i)*mm,0,360*deg);
      //G4VSolid* gas_openA   = (G4VSolid*) new G4SubtractionSolid("gas_openA",gas,openA,0,G4ThreeVector(0,0,0));
      //G4VSolid* gas_openS_1 = (G4VSolid*) new G4SubtractionSolid("gas_openS1",gas_openA,openS,0,G4ThreeVector(0,0,-0.5*thickness_A-m_Distance_AK*mm));
      //G4VSolid* gas_A       = (G4VSolid*) new G4SubtractionSolid("gas_A",gas_openS_1,openS,0,G4ThreeVector(0,0,0.5*thickness_A+m_Distance_AK*mm));
      G4VSolid* gas_A = (G4VSolid*) new G4SubtractionSolid("gas_A",gas,openA,0,G4ThreeVector(0,0,0));
      G4LogicalVolume* gas_A_vol = new G4LogicalVolume(gas_A, m_GasMaterial,"gas_A_log",0,0,0);
      gas_A_vol->SetSensitiveDetector(m_EpicScorer);
      gas_A_vol->SetVisAttributes(m_VisGasSensitive);
      gas_A_vol->SetUserLimits(user_limits);
      Tv.setZ(posZ_anode);
      m_EpicVolume->AddPlacedVolume(gas_A_vol, Tv, Rv);
      cpt_AddVolumes++; 
      m_mapping_A[cpt_AddVolumes] = i;
      
      // if segmented anode sensitive gas upstream and downstream segmented external anode i
      if (m_SegmentedAnode){
	      G4Tubs* gas_Afull = new G4Tubs("gas_around_Afull",0,37.*mm,0.5*thickness_A+m_Distance_AK*mm,0,360*deg);
      	G4Tubs* openAint  = new G4Tubs("openAint",0,m_RadiusAnode+1.*mm,thickness_A+m_Distance_AK*mm,0,360*deg);
      	G4VSolid* gas_full_Aext = (G4VSolid*) new G4SubtractionSolid("gas_full_Aext",gas_Afull,openAint,0,G4ThreeVector(0,0,0));
      	G4VSolid* gas_Aext = (G4VSolid*) new G4SubtractionSolid("gas_Aext",gas_full_Aext,openA,0,G4ThreeVector(0,0,0));
      	G4LogicalVolume* gas_Aext_vol = new G4LogicalVolume(gas_Aext, m_GasMaterial,"gas_KAext_log",0,0,0);
      	gas_Aext_vol->SetSensitiveDetector(m_EpicScorer);
      	gas_Aext_vol->SetVisAttributes(m_VisGasSensitive);
      	gas_Aext_vol->SetUserLimits(user_limits);
      	Tv.setZ(posZ_anode);
      	m_EpicVolume->AddPlacedVolume(gas_Aext_vol, Tv, Rv);
      	cpt_AddVolumes++; 
      	m_mapping_A[cpt_AddVolumes] = 100+i;
      } 

      // sample on the second cathode
      BuildSample(posZ_anode + 0.5*thickness_A + m_Distance_AK*mm - 0.5*m_SampleThickness.at(i)*mm,i);
      cpt_AddVolumes++; 
    }
    else{
        cout << "nSamplesPerAnode should be 1 or 2 only" << endl;
    }
    BuildCathode(posZ_anode + 0.5*thickness_A + m_Distance_AK*mm + 0.5*m_Thickness_K*mm);
    cpt_AddVolumes++; 
    // gas volume to the next cathode
    if(i < m_nA-1) {
      Tv.setZ(posZ_anode + 0.5*thickness_A + m_Distance_AK*mm + m_Thickness_K*mm + 0.5*m_InterDistance_KK*mm);
      m_EpicVolume->AddPlacedVolume(gas_KK_vol, Tv, Rv);
      cpt_AddVolumes++; 
    }
    posZ_anode += 2.*(m_Distance_AK*mm + m_Thickness_K*mm) + m_InterDistance_KK*mm + thickness_A;
  }
  double posZ_last_cathode = posZ_anode - 0.5*thickness_A - m_Distance_AK*mm - m_InterDistance_KK*mm - 1.5*m_Thickness_K*mm;


  // --- UNACTIVE GAS VOLUME
  //     same size than frame_open subtracted by the "active" central part
  double gas_full_width  = frame_open_width  ;
  double gas_full_length = frame_open_length ;  
  double gas_full_height = frame_open_height ; 
  G4Box*  gas_full = new G4Box("gas_full", 0.5*gas_full_width, 0.5*gas_full_height, 0.5*gas_full_length);
  G4Tubs* gas_open = new G4Tubs("gas_open_central",0,37.*mm,0.5*(posZ_last_cathode-posZ_first_cathode+2*m_Thickness_K*mm),0,360*deg);
  G4VSolid* gas_ext = (G4VSolid*) new G4SubtractionSolid("gas_ext",gas_full,gas_open,0,G4ThreeVector(0,75.*mm-0.5*gas_full_height,0));
  G4LogicalVolume* gas_ext_vol = new G4LogicalVolume(gas_ext, m_GasMaterial, "gas_ext_log", 0,0,0);
  gas_ext_vol->SetVisAttributes(m_VisGas);
  Tv.setX(0); Tv.setY(0.5*gas_full_height - 75.*mm); Tv.setZ(0);
  m_EpicVolume->AddPlacedVolume(gas_ext_vol, Tv, Rv);
  cpt_AddVolumes++; 


  // --- BUILD THE GAUSSIAN DISTRIBUTION FOR CONVOLUTION 
  // ---> time resolution (Gas+electronics)
  // ---> should be done once, depends on some parameters in config file
  double normalization_factor = 0;
  for(int bin=0; bin<m_Gauss_nbins; bin++){
    double gauss_val = TMath::Gaus((-Epic_NS::Gauss_HalfRange+bin*m_DriftStepTime)/picosecond,0,Epic_NS::Gauss_Sigma/picosecond,kTRUE); 
    m_Gauss_Distribution.push_back(gauss_val);
    normalization_factor += gauss_val;
  }
  for(int bin=0; bin<m_Gauss_nbins; bin++) m_Gauss_Distribution[bin] /= normalization_factor;


  return m_EpicVolume;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Epic::BuildCathode(double Zpos){
  // Al plate: 
  //   thickness = 12 um
  //   external radius = 37 mm
  G4Tubs* K_solid = new G4Tubs("Cathode",0,37.*mm,0.5*m_Thickness_K*mm,0,360*deg);
  G4LogicalVolume* K_vol = new G4LogicalVolume(K_solid, m_KMaterial,"Al_K_log",0,0,0);
  K_vol->SetVisAttributes(m_VisAl);

  G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(Zpos);
  m_EpicVolume->AddPlacedVolume(K_vol, Tv, Rv);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Epic::BuildSample(double Zpos, int indexA){
  // Samples: 
  //   diameter = 33 mm --> external radius = 16.5 mm
  
  G4Tubs* sample_solid = new G4Tubs(m_SampleMaterial.at(indexA),0,16.5*mm,0.5*m_SampleThickness.at(indexA)*mm,0,360*deg);
  G4Material* sample_material = MaterialManager::getInstance()->GetMaterialFromLibrary(m_SampleMaterial.at(indexA));
  G4LogicalVolume* sample_vol = new G4LogicalVolume(sample_solid, sample_material,"sample_log",0,0,0);
  if(m_SampleMaterial.at(indexA)=="235U")   sample_vol->SetVisAttributes(m_VisSample235U);
  if(m_SampleMaterial.at(indexA)=="238U")   sample_vol->SetVisAttributes(m_VisSample238U);
  if(m_SampleMaterial.at(indexA)=="252Cf")  sample_vol->SetVisAttributes(m_VisSample252Cf);
  if(m_SampleMaterial.at(indexA)=="238Pu" || 
     m_SampleMaterial.at(indexA)=="240Pu" ||
     m_SampleMaterial.at(indexA)=="242Pu") sample_vol->SetVisAttributes(m_VisSamplePu);

  G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(Zpos);
  m_EpicVolume->AddPlacedVolume(sample_vol, Tv, Rv);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Epic::BuildAnode(double Zpos){

  G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
  G4ThreeVector Tv;
  Tv.setX(0); Tv.setY(0); Tv.setZ(0);

  // Get Material
  G4Material* Cu_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
  G4Material* Kapton_material = MaterialManager::getInstance()->GetMaterialFromLibrary("Kapton");

  // Cu plate: 17 um
  G4Tubs* Cu_plate_solid = new G4Tubs("Cu_plate",0,37.*mm,0.5*Epic_NS::Cu_Thickness*mm,0,360*deg);
  G4LogicalVolume* Cu_vol = new G4LogicalVolume(Cu_plate_solid, Cu_material,"Cu_A_log",0,0,0);
  Cu_vol->SetVisAttributes(m_VisCu);

  // Kapton: 50 um
  G4Tubs* Kapton_solid = new G4Tubs("Kapton",0,37.*mm,0.5*Epic_NS::Kapton_Thickness*mm,0,360*deg);
  G4LogicalVolume* Kapton_vol = new G4LogicalVolume(Kapton_solid, Kapton_material,"Kapton_A_log",0,0,0);
  Kapton_vol->SetVisAttributes(m_VisFCWall);

  // Build
  Tv.setZ(Zpos-0.5*Epic_NS::Kapton_Thickness*mm-0.5*Epic_NS::Cu_Thickness*mm);
  m_EpicVolume->AddPlacedVolume(Cu_vol, Tv, Rv);
  Tv.setZ(Zpos);
  m_EpicVolume->AddPlacedVolume(Kapton_vol, Tv, Rv);
  Tv.setZ(Zpos+0.5*Epic_NS::Kapton_Thickness*mm+0.5*Epic_NS::Cu_Thickness*mm);
  m_EpicVolume->AddPlacedVolume(Cu_vol, Tv, Rv);
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Epic::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Epic")){
    pTree->Branch("Epic", "TEpicData", &m_Event) ;
  }
  pTree->SetBranchAddress("Epic", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAction
void Epic::ReadSensitive(const G4Event* ){
 
  //cout << "Enter in Epic::ReadSensitive" << endl; 
  m_Event->Clear();


  ///////////
  // GaseousDetector scorer
  GaseousDetectorScorers::PS_GaseousDetector* Scorer= (GaseousDetectorScorers::PS_GaseousDetector*) m_EpicScorer->GetPrimitive(0);
  /////////////
  //// Interaction scorer
  //InteractionScorers::PS_Interactions* InteractScorer= (InteractionScorers::PS_Interactions*) m_EpicScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  //unsigned int interact_size = InteractScorer->GetMult();

  //if(size != interact_size){
   // cout << "size = " << size << ", interact_size = " << interact_size << endl;
  //}

  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
    vector<string> step_name = Scorer->GetParticleName(i);
    vector<int>    step_trackID = Scorer->GetTrackID(i);
    vector<double> step_vPosZ = Scorer->GetStepPosZ(i);
    vector<double> step_vDE = Scorer->GetEnergyLossPerStep(i);
    vector<double> step_vTime = Scorer->GetStepTime(i);
    double interactionX = Scorer->GetXpos(i); 
    double interactionY = Scorer->GetYpos(i);
    double interactionZ = Scorer->GetZpos(i);


    int Anode = level[0];
    double thickness_A = 2*Epic_NS::Cu_Thickness*mm + Epic_NS::Kapton_Thickness*mm;
    double posZ_anode = (1-m_nA)*(m_Distance_AK*mm + m_Thickness_K*mm) - std::trunc(0.5*m_nA)*(m_InterDistance_KK*mm + thickness_A); 
    if(m_mapping_A[Anode]<100) 	posZ_anode +=  m_mapping_A[Anode]      * (2. * (m_Distance_AK*mm + m_Thickness_K*mm) + m_InterDistance_KK*mm + thickness_A) ;
    else                 	posZ_anode += (m_mapping_A[Anode]-100) * (2. * (m_Distance_AK*mm + m_Thickness_K*mm) + m_InterDistance_KK*mm + thickness_A) ;
      
    // === for Epic Data
    int    trackID            = step_trackID[0] ;
    int    previous_trackID   = trackID         ;
    string name_pertrackID    = step_name[0]    ;
    bool   end_of_new_trackID = false           ;
    double influence_perAnode = 0. ;

    double influence_pertrackID ;
    double detot_pertrackID     ;
    bool   valid_step           ;
    double dz_anode             ;
    vector<double> dz ;
    vector<double> de ;
    vector<double> dt ;
    int current_step = 0;
  
    vector<double> influence(m_Convo_nbins,0);

    if(name_pertrackID=="e-" || name_pertrackID=="anti_nu_e" || name_pertrackID=="gamma"){ 
      valid_step = false;
      influence_pertrackID = 0;
      detot_pertrackID = 0;
    }
    else{
      valid_step = true;
      dz_anode = TMath::Abs(step_vPosZ[0] - posZ_anode) - 0.5*thickness_A ; 
      influence_pertrackID = (step_vDE[0] * dz_anode) / m_Distance_AK*mm  ;
      detot_pertrackID = step_vDE[0] ;
      dz.push_back(dz_anode); 
      de.push_back(step_vDE[0]/eV);
      dt.push_back(step_vTime[0]);
    }


    for(int j=1; j<step_name.size(); j++){
      // FIXME
      if (m_mapping_A[Anode] == 4) {
	        cout << " j: " << j << " bad anode" << endl;
         continue;
      }
      if(step_trackID[j] != previous_trackID) 
        end_of_new_trackID = true;
      else{
        if(valid_step){
          dz_anode = TMath::Abs(step_vPosZ.at(j) - posZ_anode) - 0.5*thickness_A ; // need to remove half thickness of Anode
          if(dz_anode<m_Distance_AK){
            influence_pertrackID += (step_vDE.at(j) * dz_anode) / m_Distance_AK*mm;
	          detot_pertrackID += step_vDE.at(j); 
            dz.push_back(dz_anode); 
            de.push_back(step_vDE.at(j)/eV);
            dt.push_back(step_vTime.at(j));
          }
        }// end of fill vectors if valid_step
        current_step = j;
        if(j==step_name.size()-1) end_of_new_trackID = true;
      }

      if(end_of_new_trackID){
        if(dz.size()>0 && valid_step) {
          Propagate(dz, de, dt, detot_pertrackID,influence);
          influence_perAnode += influence_pertrackID;
        }
        // initialization for the next event [begining]
        // need to be done for any particle
        end_of_new_trackID = false;
        dz.clear(); 
        de.clear(); 
        dt.clear();
        if(current_step<step_name.size()-1){
          name_pertrackID = step_name.at(j);
	        trackID         = step_trackID.at(j);
          if(name_pertrackID=="e-" || name_pertrackID=="anti_nu_e" || name_pertrackID=="gamma") 
            valid_step = false;
          else
            valid_step = true;
          if(valid_step){
            dz_anode = TMath::Abs(step_vPosZ.at(j) - posZ_anode) - 0.5*thickness_A ; 
            if(dz_anode<m_Distance_AK){
	            // initialization for the next event [end]
              influence_pertrackID = (step_vDE.at(j) * dz_anode) / m_Distance_AK*mm;
	            detot_pertrackID     = step_vDE.at(j); 
              dz.push_back(dz_anode); 
              de.push_back(step_vDE.at(j)/eV);
              dt.push_back(step_vTime.at(j));
            }
            current_step = j;
            // in case the new track has one step only, the last one
            if((current_step==step_name.size()-1) && (dz.size()==1) ){
              Propagate(dz, de, dt, detot_pertrackID,influence);
              influence_perAnode += influence_pertrackID;
            }
          }
        }
      }// end of if end_of new trackID
      previous_trackID = step_trackID[j];
    }// end of loop over j


    // Apply time resolution, find Qmax, integrate +/- 30 ns around Qmax 
    double Qmax = 0. ;
    int    kmax = 0  ;
    vector<double> convolution(m_Convo_nbins,0.);
    for(int k=0; k<(m_Convo_nbins-m_Gauss_nbins)-1; k++){
      for(int l=0; l<m_Gauss_nbins; l++)    convolution[k+0.5*m_Gauss_nbins] += influence[k+l] * m_Gauss_Distribution[l];
      if(convolution[k+0.5*m_Gauss_nbins]>Qmax) {
	      Qmax = convolution[k+0.5*m_Gauss_nbins];
        kmax = k+0.5*m_Gauss_nbins;
      }
    }
    Qmax = 0 ;
    if (kmax > 0 ) 
    	for(int k = kmax - 0.5*m_Gauss_nbins; k < kmax + 0.5*m_Gauss_nbins; k++)
    	   Qmax += convolution[k];
    	    
    
    // Fill TEpicData
    // FIXME : need to get the interaction position in the sample
    // interactionX from the TInteractionScorers, need to find the proper index
    // interactionY
    // interactionZ
    m_Event->Set(m_mapping_A[Anode], influence_perAnode, Qmax, interactionX, interactionY, interactionZ);
    influence.clear();
  
  }// end of loop over i
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Epic::Propagate(vector<double> dz, vector<double>de, vector<double>dt, double detot, vector<double>& influence){

  // ==========================
  // for energy loss resolution
  double detot_wReso = RandGauss::shoot(detot,Epic_NS::ResoEnergyPerCent*detot);
  double weight_energy_reso = 1;
  if(detot!=0) weight_energy_reso = detot_wReso / detot;
  
  // =================================================
  // vector with fix-bin length steps along the z axis
  int num_bins = (int)(std::ceil( (m_Distance_AK-dz.at(dz.size()-1)) / (m_DriftStepLength) )) ;
  vector<double> rebinned_de(num_bins, 0.0);
  vector<double> rebinned_dt(num_bins, 0.0); 
  vector<double> rebinned_mult(num_bins, 0.0);
  vector<double> mult_electrons_per_step(num_bins, 0.0);
  double bin_Limits[num_bins+1];
  for(int bin=0; bin<=num_bins; bin++) bin_Limits[bin] = m_Distance_AK - bin * m_DriftStepLength;

  // ==================================================
  // convert vectors of G4Steps to vectors in UserSteps 

  // first G4Step is always included in the first fix_length user bin 
  rebinned_de[0] += de[0];
  rebinned_dt[0] += dt[0]; 
  rebinned_mult[0] += 1.;
          
  // loop over the G4Steps and assign the values to the fix-dz-length user bins
  int bin_current = 0;
  for(int step=1; step<dz.size()-1; step++){
    double dzmax = 0.5*(dz[step-1] + dz[step]);
    double dzmin = 0.5*(dz[step+1] + dz[step]);
    for(int bin=bin_current; bin<num_bins; bin++){
      bin_current = bin;
      if(bin_Limits[bin] >= dzmax && dzmax >  bin_Limits[bin+1] && bin_Limits[bin] >  dzmin && dzmin >= bin_Limits[bin+1]){
  	    rebinned_de[bin] += de[step];
        rebinned_dt[bin] += dt[step]; 
        rebinned_mult[bin] += 1.;
        break;
      }
      else if(bin<num_bins-1 && bin_Limits[bin] >= dzmax && dzmax >  bin_Limits[bin+1] && bin_Limits[bin+1] > dzmin && dzmin >= bin_Limits[bin+2]){
        double weight = (dzmax-bin_Limits[bin+1])/(dzmax-dzmin);
        rebinned_de[bin] += weight * de[step];   rebinned_de[bin+1] += (1.-weight) * de[step]; 
        rebinned_dt[bin] += weight * dt[step];   rebinned_dt[bin+1] += (1.-weight) * dt[step];
        rebinned_mult[bin] += weight ;           rebinned_mult[bin+1] += 1. - weight ;      
        break;
      }
	  }
	}

  // last G4Step is always included in the last fix_length user bin
  rebinned_de[num_bins-1] += de[dz.size()-1]; 
  rebinned_dt[num_bins-1] += dt[dz.size()-1]; 
  rebinned_mult[num_bins-1] += 1.;

  // =====================================================
  // apply time resolution and incremente influence vector
  for(int i=0; i<num_bins; i++){
    
    if(rebinned_mult[i]>0) rebinned_dt[i] = rebinned_dt[i] / rebinned_mult[i] ;
    mult_electrons_per_step[i] = weight_energy_reso * rebinned_de[i]/(m_IonisationEnergy/eV);
    
    int bin = rebinned_dt[i] / m_DriftStepTime;
    for(int step=i; step<ceil(m_Distance_AK/m_DriftStepLength); step++){
      if( (m_Convo_nbins-m_Gauss_nbins) < (bin+step-i) ){
        cout << " m_Convo_nbins = " << m_Convo_nbins << " (influence.size() = " << influence.size() << ") , m_Gauss_nbins = " << m_Gauss_nbins << endl;
        cout << " bin = " << bin << " (from time value) , step = " << step << " [" << i << ":"<< ceil(m_Distance_AK/m_DriftStepLength) <<"]"<< endl;
        cout << " influence["<< bin+step-i+m_Gauss_nbins << "] = " <<  influence[bin+step-i+m_Gauss_nbins] << " : rejected " << endl;
      }
      else{
        influence[bin+step-i+m_Gauss_nbins] += mult_electrons_per_step[i] * m_DriftStepLength / m_Distance_AK;
      }
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Epic::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_EpicScorer = CheckScorer("EpicScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* GaseousDetector= new GaseousDetectorScorers::PS_GaseousDetector("GaseousDetector",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  //and register it to the multifunctionnal detector
  m_EpicScorer->RegisterPrimitive(GaseousDetector);
  m_EpicScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_EpicScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Epic::Construct(){
  return  (NPS::VDetector*) new Epic();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Epic{
    public:
      proxy_nps_Epic(){
        NPS::DetectorFactory::getInstance()->AddToken("Epic","Epic");
        NPS::DetectorFactory::getInstance()->AddDetector("Epic",Epic::Construct);
      }
  };

  proxy_nps_Epic p_nps_Epic;
}
