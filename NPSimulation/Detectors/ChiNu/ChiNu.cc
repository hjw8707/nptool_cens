/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/********************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr *
 *                                                                              *
 * Creation Date  : February 2019                                               *
 * Last update    :                                                             *
 *------------------------------------------------------------------------------*
 * Decription:                                                                  *
 *  This class describes ChiNu simulation                                       *
 *                                                                              *
 *------------------------------------------------------------------------------*
 * Comment:                                                                     *
 *                                                                              *
 ********************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
#include <regex>
//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
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
#include "G4SubtractionSolid.hh"

// NPTool header
#include "ChiNu.hh"
#include "CalorimeterScorers.hh"
#include "ProcessScorers.hh"
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
namespace ChiNu_NS{
  // EJ309 Scintillator - Energy and time Resolution
  const double EnergyThreshold = 0.02*MeV;
  const double ResoTime = 0.686*ns ;
  const double ResoEnergySlope = 0.06514 ;
  const double ResoEnergyOffset = 0.09271*MeV ;
  //  const double ResoEnergy = 0.1*MeV ;
  const double Radius = 8.90*cm ; 
  const double Thickness = 5.08*cm ; 
  const string Material = "EJ309";

  // PMT
  const double PMT_Height = 392*mm; 
  const string PMT_Material = "Al"; 
 
  // Light guide
  const double LG_Rmin1 = 0*mm;
  const double LG_Rmax1 = Radius;
  const double LG_Rmin2 = 0*mm; 
  const double LG_Rmax2 = 50*mm; 
  const double LG_Thickness = 30*mm; 
  const string LG_Material = "PMMA";

  // Pyrex Window
  const double Pyrex_radius = Radius;
  const double Pyrex_thickness = 6.4*mm;
  const string Pyrex_material = "Pyrex"; 

  // Lead shield
  const double Lead_Radius = 9*cm;
  const double Lead_Thickness = 2*mm;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ChiNu Specific Method
ChiNu::ChiNu(){
  m_Event = new TChiNuData() ;
  m_ChiNuScorer = 0;
  m_PMT = 0;
  m_LightGuide = 0;
  m_CylindricalDetector = 0;
  m_LeadShield = 0;
  m_AssemblyVolume = 0;

  m_BuildLeadShield = 0;


  // RGB Color + Transparency
  m_VisCylinder = new G4VisAttributes(G4Colour(0.0, 0.5, 1, 1));   
  m_VisPMT = new G4VisAttributes(G4Colour(0.3, 0.1, 0.1, 0.3));   
  m_VisLightGuide = new G4VisAttributes(G4Colour(0.1,0.5,0.7,1));
  m_VisPyrex = new G4VisAttributes(G4Colour(0.1,0.5,0.7,0.7));
  m_VisLeadShield = new G4VisAttributes(G4Colour(0.2,0.2,0.2,1));
  m_VisFCWall = new G4VisAttributes(G4Colour(0.1,0.5,0.7,1));
  m_VisAl = new G4VisAttributes(G4Colour(0.839,0.803,0.803,1));
  m_VisTi = new G4VisAttributes(G4Colour(0.776,0.662,0.662,0.5));
  m_VisCu = new G4VisAttributes(G4Colour(0.70, 0.40, 0. ,1));
  m_VisRogers4003C = new G4VisAttributes(G4Colour(0.60, 0.60, 0.2 ,1));
}

ChiNu::~ChiNu(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChiNu::AddDetector(G4ThreeVector POS, double TimeResolution){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4 
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_ResoTime.push_back(TimeResolution);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ChiNu::AddDetector(double  R, double  Theta, double  Phi, double TimeResolution){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_ResoTime.push_back(TimeResolution);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4AssemblyVolume* ChiNu::BuildDetector(){
  if(!m_CylindricalDetector){
    m_AssemblyVolume = new G4AssemblyVolume();
    G4RotationMatrix *Rv=new G4RotationMatrix(0,0,0);
    G4ThreeVector Tv;
    Tv.setX(0); Tv.setY(0); Tv.setZ(0);

    // Scintillator
    G4Tubs* tub = new G4Tubs("ChiNu_Cyl",0,ChiNu_NS::Radius,ChiNu_NS::Thickness*0.5,0,360*deg);
    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::Material);
    m_CylindricalDetector = new G4LogicalVolume(tub,DetectorMaterial,"logic_ChiNu_tub",0,0,0);
    m_CylindricalDetector->SetVisAttributes(m_VisCylinder);
    m_CylindricalDetector->SetSensitiveDetector(m_ChiNuScorer);
    m_AssemblyVolume->AddPlacedVolume(m_CylindricalDetector, Tv, Rv);

    // Pyrex Window
    G4Tubs* Pyrex_tub = new G4Tubs("Pyrex_tub",0,ChiNu_NS::Pyrex_radius, ChiNu_NS::Pyrex_thickness*0.5, 0 , 360*deg);
    G4Material* PyrexMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::Pyrex_material);
    G4LogicalVolume* LogPyrex = new G4LogicalVolume(Pyrex_tub, PyrexMaterial,"logic_pyrex",0,0,0);
    LogPyrex->SetVisAttributes(m_VisPyrex);
    Tv.setZ(ChiNu_NS::Thickness*0.5 + ChiNu_NS::Pyrex_thickness*0.5);
    m_AssemblyVolume->AddPlacedVolume(LogPyrex, Tv, Rv);


    // Light guide
    G4Cons* LGsolid = new G4Cons("light_guide", ChiNu_NS::LG_Rmin1, ChiNu_NS::LG_Rmax1, ChiNu_NS::LG_Rmin2, ChiNu_NS::LG_Rmax2, ChiNu_NS::LG_Thickness*0.5,0,360*deg); 
    G4Material* LGMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::LG_Material);
    m_LightGuide = new G4LogicalVolume(LGsolid,LGMaterial,"logic_light_guide",0,0,0);
    m_LightGuide->SetVisAttributes(m_VisLightGuide);
    Tv.setZ(ChiNu_NS::Thickness*0.5 + ChiNu_NS::Pyrex_thickness + ChiNu_NS::LG_Thickness*0.5);
    m_AssemblyVolume->AddPlacedVolume(m_LightGuide, Tv, Rv);

    // PMT
    //G4Tubs* pmt = new G4Tubs("ChiNu_pmt",ChiNu_NS::PMT_InnerDiameter*0.5,ChiNu_NS::PMT_OuterDiameter*0.5,ChiNu_NS::PMT_Thickness*0.5,0,360*deg);
    
    double zplane[4] ={0, 18*cm, 24*cm, ChiNu_NS::PMT_Height};
    double rin[4] = {ChiNu_NS::Radius+0.2*mm, ChiNu_NS::Radius+0.2*mm, 40*mm, 40*mm};
    double rout[4] = {ChiNu_NS::Radius+2.2*mm, ChiNu_NS::Radius+2.2*mm, 42*mm, 42*mm};
    G4Polycone* pmt = new G4Polycone("ChiNu_PMT", 0, 360*deg, 4, zplane, rin, rout);
    G4Material* pmtMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ChiNu_NS::PMT_Material);
    m_PMT = new G4LogicalVolume(pmt,pmtMaterial,"logic_pmt_tub",0,0,0);
    m_PMT->SetVisAttributes(m_VisPMT);
    //Tv.setZ(ChiNu_NS::Thickness*0.5 + ChiNu_NS::LG_Thickness + ChiNu_NS::Pyrex_thickness + ChiNu_NS::PMT_Height*0.5);
    Tv.setZ(-ChiNu_NS::Thickness*0.5);
    m_AssemblyVolume->AddPlacedVolume(m_PMT, Tv, Rv);

    // Lead shield
    if(m_BuildLeadShield){
    	G4Tubs* lead = new G4Tubs("lead_shield", 0, ChiNu_NS::Lead_Radius, ChiNu_NS::Lead_Thickness*0.5, 0, 360*deg);
    	G4Material* LeadMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary("Pb");
	m_LeadShield = new G4LogicalVolume(lead, LeadMaterial, "logic_lead_shield",0,0,0);
	m_LeadShield->SetVisAttributes(m_VisLeadShield);
	Tv.setZ(-ChiNu_NS::Thickness*0.5 - ChiNu_NS::Lead_Thickness*0.5-10*mm);
        m_AssemblyVolume->AddPlacedVolume(m_LeadShield, Tv, Rv);
    }

  }
  return m_AssemblyVolume;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void ChiNu::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ChiNu");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  // In case a file of different values for different detectors is to be applied.
  // string fTimeRes=blocks[0]->GetString("fTimeResolution");
  // vector<double> TimeRes;
  // if(fTimeRes!="")
  //   {
  //     ifstream fres(fTimeRes.c_str());
  //     double T;
  //     while(fres>> T)
  // 	TimeRes.push_back(T);
  //   }
  
  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  ChiNu " << i+1 <<  endl;
    
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));

      AddDetector(Pos,ChiNu_NS::ResoTime);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  ChiNu " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      m_BuildLeadShield = blocks[i]->GetInt("LeadShield");
      double TimeResolution = ChiNu_NS::ResoTime;
      if(blocks[i]->HasToken("TimeResolution"))
	TimeResolution=blocks[i]->GetDouble("TimeResolution","ns");
      
      // if(fTimeRes!="") TimeResolution = TimeRes.at((i+1)/6+9*((i+1)%6));
      
      AddDetector(R,Theta,Phi,TimeResolution);
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
void ChiNu::ConstructDetector(G4LogicalVolume* world){
 
  // double Vm = 0.08206 * 300 * atmosphere / (0.792 * kelvin); // pressure * bar?
  // double density = ( ( 0.21*2*15.99903 + 0.79*2*14.00643 ) / Vm );

  // cout << "  //////////////////  density=" << density << " ///////////////////" << endl;
  // cout << "  atmosphere=" << atmosphere << " kelvin=" << kelvin << " Vm=" << Vm << " mg/cm3=" << mg/cm3 << endl;

  // G4Material* Air = MaterialManager::getInstance()->GetMaterialFromLibrary("Air",density * mg/cm3);
  // world->SetMaterial(Air);

  G4Material* Air = MaterialManager::getInstance()->GetMaterialFromLibrary("Air");
  world->SetMaterial(Air);

  
  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*ChiNu_NS::Thickness*0.5;
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
    BuildDetector()->MakeImprint(world,Det_pos,Rot,i);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void ChiNu::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("ChiNu")){
    pTree->Branch("ChiNu", "TChiNuData", &m_Event) ;
  }
  pTree->SetBranchAddress("ChiNu", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void ChiNu::ReadSensitive(const G4Event* ){
  m_Event->Clear();


  ///////////
  // Process scorer
  ProcessScorers::PS_Process* Process_scorer = (ProcessScorers::PS_Process*) m_ChiNuScorer->GetPrimitive(2);
  unsigned int ProcessMult = Process_scorer->GetMult();
  vector<int> IsGamma;

  int m_VID=-1;
  int g_hit=0;

  for(unsigned int i = 0 ; i < ProcessMult ; i++){
    string particle_name = Process_scorer->GetParticleName(i);
    string volume_name   = Process_scorer->GetVolumeName(i);
    string process_name   = Process_scorer->GetProcessName(i);


    regex pattern("impr_(\\d+)_logic");
    smatch match;
    int VolumeNbr;

    if(regex_search( volume_name,match,pattern))
      VolumeNbr = stoi(match[1]);
    else VolumeNbr = -1;

    if(VolumeNbr!=m_VID && particle_name=="gamma")
      {// if the first particle in the volume is a gamma, it will give a wrong time
	g_hit=1;
      }
    else if(VolumeNbr!=m_VID)
      g_hit=0;
    else if(g_hit==1){
      // Only if an electron ionizes the liquid scintillator will there be a measured energy deposit and associated time.
      if(particle_name=="e-")
	{
	  g_hit++;
	  IsGamma.push_back(m_VID);
	}
    }

    m_VID=VolumeNbr;
  }

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_ChiNuScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  vector<int> DetID;
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i);
    int DetectorNbr = level[0];
    DetID.push_back(DetectorNbr);
    // Light output equation:
    // L(Ep)=A*Ep - B*(1-exp(-C*pow(Ep,D))) from F. Pino et al. Applied Radiation and Isotopes 89 (2014) 79-84.
    double A=0.62, B=1.3, C=0.39, D=0.97;
    double Ep=RandGauss::shoot(Scorer->GetEnergy(i),ChiNu_NS::ResoEnergyOffset+ChiNu_NS::ResoEnergySlope*Scorer->GetEnergy(i));
    double Energy = A*Ep - B*(1-exp(-C*pow(Ep,D)));
    ////////////////////////
 
    
    if(Energy>ChiNu_NS::EnergyThreshold){
      //      double Time = RandGauss::shoot(Scorer->GetTime(i),ChiNu_NS::ResoTime); // One time resolution for all detectors
      double Time = RandGauss::shoot(Scorer->GetTime(i),m_ResoTime.at(DetectorNbr-1)); // Time resolution that can be set different for each detector
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time);

      // Check if the detector was first lighted by a gamma ray, which would have given a wrong time.
      auto it = find(IsGamma.begin(), IsGamma.end(), DetectorNbr);
      if(it != IsGamma.end())
	m_Event->SetIsGamma(DetectorNbr,true);
      else
	m_Event->SetIsGamma(DetectorNbr,false);
    }
  }
  IsGamma.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void ChiNu::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false; 
  m_ChiNuScorer = CheckScorer("ChiNuScorer",already_exist) ;

  if(already_exist) 
    return ;

  // Otherwise the scorer is initialised
  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
  G4VPrimitiveScorer* Process= new ProcessScorers::PS_Process("Process", 0) ;
  //and register it to the multifunctionnal detector
  m_ChiNuScorer->RegisterPrimitive(Calorimeter);
  m_ChiNuScorer->RegisterPrimitive(Interaction);
  m_ChiNuScorer->RegisterPrimitive(Process);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_ChiNuScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* ChiNu::Construct(){
  return  (NPS::VDetector*) new ChiNu();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_ChiNu{
    public:
      proxy_nps_ChiNu(){
        NPS::DetectorFactory::getInstance()->AddToken("ChiNu","ChiNu");
        NPS::DetectorFactory::getInstance()->AddDetector("ChiNu",ChiNu::Construct);
      }
  };

  proxy_nps_ChiNu p_nps_ChiNu;
}
