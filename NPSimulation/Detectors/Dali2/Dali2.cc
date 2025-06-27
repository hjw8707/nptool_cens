/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: E. Tronchin                                              *
 * contact address: elidiano.tronchin@studenti.unipd.it                      *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Dali2 simulation                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <sstream>
#include <cmath>
#include <limits>
using namespace std;

//G4 Geometry object
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4VSolid.hh"
// #ifndef G4UEXTRUDEDSOLID_hh
// #define G4UEXTRUDEDSOLID_hh
// #include "G4USolid.hh"
//  #if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

//#include "G4UExtrudedSolid.hh"
#include "G4TwoVector.hh"
#include "G4TessellatedSolid.hh"

//G4 sensitive
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"

//G4 various object
#include "G4Material.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
//#include "G4VPhysicalVolume.hh"
#include "G4PVReplica.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDParticleFilter.hh"

// NPTool header
#include "Dali2.hh"
#include "CalorimeterScorers.hh"
#include "InteractionScorers.hh"
#include "RootOutput.h"
#include "MaterialManager.hh"
#include "NPSDetectorFactory.hh"
#include "NPOptionManager.h"
#include "NPSHitsMap.hh"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Dali2_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 0*MeV;
  const double ResoTime = 0.0*ns; //4.5*ns ;
  const double ResoEnergy = 0.0311;
  //  const double ResoEnergy = 0.025;
  //Relative resolution DeltaE = 0.122*Sqrt(E) 
  //  const double ResoEnergy = 0.122;  //Relative resolution DeltaE = 0.122*Sqrt(E) 
  /* const double ResoEnergy = 1.36*MeV ; // mean Resolution(FWHM) 1.7% of 80MeV from slides 20170214-SAMURAI34-setup-DALI.pdf if 1.7% of 30MeV = 0.51 MeV // 0.001*MeV ; */
  //const double Radius = 50*mm ; 

  ////////////////////////////////////////
  // Geometry parameters
  double AlHouseWidth[3];
  double AlHouseHight[3];
  double AlHouseDepth[3];
  double AlHouseThick[3];
  double MgOCoatWidth[3];
  double MgOCoatHight[3];
  double MgOCoatDepth[3];
  double MgOCoatThick[3];
  double CrystalWidth[3];
  double CrystalHight[3];
  double CrystalDepth[3];
  ////////////////////////////////////////
  
  void SetParameters() {
    ////////////////////////////////////////
    // Geometry for Saint-Gobain
    AlHouseWidth[0] = 5.014*cm;
    AlHouseHight[0] = 8.520*cm;
    AlHouseDepth[0] = 17.44*cm;
    AlHouseThick[0] = 0.1*cm;
    MgOCoatWidth[0] = 4.780*cm;
    MgOCoatHight[0] = 8.280*cm;
    MgOCoatDepth[0] = 16.28*cm;
    MgOCoatThick[0] = 0.14*cm;
    CrystalWidth[0] = 4.500*cm;
    CrystalHight[0] = 8.000*cm;
    CrystalDepth[0] = 16.00*cm;
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    // Geometry for Scioninx
    AlHouseWidth[1] = 4.500*cm;
    AlHouseHight[1] = 8.300*cm;
    AlHouseDepth[1] = 18.00*cm;
    AlHouseThick[1] = 0.1*cm;
    MgOCoatWidth[1] = 4.300*cm;
    MgOCoatHight[1] = 8.100*cm;
    MgOCoatDepth[1] = 16.40*cm;
    MgOCoatThick[1] = 0.05*cm; // not suitable, actually
    CrystalWidth[1] = 4.000*cm;
    CrystalHight[1] = 8.000*cm;
    CrystalDepth[1] = 16.00*cm;
    ////////////////////////////////////////

    ////////////////////////////////////////
    // Geometry for DALI1
    AlHouseWidth[2] = 6.604*cm;
    AlHouseHight[2] = 6.604*cm;
    AlHouseDepth[2] = 12.93*cm;
    AlHouseThick[2] = 0.05*cm;
    MgOCoatWidth[2] = 6.500*cm;
    MgOCoatHight[2] = 6.500*cm;
    MgOCoatDepth[2] = 12.80*cm;
    MgOCoatThick[2] = 0.2*cm; // not suitable, actually
    CrystalWidth[2] = 6.096*cm;
    CrystalHight[2] = 6.096*cm;
    CrystalDepth[2] = 12.192*cm;
    ////////////////////////////////////////
  }
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Dali Specific Method
Dali2::Dali2(){
  m_Event = new TDali2Data() ;
  m_Dali2SGDetector = 0;
  m_Dali2SCDetector = 0;
  m_Dali1Detector = 0;
  m_DaliScorer = 0;

  Dali2_NS::SetParameters();
 }

Dali2::~Dali2(){
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali2::AddDetector(G4ThreeVector POS, int Type){
  // Convert the POS value to R theta Phi as Cylindrical coordinate is easier in G4 
  m_x.push_back(POS.x());
  m_y.push_back(POS.y());
  m_z.push_back(POS.z());
  m_Type.push_back(Type);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Dali2::AddDetector(G4ThreeVector POS, G4ThreeVector ANG, int Type){
  // Convert the POS value to R theta Phi as Cylindrical coordinate is easier in G4 
  m_x.    push_back(POS.x());
  m_y.    push_back(POS.y());
  m_z.    push_back(POS.z());
  m_psi.  push_back(ANG.x());
  m_theta.push_back(ANG.y());
  m_phi.  push_back(ANG.z());
  m_Type. push_back(Type);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void Dali2::AddDetector(double  R, double  Theta, double  Phi, int Type){
//  double m_r, m_alpha, m_zeta;
//  m_r = R*cos(Phi);
//  m_alpha = Theta;
//  m_zeta = R*sin(Phi);
//  m_R.push_back(m_r);
//  m_Alpha.push_back(m_alpha);
//  m_Zeta.push_back(m_zeta);
//  m_Type.push_back(Type);
//}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Dali2::AddDetector(double  R, double  Alpha, double  Zeta, int Type){
  m_x.push_back(R*cos(Alpha));
  m_y.push_back(R*sin(Alpha));
  m_z.push_back(Zeta);
  m_Type.push_back(Type);}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Definition Materials MgO and NaI(Tl)

 void Dali2::DefinitionMaterials()
 {
  
    G4Element *Tl = new G4Element("Thallium","Tl",81.,204.383*g/mole );

    NaI_Tl = new G4Material("NaI_Tl",3.6667*g/cm3, 2);
    NaI_Tl->AddMaterial(MaterialManager::getInstance()->GetMaterialFromLibrary("NaI"), 99.6*perCent);
    //    NaI_Tl->AddMaterial(MaterialManager::getInstance()->GetMaterialFromLibrary("CsI"), 99.6*perCent);
    NaI_Tl->AddElement(Tl, 0.4*perCent);

    ////////////////////////////////////////////////////////////
    // Material definition
    Air = MaterialManager::getInstance()->GetMaterialFromLibrary("Air");
    Al = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
    MgO = MaterialManager::getInstance()->GetMaterialFromLibrary("MgO");
    Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    muMetal = MaterialManager::getInstance()->GetMaterialFromLibrary("mumetal");
    BoroSili_Glass = MaterialManager::getInstance()->GetMaterialFromLibrary("Borosillicate_Glass");
    ////////////////////////////////////////////////////////////
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* Dali2::BuildDali2Detector(int type){ // Single crystal
  ///////////////////////////////////////////////////
  // type = 0: SG, type = 1: SC, type = 2: DALI1
  ///////////////////////////////////////////////////
  if      (type == 0 && m_Dali2SGDetector) return m_Dali2SGDetector;
  else if (type == 1 && m_Dali2SCDetector) return m_Dali2SCDetector;
  else if (type == 2 && m_Dali1Detector)   return m_Dali1Detector;

  
  G4Box* sDali2Volume = new G4Box("sDali2Volume",
				  Dali2_NS::AlHouseDepth[type]/2.0,
				  Dali2_NS::AlHouseHight[type]/2.0,
				  Dali2_NS::AlHouseWidth[type]/2.0);
  G4LogicalVolume* lDali2Volume = new G4LogicalVolume(sDali2Volume, 
						      Air,"lDali2Volume",0,0,0);
  // Housing:
  G4Box* sDali2AlHouseOutSG = new G4Box("sDali2AlHouseOutSG",
					Dali2_NS::AlHouseDepth[type]/2.0,
					Dali2_NS::AlHouseHight[type]/2.0,
					Dali2_NS::AlHouseWidth[type]/2.0);
  G4Box* sDali2AlHouseInSG  =
    new G4Box("sDali2AlHouseInSG",
	      Dali2_NS::AlHouseDepth[type]/2.0 - Dali2_NS::AlHouseThick[0],
	      Dali2_NS::AlHouseHight[type]/2.0 - Dali2_NS::AlHouseThick[0],
	      Dali2_NS::AlHouseWidth[type]/2.0 - Dali2_NS::AlHouseThick[0]);
  G4SubtractionSolid* sDali2AlHouseSG =
    new G4SubtractionSolid("sDali2AlHouseSG",sDali2AlHouseOutSG,sDali2AlHouseInSG);
  G4LogicalVolume* lDali2AlHouseSG =
    new G4LogicalVolume(sDali2AlHouseSG, Al,"lDali2AlHouseSG",0,0,0);

  // The Saint-Gobain Coating:
  G4Box* sDali2MgOCoatOutSG = new G4Box("sDali2MgOCoatOutSG",
					Dali2_NS::MgOCoatDepth[type]/2.0,
					Dali2_NS::MgOCoatHight[type]/2.0,
					Dali2_NS::MgOCoatWidth[type]/2.0);
  G4Box* sDali2MgOCoatInSG  =
    new G4Box("sDali2MgOCoatInSG",
	      Dali2_NS::MgOCoatDepth[type]/2.0 - Dali2_NS::MgOCoatThick[0],
	      Dali2_NS::MgOCoatHight[type]/2.0 - Dali2_NS::MgOCoatThick[0],
	      Dali2_NS::MgOCoatWidth[type]/2.0 - Dali2_NS::MgOCoatThick[0]);
  G4SubtractionSolid* sDali2MgOCoatSG =
    new G4SubtractionSolid("sDali2MgOCoatSG",sDali2MgOCoatOutSG,sDali2MgOCoatInSG);
  G4LogicalVolume* lDali2MgOCoatSG =
    new G4LogicalVolume(sDali2MgOCoatSG, MgO,"lDali2MgOCoatSG",0,0,0);

  // The Saint-Gobain:
  G4Box* sDali2CrystalSG = new G4Box("sDali2CrystalSG",
				     Dali2_NS::CrystalDepth[type]/2.0,
				     Dali2_NS::CrystalHight[type]/2.0,
				     Dali2_NS::CrystalWidth[type]/2.0);
  G4LogicalVolume* lDali2CrystalSG =
    new G4LogicalVolume(sDali2CrystalSG,NaI_Tl,"lDali2CrystalSG",0,0,0);

  /////// PMTs objets 
  G4Tubs* sAlPMT = new G4Tubs("AlPMT",19.5*mm, 20.0*mm,
			      75.0*mm,0*deg,360*deg);
  G4LogicalVolume* lAlPMT = new G4LogicalVolume(sAlPMT, Al ,"lAlPMT",0,0,0);
    

  ////////////////////////////////////////////////////////////
  // Placement (Physical)
  G4ThreeVector positionnull = G4ThreeVector(0,0,0);
  G4RotationMatrix Rot3D;

  // PMT  part -
  Rot3D.set(0,0,0);
  Rot3D.rotateY(90.*degree);
  //  new G4PVPlacement(G4Transform3D(Rot3D,G4ThreeVector(15.*cm, 0., 0.)),
  //  		    lAlPMT ,"AlPMT",lDali2Volume,false,0);
    
  // Cryst Part -
  new G4PVPlacement(0, positionnull,lDali2MgOCoatSG,"MgO_Can",
		    lDali2Volume,false,0); 
  new G4PVPlacement(0, positionnull,lDali2CrystalSG,"CrystalNaI",
		    lDali2Volume,false,0); 
     
  G4VisAttributes* MgO_Color = new G4VisAttributes(G4Colour(1,1,0, .4));
  G4VisAttributes* Al_Color[3];
  Al_Color[0] = new G4VisAttributes(G4Colour(0.5,0.5,0.5, .3));
  Al_Color[1] = new G4VisAttributes(G4Colour(0.1,0.5,0.5, .3));
  Al_Color[2] = new G4VisAttributes(G4Colour(0.5,0.5,0.1, .3));
  G4VisAttributes* Crystal_Color = new G4VisAttributes(G4Colour(0, 1, 1));   
  //G4VisAttributes* mumetal_Color = new G4VisAttributes(G4Colour(0, 0.5, 1, .3));   

  lDali2MgOCoatSG->SetVisAttributes(MgO_Color);
  lDali2CrystalSG->SetVisAttributes(Crystal_Color);
  
  lDali2AlHouseSG->SetVisAttributes(Al_Color[type]);
  lAlPMT->SetVisAttributes(Al_Color[type]);

  lDali2CrystalSG->SetSensitiveDetector(m_DaliScorer);

  if      (type == 0) return (m_Dali2SGDetector = lDali2Volume);
  else if (type == 1) return (m_Dali2SCDetector = lDali2Volume);
  else if (type == 2) return (m_Dali1Detector = lDali2Volume);
  else return NULL;
} // end BuildSquareDetector

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Dali2::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Dali2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS", "Type"};
  vector<string> cartang = {"POS", "ANG", "Type"};
  vector<string> cyli = {"R","Alpha","Zeta","Type"};

  for (unsigned int i = 0 ; i < blocks.size() ; i++){
    if (blocks[i]->HasTokenList(cartang)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Dali2 " << i+1 <<  endl;
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      G4ThreeVector Ang = NPS::ConvertVector(blocks[i]->GetTVector3("ANG","deg"));
      int    Type = blocks[i]->GetInt("Type");
      AddDetector(Pos, Ang, Type); }

    else if (blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Dali2 " << i+1 <<  endl;
      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      int    Type = blocks[i]->GetInt("Type");
      AddDetector(Pos, Type); }


    else if(blocks[i]->HasTokenList(cyli)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Dali2 " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Alpha = blocks[i]->GetDouble("Alpha","deg");
      double Zeta = blocks[i]->GetDouble("Zeta","mm");
      int    Type = blocks[i]->GetInt("Type");
      AddDetector(R,Alpha,Zeta,Type);
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
void Dali2::ConstructDetector(G4LogicalVolume* world){    

  DefinitionMaterials();
  for (unsigned short i = 0 ; i < m_x.size() ; i++) {
    G4ThreeVector Det_pos = G4ThreeVector(m_x[i], m_y[i], m_z[i]) ;

    G4RotationMatrix* Rot = new G4RotationMatrix();
    if (m_psi.size() > i) {
      Rot->rotateX(m_psi[i]);
      Rot->rotateY(m_theta[i]);
      Rot->rotateZ(m_phi[i]);}
    else {
      Rot->rotateX(180*deg);
      Rot->rotateZ(atan(m_y[i]/m_x[i]));}
    
    new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
		      BuildDali2Detector(m_Type[i]),
		      "Dali2",world,false,i+1);}}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Dali2::InitializeRootOutput(){
  RootOutput *pAnalysis = RootOutput::getInstance();
  TTree *pTree = pAnalysis->GetTree();
  if(!pTree->FindBranch("Dali2")){
    pTree->Branch("Dali2", "TDali2Data", &m_Event) ;
  }
  pTree->SetBranchAddress("Dali2", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Dali2::ReadSensitive(const G4Event* ){
  m_Event->Clear();
  ///////////

  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_DaliScorer->GetPrimitive(0);

  //  cout << m_DaliScorer->GetNumberOfPrimitives()<<endl;
  unsigned int size = Scorer->GetMult();
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i); 
     double Energy = RandGauss::shoot(Scorer->GetEnergy(i),Dali2_NS::ResoEnergy*std::sqrt(Scorer->GetEnergy(i)));
    if(Energy>Dali2_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),Dali2_NS::ResoTime);
      //      int ArrayNbr = level[0];
      //      int DetectinsArrayNbr = level[0]+1;
      //      int DetectorNbr = (ArrayNbr-1)*3+DetectinsArrayNbr;
      //      int DetectorNbr = (ArrayNbr-1)*3+DetectinsArrayNbr;
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time);
      /* m_Event->SetParticleID(Scorer->GetParticleID(i)); */
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////   
void Dali2::InitializeScorers() { 
  // This check is necessary in case the geometry is reloaded
  bool already_exist = false;
  vector<int> NestingLevel;
  NestingLevel.push_back(1);
  //  NestingLevel.push_back(4);

  m_DaliScorer = CheckScorer("DaliScorer",already_exist) ;

  if(already_exist) //Necessary?
    return ;  //Necessary?
  // Otherwise the scorer is initialised
  //  vector<int> level; level.push_back(0);
  G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter", NestingLevel) ;
  G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;

  ////////////////////////////////////////////////////////////
  // filter
  std::vector<G4String> particles;
  particles.push_back("gamma");
  particles.push_back("e-");
  particles.push_back("e+");
  particles.push_back("proton");
  G4SDParticleFilter* sdPFilter = new G4SDParticleFilter("gamma",particles);
  Calorimeter->SetFilter(sdPFilter);
  ////////////////////////////////////////////////////////////
  
  //and register it to the multifunctionnal detector
  m_DaliScorer->RegisterPrimitive(Calorimeter);
  m_DaliScorer->RegisterPrimitive(Interaction);
  G4SDManager::GetSDMpointer()->AddNewDetector(m_DaliScorer) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Dali2::Construct(){
  return  (NPS::VDetector*) new Dali2();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
  class proxy_nps_Dali2{
    public:
      proxy_nps_Dali2(){
        NPS::DetectorFactory::getInstance()->AddToken("Dali2","Dali2");
        NPS::DetectorFactory::getInstance()->AddDetector("Dali2",Dali2::Construct);
      }
  };

  proxy_nps_Dali2 p_nps_Dali2;
}
