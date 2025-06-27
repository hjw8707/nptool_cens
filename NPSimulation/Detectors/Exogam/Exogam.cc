/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Goigoux  contact address: thomas.goigoux@cea.fr   *
 *                                                                           *
 * Creation Date  : july 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Exogam simulation                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <cmath>
#include <limits>
#include <sstream>
// G4 Geometry object
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Para.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4Plane3D.hh"
#include "G4IntersectionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

// G4 sensitive
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"

// G4 various object
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4PVPlacement.hh"
#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

// NPTool header
#include "CalorimeterScorers.hh"
#include "Exogam.hh"
#include "MaterialManager.hh"
#include "NPOptionManager.h"
#include "NPSDetectorFactory.hh"
#include "NPSHitsMap.hh"
#include "RootOutput.h"
#include "InteractionScorers.hh"
#include "TExogamGeo.h"
// CLHEP header
#include "CLHEP/Random/RandGauss.h"

using namespace std;
using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace Exogam_NS {
  // Energy and time Resolution
  const double EnergyThresholdInner6MV = 10 * keV;
  const double EnergyThresholdOuter = 50 * keV;
  // const double ResoTime = 4.5*ns ;  //not used
  const double ResoEnergyInner6MV = 2. * keV;
  const double ResoEnergyOuter = 10. * keV;
} // namespace Exogam_NS
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Exogam Specific Method
Exogam::Exogam() {
  m_Event = new TExogamCalData();
  m_ExogamScorer = 0;

  InitializeMaterials();

  HalfLengthCan = 7.35 * cm;
  TaperLengthCan = 4.325 * cm;
  distCollimatorToBGOSShield = 2.95 * cm;

  rm90.rotateZ(90. * deg);
  rm90m.rotateZ(-90. * deg);
  rm180.rotateZ(180. * deg);
  rm270.rotateZ(270. * deg);
}

Exogam::~Exogam() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void Exogam::ReadConfiguration(NPL::InputParser parser) {

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> coord = {"X", "Y", "Z", "ThetaX", "ThetaY", "ThetaZ"};
  vector<string> CoordCloverXYZ = {"Flange",
    "CrystalA_Seg1","CrystalA_Seg2","CrystalA_Seg3","CrystalA_Seg4",
    "CrystalB_Seg1","CrystalB_Seg2","CrystalB_Seg3","CrystalB_Seg4",
    "CrystalC_Seg1","CrystalC_Seg2","CrystalC_Seg3","CrystalC_Seg4",
    "CrystalD_Seg1","CrystalD_Seg2","CrystalD_Seg3","CrystalD_Seg4"};
  vector<string> sphe = {"R", "Theta", "Phi"};
  // Considering the flange12 as a reference for the stanard 12 exogam config
  vector<string> flange = {"Radius", "Flange"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(coord)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Exogam " << i + 1 << endl;
      double X = blocks[i]->GetDouble("X", "mm");
      double Y = blocks[i]->GetDouble("Y", "mm");
      double Z = blocks[i]->GetDouble("Z", "mm");
      double ThetaX = blocks[i]->GetDouble("ThetaX", "deg");
      double ThetaY = blocks[i]->GetDouble("ThetaY", "deg");
      double ThetaZ = blocks[i]->GetDouble("ThetaZ", "deg");
      AddDetector(X, Y, Z, ThetaX, ThetaY, ThetaZ);
    }
    else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Exogam " << i + 1 << endl;
      double R = blocks[i]->GetDouble("R", "mm");
      double Theta = blocks[i]->GetDouble("Theta", "deg");
      double Phi = blocks[i]->GetDouble("Phi", "deg");
      AddDetector(R, Theta, Phi);
    }
    else if (blocks[i]->HasTokenList(flange)) {
      double R = blocks[i]->GetDouble("Radius", "mm");
      int Flange = blocks[i]->GetInt("Flange");
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Added Exogam " << i + 1 << " flange " << Flange << " at " << R << " mm" << endl;
      AddDetector(R,Flange);
    }
    else if (blocks[i]->HasTokenList(CoordCloverXYZ)) {
      G4ThreeVector CrystalA_Seg1 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg1", "mm"));
      G4ThreeVector CrystalA_Seg2 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg2", "mm"));
      G4ThreeVector CrystalA_Seg3 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg3", "mm"));
      G4ThreeVector CrystalA_Seg4 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg4", "mm"));
      
      G4ThreeVector CrystalB_Seg1 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg1", "mm"));
      G4ThreeVector CrystalB_Seg2 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg2", "mm"));
      G4ThreeVector CrystalB_Seg3 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg3", "mm"));
      G4ThreeVector CrystalB_Seg4 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg4", "mm"));
      
      G4ThreeVector CrystalC_Seg1 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg1", "mm"));
      G4ThreeVector CrystalC_Seg2 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg2", "mm"));
      G4ThreeVector CrystalC_Seg3 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg3", "mm"));
      G4ThreeVector CrystalC_Seg4 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg4", "mm"));
      
      G4ThreeVector CrystalD_Seg1 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg1", "mm"));
      G4ThreeVector CrystalD_Seg2 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg2", "mm"));
      G4ThreeVector CrystalD_Seg3 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg3", "mm"));
      G4ThreeVector CrystalD_Seg4 = NPS::ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg4", "mm"));
      int Flange = blocks[i]->GetInt("Flange");
      AddDetector(Flange,
        CrystalA_Seg1,CrystalA_Seg2,CrystalA_Seg3,CrystalA_Seg4,
        CrystalB_Seg1,CrystalB_Seg2,CrystalB_Seg3,CrystalB_Seg4,
        CrystalC_Seg1,CrystalC_Seg2,CrystalC_Seg3,CrystalC_Seg4,
        CrystalD_Seg1,CrystalD_Seg2,CrystalD_Seg3,CrystalD_Seg4
      ); 
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Exogam::AddDetector(double X, double Y, double Z, double ThetaX, double ThetaY, double ThetaZ) {
  m_X.push_back(X);
  m_Y.push_back(Y);
  m_Z.push_back(Z);
  m_ThetaX.push_back(ThetaX);
  m_ThetaY.push_back(ThetaY);
  m_ThetaZ.push_back(ThetaZ);

  m_R.push_back(-1000);
  m_Theta.push_back(-1000);
  m_Phi.push_back(-1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Exogam::AddDetector(double R, double Theta, double Phi) {
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);

  m_X.push_back(-1000);
  m_Y.push_back(-1000);
  m_Z.push_back(-1000);
  m_ThetaX.push_back(-1000);
  m_ThetaY.push_back(-1000);
  m_ThetaZ.push_back(-1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Exogam::AddDetector(double R, int flange) {
  Distances[flange]=R;
}

void Exogam::AddDetector(unsigned int Flange,
        G4ThreeVector CrystalA_Seg1,G4ThreeVector CrystalA_Seg2,G4ThreeVector CrystalA_Seg3,G4ThreeVector CrystalA_Seg4,
        G4ThreeVector CrystalB_Seg1,G4ThreeVector CrystalB_Seg2,G4ThreeVector CrystalB_Seg3,G4ThreeVector CrystalB_Seg4,
        G4ThreeVector CrystalC_Seg1,G4ThreeVector CrystalC_Seg2,G4ThreeVector CrystalC_Seg3,G4ThreeVector CrystalC_Seg4,
        G4ThreeVector CrystalD_Seg1,G4ThreeVector CrystalD_Seg2,G4ThreeVector CrystalD_Seg3,G4ThreeVector CrystalD_Seg4
      )
      {
        Coordinates[Flange][0][0] = CrystalA_Seg1;
        Coordinates[Flange][0][1] = CrystalA_Seg2;
        Coordinates[Flange][0][2] = CrystalA_Seg3;
        Coordinates[Flange][0][3] = CrystalA_Seg4;
        
        Coordinates[Flange][1][0] = CrystalB_Seg1;
        Coordinates[Flange][1][1] = CrystalB_Seg2;
        Coordinates[Flange][1][2] = CrystalB_Seg3;
        Coordinates[Flange][1][3] = CrystalB_Seg4;
        
        Coordinates[Flange][2][0] = CrystalC_Seg1;
        Coordinates[Flange][2][1] = CrystalC_Seg2;
        Coordinates[Flange][2][2] = CrystalC_Seg3;
        Coordinates[Flange][2][3] = CrystalC_Seg4;
        
        Coordinates[Flange][3][0] = CrystalD_Seg1;
        Coordinates[Flange][3][1] = CrystalD_Seg2;
        Coordinates[Flange][3][2] = CrystalD_Seg3;
        Coordinates[Flange][3][3] = CrystalD_Seg4;
      }

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void Exogam::ConstructDetector(G4LogicalVolume* world) {
  // G4double distBGOSShieldToGeCan = 3.2*cm;	 	// distance from the front face of the
  //  BGO Side Shield to the front face of the
  //  Ge can (theory: 3.2*cm)
  // G4double distCollimatorToGeCan=6.15*cm;		// distance from front face of the collimator
  //  to the front face of the Ge can
  
  // In this loop, CloverNbr corresponds to the Flange NUmber of Exo Structure
  for(auto Flange: Distances){
    CloverNbr = Flange.first;
    G4LogicalVolume* logicSupClover = nullptr;
    InitClover(world, logicSupClover);
    PlaceClover(world, logicSupClover, Flange);
    BuildClover(world, logicSupClover);
  }
  // If an exogam clover is not in the struct, its numerotation starts at 20 to avoid overlap with Exogam Flange Numbering
  CloverNbr = 20;
  for (unsigned i = 0; i < m_X.size(); ++i) {
    G4LogicalVolume* logicSupCLover;
    InitClover(world, logicSupCLover);
    PlaceClover(world, logicSupCLover, i);
    BuildClover(world, logicSupCLover);
    CloverNbr++;   
    // BuildSideCatcher();
    // BuildBackCatcher();

    // BuildSideShield();
    // BuildCollimator();
  }
  
  for(auto Flange: Coordinates){
    CloverNbr = Flange.first;
    G4LogicalVolume* logicSupClover = nullptr;
    InitClover(world, logicSupClover);
    PlaceClover(world, logicSupClover, Flange);
    BuildClover(world, logicSupClover);
  }
  CloverNbr = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int Exogam::InitializeMaterials() {
  m_Vacuum = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
  m_Aluminum = MaterialManager::getInstance()->GetMaterialFromLibrary("Al");
  m_Copper = MaterialManager::getInstance()->GetMaterialFromLibrary("Cu");
  m_Germanium = MaterialManager::getInstance()->GetMaterialFromLibrary("Ge");

  m_BGO = new G4Material("BGO", 7.13 * g / cm3, 3, kStateSolid); // BGO does not exist in nptool !!
  m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("Bi"), 4);
  m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("Ge"), 3);
  m_BGO->AddElement(MaterialManager::getInstance()->GetElementFromLibrary("O"), 12);

  m_CsI = MaterialManager::getInstance()->GetMaterialFromLibrary("CsI");

  return 0;
}

void Exogam::InitClover(G4LogicalVolume* world, G4LogicalVolume*& logicSupClover){
  // enveloppe of the whole Clover (i.e. suppressed Clover) called 'SupClover' including:
  //  the cryostat, dewar, side shield, back catcher, collimator
  G4double dzEnv = 40.472 * cm;
  G4double dx1Env = 3.17 * cm;
  G4double dy1Env = 3.17 * cm;
  G4double dx2Env = 2. * dzEnv * tan(22.5 * deg) + dx1Env;
  G4double dy2Env = 2. * dzEnv * tan(22.5 * deg) + dy1Env;

  G4Trd* solidSupClover = new G4Trd("SupClover", dx1Env, dx2Env, dy1Env, dy2Env, dzEnv);
  logicSupClover = new G4LogicalVolume(solidSupClover, m_Vacuum, "SupClover");
}

void Exogam::PlaceClover(G4LogicalVolume* world, G4LogicalVolume*& logicSupClover, std::pair<const unsigned int, double> Flange){
  G4double dzEnv = 40.472 * cm;
  G4double Offset = dzEnv; //-distCollimatorToGeCan;
    
  auto ExoGeo = new TExogamStructure(Distances);
  Vector3D NormalExo(0.,1.,0.);
  Vector3D SideExo1(1.,0.,0.);
  Vector3D SideExo2(0.,0.,1.);
  ExoGeo->SetCloverNbr(Flange.first);
  ExoGeo->setclover(&NormalExo);
  ExoGeo->setclover(&SideExo1);
  ExoGeo->setclover(&SideExo2);
  G4ThreeVector G4NormalExo(NormalExo.X(),NormalExo.Y(),NormalExo.Z());
  G4ThreeVector G4SideExo1(SideExo1.X(),SideExo1.Y(),SideExo1.Z());
  G4ThreeVector G4SideExo2(SideExo2.X(),SideExo2.Y(),SideExo2.Z());

  G4ThreeVector ExoPos(NormalExo.X()*(Flange.second+Offset), NormalExo.Y()*(Flange.second + Offset), NormalExo.Z()*(Flange.second + Offset));
  G4RotationMatrix ExoRot(-G4SideExo1,G4SideExo2, -G4NormalExo);
  std::string clover_name = "Clover_" + std::to_string(Flange.first);
  new G4PVPlacement(G4Transform3D(ExoRot,ExoPos), logicSupClover, clover_name, world, false, CloverNbr,
                  false); // this void overlaps the whole setup

}

void Exogam::PlaceClover(G4LogicalVolume* world, G4LogicalVolume*& logicSupClover, int i_clo){
  G4double dzEnv = 40.472 * cm;
  G4double Offset = dzEnv; //-distCollimatorToGeCan;
  
  if (m_X[i_clo] >= 0) {
    G4RotationMatrix rm;
    rm.rotateX(m_ThetaX[i_clo] / rad).rotateY(m_ThetaY[i_clo] / rad).rotateZ(m_ThetaZ[i_clo] / rad);
    new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(m_X[i_clo] * mm, m_Y[i_clo] * mm, m_Z[i_clo] * mm + Offset)),
                      logicSupClover, "Clover", world, false, CloverNbr, false); // this void overlaps the whole setup
  }
  else if (m_R[i_clo] >= 0 && m_Theta[i_clo] >= 0) {

    G4RotationMatrix* MMrot = NULL;
    G4ThreeVector MMpos = G4ThreeVector(0, 0, 0);
    G4ThreeVector MMu = G4ThreeVector(0, 0, 0);
    G4ThreeVector MMv = G4ThreeVector(0, 0, 0);
    G4ThreeVector MMw = G4ThreeVector(0, 0, 0);
    G4ThreeVector MMCenter = G4ThreeVector(0, 0, 0);
    G4double Theta = m_Theta[i_clo];
    G4double Phi = m_Phi[i_clo];

    // (u,v,w) unitary vector associated to telescope referencial
    // (u,v) // to silicon plan
    // w perpendicular to (u,v) plan and pointing ThirdStage
    // Phi is angle between X axis and projection in (X,Y) plan
    // Theta is angle between  position vector and z axis
    G4double wX = (m_R[i_clo] + Offset) * sin(Theta / rad) * cos(Phi / rad);
    G4double wY = (m_R[i_clo] + Offset) * sin(Theta / rad) * sin(Phi / rad);
    G4double wZ = (m_R[i_clo] + Offset) * cos(Theta / rad);
    MMw = G4ThreeVector(wX, wY, wZ);

    // vector corresponding to the center of the module
    G4ThreeVector CT = MMw;

    // vector parallel to one axis of silicon plane
    G4double ii = cos(Theta / rad) * cos(Phi / rad);
    G4double jj = cos(Theta / rad) * sin(Phi / rad);
    G4double kk = -sin(Theta / rad);
    G4ThreeVector Y = G4ThreeVector(ii, jj, kk);

    MMw = MMw.unit();
    MMu = MMw.cross(Y);
    MMv = MMw.cross(MMu);
    MMv = MMv.unit();
    MMu = MMu.unit();

    // Passage Matrix from Lab Referential to Telescope Referential
    // MUST2
    MMrot = new G4RotationMatrix(MMu, MMv, MMw);
    // Telescope is rotate of Beta angle around MMv axis.
    // MMrot->rotate(m_beta_u[i_clo], MMu);
    // MMrot->rotate(m_beta_v[i_clo], MMv);
    // MMrot->rotate(m_beta_w[i_clo], MMw);
    // translation to place Telescope
    double Length = 0 * cm;
    MMpos = MMw * Length * 0.5 + CT;
    std::cout << i_clo << std::endl;
    std::string clover_name = "Clover_" + std::to_string(CloverNbr);
    new G4PVPlacement(G4Transform3D(*MMrot, MMpos), logicSupClover, clover_name, world, false, CloverNbr,
                      false); // this void overlaps the whole setup
  }
}
void Exogam::PlaceClover(G4LogicalVolume* world, G4LogicalVolume*& logicSupClover, std::pair<const unsigned int, std::map<unsigned int, std::map<unsigned int, G4ThreeVector>>> Flange){
  G4double dzEnv = 40.472 * cm;
  G4double Offset = dzEnv; //-distCollimatorToGeCan;

  G4RotationMatrix* Exorot = NULL;
  G4ThreeVector Exopos = G4ThreeVector(0, 0, 0);
  G4ThreeVector ExoTest1 = G4ThreeVector(0, 0, 0);
  G4ThreeVector ExoTest2 = G4ThreeVector(0, 0, 0);
  G4ThreeVector ExoTest3 = G4ThreeVector(0, 0, 0);
  G4ThreeVector ExoTest4 = G4ThreeVector(0, 0, 0);
  G4ThreeVector Exou = G4ThreeVector(0, 0, 0);
  G4ThreeVector Exov = G4ThreeVector(0, 0, 0);
  G4ThreeVector Exow = G4ThreeVector(0, 0, 0);
  G4ThreeVector ExoCenter = G4ThreeVector(0, 0, 0);

  ExoTest1 = 0.25*(Flange.second[0][0] + Flange.second[1][0] + Flange.second[2][0] + Flange.second[3][0]);
  ExoTest2 = 0.25*(Flange.second[0][2] + Flange.second[1][2] + Flange.second[2][2] + Flange.second[3][2]);
  ExoTest3 = 0.25*(Flange.second[0][1] + Flange.second[1][3] + Flange.second[2][1] + Flange.second[3][3]);
  ExoTest4 = 0.25*(Flange.second[0][3] + Flange.second[1][1] + Flange.second[2][3] + Flange.second[3][1]);
  bool TestGeo = ((ExoTest1 - ExoTest2).mag() < 1*mm) &&((ExoTest1 - ExoTest3).mag() < 1*mm) &&((ExoTest1 - ExoTest4).mag() < 1*mm);
  if(TestGeo){
    Exou = Flange.second[1][0] - Flange.second[0][0];
    Exou.unit();
    Exov = Flange.second[3][0] - Flange.second[0][0];
    Exov.unit();
    Exow = Exov.cross(Exou);
    ExoCenter = ExoTest1;
    Exorot = new G4RotationMatrix(-Exou, Exov, -Exow);
    double Length = 0 * cm;
    Exopos = Exow * Length * 0.5 + ExoCenter.unit()*Offset + ExoCenter;
    std::string clover_name = "Clover_" + std::to_string(CloverNbr);
    new G4PVPlacement(G4Transform3D(*Exorot, Exopos), logicSupClover, clover_name, world, false, CloverNbr,
                      false); // this void overlaps the whole setup
  }
  else
    std::cout << "WARNING: EXOGAM clover number " << Flange.first << " not constructed because of cartesian geometry problem, check that your coordinates are right" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void Exogam::BuildClover(G4LogicalVolume* world, G4LogicalVolume*& logicSupClover) {
  G4double dzEnv = 40.472 * cm;
  G4double Offset = dzEnv; //-distCollimatorToGeCan;
  
  // The Cryostat
  ////////////////
  // The Aluminum Clover can ( "CloverCan" )...
  //
  G4double PhiStartCan = 45. * deg;
  G4double PhiTotCan = 360. * deg;

  G4double zPlaneCan[3];
  G4double rInnerCan[3];

  G4double zPlaneVac[3];
  G4double rInnerVac[3];
  G4double rOuterVac[3];

  zPlaneCan[0] = -HalfLengthCan;
  zPlaneCan[1] = -HalfLengthCan + TaperLengthCan;
  zPlaneCan[2] = HalfLengthCan;

  G4double rOuterCan[3]; // used to build the shield
  rOuterCan[0] = 4.4085 * cm;
  rOuterCan[1] = 6.2 * cm;
  rOuterCan[2] = 6.2 * cm;

  rInnerCan[0] = rInnerCan[1] = rInnerCan[2] = 0.1 * cm;

  G4Polyhedra* solidCloverCan =
      new G4Polyhedra("CloverCan", PhiStartCan, PhiTotCan, 4, 3, zPlaneCan, rInnerCan, rOuterCan);

  G4LogicalVolume* logicCloverCan = new G4LogicalVolume(solidCloverCan, m_Aluminum, "CloverCan");

  // The position of the Clover can in the SupClover:
  G4ThreeVector posClover(0. * cm, 0. * cm, -Offset + HalfLengthCan + 0.001 * mm); //+0.001mm to avoid roundoff errors

  new G4PVPlacement(0, posClover, logicCloverCan, "CloverCan", logicSupClover, false,CloverNbr,
                    true); // There is an overlap with vacuum SupClover

  // The vacuum clover ( "Vac" ) ...
  //
  G4double HalfLengthVac = 7.175 * cm;
  G4double TaperLengthVac = 4.0842 * cm;

  zPlaneVac[0] = -HalfLengthVac;
  zPlaneVac[1] = -HalfLengthVac + TaperLengthVac;
  zPlaneVac[2] = HalfLengthVac;
  rOuterVac[0] = 4.3083 * cm;
  rOuterVac[1] = 6.0 * cm;
  rOuterVac[2] = 6.0 * cm;

  rInnerVac[0] = rInnerVac[1] = rInnerVac[2] = 0. * cm;

  G4Polyhedra* solidVac = new G4Polyhedra("Vac", PhiStartCan, PhiTotCan, 4, 3, zPlaneVac, rInnerVac, rOuterVac);
  G4LogicalVolume* logicVac = new G4LogicalVolume(solidVac, m_Vacuum, "Vac");

  G4ThreeVector positionVac = G4ThreeVector(0. * cm, 0. * cm, -0.25 * mm);
  new G4PVPlacement(0, positionVac, logicVac, "Vac", logicCloverCan, false, CloverNbr, true);

  //
  // The enveloppe of the cold finger from the back side of the can to the Dewar
  //

  G4double zPlaneEnvColdFinger[6];
  G4double rInnerEnvColdFinger[6];
  G4double rOuterEnvColdFinger[6];

  G4double PhiStart = 0. * deg;
  G4double PhiTot = 360. * deg;
  G4double EnvColdFingerHalfLength = 7.24 * cm;

  zPlaneEnvColdFinger[0] = -EnvColdFingerHalfLength;
  zPlaneEnvColdFinger[1] = -EnvColdFingerHalfLength + 4.1 * cm;
  zPlaneEnvColdFinger[2] = -EnvColdFingerHalfLength + 4.1 * cm;
  zPlaneEnvColdFinger[3] = -EnvColdFingerHalfLength + 4.9 * cm;
  zPlaneEnvColdFinger[4] = -EnvColdFingerHalfLength + 4.9 * cm;
  zPlaneEnvColdFinger[5] = EnvColdFingerHalfLength;

  rInnerEnvColdFinger[0] = rInnerEnvColdFinger[1] = rInnerEnvColdFinger[2] = 0. * cm;
  rInnerEnvColdFinger[3] = rInnerEnvColdFinger[4] = rInnerEnvColdFinger[5] = 0. * cm;

  rOuterEnvColdFinger[0] = 2.225 * cm;
  rOuterEnvColdFinger[1] = 2.225 * cm;
  rOuterEnvColdFinger[2] = 3.1 * cm;
  rOuterEnvColdFinger[3] = 3.1 * cm;
  rOuterEnvColdFinger[4] = 2.225 * cm;
  rOuterEnvColdFinger[5] = 2.225 * cm;

  G4Polycone* solidEnvColdFinger = new G4Polycone("EnvColdFinger", PhiStart, PhiTot, 6, zPlaneEnvColdFinger,
                                                  rInnerEnvColdFinger, rOuterEnvColdFinger);

  G4LogicalVolume* logicEnvColdFinger = new G4LogicalVolume(solidEnvColdFinger, m_Aluminum, "EnvColdFinger");

  G4ThreeVector posEnvColdFinger =
      G4ThreeVector(0. * cm, 0. * cm, -Offset + 2. * HalfLengthCan + EnvColdFingerHalfLength + 0.005 * mm);

  new G4PVPlacement(0, posEnvColdFinger, logicEnvColdFinger, "EnvColdFinger", logicSupClover, false, CloverNbr, true);

  // Its internal vacuum...
  G4double minRadiusIntEnvColdFinger = 0. * cm;
  G4double maxRadiusIntEnvColdFinger = 2.025 * cm;
  G4double HalfLengthIntEnvColdFinger = 7.24 * cm;
  G4double startPhiIntEnvColdFinger = 0. * deg;
  G4double deltaPhiIntEnvColdFinger = 360. * deg;

  G4Tubs* solidIntEnvColdFinger =
      new G4Tubs("IntDewar", minRadiusIntEnvColdFinger, maxRadiusIntEnvColdFinger, HalfLengthIntEnvColdFinger,
                 startPhiIntEnvColdFinger, deltaPhiIntEnvColdFinger);

  G4LogicalVolume* logicIntEnvColdFinger = new G4LogicalVolume(solidIntEnvColdFinger, m_Vacuum, "IntEnvColdFinger");

  // and its position in the cold finger enveloppe.
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), logicIntEnvColdFinger, "IntEnvColdFinger", logicEnvColdFinger, false,
                    CloverNbr, true);

  // The cold finger and the associated plate
  //
  G4double xHalfLengthCFPlate = 5.04 * cm;
  G4double yHalfLengthCFPlate = 5.04 * cm;
  G4double zHalfLengthCFPlate = 1. * mm;

  G4Box* solidCFPlate = new G4Box("CFPlate", xHalfLengthCFPlate, yHalfLengthCFPlate, zHalfLengthCFPlate);

  G4LogicalVolume* logicCFPlate = new G4LogicalVolume(solidCFPlate, m_Copper, "CFPlate");

  G4ThreeVector posCFPlate(0. * cm, 0. * cm, -HalfLengthVac + 9.65 * cm); // 0.55(d(IntCan-Ge)
                                                                          // +9.(Ge length)+0.1(half length plate)
  new G4PVPlacement(0, posCFPlate, logicCFPlate, "CFPlate", logicVac, false, CloverNbr, true);

  // The cold finger (internal part)
  //
  G4double minRadiusIntCF = 0. * cm;
  G4double maxRadiusIntCF = 1.5 * cm;
  G4double HalfLengthIntCF = 2.30 * cm;
  G4double startPhiIntCF = 0. * deg;
  G4double deltaPhiIntCF = 360. * deg;

  G4Tubs* solidIntCF =
      new G4Tubs("IntCF", minRadiusIntCF, maxRadiusIntCF, HalfLengthIntCF, startPhiIntCF, deltaPhiIntCF);

  G4LogicalVolume* logicIntCF = new G4LogicalVolume(solidIntCF, m_Copper, "IntCF");

  // its position vs CloverCan...
  G4ThreeVector posIntCF(0. * cm, 0. * cm, 4.875 * cm); // -7.175 (halflengthcan internal)
                                                        // +0.55 (ext Can - Ge)
                                                        // +9.0 (Ge length)
                                                        // +0.2 (CF plate)
                                                        // +2.3 (IntCF length)

  new G4PVPlacement(0, posIntCF, logicIntCF, "IntCF", logicVac, false, CloverNbr, true);

  // The cold finger (external part)
  //
  G4double minRadiusExtCF = 0. * cm;
  G4double maxRadiusExtCF = 2.0 * cm;
  G4double HalfLengthExtCF = 7.2 * cm;
  G4double startPhiExtCF = 0. * deg;
  G4double deltaPhiExtCF = 360. * deg;

  G4Tubs* solidExtCF =
      new G4Tubs("IntCF", minRadiusExtCF, maxRadiusExtCF, HalfLengthExtCF, startPhiExtCF, deltaPhiExtCF);

  G4LogicalVolume* logicExtCF = new G4LogicalVolume(solidExtCF, m_Copper, "ExtCF");

  // its position vs EnvColdFinger...
  G4ThreeVector posExtCF(0. * cm, 0. * cm, 0. * cm);
  new G4PVPlacement(0, posExtCF, logicExtCF, "ExtCF", logicIntEnvColdFinger, false, CloverNbr, true);

  // The Dewar
  //
  G4double minRadiusDewar = 0. * cm;
  G4double maxRadiusDewar = 10.9 * cm;
  G4double HalfLengthDewar = 15.2 * cm;
  G4double startPhiDewar = 0. * deg;
  G4double deltaPhiDewar = 360. * deg;

  G4Tubs* solidDewar =
      new G4Tubs("Dewar", minRadiusDewar, maxRadiusDewar, HalfLengthDewar, startPhiDewar, deltaPhiDewar);

  G4LogicalVolume* logicDewar = new G4LogicalVolume(solidDewar, m_Aluminum, "Dewar");

  G4double distFrontToMidDewar = -Offset + 2. * (HalfLengthCan + EnvColdFingerHalfLength) + HalfLengthDewar + 0.01 * mm;
  //+0.01mm to avoid roundoff errors

  G4ThreeVector posDewar = G4ThreeVector(0. * cm, 0. * cm, distFrontToMidDewar);
  new G4PVPlacement(0, posDewar, logicDewar, "Dewar", logicSupClover, false, CloverNbr, true);

  /////////////////////////////////////////
  //  Construction of the active Ge volume:
  /////////////////////////////////////////
  //  A: Ge diode built from cuts subtracted from a cylinder (the "GeDiode")
  //
  G4double minRadiusGeDiode = 0. * cm;
  G4double maxRadiusGeDiode = 3.0 * cm;
  G4double HalfLengthGeDiode = 4.5 * cm;
  G4double startPhiGeDiode = 0. * deg;
  G4double deltaPhiGeDiode = 360. * deg;

  G4Tubs* solidGeDiode =
      new G4Tubs("GeDiode", minRadiusGeDiode, maxRadiusGeDiode, HalfLengthGeDiode, startPhiGeDiode, deltaPhiGeDiode);
  //
  // External Tapered volume all along the diode ( "Cut1&2" )
  //
  //
  // Cut 1 :
  //
  G4double dummy = acos(2.9 / 3.0);
  G4double xHalfLengthCut1 = 0.5 * mm;
  G4double yHalfLengthCut1 = 2.9 * tan(dummy) * cm;
  G4double zHalfLengthCut1 = 4.55 * cm;

  G4Box* solidCut1 = new G4Box("Cut1", xHalfLengthCut1, yHalfLengthCut1, zHalfLengthCut1);

  //
  //... and its position vs GeDiode
  //

  G4ThreeVector transCut1(2.95 * cm, 0. * cm, 0. * cm);
  G4SubtractionSolid* solidGeMinusCut1 = new G4SubtractionSolid("GeMinusCut1", solidGeDiode, solidCut1, 0, transCut1);

  G4ThreeVector transCut2(0., 2.95 * cm, 0.);
  G4Transform3D positionCut2(rm90, transCut2);

  G4SubtractionSolid* solidGeMinusCut12 =
      new G4SubtractionSolid("GeMinusCut12", solidGeMinusCut1, solidCut1, positionCut2);
  //
  // External Tapered volume at the front face ( "Cut3&4" )

  G4double cosTap = cos(22.5 * deg);
  G4double sinTap = sin(22.5 * deg);
  G4double tanTap = tan(22.5 * deg);

  G4double xHalfLengthCut3 = 3.0 * cm;
  G4double yHalfLengthCut3 = 1.5 * cm * sinTap;
  G4double zHalfLengthCut3 = 1.5 * cm / cosTap;

  G4Box* solidCut3 = new G4Box("Cut3", xHalfLengthCut3, yHalfLengthCut3, zHalfLengthCut3 + 0.5 * cm);

  G4double yCut3 = 2.9 * cm - 1.5 * cm * tanTap + yHalfLengthCut3 * cosTap;

  G4double temp = zHalfLengthCut3 * cosTap - yHalfLengthCut3 * sinTap;
  G4double zCut3 = -HalfLengthGeDiode + temp;

  G4RotationMatrix rmCut3;
  rmCut3.rotateX(-22.5 * deg);

  G4ThreeVector transCut3(0., yCut3, zCut3);
  G4Transform3D positionCut3(rmCut3, transCut3);

  G4SubtractionSolid* solidGeMinusCut123 =
      new G4SubtractionSolid("GeMinusCut123", solidGeMinusCut12, solidCut3, positionCut3);

  G4Box* solidCut4 = new G4Box("Cut4", yHalfLengthCut3, xHalfLengthCut3, zHalfLengthCut3);

  G4RotationMatrix rmCut4;
  rmCut4.rotateY(22.5 * deg);

  G4ThreeVector transCut4(yCut3, 0., zCut3);
  G4Transform3D positionCut4(rmCut4, transCut4);

  G4SubtractionSolid* solidGeMinusCut1234 =
      new G4SubtractionSolid("GeMinusCut1234", solidGeMinusCut123, solidCut4, positionCut4);

  dummy = acos(2.45 / 3.0);
  G4double xHalfLengthCut5 = 5.5 * mm;
  G4double yHalfLengthCut5 = 2.45 * tan(dummy) * cm;
  G4double zHalfLengthCut5 = 4.55 * cm;

  G4Box* solidCut5 = new G4Box("Cut5", xHalfLengthCut5, yHalfLengthCut5, zHalfLengthCut5);

  G4ThreeVector transCut5(-3.0 * cm, 0. * cm, 0. * cm);

  G4SubtractionSolid* solidGeMinusCut12345 =
      new G4SubtractionSolid("GeMinusCut12345", solidGeMinusCut1234, solidCut5, 0, transCut5);

  G4ThreeVector transCut6(0., -3.0 * cm, 0.);
  G4Transform3D positionCut6(rm90, transCut6);

  G4SubtractionSolid* solidGe = new G4SubtractionSolid("Ge", solidGeMinusCut12345, solidCut5, positionCut6);

  // Now the individual diode is built; create logical volumes for each of
  // the four individual diodes A, B, C and D:

  // Separating each crystal in 4 outers
  double HalfOuterBoxWidth1 = 15*mm;
  double HalfOuterBoxWidth2 = 14.5*mm;
  double HalfOuterBoxLength1 = 20*mm;
  double HalfOuterBoxLength2 = 30*mm;
  
  double CenterOuterFront = 20.54*mm;
  double CenterOuterMid = 24.5*mm;
  double DistFrontMid = 30*mm;
  
  /////////////////////// Outer 1 /////////////////////////////
  G4Box* Outer1Cut1 = new G4Box("Outer1Cut1",HalfOuterBoxWidth1, HalfOuterBoxWidth1, HalfOuterBoxLength1);
  G4Box* Outer1Cut2 = new G4Box("Outer1Cut2",HalfOuterBoxWidth2, HalfOuterBoxWidth2, HalfOuterBoxLength2);
  G4Box* Outer1Cut3 = new G4Box("Outer1Cut3",30*mm,30*mm, 15*mm);
  double RotationAngle = std::atan((CenterOuterMid-CenterOuterFront)/DistFrontMid);
 



  G4RotationMatrix* RotationOuter1Cut = new G4RotationMatrix();
  RotationOuter1Cut->rotateX(-RotationAngle);
  // G4ThreeVector Xnew(1.,0.,0.);
  // G4ThreeVector Ynew(0.,1.,0.);
  // G4ThreeVector Znew(0.,0.,1.);
  // Xnew.rotateX(RotationAngle);
  // Ynew.rotateX(RotationAngle);
  // Znew.rotateX(RotationAngle);
  // RotationOuter1Cut->rotateAxes(Xnew,Ynew,Znew);
  RotationOuter1Cut->rotateY(RotationAngle);
  RotationOuter1Cut->rotateZ(-0.009);
  G4RotationMatrix* NoRot = new G4RotationMatrix();
  
  G4Transform3D Outer1trans1(*RotationOuter1Cut, G4ThreeVector(13.2*mm, 13.2*mm,-HalfLengthGeDiode+DistFrontMid/2));
  G4Transform3D Outer1trans2(*NoRot, G4ThreeVector(HalfOuterBoxWidth2, HalfOuterBoxWidth2, -HalfLengthGeDiode+HalfOuterBoxLength2+DistFrontMid));
  G4Transform3D Outer1trans3(*NoRot, G4ThreeVector(0*mm, 0*mm, -HalfLengthGeDiode+15*mm));
  G4IntersectionSolid* Outer1part1 = new G4IntersectionSolid("Outer1Part1", solidGe, Outer1Cut1, Outer1trans1);
  Outer1part1 = new G4IntersectionSolid("Outer1Part1", Outer1part1, Outer1Cut3, Outer1trans3);
  G4IntersectionSolid* Outer1part2 = new G4IntersectionSolid("OUter1Part2", solidGe, Outer1Cut2,Outer1trans2);
  
  G4UnionSolid* Outer1union12 = new G4UnionSolid("Outer1Union12",Outer1part1,Outer1part2);
  
  
  
  ////////////////////////// Outer 2 ///////////////////////////////
  
  G4Transform3D Outer2trans1(*RotationOuter1Cut, G4ThreeVector(13.2*mm, -17.1*mm,-HalfLengthGeDiode+DistFrontMid/2));
  G4Transform3D Outer2trans2(*NoRot, G4ThreeVector(HalfOuterBoxWidth2, -HalfOuterBoxWidth2, -HalfLengthGeDiode+HalfOuterBoxLength2+DistFrontMid));
  G4Transform3D Outer2trans3(*NoRot, G4ThreeVector(0*mm, 0*mm, -HalfLengthGeDiode+15*mm));
  G4IntersectionSolid* Outer2part1 = new G4IntersectionSolid("Outer2Part1", solidGe, Outer1Cut1, Outer2trans1);
  Outer2part1 = new G4IntersectionSolid("Outer2Part1", Outer2part1, Outer1Cut3, Outer2trans3);
  G4IntersectionSolid* Outer2part2 = new G4IntersectionSolid("Outer2Part2", solidGe, Outer1Cut2,Outer2trans2);
  
  G4UnionSolid* Outer2union12 = new G4UnionSolid("Outer2Union12",Outer2part1,Outer2part2);
  
  
  ////////////////////////// Outer 3 ///////////////////////////////
  
  G4Transform3D Outer3trans1(*RotationOuter1Cut, G4ThreeVector(-17.1*mm, -17.1*mm,-HalfLengthGeDiode+DistFrontMid/2));
  G4Transform3D Outer3trans2(*NoRot, G4ThreeVector(-HalfOuterBoxWidth2, -HalfOuterBoxWidth2, -HalfLengthGeDiode+HalfOuterBoxLength2+DistFrontMid));
  G4Transform3D Outer3trans3(*NoRot, G4ThreeVector(0*mm, 0*mm, -HalfLengthGeDiode+15*mm));
  G4IntersectionSolid* Outer3part1 = new G4IntersectionSolid("Outer3Part1", solidGe, Outer1Cut1, Outer3trans1);
  Outer3part1 = new G4IntersectionSolid("Outer3Part1", Outer3part1, Outer1Cut3, Outer3trans3);
  G4IntersectionSolid* Outer3part2 = new G4IntersectionSolid("Outer3Part2", solidGe, Outer1Cut2,Outer3trans2);
  
  G4UnionSolid* Outer3union12 = new G4UnionSolid("Outer3Union12",Outer3part1,Outer3part2);
  
  
  ////////////////////////// Outer 4 ///////////////////////////////
  
  G4Transform3D Outer4trans1(*RotationOuter1Cut, G4ThreeVector(-17.1*mm, 13.2*mm,-HalfLengthGeDiode+DistFrontMid/2));
  G4Transform3D Outer4trans2(*NoRot, G4ThreeVector(-HalfOuterBoxWidth2, HalfOuterBoxWidth2, -HalfLengthGeDiode+HalfOuterBoxLength2+DistFrontMid));
  G4Transform3D Outer4trans3(*NoRot, G4ThreeVector(0*mm, 0*mm, -HalfLengthGeDiode+15*mm));
  G4IntersectionSolid* Outer4part1 = new G4IntersectionSolid("Outer4Part1", solidGe, Outer1Cut1, Outer4trans1);
  Outer4part1 = new G4IntersectionSolid("Outer4Part1", Outer4part1, Outer1Cut3, Outer4trans3);
  G4IntersectionSolid* Outer4part2 = new G4IntersectionSolid("Outer4Part2", solidGe, Outer1Cut2,Outer4trans2);
  
  G4UnionSolid* Outer4union12 = new G4UnionSolid("Outer4Union12",Outer4part1,Outer4part2);
  
  
  std::string GeA_Outer1 = "GeA_Outer1_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterA1 = new G4LogicalVolume(Outer1union12, m_Germanium, GeA_Outer1);  
  std::string GeA_Outer2 = "GeA_Outer2_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterA2 = new G4LogicalVolume(Outer2union12, m_Germanium, GeA_Outer2);
  std::string GeA_Outer3 = "GeA_Outer3_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterA3 = new G4LogicalVolume(Outer3union12, m_Germanium, GeA_Outer3);
  std::string GeA_Outer4 = "GeA_Outer4_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterA4 = new G4LogicalVolume(Outer4union12, m_Germanium, GeA_Outer4);
  
  std::string GeB_Outer1 = "GeB_Outer1_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterB1 = new G4LogicalVolume(Outer1union12, m_Germanium, GeB_Outer1);  
  std::string GeB_Outer2 = "GeB_Outer2_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterB2 = new G4LogicalVolume(Outer2union12, m_Germanium, GeB_Outer2);
  std::string GeB_Outer3 = "GeB_Outer3_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterB3 = new G4LogicalVolume(Outer3union12, m_Germanium, GeB_Outer3);
  std::string GeB_Outer4 = "GeB_Outer4_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterB4 = new G4LogicalVolume(Outer4union12, m_Germanium, GeB_Outer4);
  
  std::string GeC_Outer1 = "GeC_Outer1_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterC1 = new G4LogicalVolume(Outer1union12, m_Germanium, GeC_Outer1);  
  std::string GeC_Outer2 = "GeC_Outer2_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterC2 = new G4LogicalVolume(Outer2union12, m_Germanium, GeC_Outer2);
  std::string GeC_Outer3 = "GeC_Outer3_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterC3 = new G4LogicalVolume(Outer3union12, m_Germanium, GeC_Outer3);
  std::string GeC_Outer4 = "GeC_Outer4_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterC4 = new G4LogicalVolume(Outer4union12, m_Germanium, GeC_Outer4);
  
  std::string GeD_Outer1 = "GeD_Outer1_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterD1 = new G4LogicalVolume(Outer1union12, m_Germanium, GeD_Outer1);  
  std::string GeD_Outer2 = "GeD_Outer2_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterD2 = new G4LogicalVolume(Outer2union12, m_Germanium, GeD_Outer2);
  std::string GeD_Outer3 = "GeD_Outer3_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterD3 = new G4LogicalVolume(Outer3union12, m_Germanium, GeD_Outer3);
  std::string GeD_Outer4 = "GeD_Outer4_" + std::to_string(CloverNbr);
  G4LogicalVolume* logicOuterD4 = new G4LogicalVolume(Outer4union12, m_Germanium, GeD_Outer4);

  logicOuterA1->SetSensitiveDetector(m_ExogamScorer);
  logicOuterA2->SetSensitiveDetector(m_ExogamScorer);
  logicOuterA3->SetSensitiveDetector(m_ExogamScorer);
  logicOuterA4->SetSensitiveDetector(m_ExogamScorer);
  
  logicOuterB1->SetSensitiveDetector(m_ExogamScorer);
  logicOuterB2->SetSensitiveDetector(m_ExogamScorer);
  logicOuterB3->SetSensitiveDetector(m_ExogamScorer);
  logicOuterB4->SetSensitiveDetector(m_ExogamScorer);
  
  logicOuterC1->SetSensitiveDetector(m_ExogamScorer);
  logicOuterC2->SetSensitiveDetector(m_ExogamScorer);
  logicOuterC3->SetSensitiveDetector(m_ExogamScorer);
  logicOuterC4->SetSensitiveDetector(m_ExogamScorer);
  
  logicOuterD1->SetSensitiveDetector(m_ExogamScorer);
  logicOuterD2->SetSensitiveDetector(m_ExogamScorer);
  logicOuterD3->SetSensitiveDetector(m_ExogamScorer);
  logicOuterD4->SetSensitiveDetector(m_ExogamScorer);
  

  // positioning the tapered partial diodes (A to D)
  // into the real vacuum of the can
  G4double HalfDistanceBetweenDiodes = 0.5 * mm;

  G4double xDumVac = 2.45 * cm + HalfDistanceBetweenDiodes;
  G4double yDumVac = 2.45 * cm + HalfDistanceBetweenDiodes;
  G4double zDumVac = -HalfLengthVac + 5.05 * cm; // 5.05 = 0.55 d(int can to Ge) +4.5(half length Ge)

  G4ThreeVector positionVacD(xDumVac, yDumVac, zDumVac);

  G4ThreeVector posDumVacA(xDumVac, -yDumVac, zDumVac);
  G4Transform3D positionVacA(rm270, posDumVacA);

  G4ThreeVector posDumVacB(-xDumVac, -yDumVac, zDumVac);
  G4Transform3D positionVacB(rm180, posDumVacB);

  G4ThreeVector posDumVacC(-xDumVac, yDumVac, zDumVac);
  G4Transform3D positionVacC(rm90, posDumVacC);
  
  new G4PVPlacement(positionVacA, logicOuterA1, GeA_Outer1, logicVac, false, 1,true);
  new G4PVPlacement(positionVacA, logicOuterA2, GeA_Outer2, logicVac, false, 2,true);
  new G4PVPlacement(positionVacA, logicOuterA3, GeA_Outer3, logicVac, false, 3,true);
  new G4PVPlacement(positionVacA, logicOuterA4, GeA_Outer4, logicVac, false, 4,true);
  
  new G4PVPlacement(positionVacB, logicOuterB1, GeB_Outer1, logicVac, false, 5,true);
  new G4PVPlacement(positionVacB, logicOuterB2, GeB_Outer2, logicVac, false, 6,true);
  new G4PVPlacement(positionVacB, logicOuterB3, GeB_Outer3, logicVac, false, 7,true);
  new G4PVPlacement(positionVacB, logicOuterB4, GeB_Outer4, logicVac, false, 8,true);
  
  new G4PVPlacement(positionVacC, logicOuterC1, GeC_Outer1, logicVac, false, 9,true);
  new G4PVPlacement(positionVacC, logicOuterC2, GeC_Outer2, logicVac, false, 10,true);
  new G4PVPlacement(positionVacC, logicOuterC3, GeC_Outer3, logicVac, false, 11,true);
  new G4PVPlacement(positionVacC, logicOuterC4, GeC_Outer4, logicVac, false, 12,true);
  
  new G4PVPlacement(0,positionVacD, logicOuterD1, GeD_Outer1, logicVac, false, 13,true);
  new G4PVPlacement(0,positionVacD, logicOuterD2, GeD_Outer2, logicVac, false, 14,true);
  new G4PVPlacement(0,positionVacD, logicOuterD3, GeD_Outer3, logicVac, false, 15,true);
  new G4PVPlacement(0,positionVacD, logicOuterD4, GeD_Outer4, logicVac, false, 16,true);


  //
  // some material between the diodes to reproduce the experimental addback factor ...
  //

  G4double xAbsorb1 = 4.16 * cm;
  G4double yAbsorb1 = 200. * um; // max = HalfDistanceBetweenDiodes = 0.5*mm;
  G4double zAbsorb1 = 4.5 * cm;

  G4Box* solidAbsorb1 = new G4Box("Absorb1", xAbsorb1, yAbsorb1, zAbsorb1);

  G4double xAbsorb2 = 200 * um; // max = HalfDistanceBetweenDiodes = 0.5*mm;
  G4double yAbsorb2 = 4.16 * cm;
  G4double zAbsorb2 = 4.5 * cm;

  G4Box* solidAbsorb2 = new G4Box("Absorb2", xAbsorb2, yAbsorb2, zAbsorb2);

  // G4UnionSolid* solidAbsorb =
  // new G4UnionSolid("Absorb",solidAbsorb1,solidAbsorb2,0,0);
  G4UnionSolid* solidAbsorb = new G4UnionSolid("Absorb", solidAbsorb1, solidAbsorb2);

  G4LogicalVolume* logicAbsorb = new G4LogicalVolume(solidAbsorb, m_Copper, "Absorb");

  G4ThreeVector positionAbsorb(0., 0., zDumVac);

  new G4PVPlacement(0, positionAbsorb, logicAbsorb, "Absorb", logicVac, false, CloverNbr, true);

  //
  // Now: takes care of the holes and amorphous Ge in each diode:
  // Central hole with amorphous Ge for each diode.
  //

  G4double minRadiusAGe1 = 0. * cm;
  G4double maxRadiusAGe1 = 0.52 * cm;
  G4double HalfLengthAGe1 = 3.75 * cm;
  G4double startPhiAGe1 = 0. * deg;
  G4double deltaPhiAGe1 = 360. * deg;

  // G4Tubs* solidAGe1 = new G4Tubs("AGe1",minRadiusAGe1,maxRadiusAGe1,
  //			 HalfLengthAGe1,startPhiAGe1,deltaPhiAGe1);

  // G4LogicalVolume* logicAGe1 = new G4LogicalVolume(solidAGe1,m_Germanium,"AGe1");

  // ... and second the hole in it:

  G4Tubs* solidHole1 =
      new G4Tubs("Hole1", minRadiusAGe1, maxRadiusAGe1 - 2. * mm, HalfLengthAGe1, startPhiAGe1, deltaPhiAGe1);

  G4LogicalVolume* logicHole1 = new G4LogicalVolume(solidHole1, m_Vacuum, "Hole1");

  // Visu
  G4VisAttributes* CanVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.7)); // Grey
  G4VisAttributes* DewarVisAtt = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));    // Grey

  logicCloverCan->SetVisAttributes(CanVisAtt);
  logicEnvColdFinger->SetVisAttributes(CanVisAtt);
  logicDewar->SetVisAttributes(DewarVisAtt);
  logicSupClover->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VisAttributes* HoleVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0));      // Black
  G4VisAttributes* AbsorbVisAtt = new G4VisAttributes(G4Colour(0.5, 0.0, 0.5, 1)); // purple
  G4VisAttributes* GeAVisAtt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.1));  // Red
  G4VisAttributes* GeBVisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0, 0.6));  // Green
  G4VisAttributes* GeCVisAtt = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.6));  // Blue
  G4VisAttributes* GeDVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.6));  // White
  G4VisAttributes* VisOuter1 = new G4VisAttributes(G4Colour(0.5, 0.5, 0.0, 1.0));  // White
  G4VisAttributes* VisOuter2 = new G4VisAttributes(G4Colour(0.5, 0.7, 0.0, 1.0));  // White
  G4VisAttributes* VisOuter3 = new G4VisAttributes(G4Colour(0.2, 0.0, 0.8, 1.0));  // White
  G4VisAttributes* VisOuter4 = new G4VisAttributes(G4Colour(0.8, 0.0, 0.2, 1.0));  // White
  
  logicOuterA1->SetVisAttributes(VisOuter1);
  logicOuterA2->SetVisAttributes(VisOuter2);
  logicOuterA3->SetVisAttributes(VisOuter3);
  logicOuterA4->SetVisAttributes(VisOuter4);
  
  logicOuterB1->SetVisAttributes(VisOuter1);
  logicOuterB2->SetVisAttributes(VisOuter2);
  logicOuterB3->SetVisAttributes(VisOuter3);
  logicOuterB4->SetVisAttributes(VisOuter4);
  
  logicOuterC1->SetVisAttributes(VisOuter1);
  logicOuterC2->SetVisAttributes(VisOuter2);
  logicOuterC3->SetVisAttributes(VisOuter3);
  logicOuterC4->SetVisAttributes(VisOuter4);
  
  logicOuterD1->SetVisAttributes(VisOuter1);
  logicOuterD2->SetVisAttributes(VisOuter2);
  logicOuterD3->SetVisAttributes(VisOuter3);
  logicOuterD4->SetVisAttributes(VisOuter4);

  logicHole1->SetVisAttributes(HoleVisAtt);
  logicAbsorb->SetVisAttributes(AbsorbVisAtt);
  logicVac->SetVisAttributes(G4VisAttributes::GetInvisible());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void Exogam::InitializeRootOutput() {
  RootOutput* pAnalysis = RootOutput::getInstance();
  TTree* pTree = pAnalysis->GetTree();
  if (!pTree->FindBranch("Exogam")) {
    pTree->Branch("Exogam", "TExogamCalData", &m_Event);
  }
  pTree->SetBranchAddress("Exogam", &m_Event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void Exogam::ReadSensitive(const G4Event*) {
  m_Event->Clear();

  ///////////
  // Calorimeter scorer
  CalorimeterScorers::PS_Calorimeter* Scorer = (CalorimeterScorers::PS_Calorimeter*)m_ExogamScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  if (size > 0){
    // std::cout << "/////////: " << size << " "<< Scorer->GetName()<< std::endl;
    std::map<int,std::map<int,std::pair<double,std::map<int,double>>>> CrystalsEnergy;
    for (unsigned int i = 0; i < size; i++) {
      int CloverNbr = Scorer->GetLevel(i)[1];
      int CristalNbr = int((Scorer->GetLevel(i)[0]-1)/4);
      int OuterNbr = int((Scorer->GetLevel(i)[0]-1)%4);

      double Energy = Scorer->GetEnergy(i) * MeV;

      CrystalsEnergy[CloverNbr][CristalNbr].first = 0;
      CrystalsEnergy[CloverNbr][CristalNbr].second[OuterNbr] = Energy;
    }

// Summing energies
for(auto& Clover: CrystalsEnergy){
  for(auto& Crystal: Clover.second){
    for(auto& Outer: Crystal.second.second){
          Crystal.second.first += Outer.second;
    }
  }
}

// Applying Gaussian noise
for(auto& Clover: CrystalsEnergy){
  for(auto& Crystal: Clover.second){
    for(auto& Outer: Crystal.second.second){
          double EnergyOuter = RandGauss::shoot(Outer.second, Exogam_NS::ResoEnergyOuter);
          if(EnergyOuter > Exogam_NS::EnergyThresholdOuter)
            Outer.second = EnergyOuter;
          else
            Outer.second = 0;
        }
        double EnergyCrystal =  RandGauss::shoot(Crystal.second.first, Exogam_NS::ResoEnergyInner6MV);
        if(EnergyCrystal > Exogam_NS::EnergyThresholdInner6MV){
          Crystal.second.first = EnergyCrystal;
        }
        else{
          Crystal.second.first = 0;
        }
      }
    }

    for(auto Clover: CrystalsEnergy){
      for(auto Crystal: Clover.second){
        double Outers[4] = {0.,0.,0.,0.};
          for(auto Outer: Crystal.second.second){
            Outers[Outer.first] = Outer.second;
          }
        if(Crystal.second.first > 0)
        m_Event->SetExo(Crystal.first + Clover.first*4, Crystal.second.first, -1000, -1000, -1000, -1000, -1000, Outers[0], Outers[1], Outers[2], Outers[3]);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void Exogam::InitializeScorers() {
  bool already_exist = false;
  m_ExogamScorer = CheckScorer("ExogamScorer", already_exist);

  if (already_exist)
    return;

  // Otherwise the scorer is initialised
  vector<int> level({0, 1});

  // Calorimeter Scorer
  m_ExogamScorer->RegisterPrimitive(new CalorimeterScorers::PS_Calorimeter("Cristal", level, 0));
  //Interaction Scorer
  G4VPrimitiveScorer* Interaction = new InteractionScorers::PS_Interactions("Interaction", ms_InterCoord, 0);
  m_ExogamScorer->RegisterPrimitive(Interaction);

  G4SDManager::GetSDMpointer()->AddNewDetector(m_ExogamScorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* Exogam::Construct() { return (NPS::VDetector*)new Exogam(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_nps_Exogam {
 public:
  proxy_nps_Exogam() {
    NPS::DetectorFactory::getInstance()->AddToken("Exogam", "Exogam");
    NPS::DetectorFactory::getInstance()->AddDetector("Exogam", Exogam::Construct);
  }
};

proxy_nps_Exogam p_nps_Exogam;
}
