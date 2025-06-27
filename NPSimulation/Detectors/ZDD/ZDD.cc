/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ZDD simulation                             *
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
#include "ZDD.hh"
#include "DriftChamberScorers.hh"
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
namespace ZDD_NS{
// Energy and time Resolution
/* const double EnergyThreshold = 0.1*MeV;
  const double ResoTime = 4.5*ns ;
  const double ResoEnergy = 1.0*MeV ;
  const double Radius = 50*mm ;
  const double Width = 100*mm ;
  const double Thickness = 300*mm ;
  const string Material = "BC400";
  */
// IC E and T
const double EnergyThreshold_IC = 0.*MeV;
const double ResoTime_IC = 4.5*ns ;
const double ResoEnergy_IC = 0.; //Resolution en pourcentage

// Plastic E and T
const double EnergyThreshold_Plastic = 0.1*MeV;
const double ResoTime_Plastic = 4.5*ns ;
const double ResoEnergy_Plastic = 5; // Resolution en pourcentage

// Drift Time and position of drift chambers
const G4double ResoDriftTime   = 0.0001 * ns;
//const G4double ResoDriftTime   = 0. * ns;

// Drift features
const G4double      DriftSpeed = 450.0 * cm /microsecond; //microsecond; RQ très importante: le temps enregistré par les DC est un temps interne, ie il prend uniquement en compte le temps de dérive dans les drift
// En sortant des scorers, ce temps est en us
const G4ThreeVector DriftDir1   = G4ThreeVector(1, 0, 0);
const G4ThreeVector DriftDir2   = G4ThreeVector(0, 1, 0);


// Drift Chambers
const G4double Drift_Chamber_Width     = 180 * mm;
const G4double Drift_Chamber_Length    = 180 * mm;
const G4double Drift_Chamber_Thickness = 20 * mm;
const G4double Drift_Chamber_Mylar = 0.9 * um;


// ZDD Volume
const G4double Phi                  = 0 * deg;
const G4double ZDD_Volume_Width     = 4500 * mm;
const G4double ZDD_Volume_Length    = 1000 * mm;
const G4double ZDD_Volume_Thickness = 5000 * mm;

// Anodes and Cathodes in IC (Mylar)
const G4double AC_Width     = 400 * mm;
const G4double AC_Length    = 400 * mm;
const G4double AC_Thickness = 1.5 * um;
//const G4double AC_Angle;

// IC Entrance and Exit (Kapton)
const G4double Mylar_Foil_Width     = 400 * mm;
const G4double Mylar_Foil_Length    = 400 * mm;
const G4double Mylar_Foil_Thickness = 20. * um;
const G4double Mylar_Foil_Angle = 30 * deg;

// Ionisation Chambers
const G4double Ionisation_Chamber_Width     = 400 * mm;
const G4double Ionisation_Chamber_Length    = 400 * mm;
//const G4double Ionisation_Chamber_Thickness = 40 * mm;
const G4double Ionisation_Chamber_Angle = 30 * deg;

// Gas Gap
const G4double Gas_Gap_Width     = 400 * mm;
const G4double Gas_Gap_Length    = 400 * mm;
//const G4double Gas_Gap_Thickness = 40 * mm;
const G4double Gas_Gap_Angle = 30 * deg;

//Plastic Dimension
const G4double Plastic_Width     = 220. * mm;
const G4double Plastic_Length    = 44. * mm;
const G4double Plastic_Thickness = 20. * mm;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ZDD Specific Method
ZDD::ZDD(){
    m_Event = new TZDDData() ;
    m_ZDDScorer = 0;

    ClearGeometry();
    m_IC_Scorer = 0;
    m_Plastic_Scorer = 0;
    m_Drift_Chamber_Scorer_1 = 0;
    m_Drift_Chamber_Scorer_2 = 0;
    m_Drift_Chamber_1 = 0;
    m_Drift_Chamber_2 = 0;

    ICcounter = 0;
    ACcounter = 0;
    GGcounter = 0;
    Entry_Exit_counter = 0;
    Plasticcounter = 0;

    Mylar_mat = MaterialManager::getInstance()->GetMaterialFromLibrary("Mylar");
    m_Vis_Mylar = new G4VisAttributes(G4Colour(1, 1, 0, 1.));


    // RGB Color + Transparency
    m_VisSquare = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));
    m_VisCylinder = new G4VisAttributes(G4Colour(0, 0, 1, 0.5));
    m_Vis_ZDD = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.1));
    m_Vis_Drift_Chamber_Gas = new G4VisAttributes(G4Colour(0, 1, 0, 0.5));
    m_Vis_Ionisation_Chamber_Gas  = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0, 0.3));
    m_Vis_AC  = new G4VisAttributes(G4Colour(0.0, 0.5, 0.5, 0.5));
    m_Vis_Plastic     = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0, 0.8));

}

ZDD::~ZDD(){
}

using namespace ZDD_NS;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*void ZDD::AddDetector(G4ThreeVector POS, string  Shape){
  // Convert the POS value to R theta Phi as Spherical coordinate is easier in G4
  m_R.push_back(POS.mag());
  m_Theta.push_back(POS.theta());
  m_Phi.push_back(POS.phi());
  m_Shape.push_back(Shape);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ZDD::AddDetector(double  R, double  Theta, double  Phi, string  Shape){
  m_R.push_back(R);
  m_Theta.push_back(Theta);
  m_Phi.push_back(Phi);
  m_Shape.push_back(Shape);
}
*/
void ZDD::Add_ZDD(G4double R, double Theta) {
    m_R     = R;
    m_Theta = Theta;
}
void ZDD::Add_Drift_Chamber(G4double Z, double Thickness, string Gas, double Pressure,double Temperature) {
    m_Drift_Chamber_Z.push_back(Z);
    m_Drift_Chamber_Gas.push_back(Gas);
    m_Drift_Chamber_Pressure.push_back(Pressure);
    m_Drift_Chamber_Temperature.push_back(Temperature);
    m_Drift_Chamber_Thickness.push_back(Thickness);
}

void ZDD::Add_Ionisation_Chamber(G4int nDet, G4double Z, double Thickness, string Gas, double Pressure,
                                 double Temperature) {
    m_Ionisation_Chamber_nDet.push_back(nDet);
    m_Ionisation_Chamber_Z.push_back(Z);
    m_Ionisation_Chamber_Thickness.push_back(Thickness);
    m_Ionisation_Chamber_Gas.push_back(Gas);
    m_Ionisation_Chamber_Pressure.push_back(Pressure);
    m_Ionisation_Chamber_Temperature.push_back(Temperature);
}

void ZDD::Add_Gas_Gap(G4double Z, double Thickness, string Gas, double Pressure,
                      double Temperature) {
    m_Gas_Gap_Z.push_back(Z);
    m_Gas_Gap_Thickness.push_back(Thickness);
    m_Gas_Gap_Gas.push_back(Gas);
    m_Gas_Gap_Pressure.push_back(Pressure);
    m_Gas_Gap_Temperature.push_back(Temperature);
}

void ZDD::Add_AC(G4double Z, double Thickness, string Material, G4double Angle) {
    m_AC_Z.push_back(Z);
    m_AC_Thickness.push_back(Thickness);
    m_AC_Material.push_back(Material);
    m_AC_Angle.push_back(Angle);
}

void ZDD::Add_Entry_Exit(G4double Z, double Thickness, string Material) {
    m_Entry_Exit_Z.push_back(Z);
    m_Entry_Exit_Thickness.push_back(Thickness);
    m_Entry_Exit_Material.push_back(Material);
}

void ZDD::Add_Plastic(string Material, G4double Width, double Length,
                      double Thickness, G4ThreeVector Pos) {
    m_Plastic_Material.push_back(Material);
    m_Plastic_Width.push_back(Width);
    m_Plastic_Length.push_back(Length);
    m_Plastic_Thickness.push_back(Thickness);
    m_Plastic_Position.push_back(Pos);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*G4LogicalVolume* ZDD::BuildSquareDetector(){
  if(!m_SquareDetector){
    G4Box* box = new G4Box("ZDD_Box",ZDD_NS::Width*0.5,
        ZDD_NS::Width*0.5,ZDD_NS::Thickness*0.5);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ZDD_NS::Material);
    m_SquareDetector = new G4LogicalVolume(box,DetectorMaterial,"logic_ZDD_Box",0,0,0);
    m_SquareDetector->SetVisAttributes(m_VisSquare);
    m_SquareDetector->SetSensitiveDetector(m_ZDDScorer);
  }
  return m_SquareDetector;
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*G4LogicalVolume* ZDD::BuildCylindricalDetector(){
  if(!m_CylindricalDetector){
    G4Tubs* tub = new G4Tubs("ZDD_Cyl",0,ZDD_NS::Radius,ZDD_NS::Thickness*0.5,0,360*deg);

    G4Material* DetectorMaterial = MaterialManager::getInstance()->GetMaterialFromLibrary(ZDD_NS::Material);
    m_CylindricalDetector = new G4LogicalVolume(tub,DetectorMaterial,"logic_ZDD_tub",0,0,0);
    m_CylindricalDetector->SetVisAttributes(m_VisSquare);
    m_CylindricalDetector->SetSensitiveDetector(m_ZDDScorer);

  }
  return m_CylindricalDetector;
}*/

G4LogicalVolume* ZDD::Build_Drift_Chamber_1() {
    if (!m_Drift_Chamber_1) {
        G4Box* box = new G4Box("ZDD_Drift_Chamber_1", Drift_Chamber_Width * 0.5, Drift_Chamber_Length * 0.5,
                               m_Drift_Chamber_Thickness[0] * 0.5);

        G4Material* DetectorMaterial
                = MaterialManager::getInstance()->GetGasFromLibrary(
                    m_Drift_Chamber_Gas[0], m_Drift_Chamber_Pressure[0], m_Drift_Chamber_Temperature[0]);

        m_Drift_Chamber_1 = new G4LogicalVolume(box, DetectorMaterial, "logic_ZDD_Drift_Chamber_1", 0, 0, 0);
        m_Drift_Chamber_1->SetVisAttributes(m_Vis_Drift_Chamber_Gas);
        m_Drift_Chamber_1->SetSensitiveDetector(m_Drift_Chamber_Scorer_1);
        //        new G4PVPlacement(0, G4ThreeVector(0, 0, m_Drift_Chamber_Thickness[0] * 0.5 + m_Drift_Chamber_Z[0]),
        //                Build_Drift_Chamber_1(), "Drift_Chamber_1", world, false, 1);

        G4Box* Mylar_box = new G4Box("ZDD_Mylar_Drift_Chamber", Drift_Chamber_Width * 0.5, Drift_Chamber_Length * 0.5,
                                     Drift_Chamber_Mylar * 0.5);

        G4LogicalVolume * m_Mylar_Drift_Chamber = new G4LogicalVolume(Mylar_box, Mylar_mat, "logic_ZDD_Mylar_Drift_Chamber", 0, 0, 0);
        m_Mylar_Drift_Chamber->SetVisAttributes(m_Vis_Mylar);
        new G4PVPlacement(0, G4ThreeVector(0, 0, m_Drift_Chamber_Thickness[0] * 0.5 - Drift_Chamber_Mylar * 0.5 ),
                m_Mylar_Drift_Chamber, "ZDD_Mylar_Drift_Chamber", m_Drift_Chamber_1, false, 1);

        new G4PVPlacement(0, G4ThreeVector(0, 0, -m_Drift_Chamber_Thickness[0] * 0.5 + Drift_Chamber_Mylar * 0.5),
                m_Mylar_Drift_Chamber, "ZDD_Mylar_Drift_Chamber", m_Drift_Chamber_1, false, 2);
    }
    return m_Drift_Chamber_1;

}

G4LogicalVolume* ZDD::Build_Drift_Chamber_2() {
    if (!m_Drift_Chamber_2) {
        G4Box* box = new G4Box("ZDD_Drift_Chamber_2", Drift_Chamber_Width * 0.5, Drift_Chamber_Length * 0.5,
                               m_Drift_Chamber_Thickness[1] * 0.5);

        G4Material* DetectorMaterial
                = MaterialManager::getInstance()->GetGasFromLibrary(
                    m_Drift_Chamber_Gas[1], m_Drift_Chamber_Pressure[1], m_Drift_Chamber_Temperature[1]);

        m_Drift_Chamber_2 = new G4LogicalVolume(box, DetectorMaterial, "logic_ZDD_Drift_Chamber_2", 0, 0,
                                                0);
        m_Drift_Chamber_2->SetVisAttributes(m_Vis_Drift_Chamber_Gas);
        m_Drift_Chamber_2->SetSensitiveDetector(m_Drift_Chamber_Scorer_2);

        G4Box* Mylar_box = new G4Box("ZDD_Mylar_Drift_Chamber", Drift_Chamber_Width * 0.5, Drift_Chamber_Length * 0.5,
                                     Drift_Chamber_Mylar * 0.5);

        G4LogicalVolume * m_Mylar_Drift_Chamber = new G4LogicalVolume(Mylar_box, Mylar_mat, "logic_ZDD_Mylar_Drift_Chamber", 0, 0, 0);
        m_Mylar_Drift_Chamber->SetVisAttributes(m_Vis_Mylar);
        new G4PVPlacement(0, G4ThreeVector(0, 0, m_Drift_Chamber_Thickness[1] * 0.5 - Drift_Chamber_Mylar * 0.5),
                m_Mylar_Drift_Chamber, "ZDD_Mylar_Drift_Chamber", m_Drift_Chamber_2, false, 1);

        new G4PVPlacement(0, G4ThreeVector(0, 0, -m_Drift_Chamber_Thickness[1] * 0.5 + Drift_Chamber_Mylar * 0.5),
                m_Mylar_Drift_Chamber, "ZDD_Mylar_Drift_Chamber", m_Drift_Chamber_2, false, 2);
    }
    return m_Drift_Chamber_2;
}

void ZDD::Build_Ionisation_Chamber(G4LogicalVolume* world, G4int Det)
{
    G4Material* GasIC = MaterialManager::getInstance()->GetGasFromLibrary(
                m_Ionisation_Chamber_Gas[0], m_Ionisation_Chamber_Pressure[0], m_Ionisation_Chamber_Temperature[0]);
    G4Material* VacuumMaterial
            = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");

    if (Det == m_Ionisation_Chamber_nDet[0]-1)
    {
        G4Box* IC_MotherVolume = new G4Box("ZDD_IC_MotherVolume", Ionisation_Chamber_Width * 0.5, Ionisation_Chamber_Length * 0.5,
                               m_Ionisation_Chamber_Thickness[0] + 1.5*AC_Thickness);
        G4LogicalVolume* logic_IC_MotherVolume = new G4LogicalVolume(IC_MotherVolume, VacuumMaterial, "logic_IC_MotherVolume", 0, 0,0);

        G4Box* IC_box = new G4Box("ZDD_IC", Ionisation_Chamber_Width * 0.5, Ionisation_Chamber_Length * 0.5,
                               m_Ionisation_Chamber_Thickness[0] * 0.5);
        G4Box* Cathode_Box = new G4Box("ZDD_IC_Mylar", AC_Width * 0.5, AC_Length * 0.5, AC_Thickness * 0.5);
        G4Box* Anode_Box = new G4Box("ZDD_IC_Mylar", AC_Width * 0.5, AC_Length * 0.5, AC_Thickness * 0.5);

        G4LogicalVolume* logic_IC = new G4LogicalVolume(IC_box, GasIC, "logic_ZDD_IC", 0, 0, 0);
        logic_IC->SetVisAttributes(m_Vis_Ionisation_Chamber_Gas);
        logic_IC->SetSensitiveDetector(m_IC_Scorer);
        G4LogicalVolume* logic_Cathode_Box = new G4LogicalVolume(Cathode_Box, Mylar_mat, "logic_ZDD_Cathode_Box", 0, 0, 0);
        logic_Cathode_Box->SetVisAttributes(m_Vis_Mylar);
        G4LogicalVolume* logic_Anode_Box = new G4LogicalVolume(Anode_Box, Mylar_mat, "logic_ZDD_Anode_Box", 0, 0, 0);
        logic_Anode_Box->SetVisAttributes(m_Vis_Mylar);

        G4RotationMatrix* Rot_Ionisation_Chamber = new G4RotationMatrix();
        Rot_Ionisation_Chamber->rotateX(Ionisation_Chamber_Angle);

        new G4PVPlacement(0, G4ThreeVector(0, 0, -m_Ionisation_Chamber_Thickness[0] - AC_Thickness),
                          logic_Cathode_Box, "Cathode_Box", logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0, G4ThreeVector(0,0, -m_Ionisation_Chamber_Thickness[0]*0.5 - AC_Thickness * 0.5),
                          logic_IC, Form("IC_Box_%d1",Det), logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
                          logic_Anode_Box, "Anode_Box", logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0,G4ThreeVector(0,0, m_Ionisation_Chamber_Thickness[0]*0.5 + 0.5*AC_Thickness),
                          logic_IC, Form("IC_Box_%d2",Det), logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0, G4ThreeVector(0, 0, m_Ionisation_Chamber_Thickness[0]+ AC_Thickness),
                          logic_Cathode_Box, "Cathode_Box", logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(G4Transform3D(*Rot_Ionisation_Chamber,G4ThreeVector(0, 0, m_Ionisation_Chamber_Z[0]+ (2*(m_Ionisation_Chamber_nDet[0]-2)+1)*(m_Ionisation_Chamber_Thickness[0] + AC_Thickness)/cos(Ionisation_Chamber_Angle/rad)
                          + (2.5*AC_Thickness + 2.*m_Ionisation_Chamber_Thickness[0])/cos(Ionisation_Chamber_Angle/rad) )),
                          logic_IC_MotherVolume, "IC_MotherVolume", world, false, Det);
    }
    else
    {
        G4Box* IC_MotherVolume = new G4Box("ZDD_IC_MotherVolume", Ionisation_Chamber_Width * 0.5, Ionisation_Chamber_Length * 0.5,
                               m_Ionisation_Chamber_Thickness[0] + AC_Thickness);
        G4LogicalVolume* logic_IC_MotherVolume = new G4LogicalVolume(IC_MotherVolume, VacuumMaterial, "logic_IC_MotherVolume", 0, 0, 0);

        G4Box* IC_box = new G4Box("ZDD_IC", Ionisation_Chamber_Width * 0.5, Ionisation_Chamber_Length * 0.5,
                               m_Ionisation_Chamber_Thickness[0] * 0.5);
        G4Box* Cathode_Box = new G4Box("ZDD_IC_Mylar", AC_Width * 0.5, AC_Length * 0.5, AC_Thickness * 0.5);
        G4Box* Anode_Box = new G4Box("ZDD_IC_Mylar", AC_Width * 0.5, AC_Length * 0.5, AC_Thickness * 0.5);

        G4LogicalVolume* logic_IC = new G4LogicalVolume(IC_box, GasIC, "logic_ZDD_IC", 0, 0, 0);
        logic_IC->SetVisAttributes(m_Vis_Ionisation_Chamber_Gas);
        logic_IC->SetSensitiveDetector(m_IC_Scorer);
        G4LogicalVolume* logic_Cathode_Box = new G4LogicalVolume(Cathode_Box, Mylar_mat, "logic_ZDD_Cathode_Box", 0, 0, 0);
        logic_Cathode_Box->SetVisAttributes(m_Vis_Mylar);
        G4LogicalVolume* logic_Anode_Box = new G4LogicalVolume(Anode_Box, Mylar_mat, "logic_ZDD_Anode_Box", 0, 0, 0);
        logic_Anode_Box->SetVisAttributes(m_Vis_Mylar);

        G4RotationMatrix* Rot_Ionisation_Chamber = new G4RotationMatrix();
        Rot_Ionisation_Chamber->rotateX(Ionisation_Chamber_Angle);

        new G4PVPlacement(0, G4ThreeVector(0, 0, -m_Ionisation_Chamber_Thickness[0] - 0.5*AC_Thickness),
                          logic_Cathode_Box, "Cathode_Box", logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0, G4ThreeVector(0,0, -m_Ionisation_Chamber_Thickness[0]*0.5),
                          logic_IC, Form("IC_Box_%d1",Det), logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0, G4ThreeVector(0, 0, AC_Thickness * 0.5),
                          logic_Anode_Box, "Anode_Box", logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(0,G4ThreeVector(0,0, m_Ionisation_Chamber_Thickness[0]*0.5 + AC_Thickness),
                          logic_IC, Form("IC_Box_%d2",Det), logic_IC_MotherVolume, false, Det);
        new G4PVPlacement(G4Transform3D(*Rot_Ionisation_Chamber,G4ThreeVector(0, 0, m_Ionisation_Chamber_Z[0]+ (2*Det+1)*(m_Ionisation_Chamber_Thickness[0] + AC_Thickness)/cos(Ionisation_Chamber_Angle/rad))),
                          logic_IC_MotherVolume, "IC_MotherVolume", world, false, Det);
    }
}

void ZDD::ClearGeometry(){

    m_Drift_Chamber_Z.clear();
    m_Drift_Chamber_Gas.clear();
    m_Drift_Chamber_Pressure.clear();
    m_Drift_Chamber_Temperature.clear();
    m_Drift_Chamber_Thickness.clear();

    m_Ionisation_Chamber_Z.clear();
    m_Ionisation_Chamber_Gas.clear();
    m_Ionisation_Chamber_Pressure.clear();
    m_Ionisation_Chamber_Temperature.clear();
    m_Ionisation_Chamber_Thickness.clear();

    m_Gas_Gap_Z.clear();
    m_Gas_Gap_Gas.clear();
    m_Gas_Gap_Pressure.clear();
    m_Gas_Gap_Temperature.clear();
    m_Gas_Gap_Thickness.clear();

    m_AC_Material.clear();
    m_AC_Rotation.clear();
    m_AC_Thickness.clear();
    m_AC_Z.clear();
    m_AC_Angle.clear();

    m_Entry_Exit_Material.clear();
    m_Entry_Exit_Z.clear();
    m_Entry_Exit_Thickness.clear();

    m_Plastic_Length.clear();
    m_Plastic_Width.clear();
    m_Plastic_Thickness.clear();
    m_Plastic_Position.clear();
    m_Plastic_Material.clear();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Virtual Method of NPS::VDetector class

// Read stream at Configfile to pick-up parameters of detector (Position,...)
// Called in DetecorConstruction::ReadDetextorConfiguration Method
void ZDD::ReadConfiguration(NPL::InputParser parser){
    //ClearGeometry();

    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ZDD");
    if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    vector<string> TokenZDD  = {"R", "Theta"};
    vector<string> TokenDC = {"Z","Thickness", "Gas", "Pressure", "Temperature"};
    vector<string> TokenIC = {"NDet","ZofFirst", "Thickness", "Gas", "Pressure", "Temperature"};
    vector<string> TokenAC = {"Z", "Thickness", "Material"};
    vector<string> TokenEntryExit = {"Z", "Thickness", "Material"};
    vector<string> TokenPlastic = {"NDet","Material", "Width", "Length", "Thickness", "PosOfCentral"};

    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        /*if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  ZDD " << i+1 <<  endl;

      G4ThreeVector Pos = NPS::ConvertVector(blocks[i]->GetTVector3("POS","mm"));
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  ZDD " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,Shape);
    }*/
        if (blocks[i]->HasTokenList(TokenZDD)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  ZDD " << i + 1 << endl;
            G4double R     = blocks[i]->GetDouble("R", "mm");
            G4double Theta = blocks[i]->GetDouble("Theta", "deg");
            Add_ZDD(R, Theta);
        }
        else if (blocks[i]->GetMainValue() == "DC"
                 && blocks[i]->HasTokenList(TokenDC)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "//// DC " << i + 1 << endl;
            G4double Z           = blocks[i]->GetDouble("Z", "mm");
            G4double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
            string   Gas         = blocks[i]->GetString("Gas");
            G4double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
            G4double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
            Add_Drift_Chamber(Z, Thickness, Gas, Pressure, Temperature);
        }
        else if (blocks[i]->GetMainValue() == "IC"
                 && blocks[i]->HasTokenList(TokenIC)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "//// IC " << endl;
            G4int         nDet =  blocks[i]->GetInt("NDet");
            G4bool        checkEven = (nDet%2 ==0);
            G4double Z           = blocks[i]->GetDouble("ZofFirst", "mm");
            G4double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
            string   Gas         = blocks[i]->GetString("Gas");
            G4double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
            G4double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
            Add_Ionisation_Chamber(nDet,Z, Thickness, Gas, Pressure, Temperature);
            ICcounter++;
        }
//        else if (blocks[i]->GetMainValue() == "GasGap"
//                 && blocks[i]->HasTokenList(TokenIC)) {
//            if (NPOptionManager::getInstance()->GetVerboseLevel())
//                cout << endl << "//// GasGap " << GGcounter+1 << endl;
//            G4double Z           = blocks[i]->GetDouble("Z", "mm");
//            G4double Thickness   = blocks[i]->GetDouble("Thickness", "mm");
//            string   Gas         = blocks[i]->GetString("Gas");
//            G4double Pressure    = blocks[i]->GetDouble("Pressure", "bar");
//            G4double Temperature = blocks[i]->GetDouble("Temperature", "kelvin");
//            Add_Gas_Gap(Z, Thickness, Gas, Pressure, Temperature);
//            GGcounter++;
//        }
//        else if (blocks[i]->GetMainValue() == "AC"
//                 && blocks[i]->HasTokenList(TokenAC)) {
//            if (NPOptionManager::getInstance()->GetVerboseLevel())
//                cout << endl << "//// AC " << ACcounter+1 << endl;
//            G4double Z           = blocks[i]->GetDouble("Z", "mm");
//            G4double Thickness   = blocks[i]->GetDouble("Thickness", "um");
//            string   Material    = blocks[i]->GetString("Material");
//            G4double Angle           = blocks[i]->GetDouble("Theta","deg");
//            Add_AC(Z, Thickness, Material, Angle);
//            ACcounter++;
//        }
        else if (blocks[i]->GetMainValue() == "EntryExit"
                 && blocks[i]->HasTokenList(TokenEntryExit)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "//// AC" << Entry_Exit_counter+1 << endl;
            G4double Z           = blocks[i]->GetDouble("Z", "mm");
            G4double Thickness   = blocks[i]->GetDouble("Thickness", "um");
            string   Material    = blocks[i]->GetString("Material");
            Add_Entry_Exit(Z, Thickness, Material);
            Entry_Exit_counter++;
        }
        else if (blocks[i]->GetMainValue() == "Plastic"
                 && blocks[i]->HasTokenList(TokenPlastic))
        {
            cout << endl << "//// Plastic " << endl;
            G4int         nDet =  blocks[i]->GetInt("NDet");
            G4bool        checkEven = (nDet%2 ==0);
            G4ThreeVector Pos  = NPS::ConvertVector(blocks[i]->GetTVector3("PosOfCentral", "mm"));
            string        Material  = blocks[i]->GetString("Material");
            G4double      Width     = blocks[i]->GetDouble("Width", "mm");
            G4double      Length    = blocks[i]->GetDouble("Length", "mm");
            G4double      Thickness = blocks[i]->GetDouble("Thickness", "mm");
            for (G4int i = 0; i < nDet; i++)
            {
                if (NPOptionManager::getInstance()->GetVerboseLevel())
                    cout << endl << "       Plastic " << Plasticcounter << endl;
                if (checkEven)
                {
                    G4ThreeVector nPos(Pos.x(),(nDet/2-i-0.5)*Length+Pos.y(),Pos.z()); //Center between i/2-1 & i/2
                    Add_Plastic(Material, Width, Length, Thickness, nPos);
                    printf(" Pos : (%.4f;%.4f;%.4f)\n", nPos.x(),nPos.y(),nPos.z());
                }
                else
                {
                    G4ThreeVector nPos(Pos.x(),(nDet/2-i)*Length+Pos.y(),Pos.z()); //Center in i/2
                    Add_Plastic(Material, Width, Length, Thickness, nPos);
                    printf(" Pos : (%.4f;%.4f;%.4f)\n", nPos.x(),nPos.y(),nPos.z());
                }
                Plasticcounter++;
            }
        }
        else{
            cout << "ERROR: check your input file formatting --- Reading blocks : " << blocks[i]->GetMainValue() << endl;
            exit(1);
        }
    }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Construct detector and inialise sensitive part.
// Called After DetecorConstruction::AddDetector Method
void ZDD::ConstructDetector(G4LogicalVolume* world){
    //  for (unsigned short i = 0 ; i < m_R.size() ; i++) {

    // Mother Volume of ZDD
    G4double R  = m_R + ZDD_Volume_Thickness* 0.5;

    G4double X = R * sin(m_Theta) * cos(Phi);
    G4double Y = R * sin(m_Theta) * sin(Phi);
    G4double Z = R * cos(m_Theta);
    G4ThreeVector Det_pos = G4ThreeVector(X, Y, Z);

    // Motehr Volume de la ZDD
    G4RotationMatrix* Rot1 = new G4RotationMatrix();
    Rot1->rotateY(m_Theta);

    G4Box* MotherSolid
            = new G4Box("MotherVolume", ZDD_Volume_Width * 0.5,
                        ZDD_Volume_Length * 0.5, ZDD_Volume_Thickness* 0.5);

    G4Material* VolumeMaterial
            = MaterialManager::getInstance()->GetMaterialFromLibrary("Vacuum");
    G4LogicalVolume* MotherVolume = new G4LogicalVolume(
                MotherSolid, VolumeMaterial, "MotherVolume", 0, 0, 0);

    new G4PVPlacement(G4Transform3D(*Rot1, Det_pos), MotherVolume, "MotherVolume",
                      world, false, 0);
    //    MotherVolume->SetVisAttributes(m_Vis_ZDD);
    MotherVolume->SetVisAttributes(G4VisAttributes::GetInvisible());


    // Drift Chamber 1
    if(m_Drift_Chamber_Z.size() > 0)
    {
        //        new G4PVPlacement(0, G4ThreeVector(0, 0, (-ZDD_Volume_Thickness + m_Drift_Chamber_Thickness[0]) * 0.5 + m_Drift_Chamber_Z[0]),
        //                Build_Drift_Chamber_1(), "Drift_Chamber_1", MotherVolume, false, 1);

        //        // Drift Chamber 2
        //        new G4PVPlacement(0, G4ThreeVector(0, 0, (-ZDD_Volume_Thickness + m_Drift_Chamber_Thickness[1])* 0.5 + m_Drift_Chamber_Z[1]),
        //                Build_Drift_Chamber_2(), "Drift_Chamber_2", MotherVolume, false, 2);

        new G4PVPlacement(0, G4ThreeVector(0, 0, m_Drift_Chamber_Thickness[0] * 0.5 + m_Drift_Chamber_Z[0]),
                Build_Drift_Chamber_1(), "Drift_Chamber_1", world, false, 1);

        // Drift Chamber 2
        new G4PVPlacement(0, G4ThreeVector(0, 0, m_Drift_Chamber_Thickness[1] * 0.5 + m_Drift_Chamber_Z[1]),
                Build_Drift_Chamber_2(), "Drift_Chamber_2", world, false, 2);

    }

    // Ionisation Chambers (gas)
    if( m_Ionisation_Chamber_Z.size() > 0)
    {
        G4Box* Entrance_box = new G4Box("ZDD_Entrance_box", Mylar_Foil_Width * 0.5, Mylar_Foil_Length * 0.5,
                               Mylar_Foil_Thickness * 0.5);
        G4LogicalVolume* logic_Entrance_box
                = new G4LogicalVolume(Entrance_box, Mylar_mat, "logic_ZDD_Entrance_box", 0, 0, 0);
        logic_Entrance_box->SetVisAttributes(m_Vis_Mylar);

        G4RotationMatrix* Rot_Ionisation_Chamber = new G4RotationMatrix();
        Rot_Ionisation_Chamber->rotateX(Ionisation_Chamber_Angle);

        new G4PVPlacement(G4Transform3D(*Rot_Ionisation_Chamber,
                                        G4ThreeVector(0, 0, m_Ionisation_Chamber_Z[0]-20.*mm)),
                          logic_Entrance_box, "ZDD_Entrance_box", world, false, 1);

        for (Int_t i = 0; i < m_Ionisation_Chamber_nDet[0]; i++)
        {
            Build_Ionisation_Chamber(world,i);
        }
    }

    // AC
    if (m_AC_Z.size()>0)
    {
        for (int i = 0; i < ACcounter; i++) {
            G4RotationMatrix* Rot_AC = new G4RotationMatrix();
            Rot_AC->rotateX(m_AC_Angle[i]);
            G4Box* box = new G4Box("ZDD_AC", AC_Width * 0.5, AC_Length * 0.5,
                                   m_AC_Thickness[i] * 0.5);
            G4Material* MaterialAC = MaterialManager::getInstance()->GetMaterialFromLibrary(m_AC_Material[i]);
            G4LogicalVolume* AC
                    = new G4LogicalVolume(box, MaterialAC, "logic_ZDD_AC", 0, 0, 0);
            AC->SetVisAttributes(m_Vis_AC);
            new G4PVPlacement(G4Transform3D(*Rot_AC,
                                            G4ThreeVector(0, 0, (-ZDD_Volume_Thickness + m_AC_Thickness[i])*0.5 + m_AC_Z[i])),
                              AC, "AC", MotherVolume, false, i);
        }
    }


    // Entry_Exit
    if (m_Entry_Exit_Z.size()>0)
    {
        G4RotationMatrix* Rot_Ionisation_Chamber = new G4RotationMatrix();
        Rot_Ionisation_Chamber->rotateX(Ionisation_Chamber_Angle);
        for (int i = 0; i < Entry_Exit_counter; i++) {
            G4Box* box = new G4Box("ZDD_Entry_Exit", Mylar_Foil_Width * 0.5, Mylar_Foil_Length * 0.5,
                                   m_Entry_Exit_Thickness[i] * 0.5);
            G4Material* Material_Entry_Exit = MaterialManager::getInstance()->GetMaterialFromLibrary(m_Entry_Exit_Material[i]);
            G4LogicalVolume* Entry_Exit
                    = new G4LogicalVolume(box, Material_Entry_Exit, "logic_ZDD_Entry_Exit", 0, 0, 0);
            Entry_Exit->SetVisAttributes(m_Vis_AC);
            new G4PVPlacement(G4Transform3D(*Rot_Ionisation_Chamber,
                                            G4ThreeVector(0, 0, (-ZDD_Volume_Thickness + m_Entry_Exit_Thickness[i])*0.5 + m_Entry_Exit_Z[i])),
                              Entry_Exit, "Entry_Exit", MotherVolume, false, i);
        }

    }

    // Plastic
    if (m_Plastic_Width.size()>0 )
    {
        for (int i = 0; i < Plasticcounter; i++) {
            G4Box* box = new G4Box("ZDD_Plastic", m_Plastic_Width[i] * 0.5, m_Plastic_Length[i] * 0.5,
                                   m_Plastic_Thickness[i] * 0.5);
            G4Material* Material = MaterialManager::getInstance()->GetMaterialFromLibrary(m_Plastic_Material[i]);
            G4LogicalVolume* Plastic
                    = new G4LogicalVolume(box, Material, "logic_ZDD_Plastic", 0, 0, 0);
            Plastic->SetVisAttributes(m_Vis_Plastic);
            Plastic->SetSensitiveDetector(m_Plastic_Scorer);
            //            new G4PVPlacement(
            //                        0,
            //                        G4ThreeVector(m_Plastic_Position[i][0], m_Plastic_Position[i][1],
            //                    (-ZDD_Volume_Thickness + m_Plastic_Thickness[i])*0.5 + m_Plastic_Position[i][2]),
            //                    Plastic, "Plastic", MotherVolume, false, i);
            new G4PVPlacement(
                        0,
                        G4ThreeVector(m_Plastic_Position[i][0], m_Plastic_Position[i][1],
                    m_Plastic_Thickness[i]*0.5 + m_Plastic_Position[i][2]),
                    Plastic, "Plastic", world, false, i);
        }
    }

    /*G4double wX = m_R[i] * sin(m_Theta[i] ) * cos(m_Phi[i] ) ;
    G4double wY = m_R[i] * sin(m_Theta[i] ) * sin(m_Phi[i] ) ;
    G4double wZ = m_R[i] * cos(m_Theta[i] ) ;
    G4ThreeVector Det_pos = G4ThreeVector(wX, wY, wZ) ;
    // So the face of the detector is at R instead of the middle
    Det_pos+=Det_pos.unit()*ZDD_NS::Thickness*0.5;
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

    if(m_Shape[i] == "Cylindrical"){
      new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildCylindricalDetector(),
          "ZDD",world,false,i+1);
    }

    else if(m_Shape[i] == "Square"){
      new G4PVPlacement(G4Transform3D(*Rot,Det_pos),
          BuildSquareDetector(),
          "ZDD",world,false,i+1);
    }*/
    // }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Add Detector branch to the EventTree.
// Called After DetecorConstruction::AddDetector Method
void ZDD::InitializeRootOutput(){
    RootOutput *pAnalysis = RootOutput::getInstance();
    TTree *pTree = pAnalysis->GetTree();
    if(!pTree->FindBranch("ZDD")){
        pTree->Branch("ZDD", "TZDDData", &m_Event) ;
    }
    pTree->SetBranchAddress("ZDD", &m_Event) ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Read sensitive part and fill the Root tree.
// Called at in the EventAction::EndOfEventAvtion
void ZDD::ReadSensitive(const G4Event* ){
    m_Event->Clear();
    ///////////
    // Calorimeter scorer
    CalorimeterScorers::PS_Calorimeter* Scorer_IC
            = (CalorimeterScorers::PS_Calorimeter*)m_IC_Scorer->GetPrimitive(0);
    unsigned int size_IC = Scorer_IC->GetMult();
    for (unsigned int i = 0; i < size_IC; i++) {
        vector<unsigned int> level = Scorer_IC->GetLevel(i);
        G4double             Energy
                = RandGauss::shoot(Scorer_IC->GetEnergy(i), Scorer_IC->GetEnergy(i)*ZDD_NS::ResoEnergy_IC/(100*2.355));
        if (Energy > ZDD_NS::EnergyThreshold_IC) {
            G4double Time = RandGauss::shoot(Scorer_IC->GetTime(i), ZDD_NS::ResoTime_IC);
            int      DetectorNbr = level[0];
            //FIXME
            m_Event->SetZDDIC(DetectorNbr, Energy, Time);
        }
    }

    CalorimeterScorers::PS_Calorimeter* Scorer_Plastic
            = (CalorimeterScorers::PS_Calorimeter*)m_Plastic_Scorer->GetPrimitive(
                0);
    unsigned int size_Plastic = Scorer_Plastic->GetMult();
    for (unsigned int i = 0; i < size_Plastic; i++) {
        vector<unsigned int> level = Scorer_Plastic->GetLevel(i);
        G4double             Energy
                = RandGauss::shoot(Scorer_Plastic->GetEnergy(i), Scorer_Plastic->GetEnergy(i)*ZDD_NS::ResoEnergy_Plastic/(100*2.355));
        if (Energy > ZDD_NS::EnergyThreshold_Plastic) {
            G4double Time = RandGauss::shoot(Scorer_Plastic->GetTime(i), ZDD_NS::ResoTime_Plastic);
            int      DetectorNbr = level[0];
            m_Event->SetZDDPL(DetectorNbr,Energy,Time);
        }
    }

    ///////////
    // Calorimeter scorer
    /*CalorimeterScorers::PS_Calorimeter* Scorer= (CalorimeterScorers::PS_Calorimeter*) m_ZDDScorer->GetPrimitive(0);

  unsigned int size = Scorer->GetMult();
  for(unsigned int i = 0 ; i < size ; i++){
    vector<unsigned int> level = Scorer->GetLevel(i);
    double Energy = RandGauss::shoot(Scorer->GetEnergy(i),ZDD_NS::ResoEnergy);
    if(Energy>ZDD_NS::EnergyThreshold){
      double Time = RandGauss::shoot(Scorer->GetTime(i),ZDD_NS::ResoTime);
      int DetectorNbr = level[0];
      m_Event->SetEnergy(DetectorNbr,Energy);
      m_Event->SetTime(DetectorNbr,Time);
    }
  }*/
    DriftChamberScorers::PS_DriftChamber* Scorer_DC_1
            = (DriftChamberScorers::PS_DriftChamber*)m_Drift_Chamber_Scorer_1->GetPrimitive(0);
    unsigned int size1 = Scorer_DC_1->GetMult();

    for (unsigned int i = 0; i < size1; i++) {
        vector<unsigned int> level     = Scorer_DC_1->GetLevel(i);
        G4double               DriftTime
                = RandGauss::shoot(Scorer_DC_1->GetDriftTime(i), ZDD_NS::ResoDriftTime);
        //        printf("Ssize1 : %d \n",size1);
        //        printf("Scorer_DC_1->GetCounter(i)) : %lf \n",Scorer_DC_1->GetCounter(i));
        //        printf("DT3 : %lf \n",Scorer_DC_1->GetDriftTime(i));
        //        printf("DC : %lf \n",DriftTime);
        int DetectorNbr = level[0];
        m_Event->SetDrift(DetectorNbr,  DriftTime, 0.);
    }

    DriftChamberScorers::PS_DriftChamber* Scorer_DC_2
            = (DriftChamberScorers::PS_DriftChamber*)m_Drift_Chamber_Scorer_2->GetPrimitive(0);
    unsigned int size2 = Scorer_DC_2->GetMult();
    for (unsigned int i = 0; i < size2; i++) {
        vector<unsigned int> level     = Scorer_DC_2->GetLevel(i);
        G4double               DriftTime
                = RandGauss::shoot(Scorer_DC_2->GetDriftTime(i), ZDD_NS::ResoDriftTime);
        int DetectorNbr = level[0];
        //FIXME
        m_Event->SetDrift(DetectorNbr, DriftTime, 0.);
        //      printf("DC%d : %lf \n",DetectorNbr,Scorer_DC_1->GetDriftTime(i)*1e3);

    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////
void ZDD::InitializeScorers()
{
    // This check is necessary in case the geometry is reloaded
    bool already_exist = false;
    m_Drift_Chamber_Scorer_1  = CheckScorer("Drift_Chamber_Scorer_1", already_exist);
    m_Drift_Chamber_Scorer_2  = CheckScorer("Drift_Chamber_Scorer_2", already_exist);
    m_IC_Scorer = CheckScorer("IC_Scorer",already_exist) ;
    m_Plastic_Scorer = CheckScorer("Plastic_Scorer",already_exist) ;

    if(already_exist)
        return ;

    // Otherwise the scorer is initialised
    vector<int> level;
    level.push_back(0);

    //G4VPrimitiveScorer* Calorimeter= new CalorimeterScorers::PS_Calorimeter("Calorimeter",level, 0) ;
    //G4VPrimitiveScorer* Interaction= new InteractionScorers::PS_Interactions("Interaction",ms_InterCoord, 0) ;
    G4VPrimitiveScorer* Drift1 = new DriftChamberScorers::PS_DriftChamber("Drift", level, DriftDir1, DriftSpeed, 0);
    G4VPrimitiveScorer* Drift2 = new DriftChamberScorers::PS_DriftChamber("Drift", level, DriftDir2, DriftSpeed, 0);
    G4VPrimitiveScorer* IC = new CalorimeterScorers::PS_Calorimeter("IC", level, 0);
    G4VPrimitiveScorer* Plastic = new CalorimeterScorers::PS_Calorimeter("Plastic", level, 0);


    //and register it to the multifunctionnal detector

    //m_ZDDScorer->RegisterPrimitive(Calorimeter);
    //m_ZDDScorer->RegisterPrimitive(Interaction);
    m_Drift_Chamber_Scorer_1->RegisterPrimitive(Drift1);
    m_Drift_Chamber_Scorer_2->RegisterPrimitive(Drift2);
    m_IC_Scorer->RegisterPrimitive(IC);
    m_Plastic_Scorer->RegisterPrimitive(Plastic);

    G4SDManager::GetSDMpointer()->AddNewDetector(m_Drift_Chamber_Scorer_1) ;
    G4SDManager::GetSDMpointer()->AddNewDetector(m_Drift_Chamber_Scorer_2) ;
    G4SDManager::GetSDMpointer()->AddNewDetector(m_IC_Scorer);
    G4SDManager::GetSDMpointer()->AddNewDetector(m_Plastic_Scorer);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPS::VDetector* ZDD::Construct(){
    return  (NPS::VDetector*) new ZDD();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern"C" {
class proxy_nps_ZDD{
public:
    proxy_nps_ZDD(){
        NPS::DetectorFactory::getInstance()->AddToken("ZDD","ZDD");
        NPS::DetectorFactory::getInstance()->AddDetector("ZDD",ZDD::Construct);
    }
};

proxy_nps_ZDD p_nps_ZDD;
}
