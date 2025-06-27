#ifndef ZDD_h
#define ZDD_h 1
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

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TZDDData.h"
#include "NPInputParser.h"

class ZDD : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    ZDD() ;
    virtual ~ZDD() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS, string Shape);
    // Spherical
    void AddDetector(double R,double Theta,double Phi,string Shape);  
    
    void Add_Drift_Chamber(G4double Z,double Thickness, string Gas, double Pressure, double Temperature);       
    
    void Add_Ionisation_Chamber(G4int nDet, G4double Z, double Thickness, string Gas, double Pressure, double Temperature);
    
    void Add_Gas_Gap(G4double Z, double Thickness, string Gas, double Pressure, double Temperature);
    
    void Add_AC(G4double Z, double Thickness, string Material, G4double Angle);
    
    void Add_Entry_Exit(G4double Z, double Thickness, string Material);
    
    void Add_Plastic(string Material, G4double Width, double Length, double Thickness, G4ThreeVector Pos);
    
    void Add_ZDD(G4double R, double theta);  
    

    G4LogicalVolume* BuildSquareDetector();
    G4LogicalVolume* BuildCylindricalDetector();
    G4LogicalVolume* Build_Drift_Chamber_1();
    G4LogicalVolume* Build_Drift_Chamber_2();
    void Build_Ionisation_Chamber(G4LogicalVolume* world, G4int Det);
  
  private:
    G4LogicalVolume* m_Drift_Chamber_1;
    G4LogicalVolume* m_Drift_Chamber_2;

    G4Material* Mylar_mat;

    G4double ICcounter;
    G4double ACcounter;
    G4double GGcounter;
    G4double Entry_Exit_counter;
    G4double Plasticcounter;
    //G4LogicalVolume* m_SquareDetector;
    //G4LogicalVolume* m_CylindricalDetector;
    
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    
    void ClearGeometry();

    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetextorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_ZDDScorer ;
    G4MultiFunctionalDetector* m_Drift_Chamber_Scorer_2 ;
    G4MultiFunctionalDetector* m_Drift_Chamber_Scorer_1 ;
    G4MultiFunctionalDetector* m_IC_Scorer ;
    G4MultiFunctionalDetector* m_Plastic_Scorer ;

    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TZDDData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    G4double  m_R; 
    G4double  m_Theta;
    //vector<double>  m_Phi; 
    
    //   Shape type
    vector<string> m_Shape ;
   
    // Visualisation Attribute
    G4VisAttributes* m_VisSquare;
    G4VisAttributes* m_VisCylinder;
    G4VisAttributes* m_Vis_Mylar;


    // Drift Chambers Gas, pressure etc.
    vector<G4double> m_Drift_Chamber_Z;
    vector<string> m_Drift_Chamber_Gas;
    vector<G4double> m_Drift_Chamber_Pressure;
    vector<G4double> m_Drift_Chamber_Temperature;
    vector<G4double> m_Drift_Chamber_Thickness;

    // Ionisation Chambers Gas, pressure etc.
    vector<G4int> m_Ionisation_Chamber_nDet;
    vector<G4double> m_Ionisation_Chamber_Z;
    vector<G4double> m_Ionisation_Chamber_Thickness;
    vector<string> m_Ionisation_Chamber_Gas;
    vector<G4double> m_Ionisation_Chamber_Pressure;
    vector<G4double> m_Ionisation_Chamber_Temperature;
    
    // Gas Gap
    vector<G4double> m_Gas_Gap_Z;
    vector<G4double> m_Gas_Gap_Thickness;
    vector<string> m_Gas_Gap_Gas;
    vector<G4double> m_Gas_Gap_Pressure;
    vector<G4double> m_Gas_Gap_Temperature;

    // Anode Cathodes
    vector<G4double> m_AC_Z;
    vector<G4bool> m_AC_Rotation;
    vector<G4double> m_AC_Thickness;
    vector<string> m_AC_Material;
    vector<G4double> m_AC_Angle;
    
    // Entry and Exit (Kapton)
    vector<G4double> m_Entry_Exit_Z;
    vector<G4double> m_Entry_Exit_Thickness;
    vector<string> m_Entry_Exit_Material;
    
    // Plastic Pos, size etc.
    vector<G4double> m_Plastic_Width;
    vector<G4double> m_Plastic_Length;
    vector<G4double> m_Plastic_Thickness;
    vector<string> m_Plastic_Material;
    vector<G4ThreeVector> m_Plastic_Position;


    // Visualisation Attributes;
    G4VisAttributes* m_Vis_Drift_Chamber_Gas;
    G4VisAttributes* m_Vis_Ionisation_Chamber_Gas;
    G4VisAttributes* m_Vis_AC;
    G4VisAttributes* m_Vis_Plastic;
    G4VisAttributes* m_Vis_ZDD;
  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
