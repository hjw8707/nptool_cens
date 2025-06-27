#ifndef Vamos_h
#define Vamos_h 
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Cyril Lenain  contact address: lenain@lpccaen.in2p3.fr   *
 *                                                                           *
 * Creation Date  : octobre 2018                                             *
 * Last update    : 09/01/2019                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe The Vamos spectrometer                               *
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
#include "G4AssemblyVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4UnionSolid.hh"
#include "G4VFastSimulationModel.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TVamosData.h"
#include "NPInputParser.h"

class Vamos : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Vamos() ;
    virtual ~Vamos() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:  

    void AddVamos(G4double R,double Theta);  
    void AddBeamCatcher(string Material, G4double Width, double Length, double Thickness, G4ThreeVector Pos);
    void AddDetector(G4double Z, string Gas, double Pressure, double Temperature);       
    void AddIC(G4double Z,  double Thickness, string Gas, double Pressure, double Temperature);

    G4LogicalVolume* BuildDC1();
    G4LogicalVolume* BuildDC2();

    G4LogicalVolume* BuildQuad1();
    G4LogicalVolume* BuildQuad2();
    G4LogicalVolume* BuildDipol();

    G4LogicalVolume* BuildBeamCatcher();
    G4LogicalVolume* BuildMWPPAC();
    G4LogicalVolume* BuildDC3();
    G4LogicalVolume* BuildDC4();
    G4LogicalVolume* BuildIC(); 
    G4AssemblyVolume* BuildTMW1(); 
    G4AssemblyVolume* BuildTMW2(); 
    G4AssemblyVolume* BuildFPMW1(); 
    G4AssemblyVolume* BuildFPMW2(); 

  private:

    G4LogicalVolume* m_DC1;
    G4LogicalVolume* m_DC2;

    G4LogicalVolume* m_Quad1;
    G4LogicalVolume* m_Quad2;
    G4LogicalVolume* m_Dipol;    

    G4LogicalVolume* m_BeamCatcher;
    G4LogicalVolume* m_MWPPAC;
    G4LogicalVolume* m_DC3;
    G4LogicalVolume* m_DC4;
    G4LogicalVolume* m_IC;
    G4AssemblyVolume* m_TMW1;
    G4AssemblyVolume* m_TMW2;
    G4AssemblyVolume* m_FPMW1;
    G4AssemblyVolume* m_FPMW2;

    G4double ICcounter;

    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)

    void ClearGeometry();

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

  public: 
    // Scorer
    // Initialize all Scorer used by Vamos
    void InitializeScorers() ;

    // Associated Scorer
    G4MultiFunctionalDetector* m_CalorimeterScorer ;
    G4MultiFunctionalDetector* m_DCScorer ;
    G4MultiFunctionalDetector* m_InterScorer ;

    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////

  private:
    TVamosData* m_Event ;

  private:
    G4Region* m_ChargeStateRegion;
    vector<G4VFastSimulationModel*> m_ChargeStateModel;

    int m_StepSize;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////

  private: 

    // Geometry

    G4double R_Catcher        = 0;
    G4double CatcherWidth     = 0;
    G4double CatcherLength    = 0;
    G4double CatcherThickness = 0;
    G4double       m_R        = 0; // distance Target- Entrance of the Mother Volume
    G4double       m_Theta    = 0; 
    // Detector Coordinate  

    string CatcherMaterial;

    G4ThreeVector m_PosCatcher; 
    vector<G4double> m_Z ;
    vector<string> m_Gas;
    vector<G4double> m_Pressure;
    vector<G4double> m_Temperature;

    vector<G4double> m_ZIC;
    vector<G4double> m_ThicknessIC;
    vector<G4double> m_PressureIC;
    vector<G4double> m_TemperatureIC;
    vector<string> m_GasIC;

    // TMW1-2
    double m_Z_TMW1;
    double m_Z_TMW2;
    string m_Gas_TMW1;
    string m_Gas_TMW2;
    double m_Pressure_TMW1;
    double m_Pressure_TMW2;
    double m_Z_FPMW1;
    double m_Z_FPMW2;
    string m_Gas_FPMW1;
    string m_Gas_FPMW2;
    double m_Pressure_FPMW1;
    double m_Pressure_FPMW2;

    //   Shape type

    // Visualisation Attribute
    G4VisAttributes* m_VisQuad;
    G4VisAttributes* m_VisDC;
    G4VisAttributes* m_VisVolumeVamos;
    G4VisAttributes* m_VisCatcher;
    G4VisAttributes* m_VisGasC4H10;
    G4VisAttributes* m_VisGasCF4;
    G4VisAttributes* m_VisMylar;

  public:
    static NPS::VDetector* Construct();
};

#endif
