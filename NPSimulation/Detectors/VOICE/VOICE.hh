#ifndef VOICE_h
#define VOICE_h 1
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sunghan Bae  contact address: shbae2703@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  VOICE simulation                             *
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
#include "G4SubtractionSolid.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Material.hh"
#include "G4AssemblyVolume.hh"
#include "G4VFastSimulationModel.hh"
#include "G4FastSimulationManager.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"
#include "G4UserLimits.hh"
#include "G4PVPlacement.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TVOICEData.h"
#include "NPInputParser.h"
#include "MEventReduced.h"
#include "Decay.hh"
#include "BeamReaction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace VOICE_NS{
  // Energy and time Resolution
  const double EnergyThreshold = 1*keV;
  const double ResoTime = 0.1*ns ;
  const double ResoEnergy = 1.0*keV ;
  const double Radius = 5*micrometer ; 
  const double FGRadius = 5*micrometer ; 
  const double WireLength = 80*mm ;
  const double Thickness = 300*mm ;
  const double endthick = 5*mm ;
  const string wireMaterial = "G4_W";
  const string PCBMaterial = "PCB";

//  string detectorMaterial = "CF4";
  const G4double GasX = 100*mm;
  const G4double GasY = 100*mm;
  const G4double GasZ = 300*mm;
  const G4double ElectrodePCBX = 100*mm;
  const G4double ElectrodePCBY = 100*mm;
  const G4double ElectrodePCBZ = 1.6*mm;
  const G4double ElectrodeSubX = 80*mm;
  const G4double ElectrodeSubY = 80*mm;
  const G4double ElectrodeSubZ = 1.62*mm;

}

class VOICE : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    VOICE() ;
    virtual ~VOICE() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    void AddDetector(string Type, double zpos, double Tilt, double Rot, double radi, double pitch);
    // Cartesian
    //void AddDetector(G4ThreeVector POS, string Shape);
    // Spherical
    //void AddDetector(double R,double Theta,double Phi,string Shape);  
//    G4LogicalVolume* BuildAnodeDetector();
//    G4LogicalVolume* BuildCathodeDetector();
//    G4LogicalVolume* BuildFGDetector();
    G4AssemblyVolume* BuildAnodeDetector(double radi, double pitch);
    G4AssemblyVolume* BuildCathodeDetector(double radi, double pitch);
    G4AssemblyVolume* BuildFGDetector(double radi, double pitch);
    G4AssemblyVolume* BuildEndDetector();
//    G4LogicalVolume* BuildPCB();
    G4Material* SetGasType(string gas, double pressure);
    
    G4MultiFunctionalDetector* m_AnodeDetector;
    G4MultiFunctionalDetector* m_CathodeDetector;
    G4MultiFunctionalDetector* m_FGDetector;
    G4MultiFunctionalDetector* m_EndDetector;
    G4MultiFunctionalDetector* m_GasDetector;

    
  private:
  G4double total_length;
  G4double zoffset;

    G4AssemblyVolume* m_Anode = nullptr;
    G4AssemblyVolume* m_Cathode = nullptr;
    G4AssemblyVolume* m_FG = nullptr;
    G4AssemblyVolume* m_End = nullptr;

    G4LogicalVolume* m_logicGas;
    G4LogicalVolume* m_logicStage;
    G4Region* m_ReactionRegion;
    vector<G4VFastSimulationModel*> m_ReactionModel;
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
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

    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TVOICEData* m_Event ;
//    TVOICEGasData* m_EventGas ;

    G4int HCID_Anode;
    G4int HCID_Cathode;
    G4int HCID_FG;
    G4int HCID_END;
    G4int HCID_Gas;
    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<string>  m_Type; 
    vector<double>  m_zpos; 
    vector<double>  m_cat_zpos; 
    vector<double>  m_Tilt;
    vector<double>  m_Rot; 
    vector<double>  m_radi; 
    vector<double>  m_pitch; 
    double m_GasStepSize;
    double m_AFG_gap;
    double m_AC_gap;
    
    //   Shape type
    vector<string> m_Shape ;
   
    // Default gas type
    double m_pressure = 300./760.*atmosphere;
    string m_gas = "CF4";
    double maxZ=0.0;
    string m_mulgas[10];
    double mass_fraction[10];
    int n_gas=0;
    double mulgas_density=0.0*g/cm3;

    // Visualisation Attribute
    G4VisAttributes* m_VisGas;
    G4VisAttributes* m_VisWire;
    G4VisAttributes* m_VisPCB;
    G4VisAttributes* m_VisPCBCa;
    G4VisAttributes* m_VisPCBFG;
    G4VisAttributes* m_VisEND;
    G4VisAttributes* m_VisSquare;
   // G4VisAttributes* m_VisCylinder;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
