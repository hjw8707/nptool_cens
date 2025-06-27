#ifndef Epic_h
#define Epic_h 1
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon  contact address: audrey.chatillon@cea.fr                        *
 *                                                                           *
 * Creation Date  : décembre 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Epic simulation                             *
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
#include "G4AssemblyVolume.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TEpicData.h"
#include "NPInputParser.h"

class Epic : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Epic() ;
    virtual ~Epic() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS);

    G4AssemblyVolume* BuildEpic();
    void BuildAnode(double PosZ);
    void BuildCathode(double PosZ);
    void BuildSample(double PosZ, int indexA);

  private:
    G4AssemblyVolume* m_EpicVolume;
    
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

    void Propagate(vector<double>, vector<double>, vector<double>, double, vector<double>&);

  public:   // Scorer
    //   Initialize all Scorer 
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_EpicScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TEpicData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double> m_R; 
    vector<double> m_Theta;
    vector<double> m_Phi; 

    // Gas
    string m_GasId;
    double m_Pressure;
    G4Material* m_GasMaterial;
    double m_IonisationEnergy;
    double m_DriftStepLength;
    double m_DriftVelocity;
    double m_DriftStepTime;
    double m_DriftWindowMax;
    int    m_Gauss_nbins;
    int    m_Convo_nbins;
    vector<double> m_Gauss_Distribution;


    // Anodes and Cathodes
    int    m_SegmentedAnode; // segmented anode or full anode
    int    m_nA;   // number of anodes
    double m_Distance_AK;
    double m_InterDistance_KK;
    double m_Thickness_K;
    string m_KId;
    G4Material* m_KMaterial;
    double m_RadiusAnode;
    map<int,int> m_mapping_A; // 0-based = Interal Anode 100-based = External Anode

    // Actinide Samples
    vector<int>    m_nSamplesPerA;
    vector<string> m_SampleMaterial;
    vector<double> m_SampleThickness;

    // Visualisation Attribute
    G4VisAttributes* m_VisFCWall;
    G4VisAttributes* m_VisAl;
    G4VisAttributes* m_VisCu;
    G4VisAttributes* m_VisGasSensitive;
    G4VisAttributes* m_VisGas;
    G4VisAttributes* m_VisTi;
    G4VisAttributes* m_VisRogers4003C;
    G4VisAttributes* m_VisSample235U;
    G4VisAttributes* m_VisSample238U;
    G4VisAttributes* m_VisSamplePu;
    G4VisAttributes* m_VisSample252Cf;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
