#ifndef TOGAXSI_SI_h
#define TOGAXSI_SI_h 1
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  TOGAXSI_SI simulation                             *
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

#include "G4FastSimulationManager.hh"
#include "G4UserLimits.hh"
#include "G4VFastSimulationModel.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TTOGAXSI_SIData.h"
#include "NPInputParser.h"

#include "BeamReaction.hh"
#include "Decay.hh"

class TOGAXSI_SI : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    TOGAXSI_SI() ;
    virtual ~TOGAXSI_SI() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cylindrical Coordinates
    void AddInnerXDetector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddInnerZDetector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddOuterXDetector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddOuterZDetector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddClusterInnerDetector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddClusterX1Detector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddClusterY1Detector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddClusterX2Detector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddClusterY2Detector(double R, double Z, double Phi, G4ThreeVector Ref);
    void AddTarget(double R, double L, string Materialname, string CellMaterialname, double CellThickness, G4ThreeVector Pos);


    G4LogicalVolume* BuildInnerXDetector();
    G4LogicalVolume* BuildInnerZDetector();
    G4LogicalVolume* BuildOuterXDetector();
    G4LogicalVolume* BuildOuterZDetector();
    G4LogicalVolume* BuildClusterInnerDetector();
    G4LogicalVolume* BuildClusterX1Detector();
    G4LogicalVolume* BuildClusterY1Detector();
    G4LogicalVolume* BuildClusterX2Detector();
    G4LogicalVolume* BuildClusterY2Detector();
    G4LogicalVolume* BuildTarget(int i);
    G4LogicalVolume* BuildTargetCell(int i);
    

  private:
    G4LogicalVolume* m_InnerXDetector;
    G4LogicalVolume* m_InnerZDetector;
    G4LogicalVolume* m_OuterXDetector;
    G4LogicalVolume* m_OuterZDetector;
    G4LogicalVolume* m_ClusterInnerDetector;
    G4LogicalVolume* m_ClusterX1Detector;
    G4LogicalVolume* m_ClusterY1Detector;
    G4LogicalVolume* m_ClusterX2Detector;
    G4LogicalVolume* m_ClusterY2Detector;
    G4LogicalVolume* m_Target;
    G4LogicalVolume* m_TargetCell;

  private:
    // Initialize material used in detector definition
    void InitializeMaterial();

    // List of material
    G4Material* m_MaterialSi;
    G4Material* m_MaterialVacuum;
    G4Material* m_MaterialPCB;

    // Calculated dimensions
    double m_Active_InnerXWafer_Width;
    double m_Active_InnerXWafer_Length;
    double m_Active_InnerZWafer_Width;
    double m_Active_InnerZWafer_Length;
    double m_Active_OuterXWafer_Width;
    double m_Active_OuterXWafer_Length;
    double m_Active_OuterZWafer_Width;
    double m_Active_OuterZWafer_Length;

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

    //   Associated Scorer
    G4MultiFunctionalDetector* m_InnerXScorer;
    G4MultiFunctionalDetector* m_InnerZScorer;
    G4MultiFunctionalDetector* m_OuterXScorer1;
    G4MultiFunctionalDetector* m_OuterXScorer2;
    G4MultiFunctionalDetector* m_OuterZScorer1;
    G4MultiFunctionalDetector* m_OuterZScorer2;
    G4MultiFunctionalDetector* m_ClusterInnerScorer;
    G4MultiFunctionalDetector* m_ClusterX1Scorer1;
    G4MultiFunctionalDetector* m_ClusterX1Scorer2;
    G4MultiFunctionalDetector* m_ClusterY1Scorer1;
    G4MultiFunctionalDetector* m_ClusterY1Scorer2;
    G4MultiFunctionalDetector* m_ClusterX2Scorer1;
    G4MultiFunctionalDetector* m_ClusterX2Scorer2;
    G4MultiFunctionalDetector* m_ClusterY2Scorer1;
    G4MultiFunctionalDetector* m_ClusterY2Scorer2;
    G4MultiFunctionalDetector* m_TargetCellScorer;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TTOGAXSI_SIData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_InnerX_R; 
    vector<double>  m_InnerX_Z;
    vector<double>  m_InnerX_Phi; 
    vector<G4ThreeVector> m_InnerX_Ref;
    vector<double>  m_InnerZ_R; 
    vector<double>  m_InnerZ_Z;
    vector<double>  m_InnerZ_Phi; 
    vector<G4ThreeVector> m_InnerZ_Ref;
  
    vector<double>  m_OuterX_R; 
    vector<double>  m_OuterX_Z;
    vector<double>  m_OuterX_Phi; 
    vector<G4ThreeVector> m_OuterX_Ref;
    vector<double>  m_OuterZ_R; 
    vector<double>  m_OuterZ_Z;
    vector<double>  m_OuterZ_Phi; 
    vector<G4ThreeVector> m_OuterZ_Ref;

    vector<double>  m_ClusterInner_R; 
    vector<double>  m_ClusterInner_Z;
    vector<double>  m_ClusterInner_Phi; 
    vector<G4ThreeVector> m_ClusterInner_Ref;
  
    vector<double>  m_ClusterX1_R; 
    vector<double>  m_ClusterX1_Z;
    vector<double>  m_ClusterX1_Phi; 
    vector<G4ThreeVector> m_ClusterX1_Ref;

    vector<double>  m_ClusterY1_R; 
    vector<double>  m_ClusterY1_Z;
    vector<double>  m_ClusterY1_Phi; 
    vector<G4ThreeVector> m_ClusterY1_Ref;

    vector<double>  m_ClusterX2_R; 
    vector<double>  m_ClusterX2_Z;
    vector<double>  m_ClusterX2_Phi; 
    vector<G4ThreeVector> m_ClusterX2_Ref;

    vector<double>  m_ClusterY2_R; 
    vector<double>  m_ClusterY2_Z;
    vector<double>  m_ClusterY2_Phi; 
    vector<G4ThreeVector> m_ClusterY2_Ref;

    // Target Coordinate
    vector<double> m_Target_R;
    vector<double> m_Target_L;
    vector<string> m_Target_MaterialName;
    vector<string> m_Target_CellMaterialName;
    vector<double> m_Target_CellThickness;
    vector<G4ThreeVector> m_Target_Pos;

    //Region were reaction can occure:
    G4Region *m_ReactionRegion;
    vector<G4VFastSimulationModel*> m_ReactionModel;


    // Visualisation Attribute
    G4VisAttributes* m_VisSi;
    G4VisAttributes* m_VisPCB;
    G4VisAttributes* m_VisTarget;
    G4VisAttributes* m_VisTargetCell;
     

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
