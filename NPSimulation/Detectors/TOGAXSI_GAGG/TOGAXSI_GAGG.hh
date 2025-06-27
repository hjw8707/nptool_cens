#ifndef TOGAXSI_GAGG_h
#define TOGAXSI_GAGG_h 1
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
 *  This class describe  TOGAXSI_GAGG simulation                             *
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
#include "TTOGAXSI_GAGGData.h"
#include "NPInputParser.h"

#include "BeamReaction.hh"
#include "Decay.hh"

class TOGAXSI_GAGG : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    TOGAXSI_GAGG() ;
    virtual ~TOGAXSI_GAGG() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    void AddGAGGRecoilArray(G4ThreeVector Pos, double Phi, G4ThreeVector Ref);
    void AddGAGGClusterArray(G4ThreeVector Pos, double Phi, G4ThreeVector Ref);
    void AddTarget(double R, double L, string Materialname, string CellMaterialname, double CellThickness, G4ThreeVector Pos);

    G4LogicalVolume* BuildGAGGRecoilArray();
    G4LogicalVolume* BuildGAGGClusterArray();
    G4LogicalVolume* BuildTarget(int i);
  
  private:
    G4LogicalVolume* m_RecoilArray;
    G4LogicalVolume* m_ClusterArray;
    G4LogicalVolume* m_Target;

  private:
    // Initialize material used in detector definition
    void InitializeMaterial();
    
    //List of material
    G4Material* m_MaterialGAGG;
    G4Material* m_MaterialVacuum;


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
    void InitializeScorers();

    //   Associated Scorer
    G4MultiFunctionalDetector* m_RecoilArrayScorer ;
    G4MultiFunctionalDetector* m_ClusterArrayScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TTOGAXSI_GAGGData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<G4ThreeVector>  m_Pos_RecoilArray; 
    vector<double>  m_Phi_RecoilArray;
    vector<G4ThreeVector>  m_Ref_RecoilArray; 

    vector<G4ThreeVector>  m_Pos_ClusterArray; 
    vector<double>  m_Phi_ClusterArray;
    vector<G4ThreeVector>  m_Ref_ClusterArray; 
 
    // Target Coordinate
    vector<double> m_Target_R;
    vector<double> m_Target_L;
    vector<string> m_Target_MaterialName;    
    vector<string> m_Target_CellMaterialName;    
    vector<double> m_Target_CellThickness;
    vector<G4ThreeVector> m_Target_Pos;    
 
    //Region were reaction cann occure:
    G4Region *m_ReactionRegion;
    vector<G4VFastSimulationModel*> m_ReactionModel;

    // Visualisation Attribute
    G4VisAttributes* m_VisGAGG;
    G4VisAttributes* m_VisFrame;
    G4VisAttributes* m_VisTarget;
       

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
