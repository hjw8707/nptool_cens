#ifndef Plunger_h
#define Plunger_h 1
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Plunger simulation                             *
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
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4Region.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VFastSimulationModel.hh"

// NPTool header
#include "NPInputParser.h"
#include "NPSVDetector.hh"
#include "TPlungerData.h"

class Plunger : public NPS::VDetector {
    ////////////////////////////////////////////////////
    /////// Default Constructor and Destructor /////////
    ////////////////////////////////////////////////////
   public:
    Plunger();
    virtual ~Plunger();

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
   public:
    // Build the target, stopper, and chamber
    G4LogicalVolume* BuildTarget();
    G4LogicalVolume* BuildStopper();
    G4LogicalVolume* BuildChamber();

   private:
    G4LogicalVolume* m_TargetLV;
    G4LogicalVolume* m_StopperLV;
    G4LogicalVolume* m_ChamberLV;
    G4LogicalVolume* m_ReactionRegionLV;
    G4LogicalVolume* m_DecayRegionLV;

    G4Region* m_ReactionRegion;
    G4Region* m_DecayRegion;
    vector<G4VFastSimulationModel*> m_ReactionModel;

    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
   public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetextorConfiguration Method
    void ReadConfiguration(NPL::InputParser);

    // Construct detector and inialise sensitive part.
    // Called After DetecorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world);

    // Add Detector branch to the EventTree.
    // Called After DetecorConstruction::AddDetector Method
    void InitializeRootOutput();

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event);

    // Set the reaction region
    void SetReactionRegion(G4LogicalVolume* world);

   public:  // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers();

    //   Associated Scorer (Not used)
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
   private:
    TPlungerData* m_Event;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
   private:  // Geometry
    // Target
    bool m_TargetFound;
    double m_TargetR;
    double m_TargetThickness;
    double m_TargetPosZ;
    G4String m_TargetMaterial;

    // Stopper
    bool m_StopperFound;
    double m_StopperR;
    double m_StopperThickness;
    double m_StopperPosZ;
    G4String m_StopperMaterial;

    // Chamber
    bool m_ChamberFound;
    double m_ChamberR;
    double m_ChamberThickness;
    G4String m_ChamberMaterial;
    double m_ChamberPipeR;
    double m_ChamberPipeZ0;
    double m_ChamberPipeZ1;

    // Visualisation Attribute
    G4VisAttributes* m_VisTarget;
    G4VisAttributes* m_VisStopper;
    G4VisAttributes* m_VisChamber;

    // Needed for dynamic loading of the library
   public:
    static NPS::VDetector* Construct();
};
#endif
