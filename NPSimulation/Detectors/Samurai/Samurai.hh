#ifndef Samurai_h
#define Samurai_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elia Pilotto  contact address: pilottoelia@gmail.com     *
 *                                                                           *
 * Creation Date  : septembre 2021                                           *
 * Last update    : septembre 2021                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Samurai simulation                                  *
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
#include "G4VFastSimulationModel.hh"
#include "G4VSolid.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "NPInputParser.h"

#include "SamuraiFieldPropagation.hh"
#include "TSamuraiIdealData.h"


class Samurai : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Samurai() ;
    virtual ~Samurai() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
  
    // Cartesian
    void AddMagnet(G4ThreeVector POS, double Angle, int method,
		 string fieldmap, double n_steps);
    // Spherical
    void AddMagnet(double R, double Theta, double Phi, double Angle,
		   int method, string fieldmap, double n_steps);

    G4VSolid* BuildRYokeSolid();

    G4LogicalVolume* BuildMagnet();
    G4LogicalVolume* BuildYoke();
    G4LogicalVolume* BuildRYoke();
    G4LogicalVolume* BuildPropvol();
  
  private:
  
    //Solids used for the magnet construction
    G4VSolid* m_RYokeSolid;

    G4LogicalVolume* m_Magnet;
    G4LogicalVolume* m_Yoke;
    G4LogicalVolume* m_RYoke;
    G4LogicalVolume* m_Propvol;

    // Visualisation Attributes
    G4VisAttributes* m_VisMagnet;
    G4VisAttributes* m_VisYokes;
    G4VisAttributes* m_VisRYokes;
    G4VisAttributes* m_VisPropvol;
    
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetectorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Construct Magnet and inialise sensitive part.
    // Called After DetectorConstruction::AddDetector Method
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Magnet branch to the EventTree.
    // Called After DetectorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAction
    void ReadSensitive(const G4Event* event) ;

  public:
    // Scorer
    // Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    // Set region were magnetic field is active:
    void SetPropagationRegion();

    // Associated Scorer
    G4MultiFunctionalDetector* m_SamuraiScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    // Store data
    TSamuraiIdealData* m_Event;

    // Region were magnetic field is active:
    G4Region* m_PropagationRegion;
    vector<G4VFastSimulationModel*> m_PropagationModel;
  

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private:

    // Detector Coordinates
    double m_R;
    double m_Theta;
    double m_Phi;

    // Angle of Rotation
    double m_Angle;

    // Propagation Parameters
    NPS::PropagationMethod m_Method;
    double m_StepSize;
    string m_FieldMapFile;
    

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif







