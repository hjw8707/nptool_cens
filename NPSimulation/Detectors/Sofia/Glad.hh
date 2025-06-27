#ifndef Glad_h
#define Glad_h 1
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : November 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Glad simulation                             *
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
#include "G4VFastSimulationModel.hh"
#include "G4VSolid.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "NPInputParser.h"

#include "GladFieldPropagation.hh"

class Glad : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Glad() ;
    virtual ~Glad() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddMagnet(G4ThreeVector POS, string fieldmap);

    G4LogicalVolume* BuildGLADFromSTL();
    G4LogicalVolume* BuildMagnet();

  private:
    G4LogicalVolume* m_Magnet;
    G4LogicalVolume* m_GLAD_STL;
    
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
    void InitializeScorers() ;

    // Set region were magnetic field is active
    void SetPropagationRegion();

    //   Associated Scorer
    G4MultiFunctionalDetector* m_GladScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    //TGladData* m_Event ;
    G4Region* m_PropagationRegion;
    vector<G4VFastSimulationModel*> m_PropagationModel;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    double m_X;
    double m_Y; 
    double m_Z;
    double m_R; 
    double m_Theta;
    double m_Phi; 

    // GLAD //
    double m_GLAD_X;
    double m_GLAD_Y;
    double m_GLAD_Z;
    double m_GLAD_TiltAngle;
    double m_Current;
    int m_StepSize;
    string m_FieldMapFile;
    TVector3 m_MW3_POS;
    double m_MW3_CentralTheta;
    double m_MW3_R;

    // Visualisation Attribute
    G4VisAttributes* m_VisSquare;
    G4VisAttributes* m_VisGLAD;
    G4VisAttributes* m_VisField;
    G4VisAttributes* m_VisKapton;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
