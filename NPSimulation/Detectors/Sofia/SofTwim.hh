#ifndef SofTwim_h
#define SofTwim_h 1
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
 *  This class describe  SofTwim simulation                             *
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
#include "TSofTwimData.h"
#include "NPInputParser.h"

class SofTwim : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    SofTwim() ;
    virtual ~SofTwim() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS);
    // Spherical
    void AddDetector(double R,double Theta,double Phi);  

    G4LogicalVolume* BuildTwinMusic();
    G4LogicalVolume* BuildTwinSection();

  private:
    G4LogicalVolume* m_TwinMusic;
    G4LogicalVolume* m_TwinSection;
    G4LogicalVolume* m_AnodeDriftArea;

    double m_Pressure;
    string m_TwimGas;
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
    G4MultiFunctionalDetector* m_TwimScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TSofTwimData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_R; 
    vector<double>  m_Theta;
    vector<double>  m_Phi; 
    
    // Visualisation Attribute
    G4VisAttributes* m_VisTwimMother;
    G4VisAttributes* m_VisSquare;
    G4VisAttributes* m_VisCathode;
    G4VisAttributes* m_VisMylar;
    G4VisAttributes* m_VisKapton;

    G4Material* m_Vacuum;
    G4Material* m_Mylar;
    G4Material* m_Kapton;
    G4Material* m_Al;
    G4Material* m_Gas;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
