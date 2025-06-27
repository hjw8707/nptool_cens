#ifndef Vendeta_h
#define Vendeta_h 1
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Vendeta simulation                             *
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

// NPTool header
#include "NPSVDetector.hh"
#include "TVendetaData.h"
#include "NPInputParser.h"

class Vendeta : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    Vendeta() ;
    virtual ~Vendeta() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS);
    // Spherical
    void AddDetector(double R,double Theta,double Phi);  


    G4AssemblyVolume* BuildVendetaDetector();
    G4LogicalVolume* BuildSensitiveCell();
    G4AssemblyVolume* BuildMecanicalStructure();

  
  private:
    
    G4AssemblyVolume* m_VendetaDetector;
    G4LogicalVolume* m_SensitiveCell;
    G4AssemblyVolume* m_MecanicalStructure;
    G4LogicalVolume* m_MecanicalStructure_Al;
    G4LogicalVolume* m_MecanicalStructure_Steel;
    G4LogicalVolume* m_LeadShield;

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
    G4MultiFunctionalDetector* m_VendetaScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TVendetaData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_R; 
    vector<double>  m_Theta;
    vector<double>  m_Phi; 
    bool m_BuildLeadShield;

    int m_Build_MecanicalStructure;

    // Visualisation Attribute
    G4VisAttributes* m_VisAl;
    G4VisAttributes* m_VisInox;
    G4VisAttributes* m_VisEJ309;
    G4VisAttributes* m_VisMuMetal;
    G4VisAttributes* m_VisPyrex;
    G4VisAttributes* m_VisEJ560;
    G4VisAttributes* m_VisLeadShield;

    G4Material* m_Vacuum;
    G4Material* m_Al;
    G4Material* m_Inox;
    G4Material* m_EJ309;
    G4Material* m_EJ560;
    G4Material* m_Pyrex;
    G4Material* m_MuMetal;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
