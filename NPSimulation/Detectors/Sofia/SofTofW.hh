#ifndef SofTofW_h
#define SofTofW_h 1
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
 *  This class describe  SofTofW simulation                             *
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
#include "TSofTofWData.h"
#include "NPInputParser.h"

class SofTofW : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    SofTofW() ;
    virtual ~SofTofW() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(G4ThreeVector POS);
    // Spherical
    void AddDetector(double R,double Theta,double Phi);  


    G4AssemblyVolume* BuildTOFDetector();
    G4AssemblyVolume* BuildMagneticField();
    G4AssemblyVolume* BuildVacuumPipe();
    G4LogicalVolume* BuildGLADFromSTL();

  private:
    G4LogicalVolume* m_PlasticTof;
    G4LogicalVolume* m_GLAD_STL;
    G4AssemblyVolume* m_GLAD;
    G4AssemblyVolume* m_TofWall;
    G4AssemblyVolume* m_VacuumPipe;
    
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
    G4MultiFunctionalDetector* m_TofScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TSofTofWData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<double>  m_R; 
    vector<double>  m_Theta;
    vector<double>  m_Phi; 
    
    // GLAD //
    int m_Build_GLAD;
    int m_Build_VacuumPipe;
    int m_Build_MagneticField;
    double m_GLAD_MagField;
    double m_GLAD_DistanceFromTarget;
    double m_GLAD_Yoffset;
    double m_GLAD_TiltAngle;
    double m_VacuumPipeX;
    double m_VacuumPipeY;
    double m_VacuumPipeZ;

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
