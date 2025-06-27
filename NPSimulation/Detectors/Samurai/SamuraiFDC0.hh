#ifndef Samurai_h
#define Samurai_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Anne  contact address: anne@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : january 2024                                             *
 * Last update    : june 2024                                                *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  SamuraiFDC0 simulation                                  *
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
#include "G4VSolid.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "NPInputParser.h"
#include "NPXmlParser.h"

#include "TSamuraiFDC0Data.h"

class SamuraiFDC0 : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    SamuraiFDC0() ;
    virtual ~SamuraiFDC0() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
  
    // Cartezian FDC0
    void AddDetector(G4ThreeVector Mag_Pos, G4ThreeVector Offset);
    void ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY);

    G4LogicalVolume* BuildFDC0();
    G4LogicalVolume* BuildWire();
    G4LogicalVolume* BuildHexagon();
  
  private:
  
    //Logical Volume
    std::map<unsigned int , G4ThreeVector> m_PositionWire; 
    std::map<unsigned int , int> m_LayerNbr;
    std::map<unsigned int , int> m_WireNbr;
    std::map<unsigned int , double> m_WireAngle;
    std::map<unsigned int , bool> m_IsHorizontal;
  
    G4LogicalVolume* m_FDC0;
    G4LogicalVolume* m_Wire;
    G4LogicalVolume* m_Hexagon;
    
    ////////////////////////////////////////////////////
    //////  Inherite from NPS::VDetector class /////////
    ////////////////////////////////////////////////////
  public:
    // Read stream at Configfile to pick-up parameters of detector (Position,...)
    // Called in DetecorConstruction::ReadDetectorConfiguration Method
    void ReadConfiguration(NPL::InputParser) ;

    // Construct detector and initialise sensitive part.
    // (Called After DetecorConstruction::AddDetector Method)
    void ConstructDetector(G4LogicalVolume* world) ;

    // Add Detector branch to the EventTree.
    // (Called After DetecorConstruction::AddDetector Method)
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // (Called at in the EventAction::EndOfEventAction)
    void ReadSensitive(const G4Event* event) ;

  public:
    // Scorer
    // Initialize the scorer(s) used by the FDC0 detector
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_WireScorerFDC0;
  
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
   
    TSamuraiFDC0Data* m_Event;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private:
    //Detector coordinates
    G4ThreeVector m_Pos;

    // Visualisation Attributes
    G4VisAttributes* m_VisFDC0;
    G4VisAttributes* m_VisWire;
  
  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif







