#ifndef Samurai_h
#define Samurai_h 1
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey ANNE  contact address: anne@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Samurai BDCs simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
#include <stdlib.h>
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
#include "NPCore.h"
#include "TSamuraiBDCData.h"

class SamuraiBDC : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    SamuraiBDC() ;
    virtual ~SamuraiBDC() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
  
    // Cartezian BDC
    void AddDetector(G4ThreeVector Mag_Pos, G4ThreeVector Offset, unsigned int det);
    void ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY, unsigned int det);

    G4LogicalVolume* BuildBDC1();
    G4LogicalVolume* BuildBDC2();
    G4LogicalVolume* BuildWire();
    G4LogicalVolume* BuildHexagon();

  private:
  
    //Logical Volume

    //BDC1
    std::map<unsigned int , G4ThreeVector> m_PositionWire1;
    std::map<unsigned int , int> m_LayerNbr1;
    std::map<unsigned int , int> m_WireNbr1;
    std::map<unsigned int , bool> m_IsHorizontal1;
    std::map<unsigned int , double> m_WireAngle1; 

    //BDC2
    std::map<unsigned int , G4ThreeVector> m_PositionWire2;
    std::map<unsigned int , int> m_LayerNbr2; 
    std::map<unsigned int , int> m_WireNbr2; 
    std::map<unsigned int , bool> m_IsHorizontal2;
    std::map<unsigned int , double> m_WireAngle2;
  
    G4LogicalVolume* m_BDC1; 
    G4LogicalVolume* m_BDC2; 
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
    // Initialize the scorer(s) used by the BDC detector
    void InitializeScorers() ;

    // Associated Scorer
    G4MultiFunctionalDetector* m_WireScorerBDC ;
  
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TSamuraiBDCData* m_Event;
    ////////////////////////////////////////////////////

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private:
    //Detector coordinates
    std::map<unsigned int, G4ThreeVector> m_position;//!
    std::map<unsigned int, G4ThreeVector> m_offset;//!
    std::map<unsigned int, bool> m_invertX;//!
    std::map<unsigned int, bool> m_invertY;//!
    std::map<unsigned int, bool> m_invertD;//!

    // Visualisation Attributes
    G4VisAttributes* m_VisBDC;
    G4VisAttributes* m_VisWire;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif







