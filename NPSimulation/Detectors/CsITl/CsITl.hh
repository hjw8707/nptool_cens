#ifndef CsITl_h
#define CsITl_h 1
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com        *
 *                                                                           *
 * Creation Date  : 11월 2023                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  CsITl simulation                                   *
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

// NPTool header
#include "NPSVDetector.hh"
#include "TCsITlData.h"
#include "NPInputParser.h"

class CsITl : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    CsITl() ;
    virtual ~CsITl() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
  void AddDetector(G4ThreeVector Pos, G4RotationMatrix Rot,
		   G4double Horizontal, G4double Vertical, G4double CsIThickness,
		   G4double LeadThickness);
  G4LogicalVolume* BuildDetector(G4int i);

  void DefineMaterials();
    
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
    G4MultiFunctionalDetector* m_CsITlScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TCsITlData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  G4Material* m_matCsITl;
  G4Material* m_matLead;
  
private: // Geometry
  vector<G4ThreeVector>  m_Pos;   // Detector Coordinate
  vector<G4RotationMatrix>  m_Rot;   // Detector Direction (psi, theta, phi)

  vector<double> m_CsIHorizontal; 
  vector<double> m_CsIVertical;
  vector<double> m_CsIThickness ;
  
  vector<double> m_LeadThickness ;

  // Visualisation Attribute
  G4VisAttributes* m_VisCsITl;
  
  // Needed for dynamic loading of the library
public:
  static NPS::VDetector* Construct();
};
#endif
