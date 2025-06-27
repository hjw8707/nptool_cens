#ifndef ASGARD_h
#define ASGARD_h 1
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jason Park  contact address: jcpark@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 9월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ASGARD simulation                             *
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
#include "TASGARDData.h"
#include "NPInputParser.h"

class ASGARD : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    ASGARD() ;
    virtual ~ASGARD() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:

    // Add clover at a free position in space with coordinate
  // in spherical coordinate
  // Beta are the three angles of rotation in the Clover frame 
  void AddClover(int CloverId,
		 double R,
		 double Theta,
		 double Phi,
		 double BetaX,
		 double BetaY,
		 double BetaZ);
  // void AddClover(int CloverId,
  // 		 double R,
  // 		 double Theta,
  // 		 double phi);
  
  // Return a clover in the configuration given by option (not use a the moment)
  void ConstructClover();

  // Return a modeling of a segment
  G4LogicalVolume* ConstructSegment(G4int seg);
  
  // Return a modeling of the Crystal
  G4LogicalVolume* ConstructCrystal();
  
  // Return a modeling of the Capsule
  G4LogicalVolume* ConstructCapsule();
 
  // Return a modeling of the Dewar
  G4LogicalVolume* ConstructDewar();


  
  private:
  G4LogicalVolume* m_LogicClover;

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
    G4MultiFunctionalDetector* m_ASGARDScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TASGARDData* m_ASGARDData ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate
    // Clover Position
  vector<int>    m_CloverId;
  vector<double> m_R;
  vector<double> m_Theta;
  vector<double> m_Phi;
  vector<double> m_BetaX;
  vector<double> m_BetaY;
  vector<double> m_BetaZ;
  
  TVector3 v_xyz_seg[4][9]; // centroid xyz positions of the crystal & segments for 4 crystals in a clover
    
  private:/// Visualisation Attribute:
  G4VisAttributes* BlueVisAtt;
  G4VisAttributes* GreenVisAtt;
  G4VisAttributes* RedVisAtt;
  G4VisAttributes* TestVisAtt1;
  G4VisAttributes* TestVisAtt2;
  G4VisAttributes* TestVisAtt3;
  G4VisAttributes* WhiteVisAtt;
  G4VisAttributes* TrGreyVisAtt;

   
private:
  void InitializeMaterial();
  //  void DefineSegmentPositions(); // for generic clover
  
  //   List of material
  G4Material* m_MaterialVacuum  ;
  G4Material* m_MaterialGe; 
  G4Material* m_MaterialAl;
  G4Material* m_MaterialN2;

  

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
