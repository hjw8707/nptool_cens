#ifndef NeuLAND_h
#define NeuLAND_h 1
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  NeuLAND simulation                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
#include <map>
using namespace std;

// G4 headers
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TNeuLANDData.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"

class NeuLAND : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    NeuLAND() ;
    virtual ~NeuLAND() ;

    ////////////////////////////////////////////////////
    //////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddWall(G4ThreeVector POS, int NbrModule,bool veto,bool frame);
    void ReadXML(std::string xml_file,G4ThreeVector offset, bool InvertX,bool InvertY); 

    G4LogicalVolume* BuildModule();
    G4LogicalVolume* BuildVeto();
  
  private:
    std::map<unsigned int , G4ThreeVector> m_PositionBar;
    std::map<unsigned int , bool> m_IsVetoBar;
    std::map<unsigned int , bool> m_IsHorizontal;
    G4LogicalVolume* m_Module;
    G4LogicalVolume* m_Veto;
    double Energy;
    double Light;
    
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
    // Called After DetectorConstruction::AddDetector Method
    void InitializeRootOutput() ;

    // Read sensitive part and fill the Root tree.
    // Called at in the EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_ModuleScorer ;
    G4MultiFunctionalDetector* m_VetoScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TNeuLANDData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate 
    vector<G4ThreeVector>  m_Pos; 
    vector<int>  m_NbrModule;
    int m_TotalModule = 0;
    vector<bool> m_HasVeto; 
    vector<bool> m_HasFrame; 
    
    // Visualisation Attribute
    G4VisAttributes* m_VisModule;
    G4VisAttributes* m_VisVeto;
    G4VisAttributes* m_VisPMT;
    G4VisAttributes* m_VisFrame;

  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
