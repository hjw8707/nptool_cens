#ifndef SEASON_h
#define SEASON_h 1
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emmanuel Rey-herme                                       *
 * contact address: marine.vandebrouck@cea.fr                                *
 *                                                                           *
 * Creation Date  : septembre 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  SEASON simulation                             *
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
 #include "TSEASONData.h"
 #include "NPInputParser.h"

 class SEASON : public NPS::VDetector{
   ////////////////////////////////////////////////////
   /////// Default Constructor and Destructor /////////
   ////////////////////////////////////////////////////
   public:
     SEASON() ;
     virtual ~SEASON() ;

     ////////////////////////////////////////////////////
     /////// Specific Function of this Class ///////////
     ////////////////////////////////////////////////////
   public:
     // Cartesian
   void AddDetector(G4ThreeVector X1_Y1, G4ThreeVector X1_YMax, G4ThreeVector XMax_Y1, G4ThreeVector XMax_YMax);
   void AddDetector(double DistDSSD1, double DistTunnel, bool Deloc, bool EnableWheel, bool RotateWheel, bool EnableChamber);

     G4LogicalVolume* BuildSquareDetector();
     G4LogicalVolume* BuildWindow();
     G4LogicalVolume* BuildGrid();
     G4LogicalVolume* BuildWheel();
     G4LogicalVolume* BuildRotatedWheel();
     G4LogicalVolume* BuildChamber();
     G4LogicalVolume* BuildFoil();
   
   private:
     
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
     G4MultiFunctionalDetector* m_SEASONScorer ;
     ////////////////////////////////////////////////////
     ///////////Event class to store Data////////////////
     ////////////////////////////////////////////////////
   private:
     TSEASONData* m_Event ;

     ////////////////////////////////////////////////////
     ///////////////Private intern Data//////////////////
     ////////////////////////////////////////////////////
   private: // Geometry
     // Detector Coordinate
     vector<G4ThreeVector>   m_X1_Y1     ; // Top Left Corner Position Vector
     vector<G4ThreeVector>   m_X1_YMax   ; // Bottom Left Corner Position Vector
     vector<G4ThreeVector>   m_XMax_Y1   ; // Bottom Right Corner Position Vector??
     vector<G4ThreeVector>   m_XMax_YMax ; // Center Corner Position Vector??
     
     vector<G4ThreeVector> m_Center;
     vector<G4ThreeVector> m_Norm;
     vector<double> m_Rotation;
    
     bool m_EnableWheel;
     bool m_RotateWheel;
     bool m_EnableChamber;
     
     double m_DistDSSD1;
     double m_DistTunnel;
     
     double m_FoilRadius;
     double m_FoilThickness;
     
     double m_EnergyThreshold; //0.02*MeV;
     
     // Visualisation Attribute
     G4VisAttributes* m_VisSquare;
     G4VisAttributes* m_VisWheel;
     G4VisAttributes* m_VisFoil;
     G4VisAttributes* m_VisGe;
     G4VisAttributes* m_VisWindow;
     G4VisAttributes* m_VisGrid;
     G4VisAttributes* m_VisChamber;

   // Needed for dynamic loading of the library
   public:
     static NPS::VDetector* Construct();
 };
 #endif
