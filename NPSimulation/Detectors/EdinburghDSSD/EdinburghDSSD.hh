#ifndef EdinburghDSSD_h
#define EdinburghDSSD_h 1
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : mars 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  EdinburghDSSD simulation                             *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// C++ header
#include <string>
#include <vector>
using namespace std;

// G4 headersTEdinburghDSSDData
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "NPSVDetector.hh"
#include "TEdinburghDSSDData.h"
#include "NPInputParser.h"
 namespace EDIN_NS{
   // Resolution
   const G4double SigmaTime =  0 * ns; //
   const G4double SigmaEnergy = 5 * keV;     //

   // Threshold
   const G4double EnergyThreshold = 0 * MeV;
   const double T_sum = 0.0000000001 * ns; // if ∆T < T_sum in same strip : sum events

   // Geometry
   // Detector + Dead Layer
   const G4double SiliconThickness    = 1 * mm;
   const G4double SquareLength  = 4.95 * cm;


   const G4double DeadLayerThickness  = 0.8 * micrometer;


   //Foil (=carbon target)
   const double FoilThickness = 100. * nm;
   const double FoilRadius = 4.5 * mm; //0.5 mm in addition is covered by holder.
   const double FoilAngle = 35. * deg;

   const double InnerRadius = 4.5 * mm ;
   const double OuterRadius = 6.4 * mm;
   //Foil Holder
   const double FoilHolderThickness = 0.3 * mm ;
   const double FoilHolderInnerRadius = InnerRadius;
   const double FoilHolderOuterRadius = OuterRadius;

   //First Sledge
   const double FirstSledgeThickness = 0.37 * mm ;
   const double FirstSledgeInnerRadius = InnerRadius;
   const double FirstSledgeOuterRadius = OuterRadius;

   //Second Sledge
   const double SecondSledgeThickness = 0.73 * mm ;
   const double SecondSledgeInnerRadius = InnerRadius;
   const double SecondSledgeOuterRadius = OuterRadius;

   //Filler
   const double FillerThickness = FoilThickness ;
   const double FillerInnerRadius = InnerRadius;
   const double FillerOuterRadius = OuterRadius;


   //Bar
   const double BarLength = 50.0 * mm;
   const double BarWidth = 1.1 * mm;
   const double BarHeight = 14.5 * mm;
 }

class EdinburghDSSD : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    EdinburghDSSD() ;
    virtual ~EdinburghDSSD() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(int DetectorNumber,
      G4ThreeVector PX1_Y1 ,
      G4ThreeVector PX1_Y16=G4ThreeVector() ,
      G4ThreeVector PX16_Y1=G4ThreeVector(),
      G4ThreeVector PX16_Y16=G4ThreeVector());

    G4LogicalVolume* BuildSubWorld();
    G4LogicalVolume* BuildSquareDetector();
    G4LogicalVolume* BuildDeadLayer();
    G4LogicalVolume* BuildFoil();
    G4LogicalVolume* BuildNozzle();
    G4LogicalVolume* BuildFirstSledge();
    G4LogicalVolume* BuildSecondSledge();
    G4LogicalVolume* BuildFiller();
    G4LogicalVolume* BuildFoilHolder();
    G4LogicalVolume* BuildBarHoler();
    //G4LogicalVolume* BuildBar();


  private:
    G4LogicalVolume* m_SquareDetector;
    G4LogicalVolume* m_SquareDeadLayer;
    G4LogicalVolume* m_Foil;
    G4LogicalVolume* m_FoilHolder;
    G4LogicalVolume* m_FirstSledge;
    G4LogicalVolume* m_SecondSledge;
    G4LogicalVolume* m_Filler;
    G4LogicalVolume* m_subWorld;
    G4LogicalVolume* m_substract;
    G4LogicalVolume* m_bar;

    G4Tubs* solidTarget;

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
    void ProcessStrip(std::map<unsigned int, std::pair<double, double>>& map, int key,
                      double energy, double time);
    void UpdateMap(bool X_Y, std::map<unsigned int, std::pair<double, double>>& map,
                         std::map<unsigned int, bool>& mapInterstrip, int key, double energy,
                        double time, double t0,bool interstrip);
    void ReadSensitive(const G4Event* event) ;

  public:   // Scorer
    //   Initialize all Scorer used by the MUST2Array
    void InitializeScorers() ;

    //   Associated Scorer
    G4MultiFunctionalDetector* m_SquareScorer ;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
    private:
      TEdinburghDSSDData* m_Event;
      // Geometry
        // Detector Coordinate
        // Used for "By Point Definition"
        vector<G4ThreeVector>   m_X1_Y1     ; // Top Left Corner Position Vector
        vector<G4ThreeVector>   m_X1_Y16   ; // Bottom Left Corner Position Vector
        vector<G4ThreeVector>   m_X16_Y1   ; // Bottom Right Corner Position Vector
        vector<G4ThreeVector>   m_X16_Y16 ; // Center Corner Position Vector

        //   Shape type
        vector<string> m_Shape ;
        // DetectorNumber
        vector<int>    m_DetectorNumber;
        // Visualisation Attribute
        G4VisAttributes* m_VisSquare;
        G4VisAttributes* m_VisDeadLayer;
        /////// Default Constructor and Destructor /////////
        std::map<unsigned int, unsigned int> fMUMU_MapX;//!   // Pour éviter d'écirre dans l'abre ROOT
        std::map<unsigned int, unsigned int> fMUMU_MapY;//!

      // Needed for dynamic loading of the library
      public:
        static NPS::VDetector* Construct();
    };
#endif
