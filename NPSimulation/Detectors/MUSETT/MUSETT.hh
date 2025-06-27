#ifndef MUSETT_h
#define MUSETT_h 1
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr*
 *                                                                           *
 * Creation Date  : January 2024                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  MUSETT simulation for double alpha decay            *
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
#include "TMUSETTData.h"
#include "NPInputParser.h"
namespace MUSETT_NS {




  // Resolution
  const G4double SigmaTime =  0 * ns; //
  const G4double SigmaEnergy = 15 * keV;     //

  // Threshold
  const G4double EnergyThreshold = 0 * MeV;
  const double T_sum = 0.0000000001 * ns; // if ∆T < T_sum in same strip : sum events

  // Geometry
  // Detector + Dead Layer
  const G4double SiliconThickness    = 300 * micrometer;
  const G4double DeadLayerThickness  = 0.1 * micrometer;

  const G4double SquareLength  = 97220 * micrometer;

  //Target
  const G4double targetRadius = 1.2 * cm / 2 ;
  const G4double targetHalfThickness = 88 * nanometer /2;


  // Target Holder
  const G4double lowerPartZPosition = 17 * mm;
  const G4double holderWidth = 1 * mm;
  const G4double holderLength = 34 * mm;
  const G4double holderHeight = 135.47 * mm;

  //Aluminium lattice in front of dead layer
  const G4double AluHeight = 0.625 * um;

  const G4double XstickWidth = 660 * um;
  const G4double XstickLength = 110/3 * um;
  const G4double XstickHeight = AluHeight;

  const G4double YstickWidth = 20 * um;
  const G4double YstickLength = 97220 * um; // Adjusted Ystick length
  const G4double YstickHeight = AluHeight;


  const G4double BondPadSize = 300* um;
  const G4double BondPadHeight = 0.625*um;


  const G4double InterStrip = 60 * um;

  const G4int numXsticks = 96;
  const G4int numYsticks = 128;
} // namespace MUSETT_NS


class MUSETT : public NPS::VDetector{
  ////////////////////////////////////////////////////
  /////// Default Constructor and Destructor /////////
  ////////////////////////////////////////////////////
  public:
    MUSETT() ;
    virtual ~MUSETT() ;

    ////////////////////////////////////////////////////
    /////// Specific Function of this Class ///////////
    ////////////////////////////////////////////////////
  public:
    // Cartesian
    void AddDetector(int DetectorNumber,
      G4ThreeVector PX1_Y1 ,
      G4ThreeVector PX1_Y128=G4ThreeVector() ,
      G4ThreeVector PX128_Y1=G4ThreeVector(),
      G4ThreeVector PX128_Y128=G4ThreeVector());

    G4LogicalVolume* BuildSquareDetector();
    G4LogicalVolume* BuildDeadLayer();
    G4LogicalVolume* BuildTarget();
    G4LogicalVolume* BuildTargetHolder();
    G4LogicalVolume* BuildAlLattice();



  private:
    G4LogicalVolume* m_SquareDetector;
    G4LogicalVolume* m_SquareDeadLayer;

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
    G4MultiFunctionalDetector* m_SquareScorer;
    ////////////////////////////////////////////////////
    ///////////Event class to store Data////////////////
    ////////////////////////////////////////////////////
  private:
    TMUSETTData* m_Event ;

    ////////////////////////////////////////////////////
    ///////////////Private intern Data//////////////////
    ////////////////////////////////////////////////////
  private: // Geometry
    // Detector Coordinate
    // Used for "By Point Definition"
    vector<G4ThreeVector>   m_X1_Y1     ; // Top Left Corner Position Vector
    vector<G4ThreeVector>   m_X1_Y128   ; // Bottom Left Corner Position Vector
    vector<G4ThreeVector>   m_X128_Y1   ; // Bottom Right Corner Position Vector
    vector<G4ThreeVector>   m_X128_Y128 ; // Center Corner Position Vector

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

    ////////////////////////////////////////////////////
    //////////////////// Material //////////////////////
    ////////////////////////////////////////////////////
  private :
    void InitializeMaterial();

    G4Material* m_MaterialSilicon;
    G4Material* m_MaterialBoron;
    G4Material* m_MaterialCarbon;
    G4Material* m_MaterialAluminium;
    G4Material* m_MaterialVacuum;
  // Needed for dynamic loading of the library
  public:
    static NPS::VDetector* Construct();
};
#endif
