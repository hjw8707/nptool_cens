#ifndef AnnularS4_h
#define AnnularS4_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    : 11/10/09                                                 *
 *---------------------------------------------------------------------------*
 * Decription: Define the AnnularS4 detector from Micron                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *  + 11/10/09: Change scorer philosophy, one scorer for the detector number *
 *              added (N. de Sereville)                                      *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <vector>
using namespace std;
using namespace CLHEP;

// NPTool header
#include "NPSVDetector.hh"
#include "NPInputParser.h"
// NPTool - ROOT headers
#include "TAnnularS4Data.h"

// Geant4
#include "G4MultiFunctionalDetector.hh"
#include "G4LogicalVolume.hh"

class AnnularS4 : public NPS::VDetector{
   ////////////////////////////////////////////////////
   /////// Default Constructor and Destructor /////////
   ////////////////////////////////////////////////////
public:
   AnnularS4();
   virtual ~AnnularS4();

   ////////////////////////////////////////////////////
   //////// Specific Function of this Class ///////////
   ////////////////////////////////////////////////////
public:
   // By Position Method
   void AddModule(G4double PosZ);

   // Produce the logical volume of the detector
   G4LogicalVolume* ConstructVolume();

private:
  G4LogicalVolume* m_LogicalDetector;

   ////////////////////////////////////////////////////
   ////  Inherite from GaspardTrackerModule class /////
   ////////////////////////////////////////////////////
public:
   // Read stream at Configfile to pick-up parameters of detector (Position,...)
   // Called in DetecorConstruction::ReadDetextorConfiguration Method
   void ReadConfiguration(NPL::InputParser);

   // Construct detector and inialise sensitive part.
   // Called After DetecorConstruction::AddDetector Method
   void ConstructDetector(G4LogicalVolume* world);

   // Add Detector branch to the EventTree.
   // Called After DetecorConstruction::AddDetector Method
   void InitializeRootOutput();

   // Initialize all scorers necessary for the detector
   void InitializeScorers();

   // Read sensitive part and fill the Root tree.
   // Called at in the EventAction::EndOfEventAvtion
   void ReadSensitive(const G4Event* event);


   ////////////////////////////////////////////////////
   ///////////Event class to store Data////////////////
   ////////////////////////////////////////////////////
private:
   TAnnularS4Data* m_Event;


private:
   // Scorer
   G4MultiFunctionalDetector* m_Scorer;


   ////////////////////////////////////////////////////
   ///////////////Private intern Data//////////////////
   ////////////////////////////////////////////////////
private:
   // Used for "By Point Definition"
   vector<G4double>   m_PosZ;

  public:
    static NPS::VDetector* Construct();



};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace ANNULARS4{
   // Energy/Time resolutions for the different layers
   const G4double ResoTime        = 0.106382979    ;// = 250 ps // unchecked                //   Unit is  ns/2.35
   const G4double EnergyThreshold = 0.300          ;// = 300 keV // unchecked
   const G4double ResoEnergy      = 0.0149         ; // 0.0217 // unchecked

   // PCB
   const G4double PCBPointsX[4]={-70,70,70,-70}; // bricoled
   const G4double PCBPointsY[4]={70,70,-70,-70}; // bricoled
   const G4double PCBThickness=3.2*mm; // unchecked
   const G4double PCBInnerRadius=0*mm; // unchecked
   // Wafer
   const G4double WaferOutterRadius = 65.05*mm;
   const G4double WaferInnerRadius  = 5*mm;
   const G4double WaferThickness    = 500*micrometer; // unchecked
   const G4double WaferRCut         = 45.5*mm;  // unchecked
   const G4double ActiveWaferOutterRadius = 62.49*mm;
   const G4double ActiveWaferInnerRadius  = 7.5*mm;
   const G4double AluStripThickness   = 0.3*micrometer;
   const G4double ActiveWaferRCut     = 44.5*mm; // unchecked

   // Characteristics
   const G4int NbrRingStrips  = 1; // unsure
   const G4int NbrSectorStrips = 128; 
   const G4int NbQuadrant      = 1; // unsure

   // total WaferLength
   const G4double Length = AluStripThickness*2+WaferThickness;

   // Starting at the front and going in direction of third stage
   const G4double AluStripFront_PosZ = Length* -0.5       + 0.5*AluStripThickness; // unchecked
   const G4double Silicon_PosZ       = AluStripFront_PosZ + 0.5*AluStripThickness   + 0.5*WaferThickness; // unchecked
   const G4double AluStripBack_PosZ  = Silicon_PosZ       + 0.5*WaferThickness + 0.5*AluStripThickness; // unchecked
  
}

#endif
