#ifndef Tina2_h
#define Tina2_h 1
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describes Tina simulation                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include <string>
#include <vector>
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4AssemblyVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "NPSVDetector.hh"
#include "NPInputParser.h"
#include "TTinaData.h"

namespace TINA
{
  const G4double EnergyThreshold = 0.1*MeV;
    
  // Resolution
  // sf warning: these values were copied from MUST2Array.hh 26.2.2019
  // resolution SiLi in Must2 is 0.028 but here we take it equal to CsI
  const G4double ResoTime      = 0.213;
  const G4double ResoTTT       = 0.015;
  const G4double ResoPad       = 0.080;
  const G4double ResoYY1       = 0.015;
  const G4double ResoCsI       = 0.080;
    
  // Geometry
  // easm file gives 101.0 x 100.5 but micron data file 100.42 x 100.42
  const G4double XwidthTTT     = 100.42*mm; //100.42
  const G4double YwidthTTT     = 100.42*mm;
  const G4double XwidthActiveTTT     = 97.22*mm; //100.42 // needed for scale factor of PS_Image.
  const G4double YwidthActiveTTT     = 97.22*mm;
  const G4double ZwidthTTT     =   0.3*mm;
  const G4int Nstrips     =   128;
    
  // dimensions taken as 52x55 visible area while active detector is 50x50 within frame of 56x56
  // ZoffsetPad is distance between front face TTT and front face Pad and not fixed
  /* const G4double XwidthPad     =  55.0*mm;
     const G4double YwidthPad     =  52.0*mm;
     const G4double ZwidthPad     =  20.0*mm;
     const G4double XdistbetwPads =   2.0*mm;
     const G4double YdistbetwPads =   7.0*mm;
     const G4double ZoffsetPad    =  20.0*mm;
     const G4double XwidthBox     = 119.0*mm;
     const G4double YwidthBox     = 128.0*mm;
     const G4double ZwidthBox     =  21.0*mm;//*/
    
  const G4double XwidthPad     =  49.0*mm; //50 but cut by structure
  const G4double YwidthPad     =  47.0*mm; //50 but cut by structure
  const G4double ZwidthPad     =  25.0*mm;
  const G4double XdistbetwPads =   1.0*mm;
  const G4double YdistbetwPads =   3.0*mm; // because of central bar
  const G4double ZoffsetPad    =  11.9*mm;
  const G4double XwidthBox     = 110.0*mm;
  const G4double YwidthBox     = 110.0*mm;
  const G4double ZwidthBox     =  30.0*mm;//*/ //40

  //  const G4double XwidthTCsI    =  
      
  // dimensions taken as 50x50 and 50x17.5 full volume while active detector is 45x45 and 13.5x45.5 within that
  const G4double XwidthT1A     =  50.0*mm;
  const G4double YwidthT1A     =  50.0*mm;
  const G4double ZwidthT1A     =  25.0*mm;
  const G4double FwidthT1B     =  50.0*mm; // front plane
  const G4double BwidthT1B     =  30.0*mm; // back plane
  const G4double ZwidthT1B     =  30.0*mm;
  const G4double XwidthT3      =  50.0*mm;
  const G4double YwidthT3      =  17.5*mm;
  const G4double ZwidthT3      =  17.5*mm;
  const G4double SolidWidth   =    1*mm;//40
  const G4double SolidHeight   =  85*mm;
  const G4double SolidBaseLarge=    100*mm;
  const G4double SolidBaseSmall=    40*mm;
  const G4double CsIWidth     =    25*mm;
  const G4double CsIHeight     =    80*mm;
  const G4double CsIBaseLarge  =    30*mm;
  const G4double CsIBaseSmall  =    12*mm;
  // slight asymmetry in phi angular placement of strips within frame ignored
  // frame thickness reduced from 15.2 to 5.2 to avoid overlap with Type1A therefore YY1 floats in front of frame
  const G4double Nfrontrings   =    16   ;
  const G4double Nbacksectors  =     1   ;
  const G4double RinnerYY1     =    55*mm; // approximate 49
  const G4double RouterYY1     = 130.0*mm; // approximate 131
  const G4double ZwidthYY1     =   0.3*mm;
  const G4double phiYY1        = -110*degree;//68.81*degree; // 1.31+67.50=68.81°
  const G4double alphaYY1      = 40*degree; // 1.31+42.38+1.31=45°
  const G4double YoffsetYY1    = 108.3*mm;
  const G4double ZoffsetYY1    =   6.7*mm;
  const G4double RinnerFrame   =  39.0*mm;
  const G4double RouterFrame   = 144.0*mm;
  const G4double ZwidthFrame   =   5.2*mm;
  const G4double phiFrame      = 67.50*degree; // centered at 90°
  const G4double alphaFrame    = 45.00*degree; // sector angle
  const G4double YoffsetFrame  = 108.3*mm;
}

class Tina : public NPS::VDetector{
    
public:
    Tina();
    virtual ~Tina();
    
    // Specific functions of this class ////////////////
public:
  void AddDetector(G4ThreeVector POS, double Alpha, string Shape);
  void AddDetector(double R,double Theta,double Phi,double T, double Alpha, string Shape);
  G4AssemblyVolume* BuildTTTPadDetector();
  G4AssemblyVolume* BuildYY1CsIDetector();
  
private:
  G4AssemblyVolume* m_TTTPad;
    G4AssemblyVolume* m_YY1CsI;
    
    // Inherited from NPS::VDetector class /////////////
public:
    // Read stream at Configfile to pick-up parameters of detector
    // called in DetectorConstruction::ReadDetectorConfiguration
    void ReadConfiguration(NPL::InputParser);
    
    // Construct detector and initialise sensitive part
    // called after DetectorConstruction::AddDetector
    void ConstructDetector(G4LogicalVolume* world);
    
    // Add detector branch to the EventTree
    // called after DetectorConstruction::AddDetector
    void InitializeRootOutput();
    
    // Read sensitive part and fill the Root tree
    // called in EventAction::EndOfEventAvtion
    void ReadSensitive(const G4Event* event);
    
    // Initialise all scorers used by the detector
    void InitializeScorers();
    G4MultiFunctionalDetector* m_TTTScorer;
    G4MultiFunctionalDetector* m_PadScorer;
    G4MultiFunctionalDetector* m_YY1Scorer;
    G4MultiFunctionalDetector* m_CsIScorer;
    
    ////////////////////////////////////////////////////
private:
    // Event class to store data
    TTinaData* m_Event;
    
private:
    // Geometry
    vector<double>  m_R;
    vector<double>  m_Theta;
    vector<double>  m_Phi;
    vector<double>  m_T;
  vector<double>  m_Alpha;
    vector<string>  m_Shape;
    
    // Visualisation
    G4VisAttributes* m_VisTTT;
    G4VisAttributes* m_VisPad;
    G4VisAttributes* m_VisYY1;
    G4VisAttributes* m_VisCsI;
    
public:
    // Dynamic loading of the library
    static NPS::VDetector* Construct();
    
};
#endif
