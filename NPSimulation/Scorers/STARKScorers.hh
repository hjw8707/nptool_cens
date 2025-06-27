#ifndef STARKScorers_h
#define STARKScorers_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : February 2013                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  File old the scorer specific to the STARK Detector                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new style of scorer is aim to become the standard way of doing scorer*
 * in NPTool.                                                                *
 *The index is build using the TrackID, Detector Number and Strip Number.    *
 *The scorer Hold Energy and time together                                   *
 *Only one scorer is needed for a detector                                   *
 *****************************************************************************/
#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"
#include "NPImage.h"
#include <map>
using namespace std;
using namespace CLHEP;

namespace STARKSCORERS {

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  class PS_STARK_X6: public G4VPrimitiveScorer{

  public: // with description
    PS_STARK_X6(G4String name, G4int Level,
		G4double stripPlaneX, G4double stripPlaneY,
		G4int nFrontStrip, G4int nBackStrip,
		G4int depth=0);
    
    ~PS_STARK_X6();
    
  protected: // with description
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    
  public:
    void Initialize(G4HCofThisEvent*);
    void EndOfEvent(G4HCofThisEvent*);
    void clear();
    void DrawAll();
    void PrintAll();
    
  private: // Geometry of the detector
    G4double m_stripPlaneX;
    G4double m_stripPlaneY;
    G4int    m_nFrontStrip;
    G4int    m_nBackStrip;
    G4double m_frontStripPitch;
    G4double m_backStripPitch;
    // Level at which to find the copy number linked to the detector number
    G4int    m_Level;
    
  private: // inherited from G4VPrimitiveScorer
    G4int HCID;
    NPS::HitsMap<G4double*>* EvtMap;
    
  private: // Needed for intermediate calculation (avoid multiple instantiation in Processing Hit)
    G4ThreeVector m_Position  ;
    G4int m_detectorNumber    ;
    G4int m_frontStripNumber  ;
    G4int m_backStripNumber  ;
    G4long m_Index            ;
  };


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  class PS_STARK_BB10: public G4VPrimitiveScorer{

  public: // with description
    PS_STARK_BB10(G4String name, G4int Level,
		  G4double geoX, G4double geoY,
		  G4int nStrip,  G4int depth=0);
    
    ~PS_STARK_BB10();
    
  protected: // with description
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    
  public:
    void Initialize(G4HCofThisEvent*);
    void EndOfEvent(G4HCofThisEvent*);
    void clear();
    void DrawAll();
    void PrintAll();
    
  private: // Geometry of the detector
    G4double m_geoX;
    G4double m_geoY;
    G4int    m_nStrip;
    G4double m_stripPitch;
    // Level at which to find the copy number linked to the detector number
    G4int    m_Level;
    
  private: // inherited from G4VPrimitiveScorer
    G4int HCID;
    NPS::HitsMap<G4double*>* EvtMap;
    
  private: // Needed for intermediate calculation (avoid multiple instantiation in Processing Hit)
    G4ThreeVector m_Position  ;
    G4int m_detectorNumber    ;
    G4int m_stripNumber  ;
    G4long m_Index            ;
  };


  class PS_STARK_QQQ5: public G4VPrimitiveScorer{

  public: // with description
    PS_STARK_QQQ5(G4String name, G4int Level,
		  G4double inR, G4double outR, G4double phi0, G4double phi1,
		  G4int NRStrip,G4int NAStrip, G4int depth=0);
    ~PS_STARK_QQQ5();
    
  protected: // with description
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    
  public:
    void Initialize(G4HCofThisEvent*);
    void EndOfEvent(G4HCofThisEvent*);
    void clear();
    void DrawAll();
    void PrintAll();
    
  private: // Geometry of the detector
    G4double m_inR, m_outR, m_phi0, m_phi1;
    G4int    m_NRStrip, m_NAStrip;
    G4double m_RStripPitch, m_AStripPitch;
    // Level at which to find the copy number linked to the detector number
    G4int    m_Level;
    
  private: // inherited from G4VPrimitiveScorer
    G4int HCID;
    NPS::HitsMap<G4double*>* EvtMap;
    
  private: // Needed for intermediate calculation (avoid multiple instantiation in Processing Hit)
    G4ThreeVector m_Position  ;
    G4int m_detectorNumber    ;
    G4int m_radialStripNumber      ;
    G4int m_annularStripNumber      ;
    G4long m_Index            ;
  };
}
#endif
