#ifndef DriftChamberScorers_h
#define DriftChamberScorers_h
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sunghan Bae     contact adresse: shbae2703@gmail.com     *  
 *                                                                           *
 * Creation Date  : Feb. 2024                                                *
 * Last update    : 22/02/24                                                 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Scorer for VOICE ionization chamber					     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *****************************************************************************/
#include "G4VFastSimulationModel.hh"
#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"
#include <map>
using namespace std;
using namespace CLHEP;

namespace VOICEScorers {

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    class PS_Electrode: public G4VPrimitiveScorer{

	public: // with description
		  PS_Electrode(G4String name, G4int Level, G4double OpenX, G4double OpenY, G4double PCBZ, G4int depth);
		  ~PS_Electrode();

	protected:
		  G4bool ProcessHits(G4Step*, G4TouchableHistory*);

	public:
		void Initialize(G4HCofThisEvent*);
		void EndOfEvent(G4HCofThisEvent*);
		void clear();
		void DrawAll();
		void PrintAll();
	private:
		G4double m_OpenX;
		G4double m_OpenY;
		G4double m_PCBZ;
		G4int	m_Level;
		G4int m_detectorNumber;

	private:
		G4int HCID;
		NPS::HitsMap<G4double*>* EvtMap;
		std::map<long, long> multinessmap;
		long multiness;

	private:
		G4ThreeVector m_position;
		G4long m_index;
    };


    class PS_END: public G4VPrimitiveScorer{

	public: // with description
		  PS_END(G4String name, G4int Level, G4double ENDthick, G4int depth);
		  ~PS_END();

	protected:
		  G4bool ProcessHits(G4Step*, G4TouchableHistory*);

	public:
		void Initialize(G4HCofThisEvent*);
		void EndOfEvent(G4HCofThisEvent*);
		void clear();
		void DrawAll();
		void PrintAll();
	private:
		G4double m_ENDthick;
		G4int	m_Level;
		G4int m_detectorNumber;

	private:
		G4int HCID;
		NPS::HitsMap<G4double*>* EvtMap;
		std::map<long, long> multinessmap;
		long multiness;

	private:
		G4ThreeVector m_position;
		G4long m_index;
    };
    
    class PS_Gas: public G4VPrimitiveScorer{

	public: // with description
		  PS_Gas(G4String name, vector<double> zpos, G4double StepSize, G4int Level, G4int depth);
		  ~PS_Gas();

	protected:
		  G4bool ProcessHits(G4Step*, G4TouchableHistory*);

	public:
		void Initialize(G4HCofThisEvent*);
		void EndOfEvent(G4HCofThisEvent*);
		void clear();
		void DrawAll();
		void PrintAll();
	private:
		G4int	m_Level;
		G4int m_detectorNumber;
		vector<double> m_zpos;
		vector<int> m_checker;
		G4int m_i_stage;
		G4double m_StepSize;

	private:
		G4int HCID;
		NPS::HitsMap<G4double*>* EvtMap;
		std::map<long, long> multinessmap;
		std::map<long, int> gasmulmap;
		long multiness;
		G4bool endflag;

	private:
		G4ThreeVector m_position;
		G4long m_index;
    };
} 

#endif
