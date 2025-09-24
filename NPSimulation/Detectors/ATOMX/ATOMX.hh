#ifndef ATOMX_h
#define ATOMX_h 1
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  ATOMX simulation                             *
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
#include "G4UserLimits.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4VFastSimulationModel.hh"
#include "G4FastSimulationManager.hh"
#include "G4MultiFunctionalDetector.hh"

// NPTool header
#include "Decay.hh"
#include "TATOMXData.h"
#include "NPSVDetector.hh"
#include "NPInputParser.h"
#include "BeamReaction.hh"
class STARK;

class ATOMX : public NPS::VDetector{
    public:
        ATOMX() ;
        virtual ~ATOMX() ;

        void ReadConfiguration(NPL::InputParser) ;
        void ConstructDetector(G4LogicalVolume* world) ;
        void InitializeRootOutput() ;
        void ReadSensitive(const G4Event* event) ;
        void InitializeScorers() ;

    private:
        G4LogicalVolume* BuildDetector();
        G4MultiFunctionalDetector* m_ATOMXScorer ;

        G4LogicalVolume* m_logicChamber;
        G4LogicalVolume* m_logicGas;
        G4LogicalVolume* m_sensitive;

        TATOMXData* m_Event;

        G4ThreeVector m_Pos;
        G4ThreeVector m_ChamberSize;
        G4ThreeVector m_GasSize;
        G4ThreeVector m_PadSize;
        G4ThreeVector m_MMSSize;
        double m_ReactionBoxSize;
        double m_Mylar_Rmax;
        double m_Mylar_Thickness;

        vector<string> m_GasMaterial;
        vector<int> m_GasFraction;
        double m_Pressure;
        double m_Temperature;
        double m_ReactionZ1 = 0;
        double m_ReactionZ2 = 0;
        double m_ReactionZWidth = 0;
        double m_EProductionCut = 1000;
        double m_StepLimit = 0.5;

        G4Region* m_StepLimitRegion;
        vector<G4VFastSimulationModel*> m_ReactionModel;

        int fCopyNo = 80000;

    public:
        static NPS::VDetector* Construct();
};
#endif
