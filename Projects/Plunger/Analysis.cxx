/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe  Plunger analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include <iostream>
using namespace std;
#include "Analysis.h"
#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TChain.h"
#include "TTree.h"

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis() {}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis() {
    if (InitialConditions) delete InitialConditions;
    if (ReactionConditions) delete ReactionConditions;
    if (PlungerData) delete PlungerData;
    if (ASGARDData) delete ASGARDData;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
    InitialConditions = new TInitialConditions();
    ReactionConditions = new TReactionConditions();
    PlungerData = new TPlungerData();
    ASGARDData = new TASGARDData();

    PlungerPhysics = (TPlungerPhysics*)m_DetectorManager->GetDetector("Plunger");
    ASGARDPhysics = (TASGARDPhysics*)m_DetectorManager->GetDetector("ASGARD");

    //////////////////////////////////////////////////////////////
    // Initialize Root Input and Output: Should be put here (not automatically done)
    //////////////////////////////////////////////////////////////
    InitializeRootInput();
    InitializeRootOutput();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
    //////////////////////////////////////////////////////////////
    // Plunger Data Analysis
    //////////////////////////////////////////////////////////////
    string particleFilter = "C15";
    for (int i = 0; i < PlungerData->GetMultPlunger(); i++) {
        if (PlungerData->GetParticleName(i) == particleFilter) {
            flagVelocity = PlungerData->GetVelocity(i);
            flagKineticEnergy = PlungerData->GetKineticEnergy(i);
            break;
        }
    }

    //////////////////////////////////////////////////////////////
    // ASGARD Data Analysis
    //////////////////////////////////////////////////////////////
    // ASGARDPhysics->BuildPhysicalEvent();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}

///////////////////////////////////////////////////////////////////////////
void Analysis::InitializeRootInput() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("Plunger", true);
    inputChain->SetBranchAddress("Plunger", &PlungerData);
    inputChain->SetBranchStatus("ASGARD", true);
    inputChain->SetBranchAddress("ASGARD", &ASGARDData);
    inputChain->SetBranchStatus("InitialConditions", true);
    inputChain->SetBranchAddress("InitialConditions", &InitialConditions);
    inputChain->SetBranchStatus("ReactionConditions", true);
    inputChain->SetBranchAddress("ReactionConditions", &ReactionConditions);
}

///////////////////////////////////////////////////////////////////////////
void Analysis::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    outputTree->Branch("flagVelocity", &flagVelocity, "flagVelocity/D");
    outputTree->Branch("flagKineticEnergy", &flagKineticEnergy, "flagKineticEnergy/D");
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct() { return (NPL::VAnalysis*)new Analysis(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy {
   public:
    proxy() { NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct); }
};

proxy p;
}
