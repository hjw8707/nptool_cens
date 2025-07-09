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
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
    InitialConditions = new TInitialConditions();
    ReactionConditions = new TReactionConditions();

    PlungerPhysics = (TPlungerPhysics*)m_DetectorManager->GetDetector("Plunger");
    ASGARDPhysics = (TASGARDPhysics*)m_DetectorManager->GetDetector("ASGARD");

    //////////////////////////////////////////////////////////////
    // Initialize Root Input and Output: Should be put here (not automatically done)
    //////////////////////////////////////////////////////////////
    InitializeRootInput();
    InitializeRootOutput();
}

void Analysis::Clear() {
    flagVelocity = 0;
    flagKineticEnergy = 0;
    flagParticleName = "";
    clover_nbr.clear();
    energy.clear();
    doppler_energy.clear();
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
    Clear();

    //////////////////////////////////////////////////////////////
    // Plunger Data Analysis
    //////////////////////////////////////////////////////////////
    flagVelocity = PlungerPhysics->velocity[0];
    flagKineticEnergy = PlungerPhysics->kineticEnergy[0];
    flagParticleName = PlungerPhysics->particleName[0];

    //////////////////////////////////////////////////////////////
    // ASGARD Data Analysis
    //////////////////////////////////////////////////////////////
    // Re-analysis of the ASGARD data
    ASGARDPhysics->Clear();
    ASGARDPhysics->SetBeta(TVector3(0, 0, flagVelocity));
    ASGARDPhysics->BuildPhysicalEvent();

    for (int i = 0; i < ASGARDPhysics->AddBack_Clover.size(); i++) {
        clover_nbr.push_back(ASGARDPhysics->AddBack_Clover[i]);
        clover_theta.push_back(ASGARDPhysics->AddBack_Theta[i]);
        energy.push_back(ASGARDPhysics->AddBack_E[i]);
        doppler_energy.push_back(ASGARDPhysics->AddBack_DC[i]);
    }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}

///////////////////////////////////////////////////////////////////////////
void Analysis::InitializeRootInput() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("InitialConditions", true);
    inputChain->SetBranchAddress("InitialConditions", &InitialConditions);
    inputChain->SetBranchStatus("ReactionConditions", true);
    inputChain->SetBranchAddress("ReactionConditions", &ReactionConditions);
}

///////////////////////////////////////////////////////////////////////////
void Analysis::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    outputTree->Branch("flagParticleName", &flagParticleName, "flagParticleName/C");
    outputTree->Branch("flagVelocity", &flagVelocity, "flagVelocity/D");
    outputTree->Branch("flagKineticEnergy", &flagKineticEnergy, "flagKineticEnergy/D");

    outputTree->Branch("clover_nbr", &clover_nbr);
    outputTree->Branch("clover_theta", &clover_theta);
    outputTree->Branch("energy", &energy);
    outputTree->Branch("doppler_energy", &doppler_energy);
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
