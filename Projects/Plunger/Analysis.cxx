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
Analysis::~Analysis() {}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
    Plunger = (TPlungerPhysics*)m_DetectorManager->GetDetector("Plunger");
    ASGARD = (TASGARDPhysics*)m_DetectorManager->GetDetector("ASGARD");
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}

///////////////////////////////////////////////////////////////////////////
void Analysis::InitializeRootInput() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("Plunger", true);
    // inputChain->SetBranchAddress("Plunger", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void Analysis::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    // outputTree->Branch("Plunger", "TPlungerPhysics", &m_EventPhysics);
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
