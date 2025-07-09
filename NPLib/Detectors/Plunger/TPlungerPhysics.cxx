/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Plunger Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "TPlungerPhysics.h"

//   STL
#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
using namespace std;

//   NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "RootInput.h"
#include "RootOutput.h"

//   ROOT
#include "TChain.h"

ClassImp(TPlungerPhysics)

    ///////////////////////////////////////////////////////////////////////////
    TPlungerPhysics::TPlungerPhysics()
    : m_EventData(new TPlungerData),
      m_EventPhysics(this),
      m_TargetFound(false),
      m_TargetR(0),
      m_TargetThickness(0),
      m_TargetPosZ(0),
      m_TargetMaterial(""),
      m_StopperFound(false),
      m_StopperR(0),
      m_StopperThickness(0),
      m_StopperPosZ(0),
      m_StopperMaterial(""),
      m_ChamberFound(false),
      m_ChamberR(0),
      m_ChamberThickness(0),
      m_ChamberMaterial(""),
      m_ChamberPipeR(0),
      m_ChamberPipeZ0(0),
      m_ChamberPipeZ1(0) {}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::BuildPhysicalEvent() {
    // apply thresholds and calibration
    PreTreat();

    //////////////////////////////////////////////////////////////
    // Plunger Data Analysis
    //////////////////////////////////////////////////////////////
    for (int i = 0; i < m_EventData->GetMultPlunger(); i++) {
        if (m_EventData->GetParticleName(i) != "gamma" && m_EventData->GetParticleName(i) != "e-" &&
            m_EventData->GetVelocity(i) > 0) {
            particleName.push_back(m_EventData->GetParticleName(i));
            velocity.push_back(m_EventData->GetVelocity(i));
            kineticEnergy.push_back(m_EventData->GetKineticEnergy(i));
            break;
        }
    }
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::PreTreat() {}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::ReadAnalysisConfig() {}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::Clear() {
    particleName.clear();
    velocity.clear();
    kineticEnergy.clear();
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::ReadConfiguration(NPL::InputParser parser) {
    std::vector<NPL::InputBlock*> blocks;
    ////////////////////////////////////////////////////
    // Plunger Target
    ////////////////////////////////////////////////////
    blocks = parser.GetAllBlocksWithTokenAndValue("Plunger", "Target");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    if (blocks.size() > 0) {
        if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger Target " << endl;
        m_TargetFound = true;
        m_TargetR = blocks[0]->GetDouble("R", "mm");
        m_TargetThickness = blocks[0]->GetDouble("Thickness", "mm");
        m_TargetPosZ = blocks[0]->GetDouble("Z", "mm");
        m_TargetMaterial = blocks[0]->GetString("Material");
    }
    ////////////////////////////////////////////////////
    // Plunger Stopper
    ////////////////////////////////////////////////////
    blocks = parser.GetAllBlocksWithTokenAndValue("Plunger", "Stopper");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    if (blocks.size() > 0) {
        if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger Stopper " << endl;
        m_StopperFound = true;
        m_StopperR = blocks[0]->GetDouble("R", "mm");
        m_StopperThickness = blocks[0]->GetDouble("Thickness", "mm");
        m_StopperPosZ = blocks[0]->GetDouble("Z", "mm");
        m_StopperMaterial = blocks[0]->GetString("Material");
    }
    ////////////////////////////////////////////////////
    // Plunger Chamber
    ////////////////////////////////////////////////////
    blocks = parser.GetAllBlocksWithTokenAndValue("Plunger", "Chamber");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    if (blocks.size() > 0) {
        if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger Chamber " << endl;
        m_ChamberFound = true;
        m_ChamberR = blocks[0]->GetDouble("R", "mm");
        m_ChamberThickness = blocks[0]->GetDouble("Thickness", "mm");
        m_ChamberMaterial = blocks[0]->GetString("Material");
        m_ChamberPipeR = blocks[0]->GetDouble("PipeR", "mm");
        m_ChamberPipeZ0 = blocks[0]->GetDouble("PipeZ0", "mm");
        m_ChamberPipeZ1 = blocks[0]->GetDouble("PipeZ1", "mm");
    }
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::AddParameterToCalibrationManager() {}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::InitializeRootInputRaw() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("Plunger", true);
    inputChain->SetBranchAddress("Plunger", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::InitializeRootInputPhysics() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchAddress("Plunger", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    outputTree->Branch("Plunger", "TPlungerPhysics", &m_EventPhysics);
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TPlungerPhysics::Construct() { return (NPL::VDetector*)new TPlungerPhysics(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_Plunger {
   public:
    proxy_Plunger() {
        NPL::DetectorFactory::getInstance()->AddToken("Plunger", "Plunger");
        NPL::DetectorFactory::getInstance()->AddDetector("Plunger", TPlungerPhysics::Construct);
    }
};

proxy_Plunger p_Plunger;
}
