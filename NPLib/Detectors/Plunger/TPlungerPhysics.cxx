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
      m_PreTreatedData(new TPlungerData),
      m_EventPhysics(this),
      m_E_RAW_Threshold(0),  // adc channels
      m_E_Threshold(0),      // MeV
      m_NumberOfDetectors(0),
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
/// A usefull method to bundle all operation to add a detector
void TPlungerPhysics::AddDetector(TVector3, string) {
    // In That simple case nothing is done
    // Typically for more complex detector one would calculate the relevant
    // positions (stripped silicon) or angles (gamma array)
    m_NumberOfDetectors++;
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::AddDetector(double R, double Theta, double Phi, string shape) {
    // Compute the TVector3 corresponding
    TVector3 Pos(R * sin(Theta) * cos(Phi), R * sin(Theta) * sin(Phi), R * cos(Theta));
    // Call the cartesian method
    AddDetector(Pos, shape);
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::BuildPhysicalEvent() {
    // apply thresholds and calibration
    PreTreat();

    // match energy and time together
    // unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
    // unsigned int mysizeT = m_PreTreatedData->GetMultTime();
    // for (UShort_t e = 0; e < mysizeE; e++) {
    //     for (UShort_t t = 0; t < mysizeT; t++) {
    //         if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
    //             DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
    //             Energy.push_back(m_PreTreatedData->Get_Energy(e));
    //             Time.push_back(m_PreTreatedData->Get_Time(t));
    //         }
    //     }
    // }
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::PreTreat() {
    // This method typically applies thresholds and calibrations
    // Might test for disabled channels for more complex detector

    // clear pre-treated object
    ClearPreTreatedData();

    // instantiate CalibrationManager
    static CalibrationManager* Cal = CalibrationManager::getInstance();

    // Energy
    // unsigned int mysize = m_EventData->GetMultEnergy();
    // for (UShort_t i = 0; i < mysize; ++i) {
    //     if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
    //         Double_t Energy = Cal->ApplyCalibration("Plunger/ENERGY" + NPL::itoa(m_EventData->GetE_DetectorNbr(i)),
    //                                                 m_EventData->Get_Energy(i));
    //         if (Energy > m_E_Threshold) {
    //             m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
    //         }
    //     }
    // }

    // // Time
    // mysize = m_EventData->GetMultTime();
    // for (UShort_t i = 0; i < mysize; ++i) {
    //     Double_t Time = Cal->ApplyCalibration("Plunger/TIME" + NPL::itoa(m_EventData->GetT_DetectorNbr(i)),
    //                                           m_EventData->Get_Time(i));
    //     m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
    // }
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::ReadAnalysisConfig() {
    bool ReadingStatus = false;

    // path to file
    string FileName = "./configs/ConfigPlunger.dat";

    // open analysis config file
    ifstream AnalysisConfigFile;
    AnalysisConfigFile.open(FileName.c_str());

    if (!AnalysisConfigFile.is_open()) {
        cout << " No ConfigPlunger.dat found: Default parameter loaded for Analayis " << FileName << endl;
        return;
    }
    cout << " Loading user parameter for Analysis from ConfigPlunger.dat " << endl;

    // Save it in a TAsciiFile
    TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
    asciiConfig->AppendLine("%%% ConfigPlunger.dat %%%");
    asciiConfig->Append(FileName.c_str());
    asciiConfig->AppendLine("");
    // read analysis config file
    string LineBuffer, DataBuffer, whatToDo;
    while (!AnalysisConfigFile.eof()) {
        // Pick-up next line
        getline(AnalysisConfigFile, LineBuffer);

        // search for "header"
        string name = "ConfigPlunger";
        if (LineBuffer.compare(0, name.length(), name) == 0) ReadingStatus = true;

        // loop on tokens and data
        while (ReadingStatus) {
            whatToDo = "";
            AnalysisConfigFile >> whatToDo;

            // Search for comment symbol (%)
            if (whatToDo.compare(0, 1, "%") == 0) {
                AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
            }

            else if (whatToDo == "E_RAW_THRESHOLD") {
                AnalysisConfigFile >> DataBuffer;
                m_E_RAW_Threshold = atof(DataBuffer.c_str());
                cout << whatToDo << " " << m_E_RAW_Threshold << endl;
            }

            else if (whatToDo == "E_THRESHOLD") {
                AnalysisConfigFile >> DataBuffer;
                m_E_Threshold = atof(DataBuffer.c_str());
                cout << whatToDo << " " << m_E_Threshold << endl;
            }

            else {
                ReadingStatus = false;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
void TPlungerPhysics::Clear() {
    DetectorNumber.clear();
    Energy.clear();
    Time.clear();
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
void TPlungerPhysics::AddParameterToCalibrationManager() {
    CalibrationManager* Cal = CalibrationManager::getInstance();
    for (int i = 0; i < m_NumberOfDetectors; ++i) {
        Cal->AddParameter("Plunger", "D" + NPL::itoa(i + 1) + "_ENERGY", "Plunger_D" + NPL::itoa(i + 1) + "_ENERGY");
        Cal->AddParameter("Plunger", "D" + NPL::itoa(i + 1) + "_TIME", "Plunger_D" + NPL::itoa(i + 1) + "_TIME");
    }
}

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
