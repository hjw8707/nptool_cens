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
      m_NumberOfDetectors(0) {}

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
    unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
    unsigned int mysizeT = m_PreTreatedData->GetMultTime();
    for (UShort_t e = 0; e < mysizeE; e++) {
        for (UShort_t t = 0; t < mysizeT; t++) {
            if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
                DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
                Energy.push_back(m_PreTreatedData->Get_Energy(e));
                Time.push_back(m_PreTreatedData->Get_Time(t));
            }
        }
    }
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
    unsigned int mysize = m_EventData->GetMultEnergy();
    for (UShort_t i = 0; i < mysize; ++i) {
        if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
            Double_t Energy = Cal->ApplyCalibration("Plunger/ENERGY" + NPL::itoa(m_EventData->GetE_DetectorNbr(i)),
                                                    m_EventData->Get_Energy(i));
            if (Energy > m_E_Threshold) {
                m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
            }
        }
    }

    // Time
    mysize = m_EventData->GetMultTime();
    for (UShort_t i = 0; i < mysize; ++i) {
        Double_t Time = Cal->ApplyCalibration("Plunger/TIME" + NPL::itoa(m_EventData->GetT_DetectorNbr(i)),
                                              m_EventData->Get_Time(i));
        m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
    }
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
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Plunger");
    if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;

    vector<string> cart = {"POS", "Shape"};
    vector<string> sphe = {"R", "Theta", "Phi", "Shape"};

    for (unsigned int i = 0; i < blocks.size(); i++) {
        if (blocks[i]->HasTokenList(cart)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger " << i + 1 << endl;

            TVector3 Pos = blocks[i]->GetTVector3("POS", "mm");
            string Shape = blocks[i]->GetString("Shape");
            AddDetector(Pos, Shape);
        } else if (blocks[i]->HasTokenList(sphe)) {
            if (NPOptionManager::getInstance()->GetVerboseLevel()) cout << endl << "////  Plunger " << i + 1 << endl;
            double R = blocks[i]->GetDouble("R", "mm");
            double Theta = blocks[i]->GetDouble("Theta", "deg");
            double Phi = blocks[i]->GetDouble("Phi", "deg");
            string Shape = blocks[i]->GetString("Shape");
            AddDetector(R, Theta, Phi, Shape);
        } else {
            cout << "ERROR: check your input file formatting " << endl;
            exit(1);
        }
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
