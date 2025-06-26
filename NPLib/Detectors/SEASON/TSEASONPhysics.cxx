/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Emmanuel Rey-herme                                       *
 * contact address: marine.vandebrouck@cea.fr                                *
 *                                                                           *
 * Creation Date  : septembre 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold SEASON Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "TSEASONPhysics.h"

//   STL
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>
using namespace std;

//   NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "RootInput.h"
#include "RootOutput.h"

//   ROOT
#include "TChain.h"

ClassImp(TSEASONPhysics)

    ///////////////////////////////////////////////////////////////////////////
    TSEASONPhysics::TSEASONPhysics()
    : m_EventData(new TSEASONData), m_PreTreatedData(new TSEASONData), m_EventPhysics(this), m_Spectra(0),
      m_E_RAW_Threshold(0), // adc channels
      m_E_Threshold(0),     // MeV
      m_NumberOfDetectors(0) {}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TSEASONPhysics::AddDetector(TVector3 X1_Y1, TVector3 X1_YMax, TVector3 XMax_Y1, TVector3 XMax_YMax,
                                 int NumberOfStripsX, int NumberOfStripsY) {
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
  m_NumberOfStripsX = NumberOfStripsX;
  m_NumberOfStripsY = NumberOfStripsY;
  // Vector U on Telescope Face (parallele to 0X axis) (parallele to X Strip) (NB: remember that Y strip are allong X
  // axis)
  TVector3 U = XMax_Y1 - X1_Y1;
  double FaceX = U.Mag();
  double Ushift = (U.Mag() - FaceX) / 2.;
  U = U.Unit();

  // Vector V on Telescope Face (parallele to 0Y axis) (parallele to Y Strip)
  TVector3 V = X1_YMax - X1_Y1;
  double FaceY = V.Mag();
  double Vshift = (V.Mag() - FaceY) / 2.;
  V = V.Unit();

  // Position Vector of Strip Center
  TVector3 StripCenter = TVector3(0, 0, 0);
  // Position Vector of X=1 Y=1 Strip
  TVector3 Strip_1_1;

  // Geometry Parameter
  double StripPitchX = FaceX / NumberOfStripsX; // mm
  double StripPitchY = FaceY / NumberOfStripsY; // mm

  // buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

  // Moving StripCEnter to 1.1 corner

  Strip_1_1 = X1_Y1 + U * (StripPitchX / 2.) + V * (StripPitchY / 2.);
  Strip_1_1 += U * Ushift + V * Vshift;

  for (int i = 0; i < NumberOfStripsX; ++i) {
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < NumberOfStripsY; ++j) {
      StripCenter = Strip_1_1 + StripPitchX * i * U + StripPitchY * i * V;
      lineX.push_back(StripCenter.X());
      lineY.push_back(StripCenter.Y());
      lineZ.push_back(StripCenter.Z());
    }

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeXE = m_PreTreatedData->GetXMultEnergy();
  unsigned int mysizeYE = m_PreTreatedData->GetYMultEnergy();
  unsigned int mysizeXT = m_PreTreatedData->GetXMultTime();
  unsigned int mysizeYT = m_PreTreatedData->GetYMultTime();
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();

  for (UShort_t e = 0; e < mysizeXE; e++) {
    for (UShort_t t = 0; t < mysizeXT; t++) {
      if (m_PreTreatedData->GetXE_DetectorNbr(e) == m_PreTreatedData->GetXT_DetectorNbr(t) &&
          m_PreTreatedData->GetXE_StripNbr(e) == m_PreTreatedData->GetXT_StripNbr(t)) {
        DetectorNumberX.push_back(m_PreTreatedData->GetXE_DetectorNbr(e));
        StripNumberX.push_back(m_PreTreatedData->GetXE_StripNbr(e));
        EnergyX.push_back(m_PreTreatedData->GetX_Energy(e));
        TimeX.push_back(m_PreTreatedData->GetX_Time(t));
      }
    }
  }
  for (UShort_t e = 0; e < mysizeYE; e++) {
    for (UShort_t t = 0; t < mysizeYT; t++) {
      if (m_PreTreatedData->GetYE_DetectorNbr(e) == m_PreTreatedData->GetYT_DetectorNbr(t) &&
          m_PreTreatedData->GetYE_StripNbr(e) == m_PreTreatedData->GetYT_StripNbr(t)) {
        DetectorNumberY.push_back(m_PreTreatedData->GetYE_DetectorNbr(e));
        StripNumberY.push_back(m_PreTreatedData->GetYE_StripNbr(e));
        EnergyY.push_back(m_PreTreatedData->GetY_Energy(e));
        TimeY.push_back(m_PreTreatedData->GetY_Time(t));
      }
    }
  }

  for (UShort_t e = 0; e < mysizeE; e++) {
    for (UShort_t t = 0; t < mysizeT; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t) &&
          m_PreTreatedData->GetE_StripNbrX(e) == m_PreTreatedData->GetT_StripNbrX(t) &&
          m_PreTreatedData->GetE_StripNbrY(e) == m_PreTreatedData->GetT_StripNbrY(t)) {
        DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
        StripX.push_back(m_PreTreatedData->GetE_StripNbrX(e));
        StripY.push_back(m_PreTreatedData->GetE_StripNbrY(e));
        Energy.push_back(m_PreTreatedData->Get_Energy(e));
        Time.push_back(m_PreTreatedData->Get_Time(t));
      }
    }
  }
  // DetectorNumber = DetectorNumberX;
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();
  unsigned int mysizeX = m_EventData->GetXMultEnergy();
  unsigned int mysizeY = m_EventData->GetYMultEnergy();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  if (mysizeX == mysizeY) {
    unsigned int mysize = mysizeX;
    // Energy
    for (UShort_t i = 0; i < mysize; ++i) {
      if (m_EventData->GetXE_DetectorNbr(i) == m_EventData->GetYE_DetectorNbr(i)) {
        if (m_EventData->GetX_Energy(i) > m_E_RAW_Threshold) {
          double x = (m_EventData->GetXE_StripNbr(i) - 16.5) * 2, y = (m_EventData->GetYE_StripNbr(i) - 16.5) * 2;
          // int dist = std::round(TMath::Sqrt(x*x+y*y+4)-2);
          int dist = std::round(TMath::Sqrt(y * y + 4) - 2);
          // if(dist > 38) dist = 38;
          string name = "SEASON/D";
          name += NPL::itoa(m_EventData->GetXE_DetectorNbr(i));
          name += "_DIST";
          name += NPL::itoa(dist);
          name += "_ENERGY";
          string namestrip = "SEASON/D";
          namestrip += NPL::itoa(m_EventData->GetXE_DetectorNbr(i));
          namestrip += "_STRIPX";
          namestrip += NPL::itoa(m_EventData->GetXE_StripNbr(i));
          namestrip += "_ENERGY";
          // cout << name << endl;
          Double_t Estrip = Cal->ApplyCalibration(namestrip, m_EventData->GetX_Energy(i));
          Double_t E = Cal->ApplyCalibration(name, Estrip);
          if (E > m_XE_Threshold) {
            m_PreTreatedData->SetEnergy(m_EventData->GetXE_DetectorNbr(i), m_EventData->GetXE_StripNbr(i),
                                        m_EventData->GetYE_StripNbr(i), E);
          }
        }
      }
    }

    // Time
    for (UShort_t i = 0; i < mysize; ++i) {
      if (m_EventData->GetXT_DetectorNbr(i) == m_EventData->GetYT_DetectorNbr(i)) {
        Double_t T = Cal->ApplyCalibration("SEASON/D" + NPL::itoa(m_EventData->GetXT_DetectorNbr(i)) + "_STRIPX" +
                                               NPL::itoa(m_EventData->GetXT_StripNbr(i)) + "_STRIPY" +
                                               NPL::itoa(m_EventData->GetYT_StripNbr(i)) + "_Time",
                                           m_EventData->GetX_Time(i));
        m_PreTreatedData->SetTime(m_EventData->GetXT_DetectorNbr(i), m_EventData->GetXT_StripNbr(i),
                                  m_EventData->GetYT_StripNbr(i), T);
      }
    }
  }

  // X Energy
  for (UShort_t i = 0; i < mysizeX; ++i) {
    if (m_EventData->GetX_Energy(i) > m_E_RAW_Threshold) {
      Double_t EX = Cal->ApplyCalibration("SEASON/D" + NPL::itoa(m_EventData->GetXE_DetectorNbr(i)) + "_STRIP" +
                                              NPL::itoa(m_EventData->GetXE_StripNbr(i)) + "_Energy",
                                          m_EventData->GetX_Energy(i));
      if (EX > m_XE_Threshold) {
        m_PreTreatedData->SetXEnergy(m_EventData->GetXE_DetectorNbr(i), m_EventData->GetXE_StripNbr(i), EX);
      }
    }
  }

  // X Time
  for (UShort_t i = 0; i < mysizeX; ++i) {
    Double_t TX = Cal->ApplyCalibration("SEASON/D" + NPL::itoa(m_EventData->GetXT_DetectorNbr(i)) + "_Strip" +
                                            NPL::itoa(m_EventData->GetXT_StripNbr(i)) + "_Time",
                                        m_EventData->GetX_Time(i));
    m_PreTreatedData->SetXTime(m_EventData->GetXT_DetectorNbr(i), m_EventData->GetXT_StripNbr(i), TX);
  }

  // Y Energy
  for (UShort_t i = 0; i < mysizeY; ++i) {
    if (m_EventData->GetY_Energy(i) > m_E_RAW_Threshold) {
      Double_t EY = Cal->ApplyCalibration("SEASON/D" + NPL::itoa(m_EventData->GetYE_DetectorNbr(i)) + "_STRIP" +
                                              NPL::itoa(m_EventData->GetYE_StripNbr(i)) + "_Energy",
                                          m_EventData->GetY_Energy(i));
      if (EY > m_YE_Threshold) {
        m_PreTreatedData->SetYEnergy(m_EventData->GetYE_DetectorNbr(i), m_EventData->GetYE_StripNbr(i), EY);
      }
    }
  }

  // Y Time
  for (UShort_t i = 0; i < mysizeY; ++i) {
    Double_t TY = Cal->ApplyCalibration("SEASON/D" + NPL::itoa(m_EventData->GetYT_DetectorNbr(i)) + "_Strip" +
                                            NPL::itoa(m_EventData->GetYT_StripNbr(i)) + "_Time",
                                        m_EventData->GetY_Time(i));
    m_PreTreatedData->SetYTime(m_EventData->GetYT_DetectorNbr(i), m_EventData->GetYT_StripNbr(i), TY);
  }
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigSEASON.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigSEASON.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigSEASON.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigSEASON.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer, DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigSEASON";
    if (LineBuffer.compare(0, name.length(), name) == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {
      whatToDo = "";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }
      else {
        ReadingStatus = false;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::Clear() {
  DetectorNumber.clear();
  StripNumberX.clear();
  StripNumberY.clear();
  StripX.clear();
  StripY.clear();
  Energy.clear();
  Time.clear();

  DetectorNumberX.clear();
  DetectorNumberY.clear();
  EnergyX.clear();
  EnergyY.clear();
  TimeX.clear();
  TimeY.clear();
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("SEASON");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl;

  vector<string> cart = {"X1_Y1", "X1_YMax", "XMax_Y1", "XMax_YMax", "NStripsX", "NStripsY", "Group"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  SEASON " << i + 1 << endl;

      TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
      TVector3 B = blocks[i]->GetTVector3("X1_YMax", "mm");
      TVector3 C = blocks[i]->GetTVector3("XMax_Y1", "mm");
      TVector3 D = blocks[i]->GetTVector3("XMax_YMax", "mm");
      int NumberOfStripsX = blocks[i]->GetInt("NStripsX");
      int NumberOfStripsY = blocks[i]->GetInt("NStripsY");
      AddDetector(A, B, C, D, NumberOfStripsX, NumberOfStripsY);
    }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}
///////////////////////////////////////////////////////////////////////////
TVector3 TSEASONPhysics::GetPositionOfInteraction(const int i) const {
  TVector3 Position;
  Position = TVector3(GetStripPositionX(DetectorNumber[i], StripNumberX[i], StripNumberY[i]),
                      GetStripPositionY(DetectorNumber[i], StripNumberX[i], StripNumberY[i]),
                      GetStripPositionZ(DetectorNumber[i], StripNumberX[i], StripNumberY[i]));
  return Position;
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::InitSpectra() {
  m_Spectra = new TSEASONSpectra(m_NumberOfDetectors, m_NumberOfStripsX, m_NumberOfStripsY);
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::ClearSpectra() {
  // To be done
}

///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TSEASONPhysics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
vector<TCanvas*> TSEASONPhysics::GetCanvas() {
  /*  if (m_Spectra)
      return m_Spectra->GetCanvas();
    else {*/
  vector<TCanvas*> empty;
  return empty;
  // }
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::WriteSpectra() { m_Spectra->WriteSpectra(); }

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::AddParameterToCalibrationManager() {
  int NumberOfStrips = 32;
  CalibrationManager* Cal = CalibrationManager::getInstance();
  vector<double> standard = {0, 1};
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    for (int k = 0; k < 50; k++) {
      Cal->AddParameter("SEASON", "D" + NPL::itoa(i + 1) + "_DIST" + NPL::itoa(k) + "_ENERGY",
                        "SEASON_D" + NPL::itoa(i + 1) + "_DIST" + NPL::itoa(k) + "_ENERGY", standard);
      Cal->AddParameter("SEASON", "D" + NPL::itoa(i + 1) + "_FOILDIST" + NPL::itoa(k) + "_ENERGY",
                        "SEASON_D" + NPL::itoa(i + 1) + "_FOILDIST" + NPL::itoa(k) + "_ENERGY", standard);
    }

    for (int j = 0; j < NumberOfStrips; ++j) {
      Cal->AddParameter("SEASON", "D" + NPL::itoa(i + 1) + "_STRIPX" + NPL::itoa(j + 1) + "_ENERGY",
                        "SEASON_D" + NPL::itoa(i + 1) + "_STRIPX" + NPL::itoa(j + 1) + "_ENERGY");
      Cal->AddParameter("SEASON", "D" + NPL::itoa(i + 1) + "_STRIPX" + NPL::itoa(j + 1) + "_TIME",
                        "SEASON_D" + NPL::itoa(i + 1) + "_STRIP" + NPL::itoa(j + 1) + "_TIME");
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("SEASON", true);
  inputChain->SetBranchAddress("SEASON", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("SEASON", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TSEASONPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("SEASON", "TSEASONPhysics", &m_EventPhysics);
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSEASONPhysics::Construct() { return (NPL::VDetector*)new TSEASONPhysics(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_SEASON {
 public:
  proxy_SEASON() {
    NPL::DetectorFactory::getInstance()->AddToken("SEASON", "SEASON");
    NPL::DetectorFactory::getInstance()->AddDetector("SEASON", TSEASONPhysics::Construct);
  }
};

proxy_SEASON p_SEASON;
}

