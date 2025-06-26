/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : September 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FissionChamber Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TFissionChamberPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"

ClassImp(TFissionChamberPhysics)


  ///////////////////////////////////////////////////////////////////////////
TFissionChamberPhysics::TFissionChamberPhysics()
  : m_EventData(new TFissionChamberData),
  m_PreTreatedData(new TFissionChamberData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_E_RAW_Threshold(0), // adc channels
  m_E_Threshold(0),     // MeV
  m_NumberOfDetectors(0) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TFissionChamberPhysics::AddDetector(TVector3 pos) {
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;

  for(int i =0; i<11; i++){
    TVector3 AnodePos(pos.X() , pos.Y(), pos.Z() - 27.5 + i*5.5);
    m_AnodePosition.push_back(AnodePos);
  } 

  for(int i=0; i<12; i++){
    LastTime[i] = 0;
    CurrentTime[i] = 0;
    counter[i] = 0;
  }
} 

///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  /* m_DetectorPosition.push_back(Pos); */
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  
  PreTreat();

  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();
  for (UShort_t e = 0; e < mysizeE ; e++) {

    int A = m_EventData->GetAnodeNbr(e);
    CurrentTime[A] = m_EventData->GetTime(e);
    double DT = CurrentTime[A] - LastTime[A];

    LastTime[A] = m_EventData->GetTime(e);

    DT_FC.push_back(DT);
    AnodeNumber.push_back(m_PreTreatedData->GetAnodeNbr(e));
    Q1.push_back(m_PreTreatedData->GetQ1(e));
    Q2.push_back(m_PreTreatedData->GetQ2(e));
    Qmax.push_back(m_PreTreatedData->GetQmax(e));
    Time.push_back(m_PreTreatedData->GetTime(e));
    isFakeFission.push_back(m_PreTreatedData->GetFakeFissionStatus(e));

    Time_HF.push_back(m_PreTreatedData->GetTimeHF(e));
    ToF.push_back(m_PreTreatedData->GetToF(e));
  }

}

///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();


  unsigned int mysize = m_EventData->GetMultiplicity();
  for (UShort_t i = 0; i < mysize ; ++i) {
    Double_t Q1 = m_EventData->GetQ1(i);
    Double_t Q2 = m_EventData->GetQ2(i);
    Double_t Qmax = m_EventData->GetQmax(i);

    if (Q1 > m_E_Threshold) {
      
      int AnodeNumber = m_EventData->GetAnodeNbr(i);
      double TimeOffset = 0;
      TimeOffset = Cal->GetValue("FissionChamber/ANODE"+NPL::itoa(AnodeNumber)+"_TIMEOFFSET",0);
      double Time = m_EventData->GetTime(i);
      double tof = Time - m_EventData->GetTimeHF(i) - TimeOffset;
      if(tof < 0){
        tof += Cal->GetValue("FissionChamber/PULSE_TIMEOFFSET",0) ;
      }
      m_PreTreatedData->SetAnodeNbr(AnodeNumber);
      m_PreTreatedData->SetQ1(Q1);
      m_PreTreatedData->SetQ2(Q2);
      m_PreTreatedData->SetQmax(Qmax);
      m_PreTreatedData->SetTime(Time);
      m_PreTreatedData->SetFakeFissionStatus(m_EventData->GetFakeFissionStatus(i));
      m_PreTreatedData->SetTimeHF(m_EventData->GetTimeHF(i));
      m_PreTreatedData->SetToF(tof);
    }
  }

}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigFissionChamber.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigFissionChamber.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigFissionChamber.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigFissionChamber.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigFissionChamber";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="E_THRESHOLD") {
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
void TFissionChamberPhysics::Clear() {
  AnodeNumber.clear();
  Q1.clear();
  Q2.clear();
  Qmax.clear();
  Time.clear();
  Time_HF.clear();
  ToF.clear();
  isFakeFission.clear();
  DT_FC.clear();
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("FissionChamber");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","GasMaterial","Pressure"};
  vector<string> sphe = {"R","Theta","Phi","GasMaterial","Pressure"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FissionChamber " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      string gas = blocks[i]->GetString("GasMaterial");
      double pressure = blocks[i]->GetDouble("Pressure","bar");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FissionChamber " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string gas = blocks[i]->GetString("GasMaterial");
      double pressure = blocks[i]->GetDouble("Pressure","bar");
      AddDetector(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::InitSpectra() {
  m_Spectra = new TFissionChamberSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TFissionChamberPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i <= 11; ++i) { // 11 anodes
    Cal->AddParameter("FissionChamber", "ANODE"+ NPL::itoa(i)+"_TIMEOFFSET","FissionChamber_ANODE"+ NPL::itoa(i)+"_TIMEOFFSET");
  }
    Cal->AddParameter("FissionChamber", "PULSE_TIMEOFFSET","FissionChamber_PULSE_TIMEOFFSET");
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("FissionChamber",  true );
  inputChain->SetBranchAddress("FissionChamber", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("FissionChamber", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TFissionChamberPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("FissionChamber", "TFissionChamberPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TFissionChamberPhysics::Construct() {
  return (NPL::VDetector*) new TFissionChamberPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_FissionChamber{
    public:
      proxy_FissionChamber(){
        NPL::DetectorFactory::getInstance()->AddToken("FissionChamber","FissionChamber");
        NPL::DetectorFactory::getInstance()->AddDetector("FissionChamber",TFissionChamberPhysics::Construct);
      }
  };

  proxy_FissionChamber p_FissionChamber;
}

