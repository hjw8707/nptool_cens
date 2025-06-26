/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold quadranMSQ25 Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TquadranMSQ25Physics.h"

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

ClassImp(TquadranMSQ25Physics)


///////////////////////////////////////////////////////////////////////////
TquadranMSQ25Physics::TquadranMSQ25Physics()
   : m_EventData(new TquadranMSQ25Data),
     m_PreTreatedData(new TquadranMSQ25Data),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TquadranMSQ25Physics::AddDetector(TVector3 , double){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::AddDetector(double R, double Theta, double Phi, double Thickness){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,Thickness);
} 
  
///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    for (UShort_t t = 0; t < mysizeT ; t++) {
      if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
        DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
        Energy.push_back(m_PreTreatedData->Get_Energy(e));
        Time.push_back(m_PreTreatedData->Get_Time(t));
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  unsigned int mysize = m_EventData->GetMultEnergy();
  for (UShort_t i = 0; i < mysize ; ++i) {
    if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = Cal->ApplyCalibration("quadranMSQ25/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("quadranMSQ25/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigquadranMSQ25.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigquadranMSQ25.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigquadranMSQ25.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigquadranMSQ25.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigquadranMSQ25";
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
void TquadranMSQ25Physics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("quadranMSQ25");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Thickness"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  quadranMSQ25 " << i+1 <<  endl;
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      double Thickness = blocks[i]->GetDouble("Thickness", "um");
      AddDetector(Pos,Thickness);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  quadranMSQ25 " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double Thickness = blocks[i]->GetDouble("Thickness", "um");
      AddDetector(R,Theta,Phi,Thickness);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::InitSpectra() {
  m_Spectra = new TquadranMSQ25Spectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TquadranMSQ25Physics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("quadranMSQ25", "D"+ NPL::itoa(i+1)+"_ENERGY","quadranMSQ25_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("quadranMSQ25", "D"+ NPL::itoa(i+1)+"_TIME","quadranMSQ25_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("quadranMSQ25",  true );
  inputChain->SetBranchAddress("quadranMSQ25", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("quadranMSQ25", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TquadranMSQ25Physics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("quadranMSQ25", "TquadranMSQ25Physics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TquadranMSQ25Physics::Construct() {
  return (NPL::VDetector*) new TquadranMSQ25Physics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_quadranMSQ25{
  public:
    proxy_quadranMSQ25(){
      NPL::DetectorFactory::getInstance()->AddToken("quadranMSQ25","quadranMSQ25");
      NPL::DetectorFactory::getInstance()->AddDetector("quadranMSQ25",TquadranMSQ25Physics::Construct);
    }
};

proxy_quadranMSQ25 p_quadranMSQ25;
}

