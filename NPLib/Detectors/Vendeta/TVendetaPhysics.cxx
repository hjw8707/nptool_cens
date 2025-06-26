/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Vendeta Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TVendetaPhysics.h"

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

ClassImp(TVendetaPhysics)


  ///////////////////////////////////////////////////////////////////////////
TVendetaPhysics::TVendetaPhysics()
  : m_EventData(new TVendetaData),
  m_PreTreatedData(new TVendetaData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_E_RAW_Threshold(0), // adc channels
  m_E_Threshold(0),     // MeV
  m_AnodeNumber(1),
  m_NumberOfDetectors(0) {
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TVendetaPhysics::AddDetector(TVector3 ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
  m_DetectorPosition.push_back(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::BuildPhysicalEvent() {

  // Treat Event, only if Fission Chamber has triggered
  if(m_AnodeNumber==-1)
    return;

  // apply thresholds and calibration
  //if(m_AnodeNumber==0)
  //		return;

  PreTreat();

  // match energy and time together
  unsigned int mysizeLGE = m_PreTreatedData->GetLGMultEnergy();
  unsigned int mysizeHGE = m_PreTreatedData->GetHGMultEnergy();


  for (UShort_t e = 0; e < mysizeLGE ; e++) {

    LG_DetectorNumber.push_back(m_PreTreatedData->GetLGDetectorNbr(e));
    LG_Q1.push_back(m_PreTreatedData->GetLGQ1(e));
    LG_Q2.push_back(m_PreTreatedData->GetLGQ2(e));
    LG_Time.push_back(m_PreTreatedData->GetLGTime(e));
    LG_Qmax.push_back(m_PreTreatedData->GetLGQmax(e));
    LG_isSat.push_back(m_PreTreatedData->GetLGIsSat(e));
    LG_isGamma.push_back(m_PreTreatedData->GetLGIsGamma(e));

  }

  for (UShort_t e = 0; e < mysizeHGE ; e++) {
    HG_DetectorNumber.push_back(m_PreTreatedData->GetHGDetectorNbr(e));
    HG_Q1.push_back(m_PreTreatedData->GetHGQ1(e));
    HG_Q2.push_back(m_PreTreatedData->GetHGQ2(e));
    HG_Time.push_back(m_PreTreatedData->GetHGTime(e));
    HG_Qmax.push_back(m_PreTreatedData->GetHGQmax(e));
    HG_isSat.push_back(m_PreTreatedData->GetHGIsSat(e));
    HG_isGamma.push_back(m_PreTreatedData->GetHGIsGamma(e));
  }

  m_AnodeNumber=-1; 
}

///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();
  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  unsigned int mysizeLG = m_EventData->GetLGMultEnergy();
  unsigned int mysizeHG = m_EventData->GetHGMultEnergy();

  // LG pretreat
  for (UShort_t i = 0; i < mysizeLG ; ++i){
    int det = m_EventData->GetLGDetectorNbr(i);
    double Qmax = m_EventData->GetLGQmax(i);
    double TimeOffset=0;
    TimeOffset = Cal->GetValue("Vendeta/DET"+NPL::itoa(det)+"_LG_ANODE"+NPL::itoa(m_AnodeNumber)+"_TIMEOFFSET",0);
    if(m_AnodeNumber==0){
      // Apply calibration from Anode 6 in case of Pulser
      TimeOffset = Cal->GetValue("Vendeta/DET"+NPL::itoa(det)+"_LG_ANODE"+NPL::itoa(6)+"_TIMEOFFSET",0);
    }
    double Time = m_EventData->GetLGTime(i) + TimeOffset;
    m_PreTreatedData->SetLGDetectorNbr(det);
    m_PreTreatedData->SetLGQ1(m_EventData->GetLGQ1(i));
    m_PreTreatedData->SetLGQ2(m_EventData->GetLGQ2(i));
    m_PreTreatedData->SetLGTime(Time);
    m_PreTreatedData->SetLGQmax(Qmax);
    m_PreTreatedData->SetLGIsSat(m_EventData->GetLGIsSat(i));
    m_PreTreatedData->SetLGIsGamma(m_EventData->GetLGIsGamma(i));
  }

  // HG pretreat
  for (UShort_t i = 0; i < mysizeHG ; ++i){
    int det = m_EventData->GetHGDetectorNbr(i);
    double Qmax = m_EventData->GetHGQmax(i);
    double TimeOffset=0;
    TimeOffset = Cal->GetValue("Vendeta/DET"+NPL::itoa(det)+"_HG_ANODE"+NPL::itoa(m_AnodeNumber)+"_TIMEOFFSET",0);
    if(m_AnodeNumber==0){
      // Apply calibration from Anode 6 in case of Pulser
      TimeOffset = Cal->GetValue("Vendeta/DET"+NPL::itoa(det)+"_HG_ANODE"+NPL::itoa(6)+"_TIMEOFFSET",0);
    }

    double Time = m_EventData->GetHGTime(i) + TimeOffset;
    m_PreTreatedData->SetHGDetectorNbr(det);
    m_PreTreatedData->SetHGQ1(m_EventData->GetHGQ1(i));
    m_PreTreatedData->SetHGQ2(m_EventData->GetHGQ2(i));
    m_PreTreatedData->SetHGTime(Time);
    m_PreTreatedData->SetHGQmax(Qmax);
    m_PreTreatedData->SetHGIsSat(m_EventData->GetHGIsSat(i));
    m_PreTreatedData->SetHGIsGamma(m_EventData->GetHGIsGamma(i));
  }

}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigVendeta.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigVendeta.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigVendeta.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigVendeta.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigVendeta";
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
void TVendetaPhysics::Clear() {
  LG_DetectorNumber.clear();
  LG_Q1.clear();
  LG_Q2.clear();
  LG_Time.clear();
  LG_Qmax.clear();
  LG_isSat.clear();
  LG_isGamma.clear();

  HG_DetectorNumber.clear();
  HG_Q1.clear();
  HG_Q2.clear();
  HG_Time.clear();
  HG_Qmax.clear();
  HG_isSat.clear();
  HG_isGamma.clear();
  
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Vendeta");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Vendeta " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Vendeta " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      AddDetector(R,Theta,Phi);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::InitSpectra() {
  m_Spectra = new TVendetaSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TVendetaPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    for(int j = 0; j < 11; j++){
      Cal->AddParameter("Vendeta","DET"+NPL::itoa(i+1)+"_LG_ANODE"+NPL::itoa(j+1)+"_TIMEOFFSET","Vendeta_DET"+ NPL::itoa(i+1)+"_LG_ANODE"+NPL::itoa(j+1)+"_TIMEOFFSET");
      Cal->AddParameter("Vendeta","DET"+NPL::itoa(i+1)+"_HG_ANODE"+NPL::itoa(j+1)+"_TIMEOFFSET","Vendeta_DET"+ NPL::itoa(i+1)+"_HG_ANODE"+NPL::itoa(j+1)+"_TIMEOFFSET");
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Vendeta",  true );
  inputChain->SetBranchAddress("Vendeta", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Vendeta", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TVendetaPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Vendeta", "TVendetaPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TVendetaPhysics::Construct() {
  return (NPL::VDetector*) new TVendetaPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_Vendeta{
    public:
      proxy_Vendeta(){
        NPL::DetectorFactory::getInstance()->AddToken("Vendeta","Vendeta");
        NPL::DetectorFactory::getInstance()->AddDetector("Vendeta",TVendetaPhysics::Construct);
      }
  };

  proxy_Vendeta p_Vendeta;
}

