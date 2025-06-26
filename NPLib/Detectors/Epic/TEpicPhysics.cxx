/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon  contact address: audrey.chatillon@cea.fr                        *
 *                                                                           *
 * Creation Date  : décembre 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Epic Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TEpicPhysics.h"

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

ClassImp(TEpicPhysics)


///////////////////////////////////////////////////////////////////////////
TEpicPhysics::TEpicPhysics()
   : m_EventData(new TEpicData),
     m_PreTreatedData(new TEpicData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TEpicPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
} 
  
///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeE = m_PreTreatedData->GetMultiplicity();
  for (UShort_t e = 0; e < mysizeE ; e++) {
    int A = m_EventData->GetAnodeNbr(e);
    AnodeNumber.push_back(m_PreTreatedData->GetAnodeNbr(e));
    Q1.push_back(m_PreTreatedData->GetQ1(e));
    Qmax.push_back(m_PreTreatedData->GetQmax(e));
    Xpos.push_back(m_PreTreatedData->GetXpos(e));
    Ypos.push_back(m_PreTreatedData->GetYpos(e));
    Zpos.push_back(m_PreTreatedData->GetZpos(e));
  }
}

///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  unsigned int mysize = m_EventData->GetMultiplicity();
  for (UShort_t i = 0; i < mysize ; ++i) {
    Double_t Q1 = m_EventData->GetQ1(i);
    if (Q1 > m_E_Threshold) {
      Double_t Qmax = m_EventData->GetQmax(i);
      Double_t Xpos = m_EventData->GetXpos(i);
      Double_t Ypos = m_EventData->GetYpos(i);
      Double_t Zpos = m_EventData->GetZpos(i);
      unsigned short AnodeNumber = m_EventData->GetAnodeNbr(i);
      m_PreTreatedData->Set(AnodeNumber, Q1, Qmax, Xpos, Ypos, Zpos);
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigEpic.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigEpic.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigEpic.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigEpic.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigEpic";
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
void TEpicPhysics::Clear() {
  AnodeNumber.clear();
  Q1.clear();
  Qmax.clear();
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Epic");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Epic " << i+1 <<  endl;
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Epic " << i+1 <<  endl;
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(R,Theta,Phi,Shape);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::InitSpectra() {
  m_Spectra = new TEpicSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TEpicPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Epic", "A"+ NPL::itoa(i+1)+"_TIME_OFFSET","Epic_A"+ NPL::itoa(i+1)+"_TIMEOFFSET");
  }
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Epic",  true );
  inputChain->SetBranchAddress("Epic", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Epic", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TEpicPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Epic", "TEpicPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TEpicPhysics::Construct() {
  return (NPL::VDetector*) new TEpicPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Epic{
  public:
    proxy_Epic(){
      NPL::DetectorFactory::getInstance()->AddToken("Epic","Epic");
      NPL::DetectorFactory::getInstance()->AddDetector("Epic",TEpicPhysics::Construct);
    }
};

proxy_Epic p_Epic;
}

