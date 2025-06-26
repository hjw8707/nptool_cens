/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Léo Plagnol  contact address: leo.plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : January 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Coaxial_Germanium Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TCoaxial_GermaniumPhysics.h"

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

ClassImp(TCoaxial_GermaniumPhysics)


///////////////////////////////////////////////////////////////////////////
TCoaxial_GermaniumPhysics::TCoaxial_GermaniumPhysics()
   : m_EventData(new TCoaxial_GermaniumData),
     m_PreTreatedData(new TCoaxial_GermaniumData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TCoaxial_GermaniumPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
} 
  
///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::BuildPhysicalEvent() {
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
void TCoaxial_GermaniumPhysics::PreTreat() {
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
      Double_t Energy = Cal->ApplyCalibration("Coaxial_Germanium/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Coaxial_Germanium/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigCoaxial_Germanium.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigCoaxial_Germanium.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigCoaxial_Germanium.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigCoaxial_Germanium.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigCoaxial_Germanium";
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
void TCoaxial_GermaniumPhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Coaxial_Germanium");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Coaxial_Germanium " << i+1 <<  endl;
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Coaxial_Germanium " << i+1 <<  endl;
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
void TCoaxial_GermaniumPhysics::InitSpectra() {
  m_Spectra = new TCoaxial_GermaniumSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TCoaxial_GermaniumPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Coaxial_Germanium", "D"+ NPL::itoa(i+1)+"_ENERGY","Coaxial_Germanium_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Coaxial_Germanium", "D"+ NPL::itoa(i+1)+"_TIME","Coaxial_Germanium_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Coaxial_Germanium",  true );
  inputChain->SetBranchAddress("Coaxial_Germanium", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Coaxial_Germanium", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TCoaxial_GermaniumPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Coaxial_Germanium", "TCoaxial_GermaniumPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TCoaxial_GermaniumPhysics::Construct() {
  return (NPL::VDetector*) new TCoaxial_GermaniumPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Coaxial_Germanium{
  public:
    proxy_Coaxial_Germanium(){
      NPL::DetectorFactory::getInstance()->AddToken("Coaxial_Germanium","Coaxial_Germanium");
      NPL::DetectorFactory::getInstance()->AddDetector("Coaxial_Germanium",TCoaxial_GermaniumPhysics::Construct);
    }
};

proxy_Coaxial_Germanium p_Coaxial_Germanium;
}

