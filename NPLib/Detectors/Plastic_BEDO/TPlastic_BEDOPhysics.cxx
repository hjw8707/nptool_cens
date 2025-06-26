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
 *  This class hold Plastic_BEDO Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TPlastic_BEDOPhysics.h"

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

ClassImp(TPlastic_BEDOPhysics)


///////////////////////////////////////////////////////////////////////////
TPlastic_BEDOPhysics::TPlastic_BEDOPhysics()
   : m_EventData(new TPlastic_BEDOData),
     m_PreTreatedData(new TPlastic_BEDOData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TPlastic_BEDOPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
} 
  
///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::BuildPhysicalEvent() {
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
void TPlastic_BEDOPhysics::PreTreat() {
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
      Double_t Energy = Cal->ApplyCalibration("Plastic_BEDO/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Plastic_BEDO/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigPlastic_BEDO.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigPlastic_BEDO.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigPlastic_BEDO.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigPlastic_BEDO.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigPlastic_BEDO";
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
void TPlastic_BEDOPhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Plastic_BEDO");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS","Shape"};
  vector<string> sphe = {"R","Theta","Phi","Shape"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Plastic_BEDO " << i+1 <<  endl;
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      string Shape = blocks[i]->GetString("Shape");
      AddDetector(Pos,Shape);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Plastic_BEDO " << i+1 <<  endl;
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
void TPlastic_BEDOPhysics::InitSpectra() {
  m_Spectra = new TPlastic_BEDOSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TPlastic_BEDOPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Plastic_BEDO", "D"+ NPL::itoa(i+1)+"_ENERGY","Plastic_BEDO_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Plastic_BEDO", "D"+ NPL::itoa(i+1)+"_TIME","Plastic_BEDO_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("Plastic_BEDO",  true );
  inputChain->SetBranchAddress("Plastic_BEDO", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Plastic_BEDO", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPlastic_BEDOPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Plastic_BEDO", "TPlastic_BEDOPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TPlastic_BEDOPhysics::Construct() {
  return (NPL::VDetector*) new TPlastic_BEDOPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Plastic_BEDO{
  public:
    proxy_Plastic_BEDO(){
      NPL::DetectorFactory::getInstance()->AddToken("Plastic_BEDO","Plastic_BEDO");
      NPL::DetectorFactory::getInstance()->AddDetector("Plastic_BEDO",TPlastic_BEDOPhysics::Construct);
    }
};

proxy_Plastic_BEDO p_Plastic_BEDO;
}

