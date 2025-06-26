/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_GAGG Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TTOGAXSI_GAGGPhysics.h"

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
#include "NPSystemOfUnits.h"
using namespace NPUNITS;

//   ROOT
#include "TChain.h"

ClassImp(TTOGAXSI_GAGGPhysics)


///////////////////////////////////////////////////////////////////////////
TTOGAXSI_GAGGPhysics::TTOGAXSI_GAGGPhysics() {
/*
     m_EventData = new TTOGAXSI_GAGGData;
     m_PreTreatedData = new TTOGAXSI_GAGGData;
     m_EventPhysics = this;
     m_Spectra = NULL;
     m_E_RAW_Threshold = 0; // adc channels
     m_E_Threshold = 0;     // MeV
     m_NumberOfRecoilDetectors = 0;
     m_NumberOfClusterDetectors = 0;
     EventMultiplicity = 0;

     //Detector
     Crystal_Length = 3.5*mm;
     Crystal_Width = 3.5*mm;
     Crystal_Height = 12*mm;
*/
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
//void TTOGAXSI_GAGGPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
//  m_NumberOfDetectors++;
//} 

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::AddRecoilDetector(TVector3 Pos, double Phi, TVector3 Ref) {
//  m_NumberOfRecoilDetectors++;
}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::AddClusterDetector(TVector3 Pos, double Phi, TVector3 Ref) {
//  m_NumberOfClusterDetectors++;
}

///////////////////////////////////////////////////////////////////////////
//void TTOGAXSI_GAGGPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
//  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
//  AddDetector(Pos,shape);
//} 
  
///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::BuildPhysicalEvent() {
/*
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  unsigned int mysizeRecoilE = m_PreTreatedData->GetRecoilMultEnergy();
  unsigned int mysizeRecoilT = m_PreTreatedData->GetRecoilMultTime();

  unsigned int mysizeClusterE = m_PreTreatedData->GetClusterMultEnergy();
  unsigned int mysizeClusterT = m_PreTreatedData->GetClusterMultTime();

  for (UShort_t e = 0; e < mysizeRecoilE ; e++) {
    for (UShort_t t = 0; t < mysizeRecoilT ; t++) {
      if (m_PreTreatedData->GetRecoil_E_DetectorNbr(e) == m_PreTreatedData->GetRecoil_T_DetectorNbr(t)) {
        RecoilDetectorNumber.push_back(m_PreTreatedData->GetRecoil_E_DetectorNbr(e));
        RecoilE.push_back(m_PreTreatedData->GetRecoil_Energy(e));
        RecoilT.push_back(m_PreTreatedData->GetRecoil_Time(t));
      }
    }
  }

  for (UShort_t e = 0; e < mysizeClusterE ; e++) {
    for (UShort_t t = 0; t < mysizeClusterT ; t++) {
      if (m_PreTreatedData->GetCluster_E_DetectorNbr(e) == m_PreTreatedData->GetCluster_T_DetectorNbr(t)) {
        ClusterDetectorNumber.push_back(m_PreTreatedData->GetCluster_E_DetectorNbr(e));
        ClusterE.push_back(m_PreTreatedData->GetCluster_Energy(e));
        ClusterT.push_back(m_PreTreatedData->GetCluster_Time(t));
      }
    }
  }
*/
}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector
/*
  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Recoil Detector 
  // Energy
  unsigned int mysize = m_EventData->GetRecoilMultEnergy();
  for (UShort_t i = 0; i < mysize ; ++i) {
    if (m_EventData->GetRecoil_Energy(i) > m_E_RAW_Threshold) {
	Double_t Energy = m_EventData->GetRecoil_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("TOGAXSI_GAGG/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetRecoilEnergy(m_EventData->GetRecoil_E_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetRecoilMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time = m_EventData->GetRecoil_Time(i);
    //Double_t Time= Cal->ApplyCalibration("TOGAXSI_GAGG/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetRecoilTime(m_EventData->GetRecoil_T_DetectorNbr(i), Time);
  }

  // Cluster Detector 
  // Energy
  mysize = m_EventData->GetClusterMultEnergy();
  for (UShort_t i = 0; i < mysize ; ++i) {
    if (m_EventData->GetCluster_Energy(i) > m_E_RAW_Threshold) {
	Double_t Energy = m_EventData->GetCluster_Energy(i);
      //Double_t Energy = Cal->ApplyCalibration("TOGAXSI_GAGG/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterEnergy(m_EventData->GetCluster_E_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetClusterMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time = m_EventData->GetCluster_Time(i);
    //Double_t Time= Cal->ApplyCalibration("TOGAXSI_GAGG/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetClusterTime(m_EventData->GetCluster_T_DetectorNbr(i), Time);
  }

*/
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigTOGAXSI_GAGG.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigTOGAXSI_GAGG.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigTOGAXSI_GAGG.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigTOGAXSI_GAGG.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigTOGAXSI_GAGG";
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
void TTOGAXSI_GAGGPhysics::Clear() {
/*
  EventMultiplicity = 0;

  RecoilDetectorNumber.clear();
  RecoilE.clear();
  RecoilT.clear();

  ClusterDetectorNumber.clear();
  ClusterE.clear();
  ClusterT.clear();
*/
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::ReadConfiguration(NPL::InputParser parser) {
/*
  // Recoil Detector
  vector<NPL::InputBlock*> blocks_recoil = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_GAGG","RecoilArray");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_recoil.size() << " recoil detectors found " << endl; 

  vector<string> coord = {"Pos","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_recoil.size() ; i++){
    if(blocks_recoil[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_GAGG Recoil " << i+1 <<  endl;
      TVector3 Pos = blocks_recoil[i]->GetTVector3("Pos","mm");
      double Phi = blocks_recoil[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_recoil[i]->GetTVector3("Ref","mm");
      AddRecoilDetector(Pos,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  // Recoil Detector
  vector<NPL::InputBlock*> blocks_cluster = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_GAGG","ClusterArray");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_cluster.size() << " cluster detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_cluster.size() ; i++){
    if(blocks_cluster[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_GAGG Cluster " << i+1 <<  endl;
      TVector3 Pos = blocks_cluster[i]->GetTVector3("Pos","mm");
      double Phi = blocks_cluster[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_cluster[i]->GetTVector3("Ref","mm");
      AddClusterDetector(Pos,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
*/

}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::InitSpectra() {
//  m_Spectra = new TTOGAXSI_GAGGSpectra(m_NumberOfRecoilDetectors, m_NumberOfClusterDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TTOGAXSI_GAGGPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
/*
  for (int i = 0; i < m_NumberOfRecoilDetectors; ++i) {
    Cal->AddParameter("TOGAXSI_GAGG", "D"+ NPL::itoa(i+1)+"_ENERGY","TOGAXSI_GAGG_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("TOGAXSI_GAGG", "D"+ NPL::itoa(i+1)+"_TIME","TOGAXSI_GAGG_D"+ NPL::itoa(i+1)+"_TIME");
  }

  for (int i = 0; i < m_NumberOfClusterDetectors; ++i) {
    Cal->AddParameter("TOGAXSI_GAGG", "D"+ NPL::itoa(i+1)+"_ENERGY","TOGAXSI_GAGG_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("TOGAXSI_GAGG", "D"+ NPL::itoa(i+1)+"_TIME","TOGAXSI_GAGG_D"+ NPL::itoa(i+1)+"_TIME");
  }
*/
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("TOGAXSI_GAGG",  true );
  inputChain->SetBranchAddress("TOGAXSI_GAGG", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("TOGAXSI_GAGG", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_GAGGPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("TOGAXSI_GAGG", "TTOGAXSI_GAGGPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TTOGAXSI_GAGGPhysics::Construct() {
  return (NPL::VDetector*) new TTOGAXSI_GAGGPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_TOGAXSI_GAGG{
  public:
    proxy_TOGAXSI_GAGG(){
      NPL::DetectorFactory::getInstance()->AddToken("TOGAXSI_GAGG","TOGAXSI_GAGG");
      NPL::DetectorFactory::getInstance()->AddDetector("TOGAXSI_GAGG",TTOGAXSI_GAGGPhysics::Construct);
    }
};

proxy_TOGAXSI_GAGG p_TOGAXSI_GAGG;
}

