/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TAC Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TTACPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <chrono>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"

ClassImp(TTACPhysics)


///////////////////////////////////////////////////////////////////////////
TTACPhysics::TTACPhysics()
   : m_EventData(new TTACData),
     m_PreTreatedData(new TTACData),
     TSEvent(new TimeStamp),
     m_EventPhysics(this),
     m_Spectra(0),
     m_TAC_Time_RAW_Threshold(0), // adc channels
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TTACPhysics::AddDetector(TVector3 , string ){
} 

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::AddDetector(double R, double Theta, double Phi, string shape){
} 
  
///////////////////////////////////////////////////////////////////////////
void TTACPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::BuildPhysicalEvent() {
  
  ClaimReaderData();
 
 
  PreTreat();
  Match_TAC();
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::PreTreat() {
  ClearPreTreatedData();
  m_TAC_Mult = m_EventData->GetTAC_Mult();
  for (unsigned int i = 0; i < m_TAC_Mult; ++i) {
    if (m_EventData->GetTAC_Time(i) > m_TAC_Time_RAW_Threshold) {
      
      // If a TimeStamp reference has been added to perform the matching, prompt is asked
      if(!RefTS_Name.empty()){
        
        std::string TSName = m_EventData->GetTAC_Name(i);
        TSEvent->AddTimeStamp(TSName,m_EventData->GetTAC_TS(i));
        TSEvent->AddTimeStamp(RefTS_Name,RefTS);
        
        if(TSEvent->MatchTS(TSName)){
          m_PreTreatedData->SetTAC(m_EventData->GetTAC_N(i), m_EventData->GetTAC_Time(i), m_EventData->GetTAC_TS(i), m_EventData->GetTAC_Name(i));
        }
        TSEvent->ClearTimeStamps();
      }
      
      // else, all datas are filled
      else{
        m_PreTreatedData->SetTAC(m_EventData->GetTAC_N(i), m_EventData->GetTAC_Time(i), m_EventData->GetTAC_TS(i), m_EventData->GetTAC_Name(i));
      }
    }
  }

}

void TTACPhysics::ClaimReaderData() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }
}


void TTACPhysics::Match_TAC(){
  for(unsigned int i = 0; i < m_PreTreatedData->GetTAC_Mult(); i++){
    SortTAC[m_PreTreatedData->GetTAC_Name(i)] = std::make_pair(m_PreTreatedData->GetTAC_Time(i),m_PreTreatedData->GetTAC_TS(i));
  }
  for(auto it = SortTAC.begin(); it != SortTAC.end(); ++it){
  TAC_Name.push_back(it->first);
  TAC_Time.push_back((it->second).first);
  TAC_TS.push_back((it->second).second);
  }
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ReadAnalysisConfig() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::Clear() {
  SortTAC.clear();
  TAC_Name.clear();
  TAC_Time.clear();
  TAC_TS.clear();
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ReadConfiguration(NPL::InputParser parser) {
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ReadConfigurationTS(){
  TSEvent->ReadConfigurationFile();
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::SetRefTS(std::string TSRef_Name, unsigned long long TSRef){
  RefTS = TSRef;
  RefTS_Name = TSRef_Name; 
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::FillSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::CheckSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TTACPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TTACPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("TAC", "D"+ NPL::itoa(i+1)+"_ENERGY","TAC_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("TAC", "D"+ NPL::itoa(i+1)+"_TIME","TAC_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("TAC",  true );
  inputChain->SetBranchAddress("TAC", &m_EventData );
  }
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("TAC",  true );
  inputChain->SetBranchAddress("TAC", &m_EventPhysics);
  }
}



///////////////////////////////////////////////////////////////////////////
void TTACPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("TAC", "TTACPhysics", &m_EventPhysics);
}

void TTACPhysics::SetTreeReader(TTreeReader* TreeReader) {
   TTACPhysicsReader::r_SetTreeReader(TreeReader);
 }


////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TTACPhysics::Construct() {
  return (NPL::VDetector*) new TTACPhysics();
}

NPL::VTreeReader* TTACPhysics::ConstructReader() { return (NPL::VTreeReader*)new TTACPhysicsReader(); }



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_TAC{
  public:
    proxy_TAC(){
      NPL::DetectorFactory::getInstance()->AddToken("TAC","TAC");
      NPL::DetectorFactory::getInstance()->AddDetector("TAC",TTACPhysics::Construct);
      NPL::DetectorFactory::getInstance()->AddDetectorReader("TAC", TTACPhysics::ConstructReader);
    }
};

proxy_TAC p_TAC;
}

