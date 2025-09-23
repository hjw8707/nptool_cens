/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ATOMX Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TATOMXPhysics.h"

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

ClassImp(TATOMXPhysics)


///////////////////////////////////////////////////////////////////////////
TATOMXPhysics::TATOMXPhysics()
   : m_EventData(new TATOMXData),
     m_PreTreatedData(new TATOMXData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TATOMXPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::AddDetector(double R, double Theta, double Phi, string shape){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos,shape);
} 
  
///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  // match energy and time together
  //unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
  //unsigned int mysizeT = m_PreTreatedData->GetMultTime();
  //for (UShort_t e = 0; e < mysizeE ; e++) {
  //  for (UShort_t t = 0; t < mysizeT ; t++) {
  //    if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
  //      DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
  //      Energy.push_back(m_PreTreatedData->Get_Energy(e));
  //      Time.push_back(m_PreTreatedData->Get_Time(t));
  //    }
  //  }
  //}
}

///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  //unsigned int mysize = m_EventData->GetMultEnergy();
  //for (UShort_t i = 0; i < mysize ; ++i) {
  //  if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
  //    Double_t Energy = Cal->ApplyCalibration("ATOMX/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
  //    if (Energy > m_E_Threshold) {
  //      m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
  //    }
  //  }
  //}

  //// Time 
  //mysize = m_EventData->GetMultTime();
  //for (UShort_t i = 0; i < mysize; ++i) {
  //  Double_t Time= Cal->ApplyCalibration("ATOMX/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
  //  m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  //}
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigATOMX.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigATOMX.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigATOMX.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigATOMX.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigATOMX";
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
void TATOMXPhysics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::ReadConfiguration(NPL::InputParser parser){

    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ATOMX");
    if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl; 

    vector<string> tokens = { "TPC_Pos" };

    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        if(blocks[i]->HasTokenList(tokens))
        {
            /*
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  ATOMX " << i+1 <<  endl;
            vector<string> GasName;
            vector<int> GasFraction;
            vector<int> vReactionZ;
            if (blocks[i]->HasToken("TPC_Pos"             )) m_Pos = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Pos","mm"));
            if (blocks[i]->HasToken("TPC_Chamber_Size"    )) m_ChamberSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Chamber_Size","mm"));
            if (blocks[i]->HasToken("TPC_Gas_Size"        )) m_GasSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Gas_Size","mm"));
            if (blocks[i]->HasToken("TPC_Pad_Size"        )) m_PadSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_Pad_Size","mm"));
            if (blocks[i]->HasToken("TPC_MMS_Size"        )) m_MMSSize = NPS::ConvertVector(blocks[i]->GetTVector3("TPC_MMS_Size","mm"));
            if (blocks[i]->HasToken("TPC_Mylar_Rmax"      )) m_Mylar_Rmax = blocks[i]->GetDouble("TPC_Mylar_Rmax","cm");
            if (blocks[i]->HasToken("TPC_Mylar_Thickness" )) m_Mylar_Thickness = blocks[i]->GetDouble("TPC_Mylar_Thickness","micrometer");
            if (blocks[i]->HasToken("Gas_EProductionCut"  )) m_EProductionCut = blocks[i]->GetDouble("Gas_EProductionCut", "mm");
            if (blocks[i]->HasToken("Gas_StepLimit"       )) m_StepLimit = blocks[i]->GetDouble("Gas_StepLimit", "mm");
            if (blocks[i]->HasToken("Gas_Material"        )) GasName = blocks[i]->GetVectorString("Gas_Material");
            if (blocks[i]->HasToken("Gas_Fraction"        )) GasFraction = blocks[i]->GetVectorInt("Gas_Fraction");
            if (blocks[i]->HasToken("Gas_Temperature"     )) m_Temperature = blocks[i]->GetDouble("Gas_Temperature", "kelvin");
            if (blocks[i]->HasToken("Gas_Pressure"        )) m_Pressure = blocks[i]->GetDouble("Gas_Pressure", "bar");
            if (blocks[i]->HasToken("Gas_ReactionBox_Size")) m_ReactionBoxSize = blocks[i]->GetDouble("Gas_ReactionBox_Size","mm");
            if (GasFraction.size()==GasName.size()&&GasName.size()>0) {
                m_GasMaterial.clear();
                m_GasFraction.clear();
                for (unsigned int j = 0; j < GasName.size(); j++) {
                    m_GasMaterial.push_back(GasName[j]);
                    m_GasFraction.push_back(GasFraction[j]);
                }
            }
            vReactionZ = blocks[i]->GetVectorInt("Gas_ReactionZ");
            if (vReactionZ.size()==2) {
                m_ReactionZ1 = vReactionZ[0];
                m_ReactionZ2 = vReactionZ[1];
            }
            if (m_ReactionZ2==m_ReactionZ1) continue;
            else {
                if (m_ReactionZ2<m_ReactionZ1) {
                    m_ReactionZ1 = vReactionZ[1];
                    m_ReactionZ2 = vReactionZ[0];
                }
                m_ReactionZWidth = m_ReactionZ2 - m_ReactionZ1;
            }
            */
        }
        else{
            cout << "ERROR: check your input file formatting (TATOMXPhysics::ReadConfiguration)" << endl;
            exit(1);
        }
    }
}


///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::InitSpectra() {
  m_Spectra = new TATOMXSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::FillSpectra() {
  //m_Spectra -> FillRawSpectra(m_EventData);
  //m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  //m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TATOMXPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("ATOMX", "D"+ NPL::itoa(i+1)+"_ENERGY","ATOMX_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("ATOMX", "D"+ NPL::itoa(i+1)+"_TIME","ATOMX_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("ATOMX",  true );
  inputChain->SetBranchAddress("ATOMX", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  //inputChain->SetBranchAddress("ATOMX", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TATOMXPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  //outputTree->Branch("ATOMX", "TATOMXPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TATOMXPhysics::Construct() {
  return (NPL::VDetector*) new TATOMXPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_ATOMX{
  public:
    proxy_ATOMX(){
      NPL::DetectorFactory::getInstance()->AddToken("ATOMX","ATOMX");
      NPL::DetectorFactory::getInstance()->AddDetector("ATOMX",TATOMXPhysics::Construct);
    }
};

proxy_ATOMX p_ATOMX;
}

