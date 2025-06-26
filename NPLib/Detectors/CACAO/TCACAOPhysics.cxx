/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jwhwang@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 4월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CACAO Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TCACAOPhysics.h"

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

ClassImp(TCACAOPhysics)


///////////////////////////////////////////////////////////////////////////
TCACAOPhysics::TCACAOPhysics()
   : m_EventData(new TCACAOData),
     m_PreTreatedData(new TCACAOData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void TCACAOPhysics::AddDetector(TVector3 Pos, TRotation Rot,
				TVector3 Dim, Double_t ShieldThickness) {
  m_Pos.push_back(Pos);
  m_Rot.push_back(Rot);
  m_Dim.push_back(Dim);
  m_ShieldThickness.push_back(ShieldThickness);
  m_NumberOfDetectors++;
}


///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::BuildPhysicalEvent() {
  Clear();
  
  // apply thresholds and calibration
  // PreTreat();

  ////////////////////////////////////////////////////////////
  // Loop for All
  nhit = m_EventData->GetMult();
  for (Int_t i = 0 ; i < m_EventData->GetMult() ; i++) {
    detN[i] = m_EventData->GetDetN(i);
    E   [i] = m_EventData->GetE(i);
    T   [i] = m_EventData->GetT(i);
  }
  ////////////////////////////////////////////////////////////
}

///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigCACAO.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigCACAO.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigCACAO.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigCACAO.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigCACAO";
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
void TCACAOPhysics::Clear() {
  nhit = 0; }



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("CACAO");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> reso = {"Reso"};
  
  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};
  vector<string> cyli = {"Rho","Phi","Z"};

  vector<string> cuboid= {"DIM","ShieldThickness"};


  for(unsigned int i = 0 ; i < blocks.size() ; i++){

    ////////////////////////////////////////////////////////////
    // Resolution
    //    if (blocks[i]->HasTokenList(reso)) {
    //      CACAO_NS::ResoEnergy = blocks[i]->GetDouble("Reso","void");
    //      continue; }
    ////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////
    // checking the position items of the block
    TVector3 Pos;
    if (blocks[i]->HasTokenList(cart)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO " << i+1 <<  endl;
      Pos = blocks[i]->GetTVector3("POS","mm"); }
    else if (blocks[i]->HasTokenList(sphe)) {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO " << i+1 <<  endl;
      Pos.SetMagThetaPhi(blocks[i]->GetDouble("R","mm"),
			 blocks[i]->GetDouble("Theta","deg"),
			 blocks[i]->GetDouble("Phi","deg")); }
    else if (blocks[i]->HasTokenList(cyli)) {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CACAO " << i+1 <<  endl;
      Pos.SetPtThetaPhi(blocks[i]->GetDouble("Rho","mm"), 0,  blocks[i]->GetDouble("Phi","deg"));
      Pos.SetZ(blocks[i]->GetDouble("Z","mm")); }
    else {
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);}

    TRotation Rot;
    if (blocks[i]->HasToken("ANG")) {
      TVector3 Ang;
      Ang = blocks[i]->GetTVector3("ANG","deg");
      Rot.SetXEulerAngles(Ang.x(), Ang.y(), Ang.z()); } // (phi, theta, psi)
    ////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////
    // checking the shape items of the block
    TVector3 Dim;
    if(blocks[i]->HasTokenList(cuboid)){
      Dim = blocks[i]->GetTVector3("DIM","mm");
      double ShieldThickness = blocks[i]->GetDouble("ShieldThickness","mm");
      AddDetector(Pos, Rot, Dim, ShieldThickness); }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
    ////////////////////////////////////////////////////////////
    
  }
}

///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::InitSpectra() {
  m_Spectra = new TCACAOSpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TCACAOPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("CACAO",  true );
  inputChain->SetBranchAddress("CACAO", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("CACAO", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TCACAOPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("CACAO", "TCACAOPhysics", &m_EventPhysics);
  outputTree->Branch("nhit" ,&nhit,"nhit/I");
  outputTree->Branch("detN" ,detN, "detN[nhit]/I");
  outputTree->Branch("E",E,"E[nhit]/D");
  outputTree->Branch("T",T,"T[nhit]/D");

}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TCACAOPhysics::Construct() {
  return (NPL::VDetector*) new TCACAOPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_CACAO{
  public:
    proxy_CACAO(){
      NPL::DetectorFactory::getInstance()->AddToken("CACAO","CACAO");
      NPL::DetectorFactory::getInstance()->AddDetector("CACAO",TCACAOPhysics::Construct);
      //NPL::DetectorFactory::getInstance()->AddDetectorReader("CACAO",TCACAOPhysics::ConstructReader);
    }
};

proxy_CACAO p_CACAO;
}

