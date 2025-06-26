/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: elidiano.tronchin@studenti.unipd.it                        *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Dali Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TDali2Physics.h"

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
#include "TLorentzVector.h"

ClassImp(TDali2Physics)


///////////////////////////////////////////////////////////////////////////
TDali2Physics::TDali2Physics()
: m_EventData(new TDali2Data),
  m_PreTreatedData(new TDali2Data),
  m_EventPhysics(this),
  m_Spectra(0),
  m_E_RAW_Threshold(0), // adc channels
  m_E_Threshold(0),     // MeV
  m_NumberOfDetectors(0),
  nhit(0) {
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
//void TDaliPhysics::AddDetector(TVector3 , string ){
void TDali2Physics::AddDetector(TVector3 POS) {
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_R.push_back(POS.Perp());
  m_Alpha.push_back(POS.Phi());
  m_Zeta.push_back(POS.Z());
  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TDali2Physics::AddDetector(double R, double Theta, double Phi) { //, string shape){
  // Compute the TVector3 corresponding
  //  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  //  AddDetector(Pos);} 
  double m_r, m_alpha, m_zeta;
  m_r = R*cos(Phi);
  m_alpha = Theta;
  m_zeta = R*sin(Phi);
  m_R.push_back(m_r);
  m_Alpha.push_back(m_alpha);
  m_Zeta.push_back(m_zeta);
  m_NumberOfDetectors++;} 

///////////////////////////////////////////////////////////////////////////
void TDali2Physics::AddDetector(double  R, double  Alpha, double  Zeta, int Ring){
  m_R.push_back(R);
  m_Alpha.push_back(Alpha);
  m_Zeta.push_back(Zeta);
  m_Ring.push_back(Ring);
  m_NumberOfDetectors++;} 


///////////////////////////////////////////////////////////////////////////
void TDali2Physics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TDali2Physics::BuildPhysicalEvent() {
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

  // filling tree
  nhit = DetectorNumber.size();
  for (Int_t i = 0 ; i < nhit ; i++) {
    detN[i] = DetectorNumber[i];
    energy[i] = Energy[i];
    energyDC[i] = Energy[i];
  }
}

double TDali2Physics::GetDopplerCorrectedEnergy(double energy,
						TVector3 position,
						TVector3 beta){
  // renorm pos vector
  TLorentzVector m_GammaLV;
  position.SetMag(1);
  m_GammaLV.SetPx(energy*position.X());
  m_GammaLV.SetPy(energy*position.Y());
  m_GammaLV.SetPz(energy*position.Z());
  m_GammaLV.SetE(energy);
  m_GammaLV.Boost(-beta);
  return m_GammaLV.Energy();}


///////////////////////////////////////////////////////////////////////////
void TDali2Physics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // Energy
  unsigned int mysize = m_EventData->GetMultEnergy();
  for (UShort_t i = 0; i < mysize ; ++i) {
    
  /* ParticleID.push_back(m_EventData->GetParticleID(i)); */
    if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = Cal->ApplyCalibration("Dali2/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
      }
    }
  }

  // Time 
  mysize = m_EventData->GetMultTime();
  for (UShort_t i = 0; i < mysize; ++i) {
    Double_t Time= Cal->ApplyCalibration("Dali2/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
    m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
  }



}

///////////////////////////////////////////////////////////////////////////
void TDali2Physics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigDali2.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigDali2.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigDali.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigDali2.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigDali2";
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
void TDali2Physics::Clear() {
  DetectorNumber.clear();
  Energy.clear();
  Time.clear();

  nhit = 0;
  /* ParticleID.clear(); */
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Dali2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"}; 
  vector<string> cyli = {"R","Alpha","Zeta"}; 

  for (unsigned int i = 0 ; i < blocks.size() ; i++){ 
    if (blocks[i]->HasTokenList(cart)){ 
      if (NPOptionManager::getInstance()->GetVerboseLevel()) 
	cout << endl << "////  Dali2 " << i+1 <<  endl; 
    
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm"); 
      // string Shape = blocks[i]->GetString("Shape"); 
      AddDetector(Pos);} 
    else if(blocks[i]->HasTokenList(cyli)){ 
      if(NPOptionManager::getInstance()->GetVerboseLevel()) 
	cout << endl << "////  Dali2 " << i+1 <<  endl; 
      double R     = blocks[i]->GetDouble("R","mm"); 
      double Alpha = blocks[i]->GetDouble("Alpha","deg"); 
      double Zeta  = blocks[i]->GetDouble("Zeta","mm"); 
      int    Ring  = blocks[i]->GetInt("Ring"); 
      AddDetector(R,Alpha,Zeta,Ring); } 
    else {
      cout << "ERROR: check your input file formatting " << endl; 
      exit(1); }}}

///////////////////////////////////////////////////////////////////////////
void TDali2Physics::InitSpectra() {
  m_Spectra = new TDali2Spectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TDali2Physics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TDali2Physics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}

//TVector3 TDali2Physics::GetDALIPosition(int DALINbr) {
//  static TVector3 Pos;
//return Pos; }
  // return TVector3(m_R[DALINbr]*cos(m_Alpha[DALINbr]),
//		  m_R[DALINbr]*sin(m_Alpha[DALINbr]),
//		  m_Zeta[DALINbr]);}

//double TDali2Physics::GetDALIPPP() {
//  return 100; }


///////////////////////////////////////////////////////////////////////////
void TDali2Physics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("Dali2", "D"+ NPL::itoa(i+1)+"_ENERGY","Dali_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("Dali2", "D"+ NPL::itoa(i+1)+"_TIME","Dali_D"+ NPL::itoa(i+1)+"_TIME");
  }
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus ("Dali2",  true );
  inputChain->SetBranchAddress("Dali2", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("Dali2", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TDali2Physics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Dali2", "TDali2Physics", &m_EventPhysics);

  outputTree->Branch("nhit", &nhit, "nhit/I");
  outputTree->Branch("detN",   detN,   "detN  [nhit]/I");
  outputTree->Branch("energy", energy, "energy[nhit]/D");
  outputTree->Branch("energyDC", energyDC, "energyDC[nhit]/D");
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TDali2Physics::Construct() {
  return (NPL::VDetector*) new TDali2Physics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Dali2{
  public:
    proxy_Dali2(){
      NPL::DetectorFactory::getInstance()->AddToken("Dali2","Dali2");
      NPL::DetectorFactory::getInstance()->AddDetector("Dali2",TDali2Physics::Construct);
    }
};

proxy_Dali2 p_Dali2;
}

