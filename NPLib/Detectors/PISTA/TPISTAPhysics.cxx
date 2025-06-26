/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TPISTAPhysics.h"

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

ClassImp(TPISTAPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TPISTAPhysics::TPISTAPhysics(){
    EventMultiplicity = 0;
    m_EventData = new TPISTAData;
    m_PreTreatedData = new TPISTAData;
    m_EventPhysics = this;
    m_Spectra = NULL;
    m_E_RAW_Threshold = 0; // adc channels
    m_E_Threshold = 1;     // MeV
    m_NumberOfDetectors = 0;
    m_NumberOfStripsX = 57;
    m_NumberOfStripsY = 91;
    m_MaximumStripMultiplicityAllowed = 10;
    m_StripEnergyMatching = 0.050;
    m_DistanceBetweenDEandE = 4;
    m_Back_E_Time_min = 0;
    m_Back_E_Time_max = 100000;
    m_CALIB_BACK_E_PER_STRIP= 0;
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TPISTAPhysics::AddDetector(TVector3 A, TVector3 B, TVector3 C, TVector3 D, int version){
  // Front Face 
  //  A------------------------B
  //   *----------------------*                      
  //    *--------------------*
  //     *------------------*
  //      *----------------*
  //       *--------------*
  //        *------------*
  //         D----------C

  //double Height = 61.7; // mm
  double Height = 57.7; // mm
  //double LongBase = 78.1; // mm
  double LongBase = 72.3; // mm
  //double NumberOfStripsX = m_NumberOfStripsX;
  //double NumberOfStripsY = m_NumberOfStripsY;
 
  if(version==2){
    Height= 27;
    LongBase= 56;
    m_NumberOfStripsX= 64;
    m_NumberOfStripsY= 64;
  }
  double StripPitchY = Height/m_NumberOfStripsY; // mm
  double StripPitchX = LongBase/m_NumberOfStripsX; // mm
 
  m_NumberOfDetectors++;
  m_A.push_back(A);
  m_B.push_back(B);
  m_C.push_back(C);
  m_D.push_back(D);

  // Vector u on telescope face paralelle to Y strips
  TVector3 u = B - A;
  u = u.Unit();
  // Vector v on telescope face paralelle to X strips
  TVector3 v = (C+D)*0.5 - (A+B)*0.5;
  v = v.Unit();
  // Vector n normal to detector surface pointing to target
  TVector3 n = -0.25*(A+B+C+D);
  double norm = n.Mag();
  n = n.Unit();

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 Strip_1_1;
  double ContractedStripPitchX = StripPitchX*norm/(norm+m_DistanceBetweenDEandE);
  double ContractedLongBase = m_NumberOfStripsX*ContractedStripPitchX;
  double deltaX = LongBase/2-ContractedLongBase/2;

  //Strip_1_1 = A + u*(StripPitchX / 2.) + v*(StripPitchY / 2.);
  Strip_1_1 = A + u*deltaX + u*(ContractedStripPitchX / 2.) + v*(m_NumberOfStripsY*StripPitchY - StripPitchY / 2.);
  //Strip_1_1 = A + u*(StripPitchX / 2.) + v*(m_NumberOfStripsY*StripPitchY - StripPitchY / 2.);

  TVector3 StripPos;
  for(int i=0; i<m_NumberOfStripsX; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<m_NumberOfStripsY; j++){
      //StripPos = Strip_1_1 + i*u*StripPitchX - j*v*StripPitchY;
      StripPos = Strip_1_1 + i*u*ContractedStripPitchX - j*v*StripPitchY;
      lineX.push_back(StripPos.X());
      //lineX.push_back(StripPos.X()*norm/(norm+4*abs(sin(n.Phi()))));
      lineY.push_back(StripPos.Y());
      //lineY.push_back(StripPos.Y()*norm/(norm+4*abs(cos(n.Phi()))));
      lineZ.push_back(StripPos.Z());
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);

} 

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::AddDetector(double R, double Theta, double Phi, int version){
  m_NumberOfDetectors++;


  double Height = 61.7; // mm
  double LongBase = 78.1; // mm

  if(version==2){
    Height= 27;
    LongBase= 56;
    m_NumberOfStripsX= 64;
    m_NumberOfStripsY= 64;
  }
  
  double StripPitchY = Height/m_NumberOfStripsY; // mm
  double StripPitchX = LongBase/m_NumberOfStripsX; // mm
 

  // Vector U on detector face (parallel to Y strips) Y strips are along X axis
  TVector3 u;
  // Vector V on detector face (parallel to X strips)
  TVector3 v;
  // Vector W normal to detector face (pointing to the back)
  TVector3 W;
  // Vector C position of detector face center
  TVector3 C;
  /*C = TVector3(R*sin(Theta)*cos(Phi),
        R*sin(Theta)*sin(Phi),
        Height*0.5+R*cos(Theta));*/
  C = TVector3(R*sin(Theta)*cos(Phi),
        R*sin(Theta)*sin(Phi),
        R*cos(Theta));


  TVector3 P = TVector3(-cos(Theta)*cos(Phi),
      -cos(Theta)*sin(Phi),
      sin(Theta));

  W = C.Unit();
  u = W.Cross(P);
  v = W.Cross(u);

  u = u.Unit();
  v = v.Unit();

  double norm = C.Mag();

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  //double X, Y, Z;
  TVector3 Strip_1_1;
  double ContractedStripPitchX = StripPitchX*norm/(norm+m_DistanceBetweenDEandE);
  double ContractedLongBase = m_NumberOfStripsX*ContractedStripPitchX;
  double deltaX = LongBase/2-ContractedLongBase/2;

  TVector3 A = C - 0.5*LongBase*u - 0.5*Height*v;
  Strip_1_1 = A + u*deltaX + u*(ContractedStripPitchX / 2.) + v*(m_NumberOfStripsY*StripPitchY - StripPitchY / 2.);

  TVector3 StripPos;
  for(int i=0; i<m_NumberOfStripsX; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<m_NumberOfStripsY; j++){
      //StripPos = Strip_1_1 + i*u*StripPitchX - j*v*StripPitchY;
      StripPos = Strip_1_1 + i*u*ContractedStripPitchX - j*v*StripPitchY;
      lineX.push_back(StripPos.X());
      //lineX.push_back(StripPos.X()*norm/(norm+4*abs(sin(n.Phi()))));
      lineY.push_back(StripPos.Y());
      //lineY.push_back(StripPos.Y()*norm/(norm+4*abs(cos(n.Phi()))));
      lineZ.push_back(StripPos.Z());
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);
} 

///////////////////////////////////////////////////////////////////////////
TVector3 TPISTAPhysics::GetPositionOfInteraction(int DetectorNumber, int StripE, int StripDE){
  
  TVector3 Position = TVector3(GetStripPositionX(DetectorNumber, StripE, StripDE),
      GetStripPositionY(DetectorNumber, StripE, StripDE),
      GetStripPositionZ(DetectorNumber, StripE, StripDE));


  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TPISTAPhysics::GetDetectorNormal(const int i){
  int det = DetectorNumber[i];
  // Vector u on telescope face paralelle to Y strips
  TVector3 u = m_B[det-1] - m_A[det-1];
  u = u.Unit();
  // Vector v on telescope face paralelle to X strips
  TVector3 v = (m_C[det-1] + m_D[det-1])*0.5 - (m_A[det-1] + m_B[det-1])*0.5;
  v = v.Unit();


  TVector3 Normal = u.Cross(v);

  return (Normal.Unit());
}
///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();

  int DEMult = 0;
  int EMult = 0;
  int mult_DE_per_telescope[8];
  int mult_E_per_telescope[8];
  for(int i=0; i<8; i++){
    mult_DE_per_telescope[i] = 0;
    mult_E_per_telescope[i] = 0;
  }


  if(1 /*CheckEvent() == 1*/){
    //vector<TVector2> couple = Match_X_Y();
    //EventMultiplicity = couple.size();

    DEMult = m_PreTreatedData->GetPISTADEMult();
    EMult = m_PreTreatedData->GetPISTAEMult();
    
    /*for(unsigned int i=0; i<DEMult; i++){
      int det = m_PreTreatedData->GetPISTA_DE_DetectorNbr(i);
      mult_DE_per_telescope[det-1]++;
    }
    for(unsigned int i=0; i<EMult; i++){
      int det = m_PreTreatedData->GetPISTA_E_DetectorNbr(i);
      mult_E_per_telescope[det-1]++;
    }*/

    int DE_DetNbr = -1;
    int StripNbr_DE = -1;
    int E_DetNbr = -1;
    int StripNbr_E = -1;

    for(unsigned int i=0; i<DEMult; i++){
      DE_DetNbr = m_PreTreatedData->GetPISTA_DE_DetectorNbr(i);
      StripNbr_DE = m_PreTreatedData->GetPISTA_DE_StripNbr(i);

      mult_DE_per_telescope[DE_DetNbr-1]++;
      // *** to be removed *** //
      if(EMult==0){
        DetectorNumber.push_back(DE_DetNbr);
        DE_StripNbr.push_back(StripNbr_DE);
        DE.push_back(m_PreTreatedData->GetPISTA_DE_StripEnergy(i));
        back_DE.push_back(m_PreTreatedData->GetPISTA_DE_BackEnergy(i));
        DE_Time.push_back(m_PreTreatedData->GetPISTA_DE_StripTime(i));
        back_DE_Time.push_back(m_PreTreatedData->GetPISTA_DE_BackTime(i));

        E_StripNbr.push_back(-100);
        E.push_back(-100);
        E_Time.push_back(-100);
        back_E.push_back(-100);
        back_E_Time.push_back(-100);
      }
      // *** // 

      for(unsigned int j=0; j<EMult; j++){
        E_DetNbr = m_PreTreatedData->GetPISTA_E_DetectorNbr(j);
        StripNbr_E = m_PreTreatedData->GetPISTA_E_StripNbr(j);

        mult_E_per_telescope[E_DetNbr-1]++;
        double back_time = m_PreTreatedData->GetPISTA_E_BackTime(j);
        if(DE_DetNbr==E_DetNbr && back_time>m_Back_E_Time_min && back_time<m_Back_E_Time_max){
          // Taking Strip energy for DE
          double DE_Energy = m_PreTreatedData->GetPISTA_DE_StripEnergy(i);
          // Taking BAck Energy for E
          double E_Energy = m_PreTreatedData->GetPISTA_E_StripEnergy(j);

          DetectorNumber.push_back(DE_DetNbr);
          DE_StripNbr.push_back(StripNbr_DE);
          E_StripNbr.push_back(StripNbr_E);
          DE.push_back(DE_Energy);
          E.push_back(E_Energy);
          back_DE.push_back(m_PreTreatedData->GetPISTA_DE_BackEnergy(i));
          back_E.push_back(m_PreTreatedData->GetPISTA_E_BackEnergy(j));

          DE_Time.push_back(m_PreTreatedData->GetPISTA_DE_StripTime(i));
          back_DE_Time.push_back(m_PreTreatedData->GetPISTA_DE_BackTime(i));
          E_Time.push_back(m_PreTreatedData->GetPISTA_E_StripTime(j));
          back_E_Time.push_back(m_PreTreatedData->GetPISTA_E_BackTime(j));

          if(StripNbr_DE>0 && StripNbr_DE<m_NumberOfStripsY+1 && StripNbr_E>0 && StripNbr_E<m_NumberOfStripsX+1){
            PosX.push_back(GetPositionOfInteraction(DE_DetNbr, StripNbr_E, StripNbr_DE).X());
            PosY.push_back(GetPositionOfInteraction(DE_DetNbr, StripNbr_E, StripNbr_DE).Y());
            PosZ.push_back(GetPositionOfInteraction(DE_DetNbr, StripNbr_E, StripNbr_DE).Z());
          }
        }
      }
    }
  }
  for(int i=0; i<8; i++){
    mult_DE.push_back(mult_DE_per_telescope[i]);
    mult_E.push_back(mult_E_per_telescope[i]);
  }


  if(EMult==0)
    EventMultiplicity=0;
  else
    EventMultiplicity = DetectorNumber.size();

  if(EMult>0 && DEMult>0 && EventMultiplicity==0){
    m_DE_E_NoMatching++;
  }
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TPISTAPhysics::Match_X_Y(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_DEMult, m_EMult;
  m_DEMult = m_PreTreatedData->GetPISTADEMult();
  m_EMult = m_PreTreatedData->GetPISTAEMult();

  if(m_DEMult>m_MaximumStripMultiplicityAllowed || m_EMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  return ArrayOfGoodCouple;
}

///////////////////////////////////////////////////////////////////////////
int TPISTAPhysics::CheckEvent(){
  // Check the size of the different elements
  if(m_PreTreatedData->GetPISTADEMult() == m_PreTreatedData->GetPISTAEMult() )
    return 1;

  else
    return -1;
}


///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  //////
  // DE
  unsigned int sizeDE = m_EventData->GetPISTADEMult();
  unsigned int sizeDE_back = m_EventData->GetPISTADEBackMult();
  for (UShort_t i = 0; i < sizeDE ; ++i) {
    if (IsValidChannel(0,m_EventData->GetPISTA_DE_DetectorNbr(i),m_EventData->GetPISTA_DE_StripNbr(i)) && m_EventData->GetPISTA_DE_StripEnergy(i)<4095) {
      int DetNbr = m_EventData->GetPISTA_DE_DetectorNbr(i);
      int StripNbr = m_EventData->GetPISTA_DE_StripNbr(i);
      double StripE = m_EventData->GetPISTA_DE_StripEnergy(i);
      double StripT  = m_EventData->GetPISTA_DE_StripTime(i);

      double ped = Cal->GetValue("PISTA/T"+NPL::itoa(DetNbr)+"_STRIP"+NPL::itoa(StripNbr)+"_DE_PEDESTAL",0);
      double CalStripE = Cal->ApplyCalibration("PISTA/T"+NPL::itoa(DetNbr)+"_STRIP"+NPL::itoa(StripNbr)+"_DE_ENERGY",StripE);

      if(sizeDE_back==0 && CalStripE > m_E_Threshold){
        m_PreTreatedData->SetPISTA_DE(DetNbr, StripNbr, CalStripE, -100, StripT, -100);
      }

      for(UShort_t j = 0; j< sizeDE_back; j++){
        double BackDE = m_EventData->GetPISTA_DE_BackEnergy(j);
        double BackT  = m_EventData->GetPISTA_DE_BackTime(j);
        int BackDet   = m_EventData->GetPISTA_DE_BackDetector(j); 

        double CalBackDE = Cal->ApplyCalibration("PISTA/T"+NPL::itoa(DetNbr)+"_BACK_DE",BackDE);
        if (CalStripE > m_E_Threshold && DetNbr==BackDet) {
          m_PreTreatedData->SetPISTA_DE(DetNbr, StripNbr, CalStripE, CalBackDE, StripT, BackT);
        }
      }
    }
  }
  if(sizeDE_back>0 && sizeDE>0 && m_PreTreatedData->GetPISTADEBackMult()==0){
    m_DE_NoMatching++; 
  }

  // E
  unsigned int sizeE = m_EventData->GetPISTAEMult();
  unsigned int sizeE_back = m_EventData->GetPISTAEBackMult();
  for (UShort_t i = 0; i < sizeE ; ++i) {
    if (IsValidChannel(1,m_EventData->GetPISTA_E_DetectorNbr(i),m_EventData->GetPISTA_E_StripNbr(i)) && m_EventData->GetPISTA_E_StripEnergy(i)<4095) {
      int DetNbr = m_EventData->GetPISTA_E_DetectorNbr(i);
      int StripNbr = m_EventData->GetPISTA_E_StripNbr(i);
      double StripE = m_EventData->GetPISTA_E_StripEnergy(i);
      double StripT  = m_EventData->GetPISTA_E_StripTime(i);

      double ped = Cal->GetValue("PISTA/T"+NPL::itoa(DetNbr)+"_STRIP"+NPL::itoa(StripNbr)+"_E_PEDESTAL",0);
      double CalStripE = Cal->ApplyCalibration("PISTA/T"+NPL::itoa(DetNbr)+"_STRIP"+NPL::itoa(StripNbr)+"_E_ENERGY",StripE);
      
      for(UShort_t j = 0; j< sizeE_back; j++){
        double BackE  = m_EventData->GetPISTA_E_BackEnergy(j);
        double BackT  = m_EventData->GetPISTA_E_BackTime(j);
        int BackDet = m_EventData->GetPISTA_E_BackDetector(j);

        double CalBackE;
        if(m_CALIB_BACK_E_PER_STRIP==1){
          CalBackE = Cal->ApplyCalibration("PISTA/T"+NPL::itoa(DetNbr)+"_BACK_E_STRIP"+NPL::itoa(StripNbr),BackE);
        }
        else{
          CalBackE = Cal->ApplyCalibration("PISTA/T"+NPL::itoa(DetNbr)+"_BACK_E",BackE);
        }

        if (CalStripE > m_E_Threshold && DetNbr==BackDet) {
          m_PreTreatedData->SetPISTA_E(DetNbr, StripNbr, CalStripE, CalBackE, StripT, BackT);
        }
      }
    }
  }

  if(sizeE_back>0 &&  sizeE>0 && m_PreTreatedData->GetPISTAEBackMult()==0){
    m_E_NoMatching++; 
  }
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigPISTA.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigPISTA.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigPISTA.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigPISTA.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigPISTA";
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

      else if (whatToDo=="BACK_E_TIME_MIN") {
        AnalysisConfigFile >> DataBuffer;
        m_Back_E_Time_min = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Back_E_Time_min << endl;
      }

      else if (whatToDo=="BACK_E_TIME_MAX") {
        AnalysisConfigFile >> DataBuffer;
        m_Back_E_Time_max = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Back_E_Time_max << endl;
      }
      else if (whatToDo=="CALIB_BACK_E_PER_STRIP") {
        AnalysisConfigFile >> DataBuffer;
        m_CALIB_BACK_E_PER_STRIP = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_CALIB_BACK_E_PER_STRIP << endl;
      }
      else if(whatToDo=="DISABLE_CHANNEL"){
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << " " << DataBuffer << endl;
        int telescope = atoi(DataBuffer.substr(5,1).c_str());
        int channel = -1;
        if(DataBuffer.compare(6,3,"_DE") == 0 ){
          channel = atoi(DataBuffer.substr(9).c_str());
          *(m_DEChannelStatus[telescope -1].begin() + channel -1) = false;
        }
        else if(DataBuffer.compare(6,2,"_E") == 0 ){
          channel = atoi(DataBuffer.substr(8).c_str());
          *(m_EChannelStatus[telescope -1].begin() + channel -1) = false;
        }

      }

      else {
        ReadingStatus = false;
      }
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::Clear() {
  EventMultiplicity = 0;
  m_E_NoMatching = 0;
  m_DE_NoMatching = 0;
  m_DE_E_NoMatching = 0;


  // Position Information
  PosX.clear();
  PosY.clear();
  PosZ.clear();

  DetectorNumber.clear();
  E.clear();
  DE.clear();
  back_DE.clear();
  back_E.clear();
  DE_StripNbr.clear();
  E_StripNbr.clear();
  DE_Time.clear();
  back_DE_Time.clear();
  E_Time.clear();
  back_E_Time.clear();

  mult_DE.clear();
  mult_E.clear();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("PISTA");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS_A","POS_B","POS_C","POS_D"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;

      int version = blocks[i]->GetInt("VERSION");
      TVector3 A = blocks[i]->GetTVector3("POS_A","mm");
      TVector3 B = blocks[i]->GetTVector3("POS_B","mm");
      TVector3 C = blocks[i]->GetTVector3("POS_C","mm");
      TVector3 D = blocks[i]->GetTVector3("POS_D","mm");

      AddDetector(A,B,C,D,version);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  PISTA " << i+1 <<  endl;
      
      int version = blocks[i]->GetInt("VERSION");
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");

      AddDetector(R, Theta, Phi, version);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  InitializeStandardParameter();
  ReadAnalysisConfig();
}

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeStandardParameter() {

  // Enable all channel
  vector<bool> ChannelStatusDE;
  vector<bool> ChannelStatusE;
  m_DEChannelStatus.clear();
  m_EChannelStatus.clear();

  ChannelStatusDE.resize(m_NumberOfStripsY,true);
  ChannelStatusE.resize(m_NumberOfStripsX,true);

  for(int i=0; i<m_NumberOfStripsY; i++){
    m_DEChannelStatus[i] = ChannelStatusDE;
  }
  for(int i=0; i<m_NumberOfStripsX; i++){
    m_EChannelStatus[i] = ChannelStatusE;
  }

}

///////////////////////////////////////////////////////////////////////////
bool TPISTAPhysics::IsValidChannel(const int& DetectorType, const int& telescope, const int& channel){

  if(DetectorType==0){
    return *(m_DEChannelStatus[telescope - 1].begin() + channel -1);
  }
  else if(DetectorType==1){
    return *(m_EChannelStatus[telescope - 1].begin() + channel -1);
  }

  else 
    return false;

}
///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitSpectra() {
  m_Spectra = new TPISTASpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TPISTAPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for(int i=0; i<m_NumberOfDetectors; ++i) {
    Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_BACK_DE","PISTA_T"+ NPL::itoa(i+1)+"_BACK_DE");
    Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_BACK_E","PISTA_T"+ NPL::itoa(i+1)+"_BACK_E");

    for(int j=0; j<m_NumberOfStripsY; j++){
      Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_DE_ENERGY","PISTA_T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_DE_ENERGY");
      Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_DE_PEDESTAL","PISTA_T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_DE_PEDESTAL");
    }

    for(int j=0; j<m_NumberOfStripsX; j++){
      Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_E_ENERGY","PISTA_T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_E_ENERGY");
      Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_E_PEDESTAL","PISTA_T"+ NPL::itoa(i+1)+"_STRIP"+NPL::itoa(j+1)+"_E_PEDESTAL");
      Cal->AddParameter("PISTA", "T"+ NPL::itoa(i+1)+"_BACK_E_STRIP"+NPL::itoa(j+1),"PISTA_T"+ NPL::itoa(i+1)+"_BACK_E_STRIP"+NPL::itoa(j+1));
    }


  }
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("PISTA",  true );
  inputChain->SetBranchAddress("PISTA", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("PISTA", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TPISTAPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("PISTA", "TPISTAPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TPISTAPhysics::Construct() {
  return (NPL::VDetector*) new TPISTAPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_PISTA{
    public:
      proxy_PISTA(){
        NPL::DetectorFactory::getInstance()->AddToken("PISTA","PISTA");
        NPL::DetectorFactory::getInstance()->AddDetector("PISTA",TPISTAPhysics::Construct);
      }
  };

  proxy_PISTA p_PISTA;
}

