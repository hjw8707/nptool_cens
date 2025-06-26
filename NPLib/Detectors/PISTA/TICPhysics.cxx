/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace  contact address: pierre.morfouace@cea.fr   *
 *                                                                           *
 * Creation Date  : October 2023                                             *
 * Last update    : 22/01/24                                                 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold IC Treated  data                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TICPhysics.h"

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
#include "TKey.h"

ClassImp(TICPhysics);

///////////////////////////////////////////////////////////////////////////
TICPhysics::TICPhysics()
  : m_EventData(new TICData),
  m_PreTreatedData(new TICData),
  m_EventPhysics(this),
  m_FPMW_Section(-1),
  m_Eres_Threshold(3000),
  m_Z_SPLINE_CORRECTION(false),
  m_NumberOfDetectors(0){
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TICPhysics::AddDetector(TVector3 Pos){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
} 

///////////////////////////////////////////////////////////////////////////
void TICPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TICPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::BuildPhysicalEvent() {

  if(m_FPMW_Section<0)
    return;

  Clear();
  PreTreat();

  static CalibrationManager* Cal = CalibrationManager::getInstance();

  // fIC definition
  int size = m_PreTreatedData->GetICMult();
  for(int i=0; i<size; i++){
    int segment = m_PreTreatedData->GetIC_Section(i);
    double gain = Cal->GetValue("IC/SEC"+NPL::itoa(m_FPMW_Section)+"_SEG"+NPL::itoa(segment)+"_ALIGN",0);
    double GainInit = Cal->GetValue("IC/INIT_SEG"+NPL::itoa(segment)+"_ALIGN",0);

    fIC_raw[i] = m_PreTreatedData->GetIC_Charge(i);

    fIC[i] = gain*m_PreTreatedData->GetIC_Charge(i);

    if(i < 4){
      fIC_PID[i] = 0.5* m_PreTreatedData->GetIC_Charge(i);
    }
    else if(i >= 6 && m_Data_Year==2024){
      fIC_PID[i] = 2*m_PreTreatedData->GetIC_Charge(i);
    }
    else if((i==4 || i==5) && m_Data_Year==2024){
      fIC_PID[i] = m_PreTreatedData->GetIC_Charge(i);
    }
    else if(i >= 4 && m_Data_Year==2023){
      fIC_PID[i] = m_PreTreatedData->GetIC_Charge(i);
    }

    //fIC for different Tof in vamos
    // Order of ToF is 13 23 14 24
    for (int ToF=0 ; ToF < 4 ; ToF ++){
      double gainToF = Cal->GetValue("IC/SEC"+NPL::itoa(m_FPMW_Section)+"_SEG"+NPL::itoa(segment)+"_TOF"+ToFName[ToF]+"_ALIGN",0) ;
      fIC_AoQ[ToF][i] = gainToF * m_PreTreatedData->GetIC_Charge(i);
    }

    fIC_Init[i] = GainInit * m_PreTreatedData->GetIC_Charge(i);
    fIC_TS.push_back(m_PreTreatedData->GetIC_TS(i));
  }

  if(m_Y_SPLINE_CORRECTION && m_XY0_SPLINE_CORRECTION){
    ApplyXYCorrections(); 
  }

  if (fIC_Init[1]>0 && fIC_Init[5]>0) {
    EtotInit = 0 ;
    double ScalorInit = Cal->GetValue("IC/INIT_ETOT_SCALING",0);

    for (int Seg = 0; Seg<10; Seg++){

      if (Seg == 0 && m_Data_Year == 2024 ){
        EtotInit += 0.75 * Cal->GetValue("IC/INIT_SEG"+NPL::itoa(Seg+1)+"_ALIGN",0)* fIC_raw[1];
      }

      else { 
        EtotInit += fIC_Init[Seg];  
      }
    } 
    EtotInit = EtotInit * ScalorInit ;
  } //Condition for Init Etot

  else {
    EtotInit = -100 ;
  }

  DE = fIC_PID[0] + fIC_PID[1] + fIC_PID[2] + fIC_PID[3] + fIC_PID[4];
  //DE = 0.5*(fIC_raw[0] + fIC_raw[1] + fIC_raw[2] + fIC_raw[3]) + fIC_raw[4];
  Eres = fIC_PID[5] + fIC_PID[6] + fIC_PID[7] + fIC_PID[8] + fIC_PID[9];

  if (m_DE_SPLINE_CORRECTION && m_Y_SPLINE_CORRECTION){
    if (m_TimeData->GetMWPC13Mult() ==1 && fIC_TS.size()>=8 ){ //only mult 1 event
      UShort_t FPMW_Section = m_FPMW_Section;
      double TempY;

      //Data year sensitive loading
      if (m_Data_Year == 2024){
        TempY  =  10* (fIC_TS.at(0) - m_TimeData->GetTS_MWPC13(0)) - ((m_TimeData->GetTime_MWPC13(0)+m_TimeData->GetToff_DT13(FPMW_Section))) ;
      }
      else if (m_Data_Year == 2023){
        TempY =  m_Y ;
      }
      DE = DE * m_DEspline.at(0)->Eval(0) / m_DEspline.at(0)->Eval(TempY) ;
    } // end if mult 1
  } // end DE correction

  if(fIC[1]>0 && fIC[5]>0){
    double scalor = Cal->GetValue("IC/ETOT_SCALING_SEC"+NPL::itoa(m_FPMW_Section),0);

    for(int i=0; i<10; i++){
      if(i == 0 && m_Data_Year == 2024){
        Etot += Cal->GetValue("IC/SEC"+NPL::itoa(m_FPMW_Section)+"_SEG"+NPL::itoa(1)+"_ALIGN",0)* fIC_raw[1];
      } 
      else{
        Etot += fIC[i];
      }
    }
    Etot = scalor*Etot;


    //Etot AoQ
    for (int ToF=0 ; ToF < 4 ; ToF ++){
      double scalorToF = Cal->GetValue("IC/ETOT_SCALING_SEC"+NPL::itoa(m_FPMW_Section)+"_TOF"+ToFName[ToF],0);

      for(int i=0; i<10; i++){
        if(i == 0 && m_Data_Year == 2024){
          EtotAoQ[ToF] += Cal->GetValue("IC/SEC"+NPL::itoa(m_FPMW_Section)+"_SEG"+NPL::itoa(1)+"_TOF"+ToFName[ToF]+"_ALIGN",0)* fIC_raw[1];
        } 
        else{
          EtotAoQ[ToF] += fIC_AoQ[ToF][i];
        }
      }
      EtotAoQ[ToF] = scalorToF*EtotAoQ[ToF];
    }// End etot AoQ

    //Etot = 0.02411*(0.8686*fIC_raw[0]+0.7199*fIC_raw[1]+0.6233*fIC_raw[2]+0.4697*fIC_raw[3]+0.9787*fIC_raw[4]+0.9892*fIC_raw[5]+2.1038*fIC_raw[6]+1.9429*fIC_raw[7]+1.754*fIC_raw[8]+2.5*fIC_raw[9]);  


    if(m_Z_SPLINE_CORRECTION && Eres>m_Eres_Threshold){

      if (m_Z_THETA_CORRECTION) Chio_ZRaw = m_Z_Theta_spline.at(0)->Eval(0)/m_Z_Theta_spline.at(0)->Eval(m_Thetaf) * ApplyZSpline();
      else Chio_ZRaw = ApplyZSpline();
      Chio_Z = Cal->ApplyCalibration("IC/Z_CALIBRATION",Chio_ZRaw);
    }
    else {
      Chio_Z = -1000;
      Chio_ZRaw = -1000;
    };

  }
  else{
    DE = -100;
    Eres = -100;
    Etot = -100;
    Chio_Z = -1000;
    Chio_ZRaw = -1000;
  }

  m_FPMW_Section = -1;
}

///////////////////////////////////////////////////////////////////////////
void TICPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  //static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysize = m_EventData->GetICMult();
  for (unsigned int i = 0; i < mysize ; ++i) {
    int segment = m_EventData->GetIC_Section(i);
    //cout << section << " " << gain << endl;
    //double charge = gain*m_EventData->GetIC_Charge(i);
    double charge = m_EventData->GetIC_Charge(i);
    long TS = m_EventData->GetIC_TS(i);


    m_PreTreatedData->SetIC_Charge(charge);
    m_PreTreatedData->SetIC_Section(segment);
    m_PreTreatedData->SetIC_TS(TS);
  }

  // Find the highest section with a charge
  unsigned int secHighCharge = 0 ;
  double chargeHighest = 0 ;
  if (mysize > 1){
    for (unsigned int i = (mysize-1) ; i >= 0 ; i--){
      chargeHighest = m_EventData->GetIC_Charge(i);
      secHighCharge = i ;
      break;
    } // end if charge 
  } // end inverse loop

  if ((chargeHighest > 0) && (secHighCharge!=0)){
    // loop through all section until we get to the highest
    // if one of them is empty we clear
    for (unsigned int j = 0 ; j<secHighCharge ; j++ ){
      double charge = m_EventData->GetIC_Charge(j);
      if (charge <= 0) {
        ClearPreTreatedData();
        break;
      } // end loop 
    } // end cond on charge 


  }
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigIC.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigIC.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigIC.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigIC.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigIC";
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

      else if (whatToDo=="ERESIDUAL_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Eres_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Eres_Threshold << endl;
      }

      else if (whatToDo=="DATA_YEAR") {
        AnalysisConfigFile >> DataBuffer;
        m_Data_Year = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Data_Year << endl;
      }
      else if (whatToDo=="LOAD_Y_SPLINE") {
        AnalysisConfigFile >> DataBuffer;
        m_Y_SPLINE_PATH = DataBuffer;
        cout << "*** Loading Y spline ***" << endl;
        m_Y_SPLINE_CORRECTION = LoadSpline(m_Yspline,m_number_Y_spline,m_Y_SPLINE_PATH);
      }
      else if (whatToDo=="LOAD_XY0_PROFILE") {
        AnalysisConfigFile >> DataBuffer;
        m_XY0_PROFILE_PATH = DataBuffer;
        cout << "*** Loading XY0 profile ***" << endl;
        TString PathTemp = m_XY0_PROFILE_PATH;
        TFile* ifile = new TFile(PathTemp,"read");
        if(ifile->IsOpen() && !ifile->IsZombie()){
          m_XY0_SPLINE_CORRECTION = true;
          m_IC0_Profile.LoadProfile(PathTemp,"ICOneZeroProfile");
        }
        else {
          m_XY0_SPLINE_CORRECTION = false ;
        }
        ifile->Close();
      }
      else if (whatToDo=="LOAD_DE_SPLINE") {
        AnalysisConfigFile >> DataBuffer;
        m_DE_SPLINE_PATH = DataBuffer;
        cout << "*** Loading DE spline ***" << endl;
        m_DE_SPLINE_CORRECTION = LoadSpline(m_DEspline,m_number_DE_spline,m_DE_SPLINE_PATH);
      }
      else if (whatToDo=="LOAD_Z_SPLINE"){
        AnalysisConfigFile >> DataBuffer;
        m_Z_SPLINE_PATH = DataBuffer;
        cout << "*** Loading Z spline ***" << endl;
        m_Z_SPLINE_CORRECTION = LoadSpline(m_Zspline,m_number_zspline,m_Z_SPLINE_PATH);
      }

      else if (whatToDo=="LOAD_Z_SPLINE_EVAL"){
        AnalysisConfigFile >> DataBuffer;
        m_Z_SPLINE_EVAL_PATH = DataBuffer;
        cout << "*** Loading Z spline eval position***" << endl;
        m_Z_SPLINE_EVAL = LoadVector(m_Z_spline_eval,m_Z_SPLINE_EVAL_PATH);
      }

      else if (whatToDo=="LOAD_Z_THETA_SPLINE") {
        AnalysisConfigFile >> DataBuffer;
        m_Z_THETA_SPLINE_PATH = DataBuffer;
        cout << "*** Loading Z THETA spline ***" << endl;
        m_Z_THETA_CORRECTION = LoadSpline(m_Z_Theta_spline,m_number_Z_THETA_spline,m_Z_THETA_SPLINE_PATH);
      }
      else {
        ReadingStatus = false;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////
bool TICPhysics::LoadSpline(vector<TSpline3*> &iSpline, int &NSpline , string Path){
  TString filename = Path;
  TFile* ifile = new TFile(filename,"read");
  if(ifile->IsOpen() && !ifile->IsZombie()){

    // Get number of spline
    TIter next(ifile->GetListOfKeys());
    TKey* key;
    NSpline = 0 ;

    while ((key=(TKey*)next())){
      if (std::string(key->GetClassName()) == "TSpline3"){
        NSpline ++;
      }
    }
    cout << "This file contains  " << NSpline << " splines "  << endl;
    // Load Spline
    for(int i=0; i<NSpline; i++){
      iSpline.at(i) = (TSpline3*) ifile->FindObjectAny(Form("fspline_%d",i+1));
      iSpline.at(i)->SetName(Form("fspline_%s_%d",Path.c_str(),i+1));
      cout  << iSpline.at(i)->GetName() << " is loaded!" << endl;
    }

    ifile->Close();
    return true;
  }
  else{
    cout << "File " << filename << " not found!" << endl;
    ifile->Close();
    return false;
  }
}

///////////////////////////////////////////////////////////////////////////
template <typename T> bool TICPhysics::LoadVector(vector<T> &vec, const string &Path){    

  std::ifstream file(Path);
  if (!file.is_open()) {
    return false; // File couldn't be opened
  }

  T value;
  vec.clear();
  std::string line;
  while (std::getline(file, line)) {
    std::istringstream iss(line);
    while (iss >> value) {
      vec.push_back(value);
    }
  }
  file.close();
  return !vec.empty();
}


///////////////////////////////////////////////////////////////////////////
double TICPhysics::ApplyZSpline(){
  double DEcorr;
  double FF_DEcorr0[m_number_zspline];

  if (m_Z_SPLINE_EVAL == true){
    for(int i=0; i<m_number_zspline; i++){
      if (i < m_Z_spline_eval.size()){
        FF_DEcorr0[i] = m_Zspline[i]->Eval(m_Z_spline_eval[i]);
      }
      else if( m_Data_Year == 2024 ){
        FF_DEcorr0[i] = m_Zspline[i]->Eval(12000);
      }
      else if( m_Data_Year == 2023 ){
        FF_DEcorr0[i] = m_Zspline[i]->Eval(8500);
      }
    } // end loop spline
  } // end cond spline eval
  else if (m_Data_Year == 2023){
    for(int i=0; i<m_number_zspline; i++){
      FF_DEcorr0[i] = m_Zspline[i]->Eval(8500);
    }
  }
  else if (m_Data_Year == 2024){
    for(int i=0; i<m_number_zspline; i++){
      FF_DEcorr0[i] = m_Zspline[i]->Eval(12000);
    }
  }
  double DEspline0;
  double Eval_DEspline;
  int index=0;
  for(int i=0; i<m_number_zspline; i++){
    Eval_DEspline = m_Zspline[i]->Eval(Eres);
    if(DE<Eval_DEspline) break;
    index = i;
  }

  Eval_DEspline = m_Zspline[index]->Eval(Eres);
  DEspline0 = FF_DEcorr0[index];
  double dmin, dsup;
  if(index<(m_number_zspline-1) && DE>m_Zspline[0]->Eval(Eres)){
    dmin = DE - m_Zspline[index]->Eval(Eres);
    dsup = m_Zspline[index+1]->Eval(Eres) - DE;

    Eval_DEspline = dsup*m_Zspline[index]->Eval(Eres)/(dmin+dsup) + dmin*m_Zspline[index+1]->Eval(Eres)/(dmin+dsup);
    DEspline0     = dsup*FF_DEcorr0[index]/(dmin+dsup) + dmin*FF_DEcorr0[index+1]/(dmin+dsup);

    DEcorr        = DEspline0 * DE / Eval_DEspline;
  }
  else if(index==m_number_zspline-1){
    Eval_DEspline = m_Zspline[index]->Eval(Eres);
    DEspline0     = FF_DEcorr0[index];

    DEcorr        = DEspline0 * DE / Eval_DEspline;
  }

  return DEcorr;
}
///////////////////////////////////////////////////////////////////////////
void TICPhysics::ApplyXYCorrections(){
  vector<double> ICcorr_Y(11), ICcorr_X(11); 
  double FF_DriftTime ;
  if (m_TimeData->GetMWPC13Mult() ==1 && fIC_TS.size()>=8){ //only mult 1 event
    UShort_t FPMW_Section = m_FPMW_Section;

    // ***************************Different def of correction depending on year**********************************
    if (m_Data_Year == 2024){
      FF_DriftTime =  10* (fIC_TS.at(0) - m_TimeData->GetTS_MWPC13(0)) - ((m_TimeData->GetTime_MWPC13(0)+m_TimeData->GetToff_DT13(FPMW_Section))) ;
    }
    else if (m_Data_Year == 2023){
      FF_DriftTime = m_Y;
    }
    else {
      return ;
    }

    //**************************  Correction of section 1 to 4 ***************************************************
    for (int seg = 2; seg < fIC_TS.size() ; seg++) { // loop on multiplicity of event

      if (m_Yspline.at(seg-2)==0){
        ICcorr_Y.at(seg) = fIC_PID[seg]; 
      }

      else {
        ICcorr_Y.at(seg) = fIC_PID[seg] * m_Yspline.at(seg-2)->Eval(0)/ m_Yspline.at(seg-2)->Eval(FF_DriftTime);
        if (!(ICcorr_Y.at(seg)==ICcorr_Y.at(seg))) ICcorr_Y.at(seg) = 0;
      } //end if non empty
      if (seg == 0) break;
    }//endloop seg


    //**************************  Correction of section 0 ***************************************************

    Double_t PolX =  1.37622;

    Double_t ICRatio =  fIC_PID[1]/fIC_PID[0];
    Double_t ICRatioCorr =  ICRatio * PolX / m_IC0_Profile.Evaluate(m_X,FF_DriftTime,false);
    Double_t ICcorr_Y0 = fIC_PID[0] / PolX * m_IC0_Profile.Evaluate(m_X,FF_DriftTime,false);

    if ( ICRatioCorr<1.4 || ICRatioCorr >1.5){
      ICRatioCorr =  ICRatio * PolX / m_IC0_Profile.Evaluate(m_X,FF_DriftTime,false);
      ICcorr_Y0 =     fIC_PID[0] / PolX * m_IC0_Profile.Evaluate(m_X,FF_DriftTime,false);
      if (ICRatioCorr >100) {
        ICcorr_Y0 = fIC_PID[0];
      }
    }


    //**************************  Overwrite ICRAW ***************************************************
    //Overwrite ic_raw
    for (int i = 1 ; i<fIC_TS.size() ;i++){
      if (ICcorr_Y.at(i) != 0 ){
        fIC_PID[i] = ICcorr_Y.at(i);
      }
    }
    fIC_PID[0] = ICcorr_Y0;
  }
}

///////////////////////////////////////////////////////////////////////////
void TICPhysics::Clear() {
  DE = 0;
  Eres = 0;
  Etot = 0;
  EtotInit = 0;
  Chio_Z = 0;
  Chio_ZRaw = 0;
  for(int i=0; i<11; i++){
    fIC[i] = 0;
    fIC_raw[i] = 0;
    fIC_Init[i] = 0;
    fIC_PID[i] = 0;
  }

  fIC_TS.clear();

  for (int ToF=0; ToF<4; ToF ++){
    EtotAoQ[ToF] = 0;
    for(int i=0; i<11; i++){
      fIC_AoQ[ToF][i] = 0;
    }
  }
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("IC");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  IC " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  IC " << i+1 <<  endl;
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

  ReadAnalysisConfig();
}


///////////////////////////////////////////////////////////////////////////
void TICPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  //Calibration from online analysis
  Cal->AddParameter("IC","Z_CALIBRATION","IC_Z_CALIBRATION");
  Cal->AddParameter("IC","INIT_ETOT_SCALING","IC_INIT_ETOT_SCALING");
  for (int segment=0; segment<11;segment++){
    Cal->AddParameter("IC","INIT_SEG"+NPL::itoa(segment+1)+"_ALIGN","IC_INIT_SEG"+NPL::itoa(segment+1)+"_ALIGN");
  }


  for(int section = 0; section<20; section++){
    Cal->AddParameter("IC","ETOT_SCALING_SEC"+NPL::itoa(section),"IC_ETOT_SCALING_SEC"+NPL::itoa(section));
    for(int segment = 0; segment<11; segment++){
      Cal->AddParameter("IC","SEC"+NPL::itoa(section)+"_SEG"+NPL::itoa(segment+1)+"_ALIGN","IC_SEC"+NPL::itoa(section)+"_SEG"+NPL::itoa(segment+1)+"_ALIGN");
    }
  }

  for(string name : ToFName){
    for(int section = 0; section<20; section++){
      Cal->AddParameter("IC","ETOT_SCALING_SEC"+NPL::itoa(section) +"_TOF"+name,"IC_ETOT_SCALING_SEC"+NPL::itoa(section)+"_TOF"+name);
      for(int segment = 0; segment<11; segment++){
        Cal->AddParameter("IC","SEC"+NPL::itoa(section)+"_SEG"+NPL::itoa(segment+1)+"_TOF"+name+"_ALIGN","IC_SEC"+NPL::itoa(section)+"_SEG"+NPL::itoa(segment+1)+"_TOF"+name+"_ALIGN");
      }
    }
  }

}

///////////////////////////////////////////////////////////////////////////
void TICPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("IC",  true );
  inputChain->SetBranchAddress("IC", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("IC", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TICPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("IC", "TICPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TICPhysics::Construct() {
  return (NPL::VDetector*) new TICPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_IC{
    public:
      proxy_IC(){
        NPL::DetectorFactory::getInstance()->AddToken("IC","IC");
        NPL::DetectorFactory::getInstance()->AddDetector("IC",TICPhysics::Construct);
      }
  };

  proxy_IC p_IC;
}

