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
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FPMW Treated  data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TFPMWPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <filesystem>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

//   ROOT
#include "TChain.h"

ClassImp(TFPMWPhysics);

///////////////////////////////////////////////////////////////////////////
TFPMWPhysics::TFPMWPhysics()
  : m_EventData(new TFPMWData),
  m_PreTreatedData(new TFPMWData),
  m_EventPhysics(this),
  m_Zf(7600),
  m_GapSize(0),
  m_NumberOfDetectors(0){
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TFPMWPhysics::AddDetector(TVector3 Pos){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  DetPosX.push_back(Pos.X());
  DetPosY.push_back(Pos.Y());
  DetPosZ.push_back(Pos.Z());

  m_NumberOfDetectors++;
} 

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::AddDetector(double R, double Theta, double Phi){
  // Compute the TVector3 corresponding
  TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta));
  // Call the cartesian method
  AddDetector(Pos);
} 

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::BuildPhysicalEvent() {
  PreTreat();
  unsigned int sizeX = m_PreTreatedData->GetFPMWMultX();
  for(unsigned int i=0; i<sizeX; i++){
    DetectorHitX.insert(m_PreTreatedData->GetFPMW_DetX(i));
  }

  unsigned int sizeY = m_PreTreatedData->GetFPMWMultY();
  for(unsigned int i=0; i<sizeY; i++){
    if(DetectorHitX.find(m_PreTreatedData->GetFPMW_DetY(i)) != DetectorHitX.end()){
      DetectorHit.insert(m_PreTreatedData->GetFPMW_DetY(i));
    }
  }
  unsigned int sizeDet = DetectorHit.size();

  for(unsigned int i=0; i<sizeX; i++){
    int StrX = m_PreTreatedData->GetFPMW_StripX(i);
    int NX = m_PreTreatedData->GetFPMW_DetX(i);
    double QX = m_PreTreatedData->GetFPMW_ChargeX(i);
    if(DetectorHit.find(NX) != DetectorHit.end()){
      MapX[NX].push_back(std::make_pair(StrX,QX));
      QSumX[NX] += QX;
      ChargeX.push_back(QX);
      DetX.push_back(NX);
      if(MaxQX.find(NX)==MaxQX.end() || MaxQX[NX].second<QX){
        MaxQX[NX] = make_pair(StrX,QX);
        //QXmax[NX] = QX;
        //StripXmax[NX] = StrX;
      }
    }
  }
  for(unsigned int i=0; i<sizeY; i++){
    int StrY = m_PreTreatedData->GetFPMW_StripY(i);
    int NY = m_PreTreatedData->GetFPMW_DetY(i);
    double QY = m_PreTreatedData->GetFPMW_ChargeY(i);
    if(DetectorHit.find(NY) != DetectorHit.end()){
      MapY[NY].push_back(std::make_pair(StrY,QY));
      QSumY[NY] += QY;
      ChargeY.push_back(QY);
      DetY.push_back(NY);
      if(MaxQY.find(NY)==MaxQY.end() || MaxQY[NY].second<QY){
        MaxQY[NY] = make_pair(StrY,QY);
        //QYmax[NY] = QY;
        //StripYmax[NY] = StrY;
      }
    }
  }

  for(auto &DetN : DetectorHit){
    //double PosX = WeightedAverage(MapX[DetN]);
    //double PosY = WeightedAverage(MapY[DetN]);
    double PosX = AnalyticHyperbolicSecant(MaxQX[DetN],MapX[DetN]);
    double PosY = AnalyticHyperbolicSecant(MaxQY[DetN],MapY[DetN]);
    //double PosX = FittedHyperbolicSecant(MaxQX[DetN],MapX[DetN]);
    //double PosY = FittedHyperbolicSecant(MaxQY[DetN],MapY[DetN]);

    if(DetN==1 || DetN==2){
      PosX = PosX - DetPosX[DetN-1];
      PosY = PosY - DetPosY[DetN-1];
    }
    else if(DetN==3 || DetN==4){
      PosX = -PosX - DetPosX[DetN-1];
      PosY = PosY - DetPosY[DetN-1];
    }

    int sx0 = (int) PosX;
    int sx1 = sx0+1;
    int sy0 = (int) PosY;
    int sy1 = sy0+1;

    DetectorNbr.push_back(DetN);
    //ChargeX.push_back(QSumX[DetN]);
    //ChargeY.push_back(QSumY[DetN]);
    PositionX.push_back(PosX);
    PositionY.push_back(PosY);

    //Here we linearise the position if the yaml file is found
    if (LoadLinea == true){
      int IndexX = (DetN-1) * 2;
      double LinPosX = -1000, LinPosY = -1000;
      if (PosX > -1000){
        for (int binX = 0; binX < EdgeBinDet[IndexX].size() - 1; binX++) {
          // If X is in the correct range
          if (EdgeBinDet[IndexX][binX].size() != 0) {
            if (PosX >= EdgeBinDet[IndexX][binX][0] &&
                PosX < EdgeBinDet[IndexX][binX]
                [EdgeBinDet[IndexX][binX].size() - 1]) {
              // Apply the coorect transformation to X
              LinPosX = LinearisePos(PosX, EdgeBinDet[IndexX][binX],
                  NewEdgeBinDet[IndexX][binX]);
              LinPositionX.push_back(LinPosX);

            } // end check in range
          }   // end check empty
        }     // end loop mmX
      } // END CHECK
      else {
        LinPositionX.push_back(-1000);
      } 
      int IndexY = (DetN-1) * 2 + 1;
      if (PosY > -1000){
        for (int binY = 0; binY < EdgeBinDet[IndexY].size() - 1;
            binY++) {
          if (EdgeBinDet[IndexY][binY].size() != 0) {
            if (PosY >= EdgeBinDet[IndexY][binY][0] &&
                PosY < EdgeBinDet[IndexY][binY]
                [EdgeBinDet[IndexY][binY].size() - 1]) {
              LinPosY = LinearisePos(PosY, EdgeBinDet[IndexY][binY],
                  NewEdgeBinDet[IndexY][binY]);
              LinPositionY.push_back(LinPosY);
            } // end check in range
          }   // End non empty
        }     // End loop BINY
      }       // If real event END
      else {
        LinPositionY.push_back(-1000);
      }
    }
    else {
      LinPositionX.push_back(-1000);  
      LinPositionY.push_back(-1000);
    }
  } //End Linea

  CalculateFocalPlanePosition(m_Zf);
  CalculateTargetPosition();
}
/////////////////////////////////////////////////////////////////////
double TFPMWPhysics::LinearisePos(double X, vector<double> &EdgeBin,
    vector<double> &NewEdgeBin) {
  // ************************************************************************
  // Transform an inputed reconstructed position in order to linearise the
  // distribution of event in a FPMW.
  //
  // Input :
  //    EdgeBin : The range of the bin of the FPMW to treat, this vector is
  //    generated by LinearisationYaml.cxx
  //
  //    NewEdgeBin : The corresponding vector of new bin associated to each old
  //    bin of the FPMW. Those bins have a size made to linearise the
  //    distribution. They are also generated by LinearisationYaml.cxx
  //
  //    X : The reconstructed position before linearisation
  // Output :
  //    Xnew : The reconstructed position after linearisation.
  // ************************************************************************
  double XNew = -1000;
  bool TransfoDone = false;
  for (int i = 0; i < EdgeBin.size() - 1; i++) {
    if (X >= EdgeBin[i] && X < EdgeBin[i + 1]) {

      if (TransfoDone == true ) {
        continue;
      }

      double Xold = X;
      double Scale2 = (X - EdgeBin[i]) / (EdgeBin[i + 1] - EdgeBin[i]);
      XNew = Scale2 * (NewEdgeBin[i + 1] - NewEdgeBin[i]) + NewEdgeBin[i];
      TransfoDone = true;
    }
  }
  return XNew;
} // end transfo x


/////////////////////////////////////////////////////////////////////////

void TFPMWPhysics::LoadLineaFromFile(vector<vector<vector<double>>> &vec1, vector<vector<vector<double>>> &vec2, const string &filename) {
  ifstream inFile(filename);
  if (!inFile.is_open()) {
    cerr << "Unable to open file for reading: " << filename << endl;
    return;
  }

  vec1.clear();
  vec2.clear();
  vector<vector<vector<double>>>* currentVec = &vec1;
  string line;
  vector<vector<double>> det;
  vector<double> bin;

  while (getline(inFile, line)) {
    if (line.find("Vector 2:") != string::npos) {
      if (!bin.empty()) {
        det.push_back(bin);
        bin.clear();
      }
      if (!det.empty()) {
        vec1.push_back(det);
        det.clear();
      }
      currentVec = &vec2;
      continue;
    }

    if (line.find("Det") != string::npos) {
      if (!bin.empty()) {
        det.push_back(bin);
        bin.clear();
      }
      if (!det.empty()) {
        currentVec->push_back(det);
        det.clear();
      }
      continue;
    }

    if (line.find("Bin") != string::npos) {
      if (!bin.empty()) {
        det.push_back(bin);
        bin.clear();
      }
      istringstream iss(line.substr(line.find(":") + 1));
      double value;
      while (iss >> value) {
        bin.push_back(value);
      }
      det.push_back(bin);
      bin.clear();
      continue;
    }

    if (line.empty()) {
      if (!bin.empty()) {
        det.push_back(bin);
        bin.clear();
      }
      if (!det.empty()) {
        currentVec->push_back(det);
        det.clear();
      }
    }
  }

  if (!bin.empty()) {
    det.push_back(bin);
  }
  if (!det.empty()) {
    currentVec->push_back(det);
  }

  inFile.close();
}

/////////////////////////////////////////////////////////////////////////
// Function to calculate memory usage of a 3D vector
double TFPMWPhysics::calculate_memory_usage(vector<vector<vector<double>>>& vec) {
  double total_memory = 0;

  // Memory for the outer vector
  total_memory += vec.capacity() * sizeof(vector<vector<double>>);

  for (auto& mid_vec : vec) {
    // Memory for the middle vectors
    total_memory += mid_vec.capacity() * sizeof(vector<double>);

    for (auto& inner_vec : mid_vec) {
      // Memory for the inner vectors
      total_memory += inner_vec.capacity() * sizeof(double);
    }
  }

  return total_memory;
}
/////////////////////////////////////////////////////////////////////
void TFPMWPhysics::CalculateFocalPlanePosition(double Zf){

  if(PositionX.size()==4){
    double Z2 = DetPosZ[2]-m_GapSize;
    double Z3 = DetPosZ[3]-m_GapSize;
    double X2 = PositionX[2];
    double X3 = PositionX[3];
    double Y2 = PositionY[2];
    double Y3 = PositionY[3];

    double t = (Z3-Zf)/(Z2-Zf);

    Z2 = DetPosZ[2]+m_GapSize;
    Z3 = DetPosZ[3]+m_GapSize;

    if(Z2!=-1000 && Z3 != -1000){
      if ( (abs(X2) < 1000) && (abs(X3) < 1000)){
        Thetaf = atan((X3-X2)/(Z3-Z2));
        Xf = 1./(t-1)*(t*X2-X3);
      }
      if ((abs(Y2) < 1000) && (abs(Y3) < 1000)){
        Phif = atan((Y3-Y2)/(Z3-Z2));
        Yf = 1./(t-1)*(t*Y2-Y3);
      }
    }
    if (LoadLinea == true){
      Z2 = DetPosZ[2]-m_GapSize;
      Z3 = DetPosZ[3]-m_GapSize; 
      double X2Lin = LinPositionX[2];
      double X3Lin = LinPositionX[3];
      double Y2Lin = LinPositionY[2];
      double Y3Lin = LinPositionY[3];

      double t = (Z3-Zf)/(Z2-Zf);

      Z2 = DetPosZ[2]+m_GapSize;
      Z3 = DetPosZ[3]+m_GapSize;


      if(Z2!=-1000 && Z3 != -1000){
        if ((abs(X2) < 1000) && (abs(X3) < 1000)){
          ThetafLin = atan((X3Lin-X2Lin)/(Z3-Z2));
          XfLin = 1./(t-1)*(t*X2Lin-X3Lin);
        }
        if ((abs(Y2) < 1000) && (abs(Y3) < 1000)){
          PhifLin = atan((Y3Lin-Y2Lin)/(Z3-Z2));
          YfLin = 1./(t-1)*(t*Y2Lin-Y3Lin);
        }
      }  
    } //end load yaml true
  }
}

/////////////////////////////////////////////////////////////////////
void TFPMWPhysics::CalculateTargetPosition(){

  if(PositionX.size()>1 && PositionY.size()>1){
    double Z0 = DetPosZ[0]-m_GapSize;
    double Z1 = DetPosZ[1]-m_GapSize;
    double X0 = PositionX[0];
    double X1 = PositionX[1];
    double Y0 = PositionY[0];
    double Y1 = PositionY[1];
    Xt = (1./(Z1-Z0))*(X0*Z1-X1*Z0);

    Z0 = DetPosZ[0]+m_GapSize;
    Z1 = DetPosZ[1]+m_GapSize;
    Yt = (1./(Z1-Z0))*(Y0*Z1-Y1*Z0);
    Z0 = DetPosZ[0];
    Z1 = DetPosZ[1];
    //TVector3 vFF = TVector3(X1-X0,Y1-Y0,Z1-Z0);
    Theta_in = atan((X1-X0)/(Z1-Z0));
    Phi_in = atan((Y1-Y0)/(Z1-Z0));

    if (LoadLinea == true){
      double Z0Lin = DetPosZ[0]-m_GapSize;
      double Z1Lin = DetPosZ[1]-m_GapSize;
      double X0Lin = LinPositionX[0];
      double X1Lin = LinPositionX[1];
      double Y0Lin = LinPositionY[0];
      double Y1Lin = LinPositionY[1];
      XtLin = (1./(Z1Lin-Z0Lin))*(X0Lin*Z1Lin-X1Lin*Z0Lin);

      Z0Lin = DetPosZ[0]+m_GapSize;
      Z1Lin = DetPosZ[1]+m_GapSize;
      YtLin = (1./(Z1Lin-Z0Lin))*(Y0Lin*Z1Lin-Y1Lin*Z0Lin);

      Z0Lin = DetPosZ[0];
      Z1Lin = DetPosZ[1];
      //TVector3 vFF = TVector3(X1-X0,Y1-Y0,Z1-Z0);
      Theta_in_Lin = atan((X1Lin-X0Lin)/(Z1Lin-Z0Lin));
      Phi_in_Lin = atan((Y1Lin-Y0Lin)/(Z1Lin-Z0Lin));
    } //End Yaml
  }

}

/////////////////////////////////////////////////////////////////////
double TFPMWPhysics::WeightedAverage(std::vector<std::pair<int,double>>& Map){

  unsigned int sizeQ = Map.size(); 

  double num=0;
  double denum=0;
  for(unsigned int i = 0 ; i < sizeQ ; i++){
    num += Map[i].first * Map[i].second;
    denum += Map[i].second;
  }

  double weighted_average = num/denum;

  return weighted_average;


}

/////////////////////////////////////////////////////////////////////
double TFPMWPhysics::AnalyticHyperbolicSecant(std::pair<int,double>& MaxQ,std::vector<std::pair<int,double>>& Map){
  double sech = -1000 ;

  // std::cout << "test AnH 1" << std::endl;
  if(MaxQ.second > 0 && MaxQ.first > 2 && MaxQ.first<996){	
    //      if(Buffer_Q[MaxQ.first-1+1]==0||Buffer_Q[MaxQ.first-1-1]==0)
    //        return sech;
    // std::cout << "test AnH 2" << std::endl;
    double q2 = MaxQ.second;
    double q1 = 0,q3 = 0;
    for(auto &strip : Map){
      if(strip.first == MaxQ.first - 1){
        q1 = strip.second;
      }
      else if(strip.first == MaxQ.first + 1){
        q3 = strip.second;
      }
    }
    //std::cout << "test q " << q1 << " " << q2 << " " << q3 << std::endl;
    double vs[6];
    // std::cout << "test AnH 3" << std::endl;
    if(q1 > 0 && q3 > 0)    
    {
      // QsumSample[DetNum].push_back(QSum);
      vs[0] = sqrt(q2/q3);
      vs[1] = sqrt(q2/q1);
      vs[2] = 0.5*(vs[0] + vs[1]);
      vs[3] = log( vs[2] + sqrt(vs[2]*vs[2]-1.0) );
      vs[4] = abs((vs[0] - vs[1])/(2.0*sinh(vs[3])));	
      vs[5] = 0.5*log( (1.0+vs[4])/(1.0-vs[4]) ) ;
      // std::cout << "test AnH 4" << std::endl;

      if ( q3>q1 ) 
        sech = MaxQ.first + vs[5]/vs[3] ;
      else 
        sech = MaxQ.first - vs[5]/vs[3] ;
      //std::cout << "test sech " << sech << std::endl;
    }
  }

  return sech ;
}

///////////////////////////////////////////////////////////////////////////
double TFPMWPhysics::FittedHyperbolicSecant(std::pair<int,double>& MaxQ,std::vector<std::pair<int,double>>& Map){
  // Warning: should not delete static variable
  static TF1* f = new TF1("sechs","[0]/(cosh(TMath::Pi()*(x-[1])/[2])*cosh(TMath::Pi()*(x-[1])/[2]))",1,1000);

  // Help the fit by computing the position of the maximum by analytic method
  double StartingPoint = AnalyticHyperbolicSecant(MaxQ,Map);
  // if analytic method fails then the starting point in strip max
  if(StartingPoint==-1000) StartingPoint = MaxQ.first; 

  // Maximum is close to charge max, Mean value is close to Analytic one, typical width is 3.8 strip
  f->SetParameters(MaxQ.first,StartingPoint,1);

  static vector<double> y ;
  static vector<double> q ; 
  y.clear(); q.clear();
  double final_size = 0 ;
  unsigned int sizeQ = Map.size(); 

  for(unsigned int i = 0 ; i < sizeQ ; i++){
    if(Map[i].second > (MaxQ.second)*0.1){
      q.push_back(Map[i].second);
      y.push_back(Map[i].first);
      final_size++;
    }
  }

  // requiered at least 3 point to perfom a fit
  if(final_size<3){
    return -1000 ;
  }

  TGraph* g = new TGraph(q.size(),&y[0],&q[0]);
  g->Fit(f,"QN0");
  delete g;
  return f->GetParameter(1)  ;
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Might test for disabled channels for more complex detector

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  unsigned int mysizeX = m_EventData->GetFPMWMultX();
  for (unsigned int i = 0; i < mysizeX ; ++i) {
    int det = m_EventData->GetFPMW_DetX(i);
    int strip = m_EventData->GetFPMW_StripX(i);
    double QX = m_EventData->GetFPMW_ChargeX(i);
    double a0 = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPX"+NPL::itoa(strip),0);
    double a1 = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPX"+NPL::itoa(strip),1);
    double a2 = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPX"+NPL::itoa(strip),2);
    double scale = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPX"+NPL::itoa(strip),3);
    QX = QX - a0 + rand.Uniform();
    double Qcal = QX*a1 + pow(QX,2)*a2;
    Qcal = scale*Qcal;

    /*if(det==1){
      cout << "*******************" << endl;
      cout << "strip= " << strip << endl;
      cout << a0 << " " << a1 << " " << a2 << " " << scale << endl;
      cout << "QX= " << QX << endl;
      cout << "Qcal= " << Qcal << endl;
      }*/
    if(Qcal>100){
      m_PreTreatedData->SetFPMW_DetX(det);
      m_PreTreatedData->SetFPMW_StripX(strip);
      m_PreTreatedData->SetFPMW_ChargeX(Qcal);
    }
  }

  unsigned int mysizeY = m_EventData->GetFPMWMultY();
  for (unsigned int i = 0; i < mysizeY ; ++i) {
    int det = m_EventData->GetFPMW_DetY(i);
    int strip = m_EventData->GetFPMW_StripY(i);
    double QY = m_EventData->GetFPMW_ChargeY(i);

    double a0 = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPY"+NPL::itoa(strip),0);
    double a1 = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPY"+NPL::itoa(strip),1);
    double a2 = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPY"+NPL::itoa(strip),2);
    double scale = Cal->GetValue("FPMW/DET"+NPL::itoa(det)+"_STRIPY"+NPL::itoa(strip),3);
    QY = QY - a0 + rand.Uniform();
    double Qcal = QY*a1 + pow(QY,2)*a2;
    Qcal = scale*Qcal;

    if(Qcal>100){
      m_PreTreatedData->SetFPMW_DetY(det);
      m_PreTreatedData->SetFPMW_StripY(strip);
      m_PreTreatedData->SetFPMW_ChargeY(Qcal);
    }
  }

}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigFPMW.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigFPMW.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigFPMW.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigFPMW.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigFPMW";
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

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }

      else if (whatToDo == "LINEA_PATH") {
        AnalysisConfigFile >> DataBuffer;

        // ******************** Linearisation init******************
        string FileNameLinea = DataBuffer;
        ifstream LineaFile (FileName.c_str());	  
        if (LineaFile){
          cout << endl <<"/////////// txt Linearisation File opened ///////////" << endl;
          LoadLineaFromFile(EdgeBinDet, NewEdgeBinDet, FileNameLinea);
          LoadLinea = true ;
        }
        else if (!LineaFile){
          perror ( "Stream didnt open");	
        }
        //**********************************************************
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}


///////////////////////////////////////////////////////////////////////////
TVector3 TFPMWPhysics::GetHitPosition(int Det){
  TVector3 Hit = TVector3(PositionX[Det], PositionY[Det], GetDetectorPositionZ(Det));

  return Hit;
}

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::Clear() {
  DetectorNbr.clear();
  DetX.clear();
  DetY.clear();
  PositionX.clear();
  PositionY.clear();
  ChargeX.clear();
  ChargeY.clear();

  LinPositionX.clear();
  LinPositionY.clear();

  DetectorHitX.clear();
  DetectorHit.clear();

  MapX.clear();
  MapY.clear();
  MaxQX.clear();
  MaxQY.clear();

  Xf = -1000;
  Yf = -1000;
  XfLin = -1000;
  YfLin = -1000; 
  Thetaf = -1000;
  Phif = -1000;
  Theta_in = -1000;
  Phi_in = -1000;
  ThetafLin = -1000;
  PhifLin = -1000;
  Theta_in_Lin = -1000;
  Phi_in_Lin = -1000;

  Xt = -1000;
  Yt = -1000;
  XtLin = -1000;
  YtLin = -1000;
  /*for(int i=0; i<4; i++){
    QXmax[i] = -1;
  //QYmax[i] = -1;
  StripXmax[i] = -1;
  //StripYmax[i] = -1;
  }*/
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("FPMW");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"POS"};
  vector<string> sphe = {"R","Theta","Phi"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(cart)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FPMW " << i+1 <<  endl;

      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      AddDetector(Pos);
    }
    else if(blocks[i]->HasTokenList(sphe)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  FPMW " << i+1 <<  endl;
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
void TFPMWPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for(int i = 0; i < m_NumberOfDetectors; i++){
    for(int s = 0; s < 996; s++){
      Cal->AddParameter("FPMW","DET"+NPL::itoa(i+1)+"_STRIPX"+NPL::itoa(s+1),"FPMW_DET"+NPL::itoa(i+1)+"_STRIPX"+NPL::itoa(s+1));
      Cal->AddParameter("FPMW","DET"+NPL::itoa(i+1)+"_STRIPY"+NPL::itoa(s+1),"FPMW_DET"+NPL::itoa(i+1)+"_STRIPY"+NPL::itoa(s+1));
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("FPMW",  true );
  inputChain->SetBranchAddress("FPMW", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("FPMW", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TFPMWPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("FPMW", "TFPMWPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TFPMWPhysics::Construct() {
  return (NPL::VDetector*) new TFPMWPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_FPMW{
    public:
      proxy_FPMW(){
        NPL::DetectorFactory::getInstance()->AddToken("FPMW","FPMW");
        NPL::DetectorFactory::getInstance()->AddDetector("FPMW",TFPMWPhysics::Construct);
      }
  };

  proxy_FPMW p_FPMW;
}

