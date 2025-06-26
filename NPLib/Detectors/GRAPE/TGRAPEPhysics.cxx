/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           * 
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold GRAPE treated data                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

// STL
#include <cmath>
#include <stdlib.h>
#include <limits>
#include <fstream>
using namespace std;

#include "TGRAPEPhysics.h"

//   NPL
#include "RootInput.h"
#include "NPDetectorFactory.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "NPSystemOfUnits.h"
#include "NPOptionManager.h"
using namespace NPUNITS;

//   ROOT
#include "TChain.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TGRAPEPhysics)
///////////////////////////////////////////////////////////////////////////
TGRAPEPhysics::TGRAPEPhysics()  {
  m_EventData         = new TGRAPEData ;
  m_PreTreatedData    = new TGRAPEData ;
  m_EventPhysics      = this ;
  m_NumberOfDetectors = 0;
  m_NumberOfSegments  = 18;}

/////////////////////////////////////////////////
void TGRAPEPhysics::BuildPhysicalEvent(){
  PreTreat();

  // Addback Map
  unsigned int mysize = Gamma_Energy.size();
  for (unsigned int g = 0 ; g < mysize ; g++){
    int GRAPENum = GRAPE_Number[g];
    m_map_E[GRAPENum-1] += Gamma_Energy[g];

    if( Gamma_Energy[g] > m_map_Segment_MaxE[GRAPENum-1]) {
      m_map_Segment_MaxE[GRAPENum-1] = Gamma_Energy[g];
      m_map_Segment[GRAPENum-1] = Segment_Number[g]; }}
  
  // Final Addback and Doppler Correction 
  int zero = 0;
  for (int i = 0 ; i < 18 ; i++) {
    if (m_map_E.find(i) != m_map_E.end()) {
      int grape = i+1;
      TVector3 Pos;
      Pos = GetSegmentPosition(grape, 0, m_map_Segment[i]);

      if (Pos.Mag() != 0) {
        static TVector3 Beta = TVector3(0,0,0.125);
        double E = GetDopplerCorrectedEnergy(m_map_E[i],Pos,Beta);
        AddBack_DC.push_back(E);
        AddBack_E.push_back(m_map_E[i]);
        AddBack_Theta.push_back(Pos.Angle(Beta)*180./3.141592653589793);
        AddBack_GRAPE.push_back(grape); 
	AddBack_Segment.push_back(m_map_Segment[i]);
        AddBack_X.push_back(Pos.X());
        AddBack_Y.push_back(Pos.Y());
        AddBack_Z.push_back(Pos.Z());}}}
  return; }

/////////////////////////////////////////////////
void TGRAPEPhysics::PreTreat(){
  static CalibrationManager* cal = CalibrationManager::getInstance();
  static string name;
  unsigned int mysize = m_EventData->GetMultiplicityGe();
  double Eraw, Energy;
  int grape, crystal, segment;
  for (unsigned int i = 0 ; i < mysize ; i++) {
    Eraw = m_EventData->GetGeEnergy(i);
    if ( Eraw > 0 ) {
      grape   = m_EventData->GetGeGRAPENbr  (i);
      crystal = m_EventData->GetGeCrystalNbr(i);
      segment = m_EventData->GetGeSegmentNbr(i);
      name = "GRAPE/D"+ NPL::itoa(grape)+"_SEG"+ NPL::itoa(segment)+"_E";
      Energy =  cal->ApplyCalibration(name, Eraw);
      Gamma_Energy.  push_back(Energy);
      GRAPE_Number.  push_back(grape);
      Crystal_Number.push_back(crystal);
      Segment_Number.push_back(segment);
      Gamma_Time.    push_back(m_EventData->GetGeTimeCFD(i));}}
}

/////////////////////////////////////////////////
TVector3 TGRAPEPhysics::GetPositionOfInteraction(unsigned int& i){
  return GetSegmentPosition(GRAPE_Number[i],Crystal_Number[i],Segment_Number[i]);
}
/////////////////////////////////////////////////
// original energy, position, beta
double TGRAPEPhysics::GetDopplerCorrectedEnergy(double& energy , TVector3 position, TVector3& beta){
  // renorm pos vector
  position.SetMag(1); 
  m_GammaLV.SetPx(energy*position.X());
  m_GammaLV.SetPy(energy*position.Y());
  m_GammaLV.SetPz(energy*position.Z());
  m_GammaLV.SetE(energy);
  m_GammaLV.Boost(-beta);
  return m_GammaLV.Energy();}

/////////////////////////////////////////////////
// Add clover at the standard position of the array
// Take as argument the standard clover Id.
void TGRAPEPhysics::AddGRAPEStandard(vector<int> GRAPEId){

    //////////////////////////////////////////////////
  // standard location
  double R_std[18] = {148, 138, 138, 138, 138, 148,
		      126, 101, 101, 101, 101, 126,
		      148, 138, 138, 138, 138, 148 };
  double Theta_std[18] = { 130, 125, 125, 125, 125, 130,
			    90,  90,  90,  90,  90,  90,
			    55,  55,  55,  55,  55,  55 };
  double Phi_std[18] = { 30, 90, 150, 210, 270, 330,
			 0, 60, 120, 180, 240, 300,
			 30, 90, 150, 210, 270, 330 };
  double betaX_std[18] = { 17, 20, 20, 20, 20, 17,
			   0, 0, 0, 0, 0, 0,
			   -17, -20, -20, -20, -20, -17 };
  double betaY_std[18] = { 0. };
  double betaZ_std[18] = { 0. };
  //////////////////////////////////////////////////
  
  for (unsigned int i = 0 ;  i < GRAPEId.size(); i++) {
    if (GRAPEId[i] < 1 || GRAPEId[i] > 18) continue;
    TVector3 Pos(0,0,1);
    int ID = GRAPEId[i];
    Pos.SetMag(R_std[GRAPEId[i] - 1]);
    Pos.SetTheta(Theta_std[GRAPEId[i] - 1]*deg);
    Pos.SetPhi(Phi_std[GRAPEId[i] - 1]*deg);
    m_GRAPEPosition[ID] = Pos;
    m_NumberOfDetectors++;}
  return;}

/////////////////////////////////////////////////
void TGRAPEPhysics::AddGRAPE(unsigned int ID ,double R, double Theta, double Phi){
  TVector3 Pos(0,0,1);
  Pos.SetTheta(Theta);
  Pos.SetPhi(Phi);
  Pos.SetMag(R);
  m_GRAPEPosition[ID] = Pos;
  m_NumberOfDetectors++;
  return;}
/////////////////////////////////////////////////
TVector3 TGRAPEPhysics::GetGRAPEPosition(int GRAPENbr){
  return m_GRAPEPosition[GRAPENbr];}

/////////////////////////////////////////////////
TVector3 TGRAPEPhysics::GetSegmentPosition(int GRAPENbr, int CrystalNbr, int SegmentNbr){
    TVector3 GRAPEPos   = GetGRAPEPosition(GRAPENbr);

    TVector3 GRAPEPosW = GRAPEPos.Unit();
    TVector3 GRAPEPosV(cos(GRAPEPosW.Theta()) * cos(GRAPEPosW.Phi()),
		       cos(GRAPEPosW.Theta()) * sin(GRAPEPosW.Phi()),
		       -sin(GRAPEPosW.Theta()));
    TVector3 GRAPEPosU = GRAPEPosV.Cross(GRAPEPosW);
    GRAPEPosU = GRAPEPosU.Unit();
    GRAPEPosV = GRAPEPosW.Cross(GRAPEPosU);
    GRAPEPosV = GRAPEPosV.Unit();
    
    TRotation GRAPERot;
    GRAPERot.RotateAxes(GRAPEPosU, GRAPEPosV, GRAPEPosW);

    
    double rmaxSD = 32.;
    double offs = 20. + (rmaxSD-30.)/2.;
    TVector3 Pos;
    int i = SegmentNbr - 1;
    Pos.SetXYZ(-((i%3)-1)*offs,(1 - abs(int((i - 8.5)/3)))*offs,((i/9)*2-1)*12.50);
  
    Pos.Transform(GRAPERot);
    Pos += GRAPEPos;
    return Pos;}

/////////////////////////////////////////////////
void TGRAPEPhysics::ReadConfiguration(NPL::InputParser parser)  {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GRAPE");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " free clovers found " << endl; 

  vector<string> token = {"GRAPEID","R","Theta","Phi","Beta"}; 
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      vector<double> beta = blocks[i]->GetVectorDouble("Beta","deg"); 
      int     id = blocks[i]->GetInt("GRAPEID");
      AddGRAPE(id,R,Theta,Phi);
    }


    else{
      cout << "Warning: check your input file formatting " << endl;
    }
  }
  /*  std::ofstream fout("anapos.txt");
  //////////////////////////////////////////////////
  // just for check
  for (unsigned int j = 1 ; j <= 18 ; j++) {
    fout << " Segment Position for Det" << j << std::endl;
    for (unsigned int i = 1 ; i <= 18 ; i++) {
      TVector3 segPos = GetSegmentPosition(j, 0, i);
      fout << i << ": (" << segPos.X() << ", "<< segPos.Y() << ", "<< segPos.Z() << ")" << std::endl; }}
      fout.close();
  //////////////////////////////////////////////////*/
  blocks.clear();}

///////////////////////////////////////////////////////////////////////////
void TGRAPEPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus( "GRAPE" , true );
  if(inputChain->FindBranch( "fGRAPE_*" ))
    inputChain->SetBranchStatus( "fGRAPE_*" , true );
  inputChain->SetBranchAddress( "GRAPE" , &m_EventData );
}

///////////////////////////////////////////////////////////////////////////
void TGRAPEPhysics::InitializeRootOutput()    {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch( "GRAPE" , "TGRAPEPhysics" , &m_EventPhysics );
}
///////////////////////////////////////////////////////////////////////////  
void TGRAPEPhysics::Clear() {
  
  Gamma_Energy.clear();
  Gamma_Time.clear();
  Crystal_Number.clear();
  GRAPE_Number.clear();
  Segment_Number.clear();

  AddBack_E.clear();
  AddBack_DC.clear();
  AddBack_Theta.clear();
  AddBack_GRAPE.clear();
  AddBack_Segment.clear();
  AddBack_X.clear();
  AddBack_Y.clear();
  AddBack_Z.clear();

  m_map_E.clear();
  m_map_Segment.clear(); 
  m_map_Segment_MaxE.clear(); 

}
///////////////////////////////////////////////////////////////////////////  
void TGRAPEPhysics::ClearEventData() {
  m_EventData->Clear();
  m_PreTreatedData->Clear();
}
///////////////////////////////////////////////////////////////////////////
void TGRAPEPhysics::AddParameterToCalibrationManager(){
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for(int i = 0 ; i < 16; ++i){
    for( int j = 0 ; j < 10 ; ++j){
      Cal->AddParameter("GRAPE", "D"+ NPL::itoa(i+1)+"_SEG"+ NPL::itoa(j)+"_E",
			"GRAPE_D"+ NPL::itoa(i+1)+"_SEG"+NPL::itoa(j)+"_E");}}
  return;}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TGRAPEPhysics::Construct(){
  return (NPL::VDetector*) new TGRAPEPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_grape{
  public:
    proxy_grape(){
      NPL::DetectorFactory::getInstance()->AddToken("GRAPE","GRAPE");
      NPL::DetectorFactory::getInstance()->AddDetector("GRAPE",TGRAPEPhysics::Construct);
    }
};

proxy_grape p;
}

