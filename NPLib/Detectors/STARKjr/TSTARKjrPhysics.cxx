
/*****************************************************************************
 * Original Author: Jongwon Hwang    contact address: jwhwang@ibs.re.kr      * 
 * Creation Date  : May 2023                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold STARKjr Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include "RootInput.h"
#include "RootOutput.h"
#include "TChain.h"
#include "TAsciiFile.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPInputParser.h"
#include "NPSystemOfUnits.h"
#include "TSTARKjrPhysics.h"
#include "TRotation.h"

using namespace NPUNITS;
using namespace std;

ClassImp(TSTARKjrPhysics)

///////////////////////////////////////////////////////////////////////////
TSTARKjrPhysics::TSTARKjrPhysics(){
  m_EventData      = new TSTARKjrData;
  m_PreTreatedData = new TSTARKjrData;
  m_EventPhysics   = this;

  ////////////////////////////////////////////////////////////
  // for geometry
  X6_activeX = 40.30;  X6_activeY = 75.00; // mm
  X6_NFrontStrips = 8; X6_NBackStrips = 4;

  BB10_activeX = 39.45;  BB10_activeY = 74.15; // mm
  BB10_NFrontStrips = 8; BB10_NBackStrips = 1;

  QQQ5_activeInR = 25.25;  QQQ5_activeOutR = 81.95; // mm
  QQQ5_NRStrips = 32; QQQ5_NAStrips = 4;
  ////////////////////////////////////////////////////////////  
}

///////////////////////////////////////////////////////////////////////////
TSTARKjrPhysics::~TSTARKjrPhysics() {}

void TSTARKjrPhysics::AddDetector(string Type, TVector3 Pos, int Flip, double Beta) {
  m_Type.push_back(Type);
  m_DetPos.push_back(Pos);
  m_Flip.push_back(Flip);
  m_Beta.push_back(Beta);
  
  ////////////////////////////////////////////////////////////
  // Detector orientation (only useful for X6, but......)
  TVector3 V(0,1,0); // Y
  if (Type == "QQQ5") {
    V.RotateZ(Beta);
    if (Flip) V.RotateY(TMath::DegToRad()*180);}
  else { // for X6, BB10
    if (Flip) V.RotateZ(TMath::DegToRad()*180);
    V.RotateX(TMath::DegToRad()*90);
    V.RotateZ(TMath::DegToRad()*90 + Pos.Phi());}
  m_DetOri.push_back(V);
  ////////////////////////////////////////////////////////////

  AddStripPosition(Type, Pos, Flip, Beta);}

void TSTARKjrPhysics::AddStripPosition(string Type, TVector3 Pos, int Flip, double Beta) {
  //////////////////////////////////////////////////////
  // local coordinate, unit vector (rotation)
  TRotation rot;
  if (Type == "QQQ5") {
    rot.RotateZ(Beta);
    if (Flip) rot.RotateY(TMath::Pi()); }
  else { // for X6, BB10
    if (Flip) rot.RotateZ(TMath::Pi());
    rot.RotateX(TMath::PiOver2());
    rot.RotateZ(TMath::PiOver2() + Pos.Phi());
	rot.RotateZ(Beta);
}
  //////////////////////////////////////////////////////

  if (Type == "X6") {
    ////////////////////////////////////////////////////////////
    // Strip position calculation for X6
    ////////////////////////////////////////////////////////////
    double X6_frontStripPitch = X6_activeX/X6_NFrontStrips;
    double X6_backStripPitch  = X6_activeY/X6_NBackStrips;

    vector<vector<TVector3>> frontStripPos;
    for (int i = 0 ; i < X6_NFrontStrips ; i++) {
      vector<TVector3> backStripPos;
      for (int j = 0 ; j < X6_NBackStrips ; j++) {
	TVector3 PosStrip;
	PosStrip.SetX(-X6_frontStripPitch * (i + 0.5) + X6_activeX/2);
	PosStrip.SetY(X6_backStripPitch * (j + 0.5) - X6_activeY/2);
	PosStrip = rot * PosStrip;
	PosStrip += Pos;
	backStripPos.push_back(PosStrip);}
      frontStripPos.push_back(backStripPos);}
    m_StripPos.push_back(frontStripPos);
    ////////////////////////////////////////////////////////////

  } else if (Type == "BB10") {
    ////////////////////////////////////////////////////////////
    // Strip position calculation for BB10
    ////////////////////////////////////////////////////////////
    double BB10_frontStripPitch = BB10_activeX/BB10_NFrontStrips;
    double BB10_backStripPitch = BB10_activeY/BB10_NBackStrips;

    vector<vector<TVector3>> frontStripPos;
    for (int i = 0 ; i < BB10_NFrontStrips ; i++) {
      vector<TVector3> backStripPos;
      for (int j = 0 ; j < BB10_NBackStrips ; j++) {
	TVector3 PosStrip;
	PosStrip.SetX(-BB10_frontStripPitch * (i + 0.5) + BB10_activeX/2);
	PosStrip.SetY(BB10_backStripPitch * (j + 0.5) - BB10_activeY/2);
	PosStrip = rot * PosStrip;
	PosStrip += Pos;
	backStripPos.push_back(PosStrip);}
      frontStripPos.push_back(backStripPos);}
    m_StripPos.push_back(frontStripPos);
    ////////////////////////////////////////////////////////////
  } else if (Type == "QQQ5") {
    ////////////////////////////////////////////////////////////
    // Strip position calculation for QQQ5
    ////////////////////////////////////////////////////////////
    double QQQ5_rStripPitch = (QQQ5_activeOutR - QQQ5_activeInR)/QQQ5_NRStrips;
    double QQQ5_aStripPitch = TMath::PiOver2()/QQQ5_NAStrips;

    vector<vector<TVector3>> frontStripPos; // Radial Strip
    for (int i = 0 ; i < QQQ5_NRStrips ; i++) {
      vector<TVector3> backStripPos; // Annular strips
      for (int j = 0 ; j < QQQ5_NAStrips ; j++) {
	Double_t r = QQQ5_activeOutR - QQQ5_rStripPitch * (i + 0.5);
	Double_t a = QQQ5_aStripPitch * (i + 0.5);
	TVector3 PosStrip;
	PosStrip.SetMagThetaPhi(r, TMath::PiOver2(), a);
	PosStrip = rot * PosStrip;
	PosStrip += Pos;
	backStripPos.push_back(PosStrip);}
      frontStripPos.push_back(backStripPos);}
    m_StripPos.push_back(frontStripPos);
    ////////////////////////////////////////////////////////////

  }

  ////////////////////////////////////////////////////////////
  // for checking
  auto lastDet = m_StripPos.back();
  for (auto it2 = lastDet.begin() ;
       it2 != lastDet.end() ; it2++) {
    size_t frontIdx = it2 - lastDet.begin();
    for (auto it3 = it2->begin() ;
	 it3 != it2->end() ; it3++) {
      size_t backIdx = it3 - it2->begin();
      std::cout << "(" << frontIdx << ", " << backIdx << "): ";
      it3->Print();
    }}
  ////////////////////////////////////////////////////////////
}
////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::BuildPhysicalEvent() {    
  Clear();

  ////////////////////////////////////////////////////////////
  // Loop for All
  nhit = m_EventData->GetMult();
  for (Int_t i = 0 ; i < m_EventData->GetMult() ; i++) {
    type [i] = m_EventData->GetType(i); // 0: X6, 1: BB10, 2: QQQ5
    detN [i] = m_EventData->GetDetN(i);
    fStrN[i] = m_EventData->GetFStN(i);
    bStrN[i] = m_EventData->GetBStN(i); // only 1 for BB10
    uppE [i] = m_EventData->GetUpE(i);
    dwnE [i] = m_EventData->GetDwE(i);
    sumE [i] = m_EventData->GetFrE(i);

    sPosArr.push_back(m_StripPos[detN[i]-1][fStrN[i]-1][bStrN[i]-1]); // detector and strip number from 1 (not 0)
    
    if (type[i] == 0) { // for X6
      ////////////////////////////////////////
      // hit position with resistive strip
      //
      // 1. get the distance from the resistive strip center
      Double_t distOnStrip = (uppE[i] - dwnE[i])/sumE[i]*X6_activeY/2.;
      // 2. calculate strip center position
      // = average(strip center of ohmic #2 & strip center of ohmic #3)
      TVector3 stripCenter = m_StripPos[detN[i]-1][fStrN[i]-1][1];
      stripCenter += m_StripPos[detN[i]-1][fStrN[i]-1][2];
      stripCenter *= 0.5;
      // 3. final position = stripCenter + detOrientation(Y-axis)*distOnStrip
      TVector3 hPos = m_DetOri[detN[i]-1];
      hPos *= distOnStrip;
      hPos += stripCenter;
      hPosArr.push_back(hPos);
      ////////////////////////////////////////
    }
    else { // for BB10 & QQQ5
      hPosArr.push_back(m_StripPos[detN[i]-1][fStrN[i]-1][bStrN[i]-1]); // detector and strip number from 1 (not 0) 
    }}
  ////////////////////////////////////////////////////////////  

  return;}

///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::Clear() {
  nhit = 0;
  sPosArr.clear();
  hPosArr.clear(); }

///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("STARKjr");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 

  vector<string> cart = {"Type", "POS"};
  vector<string> sphe = {"Type", "R","Theta","Phi"};
  vector<string> cyld = {"Type", "Rho","Phi","Z"};

  for (unsigned int i = 0 ; i < blocks.size() ; i++){
    ////////////////////////////////////////////////////////////
    // Cartesian coordinate
    if (blocks[i]->HasTokenList(cart)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARKjr " << i+1 <<  endl;    
      string Type = blocks[i]->GetString("Type");
      TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
      int Flip = 0;
      if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
      double Beta = 0;
      if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");
      AddDetector(Type, Pos, Flip, Beta); }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Spherical coordinate
    else if (blocks[i]->HasTokenList(sphe)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARKjr " << i+1 <<  endl;
      string Type = blocks[i]->GetString("Type");
      double R = blocks[i]->GetDouble("R","mm");
      double Theta = blocks[i]->GetDouble("Theta","deg");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      int Flip = 0;
      if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
      double Beta = 0;
      if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");

      TVector3 Pos;
      Pos.SetMagThetaPhi(R, Theta, Phi);
      AddDetector(Type, Pos, Flip, Beta); }
    ////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////
    // Cylindrical coordinate
    else if (blocks[i]->HasTokenList(cyld)){
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  STARKjr " << i+1 <<  endl;
      string Type = blocks[i]->GetString("Type");
      double Rho = blocks[i]->GetDouble("Rho","mm");
      double Phi = blocks[i]->GetDouble("Phi","deg");
      double Z   = blocks[i]->GetDouble("Z","mm");
      int Flip = 0;
      if (blocks[i]->HasToken("Flip")) Flip = blocks[i]->GetInt("Flip");
      double Beta = 0;
      if (blocks[i]->HasToken("Beta")) Beta = blocks[i]->GetDouble("Beta", "deg");

      TVector3 Pos;
      Pos.SetXYZ(Rho*cos(Phi),Rho*sin(Phi), Z);
      AddDetector(Type, Pos, Flip, Beta); }
    ////////////////////////////////////////////////////////////
    
    else {
      cout << "Error: check your input file formatting" << endl;
      exit(1); }
  }
  std::cout << "read complete" << std::endl;}


void TSTARKjrPhysics::AddParameterToCalibrationManager() {

}


///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::InitSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::FillSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::CheckSpectra() {
}



///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::ClearSpectra() {
  // To be done
}

///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::WriteSpectra() {
}





///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::InitializeRootInputRaw() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("STARKjr",true);
    inputChain->SetBranchAddress("STARKjr",&m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::InitializeRootInputPhysics() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("STARKjr", true);
    inputChain->SetBranchAddress("STARKjr", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TSTARKjrPhysics::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    TFile* outputFile = RootOutput::getInstance()->GetFile();
    //    outputTree->Branch("STARKjr", "TSTARKjrPhysics", &m_EventPhysics);
    outputFile->WriteObjectAny(m_EventPhysics,"TSTARKjrPhysics","STARKjr");

    outputTree->Branch("nhit" ,&nhit,"nhit/I");
    outputTree->Branch("type" ,type, "type[nhit]/I");
    outputTree->Branch("detN" ,detN, "detN[nhit]/I");
    outputTree->Branch("fStrN",fStrN,"fStrN[nhit]/I");
    outputTree->Branch("bStrN",bStrN,"bStrN[nhit]/I");
    outputTree->Branch("uppE" ,uppE, "uppE[nhit]/D");
    outputTree->Branch("dwnE" ,dwnE, "dwnE[nhit]/D");
    outputTree->Branch("sumE" ,sumE, "sumE[nhit]/D");
    outputTree->Branch("sPos" ,&sPosArr);
    outputTree->Branch("hPos" ,&hPosArr);    
}

////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TSTARKjrPhysics::Construct() {
    // construct method to be passed to the DetectorFactory
    return (NPL::VDetector*) new TSTARKjrPhysics();
}

////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_STARKjr{
    // register the construct method to the factory
public:
    proxy_STARKjr(){
        NPL::DetectorFactory::getInstance()->AddToken("STARKjr","STARKjr");
        NPL::DetectorFactory::getInstance()->AddDetector("STARKjr",TSTARKjrPhysics::Construct);
    }
};

proxy_STARKjr p_STARKjr;
}
