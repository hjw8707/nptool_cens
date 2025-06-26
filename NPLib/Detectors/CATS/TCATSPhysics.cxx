/* Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sandra Giron  contact address: giron@ipno.in2p3.fr       *
 *                                                                           *
 * Creation Date  : febuary 2010                                             *
 * Last update    : modification november 2011 by Pierre Morfouace			 *
 * Contact adress : morfouac@ipno.in2p3.fr									 *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CATS treated data                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "TCATSPhysics.h"
using namespace CATS_LOCAL;

//	STL
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>
using namespace std;
//	NPL
#include "RootInput.h"
#include "NPDetectorFactory.h"
#include "RootOutput.h"
#include "NPOptionManager.h"
//	ROOT
#include "TChain.h"
#include "TF1.h"
#include "TGraph.h"
ClassImp(TCATSPhysics)

  ///////////////////////////////////////////////////////////////////////////
  TCATSPhysics::TCATSPhysics(){
    m_EventData 				= new TCATSData	;
    m_PreTreatedData    = new TCATSData ;
    m_EventPhysics 			= this			    ;
    m_NumberOfCATS      = 0             ;
    m_Spectra           = NULL          ;
    m_Zproj             = 0             ;
  }

///////////////////////////////////////////////////////////////////////////
TCATSPhysics::~TCATSPhysics(){
}

///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::PreTreat(){
  ClearPreTreatedData();   
  gRandom->SetSeed(0);
  // X
  unsigned int sizeX = m_EventData->GetCATSMultX();
  for(unsigned int i = 0 ; i < sizeX ; i++){
    // Valid Channel X
    if(IsValidChannel("X", m_EventData->GetCATSDetX(i), m_EventData->GetCATSStripX(i)) ){
      if( fCATS_Threshold_X(m_EventData , i) ){
        double QX = fCATS_X_Q(m_EventData , i);
        unsigned int stripX = -1;
        //Inversion X
        if( *(m_CATSXInversion[m_EventData->GetCATSDetX(i)-1].begin() + m_EventData->GetCATSStripX(i)-1) != m_EventData->GetCATSStripX(i) ){
          stripX = *(m_CATSXInversion[m_EventData->GetCATSDetX(i)-1].begin() + m_EventData->GetCATSStripX(i)-1); 
        }
        else {
          stripX =  m_EventData->GetCATSStripX(i);
        }
        if(QX < 0)
          QX = 0;
        m_PreTreatedData->SetStripX( m_EventData->GetCATSDetX(i) , stripX, (unsigned int)QX);
      }
    }
  }

  // Y
  unsigned int sizeY = m_EventData->GetCATSMultY();
  for(unsigned int i = 0 ; i < sizeY ; i++){
    // Valid Channel Y
    if(IsValidChannel("Y", m_EventData->GetCATSDetY(i), m_EventData->GetCATSStripY(i))){
      if( fCATS_Threshold_Y(m_EventData , i) ){
        double QY = fCATS_Y_Q(m_EventData , i);
        unsigned int stripY = -1;
        //Inversion Y
        if( *(m_CATSYInversion[m_EventData->GetCATSDetY(i)-1].begin() + m_EventData->GetCATSStripY(i)-1) != m_EventData->GetCATSStripY(i) ){
          stripY = *(m_CATSYInversion[m_EventData->GetCATSDetY(i)-1].begin() + m_EventData->GetCATSStripY(i)-1); 
        }
        else {
          stripY = m_EventData->GetCATSStripY(i) ;
        }
        if(QY < 0)
          QY = 0;
        m_PreTreatedData->SetStripY( m_EventData->GetCATSDetY(i), stripY, (unsigned int)QY );
      }
    }
  }
  return;
}

/////////////////////////////////////////////////////////////////////////////
void TCATSPhysics::BuildSimplePhysicalEvent(){
  BuildPhysicalEvent();
}

//////////////////////////////////////////////////////////////////////////////		
void TCATSPhysics::BuildPhysicalEvent(){

  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }


  PreTreat();


  sizeX = m_PreTreatedData->GetCATSMultX() ;
  for( unsigned short i = 0 ; i < sizeX; i++ ){ 
    // Insert detector number in the set, if the key already exist, do nothing
    DetectorHitX.insert(m_PreTreatedData->GetCATSDetX(i));
  }

  // std::cout << "test 3" << std::endl;
  // Correspond to CATS with both X and Y
  sizeY = m_PreTreatedData->GetCATSMultY() ;
  for( unsigned short i = 0 ; i < sizeY ; i++ ){ 
    // Insert detector number in the set, if the key already exist, do nothing
    // Only if the detector was hit on X as well
    if(DetectorHitX.find(m_PreTreatedData->GetCATSDetY(i))!=DetectorHitX.end())
      DetectorHit.insert(m_PreTreatedData->GetCATSDetY(i));
  }
  // The number of CATS hit, i.e. the number of CATS that we are going to analyse


    
  //std::cout << "StripX mult " << sizeX << std::endl;
  for( unsigned short i = 0 ; i < sizeX; i++ ){ 
    StrX					         = m_PreTreatedData->GetCATSStripX(i);
    NX						         = m_PreTreatedData->GetCATSDetX(i);
    CATS_X_Q				       = m_PreTreatedData->GetCATSChargeX(i) ;
    if(DetectorHit.find(NX)!=DetectorHit.end()){
      MapX[NX].push_back(std::make_pair(StrX,CATS_X_Q));
      QSumX[NX]+= CATS_X_Q;
      if(MaxQX.find(NX)==MaxQX.end()|| MaxQX[NX].second < CATS_X_Q ){
        MaxQX[NX] = make_pair(StrX,CATS_X_Q);
      }
    }
  }
  for( unsigned short i = 0 ; i < sizeY; i++ ){ 
    StrY					         = m_PreTreatedData->GetCATSStripY(i);
    NY						         = m_PreTreatedData->GetCATSDetY(i);
    CATS_Y_Q				       = m_PreTreatedData->GetCATSChargeY(i) ;
    if(DetectorHit.find(NY)!=DetectorHit.end()){
      MapY[NY].push_back(std::make_pair(StrY,CATS_Y_Q));
      QSumY[NY]+= CATS_Y_Q;
      if(MaxQY.find(NY)==MaxQY.end()|| MaxQY[NY].second < CATS_Y_Q ){
        MaxQY[NY] = make_pair(StrY,CATS_Y_Q);
      }
    }
  }

  // std::cout << "test 5" << std::endl;

  for(auto &DetN : DetectorHit){       
  // std::cout << "test 6" << std::endl;
    // Return the position in strip unit
    // Convention: the collected charge is atrributed to the center of the strip
    // (histogram convention) so that a reconstructed position for a single strip
    // goes from strip index -0.5 to strip index +0.5
     double PosX =  ReconstructionFunctionX[DetN](MaxQX[DetN],MapX[DetN], QSumX[DetN]);
  // std::cout << "test 7" << std::endl;
     double PosY =  ReconstructionFunctionY[DetN](MaxQY[DetN],MapY[DetN], QSumY[DetN]);
    //std::cout << "test " << std::reduce(QsumSample[DetN].begin(),QsumSample[DetN].end()) /(QsumSample[DetN]).size() << std::endl;

    // a shift - -1 is made to have PosX in between -0.5 and 27.5
    // for the following calculation of the position in the lab.
    PosX = PosX;
    PosY = PosY;

    // sx and sy are the X and Y strip number between which the PosX and PosY are
    int sx0 = (int) PosX;
    int sx1 = sx0+1; 
    int sy0 = (int) PosY;
    int sy1 = sy0+1; 
    // std::cout << StripPositionX[1][1][0] << std::endl;

    if(PosX>-1000 && PosY>-1000 && sx0 > -1 && sx1 < 28 && sy0 > -1  && sy1 < 28){
      // px and py are the x and y coordinate of strip sx and sy 
      StripNumberX.push_back(PosX);
      StripNumberY.push_back(PosY);
      DetNumber.push_back(DetN);   
      ChargeX.push_back(QSumX[DetN]);   
      ChargeY.push_back(QSumY[DetN]);
      double px0 = StripPositionX[DetN][sx0][sy0];
      double px1 = StripPositionX[DetN][sx1][sy1];

      double py0 = StripPositionY[DetN][sx0][sy0];
      double py1 = StripPositionY[DetN][sx1][sy1];

      // Positon [Detector] = <PosZ,<PosX,PosY>>
      Positions[DetN] = make_pair(StripPositionZ[DetN], make_pair(px0+(px1-px0)*(PosX-sx0),py0+(py1-py0)*(PosY-sy0)));
      // std::cout << PosX << " " << px0 << std::endl;
      PositionX.push_back(px0+(px1-px0)*(PosX-sx0));  
      PositionY.push_back(py0+(py1-py0)*(PosY-sy0));
      PositionZ.push_back(StripPositionZ[DetN]); 
    }
  }

  // At least two CATS need to gave back position in order to reconstruct on Target 
  if(Positions.size()==2){ 
      double t = (m_Zproj-Positions[1].first)/(Positions[2].first-Positions[1].first);
      PositionOnTargetX= Positions[1].second.first + (Positions[2].second.first - Positions[1].second.first)*t;
      PositionOnTargetY= Positions[1].second.second + (Positions[2].second.second - Positions[1].second.second)*t;
    if(Mask1_Z != 0 && Mask2_Z != 0)
     { 
      double tmask1 = (Positions[1].first-Mask1_Z)/(Positions[2].first - Positions[1].first);
      double tmask2 = (Positions[2].first-Mask2_Z)/(Positions[2].first - Positions[1].first);
      
      PositionOnMask1X=  Positions[1].second.first - (Positions[2].second.first  -Positions[1].second.first)*tmask1;
      PositionOnMask1Y=  Positions[1].second.second - (Positions[2].second.second -Positions[1].second.second)*tmask1;
      PositionOnMask2X=  Positions[2].second.first - (Positions[2].second.first  -Positions[1].second.first)*tmask2;
      PositionOnMask2Y=  Positions[2].second.second - (Positions[2].second.second -Positions[1].second.second)*tmask2;
     }
     else{
      PositionOnMask1X= -1000;
      PositionOnMask1Y= -1000;
      PositionOnMask2X= -1000;
      PositionOnMask2Y= -1000;
     }
    }
  else{
    BeamDirection = TVector3 (1,0,0);
    PositionOnTargetX = -1000	;
    PositionOnTargetY = -1000	;
      PositionOnMask1X= -1000;
      PositionOnMask1Y= -1000;
      PositionOnMask2X= -1000;
      PositionOnMask2Y= -1000;
  }
  BeamDirection = GetBeamDirection();

  // Does not meet the conditions for target position and beam direction 
  return;
}

void TCATSPhysics::TestCalibration(){
  // Method to be used to test calibrations
  // Thus does not match X and Y strips
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }

  PreTreat();

  sizeX = m_PreTreatedData->GetCATSMultX() ;
  sizeY = m_PreTreatedData->GetCATSMultY() ;
  
  for( unsigned short i = 0 ; i < sizeX; i++ ){ 
    StripNumberX.push_back(m_PreTreatedData->GetCATSStripX(i));
    DetNumber.push_back(m_PreTreatedData->GetCATSDetX(i));
    ChargeX.push_back(m_PreTreatedData->GetCATSChargeX(i));
    StripNumberY.push_back(-1);
    ChargeY.push_back(-1);
  }
  for( unsigned short i = 0 ; i < sizeY; i++ ){ 
    StripNumberY.push_back(m_PreTreatedData->GetCATSStripY(i));
    DetNumber.push_back(m_PreTreatedData->GetCATSDetY(i));
    ChargeY.push_back(m_PreTreatedData->GetCATSChargeY(i));
    StripNumberX.push_back(-1);
    ChargeX.push_back(-1);
  }
  BeamDirection = TVector3 (1,0,0);
  PositionOnTargetX = -1000	;
  PositionOnTargetY = -1000	;
  PositionOnMask1X= -1000;
  PositionOnMask1Y= -1000;
  PositionOnMask2X= -1000;
  PositionOnMask2Y= -1000;
}

void TCATSPhysics::SetTreeReader(TTreeReader* TreeReader) {
   TCATSPhysicsReader::r_SetTreeReader(TreeReader);
 }
///////////////////////////////////////////////////////////////////////////
//	Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
void TCATSPhysics::ReadConfiguration(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("CATSDetector");
  vector<NPL::InputBlock*> blocksMask = parser.GetAllBlocksWithToken("MASK");
  
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detectors found " << endl; 
    cout << "//// " << blocksMask.size() << " masks found " << endl; 

  vector<string> token = {"X1_Y1","X28_Y1","X1_Y28","X28_Y28","CATSNumber"};
  vector<string> tokenMask = {"Z","MaskNumber"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      TVector3 A = blocks[i]->GetTVector3("X1_Y1","mm");
      TVector3 B = blocks[i]->GetTVector3("X28_Y1","mm");
      TVector3 C = blocks[i]->GetTVector3("X1_Y28","mm");
      TVector3 D = blocks[i]->GetTVector3("X28_Y28","mm");
      UShort_t N = blocks[i]->GetInt("CATSNumber");

      AddCATS(A,B,C,D,N);
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }
  for(unsigned int i = 0 ; i < blocksMask.size() ; i++){
    if(blocksMask[i]->HasTokenList(tokenMask)){
      AddMask(blocksMask[i]->GetDouble("Z","mm"),blocksMask[i]->GetInt("MaskNumber"));
    }
    else{
      cout << "ERROR: check your input file formatting " << endl;
      exit(1);
    }
  }

  InitializeStandardParameter();
  ReadAnalysisConfig();
}

/////////////////////////////////////////////////////////////////////
//	Activated associated Branches and link it to the private member DetectorData address
//	In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
void TCATSPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain()	;
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
    std::cout << "////////////////////////////// TEST" << std::endl;
    inputChain->SetBranchStatus( "CATS" , true )			;
    inputChain->SetBranchStatus( "fCATS_*" , true )		;
    inputChain->SetBranchAddress( "CATS" , &m_EventData )           ;
  }
}

/////////////////////////////////////////////////////////////////////
//   Activated associated Branches and link it to the private member DetectorPhysics address
//   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
void TCATSPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
    inputChain->SetBranchStatus( "CATS" , true );
    inputChain->SetBranchStatus( "DetNumberX" , true );
    inputChain->SetBranchStatus( "StripX" , true );
    inputChain->SetBranchStatus( "ChargeX" , true );
    inputChain->SetBranchStatus( "StripMaxX" , true );
    inputChain->SetBranchStatus( "DetNumberY" , true );
    inputChain->SetBranchStatus( "StripY" , true );
    inputChain->SetBranchStatus( "ChargeY" , true );
    inputChain->SetBranchStatus( "StripMaxY" , true );
    inputChain->SetBranchStatus( "DetMaxX" , true );
    inputChain->SetBranchStatus( "DetMaxY" , true );
    inputChain->SetBranchStatus( "PositionX" , true );
    inputChain->SetBranchStatus( "PositionY" , true );
    inputChain->SetBranchStatus( "PositionZ" , true );
    inputChain->SetBranchStatus( "StripNumberX" , true );
    inputChain->SetBranchStatus( "StripNumberY" , true );
    inputChain->SetBranchStatus( "PositionOnTargetX" , true );
    inputChain->SetBranchStatus( "PositionOnTargetY" , true );
    inputChain->SetBranchStatus( "QsumX" , true );
    inputChain->SetBranchStatus( "QsumY" , true );
    inputChain->SetBranchAddress( "CATS" , &m_EventPhysics );
  }
}

/////////////////////////////////////////////////////////////////////
//	Create associated branches and associated private member DetectorPhysics address
void TCATSPhysics::InitializeRootOutput(){
  TTree* outputTree = RootOutput::getInstance()->GetTree()		;
  outputTree->Branch( "CATS" , "TCATSPhysics" , &m_EventPhysics )	;
}

/////////////////////////////////////////////////////////////////////
void TCATSPhysics::AddCATS(TVector3 C_X1_Y1, TVector3 C_X28_Y1, TVector3 C_X1_Y28, TVector3 C_X28_Y28, UShort_t N){
  m_NumberOfCATS++			;

  // remove warning
  C_X28_Y28 *= 1;

  //	Vector U on Telescope Face (paralelle to Y Strip) 
  TVector3 U = C_X28_Y1 - C_X1_Y1 				;	
  U = U.Unit()									;

  //	Vector V on Telescope Face (parallele to X Strip)
  TVector3 V = C_X1_Y28 - C_X1_Y1 				;
  V = V.Unit()									;

  //	Position Vector of Strip Center
  TVector3 StripCenter 					;
  //	Position Vector of X=1 Y=1 Strip 
  TVector3 Strip_1_1 					;		

  //	Geometry Parameters  
  double Face = 71.12 				; //mm
  double NumberOfStrip = 28 			;
  double StripPitch = Face / NumberOfStrip	; //mm

  //	Buffer object to fill Position Array
  vector<double> lineX ; 
  vector<double> lineY ; 

  vector< vector< double > >	OneDetectorStripPositionX	;
  vector< vector< double > >	OneDetectorStripPositionY	;
  double                 	OneDetectorStripPositionZ	;

  //	Moving StripCenter to 1.1 corner (strip center!) :
  Strip_1_1 = C_X1_Y1 + (U+V) * (StripPitch/2) 	;

  for( int i = 0 ; i < 28 ; i++ ){
    lineX.clear()	;
    lineY.clear()	;

    for( int j = 0 ; j < 28 ; j++ ){
      StripCenter  = Strip_1_1 + StripPitch*( i*U + j*V  )	;
      lineX.push_back( StripCenter.x() )	;
      lineY.push_back( StripCenter.y() )	;
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
  }

  OneDetectorStripPositionZ = C_X1_Y1.Z();

  StripPositionX[N] = OneDetectorStripPositionX	;
  StripPositionY[N] = OneDetectorStripPositionY	;
  StripPositionZ[N] = OneDetectorStripPositionZ	;

}

void TCATSPhysics::AddMask(Double_t Z, UShort_t MaskNumber){
  if(MaskNumber == 1){
    Mask1_Z = Z;
  }
  else if(MaskNumber == 2){
    Mask2_Z = Z;
  }
  else{
    std::cout << "Wrong Number for MASKS" << std::endl;
  }
}

///////////////////////////////////////////////////////////////
void TCATSPhysics::Clear(){  
  DetNumber.clear(); 
  ChargeX.clear(); 
  PositionX.clear();
  ChargeY.clear(); 
  PositionY.clear();
  PositionZ.clear();

  StripNumberX.clear();
  StripNumberY.clear();
  QSumX.clear();
  QSumY.clear();
  Positions.clear();
  MapX.clear();
  MapY.clear();
  MaxQX.clear();
  MaxQY.clear();

  DetectorHit.clear();
  DetectorHitX.clear();
}

////////////////////////////////////////////////////////////////////////////
bool TCATSPhysics :: IsValidChannel(const string& DetectorType, const int& Detector , const int& channel) {
  if(DetectorType == "X")
    return *(m_XChannelStatus[Detector-1].begin()+channel-1);

  else if(DetectorType == "Y")
    return *(m_YChannelStatus[Detector-1].begin()+channel-1);

  else return false;
}


///////////////////////////////////////////////////////////////////////////////////
void TCATSPhysics::InitializeStandardParameter(){
  //   Enable all channel and no inversion
  vector< bool > ChannelStatus;
  vector< int > InversionStatus;
  m_XChannelStatus.clear()    ;
  m_YChannelStatus.clear()    ;
  m_CATSXInversion.clear()   ;
  m_CATSYInversion.clear()   ;

  ChannelStatus.resize(28,true);
  InversionStatus.resize(28);
  for(unsigned int j = 0 ; j < InversionStatus.size() ; j++){
    InversionStatus[j] = j+1;
  }

  for(int i = 0 ; i < m_NumberOfCATS ; ++i){
    m_XChannelStatus[i] = ChannelStatus;
    m_YChannelStatus[i] = ChannelStatus;
    m_CATSXInversion[i] = InversionStatus;
    m_CATSYInversion[i] = InversionStatus;
    SetReconstructionMethod(i+1, "X", "AGAUSS");
    SetReconstructionMethod(i+1, "Y", "AGAUSS");
  }

  return;
}   

///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::ReadAnalysisConfig(){
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigCATS.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigCATS.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigCATS.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigCATS.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 10, "ConfigCATS") == 0) ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {

      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo == "DISABLE_CHANNEL") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl;
        int Detector = atoi(DataBuffer.substr(4,1).c_str());
        int channel = -1;
        if (DataBuffer.compare(5,4,"STRX") == 0) {
          channel = atoi(DataBuffer.substr(9).c_str());
          *(m_XChannelStatus[Detector-1].begin()+channel-1) = false;
        }

        else if (DataBuffer.compare(5,4,"STRY") == 0) {
          channel = atoi(DataBuffer.substr(9).c_str());
          *(m_YChannelStatus[Detector-1].begin()+channel-1) = false;
        }

        else cout << "Warning: detector type for CATS unknown!" << endl;

      }

      else if (whatToDo == "INVERSION") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer;
        int Detector = atoi(DataBuffer.substr(4,1).c_str());
        int channel1 = -1;
        int channel2 = -1;
        AnalysisConfigFile >> DataBuffer;
        cout << " " << DataBuffer;
        if (DataBuffer.compare(0,4,"STRX") == 0) {
          channel1 = atoi(DataBuffer.substr(4).c_str());
          AnalysisConfigFile >> DataBuffer;
          cout << " " << DataBuffer << endl;
          channel2 = atoi(DataBuffer.substr(4).c_str());
          *(m_CATSXInversion[Detector-1].begin()+channel1-1) = channel2;
          *(m_CATSXInversion[Detector-1].begin()+channel2-1) = channel1;
        }
      }
      
      else if (whatToDo == "INVERSIONX") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer;
        int Detector = atoi(DataBuffer.substr(4,1).c_str());
        for(unsigned int strip = 0; strip < 28; strip ++){
          *(m_CATSXInversion[Detector-1].begin()+strip) = 27-strip;
        }
      }
      
      else if (whatToDo == "INVERSIONY") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer;
        int Detector = atoi(DataBuffer.substr(4,1).c_str());
        for(unsigned int strip = 0; strip < 28; strip ++){
          *(m_CATSYInversion[Detector-1].begin()+strip) = 27-strip;
        }
      }

      else if (whatToDo == "RECONSTRUCTION_METHOD") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer ;
        // DataBuffer is of form CATSNX 
        // Look for the CATS Number removing the first 4 letters and the trailling letter
        string Duplicate = DataBuffer.substr(4); // Duplicate is of form NX
        Duplicate.resize(Duplicate.size()-1); // Duplicate is of form
        unsigned int CATSNumber =  atoi(Duplicate.c_str());

        // Look for the X or Y part of the Detector, Basically the last character
        string XorY(string(1,DataBuffer[DataBuffer.size()-1])) ; 

        // Get the Reconstruction Methods Name
        AnalysisConfigFile >> DataBuffer;       
        cout << " " << DataBuffer << endl ;
        // Set the Reconstruction Methods using above information 
        SetReconstructionMethod(CATSNumber,XorY,DataBuffer);
      }

      else if (whatToDo == "ZPROJ") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl ;
        m_Zproj = atoi(DataBuffer.c_str());
      }

      else {ReadingStatus = false;}

    }
  }
} 
///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::InitSpectra(){
  m_Spectra = new TCATSSpectra(m_NumberOfCATS);
}

///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::FillSpectra(){
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::CheckSpectra(){
  // To be done
}
///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::ClearSpectra(){
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TCATSPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();

  else{
    map< string , TH1*> empty; 
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TCATSPhysics::WriteSpectra(){
  if(m_Spectra)
    m_Spectra->WriteSpectra();
}

/////////////////////////////////////////////////////////////////////
//	Add Parameter to the CalibrationManger
void TCATSPhysics::AddParameterToCalibrationManager()	{
  CalibrationManager* Cal = CalibrationManager::getInstance();
  for(int i = 0 ; i < m_NumberOfCATS ; i++){
    for( int j = 0 ; j < 28 ; j++){
      Cal->AddParameter("CATS", "D"+NPL::itoa(i+1)+"_X"+NPL::itoa(j+1)+"_Q","CATS_D"+NPL::itoa(i+1)+"_X"+NPL::itoa(j+1)+"_Q")	;
      Cal->AddParameter("CATS", "D"+NPL::itoa(i+1)+"_Y"+NPL::itoa(j+1)+"_Q","CATS_D"+NPL::itoa(i+1)+"_Y"+NPL::itoa(j+1)+"_Q")	;
      Cal->AddParameter("CATS", "D"+NPL::itoa(i+1)+"_X"+NPL::itoa(j+1),"CATS_D"+NPL::itoa(i+1)+"_X"+NPL::itoa(j+1))	;
      Cal->AddParameter("CATS", "D"+NPL::itoa(i+1)+"_Y"+NPL::itoa(j+1),"CATS_D"+NPL::itoa(i+1)+"_Y"+NPL::itoa(j+1))	;
    } 
  }

  return;
}	

////////////////////////////////////////////////////////////////
void TCATSPhysics::SetReconstructionMethod(unsigned int CATSNumber, string XorY, string MethodName){
  if(XorY=="X"){

    if(MethodName=="ASECH") ReconstructionFunctionX[CATSNumber] = &(AnalyticHyperbolicSecant);
    else if(MethodName=="FSECH") ReconstructionFunctionX[CATSNumber] = &(FittedHyperbolicSecant);
    //else if(MethodName=="AGAUSS") ReconstructionFunctionX[CATSNumber] = &(AnalyticGaussian);
    //else if(MethodName=="CENTROIDE")  ReconstructionFunctionX[CATSNumber] = &(Centroide); 
    else cout <<"WARNING: Wrong name for reconsctuction Method, using default AGAUSS" << endl;
  }

  if(XorY=="Y"){

    if(MethodName=="ASECH") ReconstructionFunctionY[CATSNumber] = &(AnalyticHyperbolicSecant);
    else if(MethodName=="FSECH") ReconstructionFunctionY[CATSNumber] = &(FittedHyperbolicSecant);
    //else if(MethodName=="AGAUSS") ReconstructionFunctionY[CATSNumber] = &(AnalyticGaussian);
    //else if(MethodName=="CENTROIDE")  ReconstructionFunctionY[CATSNumber] = &(Centroide); 
    else cout <<"WARNING: Wrong name for reconsctuction Method, using default AGAUSS" << endl;
  }

}

///////////////////////////////////////////////////////////////
TVector3 TCATSPhysics::GetBeamDirection(){
  TVector3 Direction;
  if(Positions.size() <2)return Direction;

  else{
    Direction = TVector3 (Positions[2].second.first-Positions[1].second.first ,
    Positions[2].second.second-Positions[1].second.second,
    Positions[2].first-Positions[1].first );
    Direction.Unit();
  }


  return(Direction) ;	
}

///////////////////////////////////////////////////////////////
TVector3 TCATSPhysics::GetPositionOnTarget(){
  double Pi = 3.14159265;
  TVector3 Position = TVector3 (GetPositionOnTargetX() 	,
      GetPositionOnTargetY() 	,
      GetPositionOnTargetX()*tan(m_TargetAngle*Pi/180)); 
  return(Position) ;	
}

////////////////////////////////////////////////////////////////////////
namespace CATS_LOCAL{
  ////////////////////////////////////////////////////////////////////
  /*double AnalyticGaussian(std::pair<UShort_t,UShort_t>& MaxQ,std::vector<std::pair<UShort_t,UShort_t>>& Map, Double_t QSum){
    /*double gauss = -1000;  
    double Q[3];
    double StripPos[3];
    for(int j = 0; j<3 ; j++){
      Q[j] = 0;
      StripPos[j] = 0;
    }

    if(MaxQ.first> 3 && MaxQ.first< 26){
      // central value taken using the Strip with Max charge
      Q[0] = Buffer_Q[MaxQ.first-1] ;

      // Look at the next strip on the left
      if(Buffer_Q[MaxQ.first-2]!=-1){
        Q[1] = Buffer_Q[MaxQ.first-2];
        StripPos[1] = MaxQ.first-2;
      }

      // Look at the next next strip on the left
      else if(Buffer_Q[MaxQ.first-3]!=-1){
        Q[1] = Buffer_Q[MaxQ.first-3];
        StripPos[1] = MaxQ.first-3;
      }

      // Look at the next next next strip on the left
      else if(Buffer_Q[MaxQ.first-4]!=-1){
        Q[1] = Buffer_Q[MaxQ.first-4];
        StripPos[1] = MaxQ.first-4;
      }

      // Look at the next strip on the right
      if(Buffer_Q[MaxQ.first]!=-1){
        Q[2] = Buffer_Q[MaxQ.first];
        StripPos[2] = MaxQ.first;
      }

      // Look at the next next strip on the right
      else if(Buffer_Q[MaxQ.first+1]!=-1){
        Q[2] = Buffer_Q[MaxQ.first+1];
        StripPos[2] = MaxQ.first+1;
      }

      // Look at the next next next strip on the right 
      else if(Buffer_Q[MaxQ.first+2]!=-1){
        Q[2] = Buffer_Q[MaxQ.first+2];
        StripPos[2] = MaxQ.first+2;
      }

    }

    double Q0_Q1 = log(Q[0]/Q[1]);
    double Q0_Q2 = log(Q[0]/Q[2]);

    double num   = Q0_Q1 * (StripPos[2]*StripPos[2] - StripPos[0]*StripPos[0]) - Q0_Q2 * (StripPos[1]*StripPos[1] - StripPos[0]*StripPos[0])  ;
    double denom = Q0_Q1 * (StripPos[2] - StripPos[0]) - Q0_Q2 * (StripPos[1] - StripPos[0]) ;

    if(denom != 0){
      gauss = 0.5*num / denom;
    }

    else{
      gauss = -1000;
    }

    return gauss;
  return 0;
  }

  ///////////////////////////////////////////////////////////////
  double Centroide(std::pair<UShort_t,UShort_t>& MaxQ,std::vector<std::pair<UShort_t,UShort_t>>& Map, Double_t QSum){
    double Centroide = 0 ;
    double ChargeTotal = 0;

    unsigned int sizeQ = Buffer_Q.size(); 
    for(unsigned int i = 0 ; i < sizeQ ; i++){
      if(Buffer_Q[i]>0){ 
        Centroide += (i)*Buffer_Q[i-1] ;
        ChargeTotal+=Buffer_Q[i-1];
      }
    }

    if(ChargeTotal>0) Centroide = Centroide / ChargeTotal ;

    else {
      Centroide = -1000 ; 
    } 

    return Centroide ;
  return 0;
  }
*/
  /////////////////////////////////////////////////////////////////////
  double AnalyticHyperbolicSecant(std::pair<UShort_t,UShort_t>& MaxQ,std::vector<std::pair<UShort_t,UShort_t>>& Map,Double_t QSum){
    double sech = -1000 ;

  // std::cout << "test AnH 1" << std::endl;
    if(MaxQ.second > 0 && MaxQ.first > 2 && MaxQ.first<27){	
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

  /////////////////////////////////////////////////////////////////////
  double FittedHyperbolicSecant(std::pair<UShort_t,UShort_t>& MaxQ,std::vector<std::pair<UShort_t,UShort_t>>& Map, Double_t QSum){
    // Warning: should not delete static variable
    static TF1* f = new TF1("sechs","[0]/(cosh(TMath::Pi()*(x-[1])/[2])*cosh(TMath::Pi()*(x-[1])/[2]))",1,28);

    // Help the fit by computing the position of the maximum by analytic method
    double StartingPoint = AnalyticHyperbolicSecant(MaxQ,Map,QSum);
    // if analytic method fails then the starting point in strip max
    if(StartingPoint==-1000) StartingPoint = MaxQ.first; 

    // Maximum is close to charge max, Mean value is close to Analytic one, typical width is 3.8 strip
    f->SetParameters(MaxQ.second,StartingPoint,3.8);

    static vector<double> y ;
    static vector<double> q ; 
    y.clear(); q.clear();
    double final_size = 0 ;
    unsigned int sizeQ = Map.size(); 

    for(unsigned int i = 0 ; i < sizeQ ; i++){
      if(Map[i].second > MaxQ.second*0.2){
        q.push_back(Map[i].second);
        y.push_back(i+1);
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
    delete f;

  return 0;
  }

  ////////////////////////////////////////////////////////////////////////
  double fCATS_X_Q(const TCATSData* m_EventData , const int& i){
    static string name; 
    name = "CATS/D" ;
    name+= NPL::itoa( m_EventData->GetCATSDetX(i) ) ;
    name+= "_X" ;
    name+= NPL::itoa( m_EventData->GetCATSStripX(i) ) ;
    name+= "_Q";
    return CalibrationManager::getInstance()->ApplyCalibration( name,   
        m_EventData->GetCATSChargeX(i) + gRandom->Rndm() );
    //return CalibrationManager::getInstance()->ApplyCalibration( name,   
    //    m_EventData->GetCATSChargeX(i) + gRandom->Rndm() );
    //m_EventData->GetCATSChargeX(i) + gRandom->Rndm() - fCATS_Ped_X(m_EventData, i) );
  }
  ////////////////////////////////////////////////////////////////////////
  double fCATS_Y_Q(const TCATSData* m_EventData , const int& i){
    static string name; 
    name = "CATS/D" ;
    name+= NPL::itoa( m_EventData->GetCATSDetY(i) ) ;
    name+= "_Y" ;
    name+= NPL::itoa( m_EventData->GetCATSStripY(i) ) ;
    name+= "_Q";
    return CalibrationManager::getInstance()->ApplyCalibration( name ,   
        m_EventData->GetCATSChargeY(i) + gRandom->Rndm() );
    //return CalibrationManager::getInstance()->ApplyCalibration( name ,   
    //    m_EventData->GetCATSChargeY(i) + gRandom->Rndm() );
    //m_EventData->GetCATSChargeY(i) + gRandom->Rndm() - fCATS_Ped_Y(m_EventData, i) );
  }
  ////////////////////////////////////////////////////////////////////////
  bool fCATS_Threshold_X(const TCATSData* m_EventData , const int& i){
    static string name; 
    name = "CATS/D" ;
    name+= NPL::itoa( m_EventData->GetCATSDetX(i) ) ;
    name+= "_X" ;
    name+= NPL::itoa( m_EventData->GetCATSStripX(i) );
    return CalibrationManager::getInstance()->ApplyThreshold(name,
        m_EventData->GetCATSChargeX(i));
  }
  ////////////////////////////////////////////////////////////////////////
  bool fCATS_Threshold_Y(const TCATSData* m_EventData , const int& i){
    static string name; 
    name ="CATS/D" ;
    name+= NPL::itoa( m_EventData->GetCATSDetY(i) ) ;
    name+= "_Y" ;
    name+= NPL::itoa( m_EventData->GetCATSStripY(i) );
    return CalibrationManager::getInstance()->ApplyThreshold( name,
        m_EventData->GetCATSChargeY(i));
  }
  ////////////////////////////////////////////////////////////////////////
  double fCATS_Ped_X(const TCATSData* m_EventData, const int& i){
    static string name; 
    name =  "CATS/D" ;
    name+= NPL::itoa( m_EventData->GetCATSDetX(i) ) ;
    name+= "_X" ;
    name+= NPL::itoa( m_EventData->GetCATSStripX(i) ) ;
    return CalibrationManager::getInstance()->GetPedestal(name);
  }
  ////////////////////////////////////////////////////////////////////////
  double fCATS_Ped_Y(const TCATSData* m_EventData, const int& i){
    static string name; 
    name = "CATS/D" ;
    name+= NPL::itoa( m_EventData->GetCATSDetY(i) ) ;
    name+= "_Y" ;
    name+= NPL::itoa( m_EventData->GetCATSStripY(i) );
    return CalibrationManager::getInstance()->GetPedestal( name );
  }
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TCATSPhysics::Construct(){
  return (NPL::VDetector*) new TCATSPhysics();
}

NPL::VTreeReader* TCATSPhysics::ConstructReader() { return (NPL::VTreeReader*)new TCATSPhysicsReader(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_cats{
    public:
      proxy_cats(){
        NPL::DetectorFactory::getInstance()->AddToken("CATSDetector","CATS");
        NPL::DetectorFactory::getInstance()->AddDetector("CATSDetector",TCATSPhysics::Construct);
        NPL::DetectorFactory::getInstance()->AddDetectorReader("CATSDetector", TCATSPhysics::ConstructReader);
      }
  };

  proxy_cats p;
}

