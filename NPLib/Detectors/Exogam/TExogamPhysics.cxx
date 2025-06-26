/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sandra GIRON  contact address: giron@ipno.in2p3.fr       *
 *                  Benjamin LE CROM		   lecrom@ipno.in2p3.fr              *
 * Creation Date  : march 2014                                               *
 * Last update    : updated in 2023-2024 by H. Jacob hjacob@ijclab.in2p3.fr  *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold exogam treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TExogamPhysics.h"
using namespace EXOGAM_LOCAL;

//	STL
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <functional>
#include <iomanip>

//	NPL
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPVDetector.h"
#include "RootInput.h"
#include "RootOutput.h"
//	ROOT
#include "TChain.h"

///////////////////////////////////////////////////////////////////////////
ClassImp(TExogamPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TExogamPhysics::TExogamPhysics()
  : m_PreTreatedData(new TExogamCalData),
  m_EventData(new TExogamData),
  TSEvent(new TimeStamp),
  m_EventPhysics(this)
  {
    // m_Spectra = NULL;
    m_EXO_E_RAW_Threshold = 0;
    m_EXO_E_Threshold = 0;
    m_EXO_EHG_RAW_Threshold = 0;
    m_EXO_TDC_RAW_Threshold = 0;
    m_ExoTDC_HighThreshold = 1e6;
    m_ExoTDC_LowThreshold = 0;
    m_EXO_OuterUp_RAW_Threshold = 60000;
  DataIsCal = false;
}
  
TExogamPhysics::~TExogamPhysics() {
  delete m_PreTreatedData;
  delete m_EventData;
  // delete TSEvent; 
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::PreTreat() {
  // Clearing PreTreat TExogamData
  ClearPreTreatedData();

  //E
  m_EXO_Mult = m_EventData->GetExoMult();

  for (unsigned int i = 0; i < m_EXO_Mult; ++i) {
      bool DoPreTreat = false;
      if(!RefTS_Name.empty()){
        
        std::string TSName = "EXO_"+std::to_string(m_EventData->GetExoCrystal(i));
        TSEvent->AddTimeStamp(TSName,m_EventData->GetExoTS(i));
        TSEvent->AddTimeStamp(RefTS_Name,RefTS);
        
        if(TSEvent->MatchTS(TSName)){
          DoPreTreat = true;
        }
        TSEvent->ClearTimeStamps();
      }
      
      // else, all datas are filled
      else{
        DoPreTreat = true;
      }
    if(DoPreTreat){ 
      ResetPreTreatVariable();
    
      if (m_EventData->GetExoE(i) > m_EXO_E_RAW_Threshold)
        EXO_E = fEXO_E(m_EventData, i);

      if (m_EventData->GetExoEHG(i) > m_EXO_EHG_RAW_Threshold)
        EXO_EHG = fEXO_EHG(m_EventData, i);
    
      if (m_EventData->GetExoTDC(i) > m_EXO_TDC_RAW_Threshold)
        EXO_TDC = fEXO_T(m_EventData, i);
 
      if (m_EventData->GetExoOuter1(i) < m_EXO_OuterUp_RAW_Threshold)
        EXO_Outer1 = fEXO_Outer(m_EventData, i, 0);
      else
        EXO_Outer1 = 0;
    
      if (m_EventData->GetExoOuter2(i) < m_EXO_OuterUp_RAW_Threshold)
        EXO_Outer2 = fEXO_Outer(m_EventData, i, 1);
      else
        EXO_Outer2 = 0;
    
      if (m_EventData->GetExoOuter3(i) < m_EXO_OuterUp_RAW_Threshold)
        EXO_Outer3 = fEXO_Outer(m_EventData, i, 2);
      else
        EXO_Outer3 = 0;
    
      if (m_EventData->GetExoOuter4(i) < m_EXO_OuterUp_RAW_Threshold)
        EXO_Outer4 = fEXO_Outer(m_EventData, i, 3);
      else
        EXO_Outer4 = 0;
    
      // *1000 to convert MeV into keV
      if(EXO_E > m_EXO_E_Threshold){
        m_PreTreatedData->SetExo(m_EventData->GetExoCrystal(i), EXO_E,
        EXO_EHG, m_EventData->GetExoTS(i), EXO_TDC, 
        m_EventData->GetExoBGO(i), m_EventData->GetExoCsI(i), EXO_Outer1,
        EXO_Outer2, EXO_Outer3, EXO_Outer4);
      }
    } 
  }
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ResetPreTreatVariable(){
  EXO_E = -1000;
  EXO_EHG = -1000;
  EXO_TDC = -1000;
  EXO_Outer1 = -1000;
  EXO_Outer2 = -1000;
  EXO_Outer3 = -1000;
  EXO_Outer4 = -1000;
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::BuildPhysicalEvent() {
  ClaimReaderData();
  
  // std::cout << m_EventData << std::endl;
  if(!DataIsCal)
    PreTreat();

  // This maps stores ID of events sorted by flange number. Map key is flange nbr, vector should contain ID of events
  std::map<unsigned int,std::vector<unsigned int>> HitsID;

  for(unsigned int i = 0; i < m_PreTreatedData->GetExoMult(); i++){
    // Asking good TDC prompt, which can be ignored in EXOGAM config files (TDC not working for some good crystals in E805) 
    if(TDCMatch(i) ||  find(IgnoreTDC.begin(), IgnoreTDC.end(), m_PreTreatedData->GetExoCrystal(i)) != IgnoreTDC.end()){
      // Doing flange and crystal matching
      flange_nbr = MapCrystalFlangeCLover[m_PreTreatedData->GetExoCrystal(i)].first;
      crystal_nbr = MapCrystalFlangeCLover[m_PreTreatedData->GetExoCrystal(i)].second;
      E.push_back(m_PreTreatedData->GetExoE(i));
      EHG.push_back(m_PreTreatedData->GetExoEHG(i));
      Outer1.push_back(m_PreTreatedData->GetExoOuter1(i));
      Outer2.push_back(m_PreTreatedData->GetExoOuter2(i));
      Outer3.push_back(m_PreTreatedData->GetExoOuter3(i));
      Outer4.push_back(m_PreTreatedData->GetExoOuter4(i));
      TDC.push_back(m_PreTreatedData->GetExoTDC(i));
      TS.push_back(m_PreTreatedData->GetExoTS(i));
      Flange.push_back(flange_nbr);
      Crystal.push_back(crystal_nbr);

      HitsID[flange_nbr].push_back(i);
    }
  }

  // Now that HitsID is full, we use it to process simple AddBack of events in the same flange
  // Basically looping on all flanges, then on al events ID in each flange
  for(auto it = HitsID.begin(); it != HitsID.end(); it++){
    double E_AddBack = 0;
    double E_Max = 0;
    unsigned int Id_Max = 0;
    for(auto itvec = (*it).second.begin(); itvec !=(*it).second.end(); itvec++){
      E_AddBack+= m_PreTreatedData->GetExoE(*itvec);
      if(E_Max < m_PreTreatedData->GetExoE(*itvec)){
        E_Max = m_PreTreatedData->GetExoE(*itvec);
        Id_Max = *itvec;
      }
    }
    // Doing it again for this loop, it's a bit unhappy but didnt find a better way to do it yet
    flange_nbr = (*it).first;
    crystal_nbr = MapCrystalFlangeCLover[m_PreTreatedData->GetExoCrystal(Id_Max)].second;

    // Adding all AddBack (AB) related stuff
    E_AB.push_back(E_AddBack);
    Flange_AB.push_back(flange_nbr);
    Size_AB.push_back((*it).second.size());
    TDC_AB.push_back(m_PreTreatedData->GetExoTDC(Id_Max));
    TS_AB.push_back(m_PreTreatedData->GetExoTS(Id_Max));

    // Adding these parameters for Doppler correction purposes (D)
    Crystal_AB.push_back(crystal_nbr);
    int MaxOuterId = GetMaxOuter(Id_Max);
    Outer_AB.push_back(GetMaxOuter(Id_Max));
    
    if(MaxOuterId > -1){
      Vector3D BeamDir(0.,0.,1.);
      Vector3D BeamImpact(0.,0.,0.);
        
      ExogamGeo->SetBeam(BeamDir, BeamImpact);
      Vector3D RecoilDir = BeamDir;
      ExogamGeo->SetRecoil(RecoilDir);
      ExogamGeo->SetGammaInteractionPoint(flange_nbr, crystal_nbr, GetMaxOuter(Id_Max), E_AddBack);
        
      Vector3D GammaVector = ExogamGeo->newton_raphson(initTheta,BeamImpact,nrThreshold,nrMaxIter);
      double Theta_ = GammaVector.theta();
      double Phi_ = GammaVector.phi();
      double Angle = RecoilDir.angle(GammaVector);
      ExogamGeo->ResetGamma();

      // Doppler correction is not performed here to give the freedom to the user 
      //to use a certain way to calculate the Beta depending on the experiment
      Theta.push_back(Theta_);
      Phi.push_back(Phi_);
      }
    else{
      Theta.push_back(-1000);
      Phi.push_back(-1000);
    }
  }
}

///////////////////////////////////////////////////////////////////////////
double TExogamPhysics::DopplerCorrection(double E, double Angle, double Beta) {
  // E of gamma Doppler shifted, Angle is the angle between gamma and part, Beta is speed of part 
  double E_corr = 0;
  double Gamma = 1./sqrt(1.-Beta*Beta);

  E_corr = Gamma*E*(1.-Beta*std::cos(Angle));

  return E_corr;
}

double TExogamPhysics::DopplerCorrection(double E, double ThetaGamma, double PhiGamma, double ThetaPart, double PhiPart, double Beta){
  // E of gamma Doppler shifted, Angle is the angle between gamma and part, Beta is speed of part 
  
	double Angle =std::acos(TMath::Sin(ThetaPart)*TMath::Cos(PhiPart)*TMath::Sin(ThetaGamma)*TMath::Cos(PhiGamma)+
		  	      TMath::Sin(ThetaPart)*TMath::Sin(PhiPart)*TMath::Sin(ThetaGamma)*TMath::Sin(PhiGamma)+
			      TMath::Cos(ThetaPart)*TMath::Cos(ThetaGamma));
  
  return DopplerCorrection(E, Angle, Beta);
}

void TExogamPhysics::ClaimReaderData() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }
}

///////////////////////////////////////////////////////////////////////////
bool TExogamPhysics::TDCMatch(unsigned int event){
  return m_PreTreatedData->GetExoTDC(event) > m_ExoTDC_LowThreshold && m_PreTreatedData->GetExoTDC(event) < m_ExoTDC_HighThreshold;
}

///////////////////////////////////////////////////////////////////////////
int TExogamPhysics::GetMaxOuter(unsigned int EventId){
  // somehow starting at 50 to get something equivalent to a 50keV threshold
  double OuterMax = 0.05;
  int OuterId = -1;
  if(m_PreTreatedData->GetExoOuter1(EventId) > OuterMax){
    OuterMax =  m_PreTreatedData->GetExoOuter1(EventId);
    OuterId = 0;
  }
  if(m_PreTreatedData->GetExoOuter2(EventId) > OuterMax){
    OuterMax =  m_PreTreatedData->GetExoOuter2(EventId);
    OuterId = 1;
  }
  if(m_PreTreatedData->GetExoOuter3(EventId) > OuterMax){
    OuterMax =  m_PreTreatedData->GetExoOuter3(EventId);
    OuterId = 2;
  }
  if(m_PreTreatedData->GetExoOuter4(EventId) > OuterMax){
    OuterMax =  m_PreTreatedData->GetExoOuter4(EventId);
    OuterId = 3;
  }
  return OuterId;
}


///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::Clear() {

  E.clear();
  EHG.clear();
  Outer1.clear();
  Outer2.clear();
  Outer3.clear();
  Outer4.clear();
  Flange.clear();
  Crystal.clear();
  TDC.clear();
  TS.clear();

  E_AB.clear();
  Flange_AB.clear();
  Size_AB.clear();
  Crystal_AB.clear();
  Outer_AB.clear();
  Theta.clear();
  Phi.clear();
  TDC_AB.clear();
  TS_AB.clear();
}

///////////////////////////////////////////////////////////////////////////
//	Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " Exogam clover found " << endl;
  
  vector<string> ExoFlange = {"Flange","Radius"};
  vector<string> ExoCloverXYZ = {"Flange",
    "CrystalA_Seg1","CrystalA_Seg2","CrystalA_Seg3","CrystalA_Seg4",
    "CrystalB_Seg1","CrystalB_Seg2","CrystalB_Seg3","CrystalB_Seg4",
    "CrystalC_Seg1","CrystalC_Seg2","CrystalC_Seg3","CrystalC_Seg4",
    "CrystalD_Seg1","CrystalD_Seg2","CrystalD_Seg3","CrystalD_Seg4"};

  for(unsigned int i=0; i<blocks.size(); i++){
    if (blocks[i]->HasTokenList(ExoFlange)) {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        std::cout << std::endl << "//// EXOGAM Clover in radial geometry " << i+1 << endl;
      
        m_NumberOfClovers++;
      unsigned short flange = blocks[i]->GetInt("Flange");
      double radius = blocks[i]->GetDouble("Radius","mm");
      MapFlangeRadius[flange] = radius;
      if(blocks[i]->HasToken("CrystalA")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalA");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,0);
      }
      if(blocks[i]->HasToken("CrystalB")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalB");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,1);
      }
      if(blocks[i]->HasToken("CrystalC")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalC");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,2);
      }
      if(blocks[i]->HasToken("CrystalD")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalD");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,3);
      }
    }
    if (blocks[i]->HasTokenList(ExoCloverXYZ)) {
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "//// EXOGAM Clover in Cartesian geometry " << i+1 << endl;

      m_NumberOfClovers++;
      int flange = blocks[i]->GetInt("Flange");
      
      if(blocks[i]->HasToken("CrystalA")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalA");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,0);
      }
      if(blocks[i]->HasToken("CrystalB")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalB");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,1);
      }
      if(blocks[i]->HasToken("CrystalC")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalC");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,2);
      }
      if(blocks[i]->HasToken("CrystalD")){
        int RawCrystalNumber = blocks[i]->GetInt("CrystalD");
        MapCrystalFlangeCLover[RawCrystalNumber] = std::make_pair(flange,3);
      }
      
      MapFlangeCoordinates[flange][0][0] = ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg1","mm"));
      MapFlangeCoordinates[flange][0][1] = ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg2","mm"));
      MapFlangeCoordinates[flange][0][2] = ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg3","mm"));
      MapFlangeCoordinates[flange][0][3] = ConvertVector(blocks[i]->GetTVector3("CrystalA_Seg4","mm"));
      MapFlangeCoordinates[flange][1][0] = ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg1","mm"));
      MapFlangeCoordinates[flange][1][1] = ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg2","mm"));
      MapFlangeCoordinates[flange][1][2] = ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg3","mm"));
      MapFlangeCoordinates[flange][1][3] = ConvertVector(blocks[i]->GetTVector3("CrystalB_Seg4","mm"));
      MapFlangeCoordinates[flange][2][0] = ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg1","mm"));
      MapFlangeCoordinates[flange][2][1] = ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg2","mm"));
      MapFlangeCoordinates[flange][2][2] = ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg3","mm"));
      MapFlangeCoordinates[flange][2][3] = ConvertVector(blocks[i]->GetTVector3("CrystalC_Seg4","mm"));
      MapFlangeCoordinates[flange][3][0] = ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg1","mm"));
      MapFlangeCoordinates[flange][3][1] = ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg2","mm"));
      MapFlangeCoordinates[flange][3][2] = ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg3","mm"));
      MapFlangeCoordinates[flange][3][3] = ConvertVector(blocks[i]->GetTVector3("CrystalD_Seg4","mm"));

    } 
  }
  if(MapFlangeRadius.size() > 0)
    ExogamGeo = new TExogamStructure(MapFlangeRadius);
  if(MapFlangeCoordinates.size() > 0)
    ExogamGeo = new TExogamCartesian(MapFlangeCoordinates);

  ReadAnalysisConfig();
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;


  // path to photon cross section
  string CSFilename = string(getenv("NPTOOL")) + "/Inputs/PhotonCrossSection/CoherentGe.xcom";
  string LineBuffer;
  // path to file
  string FileName = "./configs/ConfigExogam.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigExogam.dat found: Default parameters loaded for "
      "Analysis "
      << FileName << endl;
    return;
  }


  string DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 12, "ConfigExogam") == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {

      whatToDo = "";
      AnalysisConfigFile >> whatToDo;
      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }

      else if (whatToDo == "EXO_Threshold") {
        //AnalysisConfigFile >> DataBuffer;
        //m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
        //cout << "MAXIMUN STRIP MULTIPLICITY " << m_MaximumStripMultiplicityAllowed << endl;
      }
      else if (whatToDo == "TDC_THRESHOLDS") {
        AnalysisConfigFile >> DataBuffer;
        m_ExoTDC_LowThreshold = stoi(DataBuffer);
        AnalysisConfigFile >> DataBuffer;
        m_ExoTDC_HighThreshold = stoi(DataBuffer);
        cout << "TDC Thresholds " << m_ExoTDC_LowThreshold << " " <<m_ExoTDC_HighThreshold << endl;
      }
      else if (whatToDo == "IGNORE_TDC") {
        AnalysisConfigFile >> DataBuffer;
        IgnoreTDC.push_back(stoi(DataBuffer));
        cout << "TDC Ignored for Crystals : " << DataBuffer << endl;
      }
      else if (whatToDo == "DATA_IS_CAL") {
        AnalysisConfigFile >> DataBuffer;
        DataIsCal = (stoi(DataBuffer) == 1);
        if(DataIsCal)
          cout << "Using Calibrated Data for Exogam" << endl;
      }
      else{
        ReadingStatus = false;
      }
    }
  }
}
void TExogamPhysics::InitSpectra() {
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::FillSpectra() {
}
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }
///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ClearSpectra() {
  // To be done
}
///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TExogamPhysics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}


//////////////////////////////////////////////////////////////////////////
//	Add Parameter to the CalibrationManger
//////////////////////////////////////////////////////////////////////////
void TExogamPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  for (auto it = MapCrystalFlangeCLover.begin(); it != MapCrystalFlangeCLover.end(); it++)
  {  unsigned int i = it->first;
    Cal->AddParameter("EXO", "E" + NPL::itoa(i),
        "EXO_E" + NPL::itoa(i));
    Cal->AddParameter("EXO", "EHG" + NPL::itoa(i),
        "EXO_EHG" + NPL::itoa(i));
    // Cal->AddParameter("EXOGAM", "Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_T",
    // "EXOGAM_Cl" + NPL::itoa(i) + "_Cr" + NPL::itoa(j) + "_T");

    for (int j = 0; j < 4; j++) {
      Cal->AddParameter("EXO", "Outer" + NPL::itoa(i) + "_" + NPL::itoa(j),
          "EXO_Outer" + NPL::itoa(i) + "_" + NPL::itoa(j));
    }
  }
}


//////////////////////////////////////////////////////////////////////////
//	Activated associated Branches and link it to the private member DetectorData address
//	In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
//////////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use pre calibrated NPTOOL data
  else if (DataIsCal) {
    inputChain->SetBranchStatus("Exogam", true);
    inputChain->SetBranchStatus("cEXO_*", true);
    inputChain->SetBranchAddress("Exogam", &m_PreTreatedData); 
  }
  // Option to use the standard npanalysis
  else{
    inputChain->SetBranchStatus("Exogam", true);
    inputChain->SetBranchStatus("fEXO_*", true);
    inputChain->SetBranchAddress("Exogam", &m_EventData);

  }
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::ReadConfigurationTS(){
  TSEvent->ReadConfigurationFile();
}

///////////////////////////////////////////////////////////////////////////
void TExogamPhysics::SetRefTS(std::string TSRef_Name, unsigned long long TSRef){
  RefTS = TSRef;
  RefTS_Name = TSRef_Name; 
}


/////////////////////////////////////////////////////////////////////
//   Activated associated Branches and link it to the private member DetectorPhysics address
//   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
void TExogamPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else{
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchAddress("Exogam", &m_EventPhysics);
  }
}

/////////////////////////////////////////////////////////////////////
//	Create associated branches and associated private member DetectorPhysics address
/////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("Exogam", "TExogamPhysics", &m_EventPhysics);

  // control histograms if needed
  /*
     TList* outputList = RootOutput::getInstance()->GetList();
     controle = new TH1F("controle","histo de controle",20,0,20);
     outputList->Add(controle);
     */
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::SetTreeReader(TTreeReader* TreeReader) {
  TExogamPhysicsReader::r_SetTreeReader(TreeReader);
}

/////////////////////////////// DoCalibration Part //////////////////////////:
void TExogamPhysics::InitializeRootHistogramsCalib() {
  std::cout << "Initialize Exogam Histograms" << std::endl;
  map<int, bool>::iterator it;
  map<int, map<int,bool>>::iterator it2;
  for (it = DoCalibrationE.begin(); it != DoCalibrationE.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsE_F(it->first);
    }
  }
  for (it = DoCalibrationEHG.begin(); it != DoCalibrationEHG.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsEHG_F(it->first);
    }
  }
  //for (it = DoCalibrationT.begin(); it != DoCalibrationT.end(); it++) {
  //  if (it->second) {
  //    InitializeRootHistogramsT_F(it->first);
  //  }
  //}
  for (it2 = DoCalibrationOuter.begin(); it2 != DoCalibrationOuter.end(); it2++) {
    for (it = (it2->second).begin(); it != (it2->second).end(); it++) {
      if (it->second) {
        InitializeRootHistogramsOuter_F(it2->first,it->first);
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::FillHistogramsCalib() {
  if (NPOptionManager::getInstance()->IsReader())
    m_EventData = &(**r_ReaderEventData);

  FillRootHistogramsCalib_F();
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitializeRootHistogramsE_F(unsigned int DetectorNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();

  TString hnameEXOE = Form("EXO_E%d", DetectorNumber);
  TString htitleEXOE = Form("EXO_E%d", DetectorNumber);
  (*TH1Map)["Exogam"][hnameEXOE] = new TH1F(hnameEXOE, htitleEXOE, 65536, 0, 65536);
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitializeRootHistogramsOuter_F(unsigned int DetectorNumber, unsigned int OuterNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();

  TString hnameEXOOuter = Form("EXO_Outer%d_%d", DetectorNumber, OuterNumber);
  TString htitleEXOOuter = Form("EXO_Outer%d_%d", DetectorNumber, OuterNumber);
  (*TH1Map)["Exogam"][hnameEXOOuter] = new TH1F(hnameEXOOuter, htitleEXOOuter, 65536, 0, 65536);
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::InitializeRootHistogramsEHG_F(unsigned int DetectorNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();

  TString hnameEXOEHG = Form("EXO_EHG%d", DetectorNumber);
  TString htitleEXOEHG = Form("EXO_EHG%d", DetectorNumber);
  (*TH1Map)["Exogam"][hnameEXOEHG] = new TH1F(hnameEXOEHG, htitleEXOEHG, 65536, 0, 65536);

}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::FillRootHistogramsCalib_F(){
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  TString hname;

  for (UShort_t i = 0; i < m_EventData->GetExoMult(); i++) {
    unsigned int DetectorNbr = m_EventData->GetExoCrystal(i);

    if(DoCalibrationE[DetectorNbr] && m_EventData->GetExoE(i) >0){
      hname = Form("EXO_E%d", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoE(i));
    }
    if(DoCalibrationEHG[DetectorNbr] && m_EventData->GetExoEHG(i) >0){
      hname = Form("EXO_EHG%d", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoEHG(i));
    }
    if(DoCalibrationOuter[DetectorNbr][0] && m_EventData->GetExoOuter1(i) >0){
      hname = Form("EXO_Outer%d_0", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter1(i));
    }
    if(DoCalibrationOuter[DetectorNbr][1] && m_EventData->GetExoOuter2(i) >0){
      hname = Form("EXO_Outer%d_1", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter2(i));
    }
    if(DoCalibrationOuter[DetectorNbr][2] && m_EventData->GetExoOuter3(i) >0){
      hname = Form("EXO_Outer%d_2", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter3(i));
    }
    if(DoCalibrationOuter[DetectorNbr][3] && m_EventData->GetExoOuter4(i) >0){
      hname = Form("EXO_Outer%d_3", DetectorNbr);
      (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoOuter4(i));
    }
  }
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::DoCalibration() {
  std::cout << "Do Calibration Exogam" << std::endl;
  DefineCalibrationSource(Source_name);
  map<int, bool>::iterator it;
  map<int, map<int,bool>>::iterator it2;

  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  std::string make_folder = "mkdir " + Path + OutputName;

  MakeFolder(make_folder);

  ofstream* calib_file = new ofstream;
  ofstream* dispersion_file = new ofstream;

  if(!DoCalibrationE.empty()){
    MakeECalibFolders(make_folder);
    CreateCalibrationEFiles(calib_file, dispersion_file);
  }
  for (it = DoCalibrationE.begin(); it != DoCalibrationE.end(); it++) {
    if (it->second) {
      DoCalibrationE_F(it->first,"E", calib_file, dispersion_file, Threshold_E_Cal);
    }
  }
  calib_file->close();
  dispersion_file->close();

  if(!DoCalibrationEHG.empty()){
    MakeEHGCalibFolders(make_folder);
    CreateCalibrationEHGFiles(calib_file, dispersion_file);
  }
  for (it = DoCalibrationEHG.begin(); it != DoCalibrationEHG.end(); it++) {
    if (it->second) {
      DoCalibrationE_F(it->first,"EHG", calib_file, dispersion_file, Threshold_EHG_Cal);
    }
  }
  calib_file->close();
  dispersion_file->close();

  if(!DoCalibrationOuter.empty()){
    MakeOuterCalibFolders(make_folder);
    CreateCalibrationOuterFiles(calib_file, dispersion_file);
  }
  for (it2 = DoCalibrationOuter.begin(); it2 != DoCalibrationOuter.end(); it2++) {
    for (it = (it2->second).begin(); it != (it2->second).end(); it++) {
      if (it->second) {
        DoCalibrationE_F(it->first,Form("Outer%d_",it2->first), calib_file, dispersion_file, Threshold_Outers_Cal);
      }
    }
  }
  calib_file->close();
  dispersion_file->close();
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::MakeFolder(std::string make_folder) {
  int sys = system(make_folder.c_str());
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::MakeECalibFolders(std::string make_folder) {
  int sys =system((make_folder+"/Exogam_E").c_str()); 
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::MakeEHGCalibFolders(std::string make_folder) {
  int sys =system((make_folder+"/Exogam_EHG").c_str());
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::MakeOuterCalibFolders(std::string make_folder) {
  int sys =system((make_folder+"/Exogam_Outer").c_str());
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::DoCalibrationE_F(unsigned int DetectorNumber,std::string CalibType, ofstream* calib_file, ofstream* dispersion_file, unsigned int Threshold) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();

#if CUBIX
  CubixEnergyCal->Reset();
  std::string hnameEXOE = Form("EXO_%s%d",CalibType.c_str(), DetectorNumber);
  std::string htitleEXOE = Form("EXO_%s%d",CalibType.c_str(), DetectorNumber);

  auto hist = ((*TH1Map)["Exogam"][hnameEXOE]);

  CubixEnergyCal->SetDataFromHistTH1(hist,0);

  for (auto ie : Source_E)
    CubixEnergyCal->AddPeak(ie);

  CubixEnergyCal->SetGain(1.);
  CubixEnergyCal->SetVerbosityLevel(1);

  CubixEnergyCal->SetFitPlynomialOrder(FitPolOrder);
  CubixEnergyCal->SetNoOffset(false);
  CubixEnergyCal->UseLeftTail(true);
  CubixEnergyCal->UseRightTail(true);

  CubixEnergyCal->UseFirstDerivativeSearch();

  CubixEnergyCal->SetGlobalChannelLimits(hist->GetXaxis()->GetBinLowEdge(1)+Threshold,hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetNbins()));      // limit the search to this range in channels
  CubixEnergyCal->SetGlobalPeaksLimits(15,5);   // default fwhm and minmum amplitude for the peaksearch [15 5]

  CubixEnergyCal->StartCalib();
  vector < Fitted > FitResults = CubixEnergyCal->GetFitResults();



  std:: cout << calib_file << " " << (*calib_file).is_open() << std::endl;
  std:: cout << hnameEXOE << " ";
  (*calib_file) << hnameEXOE << " ";
  if(FitResults.size() > 1)
  {
    for(unsigned int i = 0; i <= FitPolOrder; i++){
      (*calib_file) << scientific << setprecision(6) << setw(14) << CubixEnergyCal->fCalibFunction->GetParameter(i) << " ";
      std::cout << scientific << setprecision(6) << setw(14) << CubixEnergyCal->fCalibFunction->GetParameter(i) << " ";
    }
  }
  else
  {
    for(unsigned int i = 0; i <= FitPolOrder; i++){
      (*calib_file) << scientific << setprecision(6) << setw(14) << 0. << " ";
      std::cout << scientific << setprecision(6) << setw(14) << 0. << " ";
    }
  }
  (*calib_file) << "\n";
  std::cout << "\n";


  if(FitResults.size()>1 && CubixEnergyCal->fCalibFunction) {
    auto c = new TCanvas;
    c->SetName("CalibrationResults");
    c->SetTitle("Calibration Results");
    c->Divide(1,2,0.0001,0.0001);
    c->cd(1);
    CubixEnergyCal->fCalibGraph->Draw("ap");
    CubixEnergyCal->fCalibFunction->Draw("same");
    c->cd(2);
    CubixEnergyCal->fResidueGraph->Draw("ape");
    c->Update();
    c->Modified();
  }
  if(FitResults.size() > 1)
  {
    (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())] = (TGraphErrors*)(CubixEnergyCal->fCalibGraph->Clone());
    (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->GetYaxis()->SetTitle("Energy (MeV)");
    (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->SetTitle(Form("Calibration_Graph_%s",hnameEXOE.c_str()));

    (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())] = (TGraphErrors*)(CubixEnergyCal->fResidueGraph->Clone());
    (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->GetXaxis()->SetTitle("Energy (MeV)");
    (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->GetYaxis()->SetTitle("Residue (MeV)");
    (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->SetTitle(Form("Residue_Graph_%s",hnameEXOE.c_str()));
  }
#else
  std::cout << "Exogam calibration currently not supported without CUBIX. Download CUBIX and set -DCUBIX=1 to use EXOGAM calibration\n";
  exit(1);

#endif

}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::DefineCalibrationSource(std::string source) {
  // 239Pu
  if(source == "60Co"){
    Source_isotope.push_back("$^{60}$Co");
    Source_E.push_back(1.17322);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(99.85);

    Source_isotope.push_back("$^{60}$Co");
    Source_E.push_back(1.33249);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(99.98);
  }
  else if(source == "152Eu"){
    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.121782);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(28.58);
    Source_branching_ratio_err.push_back(0.16);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.344279);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(26.59);
    Source_branching_ratio_err.push_back(0.20);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.40801);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(20.87);
    Source_branching_ratio_err.push_back(0.09);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.964079);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(14.51);
    Source_branching_ratio_err.push_back(0.07);

    // FIXME Fit does not work properly for this peak (no idea why)
    // Source_isotope.push_back("$^{152}$Eu");
    // Source_E.push_back(1.11207);
    // Source_Sig.push_back(0.0001);
    // Source_branching_ratio.push_back(13.67);
    // Source_branching_ratio_err.push_back(0.08);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.778904);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(12.93);
    Source_branching_ratio_err.push_back(0.08);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.08587);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(10.11);
    Source_branching_ratio_err.push_back(0.05);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.244698);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(7.58);
    Source_branching_ratio_err.push_back(0.04);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.867378);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(4.23);
    Source_branching_ratio_err.push_back(0.03);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.443965);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(2.827);
    Source_branching_ratio_err.push_back(0.014);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(0.411116);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(2.237);
    Source_branching_ratio_err.push_back(0.013);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.08974);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(1.734);
    Source_branching_ratio_err.push_back(0.011);

    // FIXME Fit does not work properly for this peak (no idea why)
    // Source_isotope.push_back("$^{152}$Eu");
    // Source_E.push_back(1.29914);
    // Source_Sig.push_back(0.0001);
    // Source_branching_ratio.push_back(1.633);
    // Source_branching_ratio_err.push_back(0.011);

    Source_isotope.push_back("$^{152}$Eu");
    Source_E.push_back(1.21295);
    Source_Sig.push_back(0.0001);
    Source_branching_ratio.push_back(1.415);
    Source_branching_ratio_err.push_back(0.008);

  }
  else{
    std::cout << "Please enter a valid source for gamma ray calibration\nCurrently supported sources are 60Co and 152Eu\n";
    exit(1);
  }
  std::cout << "///////////////////////////////// " << Source_E.size() << " " << Source_branching_ratio.size() << std::endl;
}


/////////////////////////////////////////////////////////////////////
// FIXME Probably could be done better, currently a but inelegant
/////////////////////////////////////////////////////////////////////
void TExogamPhysics::CreateCalibrationEFiles(ofstream* calib_file,
    ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Exogam_E";
  (*calib_file).open(((string)(Path + OutputName + "/Exogam_E/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/Exogam_E/" +Filename + ".dispersion")).c_str());
}

void TExogamPhysics::CreateCalibrationEHGFiles(ofstream* calib_file,
    ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Exogam_EHG";
  (*calib_file).open(((string)(Path + OutputName + "/Exogam_EHG/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/Exogam_EHG/" +Filename + ".dispersion")).c_str());
}

void TExogamPhysics::CreateCalibrationOuterFiles(ofstream* calib_file,
    ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Exogam_Outer";
  (*calib_file).open(((string)(Path + OutputName + "/Exogam_Outer/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/Exogam_Outer/" +Filename + ".dispersion")).c_str());
}

void TExogamPhysics::ReadDoCalibration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");

  vector<string> calibs = {"Threshold_E","Threshold_EHG","Threshold_Outers","FirstCr","LastCr","FirstOuter","LastOuter","FitOrder","Source"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(calibs)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Exogam Calibration" << endl;
      unsigned int FirstCr = blocks[i]->GetInt("FirstCr");
      unsigned int LastCr = blocks[i]->GetInt("LastCr");
      unsigned int FirstOuter = blocks[i]->GetInt("FirstOuter");
      unsigned int LastOuter = blocks[i]->GetInt("LastOuter");
      FitPolOrder = blocks[i]->GetInt("FitOrder");
      Source_name = blocks[i]->GetString("Source");
      Threshold_E_Cal = blocks[i]->GetInt("Threshold_E");
      Threshold_EHG_Cal = blocks[i]->GetInt("Threshold_EHG");
      Threshold_Outers_Cal = blocks[i]->GetInt("Threshold_Outers");
      for(unsigned int k = FirstCr; k <= LastCr; k++){
        DoCalibrationE[k] = true;
        DoCalibrationEHG[k] = true; 
        for(unsigned int p = FirstOuter; p <= LastOuter; p++){
          DoCalibrationOuter[k][p] = true; 
        }
      }
    }
    else {
      cout << "ERROR: Missing token for Exogam DoCalibration blocks, check your "
        "input "
        "file"
        << endl;
      exit(1);
    }
  }
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::WriteHistogramsCalib() {
  std::cout << "Writing Exogam Histograms\n";
  WriteHistogramsE();
  RootHistogramsCalib::getInstance()->GetFile()->Close();
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::WriteHistogramsEfficiency() {
  std::cout << "Writing Exogam Efficiency\n";
  auto File = RootHistogramsCalib::getInstance()->GetFile();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  auto TCanvasMap = RootHistogramsCalib::getInstance()->GetTCanvasMap();
  
  std::string Path = NPOptionManager::getInstance()->GetEfficiencyOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();

  map<int, bool>::iterator it;
  std::string hnameEXOE;

  if (!File->GetDirectory("Exogam"))
    File->mkdir("Exogam");
  File->cd("Exogam");

  for (it = DoEfficiency.begin(); it != DoEfficiency.end(); it++) {
    if (it->second) {

      if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it->first)))
        gDirectory->mkdir(Form("EXO_Cr%d", it->first));
      gDirectory->cd(Form("EXO_Cr%d", it->first));

      hnameEXOE = Form("EXO_E%d", it->first);
      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      (*TCanvasMap)["Exogam"][hnameEXOE]->SaveAs((Path+OutputName+"/"+ hnameEXOE+".png").c_str());
      hnameEXOE = Form("EXO_TS%d", it->first);
      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      hnameEXOE = Form("EXO_Eff%d", it->first);
      (*TGraphMap)["Exogam"][hnameEXOE]->Write();
      (*TCanvasMap)["Exogam"][hnameEXOE]->SaveAs((Path+OutputName+"/"+ hnameEXOE+".png").c_str());
      //if((*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_E%d",it->first)]!= nullptr)
      //  (*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_E%d",it->first)]->Write();
      //if((*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_E%d",it->first)]!= nullptr)
      //  (*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_E%d",it->first)]->Write();
    }
    File->cd("Exogam");
  }
  hnameEXOE = "EXO_E_all";
  (*TH1Map)["Exogam"][hnameEXOE]->Write();
  (*TCanvasMap)["Exogam"][hnameEXOE]->SaveAs((Path+OutputName+"/"+ hnameEXOE+".png").c_str());
  hnameEXOE = "EXO_TS_all";
  (*TH1Map)["Exogam"][hnameEXOE]->Write();
  hnameEXOE = "EXO_Eff_all";
  (*TGraphMap)["Exogam"][hnameEXOE]->Write();
  (*TCanvasMap)["Exogam"][hnameEXOE]->SaveAs((Path+OutputName+"/"+ hnameEXOE+".png").c_str());
  for(auto PeakE: Source_E)
    (*TGraphMap)["Exogam"]["Exo_Eff_crystal "+NPL::itoa(int(PeakE*1000.))+" keV"]->Write();

  
  
  RootHistogramsCalib::getInstance()->GetFile()->Close();
}

/////////////////////////////////////////////////////////////////////
void TExogamPhysics::WriteHistogramsE() {
  auto File = RootHistogramsCalib::getInstance()->GetFile();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();

  map<int, bool>::iterator it;
  map<int, map<int,bool>>::iterator it2;
  std::string hnameEXOE;

  if (!File->GetDirectory("Exogam"))
    File->mkdir("Exogam");
  File->cd("Exogam");

  for (it = DoCalibrationE.begin(); it != DoCalibrationE.end(); it++) {
    if (it->second) {
      hnameEXOE = Form("EXO_E%d", it->first);

      if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it->first)))
        gDirectory->mkdir(Form("EXO_Cr%d", it->first));
      gDirectory->cd(Form("EXO_Cr%d", it->first));

      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      if((*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_E%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_E%d",it->first)]->Write();
      if((*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_E%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_E%d",it->first)]->Write();
    }
    File->cd("Exogam");
  }

  for (it = DoCalibrationEHG.begin(); it != DoCalibrationEHG.end(); it++) {
    if (it->second) {
      hnameEXOE = Form("EXO_EHG%d", it->first);

      if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it->first)))
        gDirectory->mkdir(Form("EXO_Cr%d", it->first));
      gDirectory->cd(Form("EXO_Cr%d", it->first));

      (*TH1Map)["Exogam"][hnameEXOE]->Write();
      if((*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_EHG%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Calib_Graph_EXO_EHG%d",it->first)]->Write();
      if((*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_EHG%d",it->first)]!= nullptr)
        (*TGraphMap)["Exogam"][Form("Residue_Graph_EXO_EHG%d",it->first)]->Write();
    }
    File->cd("Exogam");
  }
  for (it2 = DoCalibrationOuter.begin(); it2 != DoCalibrationOuter.end(); it2++) {
    for (it = (it2->second).begin(); it != (it2->second).end(); it++) {
      if (it->second) {
        hnameEXOE = Form("EXO_Outer%d_%d",it2->first,it->first);

        if (!gDirectory->GetDirectory(Form("EXO_Cr%d", it2->first)))
          gDirectory->mkdir(Form("EXO_Cr%d", it2->first));
        gDirectory->cd(Form("EXO_Cr%d", it2->first));

        if (!gDirectory->GetDirectory("Outers"))
          gDirectory->mkdir("Outer");
        gDirectory->cd("Outer");

        (*TH1Map)["Exogam"][hnameEXOE]->Write();
        if((*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]!= nullptr)
          (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->Write();
        if((*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]!= nullptr)
          (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->Write();
      }
      File->cd("Exogam");
    }
  }
}

void TExogamPhysics::InitializeRootHistogramsEfficiency() {
  DefineCalibrationSource(Source_name);
  
  std::cout << "Initialize Exogam Efficiency Histograms" << std::endl;
  map<int, bool>::iterator it;
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();

  // unsigned long long FirstTS = FindFirstExoTS();
  // unsigned long long LastTS = FindLastExoTS();
  // std::cout << "Les TS: " << FirstTS << " " << LastTS << std::endl;
  //FIXME First TS and Last TS work byut idk how to handle runs with bad TS indexin//FIXME First TS and Last TS work byut idk how to handle runs with bad TS indexingg
  for (it = DoEfficiency.begin(); it != DoEfficiency.end(); it++) {
    if (it->second) {
      TString hnameEXOE = Form("EXO_E%d", it->first);
      TString htitleEXOE = Form("EXO_E%d", it->first);
      (*TH1Map)["Exogam"][hnameEXOE] = new TH1F(hnameEXOE, htitleEXOE, 5000, 0, 5000);
      hnameEXOE = Form("EXO_TS%d", it->first);
      htitleEXOE = Form("EXO_TS%d", it->first);
      (*TH1Map)["Exogam"][hnameEXOE] = new TH1F(hnameEXOE, htitleEXOE,10000,0,pow(2,48)); 
    }
  }
  (*TH1Map)["Exogam"]["EXO_E_all"] = new TH1F("EXO_E_all", "EXO_E_all", 5000, 0, 5000);
  (*TH1Map)["Exogam"]["EXO_TS_all"] = new TH1F("EXO_TS_all", "EXO_TS_all", 10000,0,pow(2,48));
  for(auto PeakE: Source_E){
    (*TGraphMap)["Exogam"]["Exo_Eff_crystal "+NPL::itoa(int(PeakE*1000.))+" keV"] = new TGraphErrors;
    std::cout << "Exo_Eff_crystal "+NPL::itoa(int(PeakE*1000.))+" keV" << std::endl;
    (*TGraphMap)["Exogam"]["Exo_Eff_crystal "+NPL::itoa(int(PeakE*1000.))+" keV"]->SetName(("EfficiencyByCrystal"+NPL::itoa(int(PeakE*1000.))+"keV").c_str());
    (*TGraphMap)["Exogam"]["Exo_Eff_crystal "+NPL::itoa(int(PeakE*1000.))+" keV"]->SetTitle(("Efficiency By Crystal "+NPL::itoa(int(PeakE*1000.))+" keV").c_str());
  }
}


unsigned long long TExogamPhysics::FindFirstExoTS(){
  TTreeReader* TreeReader = RootInput::getInstance()->GetTreeReader();
  unsigned long long FirstTSValue = 0;
  bool FirstTSFound = false;
  TreeReader->Restart();
  while(TreeReader->Next() && !FirstTSFound){
    ClaimReaderData();
    for(auto TSvalue: m_EventData->fExo_TS){
      if(TSvalue > 0){
        FirstTSValue = TSvalue;
        FirstTSFound = true; 
      }
    }
  }
  TreeReader->Restart();
  return FirstTSValue;
}

unsigned long long TExogamPhysics::FindLastExoTS(){
  TTreeReader* TreeReader = RootInput::getInstance()->GetTreeReader();
  unsigned long long LastTSValue = 0;
  bool LastTSFound = false;
  TreeReader->Restart();
  for(unsigned long long i = TreeReader->GetTree()->GetEntries() - 1; i >= 0 && !LastTSFound; i--){
    TreeReader->SetEntry(i);
    ClaimReaderData();
    for(auto TSvalue: m_EventData->fExo_TS){
      if(TSvalue > 0){
        LastTSValue = TSvalue;
        std::cout << LastTSValue << " " << m_EventData->fExo_Crystal[0] << " " << m_EventData->fExo_Crystal.size() <<  " " << i <<  std::endl;
        LastTSFound = true; 
      }
    }
  }
  TreeReader->Restart();
  return LastTSValue;
}

void TExogamPhysics::Efficiency() {
  std::cout << "Evaluate Efficiency Exogam" << std::endl;
  
  map<int, bool>::iterator it;
  
  std::string Path = NPOptionManager::getInstance()->GetEfficiencyOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  std::string make_folder = "mkdir " + Path + OutputName;

  MakeFolder(make_folder);

  ofstream* efficiency_file = new ofstream;
  if(!DoEfficiency.empty()){
    CreateEfficiencyFile(efficiency_file);
  }
  for (it = DoEfficiency.begin(); it != DoEfficiency.end(); it++) {
    if (it->second) {
      EvaluateEfficiency_F(NPL::itoa((int)it->first), efficiency_file, Source_activity);
    }
  }
  EvaluateEfficiency_F("_all", efficiency_file, Source_activity);
  efficiency_file->close();
}

void TExogamPhysics::CreateEfficiencyFile(ofstream* efficiency_file) {
  std::string Path = NPOptionManager::getInstance()->GetEfficiencyOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Efficiency_EXOGAM";
  (*efficiency_file).open(((string)(Path + OutputName + Filename + ".ef")).c_str());
}

void TExogamPhysics::ReadEfficiency(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Exogam");

  vector<string> eff = {"FirstCr","LastCr","Source","SourceActivity","minmax","Time"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(eff)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Exogam Efficiency" << endl;
      unsigned int FirstCr = blocks[i]->GetInt("FirstCr");
      unsigned int LastCr = blocks[i]->GetInt("LastCr");
      Source_name = blocks[i]->GetString("Source");
      Source_activity = blocks[i]->GetInt("SourceActivity");
      minmax = blocks[i]->GetInt("minmax");
      Time = blocks[i]->GetInt("Time");
      for(unsigned int k = FirstCr; k <= LastCr; k++){
        DoEfficiency[k] = true;
      }
    }
    else {
      cout << "ERROR: Missing token for Exogam Efficiency blocks, check your "
        "input "
        "file"
        << endl;
      exit(1);
    }
  }
}

void TExogamPhysics::FillHistogramsEfficiency() {
  Clear();
  if (NPOptionManager::getInstance()->IsReader())
    m_EventData = &(**r_ReaderEventData);

  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  TString hname;
  
  std::map<unsigned int,std::vector<unsigned int>> HitsID;

  for(unsigned int i = 0; i < m_EventData->GetExoMult(); i++){
      // Doing flange and crystal matching
      if (m_EventData->GetExoE(i) > m_EXO_E_RAW_Threshold){
        E.push_back(fEXO_E(m_EventData, i));
        flange_nbr = MapCrystalFlangeCLover[m_EventData->GetExoCrystal(i)].first;
        crystal_nbr = MapCrystalFlangeCLover[m_EventData->GetExoCrystal(i)].second;
        Flange.push_back(flange_nbr);
        Crystal.push_back(crystal_nbr);
        TS.push_back(m_EventData->GetExoTS(i));
        HitsID[flange_nbr].push_back(i);
      }
    }
  for(auto it = HitsID.begin(); it != HitsID.end(); it++){
    double E_AddBack = 0;
    double E_Max = 0;
    unsigned int Id_Max = 0;
    for(auto itvec = (*it).second.begin(); itvec !=(*it).second.end(); itvec++){
      double E_Crystal = fEXO_E(m_EventData, *itvec);
      // std::cout << E_Crystal << " " << m_EventData->GetExoE(*itvec) << std::endl;
      E_AddBack+= E_Crystal;
      if(E_Max < E_Crystal){
        E_Max = E_Crystal;
        Id_Max = *itvec;
        
      }
    }
    hname = Form("EXO_E%d", m_EventData->GetExoCrystal(Id_Max));
    (*TH1Map)["Exogam"][hname]->Fill(1000.*E_AddBack);
    hname = Form("EXO_TS%d",m_EventData->GetExoCrystal(Id_Max));
    (*TH1Map)["Exogam"][hname]->Fill(m_EventData->GetExoTS(Id_Max));
    (*TH1Map)["Exogam"]["EXO_E_all"]->Fill(1000.*E_AddBack);
    (*TH1Map)["Exogam"]["EXO_TS_all"]->Fill(m_EventData->GetExoTS(Id_Max));
    // if(m_EventData->GetExoTS(Id_Max) > 15657103568298)
    // std::cout << E_AddBack << " " << m_EventData->GetExoTS(Id_Max) << " " <<m_EventData->GetExoCrystal(Id_Max)  << std::endl;
  }
}

void TExogamPhysics::EvaluateEfficiency_F(std::string DetectorID, ofstream* efficiency_file, double Source_activity) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TCanvasMap = RootHistogramsCalib::getInstance()->GetTCanvasMap();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  
  
  TString hname = "EXO_E"+DetectorID;
  
  auto c1 = new TCanvas("c1","c1",800,600);
  c1->SetName("Peaks Fitted "+hname);
  c1->SetTitle("Peaks Fitted "+hname);
  FindPeaks((*TH1Map)["Exogam"][hname]);
  (*TH1Map)["Exogam"][hname]->Draw();
  (*TH1Map)["Exogam"][hname]->GetXaxis()->SetRangeUser(0.,*std::max_element(Source_E.begin(),Source_E.end())*1.1*1000);
  c1->Modified();
  c1->Update();
  (*TCanvasMap)["Exogam"][hname] = c1;
  std::vector<std::pair<double, double>> EfficiencyPoints;
  std::vector<std::pair<double, double>> EfficiencyPointsErr;


  // Get number of counts for each energy
  for(int j = 0; j < Energies.size(); j++){
    FitResult= nullptr;
    double Area[2];
    FitFunction((*TH1Map)["Exogam"][hname],j,Area);
    std::cout << "TEST HUGO " << Area[0] << " " << Area[1] << std::endl;
    double FWHM = fFitFunction->GetParameter(6);
    double FWHMerr = fFitFunction->GetParError(6);
    double mean = fFitFunction->GetParameter(5);
    double meanerr = fFitFunction->GetParError(5);
    EfficiencyPoints.push_back(std::make_pair(mean,Area[0]));
    EfficiencyPointsErr.push_back(std::make_pair(meanerr,Area[1]));
  }
  // Get TS difference
  hname = "EXO_TS"+DetectorID;
  // double TimeDifference = ((*TH1Map)["Exogam"][hname]->GetXaxis()->GetXmax() - (*TH1Map)["Exogam"][hname]->GetXaxis()->GetXmin())/(pow(10,8));
  double TimeDifference = Time;


  auto g = new TGraphErrors;
  std:: cout << efficiency_file << " " << (*efficiency_file).is_open() << std::endl;
  hname = "Detector "+DetectorID;
  std:: cout << hname << " ";
  g->SetName("EfficiencyResults "+hname);
  g->SetTitle("Efficiency Results "+hname);
  (*efficiency_file) << hname << " ";
    for(unsigned int i = 0; i < EfficiencyPoints.size(); i++){
  
      double Efficiency = EfficiencyPoints[i].second/(Source_activity*Source_branching_ratio[SourceID[i]]/100.*TimeDifference);
      g->SetPoint(i,EfficiencyPoints[i].first, Efficiency);
      double EfficiencyErr = Efficiency*sqrt(pow(1e-8/TimeDifference,2) + pow(Source_branching_ratio_err[SourceID[i]]/Source_branching_ratio[SourceID[i]],2) + pow(EfficiencyPointsErr[i].second/EfficiencyPoints[i].second,2));      
      g->SetPointError(i,EfficiencyPointsErr[i].first, EfficiencyErr);
      if(DetectorID.compare("_all")!=0)
      { 
        auto crystalgraph = (*TGraphMap)["Exogam"]["Exo_Eff_crystal "+NPL::itoa(int(Source_E[SourceID[i]]*1000.))+" keV"];
        std::cout << "Exo_Eff_crystal "+NPL::itoa(int(Source_E[SourceID[i]]*1000.))+" keV" << std::endl;
        unsigned int N = crystalgraph->GetN();
        crystalgraph->SetPoint(N, std::atoi(DetectorID.c_str()), Efficiency);
        crystalgraph->SetPointError(N, 0., EfficiencyErr);
      }
      (*efficiency_file) << scientific << setprecision(6) << setw(14) << EfficiencyPoints[i].first << " " 
      << EfficiencyPointsErr[i].first << " " << Efficiency << " " << EfficiencyErr;
  }
  (*efficiency_file) << "\n";
  hname = "EXO_Eff"+DetectorID;
  (*TGraphMap)["Exogam"][hname] = g;

  auto c2 = new TCanvas("c2","c2",800,600);
  c2->SetName("EfficiencyResults "+hname);
  c2->SetTitle("Efficiency Results "+hname);
  (*TGraphMap)["Exogam"][hname]->Draw("ap");
  // gDirectory->Add((*TGraphMap)["Exogam"][hname]);
  c2->Modified();
  c2->Update();
  (*TCanvasMap)["Exogam"][hname] = c2;

  //if(FitResults.size() > 1)
  //{
  //  (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())] = (TGraphErrors*)(CubixEnergyCal->fCalibGraph->Clone());
  //  (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->GetYaxis()->SetTitle("Energy (MeV)");
  //  (*TGraphMap)["Exogam"][Form("Calib_Graph_%s",hnameEXOE.c_str())]->SetTitle(Form("Calibration_Graph_%s",hnameEXOE.c_str()));

  //  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())] = (TGraphErrors*)(CubixEnergyCal->fResidueGraph->Clone());
  //  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->GetXaxis()->SetTitle("Energy (MeV)");
  //  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->GetYaxis()->SetTitle("Residue (MeV)");
  //  (*TGraphMap)["Exogam"][Form("Residue_Graph_%s",hnameEXOE.c_str())]->SetTitle(Form("Residue_Graph_%s",hnameEXOE.c_str()));
  //}

}

int TExogamPhysics::FindMatchingEnergy(double Energy){
  for(unsigned int i = 0; i < Source_E.size(); i++){
    if(abs(Source_E[i]*1000 - Energy) < minmax)
      return i;
  }
  return -1;
}

void TExogamPhysics::FindPeaks(TH1* Hist){
  Energies.clear();
  MinMax.clear();
  SourceID.clear();
  TSpectrum *tspec = new TSpectrum;
  tspec->Search(Hist,sigma,"goff",threshold);

  int NPeaks = tspec->GetNPeaks();

  for(int i=0 ; i<NPeaks ; i++) {
    double Energy = tspec->GetPositionX()[i];
    double Value = tspec->GetPositionY()[i];
    double MaxGlob = Hist->GetMaximum();
    // Next line finds is a peak in the source is recognized
    std::cout << "Energy 1 " << Energy << std::endl;
    int SourceIndex = FindMatchingEnergy(Energy);
    std::cout << Energy << " " << SourceIndex << std::endl;
      if(SourceIndex > -1){
        Energies.push_back(Energy);
        MinMax.push_back(std::make_pair(Energy-minmax,Energy+minmax));
        SourceID.push_back(SourceIndex);
      }
  }
  
  Hist->GetYaxis()->SetRangeUser(0,Hist->GetMaximum()*1.1);

  std::vector<int> Indexes(Energies.size());
  for(unsigned int i = 0; i < Indexes.size();i++){
    Indexes[i] = i;
  }

  std::sort(Indexes.begin(), Indexes.end(), [this](int i1, int i2) {
        return this->Energies[i1] < this->Energies[i2];});
  std::vector<double> Energies_sorted(Energies.size());
  std::vector<int> SourceID_sorted(SourceID.size());
  std::vector<std::pair<double,double>> MinMax_sorted(MinMax.size());

  for(unsigned int i = 0; i < Indexes.size();i++){
    Energies_sorted[i] = Energies[Indexes[i]];
    SourceID_sorted[i] = SourceID[Indexes[i]];
    MinMax_sorted[i] = MinMax[Indexes[i]];
  }

  Energies = Energies_sorted;
  SourceID = SourceID_sorted;
  MinMax = MinMax_sorted;
  
  for(auto v: Indexes)
    std::cout << v << std::endl;
  
  for(auto v: Energies)
    std::cout << v << std::endl;
  for(auto v: MinMax)
    std::cout << v.first << " " << v.second << std::endl;
  std::cout << std::endl;
}

void TExogamPhysics::InitFitParameters(){
  Minimizer = "Minuit2";
  Algorithm = "Migrad";
  Tolerance = 0.1;
  PrintLevel = 0;


  // Fit Parameters
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer(Minimizer.c_str(),Algorithm.c_str());
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(Tolerance);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(PrintLevel);
  DefFWHM = 2;
  DefFWHM_min = 1;
  DefFWHM_max = 5;

  UseLT = false;
  LeftTailVal = 5;
  LeftTailValMin = 2;
  LeftTailValMax = 10;

  UseRT = false;
  RightTailVal = 2;
  RightTailValMin = 0.1;
  RightTailValMax = 5;

  UseStep = false;
  StepVal = 0.01;
  StepValMin = -1.;
  StepValMax = 1.;

  BackgroundType = "exp";
}

void TExogamPhysics::FitFunction(TH1* Hist, int PeakID, double* Area){
    // Somehow a copy of cubix fit function, but removing all the parts of the interface 
  
    InitFitParameters();

    if(Energies.size() != MinMax.size())
      {
        std::cout << "Error: Energies and MinMax should have the same size !" << std::endl;
        return;
      }
    // Int_t NPars = 4+6*Energies.size();
      int NPars = 10;

      delete fFitFunction;
      fFitFunction = new TF1("MyFit", this, &TExogamPhysics::DoubleTailedStepedGaussian, MinMax[PeakID].first, MinMax[PeakID].second, NPars, "TExogamPhysicsFit", "DoubleTailedStepedGaussian");

      // Changing parameter 0 so that its the peakID (nb of peaks always 1)
      fFitFunction->SetParName(0, "PeakID");
      fFitFunction->SetParName(1, "BkgConst");
      fFitFunction->SetParName(2, "BkgSlope");
      fFitFunction->SetParName(3, "BkgExp");

      fFitFunction->SetParName(4+0, Form("Height"));
      fFitFunction->SetParName(4+1, Form("Position"));
      fFitFunction->SetParName(4+2, Form("FWHM"));
      fFitFunction->SetParName(4+3, Form("LeftTail"));
      fFitFunction->SetParName(4+4, Form("RightTail"));
      fFitFunction->SetParName(4+5, Form("AmplitudeStep"));

      fFitFunction->SetNpx(1000);
      fFitFunction->SetLineColor(kRed);

      fFitFunction->FixParameter(0, PeakID); // 1 peak

      //Double_t x,y;
      //x = fPad->GetUxmin();
      //y = fPad->GetUxmax();

      //Calc Bckd
      fFitFunction->SetParameter(1, Hist->GetBinContent(Hist->FindBin(MinMax[PeakID].first)));
      fFitFunction->SetParLimits(1, 0.,Hist->GetMaximum());

      fFitFunction->SetParameter(2, 0);
      fFitFunction->SetParLimits(2, -50., 50.);
      fFitFunction->SetParameter(3, 0.);
      fFitFunction->SetParLimits(3, -50., 50.);

      if(BackgroundType == "pol1")
          fFitFunction->FixParameter(3,0);
      else if(BackgroundType == "exp")
          fFitFunction->FixParameter(2,0);

          //Height
          fFitFunction->SetParameter(4+0, Hist->GetBinContent(Hist->FindBin(Energies[PeakID])) - (Hist->GetBinContent(Hist->FindBin(MinMax[PeakID].first))+Hist->GetBinContent(Hist->FindBin(MinMax[PeakID].first)))*0.5 );
          fFitFunction->SetParLimits(4+0, fFitFunction->GetParameters()[4+0]*0.5, fFitFunction->GetParameters()[4+0]*1.5);
          
          //Position
          fFitFunction->SetParameter(4+1, Energies[PeakID]);
          // fFitFunction->SetParLimits(4+1, Energies[PeakID]-DefFWHM, Energies[PeakID]+DefFWHM);
          fFitFunction->SetParLimits(4+1, MinMax[PeakID].first, MinMax[PeakID].second);
          // In case there is a need to fix mean pos (not needed her)
          //fFitFunction->FixParameter(4+i*6+1,fEnergies[i]);
          
          //FWHM
          fFitFunction->SetParameter(4+2, DefFWHM);
          fFitFunction->SetParLimits(4+2, DefFWHM_min, DefFWHM_max);
          // In case there is a need to fix FWHM (not needed her)
          //fFitFunction->FixParameter(4+i*6+2,DefFWHM);
          
          //LeftTail
          fFitFunction->SetParameter(4+3, LeftTailVal);
          fFitFunction->SetParLimits(4+3, LeftTailValMin, LeftTailValMax);
          if(!UseLT)
              fFitFunction->FixParameter(4+3,-5);
          // In case there is a need to fix LT (not needed her)
          //fFitFunction->FixParameter(4+3,LeftTailVal);
          
          //RightTail
          fFitFunction->SetParameter(4+4, RightTailVal);
          fFitFunction->SetParLimits(4+4, RightTailValMin, RightTailValMax);
          if(!UseRT)
              fFitFunction->FixParameter(4+4,5);
          // In case there is a need to fix RT (not needed her)
          // fFitFunction->FixParameter(4+4,RightTailVal);

          //AmplitudeStep
          fFitFunction->SetParameter(4+5, StepVal);
          fFitFunction->SetParLimits(4+5, StepValMin, StepValMax);
          if(!UseStep)
              fFitFunction->FixParameter(4+5,0);
     
      // Hist->GetXaxis()->SetRangeUser(x,y);

      TString FitOpt = "QR0S";
      if(PrintLevel>0) FitOpt +="V";
      FitOpt += FitOptions;
      FitResult = Hist->Fit(fFitFunction,FitOpt.Data());
      //auto fout = new TFile("./ssd/testbetamin.root","RECREATE");
      //Hist->Write();
      //exit(1);

      // Not necessary for the moment, might as well comment it
      //if(fPlayer->fFixAmpli->GetState() == kButtonDown) {

      //    //Extract Background
      //    delete fBackFunction;
      //    fBackFunction = new TF1("Background", this, &CXFit::StepedBackground, fBackgd[0], fBackgd[1], NPars, "CXFit", "StepedBackground");
      //    fBackFunction->SetParameters(fFitFunction->GetParameters());

      //    for(auto i=0U ; i<fEnergies.size() ; i++) {
      //        Double_t Ampli = fHistogram->GetBinContent(fHistogram->FindBin(fEnergies[i]))- fBackFunction->Eval(fEnergies[i]);
      //        fFitFunction->FixParameter(4+i*6+0,Ampli);
      //    }

      //    r = fHistogram->Fit(fFitFunction,FitOpt.Data());
      //}

      if(FitResult.Get() == nullptr)
      {
        std::cout << "Fit result empty, exitting" << std::endl;
        return;
      }

      //Extract Background
      delete fBackFunction;
      fBackFunction = new TF1("Background", this, &TExogamPhysics::StepedBackground, MinMax[PeakID].first, MinMax[PeakID].second, NPars, "TExogamPhysicsFit", "StepedBackground");
      fBackFunction->SetParameters(fFitFunction->GetParameters());
      fBackFunction->SetNpx(1000);

      fBackFunction->SetLineColor(kBlue);
      fBackFunction->DrawCopy("same");

      //Extract Residue
      fHistogram = Hist;
      delete fResidue;
      fResidue = new TF1("Residue", this, &TExogamPhysics::Residue, MinMax[PeakID].first, MinMax[PeakID].second, NPars, "TExogamPhysicsFit", "Residue");
      fResidue->SetParameters(fFitFunction->GetParameters());
      fResidue->SetNpx(1000);

      fResidue->SetLineWidth(1);
      fResidue->SetLineColor(kBlack);
      fResidue->DrawCopy("same");

      ostringstream text;

      text << "Fit results :";
      text << "Status: ";
      if(FitResult->Status()==0)
          text << " Successeful" << endl;
      else
          text << " Failed" << endl;
      std::cout<<text.str();
      text.str("");

      text << "Chi2  = "<< FitResult->Chi2();
      text << "Ndf   = "<< FitResult->Ndf();
      text << "P val = "<< FitResult->Prob();
      std::cout<<text.str()<<endl;

      // fListOfPeaks->Clear();

          text<<"Peak :";
          std::cout<<text.str()<<endl;

          TF1 *peak = new TF1("Peak", this, &TExogamPhysics::PeakFunction, MinMax[PeakID].first, MinMax[PeakID].second, NPars, "TExogamPhysicsFit", "PeakFunction");
          peak->SetParameters(fFitFunction->GetParameters());
          peak->SetParErrors(fFitFunction->GetParErrors());
          peak->SetParameter(1,1);//with backgroud
          peak->SetParameter(0,PeakID);
          peak->SetNpx(1000);

          peak->SetLineColor(kGreen);
          peak->SetLineStyle(kDashed);
          peak->DrawCopy("same");


          Area[0]     = (peak->Integral(MinMax[PeakID].first,MinMax[PeakID].second,1e-6)-fBackFunction->Integral(MinMax[PeakID].first,MinMax[PeakID].second,1e-6))/fHistogram->GetBinWidth(1);
          Area[1]  = 2*sqrt(Area[0]);
          Double_t Mean     = peak->GetParameter(4+1);
          Double_t MeanErr  = peak->GetParError(4+1);
          Double_t FWHM     = peak->GetParameter(4+2);
          Double_t FWHMErr  = peak->GetParError(4+2);
//          Double_t LeftT    = TMath::Abs(peak->GetParameter(4+i*6+3));
//          Double_t LeftTErr = peak->GetParError(4+i*6+3);
//          Double_t Right    = peak->GetParameter(4+i*6+4);
//          Double_t RightErr = peak->GetParError(4+i*6+4);

          peak->SetParameter(1,0);//without backgroud
          Double_t Max      = peak->GetParameter(4+0);
          Double_t MaxErr   = peak->GetParError(4+0);

          Double_t FWHM_L     = peak->GetX(Max/2,MinMax[PeakID].first,Mean,1e-6);
          Double_t FWHM_L_err = peak->GetX((Max-MaxErr)/2,MinMax[PeakID].first,Mean,1e-6);

          Double_t F01_L     = Mean-peak->GetX(Max/10.,MinMax[PeakID].first,Mean,1e-6);
          Double_t F01_R     = peak->GetX(Max/10.,Mean,MinMax[PeakID].second,1e-6)-Mean;

          Double_t FWHM_R     = peak->GetX(Max/2,Mean,MinMax[PeakID].second,1e-6);
          Double_t FWHM_R_err = peak->GetX((Max-MaxErr)/2,Mean,MinMax[PeakID].second,1e-6);


          Double_t LeftTailParam = F01_L/(FWHM*0.5);
          Double_t RightTailParam = F01_R/(FWHM*0.5);
          std::cout<<F01_L<<" "<<FWHM<<" "<<endl;
          Double_t FWHM_Real     = FWHM_R-FWHM_L;
          Double_t FWHM_Real_err = (FWHM_R_err-FWHM_L_err)-FWHM_Real;
          peak->SetParameter(1,1);//with backgroud

          text<<left<<setw(11)<<"Mean"<<": "<<setprecision(7)<<setw(10)<<Mean<<" ("<<setprecision(7)<<setw(10)<<MeanErr<<")";
          std::cout<<text.str()<<endl;
          text.str("");

          text<<left<<setw(11)<<"Amplitude"<<": "<<setprecision(7)<<setw(10)<<Max<<" ("<<setprecision(7)<<setw(10)<<MaxErr<<")";
          std::cout<<text.str()<<endl;
          text.str("");

          text<<left<<setw(11)<<"FWHM (gaus)"<<": "<<setprecision(7)<<setw(10)<<FWHM<<" ("<<setprecision(7)<<setw(10)<<FWHMErr<<")";
          std::cout<<text.str()<<endl;
          text.str("");

          text<<left<<setw(11)<<"FWHM (real)"<<": "<<setprecision(7)<<setw(10)<<FWHM_Real<<" ("<<setprecision(7)<<setw(10)<<FWHM_Real_err<<")";
          std::cout<<text.str()<<endl;
          text.str("");

          text<<left<<setw(11)<<"L Tail"<<": "<<setprecision(7)<<setw(10)<<LeftTailParam;
          std::cout<<text.str()<<endl;
          text.str("");

          text<<left<<setw(11)<<"R Tail"<<": "<<setprecision(7)<<setw(10)<<RightTailParam;
          std::cout<<text.str()<<endl;
          text.str("");

          text<<left<<setw(11)<<"Area"<<": "<<setprecision(7)<<setw(10)<<Area[0]<<" ("<<setprecision(7)<<setw(10)<<Area[1]<<")";
          std::cout<<text.str()<<endl;
          text.str("");

      fFitFunction->DrawCopy("same");

}

double TExogamPhysics::DoubleTailedStepedGaussian(Double_t*xx,Double_t*pp)
{
    Double_t x   = xx[0];

    int    PeakID = (int)pp[0]; //Number of subpeaks in the peak range

    Double_t f_tot = 0.;

    Int_t Npar = 6;

    Double_t Back_const = pp[1];
    Double_t Back_slope = pp[2];
    Double_t Back_Exp = pp[3];

    f_tot += (Back_const + (x-MinMax[PeakID].first)*Back_slope)*exp((x-MinMax[PeakID].first)*Back_Exp);
      Double_t Ampli     = pp[4+0];
      Double_t Mean      = pp[4+1];
      Double_t Sigma     = pp[4+2]*1./sqrt(8.*log(2.));
      Double_t Lambda    = pp[4+3];
      Double_t Rho       = pp[4+4];
      Double_t S         = pp[4+5];
    
    // std::cout << Back_const << " " << Back_slope << " " << Back_Exp << " " << Ampli <<  std::endl;


      Double_t U         = (x-Mean)/Sigma;
      Double_t f_g       = Ampli*TMath::Exp(-U*U*0.5);
      Double_t f_lambda  = Ampli*TMath::Exp(-0.5*Lambda*(2.*U-Lambda));
      Double_t f_rho     = Ampli*TMath::Exp(-0.5*Rho*(2.*U-Rho));
      Double_t f_S       = Ampli*S*1./((1+TMath::Exp(U))*(1+TMath::Exp(U)));

      if(U<Lambda) f_tot += f_lambda;
      else if(U>Rho) f_tot += f_rho;
      else f_tot += f_g;

      f_tot += f_S;

    return f_tot;
}

double TExogamPhysics::StepedBackground(Double_t*xx,Double_t*pp)
{
    Double_t x   = xx[0];

    Double_t f_tot = 0.;

    int    PeakID = (int)pp[0]; //Number of subpeaks in the peak range

    Double_t Back_const = pp[1];
    Double_t Back_slope = pp[2];
    Double_t Back_Exp = pp[3];

    Int_t Npar = 6;

    f_tot += (Back_const + (x-MinMax[PeakID].first)*Back_slope)*exp((x-MinMax[PeakID].first)*Back_Exp);

        Double_t Ampli = pp[4+0];
        Double_t Mean  = pp[4+1];
        Double_t Sigma = pp[4+2]*1./sqrt(8.*log(2.));
        Double_t S     = pp[4+5];

        f_tot += Ampli*S*1./((1+TMath::Exp((x-Mean)/Sigma))*(1+TMath::Exp((x-Mean)/Sigma)));

    return f_tot;
}

double TExogamPhysics::Residue(Double_t*xx,Double_t*/*pp*/)
{
    Double_t x   = xx[0];

    return fFitFunction->Eval(fHistogram->GetBinCenter(fHistogram->FindBin(x))) - fHistogram->GetBinContent(fHistogram->FindBin(x));
}

double TExogamPhysics::PeakFunction(Double_t*xx,Double_t*pp)
{
    Double_t x   = xx[0];


    Double_t f_tot = 0.;

    Int_t Npar = 6;

    Bool_t WithBackground = true;
    if(pp[1]==0)
        WithBackground = false;

    Double_t BackGround = fBackFunction->Eval(x);

    if(WithBackground)
        f_tot += BackGround;

    Double_t Ampli     = pp[4+0];
    Double_t Mean      = pp[4+1];
    Double_t Sigma     = pp[4+2]*1./sqrt(8.*log(2.));
    Double_t Lambda    = pp[4+3];
    Double_t Rho       = pp[4+4];
    Double_t S         = 0;

    Double_t U         = (x-Mean)/Sigma;
    Double_t f_g       = Ampli*TMath::Exp(-U*U*0.5);
    Double_t f_lambda  = Ampli*TMath::Exp(-0.5*Lambda*(2.*U-Lambda));
    Double_t f_rho     = Ampli*TMath::Exp(-0.5*Rho*(2.*U-Rho));
    Double_t f_S       = Ampli*S*1./((1+TMath::Exp(U))*(1+TMath::Exp(U)));

    if(U<Lambda) f_tot += f_lambda;
    else if(U>Rho) f_tot += f_rho;
    else f_tot += f_g;

    f_tot += f_S;

    return f_tot;
}

///////////////////////////////////////////////////////////////////////////
namespace EXOGAM_LOCAL {
  //	tranform an integer to a string
  double fEXO_E(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXO/E";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoE(i),1);
  }

  double fEXO_EHG(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXO/EHG";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoEHG(i),1);
  }

  double fEXO_T(const TExogamData* m_EventData, const unsigned int& i) {
    static string name;
    name = "EXOGAM/Cr_";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    name += "_TDC";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoTDC(i),1);
  }

  double fEXO_Outer(const TExogamData* m_EventData, const unsigned int& i, const unsigned int OuterNumber) {
    static string name;
    name = "EXO/Outer";
    name += NPL::itoa(m_EventData->GetExoCrystal(i));
    name += "_";
    name += NPL::itoa(OuterNumber);
    if(OuterNumber == 0)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter1(i),1);
    else if(OuterNumber == 1)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter2(i),1);
    else if(OuterNumber == 2)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter3(i),1);
    else if(OuterNumber == 3)
      return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetExoOuter4(i),1);
    else{
      std::cout << "WARNING: Outer number != 0-3, something is wrong\n";
      return 0;
    };
  }

  string itoa(int value) {
    std::ostringstream o;

    if (!(o << value))
      return "";

    return o.str();
  }
} // namespace EXOGAM_LOCAL

/////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TExogamPhysics::Construct() { return (NPL::VDetector*)new TExogamPhysics(); }

NPL::VTreeReader* TExogamPhysics::ConstructReader() { return (NPL::VTreeReader*)new TExogamPhysicsReader(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  class proxy_exogam {
    public:
      proxy_exogam() {
        NPL::DetectorFactory::getInstance()->AddToken("Exogam", "Exogam");
        NPL::DetectorFactory::getInstance()->AddDetector("Exogam", TExogamPhysics::Construct);
        NPL::DetectorFactory::getInstance()->AddDetectorReader("Exogam", TExogamPhysics::ConstructReader);
      }
  };

  proxy_exogam p;
}
