/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Freddy Flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : December 2019                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NebulaPlus Treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNebulaPlusPhysics.h"
using namespace NEBULAPLUS_LOCAL;

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
#include "NPSystemOfUnits.h"
//   ROOT
#include "TChain.h"

ClassImp(TNebulaPlusPhysics)


  ///////////////////////////////////////////////////////////////////////////
TNebulaPlusPhysics::TNebulaPlusPhysics()
  : m_EventData(new TNebulaPlusData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_Q_RAW_Threshold(15000), // adc channels
  m_Q_Threshold(0),     // normal bars in MeV
  m_V_Threshold(0),     // veto bar in MeV
  m_TotalNbrModules(0),  
  m_BarWidth(12),     // bar width in cm
  m_BarLength(180),     // bar length in cm
  m_TdiffRange(12),     // +- range in (Td-Tu) in ns
  m_NumberOfBars(0) {
  }

/////////////////////////////////////////////
////   Innherited from VDetector Class   ////
/////////////////////////////////////////////
void TNebulaPlusPhysics::ReadConfiguration(NPL::InputParser parser) {

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEBULAPLUS");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " Nebula layers found " << endl;

  unsigned int det=0;
  // Cartesian Case
  vector<string> info
    = {"POS","NumberOfModules", "Veto", "Frame"};
  string Type; 

  for (unsigned int i = 0; i < blocks.size(); i++) {

    if (blocks[i]->HasTokenList(info)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())

      cout << endl << "//// Nebula Layer " << Type << " " << i + 1 << endl;
      int NbrOfModulesInLayer = blocks[i]->GetInt("NumberOfModules");
      //det = i+1;
      //m_DetectorNumberIndex[detectorNbr]=det;
      TVector3 Pos = blocks[i]->GetTVector3("Pos", "cm");
      AddLayer(Pos, NbrOfModulesInLayer);
    }

    else {
      cout << "ERROR: Missing token for NebulaPlus, check your input "
        "file"
        << endl;
      exit(1);
    }

  }

  //InitializeStandardParameter();
  //ReadAnalysisConfig();
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::BuildPhysicalEvent() {

  // apply thresholds and calibration
  //std::cout << "----------------" << std::endl;
  //std::cout << "in BuildPhysical" << std::endl;

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  //Cal->Dump();
  //double offset = Cal->GetValue("NebulaPlus/Tdiff_Offset_58",0); 
  //double offset = Cal->GetValue("NebulaPlus/U3_ENERGY",1); 


  static double rawQup,calQup,rawQdown,calQdown,rawTup,calTup,rawTdown,calTdown,calQ,calT,Y;
  static unsigned int ID;
  // All vector size 
  static unsigned int Usize, Dsize ; 
  Usize = m_EventData->GetMultUp();
  Dsize = m_EventData->GetMultDown();
  static double threshold;
  bool matchUD;
  int unmatched=0;

  // loop on up
  //std::cout << "size up" << Usize<< std::endl;
  //std::cout << "size down" << Dsize<< std::endl;

  for (unsigned int up = 0; up < Usize ; up++) {

    rawQup = m_EventData->GetQUp(up);
    rawTup= m_EventData->GetTUp(up);
    rawQdown=-1;
    rawTdown=-1;
    matchUD = false;
    if (rawQup > m_Q_RAW_Threshold) {

      ID = m_EventData->GetIDUp(up);

      // look for associated down
      for(unsigned int down = 0 ; down < Dsize ; down++){
        if(m_EventData->GetIDDown(down)==ID){

      	  //std::cout << "-------" << std::endl;
      	  //std::cout << "same ID:\t" << ID<< std::endl;
          if( m_EventData->GetQDown(down) > m_Q_RAW_Threshold){
              rawQdown = m_EventData->GetQDown(down); 
              rawTdown = m_EventData->GetTDown(down); 
      	      calQup   = fCalQUp(m_EventData,up);
      	      calQdown = fCalQDown(m_EventData,down);
              matchUD = true;
              break;
          }//if raw threshold down
        } //if match ID 
      }// down 

      // Got a Down-Up match, do the math
      if(matchUD){

        // cal Q Up and Down
        //calQup = aQu[ID] * (rawQup-bQu[ID]);
        //calQdown= aQd[ID] * (rawQdown-bQd[ID]);
        
        // average value of Up and Down
        calQ=sqrt(calQup*calQdown); 

        // cal T  Up
        //calTup=aTu[ID]*rawTup+bTu[ID];
        calTup=rawTup;

        // cal T Down
        //calTdown=aTd[ID]*rawTdown+bTd[ID];
        calTdown=rawTdown;

        if(calQ>threshold){
          calT= (calTdown+calTup)*0.5; 
          Y=(calTdown-calTup) - Cal->GetValue("NebulaPlus/Tdiff_Offset_"+NPL::itoa(ID),0);

          DetectorNumber.push_back(ID);
          Charge.push_back(calQ);
          TOF.push_back(calT);
          PosX.push_back(m_BarPositionX[ID-1]);
          PosY.push_back(Y/m_TdiffRange*m_BarLength);
          PosZ.push_back(m_BarPositionZ[ID-1]);
/*
          if(ID<121)
            IsVeto.push_back(0);
          else
            IsVeto.push_back(1);
*/
        }
      } else {
          unmatched++;/*std::cout << "unmatched=" <<unmatched<< std::endl;*/
      }
    }// if raw threshold up
  } // Qup

}

///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::PreTreat() {

}

///////////////////////////////////////////////////////////////////////////

void TNebulaPlusPhysics::AddLayer(TVector3 Pos, int n_modules) {

    for(int i=0; i<n_modules; i++){
        m_BarPositionX.push_back(Pos.X() + i * m_BarWidth);
        m_BarPositionY.push_back(Pos.Y());
        m_BarPositionZ.push_back(Pos.Z());
        m_TotalNbrModules++;
    }
}
///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::ReadAnalysisConfig() {
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::Clear() {
  DetectorNumber.clear();
  Charge.clear();
  TOF.clear();
  PosY.clear();
  PosX.clear();
  PosZ.clear();
  IsVeto.clear();
}



///////////////////////////////////////////////////////////////////////////
/*
void TNebulaPlusPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEBULAPLUS");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","Offset","InvertX","InvertY"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      cout << endl << "////  NebulaPlus (" << i+1 << ")" << endl;
      unsigned int det = std::atoi(blocks[i]->GetMainValue().c_str());
      string xmlpath = blocks[i]->GetString("XML");
      NPL::XmlParser xml;
      xml.LoadFile(xmlpath);
      ReadXML(xml);
      TVector3 offset = blocks[i]->GetTVector3("Offset","mm"); 
      bool invertX = blocks[i]->GetInt("InvertX"); 
      bool invertY = blocks[i]->GetInt("InvertY"); 
      m_offset[det] = offset;
      m_invertX[det] = invertX;
      m_invertY[det] = invertY;
    }
  }
}
*/
///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::InitSpectra() {
  m_Spectra = new TNebulaPlusSpectra(m_NumberOfBars);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TNebulaPlusPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  vector<double> standard={8192, 6};
  for (int i = 0; i < m_TotalNbrModules; ++i) {
    Cal->AddParameter("NebulaPlus","U"+ NPL::itoa(i+1)+"_ENERGY","NebulaPlus_U"+ NPL::itoa(i+1)+"_ENERGY",standard);
    Cal->AddParameter("NebulaPlus","D"+ NPL::itoa(i+1)+"_ENERGY","NebulaPlus_D"+ NPL::itoa(i+1)+"_ENERGY",standard);
    Cal->AddParameter("NebulaPlus","Tdiff_Offset_"+ NPL::itoa(i+1),"NebulaPlus_Tdiff_Offset_"+ NPL::itoa(i+1),standard);
  }
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("NebulaPlus",  true );
  inputChain->SetBranchAddress("NebulaPlus", &m_EventData );
}

///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("NebulaPlus", &m_EventPhysics);
}


///////////////////////////////////////////////////////////////////////////
void TNebulaPlusPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree("NebulaPlus");
  outputTree->Branch("NebulaPlus", "TNebulaPlusPhysics", &m_EventPhysics);
}

////////////////////////////////////////////////////////////////////////////////
namespace NEBULAPLUS_LOCAL{

double fCalQUp(const TNebulaPlusData* m_EventData, const int& i){
    static string name;
    name = "NebulaPlus/U";
    name += NPL::itoa(m_EventData->GetIDUp(i));
    name += "_ENERGY";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetQUp(i), 1);
}

double fCalQDown(const TNebulaPlusData* m_EventData, const int& i){
    static string name;
    name = "NebulaPlus/D";
    name += NPL::itoa(m_EventData->GetIDDown(i));
    name += "_ENERGY";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetQDown(i), 1);
}

}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TNebulaPlusPhysics::Construct() {
  return (NPL::VDetector*) new TNebulaPlusPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_NebulaPlus{
    public:
      proxy_NebulaPlus(){
        NPL::DetectorFactory::getInstance()->AddToken("NEBULAPLUS","NebulaPlus");
        NPL::DetectorFactory::getInstance()->AddDetector("NEBULAPLUS",TNebulaPlusPhysics::Construct);
      }
  };

  proxy_NebulaPlus p_NebulaPlus;
}

