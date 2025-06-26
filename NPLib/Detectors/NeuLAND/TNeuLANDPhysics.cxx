/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey ANNE  contact address: anne@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : November 2023                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NeuLAND Treated data                                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TNeuLANDPhysics.h"

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

ClassImp(TNeuLANDPhysics)


  ///////////////////////////////////////////////////////////////////////////
TNeuLANDPhysics::TNeuLANDPhysics()
  : m_EventData(new TNeuLANDData),
  m_EventPhysics(this),
  m_Spectra(0),
  m_Q_RAW_Threshold(0), // adc channels
  m_Q_Threshold(7),     // normal bars in MeV
  m_V_Threshold(1),     // veto bars in MeV
  m_NumberOfBars(0)
  {
  m_Material_Index = 1.58;
  m_Signal_Effective_Velocity = (3.0e8/m_Material_Index)*NPUNITS::m/NPUNITS::s;  
  }

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
void TNeuLANDPhysics::ReadXML(NPL::XmlParser xml){ 
  std::vector<NPL::XML::block*> b = xml.GetAllBlocksWithName("NEULAND");  

  for(unsigned int i = 0 ; i < b.size() ; i++){
    m_NumberOfBars++;
    unsigned int id = b[i]->AsInt("ID");

    // position
    PositionX[id] = b[i]->AsDouble("xpos");
    PositionY[id] = b[i]->AsDouble("ypos"); 
    PositionZ[id] = b[i]->AsDouble("zpos");

    // bar direction
    DirectionBar[id] = b[i]->AsString("direction");
    
    std::string tree_name = RootInput::getInstance()->GetChain()->GetName();

    if(tree_name!="SimulatedTree"){
      // linear cal
      aQu[id] = b[i]->AsDouble("QUCal");
      bQu[id] = b[i]->AsDouble("QUPed");
      aQd[id] = b[i]->AsDouble("QDCal");
      bQd[id] = b[i]->AsDouble("QDPed");
      aTu[id] = b[i]->AsDouble("TUCal");
      bTu[id] = b[i]->AsDouble("TUOff");
      aTd[id] = b[i]->AsDouble("TDCal");
      bTd[id] = b[i]->AsDouble("TDOff");

      // T average offset
      avgT0[id] = b[i]->AsDouble("TAveOff");

      // slew correction T= tcal +slwT/sqrt(Qcal)
      slwTu[id] = b[i]->AsDouble("TUSlw");
      slwTd[id] = b[i]->AsDouble("TDSlw");

      // DT position cal
      DTa[id] = b[i]->AsDouble("DTCal");//!
      DTb[id] = b[i]->AsDouble("DTOff");//!
    }
    else{
      // linear cal
      aQu[id] =1;
      bQu[id] =0; 
      aQd[id] =1;
      bQd[id] =0; 
      aTu[id] =1; 
      bTu[id] =0; 
      aTd[id] =1; 
      bTd[id] =0; 

      // T average offset
      //avgT0[id]= 0.0;
      avgT0[id]= -6.5833;

      // slew correction T= tcal +slwT/sqrt(Qcal)
      slwTu[id] = 0;
      slwTd[id] = 0;

      // DT position cal
      // DTa[id] = 1;//!
      DTb[id] = 0;//!

      DTa[id] = 0.5*m_Signal_Effective_Velocity;//!
      //DTa[id] = 94.87;//!


      }

  } 
  cout << " -> " << m_NumberOfBars << " bars found" << endl;
} 

///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::BuildPhysicalEvent() {
/*  // apply thresholds and calibration

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();
  static double rawQup,calQup,rawQdown,calQdown,rawTup,calTup,rawTdown,calTdown,calQ,calT,X,Y;
  static unsigned int ID;

  static double neutron_mass = 939.57; //  MeV/c^2
  static double c_light = 299792458;    //m.s-1

  // All vector size 
  static unsigned int QUsize, QDsize, TUsize, TDsize ; 
  QUsize = m_EventData->GetChargeUpMult();
  QDsize = m_EventData->GetChargeDownMult();
  TUsize = m_EventData->GetTimeUpMult();
  TDsize = m_EventData->GetTimeDownMult();
  static double threshold;
  // loop on Qup
  for (unsigned int qup = 0; qup < QUsize ; qup++) {
    rawQup = m_EventData->GetChargeUp(qup);
    rawTup=-1;
    rawQdown=-1;
    rawTdown=-1;
    if (rawQup > m_Q_RAW_Threshold) {
      ID = m_EventData->GetChargeUpID(qup);
      if(ID<401)
        threshold=m_Q_Threshold;
      else
        threshold=m_V_Threshold;

      // look for associated Charge down
      for(unsigned int qdown = 0 ; qdown < QDsize ; qdown++){
        if(m_EventData->GetChargeDownID(qdown)==ID){
          rawQdown=m_EventData->GetChargeDown(qdown); 
          if(rawQdown > m_Q_RAW_Threshold)
	    {
	      // Look for the associate time 
	      for(unsigned int tdown = 0 ; tdown < TDsize; tdown++)
		{
		  if(m_EventData->GetTimeDownID(qdown)==ID)
		    {
		      rawTdown=m_EventData->GetTimeDown(qdown);
		      break;
		    }
		}// TDown
	    }//if raw threshold down
          break;
        } //if match ID 

      }// Qdwown 

      if(rawTdown>0){ // Tdown is found, means Qdown as well
        // look for Tup  
        for(unsigned int tup = 0 ; tup < TUsize ; tup++){
          if(m_EventData->GetTimeUpID(tup)==ID){
            rawTup = m_EventData->GetTimeUp(tup);
            break;
          }
        }
      }
      // Got everything, do the math
      if(rawTup>0){
        // cal Q Up and Down
        calQup=aQu[ID]*(rawQup-bQu[ID]);
        calQdown=aQd[ID]*(rawQdown-bQd[ID]);

        // average value of Up and Down
        calQ=sqrt(calQup*calQdown); 

        // cal T  Up
        calTup=aTu[ID]*rawTup+bTu[ID];
        // slew correction
        calTup -= slwTu[ID]/sqrt(rawQup-bQu[ID]);

        // cal T Down
        calTdown=aTd[ID]*rawTdown+bTd[ID];
        // slew correction
        calTdown -= slwTd[ID]/sqrt(rawQdown-bQd[ID]);
    
        if(calQ>threshold){
          calT= (calTdown+calTup)*0.5+avgT0[ID]+Cal->GetPedestal("NEULAND_T_ID"+NPL::itoa(ID));

	  if(DirectionBar[ID]=='H')
	    {
	      X = (calTdown-calTup)*DTa[ID]+DTb[ID]+Cal->GetPedestal("NEULAND_X_ID"+NPL::itoa(ID));
	      Y = 0.0;
	    }
	      
	  if(DirectionBar[ID]=='V')
	    {
	      X = 0.0;
	      Y=(calTdown-calTup)*DTa[ID]+DTb[ID]+Cal->GetPedestal("NEULAND_Y_ID"+NPL::itoa(ID));     
	    }
	   
          DetectorNumber.push_back(ID);
          Charge.push_back(calQ);
          TOF.push_back(calT);
          PosY.push_back(Y+PositionY[ID]);
          PosX.push_back(X+PositionX[ID]);
          PosZ.push_back(PositionZ[ID] + m_offset[0][2]); // [0]-> only one Nebula detector    [2]-> z position

          if(ID<401)
	    {
            IsVeto.push_back(0);
	    }
          else
            IsVeto.push_back(1);

	  X = 0.0;
	  Y = 0.0;

        }
      }

    }// if raw threshold up
  } // Qup
*/
}

///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::PreTreat() {

}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::ReadAnalysisConfig() {
}

///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::Clear() {
  DetectorNumber.clear();
  Charge.clear();
  TOF.clear();
  PosY.clear();
  PosX.clear();
  PosZ.clear();
  IsVeto.clear();
}

///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("NEULAND");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " detector(s) found " << endl; 

  vector<string> token= {"XML","Offset","InvertX","InvertY"};

  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      cout << endl << "////  NeuLAND (" << i+1 << ")" << endl;
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

///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::InitSpectra() {
  m_Spectra = new TNeuLANDSpectra(m_NumberOfBars);
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::ClearSpectra() {
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TNeuLANDPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();

  vector<double> standardO={0};
  for (int i = 0; i < m_NumberOfBars; ++i) {
    Cal->AddParameter("NEULAND_T_ID"+ NPL::itoa(i+1),standardO);
    Cal->AddParameter("NEULAND_X_ID"+ NPL::itoa(i+1),standardO);
    Cal->AddParameter("NEULAND_Y_ID"+ NPL::itoa(i+1),standardO);
  }
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("NeuLAND",  true );
  inputChain->SetBranchAddress("NeuLAND", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("NeuLAND", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TNeuLANDPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree("NeuLAND");
  outputTree->Branch("NeuLAND", "TNeuLANDPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TNeuLANDPhysics::Construct() {
  return (NPL::VDetector*) new TNeuLANDPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_NeuLAND{
    public:
      proxy_NeuLAND(){
        NPL::DetectorFactory::getInstance()->AddToken("NEULAND","NeuLAND");
        NPL::DetectorFactory::getInstance()->AddDetector("NEULAND",TNeuLANDPhysics::Construct);
      }
  };

  proxy_NeuLAND p_NeuLAND;
}

