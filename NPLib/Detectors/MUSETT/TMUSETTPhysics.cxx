/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : June 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold MUSETT Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include "TMUSETTPhysics.h"
using namespace MUSETT_LOCAL;

//   STL
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <cstdlib>
#include <time.h>

//   NPL
#include "NPDetectorFactory.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
using namespace NPUNITS;

//   ROOT
#include "TChain.h"
#include "TMath.h"

///////////////////////////////////////////////////////////////////////////

ClassImp(TMUSETTPhysics)
  ///////////////////////////////////////////////////////////////////////////
  TMUSETTPhysics::TMUSETTPhysics() {
    EventMultiplicity                  = 0;
    EventTrigger                       = 0;
    m_EventData                        = new TMUSETTData;
    m_PreTreatedData                   = new TMUSETTData;
    m_random                           = new TRandom3();
    m_EventPhysics                     = this;
    // m_Spectra                          = NULL;
    m_NumberOfDetector                = 0;
    m_MaximumStripMultiplicityAllowed  = 10;
    m_StripEnergyMatching = 0.050;
    // Raw Threshold
    m_DSSD_X_E_RAW_Threshold = 8200;
    m_DSSD_Y_E_RAW_Threshold = 8200;
    m_SecondLayer_E_RAW_Threshold = 8200;
    // Calibrated Threshold
    m_DSSD_X_E_Threshold = 0;
    m_DSSD_Y_E_Threshold = 0;
    m_SecondLayer_E_Threshold  = 0;

    m_Take_E_Y = false;
    m_Take_T_Y = true;

  }

///////////////////////////////////////////////////////////////////////////
TMUSETTPhysics::~TMUSETTPhysics() {}
///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

//////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::PreTreat() {
  ClearPreTreatedData();
  static unsigned int DSSDX_EMult, DSSDY_EMult, SecondLayer_EMult;
  static unsigned int DSSDX_TMult, DSSDY_TMult, SecondLayer_TMult;
  DSSDX_EMult = m_EventData->GetDSSDXEMult();
  DSSDY_EMult = m_EventData->GetDSSDYEMult();
  DSSDX_TMult = m_EventData->GetDSSDXTMult();
  DSSDY_TMult = m_EventData->GetDSSDYTMult();
  m_PreTreatedData->SetTRIGGER(m_EventData->GetDSSDTRIG());
  //   X
  //   E

  for (unsigned int i = 0; i < DSSDX_EMult; ++i) {
    if (m_EventData->GetDSSDXEEnergy(i) > m_DSSD_X_E_RAW_Threshold){
        //&& IsValidChannel(0, m_EventData->GetDSSDXEDetectorNbr(i),
          //m_EventData->GetDSSDXEStripNbr(i))) {
      double EX = fDSSD_X_E(m_EventData, i);
      if (EX > m_DSSD_X_E_Threshold)
      {

        //Fo r DATA
        /*
        m_PreTreatedData->SetDSSDXE(true,
            m_EventData->GetDSSDXEDetectorNbr(i),
            m_EventData->GetDSSDXEStripNbr(i), EX);
        */



        //For Simu
        m_PreTreatedData->SetDSSDXE(false,
            m_EventData->GetDSSDXEDetectorNbr(i),
            m_EventData->GetDSSDXEStripNbr(i)-1, EX);

      }
    }
  }

  //   T
  for (unsigned int i = 0; i < DSSDX_TMult; ++i) {
    //if (IsValidChannel(0, m_EventData->GetDSSDXTDetectorNbr(i),
    //      m_EventData->GetDSSDXTStripNbr(i)))

    //For Data
    /*
      m_PreTreatedData->SetDSSDXT(true,
          m_EventData->GetDSSDXTDetectorNbr(i),
          m_EventData->GetDSSDXTStripNbr(i),
          fDSSD_X_T(m_EventData, i));
    */
    //For Simu

      m_PreTreatedData->SetDSSDXT(false,
          m_EventData->GetDSSDXTDetectorNbr(i),
          m_EventData->GetDSSDXTStripNbr(i)-1,
          fDSSD_X_T(m_EventData, i));

  }

  //   Y
  //   E
  for (unsigned int i = 0; i < DSSDY_EMult; ++i) {
    //std::cout << "EY = " << m_EventData->GetDSSDYEEnergy(i) ;
    if (m_EventData->GetDSSDYEEnergy(i) < m_DSSD_Y_E_RAW_Threshold){
        //&& IsValidChannel(1, m_EventData->GetDSSDYEDetectorNbr(i),
          //m_EventData->GetDSSDYEStripNbr(i))) {
          //std::cout << "-> OK  " << std::endl;
      double EY = fDSSD_Y_E(m_EventData, i);
      if (EY > m_DSSD_Y_E_Threshold)

      //For DATA
      /*
      m_PreTreatedData->SetDSSDYE(true,
          m_EventData->GetDSSDYEDetectorNbr(i),
          m_EventData->GetDSSDYEStripNbr(i), EY);

      */
      //For Simu
        m_PreTreatedData->SetDSSDYE(false,
            m_EventData->GetDSSDYEDetectorNbr(i),
            m_EventData->GetDSSDYEStripNbr(i)-1, EY);

    }
  }

  //   T
  for (unsigned int i = 0; i < DSSDY_TMult; ++i) {
    //if (IsValidChannel(1, m_EventData->GetDSSDYTDetectorNbr(i),
    //      m_EventData->GetDSSDYTStripNbr(i)))
    double TY = fDSSD_Y_T(m_EventData, i);

    //For Data
    /*
      m_PreTreatedData->SetDSSDYT(true,
          m_EventData->GetDSSDYTDetectorNbr(i),
          m_EventData->GetDSSDYTStripNbr(i),TY);


    */
    //For Simu
    m_PreTreatedData->SetDSSDYT(false,
        m_EventData->GetDSSDYTDetectorNbr(i),
        m_EventData->GetDSSDYTStripNbr(i)-1,TY);
    

  }

  return;
}


///////////////////////////////////////////////////////////////////////////

void TMUSETTPhysics::BuildPhysicalEvent() {
  PreTreat();
  bool check_SecondLayer  = false;
  static unsigned int DSSDXEMult, DSSDYEMult, DSSDXTMult, DSSDYTMult,SecondLayerEMult,SecondLayerTMult;
  DSSDXEMult = m_PreTreatedData->GetDSSDXEMult();
  DSSDYEMult = m_PreTreatedData->GetDSSDYEMult();
  DSSDXTMult = m_PreTreatedData->GetDSSDXTMult();
  DSSDYTMult = m_PreTreatedData->GetDSSDYTMult();
  // random->SetSeed(42);

  // srand(time(NULL));

  if (1 /*CheckEvent() == 1*/) {
    vector<TVector2> couple = Match_X_Y();
    EventTrigger = m_PreTreatedData->GetDSSDTRIG();
    EventMultiplicity = couple.size();
    for (unsigned int i = 0; i < couple.size(); ++i) {
      check_SecondLayer  = false;

      int N = m_PreTreatedData->GetDSSDXEDetectorNbr(couple[i].X());

      int X = m_PreTreatedData->GetDSSDXEStripNbr(couple[i].X());
      int Y = m_PreTreatedData->GetDSSDYEStripNbr(couple[i].Y());
      double DSSD_X_E = m_PreTreatedData->GetDSSDXEEnergy(couple[i].X());
      double DSSD_Y_E = m_PreTreatedData->GetDSSDYEEnergy(couple[i].Y());
      /*
      int TRIG_X_E = m_PreTreatedData->GETDSSDXETrig(couple[i].X());
      int TRIG_X_T = m_PreTreatedData->GETDSSDXTTrig(couple[i].X());
      int TRIG_Y_E = m_PreTreatedData->GETDSSDYETrig(couple[i].Y());
      int TRIG_Y_T = m_PreTreatedData->GETDSSDYTTrig(couple[i].Y());
      */

      //  Search for associate Time
      double DSSD_X_T = -1000;
      for (unsigned int t = 0; t < DSSDXTMult; ++t) {
        if (m_PreTreatedData->GetDSSDXTStripNbr(couple[i].X())
            == m_PreTreatedData->GetDSSDXTStripNbr(t)
            && m_PreTreatedData->GetDSSDXTDetectorNbr(couple[i].X())
            == m_PreTreatedData->GetDSSDXTDetectorNbr(t)) {
          DSSD_X_T = m_PreTreatedData->GetDSSDXTTime(t);
          //  std::cout << "X :  " << DSSD_X_T << std::endl;
          break;
        }
      }
      double DSSD_Y_T = -1000;
      for (unsigned int t = 0; t < DSSDYTMult; ++t) {
        if (m_PreTreatedData->GetDSSDYTStripNbr(couple[i].Y())
            == m_PreTreatedData->GetDSSDYTStripNbr(t)
            && m_PreTreatedData->GetDSSDYTDetectorNbr(couple[i].Y())
            == m_PreTreatedData->GetDSSDYTDetectorNbr(t)) {
          DSSD_Y_T = m_PreTreatedData->GetDSSDYTTime(t);
          //cout << DSSD_Y_T << endl;
          break;
        }
      }
      DSSD_X.push_back(X);
      DSSD_Y.push_back(Y);
      DetectorNumber.push_back(N);


      PosX.push_back(GetPositionOfInteraction(i).x());
      PosY.push_back(GetPositionOfInteraction(i).y());
      PosZ.push_back(GetPositionOfInteraction(i).z());


      DSSD_E_X.push_back(DSSD_X_E);
      DSSD_E_Y.push_back(DSSD_Y_E);
      DSSD_T_X.push_back(DSSD_X_T);
      DSSD_T_Y.push_back(DSSD_Y_T);

      /*
      DSSD_E_X_TRIG.push_back(TRIG_X_E);
      DSSD_E_Y_TRIG.push_back(TRIG_Y_E);
      DSSD_T_X_TRIG.push_back(TRIG_X_T);
      DSSD_T_Y_TRIG.push_back(TRIG_Y_T);
      */

      if (m_Take_E_Y){
        DSSD_E.push_back(DSSD_Y_E);
      }
      else{
        DSSD_E.push_back(DSSD_X_E);
      }

      if (m_Take_T_Y)
        DSSD_T.push_back(DSSD_Y_T);
      else{
        DSSD_T.push_back(DSSD_X_T);
      }

    } // loop on couples

  } // if (CheckEvent)
  return;
}

///////////////////////////////////////////////////////////////////////////
int TMUSETTPhysics::CheckEvent() {
  // Check the size of the different elements
  if (m_PreTreatedData->GetDSSDXEMult()
      == m_PreTreatedData->GetDSSDYEMult())
    return 1; // Regular Event

  // INterstrip management is not coded, so waste of time to make this test
  /*  else if(   m_PreTreatedData->GetMMStripXEMult() ==
      m_PreTreatedData->GetMMStripYEMult()+1
      || m_PreTreatedData->GetMMStripXEMult() ==
      m_PreTreatedData->GetMMStripYEMult()-1  )
      return 2 ; // Pseudo Event, potentially interstrip*/

  else
    return -1; // Rejected Event
}


///////////////////////////////////////////////////////////////////////////
vector<TVector2> TMUSETTPhysics::Match_X_Y() {
  vector<TVector2> ArrayOfGoodCouple;
  static unsigned int m_DSSDXEMult,m_DSSDYEMult;
  m_DSSDXEMult = m_PreTreatedData->GetDSSDXEMult();
  m_DSSDYEMult = m_PreTreatedData->GetDSSDYEMult();

  // Prevent code from treating very high multiplicity Event
  // Those event are not physical anyway and that improve speed.
  if (m_DSSDXEMult > m_MaximumStripMultiplicityAllowed
      || m_DSSDYEMult > m_MaximumStripMultiplicityAllowed) {
    return ArrayOfGoodCouple;
  }


  //std::cout << "multX = " << m_DSSDXEMult << ", multY = " << m_DSSDYEMult << std::endl;
  for (unsigned int i = 0; i < m_DSSDXEMult; i++) {
    for (unsigned int j = 0; j < m_DSSDYEMult; j++) {

      // Declaration of variable for clarity
      double DSSDXDetNbr = m_PreTreatedData->GetDSSDXEDetectorNbr(i);
      double DSSDYDetNbr = m_PreTreatedData->GetDSSDYEDetectorNbr(j);

      //   if same detector check energy
      if (DSSDXDetNbr == DSSDYDetNbr) {

        // Declaration of variable for clarity
        double DSSDXEnergy = m_PreTreatedData->GetDSSDXEEnergy(i);
        double DSSDXNbr    = m_PreTreatedData->GetDSSDXEStripNbr(i);
        double DSSDYEnergy = m_PreTreatedData->GetDSSDYEEnergy(j);
        double DSSDYNbr    = m_PreTreatedData->GetDSSDYEStripNbr(j);


        //   Look if energy match
        if (abs((DSSDXEnergy - DSSDYEnergy) / 2.)
            < m_StripEnergyMatching) {
          // Gives a unique ID for every telescope and strip combination
          int IDX = m_NumberOfDetector * DSSDXNbr + DSSDXDetNbr;
          int IDY = m_NumberOfDetector * DSSDYNbr + DSSDYDetNbr;

          m_HitDSSDX[IDX]++;
          m_HitDSSDY[IDY]++;

          ArrayOfGoodCouple.push_back(TVector2(i, j));
        }
      }
    }
  }

  //std::cout << "Size before treatment : " << ArrayOfGoodCouple.size() << std::endl;
  //m_all+=ArrayOfGoodCouple.size();
  // Prevent to treat event with ambiguous matching beetween X and Y
  /*
  map<int, int>::iterator itX = m_HitDSSDX.begin();
  for (; itX != m_HitDSSDX.end(); itX++) {
    if (itX->second > 1) {
      m_rejected+=ArrayOfGoodCouple.size();
      ArrayOfGoodCouple.clear();
    }
  }

  //std::cout << "    Size after X : " << ArrayOfGoodCouple.size() << std::endl;


  map<int, int>::iterator itY = m_HitDSSDY.begin();
  for (; itY != m_HitDSSDY.end(); itY++) {
    if (itY->second > 1) {
      m_rejected+=ArrayOfGoodCouple.size();
      ArrayOfGoodCouple.clear();
    }
  }
  */

  //std::cout << "    Size after Y : " << ArrayOfGoodCouple.size() << std::endl;
  //std::cout << "-> rej/all = " << (double) m_rejected/m_all << std::endl;
  //std::cout << std::endl;

  m_HitDSSDX.clear();
  m_HitDSSDY.clear();

  return ArrayOfGoodCouple;
}

////////////////////////////////////////////////////////////////////////////
bool TMUSETTPhysics::IsValidChannel(const int& Type, const int& telescope, const int& channel) {

  if (Type == 0 && channel >= StripLimit && channel <= 127 - StripLimit){
    return *(m_XChannelStatus[telescope].begin() + channel);
  }
  else if (Type == 1&& channel >= StripLimit && channel <= 127 - StripLimit){
    return *(m_YChannelStatus[telescope].begin() + channel);
  }

  else
    return false;
}

///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::ReadAnalysisConfig() {

  NPL::InputParser parser("./configs/ConfigMUSETT.dat",false);
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("ConfigMUSETT");

  cout << endl << "//// Read MUSETT analysis configuration" <<endl;

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if(blocks[i]->HasToken("MAX_STRIP_MULTIPLICITY"))
      m_MaximumStripMultiplicityAllowed = blocks[i]->GetInt("MAX_STRIP_MULTIPLICITY");

    if(blocks[i]->HasToken("STRIP_ENERGY_MATCHING"))
      m_StripEnergyMatching = blocks[i]->GetDouble("STRIP_ENERGY_MATCHING","MeV");

    if(blocks[i]->HasToken("DISABLE_CHANNEL_X")){
      vector<int> v = blocks[i]->GetVectorInt("DISABLE_CHANNEL_X");
      *(m_XChannelStatus[v[0]].begin() + v[1] - 1) = false;
    }

    if(blocks[i]->HasToken("DISABLE_CHANNEL_Y")){
      vector<int> v = blocks[i]->GetVectorInt("DISABLE_CHANNEL_Y");
      *(m_YChannelStatus[v[0]].begin() + v[1] - 1) = false;
    }

    if(blocks[i]->HasToken("DISABLE_ALL")){
      int telescope = blocks[i]->GetInt("DISABLE_ALL");
      vector<bool> ChannelStatus;
      ChannelStatus.resize(128, false);
      m_XChannelStatus[telescope] = ChannelStatus;
      m_YChannelStatus[telescope] = ChannelStatus;
      ChannelStatus.resize(16, false);
    }

    if (blocks[i]->HasToken("TAKE_E_Y"))
      m_Take_E_Y = blocks[i]->GetInt("TAKE_E_Y");

    if (blocks[i]->HasToken("TAKE_T_Y"))
      m_Take_T_Y = blocks[i]->GetInt("TAKE_T_Y");

    if (blocks[i]->HasToken("TAKE_E_X"))
      m_Take_E_Y = !(blocks[i]->GetInt("TAKE_E_X"));

    if (blocks[i]->HasToken("TAKE_T_X"))
      m_Take_T_Y = !(blocks[i]->GetInt("TAKE_T_X"));

    if (blocks[i]->HasToken("DSSD_X_E_RAW_THRESHOLD"))
      m_DSSD_X_E_RAW_Threshold = blocks[i]->GetInt("DSSD_X_E_RAW_THRESHOLD");

    if (blocks[i]->HasToken("DSSD_Y_E_RAW_THRESHOLD"))
      m_DSSD_Y_E_RAW_Threshold = blocks[i]->GetInt("DSSD_Y_E_RAW_THRESHOLD");

    if (blocks[i]->HasToken("SECONDLAYER_E_RAW_THRESHOLD"))
      m_SecondLayer_E_RAW_Threshold = blocks[i]->GetInt("SECONDLAYER_E_RAW_THRESHOLD");

    if (blocks[i]->HasToken("DSSD_X_E_THRESHOLD"))
      m_DSSD_X_E_Threshold = blocks[i]->GetDouble("DSSD_X_E_THRESHOLD","MeV");

    if (blocks[i]->HasToken("DSSD_Y_E_THRESHOLD"))
      m_DSSD_Y_E_Threshold = blocks[i]->GetDouble("DSSD_Y_E_THRESHOLD","MeV");

  }
}


///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::Clear() {
  EventMultiplicity = 0;
  EventTrigger = 0;
  DetectorNumber.clear();

  PosX.clear();
  PosY.clear();
  PosZ.clear();

  // DSSD
  DSSD_E.clear();
  DSSD_E_X.clear();
  DSSD_E_Y.clear();
  DSSD_T.clear();
  DSSD_T_X.clear();
  DSSD_T_Y.clear();
  DSSD_X.clear();
  DSSD_Y.clear();
  /*
  DSSD_E_X_TRIG.clear();
  DSSD_E_Y_TRIG.clear();
  DSSD_T_X_TRIG.clear();
  DSSD_T_Y_TRIG.clear();
  */

}

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("MUSETT");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " MUSETT found " << endl;

  // Cartesian Case
  vector<string> cart = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128"};
  // Spherical Case
  vector<string> sphe = {"R", "THETA", "PHI", "BETA"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  MUSETT Detector " << i + 1 << endl;
      TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
      TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
      TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
      TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");
      AddDetector(A, B, C, D);
    }

    else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  MUSETT Detector " << i + 1 << endl;

      double Theta = blocks[i]->GetDouble("THETA", "deg");
      double Phi = blocks[i]->GetDouble("PHI", "deg");
      double R = blocks[i]->GetDouble("R", "mm");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA", "deg");
      cout << Theta << " " << Phi << endl;
      AddDetector(Theta, Phi, R, beta[0], beta[1], beta[2]);
    }

    else {
      cout << "ERROR: Missing token for MUSETT blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }

  InitializeStandardParameter();
  ReadAnalysisConfig();
}
//////////////////////////////////////////////////////////////////////////
//void TMUSETTPhysics::InitSpectra() {
//  m_Spectra = new TMUSETTSpectra(m_DetectorNumberIndex);
//}

///////////////////////////////////////////////////////////////////////////
/*void TMUSETTPhysics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::CheckSpectra() { m_Spectra->CheckSpectra(); }
///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::ClearSpectra() {
  // To be done
}

///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::WriteSpectra() {
  if (m_Spectra)
    m_Spectra->WriteSpectra();
}

///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TMUSETTPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}
*/
///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  // Good for simulation, close to typical values
  vector<double> standardX    = {-63, 63. / 8192.};
  vector<double> standardY    = {63, -63. / 8192.};
  vector<double> standardSecondLayer  = {-63, 63. / 8192.};
  vector<double> standardT    = {-1000, 1000. / 8192.};



  for (int i = 0; i < m_NumberOfDetector; i++) {

    for (int j = 0; j < 128; j++) {
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_E",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_E",
          standardX);
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_E",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_E",
          standardY);
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_T",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_X" + NPL::itoa(j ) + "_T",
          standardT);
      Cal->AddParameter(
          "MUSETT", "T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_T",
          "MUSETT_T" + NPL::itoa(i) + "_DSSD_Y" + NPL::itoa(j ) + "_T",
          standardT);
    }

  }

  return;
}

///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("MUSETT", true);
  inputChain->SetBranchStatus("fMM_*", true);
  inputChain->SetBranchAddress("MUSETT", &m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("MUSETT", true);
  inputChain->SetBranchStatus("EventMultiplicity", true);
  inputChain->SetBranchStatus("EventTrigger", true);
  inputChain->SetBranchStatus("EventType", true);
  inputChain->SetBranchStatus("DetectorNumber", true);
  inputChain->SetBranchStatus("Si_E", true);
  inputChain->SetBranchStatus("Si_T", true);
  inputChain->SetBranchStatus("Si_X", true);
  inputChain->SetBranchStatus("Si_Y", true);
  inputChain->SetBranchStatus("Si_EX", true);
  inputChain->SetBranchStatus("Si_TX", true);
  inputChain->SetBranchStatus("Si_EY", true);
  inputChain->SetBranchStatus("Si_TY", true);
  inputChain->SetBranchStatus("DetectorNumber_X", true);
  inputChain->SetBranchStatus("DetectorNumber_Y", true);
  inputChain->SetBranchStatus("TotalEnergy", true);
  inputChain->SetBranchAddress("MUSETT", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("MUSETT", "TMUSETTPhysics", &m_EventPhysics);
}

/////   Specific to MUSETTArray   ////
void TMUSETTPhysics::AddDetector(TVector3 C_X1_Y1, TVector3 C_X128_Y1,
    TVector3 C_X1_Y128, TVector3 C_X128_Y128) {
  // To avoid warning
  C_X128_Y128 *= 1;

  m_NumberOfDetector++;

  // Vector U parallel to BaseLarge
  TVector3 U = C_X128_Y1 - C_X1_Y1;
  U = U.Unit();

  // Vector V parallel to height
  TVector3 V = 0.5 * (C_X1_Y128 + C_X128_Y128 - C_X1_Y1 - C_X128_Y1);
  V = V.Unit();

  //   Position Vector of Strip Center
  TVector3 StripCenter = TVector3(0, 0, 0);
  //   Position Vector of X=1 Y=1 Strip
  TVector3 Strip_1_1;

  //   Geometry Parameter
  double Base,Height;
  Base          = 100.00; // mm
  Height        = 100.00; // mm

  //double Face          = 98; // mm
  double NumberOfStrip = 128;
  double StripPitchBase    = Base / NumberOfStrip; // mm
  double StripPitchHeight  = Height / NumberOfStrip; // mm
  //   Buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  //   Moving StripCenter to 1.1 corner:
 // Strip_1_1 = C_X1_Y1 + U  * (StripPitchBase / 2.) + V * (StripPitchHeight / 2.);
  // This calculation recenter the strip around the detector center.
  // This account for cases where the provided corner coordinates
  // does not match the detector size
  TVector3 Center = 0.25*(C_X1_Y128 + C_X128_Y128 + C_X1_Y1 + C_X128_Y1);
  Strip_1_1 = Center-(0.5*Base*U+0.5*Height*V) + U  * (StripPitchBase / 2.) + V * (StripPitchHeight / 2.);

  for (int i = 0; i < 128; ++i) {
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      //StripCenter = Strip_1_1 + StripPitch * (i * U + j * V);
      StripCenter = Strip_1_1 + i*U*StripPitchBase  + j*V*StripPitchHeight;
      lineX.push_back(StripCenter.X());
      lineY.push_back(StripCenter.Y());
      lineZ.push_back(StripCenter.Z());
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);
}

void TMUSETTPhysics::InitializeStandardParameter() {
  //   Enable all channel
  vector<bool> ChannelStatus;
  m_XChannelStatus.clear();
  m_YChannelStatus.clear();

  ChannelStatus.resize(128, true);
  for (int i = 0; i < m_NumberOfDetector; ++i) {
    m_XChannelStatus[i] = ChannelStatus;
    m_YChannelStatus[i] = ChannelStatus;
  }


  m_MaximumStripMultiplicityAllowed = m_NumberOfDetector;

  return;
}

////////////////////////////////////////////////////////////////////////////////
void TMUSETTPhysics::AddDetector(double theta, double phi, double distance,
    double beta_u, double beta_v, double beta_w) {

  m_NumberOfDetector++;

  double Pi = 3.141592654;

  // convert from degree to radian:
  theta = theta * Pi / 180.;
  phi   = phi * Pi / 180.;

  // Vector U on Detector Face (paralelle to Y Strip) (NB: remember that Y
  // strip are allong X axis)
  TVector3 U;
  // Vector V on Detector Face (parallele to X Strip)
  TVector3 V;
  // Vector W normal to Detector Face (pointing CsI)
  TVector3 W;
  // Vector position of Detector Face center
  TVector3 C;

  C = TVector3(distance * sin(theta) * cos(phi),
      distance * sin(theta) * sin(phi), distance * cos(theta));

  TVector3 P
    = TVector3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));

  W = C.Unit();
  U = W.Cross(P);
  V = W.Cross(U);

  U = U.Unit();
  V = V.Unit();

  U.Rotate(beta_u * Pi / 180., U);
  V.Rotate(beta_u * Pi / 180., U);

  U.Rotate(beta_v * Pi / 180., V);
  V.Rotate(beta_v * Pi / 180., V);

  U.Rotate(beta_w * Pi / 180., W);
  V.Rotate(beta_w * Pi / 180., W);

  double Face          = 98; // mm
  double NumberOfStrip = 128;
  double StripPitch    = Face / NumberOfStrip; // mm

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  double X, Y, Z;

  // Moving C to the 1.1 corner:
  C.SetX(C.X() - (Face / 2 - StripPitch / 2) * (V.X() + U.X()));
  C.SetY(C.Y() - (Face / 2 - StripPitch / 2) * (V.Y() + U.Y()));
  C.SetZ(C.Z() - (Face / 2 - StripPitch / 2) * (V.Z() + U.Z()));

  for (int i = 0; i < 128; ++i) {

    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      X = C.X() + StripPitch * (U.X() * i + V.X() * j);
      Y = C.Y() + StripPitch * (U.Y() * i + V.Y() * j);
      Z = C.Z() + StripPitch * (U.Z() * i + V.Z() * j);

      lineX.push_back(X);
      lineY.push_back(Y);
      lineZ.push_back(Z);
    }

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneDetectorStripPositionX);
  m_StripPositionY.push_back(OneDetectorStripPositionY);
  m_StripPositionZ.push_back(OneDetectorStripPositionZ);
}

///////////////////////////////////////////////////////////////////////////////
TVector3 TMUSETTPhysics::GetPositionOfInteraction(const int i,bool random) {
  TVector3 Position
    = TVector3(GetStripPositionX(DetectorNumber[i], DSSD_X[i], DSSD_Y[i]),
        GetStripPositionY(DetectorNumber[i], DSSD_X[i], DSSD_Y[i]),
        GetStripPositionZ(DetectorNumber[i], DSSD_X[i], DSSD_Y[i]));
  return Position;
}

///////////////////////////////////////////////////////////////////////////////
/*double TMUSETTPhysics::GetRelativeAngle() {
  TVector3 Position0= TVector3(PosX[0],PosY[0],PosZ[0]);
  TVector3 Position1= TVector3(PosX[1],PosY[1],PosZ[1]);

  //std::cout << Position0.Angle(Position1)/(TMath::Pi())*180 << std::endl;
  //std::cout << DetectorNumber[0] << " " << DetectorNumber[0] << std::endl;
  //std::cout << DetectorNumber[1] << " " << DetectorNumber[1] << std::endl;

  //std::cout << PosX[0] << " " << PosY[0] << " " << PosZ[0] << std::endl;
  //std::cout << PosX[1] << " " << PosY[1] << " " << PosZ[1] << std::endl;
  //std::cout << std::endl;
  return Position0.Angle(Position1)/(TMath::Pi())*180;
}
*/
/*
double TMUSETTPhysics::GetS(const unsigned int A, const unsigned int Q2A) {
  double eta = DSSD_E[0]/DSSD_E[1];
  if(eta > 1) eta = 1./eta;
  double epsilon = 4./A;
  return Q2A/(1.+epsilon+ (2.*epsilon*pow(eta,1./2.)*cos(RelativeAngle*TMath::Pi()/180.))/(eta + 1.));

}
*/
///////////////////////////////////////////////////////////////////////////////
TVector3 TMUSETTPhysics::GetDetectorNormal(const int i) {
  TVector3 U = TVector3(GetStripPositionX(DetectorNumber[i], 128, 1),
      GetStripPositionY(DetectorNumber[i], 128, 1),
      GetStripPositionZ(DetectorNumber[i], 128, 1))

    - TVector3(GetStripPositionX(DetectorNumber[i], 1, 1),
        GetStripPositionY(DetectorNumber[i], 1, 1),
        GetStripPositionZ(DetectorNumber[i], 1, 1));

  TVector3 V = TVector3(GetStripPositionX(DetectorNumber[i], 128, 128),
      GetStripPositionY(DetectorNumber[i], 128, 128),
      GetStripPositionZ(DetectorNumber[i], 128, 128))

    - TVector3(GetStripPositionX(DetectorNumber[i], 128, 1),
        GetStripPositionY(DetectorNumber[i], 128, 1),
        GetStripPositionZ(DetectorNumber[i], 128, 1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}

///////////////////////////////////////////////////////////////////////////
namespace MUSETT_LOCAL {
  //   DSSD
  //   X
  double fDSSD_X_E(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDXEDetectorNbr(i));
    name += "_DSSD_X";
    name += NPL::itoa(m_EventData->GetDSSDXEStripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDXEEnergy(i),0);
  }

  double fDSSD_X_T(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDXTDetectorNbr(i));
    name += "_DSSD_X";
    name += NPL::itoa(m_EventData->GetDSSDXTStripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDXTTime(i),0);
  }

  //   Y
  double fDSSD_Y_E(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDYEDetectorNbr(i));
    name += "_DSSD_Y";
    name += NPL::itoa(m_EventData->GetDSSDYEStripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDYEEnergy(i),0);
  }

  double fDSSD_Y_T(const TMUSETTData* m_EventData, const int& i) {
    static string name;
    name = "MUSETT/T";
    name += NPL::itoa(m_EventData->GetDSSDYTDetectorNbr(i));
    name += "_DSSD_Y";
    name += NPL::itoa(m_EventData->GetDSSDYTStripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(
        name, m_EventData->GetDSSDYTTime(i),0);
  }
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TMUSETTPhysics::Construct() {
  return (NPL::VDetector*)new TMUSETTPhysics();
}

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
  class proxy_MUSETT {
    public:
      proxy_MUSETT() {
        NPL::DetectorFactory::getInstance()->AddToken("MUSETT", "MUSETT");
        NPL::DetectorFactory::getInstance()->AddDetector("MUSETT",
            TMUSETTPhysics::Construct);
      }
  };

  proxy_MUSETT p_MUSETT;
}
