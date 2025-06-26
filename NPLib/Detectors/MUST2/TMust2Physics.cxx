/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : febuary 2009                                             *
 * Last update    : July 2021
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold must2 treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "TMust2Physics.h"
using namespace MUST2_LOCAL;

//   STL
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdlib.h>

//   NPL
// #include "NPDetectorFactory.h"
#include "NPInputParser.h"
#include "NPOptionManager.h"
#include "NPSystemOfUnits.h"
#include "RootHistogramsCalib.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TAsciiFile.h"
#include "TF1.h"
#include "TSpectrum.h"
using namespace NPUNITS;
//   ROOT
#include "TFitResult.h"
#include "TChain.h"
///////////////////////////////////////////////////////////////////////////

ClassImp(TMust2Physics)

    ///////////////////////////////////////////////////////////////////////////
    TMust2Physics::TMust2Physics() {
  EventMultiplicity = 0;
  m_EventData = new TMust2Data;
  m_PreTreatedData = new TMust2Data;
  m_EventPhysics = this;
  m_Spectra = NULL;
  m_NumberOfTelescope = 0;
  m_MaximumStripMultiplicityAllowed = 10;
  m_StripEnergyMatchingSigma = 0.020;
  m_StripEnergyMatchingNumberOfSigma = 3;
  // Raw Threshold
  m_Si_X_E_RAW_Threshold = 8192;
  m_Si_Y_E_RAW_Threshold = 8192;
  m_SiLi_E_RAW_Threshold = 8192;
  m_CsI_E_RAW_Threshold = 8192;
  m_CsI_E_RAW_Threshold = 0;
  // Calibrated Threshold
  m_Si_X_E_Threshold = 0;
  m_Si_Y_E_Threshold = 0;
  m_SiLi_E_Threshold = 0;
  m_CsI_E_Threshold = 0;

  m_Ignore_not_matching_SiLi = false;
  m_Ignore_not_matching_CsI = false;

  m_Take_E_Y = false;
  m_Take_T_Y = true;

  //////////////////
  // SiLi matching //
  //////////////////

  m_SiLi_Size = 32;
  m_SiLi_MatchingX.resize(16, 0);
  m_SiLi_MatchingY.resize(16, 0);

  m_SiLi_MatchingX[0] = 112;
  m_SiLi_MatchingY[0] = 112;

  m_SiLi_MatchingX[1] = 112;
  m_SiLi_MatchingY[1] = 80;

  m_SiLi_MatchingX[2] = 80;
  m_SiLi_MatchingY[2] = 112;

  m_SiLi_MatchingX[3] = 80;
  m_SiLi_MatchingY[3] = 80;
  //
  m_SiLi_MatchingX[4] = 48;
  m_SiLi_MatchingY[4] = 80;

  m_SiLi_MatchingX[5] = 48;
  m_SiLi_MatchingY[5] = 112;

  m_SiLi_MatchingX[6] = 16;
  m_SiLi_MatchingY[6] = 80;

  m_SiLi_MatchingX[7] = 16;
  m_SiLi_MatchingY[7] = 112;
  //
  m_SiLi_MatchingX[8] = 112;
  m_SiLi_MatchingY[8] = 16;

  m_SiLi_MatchingX[9] = 112;
  m_SiLi_MatchingY[9] = 48;

  m_SiLi_MatchingX[10] = 80;
  m_SiLi_MatchingY[10] = 16;

  m_SiLi_MatchingX[11] = 80;
  m_SiLi_MatchingY[11] = 48;
  //
  m_SiLi_MatchingX[12] = 48;
  m_SiLi_MatchingY[12] = 48;

  m_SiLi_MatchingX[13] = 48;
  m_SiLi_MatchingY[13] = 16;

  m_SiLi_MatchingX[14] = 16;
  m_SiLi_MatchingY[14] = 48;

  m_SiLi_MatchingX[15] = 16;
  m_SiLi_MatchingY[15] = 16;

  //////////////////
  // CsI matching //
  //////////////////

  m_CsI_Size = 32;
  m_CsI_MatchingX.resize(16, 0);
  m_CsI_MatchingY.resize(16, 0);

  m_CsI_MatchingX[0] = 112;
  m_CsI_MatchingY[0] = 112;

  m_CsI_MatchingX[1] = 112;
  m_CsI_MatchingY[1] = 80;

  m_CsI_MatchingX[2] = 112;
  m_CsI_MatchingY[2] = 48;

  m_CsI_MatchingX[3] = 112;
  m_CsI_MatchingY[3] = 16;
  //
  m_CsI_MatchingX[4] = 80;
  m_CsI_MatchingY[4] = 16;

  m_CsI_MatchingX[5] = 80;
  m_CsI_MatchingY[5] = 48;

  m_CsI_MatchingX[6] = 80;
  m_CsI_MatchingY[6] = 80;

  m_CsI_MatchingX[7] = 80;
  m_CsI_MatchingY[7] = 112;
  //
  m_CsI_MatchingX[8] = 48;
  m_CsI_MatchingY[8] = 16;

  m_CsI_MatchingX[9] = 48;
  m_CsI_MatchingY[9] = 48;

  m_CsI_MatchingX[10] = 48;
  m_CsI_MatchingY[10] = 80;

  m_CsI_MatchingX[11] = 48;
  m_CsI_MatchingY[11] = 112;
  //
  m_CsI_MatchingX[12] = 16;
  m_CsI_MatchingY[12] = 16;

  m_CsI_MatchingX[13] = 16;
  m_CsI_MatchingY[13] = 48;

  m_CsI_MatchingX[14] = 16;
  m_CsI_MatchingY[14] = 80;

  m_CsI_MatchingX[15] = 16;
  m_CsI_MatchingY[15] = 112;


  // FIXME: Temporary fix... for zero degree
  // Particular strip matching required for zero degree telescope...
  // Might be a geometrical problem, only relevant if a telescope is placed at
  // zero degree
  // m_ZeroDegree_CsI_MatchingX[0].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[0].second = 128;
  // m_ZeroDegree_CsI_MatchingY[0].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[0].second = 128;

  // m_ZeroDegree_CsI_MatchingX[1].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[1].second = 128;
  // m_ZeroDegree_CsI_MatchingY[1].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[1].second = 102;

  // m_ZeroDegree_CsI_MatchingX[2].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[2].second = 128;
  // m_ZeroDegree_CsI_MatchingY[2].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[2].second = 64;

  // m_ZeroDegree_CsI_MatchingX[3].first  = 103;
  // m_ZeroDegree_CsI_MatchingX[3].second = 128;
  // m_ZeroDegree_CsI_MatchingY[3].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[3].second = 26;

  // m_ZeroDegree_CsI_MatchingX[4].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[4].second = 102;
  // m_ZeroDegree_CsI_MatchingY[4].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[4].second = 26;

  // m_ZeroDegree_CsI_MatchingX[5].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[5].second = 102;
  // m_ZeroDegree_CsI_MatchingY[5].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[5].second = 64;

  // m_ZeroDegree_CsI_MatchingX[6].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[6].second = 102;
  // m_ZeroDegree_CsI_MatchingY[6].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[6].second = 102;

  // m_ZeroDegree_CsI_MatchingX[7].first  = 65;
  // m_ZeroDegree_CsI_MatchingX[7].second = 102;
  // m_ZeroDegree_CsI_MatchingY[7].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[7].second = 128;

  // m_ZeroDegree_CsI_MatchingX[8].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[8].second = 64;
  // m_ZeroDegree_CsI_MatchingY[8].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[8].second = 26;

  // m_ZeroDegree_CsI_MatchingX[9].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[9].second = 64;
  // m_ZeroDegree_CsI_MatchingY[9].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[9].second = 64;

  // m_ZeroDegree_CsI_MatchingX[10].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[10].second = 64;
  // m_ZeroDegree_CsI_MatchingY[10].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[10].second = 102;

  // m_ZeroDegree_CsI_MatchingX[11].first  = 27;
  // m_ZeroDegree_CsI_MatchingX[11].second = 64;
  // m_ZeroDegree_CsI_MatchingY[11].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[11].second = 128;

  // m_ZeroDegree_CsI_MatchingX[12].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[12].second = 26;
  // m_ZeroDegree_CsI_MatchingY[12].first  = 1;
  // m_ZeroDegree_CsI_MatchingY[12].second = 26;

  // m_ZeroDegree_CsI_MatchingX[13].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[13].second = 26;
  // m_ZeroDegree_CsI_MatchingY[13].first  = 27;
  // m_ZeroDegree_CsI_MatchingY[13].second = 64;

  // m_ZeroDegree_CsI_MatchingX[14].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[14].second = 26;
  // m_ZeroDegree_CsI_MatchingY[14].first  = 65;
  // m_ZeroDegree_CsI_MatchingY[14].second = 102;

  // m_ZeroDegree_CsI_MatchingX[15].first  = 1;
  // m_ZeroDegree_CsI_MatchingX[15].second = 26;
  // m_ZeroDegree_CsI_MatchingY[15].first  = 103;
  // m_ZeroDegree_CsI_MatchingY[15].second = 128;
}

///////////////////////////////////////////////////////////////////////////
TMust2Physics::~TMust2Physics() {}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::BuildSimplePhysicalEvent() { BuildPhysicalEvent(); }

///////////////////////////////////////////////////////////////////////////

void TMust2Physics::BuildPhysicalEvent() {
  if (NPOptionManager::getInstance()->IsReader() == true) {
    m_EventData = &(**r_ReaderEventData);
  }
  PreTreat();
  
  std::vector<unsigned int> Raw_Si_X_E;
  std::vector<unsigned int> Raw_Si_Y_E;
  Raw_Si_X_E.clear();
  Raw_Si_Y_E.clear();

      // Treating again raw data to get raw Si energy in the data: should be removed at some point
      // Added because I did cut on wrong calib and now need raw Si to apply CUTS
  for (unsigned int i = 0; i < m_EventData->GetMMStripXEMult(); ++i) {
    if (m_EventData->GetMMStripXEEnergy(i) > m_Si_X_E_RAW_Threshold &&
        IsValidChannel(0, m_EventData->GetMMStripXEDetectorNbr(i), m_EventData->GetMMStripXEStripNbr(i))) {
      double EX = fSi_X_E(m_EventData, i);
      if (EX > m_Si_X_E_Threshold)
        Raw_Si_X_E.push_back(m_EventData->GetMMStripXEEnergy(i));
    }
  }
  
  for (unsigned int i = 0; i < m_EventData->GetMMStripYEMult(); ++i) {
    if (m_EventData->GetMMStripYEEnergy(i) < m_Si_Y_E_RAW_Threshold &&
        IsValidChannel(1, m_EventData->GetMMStripYEDetectorNbr(i), m_EventData->GetMMStripYEStripNbr(i))) {
      double EY = fSi_Y_E(m_EventData, i);
      if (EY > m_Si_Y_E_Threshold)
        Raw_Si_Y_E.push_back(m_EventData->GetMMStripYEEnergy(i));
    }
  }

  m_multimatch = false;

  bool check_SILI = false;
  bool check_CSI = false;

  m_StripXEMult = m_PreTreatedData->GetMMStripXEMult();
  m_StripYEMult = m_PreTreatedData->GetMMStripYEMult();
  m_StripXTMult = m_PreTreatedData->GetMMStripXTMult();
  m_StripYTMult = m_PreTreatedData->GetMMStripYTMult();
  m_SiLiEMult = m_PreTreatedData->GetMMSiLiEMult();
  m_SiLiTMult = m_PreTreatedData->GetMMSiLiTMult();
  m_CsIEMult = m_PreTreatedData->GetMMCsIEMult();
  m_CsITMult = m_PreTreatedData->GetMMCsITMult();

  // Just a sanity check, to be sure that mult is consistent
  if(m_StripXEMult != Raw_Si_X_E.size())
    std::cout << "Si_X_E mult not consistent" << std::endl;
  if(m_StripYEMult != Raw_Si_Y_E.size())
    std::cout << "Si_Y_E mult not consistent" << std::endl;

  // Returns the matching couples and set the type
  // of matching. In case of a multiple match in one
  // detector, m_match_type is set to 1, else it is
  // set to 2

  auto couple = Match_X_Y();
  auto match = Match_Si_CsI(couple);

      /////////////////////////////////////////////////
      for(auto event: match){
      unsigned int Xindex = event.second.first;
      unsigned int Yindex = event.second.second;
      int CsIindex = event.first;
      
      int N = m_PreTreatedData->GetMMStripXEDetectorNbr(Xindex);
      int X = m_PreTreatedData->GetMMStripXEStripNbr(Xindex);
      int Y = m_PreTreatedData->GetMMStripYEStripNbr(Yindex);
      int PixelN;
      if(CsIindex >= 0){
        CsI_N.push_back(m_PreTreatedData->GetMMCsIECristalNbr(CsIindex));
        CsI_E_Raw.push_back(m_PreTreatedData->GetMMCsIEEnergy(CsIindex));
        CsI_E.push_back(fCsI_E(m_PreTreatedData,CsIindex));
        CsI_T.push_back(-1000);
        PixelN = GetPixel(event).second;
      
         //Look for associate Time
        for (unsigned int k = 0; k < m_CsITMult; ++k) {
          // Same Cristal, Same Detector
          if (N == m_PreTreatedData->GetMMCsITDetectorNbr(k) && 
            m_PreTreatedData->GetMMCsIECristalNbr(CsIindex) == m_PreTreatedData->GetMMCsITCristalNbr(k)) {
            CsI_T[CsI_T.size() - 1] = m_PreTreatedData->GetMMCsITTime(CsIindex);
            break;}
        }
      }
      else{
        CsI_N.push_back(-1000);
        CsI_E_Raw.push_back(-1000);
        CsI_E.push_back(-1000);
        CsI_T.push_back(-1000);
        PixelN = -1;
      }
      double Si_X_T = -1000;
      for (unsigned int t = 0; t < m_StripXTMult; ++t) {
        if (N == m_PreTreatedData->GetMMStripXTDetectorNbr(t) && X == m_PreTreatedData->GetMMStripXTStripNbr(t)) {
          Si_X_T = m_PreTreatedData->GetMMStripXTTime(t);
          break;
        }
      }

      double Si_Y_T = -1000;
      for (unsigned int t = 0; t < m_StripYTMult; ++t) {
        if (N == m_PreTreatedData->GetMMStripYTDetectorNbr(t) && Y == m_PreTreatedData->GetMMStripYTStripNbr(t)) {
          Si_Y_T = m_PreTreatedData->GetMMStripYTTime(t);
          break;
        }
      }



       double Si_X_E = m_PreTreatedData->GetMMStripXEEnergy(Xindex);
       double Si_Y_E = m_PreTreatedData->GetMMStripYEEnergy(Yindex);
       double Si_X_E_Raw = Raw_Si_X_E[Xindex];
       double Si_Y_E_Raw = Raw_Si_Y_E[Yindex];

      TelescopeNumber.push_back(N);
      Si_X.push_back(X);
      Si_Y.push_back(Y);

      // Store the other value for checking purpose
      Si_EX.push_back(Si_X_E);
      Si_TX.push_back(Si_X_T);

      Si_EY.push_back(Si_Y_E);
      Si_TY.push_back(Si_Y_T);
      if(Cal_Pixel)
        Pixel.push_back(PixelN);

      if (m_Take_E_Y){
        Si_E.push_back(Si_Y_E);
        Si_E_Raw.push_back(Si_Y_E_Raw);
      }
      else{
        Si_E.push_back(Si_X_E);
        Si_E_Raw.push_back(Si_X_E_Raw);
      }
      if (m_Take_T_Y)
        Si_T.push_back(Si_Y_T);
      else
        Si_T.push_back(Si_X_T);

    }
  EventMultiplicity = TelescopeNumber.size();
  return;
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::PreTreat() {
  ClearPreTreatedData();
  m_StripXEMult = m_EventData->GetMMStripXEMult();
  m_StripYEMult = m_EventData->GetMMStripYEMult();
  m_StripXTMult = m_EventData->GetMMStripXTMult();
  m_StripYTMult = m_EventData->GetMMStripYTMult();
  m_SiLiEMult = m_EventData->GetMMSiLiEMult();
  m_SiLiTMult = m_EventData->GetMMSiLiTMult();
  m_CsIEMult = m_EventData->GetMMCsIEMult();
  m_CsITMult = m_EventData->GetMMCsITMult();
  map<int, int> hitX;
  map<int, int> hitY;

  //   X
  //   E
  for (unsigned int i = 0; i < m_StripXEMult; ++i) {
    if (m_EventData->GetMMStripXEEnergy(i) > m_Si_X_E_RAW_Threshold &&
        IsValidChannel(0, m_EventData->GetMMStripXEDetectorNbr(i), m_EventData->GetMMStripXEStripNbr(i))) {
      double EX = fSi_X_E(m_EventData, i);
      if (EX > m_Si_X_E_Threshold)
        m_PreTreatedData->SetStripXE(m_EventData->GetMMStripXEDetectorNbr(i), m_EventData->GetMMStripXEStripNbr(i), EX);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_StripXTMult; ++i) {
    if (IsValidChannel(0, m_EventData->GetMMStripXTDetectorNbr(i), m_EventData->GetMMStripXTStripNbr(i)))
      m_PreTreatedData->SetStripXT(m_EventData->GetMMStripXTDetectorNbr(i), m_EventData->GetMMStripXTStripNbr(i),
                                   fSi_X_T(m_EventData, i));
  }

  //   Y
  //   E
  for (unsigned int i = 0; i < m_StripYEMult; ++i) {
    if (m_EventData->GetMMStripYEEnergy(i) < m_Si_Y_E_RAW_Threshold &&
        IsValidChannel(1, m_EventData->GetMMStripYEDetectorNbr(i), m_EventData->GetMMStripYEStripNbr(i))) {
      double EY = fSi_Y_E(m_EventData, i);
      if (EY > m_Si_Y_E_Threshold)
        m_PreTreatedData->SetStripYE(m_EventData->GetMMStripYEDetectorNbr(i), m_EventData->GetMMStripYEStripNbr(i), EY);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_StripYTMult; ++i) {
    if (IsValidChannel(1, m_EventData->GetMMStripYTDetectorNbr(i), m_EventData->GetMMStripYTStripNbr(i)))
      m_PreTreatedData->SetStripYT(m_EventData->GetMMStripYTDetectorNbr(i), m_EventData->GetMMStripYTStripNbr(i),
                                   fSi_Y_T(m_EventData, i));
  }

  //   CsI
  //   E
  for (unsigned int i = 0; i < m_CsIEMult; ++i) {
    if (m_EventData->GetMMCsIEEnergy(i) > m_CsI_E_RAW_Threshold &&
        IsValidChannel(3, m_EventData->GetMMCsIEDetectorNbr(i), m_EventData->GetMMCsIECristalNbr(i))) {
      // Implementing special CSI E treatment: to get a calibrated and non calibrated branch in the analysis, the calibration is applied later,
      // but the threshold is still checked here
      double ECsI = fCsI_E(m_EventData, i);
      if (ECsI > m_CsI_E_Threshold) {
      m_PreTreatedData->SetCsIE(m_EventData->GetMMCsIEDetectorNbr(i), m_EventData->GetMMCsIECristalNbr(i), m_EventData->GetMMCsIEEnergy(i));
      }
    }
  }

  //   T
  for (unsigned int i = 0; i < m_CsITMult; ++i) {
    if (IsValidChannel(3, m_EventData->GetMMCsITDetectorNbr(i), m_EventData->GetMMCsITCristalNbr(i)))
      m_PreTreatedData->SetCsIT(m_EventData->GetMMCsITDetectorNbr(i), m_EventData->GetMMCsITCristalNbr(i),
                                fCsI_T(m_EventData, i));
  }

  //   SiLi
  //   E
  for (unsigned int i = 0; i < m_SiLiEMult; ++i) {
    if (m_EventData->GetMMSiLiEEnergy(i) > m_SiLi_E_RAW_Threshold &&
        IsValidChannel(2, m_EventData->GetMMSiLiEDetectorNbr(i), m_EventData->GetMMSiLiEPadNbr(i))) {
      double ESiLi = fSiLi_E(m_EventData, i);
      if (ESiLi > m_SiLi_E_Threshold)
        m_PreTreatedData->SetSiLiE(m_EventData->GetMMSiLiEDetectorNbr(i), m_EventData->GetMMSiLiEPadNbr(i), ESiLi);
    }
  }

  //   T
  for (unsigned int i = 0; i < m_SiLiTMult; ++i) {
    if (IsValidChannel(2, m_EventData->GetMMSiLiTDetectorNbr(i), m_EventData->GetMMSiLiTPadNbr(i)))
      m_PreTreatedData->SetSiLiT(m_EventData->GetMMSiLiTDetectorNbr(i), m_EventData->GetMMSiLiTPadNbr(i),
                                 fSiLi_T(m_EventData, i));
  }

  return;
}

///////////////////////////////////////////////////////////////////////////
bool TMust2Physics::ResolvePseudoEvent() { return false; }

///////////////////////////////////////////////////////////////////////////

std::vector<std::pair<unsigned int,unsigned int>> TMust2Physics::Match_X_Y() {
  // Simplified Match_X_Y to get all possible matching energy couples in same detector
  // The ony dubious case is 2 events in the same detector with the same energy, giving 4 possiblilities
  // If one of the matching couples finds a corresponding CsI in Matchi_Si_CsI, then the problem is solved
  // If none of the 4 couples finds a matching CsI, 2 Si only events are added to preserve good multiplicity
  // Thus, the multiplicity given by Match_X_Y is not a physical multiplicity, but a vector of all possible matches
  // The good physical multiplicity is given by Match_Si_CsI
  std::vector<std::pair<unsigned int,unsigned int>> ArrayOfGoodCouple;
  ArrayOfGoodCouple.clear();

  m_StripXEMult = m_PreTreatedData->GetMMStripXEMult();
  m_StripYEMult = m_PreTreatedData->GetMMStripYEMult();

  double matchSigma = m_StripEnergyMatchingSigma;
  double NmatchSigma = m_StripEnergyMatchingNumberOfSigma;

  // Prevent code from treating very high multiplicity Event
  // Those event are not physical anyway and that improve speed.
  if (m_StripXEMult > m_MaximumStripMultiplicityAllowed || m_StripYEMult > m_MaximumStripMultiplicityAllowed) {
    return ArrayOfGoodCouple;
  }

  for (unsigned int i = 0; i < m_StripXEMult; i++) {
    for (unsigned int j = 0; j < m_StripYEMult; j++) {

      // Declaration of variable for clarity
      int StripXDetNbr = m_PreTreatedData->GetMMStripXEDetectorNbr(i);
      int StripYDetNbr = m_PreTreatedData->GetMMStripYEDetectorNbr(j);

      //   if same detector check energy
      if (StripXDetNbr == StripYDetNbr) {

        int DetNbr = StripXDetNbr;

        // Declaration of variable for clarity
        double StripXEnergy = m_PreTreatedData->GetMMStripXEEnergy(i);

        double StripYEnergy = m_PreTreatedData->GetMMStripYEEnergy(j);

        //   Look if energy match
        if (abs((StripXEnergy - StripYEnergy) / 2.) < NmatchSigma * matchSigma) {

            ArrayOfGoodCouple.push_back(std::make_pair(i,j));

        }
      } // if same detector
    }   // loop on StripY Mult
  }     // loop on StripX Mult


  return ArrayOfGoodCouple;
}


////////////////////////////////////////////////////////////////////////////
bool TMust2Physics::IsValidChannel(const int& DetectorType, const int& telescope, const int& channel) {
  if (DetectorType == 0)
    return *(m_XChannelStatus[telescope - 1].begin() + channel - 1);

  else if (DetectorType == 1)
    return *(m_YChannelStatus[telescope - 1].begin() + channel - 1);

  else if (DetectorType == 2)
    return *(m_SiLiChannelStatus[telescope - 1].begin() + channel - 1);

  if (DetectorType == 3)
    return *(m_CsIChannelStatus[telescope - 1].begin() + channel - 1);

  else
    return false;
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigMust2.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigMust2.dat found: Default parameters loaded for "
            "Analysis "
         << FileName << endl;
    return;
  }
  cout << " Loading user parameters for Analysis from ConfigMust2.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigMust2.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer, DataBuffer, whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    if (LineBuffer.compare(0, 11, "ConfigMust2") == 0)
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus) {

      whatToDo = "";
      AnalysisConfigFile >> whatToDo;
      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n');
      }

      else if (whatToDo == "MAX_STRIP_MULTIPLICITY") {
        AnalysisConfigFile >> DataBuffer;
        m_MaximumStripMultiplicityAllowed = atoi(DataBuffer.c_str());
        cout << "MAXIMUN STRIP MULTIPLICITY " << m_MaximumStripMultiplicityAllowed << endl;
      }

      else if (whatToDo == "STRIP_ENERGY_MATCHING_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingSigma = atof(DataBuffer.c_str());
        cout << "STRIP ENERGY MATCHING SIGMA " << m_StripEnergyMatchingSigma << endl;
      }

      else if (whatToDo == "STRIP_ENERGY_MATCHING_NUMBER_OF_SIGMA") {
        AnalysisConfigFile >> DataBuffer;
        m_StripEnergyMatchingNumberOfSigma = atoi(DataBuffer.c_str());
        cout << "STRIP ENERGY MATCHING NUMBER OF SIGMA " << m_StripEnergyMatchingNumberOfSigma << endl;
      }

      else if (whatToDo == "DISABLE_ALL") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl;
        int telescope = atoi(DataBuffer.substr(2, 1).c_str());
        vector<bool> ChannelStatus;
        ChannelStatus.resize(128, false);
        m_XChannelStatus[telescope - 1] = ChannelStatus;
        m_YChannelStatus[telescope - 1] = ChannelStatus;
        ChannelStatus.resize(16, false);
        m_SiLiChannelStatus[telescope - 1] = ChannelStatus;
        m_CsIChannelStatus[telescope - 1] = ChannelStatus;
      }

      else if (whatToDo == "DISABLE_CHANNEL") {
        AnalysisConfigFile >> DataBuffer;
        cout << whatToDo << "  " << DataBuffer << endl;
        int telescope = atoi(DataBuffer.substr(2, 1).c_str());
        int channel = -1;
        if (DataBuffer.compare(3, 4, "STRX") == 0) {
          channel = atoi(DataBuffer.substr(7).c_str());
          *(m_XChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else if (DataBuffer.compare(3, 4, "STRY") == 0) {
          channel = atoi(DataBuffer.substr(7).c_str());
          *(m_YChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else if (DataBuffer.compare(3, 4, "SILI") == 0) {
          channel = atoi(DataBuffer.substr(7).c_str());
          *(m_SiLiChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else if (DataBuffer.compare(3, 3, "CSI") == 0) {
          channel = atoi(DataBuffer.substr(6).c_str());
          *(m_CsIChannelStatus[telescope - 1].begin() + channel - 1) = false;
        }

        else
          cout << "Warning: detector type for Must2 unknown!" << endl;
      }

      else if (whatToDo == "TAKE_E_Y") {
        m_Take_E_Y = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "TAKE_T_Y") {
        m_Take_T_Y = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "TAKE_E_X") {
        m_Take_E_Y = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "TAKE_T_X") {
        m_Take_T_Y = false;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "IGNORE_NOT_MATCHING_CSI") {
        m_Ignore_not_matching_CsI = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "CSI_SIZE") {
        AnalysisConfigFile >> DataBuffer;
        m_CsI_Size = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_CsI_Size << endl;
      }

      else if (whatToDo == "IGNORE_NOT_MATCHING_SILI") {
        m_Ignore_not_matching_SiLi = true;
        cout << whatToDo << endl;
      }

      else if (whatToDo == "SILI_SIZE") {
        AnalysisConfigFile >> DataBuffer;
        m_SiLi_Size = atoi(DataBuffer.c_str());
        cout << whatToDo << " " << m_SiLi_Size << endl;
      }

      else if (whatToDo == "SI_X_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_X_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_X_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "SI_Y_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_Y_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_Y_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "SILI_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_SiLi_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_SiLi_E_RAW_Threshold << endl;
      }

      else if (whatToDo == "CSI_E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_CsI_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_CsI_E_Threshold << endl;
      }

      else if (whatToDo == "SI_X_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_X_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_X_E_Threshold << endl;
      }

      else if (whatToDo == "SI_Y_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_Si_Y_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_Si_Y_E_Threshold << endl;
      }

      else if (whatToDo == "SILI_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_SiLi_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_SiLi_E_Threshold << endl;
      }

      else if (whatToDo == "CSI_E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_CsI_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_CsI_E_Threshold << endl;
      }
      else if (whatToDo == "CAL_PIXEL") {
        Cal_Pixel = true;
        cout << whatToDo << endl;
      }
      else if (whatToDo == "PIXEL_SIZE") {
        AnalysisConfigFile >> DataBuffer;
        PixelSize = (atoi(DataBuffer.c_str()));
        cout << "PIXEL SIZE " << PixelSize << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
bool TMust2Physics::Match_Si_SiLi(int X, int Y, int PadNbr) {

  // remove the central part and surrounding
  if (X < 8 || X > 120 || (Y < 68 && Y > 60))
    return false;

  if (abs(m_SiLi_MatchingX[PadNbr - 1] - X) < m_SiLi_Size / 2. &&
      abs(m_SiLi_MatchingY[PadNbr - 1] - Y) < m_SiLi_Size / 2.)
    return true;

  else
    return false;
}

///////////////////////////////////////////////////////////////////////////
bool TMust2Physics::Match_Si_CsI(int X, int Y, int CristalNbr, int TelescopeNumber) {
  // Legacy Match Si_CsI
  // FIXME: Temporary fix... for zero degree
  // Added to correct visible gaps at zero degree. This seem to be a geometrical
  // issue, this is not required in most cases (ZeroDegreeTelescopeNumber should
  // be added to TMust2Physics.h or to the input file if more experiment use a
  // telescope at zero degree
  // if (TelescopeNumber == ZeroDegreeTelescopeNumber) {
  //   if ((X >= m_ZeroDegree_CsI_MatchingX[CristalNbr - 1].first - 1
  //        && X <= m_ZeroDegree_CsI_MatchingX[CristalNbr - 1].second + 1)
  //       && (Y >= m_ZeroDegree_CsI_MatchingY[CristalNbr - 1].first - 1
  //           && Y <= m_ZeroDegree_CsI_MatchingY[CristalNbr - 1].second + 1)) {
  //     return true;
  //   } else {
  //     return false;
  //   }
  // } else {
  if (abs((double)m_CsI_MatchingX[CristalNbr - 1] +0.5 - X) <= (double)m_CsI_Size / 2. &&
      abs((double)m_CsI_MatchingY[CristalNbr - 1] +0.5 - Y) <= (double)m_CsI_Size / 2.) {
    return true;
  }
  else {
    return false;
  }
  // }
}

std::vector<std::pair<int, std::pair<unsigned int, unsigned int>>> 
TMust2Physics::Match_Si_CsI(std::vector<std::pair<unsigned int, unsigned int>> Array_of_good_XY_couples) {
  // New Match Si_CsI->Takes into account possible events with an overlap in 2 CsI
  // vector of matching pairs: first is CsI index, second is the pair of XY indexes
  std::vector<std::pair<int, std::pair<unsigned int, unsigned int>>> Array_of_matches;
  Array_of_matches.reserve(Array_of_good_XY_couples.size());
  std::vector<std::vector<std::pair<int, std::pair<unsigned int,unsigned int>>>::iterator> Wrong_matches;
  for(unsigned int i = 0; i < m_PreTreatedData->GetMMCsIEMult(); i++){
    for(auto couple: Array_of_good_XY_couples){
      unsigned int Xindex = couple.first;
      unsigned int Yindex = couple.second;

      // a couple is already matched in detector 
      unsigned int XYdet = m_PreTreatedData->GetMMStripXEDetectorNbr(Xindex);
      unsigned int Xstrip = m_PreTreatedData->GetMMStripXEStripNbr(Xindex);
      unsigned int Ystrip = m_PreTreatedData->GetMMStripYEStripNbr(Yindex);
      
      unsigned int CsIdet = m_PreTreatedData->GetMMCsIEDetectorNbr(i);
      
      // Checking if detnbr matches
      if(XYdet == CsIdet){
        unsigned int CsINbr = m_PreTreatedData->GetMMCsIECristalNbr(i);
        
        //Checking if the projection of XY position on CsI matches
        if(abs((double)m_CsI_MatchingX[CsINbr - 1] + 0.5 - (double)Xstrip) <= (double)m_CsI_Size / 2. &&
        abs((double)m_CsI_MatchingY[CsINbr - 1] + 0.5 - (double)Ystrip) <= (double)m_CsI_Size / 2.){
          
          // It is a priori possible that 2 XY pairs match with the same CsI (in particular with large CsI sizes)
          // We need to check if 2 XY pairs match with the same CsI, and we then reject the 2 XY pairs AND the CsI event 
        
          // iterator that points to an element of Array_of_matches if, in the same detector, this CsI or this couple has already been matched. Else it points to end. 
          auto it = std::find_if(Array_of_matches.begin(),Array_of_matches.end(),
          [i,couple](const std::pair<int,std::pair<unsigned int,unsigned int>>& match){
          return (match.first == i || match.second.first == couple.first || match.second.second == couple.second);});
          
          // if this CsI has not been found, we can add the CsI XY match to Array_of_matches
          if(it == Array_of_matches.end()){
            Array_of_matches.push_back(std::make_pair(i,couple));
          }
          // if this CsI has been found, the previous match is kept but the associated iterator is saved in a vector to be deleted at the end of the loop
          else{
          // Creating a vector of wrong matches to delete them from the vector after looping over all events
            if(std::find(Wrong_matches.begin(), Wrong_matches.end(),it)==Wrong_matches.end())
              Wrong_matches.push_back(it);
          }
        
        }  
      }
    }
  }
  for(auto it: Wrong_matches){
      Array_of_matches.erase(it);
  }
  // However, when we deleted all these wrong matches, we removed events that we cant match well to a CsI (double matches) or events that stop in the Si
  // Removing all these events will make false multiplicities. We add these events again with index -1 to indicate that these events are 
  // physical events in DSSD that we cant match to a CsI
  for(auto couple: Array_of_good_XY_couples){
    // iterator that points to an element of Array_of_matches if an X Y couple has been matched, else it points to end 
    auto it = std::find_if(Array_of_matches.begin(),Array_of_matches.end(),
    [couple](const std::pair<int,std::pair<unsigned int,unsigned int>>& match){
    return (match.second.first == couple.first || match.second.second == couple.second);});
    
    if(it == Array_of_matches.end()){
      Array_of_matches.push_back(std::make_pair(-1,couple));
    }
  }

  return Array_of_matches;
}

std::pair<unsigned int, unsigned int> TMust2Physics::GetPixel(std::pair<int,std::pair <unsigned int, unsigned int>> match){
  if(m_CsI_Size%PixelSize != 0)
    std::cout << "WARNING: Your current Pixel Size: " << PixelSize << " does not divide your CsI size: " << m_CsI_Size << ". The actual size of the pixels may be dubious." << std::endl;
  unsigned int CsINbr = m_PreTreatedData->GetMMCsIECristalNbr(match.first);
  double CsIPosX = m_CsI_MatchingX[CsINbr-1] + 0.5;
  double CsIPosY = m_CsI_MatchingY[CsINbr-1] + 0.5;
  double StripX = m_PreTreatedData->GetMMStripXEStripNbr(match.second.first);
  double StripY = m_PreTreatedData->GetMMStripYEStripNbr(match.second.second);
  double FirstPixelX = CsIPosX - m_CsI_Size/2.;
  double FirstPixelY = CsIPosY - m_CsI_Size/2.;
  unsigned int PixelNumber = (int)((StripX - FirstPixelX)/((double)PixelSize)) + (m_CsI_Size/PixelSize)*(int)((StripY - FirstPixelY)/((double)PixelSize));
  return std::make_pair(CsINbr,PixelNumber);
}


///////////////////////////////////////////////////////////////////////////
void TMust2Physics::Clear() {
  EventMultiplicity = 0;

  m_match_type.clear();

  m_StripXMultDet.clear();
  m_StripYMultDet.clear();
  m_NMatchDet.clear();
  m_NMatchX.clear();
  m_NMatchY.clear();

  TelescopeNumber.clear();
  EventType.clear();
  TotalEnergy.clear();

  // Si X
  Si_E.clear();
  Si_E_Raw.clear();
  Si_T.clear();
  Si_X.clear();
  Si_Y.clear();

  // Si(Li)
  SiLi_E.clear();
  SiLi_T.clear();
  SiLi_N.clear();

  // CsI
  CsI_E.clear();
  CsI_E_Raw.clear();
  CsI_T.clear();
  CsI_N.clear();

  Si_EX.clear();
  Si_TX.clear();
  Si_EY.clear();
  Si_TY.clear();
  TelescopeNumber_X.clear();
  TelescopeNumber_Y.clear();
  Pixel.clear();
}
///////////////////////////////////////////////////////////////////////////

void TMust2Physics::ReadCalibrationRun() {
  m_StripXEMult = m_EventData->GetMMStripXEMult();
  m_StripYEMult = m_EventData->GetMMStripYEMult();
  m_StripXTMult = m_EventData->GetMMStripXTMult();
  m_StripYTMult = m_EventData->GetMMStripYTMult();
  m_SiLiEMult = m_EventData->GetMMSiLiEMult();
  m_SiLiTMult = m_EventData->GetMMSiLiTMult();
  m_CsIEMult = m_EventData->GetMMCsIEMult();
  m_CsITMult = m_EventData->GetMMCsITMult();

  //   X
  //   E
  for (unsigned int i = 0; i < m_StripXEMult; ++i) {
    TelescopeNumber_X.push_back(m_EventData->GetMMStripXEDetectorNbr(i));
    Si_EX.push_back(fSi_X_E(m_EventData, i));
    Si_X.push_back(m_EventData->GetMMStripXEStripNbr(i));
  }
  //   T
  for (unsigned int i = 0; i < m_StripXTMult; ++i) {
    TelescopeNumber_X.push_back(m_EventData->GetMMStripXTDetectorNbr(i));
    Si_TX.push_back(fSi_X_T(m_EventData, i));
    Si_X.push_back(m_EventData->GetMMStripXTStripNbr(i));
  }

  //   Y
  //   E
  for (unsigned int i = 0; i < m_StripYEMult; ++i) {
    TelescopeNumber_Y.push_back(m_EventData->GetMMStripYEDetectorNbr(i));
    Si_EY.push_back(fSi_Y_E(m_EventData, i));
    Si_Y.push_back(m_EventData->GetMMStripYEStripNbr(i));
  }

  //   T
  for (unsigned int i = 0; i < m_StripYTMult; ++i) {
    TelescopeNumber_Y.push_back(m_EventData->GetMMStripYTDetectorNbr(i));
    Si_TY.push_back(fSi_Y_T(m_EventData, i));
    Si_Y.push_back(m_EventData->GetMMStripYTStripNbr(i));
  }

  // CsI Energy
  for (unsigned int j = 0; j < m_CsIEMult; ++j) {
    CsI_N.push_back(m_EventData->GetMMCsIECristalNbr(j));
    CsI_E.push_back(fCsI_E(m_EventData, j));
  }

  // CsI Time
  for (unsigned int j = 0; j < m_CsITMult; ++j) {
    CsI_T.push_back(fCsI_T(m_EventData, j));
  }
}

////   Innherited from VDetector Class   ////

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::ReadConfiguration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("M2Telescope");
  if (NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks.size() << " Telescope found " << endl;

  // Cartesian Case
  vector<string> cart = {"X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128", "SI", "SILI", "CSI"};
  // Spherical Case
  vector<string> sphe = {"R", "THETA", "PHI", "BETA", "SI", "SILI", "CSI"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(cart)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Must2 Telescope " << i + 1 << endl;
      TVector3 A = blocks[i]->GetTVector3("X1_Y1", "mm");
      TVector3 B = blocks[i]->GetTVector3("X128_Y1", "mm");
      TVector3 C = blocks[i]->GetTVector3("X1_Y128", "mm");
      TVector3 D = blocks[i]->GetTVector3("X128_Y128", "mm");
      AddTelescope(A, B, C, D);
      m_CsIPresent[i + 1] = blocks[i]->GetInt("CSI");
      m_SiLiPresent[i + 1] = blocks[i]->GetInt("SILI");
      m_CsIOffset[i + 1] = blocks[i]->GetInt("CsIOffset");
    }

    else if (blocks[i]->HasTokenList(sphe)) {
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Must2 Telescope " << i + 1 << endl;

      double Theta = blocks[i]->GetDouble("THETA", "deg");
      double Phi = blocks[i]->GetDouble("PHI", "deg");
      double R = blocks[i]->GetDouble("R", "mm");
      vector<double> beta = blocks[i]->GetVectorDouble("BETA", "deg");
      cout << Theta << " " << Phi << endl;
      AddTelescope(Theta, Phi, R, beta[0], beta[1], beta[2]);

      m_CsIPresent[i + 1] = blocks[i]->GetInt("CSI");
      m_SiLiPresent[i + 1] = blocks[i]->GetInt("SILI");
      m_CsIOffset[i + 1] = blocks[i]->GetInt("CsIOffset");
    }

    else {
      cout << "ERROR: Missing token for M2Telescope blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
  // m_MultEvt = new TMust2MultTelescope[m_NumberOfTelescope];

  InitializeStandardParameter();
  ReadAnalysisConfig();
}

void TMust2Physics::ReadDoCalibration(NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("M2Telescope");
  vector<NPL::InputBlock*> Energyblocks = parser.GetAllBlocksWithToken("EnergyCalibrationParameters");
  vector<NPL::InputBlock*> CSIblocks = parser.GetAllBlocksWithToken("CSICalibrationParameters");
  // if (NPOptionManager::getInstance()->GetVerboseLevel())
  // cout << "//// " << blocks.size() << " Telescope found " << endl;

  vector<string> calibs = {"TelescopeNumber", "Time", "Energy", "CSI"};
  vector<string> EnergyParameters = {"TelescopeNumber", "XThreshold", "YThreshold", "AlphaFitType"};
  vector<string> CSIParameters = {"TelescopeNumber", "CsIEnergyXThreshold", "CsIEnergyYThreshold", "CSIEThreshold","SiThickness","AlThickness","X1_Y1", "X1_Y128", "X128_Y1", "X128_Y128"};

  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(calibs)) {
      unsigned int TelescopeNumber = blocks[i]->GetInt("TelescopeNumber");
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Must2 Telescope " << TelescopeNumber << endl;
      DoCalibrationTime[TelescopeNumber] = blocks[i]->GetInt("Time");
      DoCalibrationEnergy[TelescopeNumber] = blocks[i]->GetInt("Energy");
      DoCalibrationCsI[TelescopeNumber] = blocks[i]->GetInt("CSI");
      if (DoCalibrationEnergy[TelescopeNumber] == 1) {
        IsCalibEnergy = true;
        if (EnergyXThreshold.count(TelescopeNumber) == 0) {
          EnergyXThreshold[TelescopeNumber] = 8192;
        }
        if (EnergyYThreshold.count(TelescopeNumber) == 0) {
          EnergyYThreshold[TelescopeNumber] = 8192;
        }
        if (AlphaFitType.count(TelescopeNumber) == 0) {
          AlphaFitType[TelescopeNumber] = "NoSatellite";
        }
      }
      if (DoCalibrationCsI[TelescopeNumber] == 1) {
        IsCalibCSI = true;
        if (CSIEnergyXThreshold.count(TelescopeNumber) == 0) {
          CSIEnergyXThreshold[TelescopeNumber] = 0;
        }
        if (CSIEnergyYThreshold.count(TelescopeNumber) == 0) {
          CSIEnergyYThreshold[TelescopeNumber] = 0;
        }
        if (CSIEThreshold.count(TelescopeNumber) == 0) {
          CSIEThreshold[TelescopeNumber] = 8192;
        }
      }
    }
    else {
      cout << "ERROR: Missing token for M2Telescope DoCalibration blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
  for (unsigned int i = 0; i < Energyblocks.size(); i++) {
    if (Energyblocks[i]->HasTokenList(EnergyParameters)) {
      unsigned int TelescopeNumber = Energyblocks[i]->GetInt("TelescopeNumber");
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  Energy Calibration parameters for MUST2 Telescope " << TelescopeNumber << endl;
      EnergyXThreshold[TelescopeNumber] = Energyblocks[i]->GetInt("XThreshold");
      EnergyYThreshold[TelescopeNumber] = Energyblocks[i]->GetInt("YThreshold");
      AlphaFitType[TelescopeNumber] = Energyblocks[i]->GetString("AlphaFitType");
      m_NumberOfTelescope++;
    }
    else {
      cout << "ERROR: Missing token for EnergyParam DoCalibration blocks, check your "
              "input "
              "file"
           << endl;
      exit(1);
    }
  }
  for (unsigned int i = 0; i < CSIblocks.size(); i++) {
    if (CSIblocks[i]->HasTokenList(CSIParameters)) {
      unsigned int TelescopeNumber = CSIblocks[i]->GetInt("TelescopeNumber");
      if (NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  CSI Calibration parameters for MUST2 Telescope " << TelescopeNumber << endl;
      PixelSize = CSIblocks[i]->GetInt("PixelSize"); 
      CSIEnergyXThreshold[TelescopeNumber] = CSIblocks[i]->GetInt("CsIEnergyXThreshold");
      CSIEnergyYThreshold[TelescopeNumber] = CSIblocks[i]->GetInt("CsIEnergyYThreshold");
      CSIEThreshold[TelescopeNumber] = CSIblocks[i]->GetInt("CSIEThreshold");
      SiThickness[TelescopeNumber] = CSIblocks[i]->GetDouble("SiThickness","um"); 
      AlThickness[TelescopeNumber] = CSIblocks[i]->GetDouble("AlThickness","um"); 
      TVector3 A = CSIblocks[i]->GetTVector3("X1_Y1", "mm");
      TVector3 B = CSIblocks[i]->GetTVector3("X128_Y1", "mm");
      TVector3 C = CSIblocks[i]->GetTVector3("X1_Y128", "mm");
      TVector3 D = CSIblocks[i]->GetTVector3("X128_Y128", "mm");
      AddTelescope(A, B, C, D);
    }
    else {
      cout << "ERROR: Missing token for CSI DoCalibration blocks, check your "
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
void TMust2Physics::InitSpectra() { m_Spectra = new TMust2Spectra(m_NumberOfTelescope); }

void TMust2Physics::SetTreeReader(TTreeReader* TreeReader) { TMust2PhysicsReader::r_SetTreeReader(TreeReader); }

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::FillSpectra() {
  m_Spectra->FillRawSpectra(m_EventData);
  m_Spectra->FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra->FillPhysicsSpectra(m_EventPhysics);
}
///////////////////////////////////////////////////////////////////////////
void TMust2Physics::CheckSpectra() { m_Spectra->CheckSpectra(); }
///////////////////////////////////////////////////////////////////////////
void TMust2Physics::ClearSpectra() {
  // To be done
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::WriteSpectra() {
  if (m_Spectra)
    m_Spectra->WriteSpectra();
}

///////////////////////////////////////////////////////////////////////////
map<string, TH1*> TMust2Physics::GetSpectra() {
  if (m_Spectra)
    return m_Spectra->GetMapHisto();
  else {
    map<string, TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
  // Good for simulation, close to typical values
  vector<double> standardX = {-63, 63. / 8192.};
  vector<double> standardY = {63, -63. / 8192.};
  vector<double> standardSiLi = {-250, 250. / 8192.};
  vector<double> standardT = {-1000, 1000. / 8192.};

  vector<double> standardCsI;
  for (int i = 0; i < m_NumberOfTelescope; ++i) {

    if (m_CsIOffset[i+1] == 1) {
      standardCsI = {0, 500. / 16384.};
    }
    else{
      standardCsI = {-250, 250. / 8192.};
    }

    for (int j = 0; j < 128; ++j) {
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_E",
                        "MUST2_T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_E", standardX);
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_E",
                        "MUST2_T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_E", standardY);
      //FIXME This line should be removed, only here because I made a mistake while doing my CUTS for CSI calib
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_Si_X_O" + NPL::itoa(j + 1) + "_E",
                        "MUST2_T" + NPL::itoa(i + 1) + "_Si_X_O" + NPL::itoa(j + 1) + "_E", standardX);
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_Si_Y_O" + NPL::itoa(j + 1) + "_E",
                        "MUST2_T" + NPL::itoa(i + 1) + "_Si_Y_O" + NPL::itoa(j + 1) + "_E", standardY);
      
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_T",
                        "MUST2_T" + NPL::itoa(i + 1) + "_Si_X" + NPL::itoa(j + 1) + "_T", standardT);
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_T",
                        "MUST2_T" + NPL::itoa(i + 1) + "_Si_Y" + NPL::itoa(j + 1) + "_T", standardT);
    }

    for (int j = 0; j < 16; ++j) {
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_E",
                        "MUST2_T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_E", standardSiLi);
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_T",
                        "MUST2_T" + NPL::itoa(i + 1) + "_SiLi" + NPL::itoa(j + 1) + "_T", standardT);
    }

    for (int j = 0; j < 16; ++j) {
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
                        "MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E", standardCsI);
      Cal->AddParameter("MUST2", "T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
                        "MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T", standardT);
      for(auto p: ParticleType){
        Cal->AddParameter("MUST2", p+"_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
                          p+"_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E", standardCsI);
        Cal->AddParameter("MUST2", p+"_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
                          p+"_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T", standardT);
        if(Cal_Pixel && std::find(ParticleTypePixel.begin(),ParticleTypePixel.end(),p)!= ParticleTypePixel.end()){
          unsigned int PixelNumber = (m_CsI_Size/PixelSize)*(m_CsI_Size/PixelSize);
          for(unsigned int k = 0; k < PixelNumber;k++ ){
            Cal->AddParameter("MUST2", p+"_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_Pixel" + NPL::itoa(k) + "_E",
                              p+"_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_Pixel" + NPL::itoa(k) + "_E", standardCsI);
            Cal->AddParameter("MUST2", p+"_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_Pixel" + NPL::itoa(k) + "_T",
                              p+"_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_Pixel" + NPL::itoa(k) + "_T", standardT);
          }
        }
      }


      //Cal->AddParameter("MUST2", "deuteron_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
      //                  "deuteron_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E", standardCsI);
      //Cal->AddParameter("MUST2", "deuteron_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
      //                  "deuteron_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T", standardT);

      //Cal->AddParameter("MUST2", "triton_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
      //                  "triton_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E", standardCsI);
      //Cal->AddParameter("MUST2", "triton_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
      //                  "triton_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T", standardT);

      //Cal->AddParameter("MUST2", "alpha_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E",
      //                  "alpha_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_E", standardCsI);
      //Cal->AddParameter("MUST2", "alpha_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T",
      //                  "alpha_MUST2_T" + NPL::itoa(i + 1) + "_CsI" + NPL::itoa(j + 1) + "_T", standardT);
    }
  }

  return;
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else {
    inputChain->SetBranchStatus("MUST2", true);
    inputChain->SetBranchStatus("fMM_*", true);
    inputChain->SetBranchAddress("MUST2", &m_EventData);
  }
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  // Option to use the nptreereader anaysis
  if (NPOptionManager::getInstance()->IsReader() == true) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    inputTreeReader->SetTree(inputChain);
  }
  // Option to use the standard npanalysis
  else {
    inputChain->SetBranchStatus("MUST2", true);
    inputChain->SetBranchStatus("EventMultiplicity", true);
    inputChain->SetBranchStatus("EventType", true);
    inputChain->SetBranchStatus("TelescopeNumber", true);
    inputChain->SetBranchStatus("Si_E", true);
    inputChain->SetBranchStatus("Si_T", true);
    inputChain->SetBranchStatus("Si_X", true);
    inputChain->SetBranchStatus("Si_Y", true);
    inputChain->SetBranchStatus("Si_EX", true);
    inputChain->SetBranchStatus("Si_TX", true);
    inputChain->SetBranchStatus("Si_EY", true);
    inputChain->SetBranchStatus("Si_TY", true);
    inputChain->SetBranchStatus("TelescopeNumber_X", true);
    inputChain->SetBranchStatus("TelescopeNumber_Y", true);
    inputChain->SetBranchStatus("SiLi_E", true);
    inputChain->SetBranchStatus("SiLi_T", true);
    inputChain->SetBranchStatus("SiLi_N", true);
    inputChain->SetBranchStatus("CsI_E", true);
    inputChain->SetBranchStatus("CsI_E_Raw", true);
    inputChain->SetBranchStatus("CsI_T", true);
    inputChain->SetBranchStatus("CsI_N", true);
    inputChain->SetBranchStatus("TotalEnergy", true);
    inputChain->SetBranchStatus("Pixel", true);
    inputChain->SetBranchAddress("MUST2", &m_EventPhysics);
  }
}

///////////////////////////////////////////////////////////////////////////
void TMust2Physics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("MUST2", "TMust2Physics", &m_EventPhysics);
}

/////   Specific to MUST2Array   ////

void TMust2Physics::AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1, TVector3 C_X1_Y128, TVector3 C_X128_Y128) {
  // To avoid warning
  C_X128_Y128 *= 1;

  m_NumberOfTelescope++;

  //   Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that
  //   Y strip are allong X axis)
  TVector3 U = C_X128_Y1 - C_X1_Y1;
  double Ushift = (U.Mag() - 98) / 2.;
  U = U.Unit();
  //   Vector V on Telescope Face (parallele to X Strip)
  TVector3 V = C_X1_Y128 - C_X1_Y1;
  double Vshift = (V.Mag() - 98) / 2.;
  V = V.Unit();

  //   Position Vector of Strip Center
  TVector3 StripCenter = TVector3(0, 0, 0);
  //   Position Vector of X=1 Y=1 Strip
  TVector3 Strip_1_1;

  //   Geometry Parameter
  double Face = 98; // mm
  double NumberOfStrip = 128;
  double StripPitch = Face / NumberOfStrip; // mm
  //   Buffer object to fill Position Array
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

  //   Moving StripCenter to 1.1 corner:
  Strip_1_1 = C_X1_Y1 + (U + V) * (StripPitch / 2.);
  Strip_1_1 += U * Ushift + V * Vshift;

  for (int i = 0; i < 128; ++i) {
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    for (int j = 0; j < 128; ++j) {
      StripCenter = Strip_1_1 + StripPitch * (i * U + j * V);
      lineX.push_back(StripCenter.X());
      lineY.push_back(StripCenter.Y());
      lineZ.push_back(StripCenter.Z());
    }

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}

void TMust2Physics::InitializeStandardParameter() {
  //   Enable all channel
  vector<bool> ChannelStatus;
  m_XChannelStatus.clear();
  m_YChannelStatus.clear();
  m_SiLiChannelStatus.clear();
  m_CsIChannelStatus.clear();

  ChannelStatus.resize(128, true);
  for (int i = 0; i < m_NumberOfTelescope; ++i) {
    m_XChannelStatus[i] = ChannelStatus;
    m_YChannelStatus[i] = ChannelStatus;
  }

  ChannelStatus.resize(16, true);
  for (int i = 0; i < m_NumberOfTelescope; ++i) {
    m_SiLiChannelStatus[i] = ChannelStatus;
    m_CsIChannelStatus[i] = ChannelStatus;
  }

  m_MaximumStripMultiplicityAllowed = m_NumberOfTelescope;

  return;
}

void TMust2Physics::AddTelescope(double theta, double phi, double distance, double beta_u, double beta_v,
                                 double beta_w) {

  m_NumberOfTelescope++;

  double Pi = 3.141592654;

  cout << theta << " " << phi << endl;
  // convert from degree to radian:
  // theta = theta * Pi / 180.;
  // phi   = phi * Pi / 180.;

  // Vector U on Telescope Face (paralelle to Y Strip) (NB: remember that Y
  // strip are allong X axis)
  TVector3 U;
  // Vector V on Telescope Face (parallele to X Strip)
  TVector3 V;
  // Vector W normal to Telescope Face (pointing CsI)
  TVector3 W;
  // Vector position of Telescope Face center
  TVector3 C;

  C = TVector3(distance * sin(theta) * cos(phi), distance * sin(theta) * sin(phi), distance * cos(theta));

  //   std::cout << C.Theta() * 180. / Pi << " " << C.Phi() * 180. / Pi << " "
  //             << C.Mag() << std::endl;

  std::cout << C.X() << " " << C.Y() << std::endl;

  TVector3 P = TVector3(cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta));
  // = TVector3(1, 1, 1);

  W = C.Unit();
  // U = W.Cross(P);
  V = W.Cross(P);
  // U = TVector3(0, -1, 0);
  // V = W.Cross(U);
  U = W.Cross(V);
  // V = TVector3(-1, 0, 0);

  U = U.Unit();
  // U.Rotate(90. * Pi / 180.,W);
  V = V.Unit();
  // V.Rotate(90. * Pi / 180.,W);

  std::cout << "U: " << U.X() << " " << U.Y() << " " << U.Z() << " "
            << "V: " << V.X() << " " << V.Y() << " " << V.Z() << " "
            << "W: " << W.X() << " " << W.Y() << " " << W.Z() << std::endl;

  // exit(1);

  U.Rotate(beta_u * Pi / 180., U);
  V.Rotate(beta_u * Pi / 180., U);

  U.Rotate(beta_v * Pi / 180., V);
  V.Rotate(beta_v * Pi / 180., V);

  U.Rotate(beta_w * Pi / 180., W);
  V.Rotate(beta_w * Pi / 180., W);

  double Face = 98; // mm
  double NumberOfStrip = 128;
  double StripPitch = Face / NumberOfStrip; // mm

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneTelescopeStripPositionX;
  vector<vector<double>> OneTelescopeStripPositionY;
  vector<vector<double>> OneTelescopeStripPositionZ;

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

    OneTelescopeStripPositionX.push_back(lineX);
    OneTelescopeStripPositionY.push_back(lineY);
    OneTelescopeStripPositionZ.push_back(lineZ);
  }
  m_StripPositionX.push_back(OneTelescopeStripPositionX);
  m_StripPositionY.push_back(OneTelescopeStripPositionY);
  m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
}

TVector3 TMust2Physics::GetPositionOfInteraction(const int i) const {
  TVector3 Position = TVector3(GetStripPositionX(TelescopeNumber[i], Si_X[i], Si_Y[i]),
                               GetStripPositionY(TelescopeNumber[i], Si_X[i], Si_Y[i]),
                               GetStripPositionZ(TelescopeNumber[i], Si_X[i], Si_Y[i]));

  return (Position);
}

TVector3 TMust2Physics::GetTelescopeNormal(const int i) const {
  TVector3 U = TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1), GetStripPositionY(TelescopeNumber[i], 128, 1),
                        GetStripPositionZ(TelescopeNumber[i], 128, 1))

               - TVector3(GetStripPositionX(TelescopeNumber[i], 1, 1), GetStripPositionY(TelescopeNumber[i], 1, 1),
                          GetStripPositionZ(TelescopeNumber[i], 1, 1));

  TVector3 V =
      TVector3(GetStripPositionX(TelescopeNumber[i], 128, 128), GetStripPositionY(TelescopeNumber[i], 128, 128),
               GetStripPositionZ(TelescopeNumber[i], 128, 128))

      - TVector3(GetStripPositionX(TelescopeNumber[i], 128, 1), GetStripPositionY(TelescopeNumber[i], 128, 1),
                 GetStripPositionZ(TelescopeNumber[i], 128, 1));

  TVector3 Normal = U.Cross(V);

  return (Normal.Unit());
}
void TMust2Physics::InitializeRootHistogramsCalib() {
  std::cout << "Initialize Histograms" << std::endl;
  map<int, bool>::iterator it;
  for (it = DoCalibrationEnergy.begin(); it != DoCalibrationEnergy.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsEnergyF(it->first);
    }
  }
  for (it = DoCalibrationTime.begin(); it != DoCalibrationTime.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsTimeF(it->first);
    }
  }
  for (it = DoCalibrationCsI.begin(); it != DoCalibrationCsI.end(); it++) {
    if (it->second) {
      InitializeRootHistogramsCSIF(it->first);
    }
  }
  if (NPOptionManager::getInstance()->IsReader() == true && IsCalibCSI) {
    TTreeReader* inputTreeReader = RootInput::getInstance()->GetTreeReader();
    GATCONF_ = new TTreeReaderValue<std::vector<unsigned int>>(*inputTreeReader, "GATCONF");
  }
}

void TMust2Physics::InitializeRootHistogramsCSIF(Int_t DetectorNumber) {
  string EnergyLossPath = NPOptionManager::getInstance()->GetEnergyLossPath();
  string CutsPath = NPOptionManager::getInstance()->GetCutsPath();
  DoCSIFit = true;

  if (DoCSIFit) {
    for (unsigned int i = 0; i < ParticleType.size(); i++) {
      ParticleSi[ParticleType[i].c_str()] =
          new NPL::EnergyLoss(EnergyLossPath + ParticleType[i].c_str() + "_Si.G4table", "G4Table", 100);
      ParticleAl[ParticleType[i].c_str()] =
          new NPL::EnergyLoss(EnergyLossPath + ParticleType[i].c_str() + "_Al.G4table", "G4Table", 100);
    }
  }
  unsigned int NbCSI = 16;
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TH2Map = RootHistogramsCalib::getInstance()->GetTH2Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  auto TCutGMap = RootHistogramsCalib::getInstance()->GetTCutGMap();
  auto TFileMap = RootHistogramsCalib::getInstance()->GetTFileMap();
  
  TString hnameCSIE = Form("hMM%d_SiThickness", DetectorNumber);
  TString htitleCSIE = Form("MM%d_SiThickness", DetectorNumber);
  (*TH1Map)["MUST2"][hnameCSIE] = new TH1F(hnameCSIE, htitleCSIE, 1000, 0, 400);
  
  hnameCSIE = Form("hMM%d_SiE_SiStop", DetectorNumber);
  htitleCSIE = Form("MM%d_SiE_SiStop", DetectorNumber);
  (*TH1Map)["MUST2"][hnameCSIE] = new TH1F(hnameCSIE, htitleCSIE, 500, 0, 50);
  
  hnameCSIE = Form("hMM%d_AlThickness_CsIStop", DetectorNumber);
  htitleCSIE = Form("MM%d_AlThickness_CsIStop", DetectorNumber);
  (*TH1Map)["MUST2"][hnameCSIE] = new TH1F(hnameCSIE, htitleCSIE, 500, 0, 50);

  for (Int_t j = 0; j < NbCSI; j++) {
    TString hnameCSIE = Form("hMM%d_CSI_E%d", DetectorNumber, j + 1);
    TString htitleCSIE = Form("MM%d_CSI_E%d", DetectorNumber, j + 1);
    (*TH2Map)["MUST2"][hnameCSIE] = new TH2S(hnameCSIE, htitleCSIE, 1024, 8192, 16384, 1000, 0, 60);

    if (DoCSIFit) {
      string EnergyLossPath = NPOptionManager::getInstance()->GetEnergyLossPath();
      string CutsPath = NPOptionManager::getInstance()->GetCutsPath();

      for (unsigned int i = 0; i < ParticleType.size(); i++) {
        std::cout << ParticleType[i] << "\n";
        TString CutName = Form("%s_hMM%u_CSI%u", ParticleType[i].c_str(), DetectorNumber, j + 1);
        TString cFileName = CutName + ".root";
        std::cout << CutName << "  " << cFileName << " " << CutsPath + cFileName << "\n";

        htitleCSIE = Form("%s_MM%u_CSI%u", ParticleType[i].c_str(), DetectorNumber, j + 1);
        (*TH2Map)["MUST2"][CutName] = new TH2S(CutName, htitleCSIE, 1024, 8192, 16384, 800, 0, 200);

        if((*TFileMap)["MUST2"][CutName] = new TFile(CutsPath+cFileName))
          (*TCutGMap)["MUST2"][CutName] = (TCutG*)(*TFileMap)["MUST2"][CutName]->FindObjectAny(CutName);
      }
      if(Cal_Pixel){
        std::cout << "Initializing Pixel Calib for detector " << DetectorNumber << std::endl;
        for (unsigned int i = 0; i < ParticleTypePixel.size(); i++) {
          if(m_CsI_Size%PixelSize != 0)
            std::cout << "WARNING: Your current Pixel Size: " << PixelSize << " does not divide your CsI size: " << m_CsI_Size << ". The actual size of the pixels may be dubious." << std::endl;
            unsigned int PixelNumber = (m_CsI_Size/PixelSize)*(m_CsI_Size/PixelSize);
              for (unsigned int k = 0; k < PixelNumber; k++) {
          
                // std::cout << ParticleTypePixel[i] << "\n";
                TString HistName = Form("%s_hMM%u_CSI%u_Pixel%u", ParticleTypePixel[i].c_str(), DetectorNumber, j + 1, k);

                htitleCSIE = Form("%s_MM%u_CSI%u_Pixel%u", ParticleTypePixel[i].c_str(), DetectorNumber, j + 1, k);
                (*TH2Map)["MUST2"][HistName] = new TH2S(HistName, htitleCSIE, 1024, 8192, 16384, 800, 0, 200);

          }
        }
      }
    }
  }
}

void TMust2Physics::InitializeRootHistogramsEnergyF(Int_t DetectorNumber) {
  unsigned int NbStrips = 128;
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  string Path = "../../Inputs/EnergyLoss/";
  if(Alpha_Al == nullptr)
    Alpha_Al = new NPL::EnergyLoss(Path+"alpha_Al.G4table","G4table",100);

  for (Int_t j = 0; j < NbStrips; j++) {
    TString hnameXE = Form("hMM%d_STRX_E%d", DetectorNumber, j + 1);
    TString htitleXE = Form("MM%d_STRX_E%d", DetectorNumber, j + 1);
    (*TH1Map)["MUST2"][hnameXE] = new TH1F(hnameXE, htitleXE, 16384, 0, 16384);
    // strips YE
    TString hnameYE = Form("hMM%d_STRY_E%d", DetectorNumber, j + 1);
    TString htitleYE = Form("MM%d_STRY_E%d", DetectorNumber, j + 1);
    (*TH1Map)["MUST2"][hnameYE] = new TH1F(hnameYE, htitleYE, 16384, 0, 16384);

    TString hnameFITXE = Form("hMM%d_FITX_E%d", DetectorNumber, j + 1);
    TString htitleFITXE = Form("MM%d_FITX_E%d", DetectorNumber, j + 1);
    (*TGraphMap)["MUST2"][hnameFITXE] = new TGraphErrors(3);
    (*TGraphMap)["MUST2"][hnameFITXE]->SetTitle(htitleFITXE);
    (*TGraphMap)["MUST2"][hnameFITXE]->SetName(hnameFITXE);

    TString hnameFITYE = Form("hMM%d_FITY_E%d", DetectorNumber, j + 1);
    TString htitleFITYE = Form("MM%d_FITY_E%d", DetectorNumber, j + 1);
    (*TGraphMap)["MUST2"][hnameFITYE] = new TGraphErrors(3);
    (*TGraphMap)["MUST2"][hnameFITYE]->SetTitle(htitleFITYE);
    (*TGraphMap)["MUST2"][hnameFITYE]->SetName(hnameFITYE);
  }
  TString hname = Form("SigmaFit_T%d", DetectorNumber);
  (*TH1Map)["MUST2"][hname] = new TH1F("Sigma", "Sigma from fit (channel)", 80, 0, 10);

  hname = Form("DispersionX_T%d", DetectorNumber);
  (*TGraphMap)["MUST2"][hname] = new TGraphErrors(NbStrips);
  (*TGraphMap)["MUST2"][hname]->SetTitle(hname);
  (*TGraphMap)["MUST2"][hname]->SetName(hname);
  (*TGraphMap)["MUST2"][hname]->SetMarkerStyle(2);
  (*TGraphMap)["MUST2"][hname]->SetLineColorAlpha(0,0);
  
  hname = Form("DispersionY_T%d", DetectorNumber);
  (*TGraphMap)["MUST2"][hname] = new TGraphErrors(NbStrips);
  (*TGraphMap)["MUST2"][hname]->SetTitle(hname);
  (*TGraphMap)["MUST2"][hname]->SetName(hname);
  (*TGraphMap)["MUST2"][hname]->SetMarkerStyle(2);
  (*TGraphMap)["MUST2"][hname]->SetLineColorAlpha(0,0);

  hname = Form("coeffX_a_T%d", DetectorNumber);
  (*TGraphMap)["MUST2"][hname] = new TGraphErrors(NbStrips);
  (*TGraphMap)["MUST2"][hname]->SetTitle(hname);
  (*TGraphMap)["MUST2"][hname]->SetName(hname);
  (*TGraphMap)["MUST2"][hname]->SetMarkerStyle(2);
  (*TGraphMap)["MUST2"][hname]->SetLineColorAlpha(0,0);

  hname = Form("coeffX_b_T%d", DetectorNumber);
  (*TGraphMap)["MUST2"][hname] = new TGraphErrors(NbStrips);
  (*TGraphMap)["MUST2"][hname]->SetTitle(hname);
  (*TGraphMap)["MUST2"][hname]->SetName(hname);
  (*TGraphMap)["MUST2"][hname]->SetMarkerStyle(2);
  (*TGraphMap)["MUST2"][hname]->SetLineColorAlpha(0,0);

  hname = Form("coeffY_a_T%d", DetectorNumber);
  (*TGraphMap)["MUST2"][hname] = new TGraphErrors(NbStrips);
  (*TGraphMap)["MUST2"][hname]->SetTitle(hname);
  (*TGraphMap)["MUST2"][hname]->SetName(hname);
  (*TGraphMap)["MUST2"][hname]->SetMarkerStyle(2);
  (*TGraphMap)["MUST2"][hname]->SetLineColorAlpha(0,0);

  hname = Form("coeffY_b_T%d", DetectorNumber);
  (*TGraphMap)["MUST2"][hname] = new TGraphErrors(NbStrips);
  (*TGraphMap)["MUST2"][hname]->SetTitle(hname);
  (*TGraphMap)["MUST2"][hname]->SetName(hname);
  (*TGraphMap)["MUST2"][hname]->SetMarkerStyle(2);
  (*TGraphMap)["MUST2"][hname]->SetLineColorAlpha(0,0);
}

void TMust2Physics::FillHistogramsCalib() {
  map<int, bool>::iterator it;
  if (IsCalibEnergy) {
    if (NPOptionManager::getInstance()->IsReader())
      m_EventData = &(**r_ReaderEventData);
    FillHistogramsCalibEnergyF();
  }

  if (IsCalibCSI) {
    if (!IsCalibEnergy && NPOptionManager::getInstance()->IsReader() && (**GATCONF_).size() > 0) {
      m_EventData = &(**r_ReaderEventData);
      FillHistogramsCalibCSIF();
    }
  }
  for (it = DoCalibrationTime.begin(); it != DoCalibrationTime.end(); it++)
  {
    if(it->second)
    {
      FillHistogramsCalibTimeF();
    }
  }
}

void TMust2Physics::FillHistogramsCalibEnergyF() {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  TString hnameXE, hnameYE;
  for (UShort_t i = 0; i < m_EventData->GetMMStripXEMult(); i++) {
    unsigned int DetectorNbr = m_EventData->GetMMStripXEDetectorNbr(i);
    unsigned int StripNbr = m_EventData->GetMMStripXEStripNbr(i);
    unsigned int Energy = m_EventData->GetMMStripXEEnergy(i);
    if (DoCalibrationEnergy[DetectorNbr] && EnergyXThreshold[DetectorNbr] < Energy) {
      hnameXE = Form("hMM%d_STRX_E%d", DetectorNbr, StripNbr);
      (*TH1Map)["MUST2"][hnameXE]->Fill(Energy);
    }
  }
  for (UShort_t i = 0; i < m_EventData->GetMMStripYEMult(); i++) {
    unsigned int DetectorNbr = m_EventData->GetMMStripYEDetectorNbr(i);
    unsigned int StripNbr = m_EventData->GetMMStripYEStripNbr(i);
    unsigned int Energy = m_EventData->GetMMStripYEEnergy(i);
    if (DoCalibrationEnergy[DetectorNbr] && EnergyYThreshold[DetectorNbr] > Energy) {
      hnameYE = Form("hMM%d_STRY_E%d", DetectorNbr, StripNbr);
      (*TH1Map)["MUST2"][hnameYE]->Fill(Energy);
    }
  }
}

void TMust2Physics::FillHistogramsCalibCSIF() {
  DoCalibrationCSIPreTreat();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TH2Map = RootHistogramsCalib::getInstance()->GetTH2Map();
  auto TCutGMap = RootHistogramsCalib::getInstance()->GetTCutGMap();
  auto Cal = CalibrationManager::getInstance();



  auto couple = Match_X_Y();
  auto match = Match_Si_CsI(couple);
  for(auto event: match){
    unsigned int Xindex = event.second.first;
    unsigned int Yindex = event.second.second;
    // Very important that this is an int and not unsigned int for the -1 case
    int CsIindex = event.first;
  
    unsigned int DetN = m_PreTreatedData->GetMMStripXEDetectorNbr(Xindex);
    unsigned int StripXN = m_PreTreatedData->GetMMStripXEStripNbr(Xindex);
    unsigned int StripYN = m_PreTreatedData->GetMMStripYEStripNbr(Yindex);
    
    // This part is there only to get pos of interaction (actually a bit messy, maybe could be fixed)
    Si_X.clear();
    Si_Y.clear();
    TelescopeNumber.clear();
    Si_X.push_back(StripXN);
    Si_Y.push_back(StripYN);
    TelescopeNumber.push_back(DetN);
    TVector3 BeamImpact = TVector3(0, 0, 0);
    TVector3 HitDirection = GetPositionOfInteraction(0) - BeamImpact;
    double ThetaM2Surface = HitDirection.Angle(-GetTelescopeNormal(0));
    
    if (DoCalibrationCsI.find(DetN) != DoCalibrationCsI.end() && DoCalibrationCsI[DetN] == 1) {
      
      
      double SiXE = m_PreTreatedData->GetMMStripXEEnergy(Xindex);
      double SiYE = m_PreTreatedData->GetMMStripYEEnergy(Yindex);
  
      if(CsIindex >= 0){
              
      unsigned int CSIE = m_PreTreatedData->GetMMCsIEEnergy(CsIindex);
      unsigned int CSIN = m_PreTreatedData->GetMMCsIECristalNbr(CsIindex);
      //FIXME This line should be removed, only here because I made a mistake while doing my CUTS for CSI calib
      double StripXEnergy_O = 0;
      double StripYEnergy_O = 0;
        if (m_EventData->GetMMStripXEEnergy(Xindex) > CSIEnergyXThreshold[m_EventData->GetMMStripXEDetectorNbr(Xindex)] &&
            IsValidChannel(0, m_EventData->GetMMStripXEDetectorNbr(Xindex), m_EventData->GetMMStripXEStripNbr(Xindex))) {
              static string name;
              name = "MUST2/T";
              name += NPL::itoa(m_EventData->GetMMStripXEDetectorNbr(Xindex));
              name += "_Si_X_O";
              name += NPL::itoa(m_EventData->GetMMStripXEStripNbr(Xindex));
              name += "_E";
              StripXEnergy_O = Cal->ApplyCalibration(name,m_EventData->GetMMStripXEEnergy(Xindex));
        }
        if (m_EventData->GetMMStripYEEnergy(Yindex) < CSIEnergyYThreshold[m_EventData->GetMMStripYEDetectorNbr(Yindex)] &&
            IsValidChannel(0, m_EventData->GetMMStripYEDetectorNbr(Yindex), m_EventData->GetMMStripYEStripNbr(Yindex))) {
              static string name;
              name = "MUST2/T";
              name += NPL::itoa(m_EventData->GetMMStripYEDetectorNbr(Yindex));
              name += "_Si_Y_O";
              name += NPL::itoa(m_EventData->GetMMStripYEStripNbr(Yindex));
              name += "_E";
              StripYEnergy_O = Cal->ApplyCalibration(name,m_EventData->GetMMStripYEEnergy(Yindex));
        }
              
        TString hnameCSIE = Form("hMM%d_CSI_E%d", DetN, CSIN);
        (*TH2Map)["MUST2"][hnameCSIE]->Fill(CSIE, StripXEnergy_O);
        if (DoCSIFit) {
          for (unsigned int i = 0; i < ParticleType.size(); i++) {
            TString CutName = Form("%s_hMM%u_CSI%u", ParticleType[i].c_str(), DetN, CSIN);

            if ((*TCutGMap)["MUST2"][CutName] != 0 &&
                (*TCutGMap)["MUST2"][CutName]->IsInside(CSIE, StripXEnergy_O)) {
                    double E_from_delta_E = ParticleSi[ParticleType[i].c_str()]->EvaluateEnergyFromDeltaE(
                              SiXE, SiThickness[DetN], ThetaM2Surface, 6.0 * MeV, 300.0 * MeV, 0.001 * MeV, 10000);
              // if(!(ParticleType[i] == "proton" && )
              (*TH2Map)["MUST2"][CutName]->Fill(CSIE, E_from_delta_E- SiXE);
            }
          }
          if(Cal_Pixel){
          for (unsigned int i = 0; i < ParticleTypePixel.size(); i++) {
            TString CutName = Form("%s_hMM%u_CSI%u", ParticleTypePixel[i].c_str(), DetN, CSIN);

            if ((*TCutGMap)["MUST2"][CutName] != 0 &&
                (*TCutGMap)["MUST2"][CutName]->IsInside(CSIE, StripXEnergy_O)) {
                    double E_from_delta_E = ParticleSi[ParticleTypePixel[i].c_str()]->EvaluateEnergyFromDeltaE(
                              SiXE, SiThickness[DetN], ThetaM2Surface, 6.0 * MeV, 300.0 * MeV, 0.001 * MeV, 10000);

                    auto pixel = GetPixel(event);
                    TString HistName = Form("%s_hMM%u_CSI%u_Pixel%u", ParticleTypePixel[i].c_str(), DetN, CSIN, pixel.second);
                    (*TH2Map)["MUST2"][HistName]->Fill(CSIE, E_from_delta_E- SiXE);
                    
                }
          }
        }
      }
      }
    else{
      TString hnameCSIE = Form("hMM%d_SiThickness", DetN);
      (*TH1Map)["MUST2"][hnameCSIE]->Fill(2.97*pow(SiXE,1.45)*cos(ThetaM2Surface));
      hnameCSIE = Form("hMM%d_SiE_SiStop", DetN);
      (*TH1Map)["MUST2"][hnameCSIE]->Fill(SiXE*pow(cos(ThetaM2Surface),1./1.45)); // Accounting for the angle to exptrapolate real distances
      // std::cout << ThetaM2Surface << " " <<2.97*pow(StripXEnergy,1.45)*cos(ThetaM2Surface) << std::endl; 
      }
    }
  }
}
  



  

void TMust2Physics::DoCalibrationCSIPreTreat() {
  ClearPreTreatedData();
  m_StripXEMult = m_EventData->GetMMStripXEMult();
  m_StripYEMult = m_EventData->GetMMStripYEMult();
  m_CsIEMult = m_EventData->GetMMCsIEMult();

  //   X
  //   E
  for (unsigned int i = 0; i < m_StripXEMult; ++i) {
    if (m_EventData->GetMMStripXEEnergy(i) > CSIEnergyXThreshold[m_EventData->GetMMStripXEDetectorNbr(i)] &&
        IsValidChannel(0, m_EventData->GetMMStripXEDetectorNbr(i), m_EventData->GetMMStripXEStripNbr(i))) {
      double EX = fSi_X_E(m_EventData, i);
      if(EX > 0)
        m_PreTreatedData->SetStripXE(m_EventData->GetMMStripXEDetectorNbr(i), m_EventData->GetMMStripXEStripNbr(i), EX);
    }
  }
  //   Y
  //   E
  for (unsigned int i = 0; i < m_StripYEMult; ++i) {
    if (m_EventData->GetMMStripYEEnergy(i) < CSIEnergyYThreshold[m_EventData->GetMMStripYEDetectorNbr(i)] &&
        IsValidChannel(1, m_EventData->GetMMStripYEDetectorNbr(i), m_EventData->GetMMStripYEStripNbr(i))) {
      double EY = fSi_Y_E(m_EventData, i);
      if(EY > 0)
        m_PreTreatedData->SetStripYE(m_EventData->GetMMStripYEDetectorNbr(i), m_EventData->GetMMStripYEStripNbr(i), EY);
    }
  }
  for (unsigned int i = 0; i < m_CsIEMult; ++i) {
    if (m_EventData->GetMMCsIEEnergy(i) > CSIEThreshold[m_EventData->GetMMCsIEDetectorNbr(i)] &&
        IsValidChannel(3, m_EventData->GetMMCsIEDetectorNbr(i), m_EventData->GetMMCsIECristalNbr(i))) {
      m_PreTreatedData->SetCsIE(m_EventData->GetMMCsIEDetectorNbr(i), m_EventData->GetMMCsIECristalNbr(i),
                                m_EventData->GetMMCsIEEnergy(i));
    }
  }
  return;
}

void TMust2Physics::WriteHistogramsCalib() {
  std::cout << "Writing Histograms\n";
  // map<int, bool>::iterator it;
  // for (it = DoCalibrationTime.begin(); it != DoCalibrationTime.end(); it++)
  //  {
  WriteHistogramsEnergyF();
  //  }
  //  for (it = DoCalibrationTime.begin(); it != DoCalibrationTime.end(); it++)
  //  {
  WriteHistogramsTimeF();
  //  }
  //  for (it = DoCalibrationCsI.begin(); it != DoCalibrationCsI.end(); it++)
  //  {
  WriteHistogramsCSIF();
  //  }
}

void TMust2Physics::WriteHistogramsCSIF() {
  auto File = RootHistogramsCalib::getInstance()->GetFile();
  auto TH2Map = RootHistogramsCalib::getInstance()->GetTH2Map();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  unsigned int NbCSI = 16;

  if (!File->GetDirectory("MUST2"))
    File->mkdir("MUST2");
  File->cd("MUST2");
  map<int, bool>::iterator it;
  for (it = DoCalibrationCsI.begin(); it != DoCalibrationCsI.end(); it++) {
    std::cout << "Writing Calibs for MUST2 : " << it->first << "\n";
    if (it->second) {
      if (!gDirectory->GetDirectory(Form("M2_Telescope%d", it->first))) {
        gDirectory->mkdir(Form("M2_Telescope%d", it->first));
      }
      gDirectory->cd(Form("M2_Telescope%d", it->first));
      gDirectory->mkdir("CSI");
      // gDirectory->mkdir("Time");
      gDirectory->cd("CSI");
      TString hnameCSIE = Form("hMM%d_SiThickness", it->first);
      (*TH1Map)["MUST2"][hnameCSIE]->Write();
      hnameCSIE = Form("hMM%d_SiE_SiStop", it->first);
      (*TH1Map)["MUST2"][hnameCSIE]->Write();
      for (Int_t j = 1; j <= NbCSI; j++) {
        TString hnameCSIE = Form("hMM%d_CSI_E%d", it->first, j);
        (*TH2Map)["MUST2"][hnameCSIE]->Write();
        if (DoCSIFit) {
          for (unsigned int i = 0; i < ParticleType.size(); i++) {
            TString CutName = Form("%s_hMM%u_CSI%u", ParticleType[i].c_str(), it->first, j);
            (*TH2Map)["MUST2"][CutName]->Write();
          }
          for (unsigned int i = 0; i < ParticleTypePixel.size(); i++) {
            unsigned int PixelNumber = (m_CsI_Size/PixelSize)*(m_CsI_Size/PixelSize);
            for (unsigned int k = 0; k < PixelNumber; k++) {
              TString HistName = Form("%s_hMM%u_CSI%u_Pixel%u", ParticleTypePixel[i].c_str(), it->first, j, k);
              (*TH2Map)["MUST2"][HistName]->Write();
            }
          }
        }
      }
    }
    File->cd("MUST2");
  }
  File->Close();
}

void TMust2Physics::WriteHistogramsEnergyF() {
  auto File = RootHistogramsCalib::getInstance()->GetFile();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  TString hnameXE, hnameYE;
  unsigned int NbStrips = 128;

  if (!File->GetDirectory("MUST2"))
    File->mkdir("MUST2");
  File->cd("MUST2");
  map<int, bool>::iterator it;
  for (it = DoCalibrationEnergy.begin(); it != DoCalibrationEnergy.end(); it++) {
    if (it->second) {
      if (!gDirectory->GetDirectory(Form("M2_Telescope%d", it->first)))
        gDirectory->mkdir(Form("M2_Telescope%d", it->first));
      gDirectory->cd(Form("M2_Telescope%d", it->first));
      gDirectory->mkdir("Energy");
      // gDirectory->mkdir("Time");
      gDirectory->cd("Energy");
      (*TGraphMap)["MUST2"][Form("DispersionX_T%d", it->first)]->Write();
      (*TGraphMap)["MUST2"][Form("DispersionY_T%d", it->first)]->Write();
      (*TGraphMap)["MUST2"][Form("coeffX_a_T%d", it->first)]->Write();
      (*TGraphMap)["MUST2"][Form("coeffY_a_T%d", it->first)]->Write();
      (*TGraphMap)["MUST2"][Form("coeffX_b_T%d", it->first)]->Write();
      (*TGraphMap)["MUST2"][Form("coeffY_b_T%d", it->first)]->Write();
      for (Int_t j = 0; j < NbStrips; j++) {
        TString hnameXE = Form("hMM%d_STRX_E%d", it->first, j + 1);
        TString hnameYE = Form("hMM%d_STRY_E%d", it->first, j + 1);
        (*TH1Map)["MUST2"][hnameXE]->Write();
        (*TH1Map)["MUST2"][hnameYE]->Write();
        TString hnameFITXE = Form("hMM%d_FITX_E%d", it->first, j + 1);
        TString hnameFITYE = Form("hMM%d_FITY_E%d", it->first, j + 1);
        (*TGraphMap)["MUST2"][hnameFITXE]->Write();
        (*TGraphMap)["MUST2"][hnameFITYE]->Write();
      }
    }
    File->cd("MUST2");
  }
}

void TMust2Physics::DoCalibration() {
  std::cout << "Do Calibration" << std::endl;
  map<int, bool>::iterator it;
  for (it = DoCalibrationEnergy.begin(); it != DoCalibrationEnergy.end(); it++) {
    if (it->second) {
      MakeEnergyCalibFolders();
      DoCalibrationEnergyF(it->first);
    }
  }
  for (it = DoCalibrationTime.begin(); it != DoCalibrationTime.end(); it++) {
    if (it->second) {
      DoCalibrationTimeF(it->first);
    }
  }
  for (it = DoCalibrationCsI.begin(); it != DoCalibrationCsI.end(); it++) {
    if (it->second) {
      MakeCSICalibFolders();
      DoCalibrationCsIF(it->first);
    }
  }
}

void TMust2Physics::DoCalibrationEnergyF(Int_t DetectorNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  unsigned int NbStrips = 128;
  // We assume that the pedestals are aligned at 8192 channels
  unsigned int PedestalValue = 8192;
  ofstream* calib_file = new ofstream;
  ofstream* dispersion_file = new ofstream;

  DefineCalibrationSource();

  unsigned int max_steps = 100;
  double AlThickness = 0.3*um;
  double Al_step = 0.1*um;
  double mean_extrapolation;
  std::map<unsigned int, double> a;
  std::map<unsigned int, double> b;
  std::map<unsigned int, double> dispersion;
  CreateCalibrationEnergyFiles(DetectorNumber, "X", calib_file, dispersion_file);
  
  // Prepare the calibration : find all peaks in channels for each strip and do a first calibration assuming 0.3 um of aluminum
  std::vector<double> Source_E_Slowed;
  Source_E_Slowed = SlowSource(AlThickness);
  // (*TGraphMap)["MUST2"][Form("coeffX_a_T%d", DetectorNumber)]->Clear();
  // (*TGraphMap)["MUST2"][Form("coeffX_b_T%d", DetectorNumber)]->Clear();
  // (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->Clear();
  
  for (unsigned int StripNb = 1; StripNb < NbStrips + 1; StripNb++) {
    if(FindAlphas(((*TH1Map)["MUST2"][Form("hMM%d_STRX_E%d", DetectorNumber, StripNb)]), "X", StripNb, DetectorNumber)){
      double a_temp, b_temp;
      FitLinearEnergy(((*TGraphMap)["MUST2"][Form("hMM%d_FITX_E%d", DetectorNumber, StripNb)]), "X", StripNb,
                      DetectorNumber, &a_temp, &b_temp,Source_E_Slowed);
      (*TGraphMap)["MUST2"][Form("coeffX_a_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, a_temp);
      (*TGraphMap)["MUST2"][Form("coeffX_b_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, b_temp);
      dispersion[StripNb] = PedestalValue + b_temp / a_temp;
      if(abs(dispersion[StripNb]) < 40){
        (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, dispersion[StripNb]);
        a[StripNb] = a_temp; b[StripNb] = b_temp;
      }
      else{
        a[StripNb] = 0.0; b[StripNb] = 0.0;
        dispersion[StripNb] = -1000;
      }
    }
    else{
      a[StripNb] = 0.0; b[StripNb] = 0.0;
      dispersion[StripNb] = -1000;
    }
  }
  mean_extrapolation = FindMeanExtrapolation((*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]);
  (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->SetMaximum(mean_extrapolation+30);
  (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->SetMinimum(mean_extrapolation-30);
  std::cout << "extrapo value " << mean_extrapolation << std::endl;
  double step = 0;
  bool check1 = false; bool check2 = false;
  while(abs(mean_extrapolation) > 0.1 && step < max_steps){
    std::cout << "step " << step << " " << mean_extrapolation << std::endl;
    // (*TGraphMap)["MUST2"][Form("coeffX_a_T%d", DetectorNumber)]->Clear();
    // (*TGraphMap)["MUST2"][Form("coeffX_b_T%d", DetectorNumber)]->Clear();
    // (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->Clear();
    // Here we change Al thickness, if the step is too big, we divide it by 10
    if(mean_extrapolation < 0)
    {
      AlThickness -= Al_step;
      check1=true;
    }
    else if (mean_extrapolation > 0)
    {
      AlThickness += Al_step;
      check2=true;
    }
    if(check1&&check2)
    {
      Al_step/=10.;
      check1=false;check2=false;
    }
    Source_E_Slowed = SlowSource(AlThickness);
    for (unsigned int StripNb = 1; StripNb < NbStrips + 1; StripNb++) {    
      if (AlphaMean[StripNb].size() > 0 && dispersion[StripNb] > -1000) {
        double a_temp, b_temp;
        FitLinearEnergy(((*TGraphMap)["MUST2"][Form("hMM%d_FITX_E%d", DetectorNumber, StripNb)]), "X", StripNb,
                        DetectorNumber, &a_temp, &b_temp,Source_E_Slowed);
        (*TGraphMap)["MUST2"][Form("coeffX_a_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, a_temp);
        (*TGraphMap)["MUST2"][Form("coeffX_b_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, b_temp);
        dispersion[StripNb] = PedestalValue + b_temp / a_temp;
        if(abs(dispersion[StripNb]) < 40){
          (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, dispersion[StripNb]);
          a[StripNb] = a_temp; b[StripNb] = b_temp;
        }
        else{
          a[StripNb] = 0.0; b[StripNb] = 0.0;
          dispersion[StripNb] = -1000;
        }
      }
    }
    mean_extrapolation = FindMeanExtrapolation((*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]);
    (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->SetMaximum(mean_extrapolation+30);
    (*TGraphMap)["MUST2"][Form("DispersionX_T%d", DetectorNumber)]->SetMinimum(mean_extrapolation-30);
    step++;
  }
  for (unsigned int StripNb = 1; StripNb < NbStrips + 1; StripNb++) {    
    *dispersion_file << "MUST2_T" << DetectorNumber << "_Si_X" << StripNb << "_E_Zero_Dispersion " << dispersion[StripNb] << endl;
    *calib_file << "MUST2_T" << DetectorNumber << "_Si_X" << StripNb << "_E " << b[StripNb] << " " << a[StripNb] << endl;
  }
  *dispersion_file << "MUST2_T" << DetectorNumber << " Dead Layer: " << AlThickness/um << " um of Aluminum" << endl;
  
  
  AlphaMean.clear();
  AlphaSigma.clear();
  CloseCalibrationEnergyFiles(calib_file, dispersion_file);
  
  max_steps = 100;
  AlThickness = 0.3*um;
  Al_step = 0.1*um;
  mean_extrapolation = 1000;
  a.clear();
  b.clear();
  dispersion.clear();
  
  CreateCalibrationEnergyFiles(DetectorNumber, "Y", calib_file, dispersion_file);
  
  // Prepare the calibration : find all peaks in channels for each strip and do a first calibration assuming 0.3 um of aluminum
  Source_E_Slowed.clear();
  Source_E_Slowed = SlowSource(AlThickness);
  // (*TGraphMap)["MUST2"][Form("coeffY_a_T%d", DetectorNumber)]->Clear();
  // (*TGraphMap)["MUST2"][Form("coeffY_b_T%d", DetectorNumber)]->Clear();
  // (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->Clear();
  
  for (unsigned int StripNb = 1; StripNb < NbStrips + 1; StripNb++) {
    if(FindAlphas(((*TH1Map)["MUST2"][Form("hMM%d_STRY_E%d", DetectorNumber, StripNb)]), "Y", StripNb, DetectorNumber)){
      double a_temp, b_temp;
        FitLinearEnergy(((*TGraphMap)["MUST2"][Form("hMM%d_FITY_E%d", DetectorNumber, StripNb)]), "Y", StripNb,
                        DetectorNumber, &a_temp, &b_temp,Source_E_Slowed);
        (*TGraphMap)["MUST2"][Form("coeffY_a_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, a_temp);
        (*TGraphMap)["MUST2"][Form("coeffY_b_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, b_temp);
        dispersion[StripNb] = PedestalValue + b_temp / a_temp;
        if(abs(dispersion[StripNb]) < 40){
          (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, dispersion[StripNb]);
          a[StripNb] = a_temp; b[StripNb] = b_temp;
        }
        else{
          a[StripNb] = 0.0; b[StripNb] = 0.0;
          dispersion[StripNb] = -1000;
        }
    }
    else{
      a[StripNb] = 0.0; b[StripNb] = 0.0;
      dispersion[StripNb] = -1000;
    }
  }
  mean_extrapolation = FindMeanExtrapolation((*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]);
  (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->SetMaximum(mean_extrapolation+30);
  (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->SetMinimum(mean_extrapolation-30);
  step = 0;
  check1 = false;check2 = false;
  while(abs(mean_extrapolation) > 0.1 && step < max_steps){
    // (*TGraphMap)["MUST2"][Form("coeffY_a_T%d", DetectorNumber)]->Clear();
    // (*TGraphMap)["MUST2"][Form("coeffY_b_T%d", DetectorNumber)]->Clear();
    // (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->Clear();
    // Here we change Al thickness, if the step is too big, we divide it by 10
    if(mean_extrapolation < 0)
    {
      AlThickness += Al_step;
      check1=true;
    }
    else if (mean_extrapolation > 0)
    {
      AlThickness -= Al_step;
      check2=true;
    }
    if(check1&&check2)
    {
      Al_step/=10.;
      check1=false;check2=false;
    }
    Source_E_Slowed = SlowSource(AlThickness);
    for (unsigned int StripNb = 1; StripNb < NbStrips + 1; StripNb++) {    
      if (AlphaMean[StripNb].size() > 0&& dispersion[StripNb] > -1000) {
      double a_temp, b_temp;
      FitLinearEnergy(((*TGraphMap)["MUST2"][Form("hMM%d_FITY_E%d", DetectorNumber, StripNb)]), "Y", StripNb,
                      DetectorNumber, &a_temp, &b_temp,Source_E_Slowed);
      (*TGraphMap)["MUST2"][Form("coeffY_a_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, a_temp);
      (*TGraphMap)["MUST2"][Form("coeffY_b_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, b_temp);
      dispersion[StripNb] = PedestalValue + b_temp / a_temp;
        if(abs(dispersion[StripNb]) < 40){
          (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->SetPoint(StripNb, StripNb, dispersion[StripNb]);
          a[StripNb] = a_temp; b[StripNb] = b_temp;
        }
        else{
          a[StripNb] = 0.0; b[StripNb] = 0.0;
          dispersion[StripNb] = -1000;
        }
      }
    }
    mean_extrapolation = FindMeanExtrapolation((*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]);
    (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->SetMaximum(mean_extrapolation+30);
    (*TGraphMap)["MUST2"][Form("DispersionY_T%d", DetectorNumber)]->SetMinimum(mean_extrapolation-30);
    step++;
  }
  for (unsigned int StripNb = 1; StripNb < NbStrips + 1; StripNb++) {    
    *dispersion_file << "MUST2_T" << DetectorNumber << "_Si_Y" << StripNb << "_E_Zero_Dispersion " << dispersion[StripNb] << endl;
    *calib_file << "MUST2_T" << DetectorNumber << "_Si_Y" << StripNb << "_E " << b[StripNb] << " " << a[StripNb] << endl;
  }
  *dispersion_file << "MUST2_T" << DetectorNumber << " Dead Layer: " << AlThickness/um << " um of Aluminum" << endl;
  
  
  AlphaMean.clear();
  AlphaSigma.clear();
  CloseCalibrationEnergyFiles(calib_file, dispersion_file);

}

std::vector<double> TMust2Physics::SlowSource(double AlThickness){
  std::vector<double> Slowed_Source;
  for(auto AlphaE: Source_E){
    Slowed_Source.push_back(Alpha_Al->Slow(AlphaE,AlThickness,0));
  }
  return Slowed_Source;
}

double TMust2Physics::FindMeanExtrapolation(TGraphErrors* DispersionGraph){
  return DispersionGraph->GetMean(2);
}

void TMust2Physics::DefineCalibrationSource() {
  // 239Pu
  Source_isotope.push_back("$^{239}$Pu");
  Source_E.push_back(5.15659);
  Source_Sig.push_back(0.00014);
  Source_branching_ratio.push_back(70.77);
  Source_isotope.push_back("$^{239}$Pu");
  Source_E.push_back(5.14438);
  Source_Sig.push_back(0.00014);
  Source_branching_ratio.push_back(17.11);
  Source_isotope.push_back("$^{239}$Pu");
  Source_E.push_back(5.1055);
  Source_Sig.push_back(0.00014);
  Source_branching_ratio.push_back(11.94);
  // 241Am
  Source_isotope.push_back("$^{241}$Am");
  Source_E.push_back(5.48556);
  Source_Sig.push_back(0.00012);
  Source_branching_ratio.push_back(84.8);
  Source_isotope.push_back("$^{241}$Am");
  Source_E.push_back(5.44280);
  Source_Sig.push_back(0.00012);
  Source_branching_ratio.push_back(13.1);
  Source_isotope.push_back("$^{241}$Am");
  Source_E.push_back(5.388);
  Source_Sig.push_back(0.00012);
  Source_branching_ratio.push_back(1.66);
  // 244Cm
  Source_isotope.push_back("$^{244}$Cm");
  Source_E.push_back(5.80477);
  Source_Sig.push_back(0.00005);
  Source_branching_ratio.push_back(76.40);
  Source_isotope.push_back("$^{244}$Cm");
  Source_E.push_back(5.76264);
  Source_Sig.push_back(0.00005);
  Source_branching_ratio.push_back(23.60);
}

void TMust2Physics::FitLinearEnergy(TGraphErrors* FitHist, TString side, unsigned int StripNb,
                                    unsigned int DetectorNumber, double* a, double* b,std::vector<double> Source_E_slowed){                               
  // FitHist->Reset();                
  if (AlphaMean[StripNb].size() == 3 && AlphaSigma[StripNb].size() == 3) {
    double AlphaSourceEnergy[3];
    double AlphaSourceSigma[3];
    // double AlphaMeanP[3];
    // double AlphaSigmaP[3];
    if (side == "X") {
      for (unsigned int i = 0; i < 3; i++) {
        AlphaSourceEnergy[i] = Source_E_slowed[3 * i];
        AlphaSourceSigma[i] = Source_Sig[3 * i];
        // AlphaMeanP[i] = AlphaMean[i];
        // AlphaSigmaP[i] = AlphaSigma[i];
      }
    }
    else if (side == "Y") {
      for (unsigned int i = 0; i < 3; i++) {
        AlphaSourceEnergy[i] = Source_E_slowed[6 - 3 * i];
        AlphaSourceSigma[i] = Source_Sig[6 - 3 * i];
        // AlphaMeanP[i] = AlphaMean[i];
        // AlphaSigmaP[i] = AlphaSigma[i];
      }
    }
    else {
      std::cout << "Check side string formatting" << std::endl;
    }
    for (Int_t p = 0; p < 3; p++) {
      FitHist->SetPoint(p, AlphaMean[StripNb][p], AlphaSourceEnergy[p]);
      FitHist->SetPointError(p, AlphaSigma[StripNb][p], AlphaSourceSigma[p]);
    }

    TF1* f1 = new TF1("f1", "[1]+[0]*x");
    if (side == "X") {
      f1->SetParameter(0, 0.007);
      f1->SetParameter(1, -60);
    }
    else if (side == "Y") {
      f1->SetParameter(0, -0.007);
      f1->SetParameter(1, 60);
    }
    FitHist->Fit("f1", "Q");

    *a = f1->GetParameter(0);
    *b = f1->GetParameter(1);
  }
}

bool TMust2Physics::FindAlphas(TH1* CalibHist, TString side, unsigned int StripNb, unsigned int DetectorNumber) {
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  Double_t ResSigma = 5;
  Double_t ResSigmaTSpec = 1;
  Double_t Threshold = 0.05;
  Int_t Npeaks = 3; // maximum number of peaks that can be found

  //////// Peak finder
  TSpectrum* s = new TSpectrum(Npeaks, ResSigmaTSpec);

  Int_t Nfound = s->Search(CalibHist, ResSigma, "new", Threshold);
  Double_t* Xpeaks = s->GetPositionX();

  // If 3 peaks found
  if (Nfound == 3) {
    for (Int_t p = 0; p < Nfound; p++) {
      for (Int_t i = 0; i < Nfound - 1; i++) {
        if (Xpeaks[i] > Xpeaks[i + 1]) {
          Float_t varia = Xpeaks[i];
          Xpeaks[i] = Xpeaks[i + 1];
          Xpeaks[i + 1] = varia;
        }
      }
    }
    if ((AlphaFitType[DetectorNumber] == "NoSatellite") || (AlphaFitType[DetectorNumber] == "")) {
      Float_t linf = 0, lsup = 0;
      for (Int_t p = 0; p < Nfound; p++) {
        if (side == "X") {
          // linf = Xpeaks[p]-2;
          // lsup = Xpeaks[p]+8;
          linf = Xpeaks[p] - 2;
          lsup = Xpeaks[p] + 8;
        }

        else if (side == "Y") {
          linf = Xpeaks[p] - 8;
          lsup = Xpeaks[p] + 2;
          // linf = Xpeaks[p]-8;
          // lsup = Xpeaks[p]+2;
        }
        TF1* gauss = new TF1("gauss", "gaus", linf, lsup);
        CalibHist->Fit(gauss, "RQ");
        AlphaMean[StripNb].push_back(gauss->GetParameter(1));
        AlphaSigma[StripNb].push_back(gauss->GetParameter(2));
        (*TH1Map)["MUST2"][Form("SigmaFit_T%d", DetectorNumber)]->Fill(gauss->GetParameter(2));
      }
    }
    else if (AlphaFitType[DetectorNumber] == "WithSatellite") {
      Float_t linf[3];
      Float_t lsup[3];
      for (Int_t p = 0; p < Nfound; p++) {
        if (side == "X") {
          linf[p] = Xpeaks[p] - 20;
          lsup[p] = Xpeaks[p] + 20;
        }

        else if (side == "Y") {
          linf[p] = Xpeaks[p] + 20;
          lsup[p] = Xpeaks[p] - 20;
        }
      }
      TF1* SatellitePu = new TF1("gauss", source_Pu, linf[0], lsup[0], 4);
      SatellitePu->SetParameters(150, Xpeaks[0], Xpeaks[0] - 1, 0.1);
      CalibHist->Fit(SatellitePu, "RQ+");
      AlphaMean[StripNb].push_back(SatellitePu->GetParameter(1));
      AlphaSigma[StripNb].push_back(SatellitePu->GetParameter(3));
      (*TH1Map)["MUST2"][Form("SigmaFit_T%d", DetectorNumber)]->Fill(SatellitePu->GetParameter(3));

      TF1* SatelliteAm = new TF1("gauss", source_Am, linf[1], lsup[1], 4);
      SatelliteAm->SetParameters(150, Xpeaks[1], Xpeaks[1] - 1, 0.1);
      CalibHist->Fit(SatelliteAm, "RQ+");
      AlphaMean[StripNb].push_back(SatelliteAm->GetParameter(1));
      AlphaSigma[StripNb].push_back(SatelliteAm->GetParameter(3));
      (*TH1Map)["MUST2"][Form("SigmaFit_T%d", DetectorNumber)]->Fill(SatelliteAm->GetParameter(3));

      TF1* SatelliteCm = new TF1("gauss", source_Cm, linf[2], lsup[2], 4);
      SatelliteCm->SetParameters(150, Xpeaks[2], Xpeaks[2] - 1, 0.1);
      CalibHist->Fit(SatelliteCm, "RQ+");
      AlphaMean[StripNb].push_back(SatelliteCm->GetParameter(1));
      AlphaSigma[StripNb].push_back(SatelliteCm->GetParameter(3));
      (*TH1Map)["MUST2"][Form("SigmaFit_T%d", DetectorNumber)]->Fill(SatelliteCm->GetParameter(3));
    }
  }

  if (Nfound != 3) {
    BadStrip[side][StripNb] = Nfound;
    return false;
  }

  return true;
}

Double_t TMust2Physics::source_Pu(Double_t* x, Double_t* par) {
  // [0] : constant
  // [1] : position peak1
  // [2] : position peak2
  // [3] : sigma

  Double_t arg1 = 0;
  Double_t arg2 = 0;
  Double_t arg3 = 0;

  if (par[4] != 0) {
    arg1 = (x[0] - par[1]) / par[3];
    arg2 = (x[0] - par[2]) / par[3];
    arg3 = (x[0] - par[1] - (par[1] - par[2]) * (5.15659 - 5.1055) / (5.15659 - 5.14438)) / par[3];
  }

  else
    cout << " Attention, sigma est nul !" << endl;

  Double_t gaus1 = par[0] * exp(-0.5 * arg1 * arg1);
  Double_t gaus2 = 15.1 / 73.8 * par[0] * exp(-0.5 * arg2 * arg2);
  Double_t gaus3 = 11.5 / 73.8 * par[0] * exp(-0.5 * arg3 * arg3);
  Double_t fitval = gaus1 + gaus2 + gaus3;

  return fitval;
}

Double_t TMust2Physics::source_Am(Double_t* x, Double_t* par) {
  // [0] : constant
  // [1] : position peak1
  // [2] : position peak2
  // [3] : sigma

  Double_t arg1 = 0;
  Double_t arg2 = 0;
  Double_t arg3 = 0;

  if (par[4] != 0) {
    arg1 = (x[0] - par[1]) / par[3];
    arg2 = (x[0] - par[2]) / par[3];
    arg3 = (x[0] - par[1] - (par[1] - par[2]) * (5.48556 - 5.388) / (5.48556 - 5.44280)) / par[3];
  }

  else
    cout << " Attention, sigma est nul !" << endl;

  Double_t gaus1 = par[0] * exp(-0.5 * arg1 * arg1);
  Double_t gaus2 = 13.0 / 84.5 * par[0] * exp(-0.5 * arg2 * arg2);
  Double_t gaus3 = 1.6 / 84.5 * par[0] * exp(-0.5 * arg3 * arg3);
  Double_t fitval = gaus1 + gaus2 + gaus3;

  return fitval;
}

Double_t TMust2Physics::source_Cm(Double_t* x, Double_t* par) {
  // [0] : constante
  // [1] : position peak1
  // [2] : position peak2
  // [3] : sigma

  Double_t arg1 = 0;
  Double_t arg2 = 0;

  if (par[3] != 0) {
    arg1 = (x[0] - par[1]) / par[3];
    arg2 = (x[0] - par[2]) / par[3];
  }

  else
    cout << " Attention, sigma est nul !" << endl;

  Double_t gaus1 = par[0] * exp(-0.5 * arg1 * arg1);
  Double_t gaus2 = 23.6 / 76.4 * par[0] * exp(-0.5 * arg2 * arg2);
  Double_t fitval = gaus1 + gaus2;

  return fitval;
}

void TMust2Physics::MakeEnergyCalibFolders() {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString test_folder = "test -f " + Path + OutputName;
  TString make_folder = "mkdir " + Path + OutputName;

  int sys =system(make_folder);
  sys =system(make_folder + "/peaks");
  sys =system(make_folder + "/dispersion");
  sys =system(make_folder + "/latex");
  sys =system(make_folder + "/latex/pictures");
}

void TMust2Physics::MakeCSICalibFolders() {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString test_folder = "test -f " + Path + OutputName;
  TString make_folder = "mkdir " + Path + OutputName;

  int sys =system(make_folder);
}

void TMust2Physics::CreateCalibrationEnergyFiles(unsigned int DetectorNumber, TString side, ofstream* calib_file,
                                                 ofstream* dispersion_file) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "Cal_Str_" + side + "_E_MM" + std::to_string(DetectorNumber);
  //  fname =  Path+OutputName + "/peaks/" + str1 + ".peak";
  // peaks_file.open( ( (string)(Path+OutputName+"/peaks/"+Filename+".peak") ).c_str() );
  (*calib_file).open(((string)(Path + OutputName + "/" + Filename + ".cal")).c_str());
  (*dispersion_file).open(((string)(Path + OutputName + "/dispersion/" + Filename + ".dispersion")).c_str());
}

void TMust2Physics::CreateCalibrationCSIFiles(unsigned int DetectorNumber, ofstream* calib_file, TString ParticleType) {
  std::string Path = NPOptionManager::getInstance()->GetCalibrationOutputPath();
  std::string OutputName = NPOptionManager::getInstance()->GetOutputFile();
  if (OutputName.size() > 5) {
    if (OutputName.substr(OutputName.size() - 5, OutputName.size()) == ".root") {
      OutputName = OutputName.substr(0, OutputName.size() - 5);
    }
  }
  TString Filename = "CsI_MM" + std::to_string(DetectorNumber) + "_" + ParticleType;
  std::cout << Filename << std::endl;
  (*calib_file).open(((string)(Path + OutputName + "/" + Filename + ".txt")).c_str());
  std::cout << (string)(Path + OutputName + "/" + Filename + ".txt") << std::endl;
}

void TMust2Physics::CloseCalibrationCSIFiles(ofstream* calib_file) { calib_file->close(); }

void TMust2Physics::CloseCalibrationEnergyFiles(ofstream* calib_file, ofstream* dispersion_file) {
  // peaks_file.close();
  calib_file->close();
  dispersion_file->close();
}

void TMust2Physics::DoCalibrationTimeF(Int_t DetectorNumber) {}

void TMust2Physics::DoCalibrationCsIF(Int_t DetectorNumber) {
  gErrorIgnoreLevel = kWarning;
  TF1* Gaus = new TF1("Gaus", "gaus", 0, 200);
  // TF1* f3 = new TF1("f3","[0] + [1]*(x-8192) + [2]*(x-[3])*(x-[3])/(1+exp([3] - x)/[4])", 8192, 16384);
  TF1* f4 = new TF1("f4", "pol4", 8192, 16384);
  
  auto File = new TFile("./FitSlices.root", "RECREATE");
  auto TH2Map = RootHistogramsCalib::getInstance()->GetTH2Map();
  auto TH1Map = RootHistogramsCalib::getInstance()->GetTH1Map();
  auto TGraphMap = RootHistogramsCalib::getInstance()->GetTGraphMap();
  unsigned int NbCSI = 16;

  ofstream* calib_file = new ofstream;

  for (unsigned int j = 0; j < ParticleType.size(); j++) {
    CreateCalibrationCSIFiles(DetectorNumber, calib_file, ParticleType[j]);
    for (unsigned int i = 1; i <= NbCSI; i++) {
      TString CutName = Form("%s_hMM%u_CSI%u", ParticleType[j].c_str(), DetectorNumber, i);
      TString htitleCSIE = Form("%s_hMM%u_CSI%u", ParticleType[j].c_str(), DetectorNumber, i);
      double a = 0;
      double b = 0;
      double c = 0;
      double d = 0;
      double e = 0;
      // double f = 0;
      if ((*TH2Map)["MUST2"][CutName] != 0&& (*TH2Map)["MUST2"][CutName]->Integral() > 1000) {
        // f4->SetParameters()
        std::cout << "Fit on CSI " << i << " Detector " << DetectorNumber << " with particle " << ParticleType[j] << std::endl;
        int Res = -1;
        Res = (*TH2Map)["MUST2"][CutName]->Fit(f4,"BFM");
        File->Write();
        std::cout << Res << std::endl;
        if(Res != -1){

        a = f4->GetParameter(0);
        b = f4->GetParameter(1);
        c = f4->GetParameter(2);
        d = f4->GetParameter(3);
        e = f4->GetParameter(4);
        // f = f4->GetParameter(5);
        }
      }
      *calib_file << ParticleType[j].c_str() << "_MUST2_T" << DetectorNumber << "_CsI" << i << "_E " << a << " " << b << " " << c << " " << d << " " << e << " " << endl;
      if(Cal_Pixel && std::find(ParticleTypePixel.begin(),ParticleTypePixel.end(),ParticleType[j])!=ParticleTypePixel.end()){
        unsigned int PixelNumber = (m_CsI_Size/PixelSize)*(m_CsI_Size/PixelSize);
        for (unsigned int k = 0; k < PixelNumber; k++) {
          CutName = Form("%s_hMM%u_CSI%u_Pixel%u", ParticleType[j].c_str(), DetectorNumber, i,k);
          htitleCSIE = Form("%s_hMM%u_CSI%u_Pixel%u", ParticleType[j].c_str(), DetectorNumber, i,k);
          a = 0;
          b = 0;
          c = 0;
          d = 0;
          e = 0;
          // f = 0;
          if ((*TH2Map)["MUST2"][CutName] != 0 && (*TH2Map)["MUST2"][CutName]->Integral() > 1000) {
            std::cout << "Fit on CSI " << i << " Pixel " << k << " Detector " << DetectorNumber << " with particle " << ParticleType[j] << std::endl;
            int Res = -1;
            Res = (*TH2Map)["MUST2"][CutName]->Fit(f4,"BFM");
            File->Write();
            std::cout << Res << std::endl;
            if(Res != -1){

            a = f4->GetParameter(0);
            b = f4->GetParameter(1);
            c = f4->GetParameter(2);
            d = f4->GetParameter(3);
            e = f4->GetParameter(4);
            // f = f4->GetParameter(5);
            }
          }
          *calib_file << ParticleType[j].c_str() << "_MUST2_T" << DetectorNumber << "_CsI" << i << "_Pixel" << k <<  "_E " << a << " " << b << " " << c << " " << d << " " << e << " " << endl;
          
        }
      }
    }
    if(Cal_Pixel && std::find(ParticleTypePixel.begin(),ParticleTypePixel.end(),ParticleType[j])!=ParticleTypePixel.end()){
      *calib_file << "CsI size used for pixel cal: " << m_CsI_Size<< endl;
      *calib_file << "Pixel size used for pixel cal: " << PixelSize<< endl;
    }
    CloseCalibrationCSIFiles(calib_file);
  }
}

///////////////////////////////////////////////////////////////////////////
namespace MUST2_LOCAL {
  //   DSSD
  //   X
  double fSi_X_E(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMStripXEDetectorNbr(i));
    name += "_Si_X";
    name += NPL::itoa(m_EventData->GetMMStripXEStripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMStripXEEnergy(i));
  }

  double fSi_X_T(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMStripXTDetectorNbr(i));
    name += "_Si_X";
    name += NPL::itoa(m_EventData->GetMMStripXTStripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMStripXTTime(i));
  }

  //   Y
  double fSi_Y_E(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMStripYEDetectorNbr(i));
    name += "_Si_Y";
    name += NPL::itoa(m_EventData->GetMMStripYEStripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMStripYEEnergy(i));
  }

  double fSi_Y_T(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMStripYTDetectorNbr(i));
    name += "_Si_Y";
    name += NPL::itoa(m_EventData->GetMMStripYTStripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMStripYTTime(i));
  }

  //   SiLi
  double fSiLi_E(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMSiLiEDetectorNbr(i));
    name += "_SiLi";
    name += NPL::itoa(m_EventData->GetMMSiLiEPadNbr(i));
    name += "_E";

    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMSiLiEEnergy(i));
  }

  double fSiLi_T(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMSiLiTDetectorNbr(i));
    name += "_SiLi";
    name += NPL::itoa(m_EventData->GetMMSiLiTPadNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMSiLiTTime(i));
  }

  //   CsI
  double fCsI_E(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMCsIEDetectorNbr(i));
    name += "_CsI";
    name += NPL::itoa(m_EventData->GetMMCsIECristalNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMCsIEEnergy(i));
  }

  double fCsI_T(const TMust2Data* m_EventData, const int& i) {
    static string name;
    name = "MUST2/T";
    name += NPL::itoa(m_EventData->GetMMCsITDetectorNbr(i));
    name += "_CsI";
    name += NPL::itoa(m_EventData->GetMMCsITCristalNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetMMCsITTime(i));
  }
} // namespace MUST2_LOCAL

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TMust2Physics::Construct() { return (NPL::VDetector*)new TMust2Physics(); }

NPL::VTreeReader* TMust2Physics::ConstructReader() { return (NPL::VTreeReader*)new TMust2PhysicsReader(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_must2 {
 public:
  proxy_must2() {
    NPL::DetectorFactory::getInstance()->AddToken("M2Telescope", "MUST2");
    NPL::DetectorFactory::getInstance()->AddDetector("M2Telescope", TMust2Physics::Construct);
    NPL::DetectorFactory::getInstance()->AddDetectorReader("M2Telescope", TMust2Physics::ConstructReader);
  }
};

proxy_must2 p_must2;
}
