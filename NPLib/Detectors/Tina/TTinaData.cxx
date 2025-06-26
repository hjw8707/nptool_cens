/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Tina Raw data                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TTinaData.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std; 

ClassImp(TTinaData)

//////////////////////////////////////////////////////////////////////
TTinaData::TTinaData() {
}

//////////////////////////////////////////////////////////////////////
TTinaData::~TTinaData() {
}

//////////////////////////////////////////////////////////////////////
void TTinaData::Clear() {
  fKE = 0;
  // TTT front
  fTTTfront_E_DetectorNbr.clear();
  fTTTfront_E_StripNbr.clear();
  fTTTfront_Energy.clear();
  fTTTfront_T_DetectorNbr.clear();
  fTTTfront_T_StripNbr.clear();
  fTTTfront_Time.clear();
  // TTT back
  fTTTback_E_DetectorNbr.clear();
  fTTTback_E_StripNbr.clear();
  fTTTback_Energy.clear();
  fTTTback_T_DetectorNbr.clear();
  fTTTback_T_StripNbr.clear();
  fTTTback_Time.clear();
  // Pad
  fPad_E_DetectorNbr.clear();
  fPad_Energy.clear();
  fPad_T_DetectorNbr.clear();
  fPad_Time.clear();
  // YY1 ring
  fYY1ring_E_DetectorNbr.clear();
  fYY1ring_E_StripNbr.clear();
  fYY1ring_Energy.clear();
  fYY1ring_T_DetectorNbr.clear();
  fYY1ring_T_StripNbr.clear();
  fYY1ring_Time.clear();
  // YY1 sector
  fYY1sector_E_DetectorNbr.clear();
  fYY1sector_E_StripNbr.clear();
  fYY1sector_Energy.clear();
  fYY1sector_T_DetectorNbr.clear();
  fYY1sector_T_StripNbr.clear();
  fYY1sector_Time.clear();
  // CsI
  fCsI_E_DetectorNbr.clear();
  fCsI_Energy.clear();
  fCsI_T_DetectorNbr.clear();
  fCsI_Time.clear();
}

//////////////////////////////////////////////////////////////////////
void TTinaData::Dump() const {
  // this method is very useful for debugging and worth the development
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TTinaData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // TTT front
  size_t mysize = fTTTfront_E_DetectorNbr.size();
  cout << "TinaTTTfrontMultEnergy: " << mysize << endl;
  for (size_t i = 0 ; i < mysize ; i++){
    cout << " DetNbr: " << fTTTfront_E_DetectorNbr[i]
         << " Energy: " << fTTTfront_Energy[i];
  }
  mysize = fTTTfront_T_DetectorNbr.size();
  cout << "TinaTTTfront_T_MultTime: " << mysize << endl;
  for (size_t i = 0 ; i < mysize ; i++){
    cout << " DetNbr: " << fTTTfront_T_DetectorNbr[i]
         << " Time  : " << fTTTfront_Time[i];
  }

  // TTT back
  mysize = fTTTback_E_DetectorNbr.size();
  cout << "TinaTTTbackMultEnergy: " << mysize << endl;
  for (size_t i = 0 ; i < mysize ; i++){
    cout << " DetNbr: " << fTTTback_E_DetectorNbr[i]
         << " Energy: " << fTTTback_Energy[i];
  }
  mysize = fTTTback_T_DetectorNbr.size();
  cout << "TinaTTTback_T_MultTime: " << mysize << endl;
  for (size_t i = 0 ; i < mysize ; i++){
    cout << " DetNbr: " << fTTTback_T_DetectorNbr[i]
         << " Time  : " << fTTTback_Time[i];
  }

  // sf warning: need to add similar routine for Pad, YY1, and CsI

}
