/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold QQQ3 Raw data                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

#include "TQQQ3Data.h"

ClassImp(TQQQ3Data)

    /////////////////////////
    TQQQ3Data::TQQQ3Data() {}

/////////////////////////
TQQQ3Data::~TQQQ3Data() {}

/////////////////////////
void TQQQ3Data::Clear() {
  fQQQ3_StripFront_DetectorNbr.clear();
  fQQQ3_StripFront_StripNbr.clear();
  fQQQ3_StripFront_Energy.clear();
  fQQQ3_StripFront_TimeCFD.clear();
  fQQQ3_StripFront_TimeLED.clear();
  fQQQ3_StripFront_Time.clear();

  fQQQ3_StripBack_DetectorNbr.clear();
  fQQQ3_StripBack_StripNbr.clear();
  fQQQ3_StripBack_Energy.clear();
  fQQQ3_StripBack_TimeCFD.clear();
  fQQQ3_StripBack_TimeLED.clear();
  fQQQ3_StripBack_Time.clear();

  fQQQ3_PAD_DetectorNbr.clear();
  fQQQ3_PAD_Energy.clear();
  fQQQ3_PAD_TimeCFD.clear();
  fQQQ3_PAD_TimeLED.clear();
  fQQQ3_PAD_Time.clear();
}

/////////////////////////
void TQQQ3Data::Dump() const {

  // Front
  cout << "QQQ3 Strip Front Mult = " << fQQQ3_StripFront_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fQQQ3_StripFront_DetectorNbr.size(); i++) {
    cout << "DetNbr (Front): " << fQQQ3_StripFront_DetectorNbr[i] << "   Strip: " << fQQQ3_StripFront_StripNbr[i]
         << "   Energy: " << fQQQ3_StripFront_Energy[i] << "   Time CFD: " << fQQQ3_StripFront_TimeCFD[i]
         << "   Time LED: " << fQQQ3_StripFront_TimeLED[i] << endl;
  }

  // Back
  cout << "QQQ3 Strip Back Mult  = " << fQQQ3_StripBack_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fQQQ3_StripBack_DetectorNbr.size(); i++) {
    cout << "DetNbr (Back): " << fQQQ3_StripBack_DetectorNbr[i] << "   Strip: " << fQQQ3_StripBack_StripNbr[i]
         << "   Energy: " << fQQQ3_StripBack_Energy[i] << "   Time CFD: " << fQQQ3_StripBack_TimeCFD[i]
         << "   Time LED: " << fQQQ3_StripBack_TimeLED[i] << endl;
  }

  // PAD
  cout << "QQQ3 Strip PAD Mult  = " << fQQQ3_PAD_DetectorNbr.size() << endl;
  for (UShort_t i = 0; i < fQQQ3_PAD_DetectorNbr.size(); i++) {
    cout << "DetNbr (PAD): " << fQQQ3_PAD_DetectorNbr[i] << "   Energy: " << fQQQ3_PAD_Energy[i]
         << "   Time CFD: " << fQQQ3_PAD_TimeCFD[i] << "   Time LED: " << fQQQ3_PAD_TimeLED[i] << endl;
  }
}
