/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jw.hwang@ibs.re.kr       *
 *                                                                           *
 * Creation Date  : February 2021                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the GRAPE  raw data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

#include "TGRAPEData.h"

ClassImp(TGRAPEData)

/////////////////////////
TGRAPEData::TGRAPEData(){
}

/////////////////////////
TGRAPEData::~TGRAPEData(){
}

/////////////////////////
void TGRAPEData::Clear(){
  fGRAPE_Ge_GRAPENbr.clear();
  fGRAPE_Ge_CrystalNbr.clear();
  fGRAPE_Ge_SegmentNbr.clear();
  fGRAPE_Ge_Energy.clear();
  fGRAPE_Ge_TimeCFD.clear();
  fGRAPE_Ge_TimeLED.clear();
}

/////////////////////////
void TGRAPEData::Dump() const{
  // Energy
 // cout << "GRAPE_Mult = " << fGRAPE_GRAPENbr.size() << endl;
  
  // Front
 // for (UShort_t i = 0; i < fTIG_GRAPENbr.size(); i++){
 //   cout << "GRAPE: " << fTIG_GRAPENbr[i]
 //        << " Crystal: " << fTIG_CrystalNbr[i]
 //        << " Energy: " << fTIG_Energy[i]
 //        << " Time: " << fTIG_Time[i] << endl;
 // }
}
