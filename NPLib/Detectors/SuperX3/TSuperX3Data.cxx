/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : january 2011                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *     This class holds the raw data storage for the SuperX3 detector from Micron *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TSuperX3Data.h"

#include <iostream>
using namespace std;

ClassImp(TSuperX3Data)

    TSuperX3Data::TSuperX3Data() {
  // Default constructor
  Clear();
}

TSuperX3Data::~TSuperX3Data() {}

void TSuperX3Data::Clear() {
  // (Up, E)
  fSuperX3_UpE_DetectorNbr.clear();
  fSuperX3_UpE_StripNbr.clear();
  fSuperX3_UpE_Energy.clear();
  // (Up, T)
  fSuperX3_UpT_DetectorNbr.clear();
  fSuperX3_UpT_StripNbr.clear();
  fSuperX3_UpT_Time.clear();
  // (Down, E)
  fSuperX3_DownE_DetectorNbr.clear();
  fSuperX3_DownE_StripNbr.clear();
  fSuperX3_DownE_Energy.clear();
  // (Down, T)
  fSuperX3_DownT_DetectorNbr.clear();
  fSuperX3_DownT_StripNbr.clear();
  fSuperX3_DownT_Time.clear();

  // (Back, E)
  fSuperX3_BackE_DetectorNbr.clear();
  fSuperX3_BackE_StripNbr.clear();
  fSuperX3_BackE_Energy.clear();
  // (Back, T)
  fSuperX3_BackT_DetectorNbr.clear();
  fSuperX3_BackT_StripNbr.clear();
  fSuperX3_BackT_Time.clear();
}

void TSuperX3Data::Dump() const {}
