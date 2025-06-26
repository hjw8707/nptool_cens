/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2009                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Exogam Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
using namespace std;

#include "TExogamCalData.h"


ClassImp(TExogamCalData)

TExogamCalData::TExogamCalData() {
  // Default constructor
  Clear();
}

TExogamCalData::~TExogamCalData() {}

void TExogamCalData::Clear() {
  cExo_Crystal.clear();
  cExo_E.clear();
  cExo_E_HG.clear(); // High gain x20
  cExo_TS.clear();
  cExo_TDC.clear();
  cExo_BGO.clear();
  cExo_CsI.clear();
  cExo_Outer1.clear();
  cExo_Outer2.clear();
  cExo_Outer3.clear();
  cExo_Outer4.clear();
}

void TExogamCalData::Dump() const {}
