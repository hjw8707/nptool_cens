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

#include "TExogamData.h"

ClassImp(TExogamData)

    TExogamData::TExogamData() {
  // Default constructor
  Clear();
}

TExogamData::~TExogamData() {}

void TExogamData::Clear() {
  fExo_Crystal.clear();
  fExo_E.clear();
  fExo_E_HG.clear(); // High gain x20
  fExo_TS.clear();
  fExo_TDC.clear();
  fExo_BGO.clear();
  fExo_CsI.clear();
  fExo_Outer1.clear();
  fExo_Outer2.clear();
  fExo_Outer3.clear();
  fExo_Outer4.clear();
}

void TExogamData::Dump() const {}

