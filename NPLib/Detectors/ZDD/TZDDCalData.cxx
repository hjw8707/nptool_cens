/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Cal data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TZDDCalData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TZDDCalData)


//////////////////////////////////////////////////////////////////////
TZDDCalData::TZDDCalData() {
}



//////////////////////////////////////////////////////////////////////
TZDDCalData::~TZDDCalData() {
}



//////////////////////////////////////////////////////////////////////
void TZDDCalData::Clear() {
    cZDD_IC_E.clear();
    cZDD_IC_N.clear();
    cZDD_IC_TS.clear();

    cZDD_PL_E.clear();
    cZDD_PL_N.clear();
    cZDD_PL_TS.clear();

    cZDD_DC_E.clear();
    cZDD_DC_N.clear();
    cZDD_DC_TS.clear();

    cZDD_EXO_E.clear();
    cZDD_EXO_N.clear();
    cZDD_EXO_TS.clear();
}



//////////////////////////////////////////////////////////////////////
void TZDDCalData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TZDDCalData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  cout << "ZDD_ICE Mult: " << cZDD_IC_E.size() << endl;
  for (size_t i = 0 ; i < cZDD_IC_E.size(); i++){
    cout << "DetNbr: " << cZDD_IC_N[i] << " Energy: " << cZDD_IC_E[i] << " TimeStamp: " << cZDD_IC_TS[i] << endl;
  }
  cout << "ZDD_PLE Mult: " << cZDD_PL_E.size() << endl;
  for (size_t i = 0 ; i < cZDD_PL_E.size(); i++){
    cout << "DetNbr: " << cZDD_PL_N[i] << " Energy: " << cZDD_PL_E[i] << " TimeStamp: " << cZDD_PL_TS[i] << endl;
  }
  cout << "ZDD_DCE Mult: " << cZDD_DC_E.size() << endl;
  for (size_t i = 0 ; i < cZDD_DC_E.size(); i++){
    cout << "DetNbr: " << cZDD_DC_N[i] << " Energy: " << cZDD_DC_E[i] << " TimeStamp: " << cZDD_DC_TS[i] << endl;
  }
  cout << "ZDD_EXOE Mult: " << cZDD_EXO_E.size() << endl;
  for (size_t i = 0 ; i < cZDD_EXO_E.size(); i++){
    cout << "DetNbr: " << cZDD_EXO_N[i] << " Energy: " << cZDD_EXO_E[i] << " TimeStamp: " << cZDD_EXO_TS[i] << endl;
  }  
}
