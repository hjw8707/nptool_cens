/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TZDDData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TZDDData)


//////////////////////////////////////////////////////////////////////
TZDDData::TZDDData() {
}



//////////////////////////////////////////////////////////////////////
TZDDData::~TZDDData() {
}



//////////////////////////////////////////////////////////////////////
void TZDDData::Clear() {
    fZDD_IC_E.clear();
    fZDD_IC_N.clear();
    fZDD_IC_TS.clear();

    fZDD_PL_E.clear();
    fZDD_PL_N.clear();
    fZDD_PL_TS.clear();

    fZDD_DC_E.clear();
    fZDD_DC_N.clear();
    fZDD_DC_TS.clear();

    fZDD_Drift_DetectorNbr.clear();
    fZDD_DriftTime.clear();
    fZDD_X.clear();

    fZDD_EXO_E.clear();
    fZDD_EXO_N.clear();
    fZDD_EXO_TS.clear();
}



//////////////////////////////////////////////////////////////////////
void TZDDData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TZDDData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  cout << "ZDD_ICE Mult: " << fZDD_IC_E.size() << endl;
  for (size_t i = 0 ; i < fZDD_IC_E.size(); i++){
    cout << "DetNbr: " << fZDD_IC_N[i] << " Energy: " << fZDD_IC_E[i] << " TimeStamp: " << fZDD_IC_TS[i] << endl;
  }
  cout << "ZDD_PLE Mult: " << fZDD_PL_E.size() << endl;
  for (size_t i = 0 ; i < fZDD_PL_E.size(); i++){
    cout << "DetNbr: " << fZDD_PL_N[i] << " Energy: " << fZDD_PL_E[i] << " TimeStamp: " << fZDD_PL_TS[i] << endl;
  }
  cout << "ZDD_DCE Mult: " << fZDD_DC_E.size() << endl;
  for (size_t i = 0 ; i < fZDD_DC_E.size(); i++){
    cout << "DetNbr: " << fZDD_DC_N[i] << " Energy: " << fZDD_DC_E[i] << " TimeStamp: " << fZDD_DC_TS[i] << endl;
  }
  cout << "ZDD_EXOE Mult: " << fZDD_EXO_E.size() << endl;
  for (size_t i = 0 ; i < fZDD_EXO_E.size(); i++){
    cout << "DetNbr: " << fZDD_EXO_N[i] << " Energy: " << fZDD_EXO_E[i] << " TimeStamp: " << fZDD_EXO_TS[i] << endl;
  }  
}
