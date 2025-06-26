/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Théodore Efremov contact address: theodore.efremov@cea.fr*
 *                                                                           *
 * Creation Date  : Nov 2024                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Time Raw data                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TTimeData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TTimeData)


//////////////////////////////////////////////////////////////////////
TTimeData::TTimeData() {
}



//////////////////////////////////////////////////////////////////////
TTimeData::~TTimeData() {
}



//////////////////////////////////////////////////////////////////////
void TTimeData::Clear() {
  fTS_MWPC13.clear();
  fTS_MWPC14.clear();
  fTS_MWPC23.clear();
  fTS_MWPC24.clear();

  fTime_MWPC13.clear();
  fTime_MWPC14.clear();
  fTime_MWPC23.clear();
  fTime_MWPC24.clear();

  fToff_DT13.clear();
  fToff_DT14.clear();

  fSection_MWPC3.clear();
  fSection_MWPC4.clear();

}



//////////////////////////////////////////////////////////////////////
void TTimeData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TTimeData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  size_t mysize = fTime_MWPC13.size();
  cout << "MWPC_Mult: " << mysize << endl;
  for(unsigned int i=0; i<mysize; i++){
    cout << "Time 13 mult " << i+1 << " / T= " << fTime_MWPC13[i] << endl;
    cout << "Time 14 mult " << i+1 << " / T= " << fTime_MWPC14[i] << endl;
    cout << "Time 23 mult " << i+1 << " / T= " << fTime_MWPC23[i] << endl;
    cout << "Time 24 mult " << i+1 << " / T= " << fTime_MWPC24[i] << endl;
  }
}
