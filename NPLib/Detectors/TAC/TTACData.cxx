/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TAC Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TTACData.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

ClassImp(TTACData)

    //////////////////////////////////////////////////////////////////////
    TTACData::TTACData() {}

//////////////////////////////////////////////////////////////////////
TTACData::~TTACData() {}

//////////////////////////////////////////////////////////////////////
void TTACData::Clear() {
  fTAC_Time.clear();
  fTAC_N.clear();
  fTAC_Name.clear();
  fTAC_TS.clear();
}

//////////////////////////////////////////////////////////////////////
void TTACData::Dump() const {
  // // This method is very useful for debuging and worth the dev.
  // cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TTACData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // size_t mysize = fTAC_Channel.size();
  // cout << "TAC_Mult: " << mysize << endl;
  // for (size_t i = 0 ; i < mysize ; i++){
  //   cout << "TACNbr: " << fTAC_Nbr[i]  << " Channel: " << fTAC_Channel[i] << " TS: " << fTAC_TS[i] << endl;
  // }
}
