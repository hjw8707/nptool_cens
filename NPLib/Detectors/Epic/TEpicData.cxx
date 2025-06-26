/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon  contact address: audrey.chatillon@cea.fr                        *
 *                                                                           *
 * Creation Date  : décembre 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Epic Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TEpicData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TEpicData)

//////////////////////////////////////////////////////////////////////
TEpicData::TEpicData() {
}
//////////////////////////////////////////////////////////////////////
TEpicData::~TEpicData() {
}

//////////////////////////////////////////////////////////////////////
void TEpicData::Clear() {
  fEpic_Data.clear();
}

//////////////////////////////////////////////////////////////////////
void TEpicData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TEpicData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  size_t mysize = GetMultiplicity();
  cout << "MultAnodes = : " << mysize << endl;
 
  for (size_t i = 0 ; i < mysize ; i++){
    cout << "==================================" << endl;
    cout << "ANODE #" << GetAnodeNbr(i) << " Q1: " << GetQ1(i) << ", Qmax = " << GetQmax(i) << endl;
    cout << "POSITION : X = " << GetXpos(i) << ", Y =  " << GetYpos(i) << ", X = " << GetZpos(i) << endl;
  }
}
