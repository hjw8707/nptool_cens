/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FPMW Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TFPMWData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TFPMWData)


//////////////////////////////////////////////////////////////////////
TFPMWData::TFPMWData() {
}



//////////////////////////////////////////////////////////////////////
TFPMWData::~TFPMWData() {
}



//////////////////////////////////////////////////////////////////////
void TFPMWData::Clear() {
  // X //
  fFPMW_DetX.clear();
  fFPMW_StripX.clear();
  fFPMW_ChargeX.clear();
 
  // Y //
  fFPMW_DetY.clear();
  fFPMW_StripY.clear();
  fFPMW_ChargeY.clear();
 
}



//////////////////////////////////////////////////////////////////////
void TFPMWData::Dump() const {
  // This method is very useful for debuging and worth the dev.
  cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event [TFPMWData::Dump()] XXXXXXXXXXXXXXXXX" << endl;

  // X
  size_t mysizeX = fFPMW_DetX.size();
  cout << "FPMW_MultX: " << mysizeX << endl;

  // Y
  size_t mysizeY = fFPMW_DetY.size();
  cout << "FPMW_MultY: " << mysizeY << endl;

}
