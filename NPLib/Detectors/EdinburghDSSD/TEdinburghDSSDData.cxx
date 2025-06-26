/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Louis Heitz  contact address: louis.heitz@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : mars 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold EdinburghDSSD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TEdinburghDSSDData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

ClassImp(TEdinburghDSSDData)


//////////////////////////////////////////////////////////////////////
TEdinburghDSSDData::TEdinburghDSSDData() {
}



//////////////////////////////////////////////////////////////////////
TEdinburghDSSDData::~TEdinburghDSSDData() {
}



//////////////////////////////////////////////////////////////////////
void TEdinburghDSSDData::Clear() {

  fEdin_DSSDXE_DetectorNbr.clear();
  fEdin_DSSDXE_StripNbr.clear();
  fEdin_DSSDXE_Energy.clear();



  fEdin_DSSDXT_DetectorNbr.clear();
  fEdin_DSSDXT_StripNbr.clear();
  fEdin_DSSDXT_Time.clear();


  fEdin_DSSDX_TimeStamp.clear();
  fEdin_DSSDX_IsInterstrip.clear();


  fEdin_DSSDYE_DetectorNbr.clear();
  fEdin_DSSDYE_StripNbr.clear();
  fEdin_DSSDYE_Energy.clear();


  fEdin_DSSDYT_DetectorNbr.clear();
  fEdin_DSSDYT_StripNbr.clear();
  fEdin_DSSDYT_Time.clear();

  fEdin_DSSDY_TimeStamp.clear();
  fEdin_DSSDY_IsInterstrip.clear();
}



//////////////////////////////////////////////////////////////////////
void TEdinburghDSSDData::Dump() const {
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXX MUSETT Event XXXXXXXXXXXXXXXXX" << std::endl;

  std::cout << "// First Layer " << std::endl;
  // (X,E)
  std::cout << " DSSDXE_Mult = " << fEdin_DSSDXE_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fEdin_DSSDXE_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fEdin_DSSDXE_DetectorNbr[i] << " strip: " << fEdin_DSSDXE_StripNbr[i]
         << " Energy: " << fEdin_DSSDXE_Energy[i] << std::endl;
  // (X,T)
  std::cout << " DSSDXT_Mult = " << fEdin_DSSDXT_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fEdin_DSSDXT_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fEdin_DSSDXT_DetectorNbr[i] << " DSSD: " << fEdin_DSSDXT_StripNbr[i]
         << " Time: " << fEdin_DSSDXT_Time[i] << std::endl;
  // (Y,E)
  std::cout << " DSSDYE_Mult = " << fEdin_DSSDYE_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fEdin_DSSDYE_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fEdin_DSSDYE_DetectorNbr[i] << " DSSD: " << fEdin_DSSDYE_StripNbr[i]
         << " Energy: " << fEdin_DSSDYE_Energy[i] << std::endl;
  // (Y,T)
  std::cout << " DSSDYT_Mult = " << fEdin_DSSDYT_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fEdin_DSSDYT_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fEdin_DSSDYT_DetectorNbr[i] << " DSSD: " << fEdin_DSSDYT_StripNbr[i]
         << " Time: " << fEdin_DSSDYT_Time[i] << std::endl;
}
