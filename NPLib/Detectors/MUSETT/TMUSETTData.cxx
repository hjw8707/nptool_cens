/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : June 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold MUSETT Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include "TMUSETTData.h"
#include <iostream>

#include "MUSETTMap.h"
#include "TMUSETTData.h"
ClassImp(TMUSETTData)
    ////////////////////////////////////////////////////////////////////////////////
    TMUSETTData::TMUSETTData() {
  // Init the correspondace table
  for (unsigned int i = 0; i < 128; i++) {
    fMUMU_MapX[i] = MUSETT_MAP::MapX[i];
    fMUMU_MapY[i] = MUSETT_MAP::MapY[i];
  }
}
////////////////////////////////////////////////////////////////////////////////
TMUSETTData::~TMUSETTData() {}
////////////////////////////////////////////////////////////////////////////////
void TMUSETTData::Clear() {
  fMUMU_TRIGGER = 137;
  fMUMU_DSSDXE_DetectorNbr.clear();
  fMUMU_DSSDXE_StripNbr.clear();
  fMUMU_DSSDXE_Energy.clear();



  fMUMU_DSSDXT_DetectorNbr.clear();
  fMUMU_DSSDXT_StripNbr.clear();
  fMUMU_DSSDXT_Time.clear();


  fMUMU_DSSDX_TimeStamp.clear();
  fMUMU_DSSDX_IsInterstrip.clear();


  fMUMU_DSSDYE_DetectorNbr.clear();
  fMUMU_DSSDYE_StripNbr.clear();
  fMUMU_DSSDYE_Energy.clear();


  fMUMU_DSSDYT_DetectorNbr.clear();
  fMUMU_DSSDYT_StripNbr.clear();
  fMUMU_DSSDYT_Time.clear();

  fMUMU_DSSDY_TimeStamp.clear();
  fMUMU_DSSDY_IsInterstrip.clear();
}
////////////////////////////////////////////////////////////////////////////////
void TMUSETTData::Dump() const {
  std::cout << "XXXXXXXXXXXXXXXXXXXXXXXX Hey MUSETT Event XXXXXXXXXXXXXXXXX" << std::endl;

  std::cout << "// First Layer " << std::endl;
  // (X,E)
  std::cout << " DSSDXE_Mult = " << fMUMU_DSSDXE_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fMUMU_DSSDXE_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fMUMU_DSSDXE_DetectorNbr[i] << " strip: " << fMUMU_DSSDXE_StripNbr[i]
         << " Energy: " << fMUMU_DSSDXE_Energy[i] << std::endl;
  // (X,T)
  std::cout << " DSSDXT_Mult = " << fMUMU_DSSDXT_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fMUMU_DSSDXT_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fMUMU_DSSDXT_DetectorNbr[i] << " DSSD: " << fMUMU_DSSDXT_StripNbr[i]
         << " Time: " << fMUMU_DSSDXT_Time[i] << std::endl;
  // (Y,E)
  std::cout << " DSSDYE_Mult = " << fMUMU_DSSDYE_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fMUMU_DSSDYE_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fMUMU_DSSDYE_DetectorNbr[i] << " DSSD: " << fMUMU_DSSDYE_StripNbr[i]
         << " Energy: " << fMUMU_DSSDYE_Energy[i] << std::endl;
  // (Y,T)
  std::cout << " DSSDYT_Mult = " << fMUMU_DSSDYT_DetectorNbr.size() << std::endl;
  for (UShort_t i = 0; i < fMUMU_DSSDYT_DetectorNbr.size(); i++)
    std::cout << "  DetNbr: " << fMUMU_DSSDYT_DetectorNbr[i] << " DSSD: " << fMUMU_DSSDYT_StripNbr[i]
         << " Time: " << fMUMU_DSSDYT_Time[i] << std::endl;
}
