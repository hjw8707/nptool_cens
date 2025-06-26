/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jwhwang@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 4월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CACAO Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/
#include "TCACAOData.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std; 

ClassImp(TCACAOData)


//////////////////////////////////////////////////////////////////////
TCACAOData::TCACAOData() {
}

//////////////////////////////////////////////////////////////////////
TCACAOData::~TCACAOData() {
}

//////////////////////////////////////////////////////////////////////
void TCACAOData::Clear()
{
  fDetN.clear();
  fE.clear();
  fT.clear();
}

//////////////////////////////////////////////////////////////////////
void TCACAOData::Dump() const
{
  std::cout << "========== Check CACAO Data ==============" << std::endl;
  std::cout << "  Total Size = " << GetMult() << std::endl;
  for (Int_t i = 0 ; i < GetMult() ; i++) {
    std::cout << " DetN = " << fDetN[i] << ", ";
    std::cout << " E = " << fE[i] << ", ";
    std::cout << " T = " << fT[i] << std::endl;
  }
  std::cout << "=======================================" << std::endl;}
