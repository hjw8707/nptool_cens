/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien MATTA  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   :                                                         *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include "TCsITlData.h"


ClassImp(TCsITlData)

TCsITlData::TCsITlData()
{
}



TCsITlData::~TCsITlData()
{
}



void TCsITlData::Clear()
{
  fDetN.clear();
  fE.clear();
  fT.clear();
}



void TCsITlData::Dump() const
{
  std::cout << "========== Check CsITl Data ==============" << std::endl;
  std::cout << "  Total Size = " << GetMult() << std::endl;
  for (Int_t i = 0 ; i < GetMult() ; i++) {
    std::cout << " DetN = " << fDetN[i] << ", ";
    std::cout << " E = " << fE[i] << ", ";
    std::cout << " T = " << fT[i] << std::endl;
  }
  std::cout << "=======================================" << std::endl;}
