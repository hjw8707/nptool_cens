/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold STARK Raw data                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TSTARKRaw.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std; 

ClassImp(TSTARKRaw)

//////////////////////////////////////////////////////////////////////
TSTARKRaw::TSTARKRaw() {}

//////////////////////////////////////////////////////////////////////
TSTARKRaw::~TSTARKRaw() {}

//////////////////////////////////////////////////////////////////////
void TSTARKRaw::Clear() {
  // X6
  fX6_Det.clear();
  fX6_Ch.clear();
  fX6_E.clear();
  fX6_T.clear();}

void TSTARKRaw::SetRaw(const TSTARKData* sData){
  /*  for (Int_t i = 0 ; i < sData->GetX6_Mult() ; i++) {
    fX6_Det.push_back(sData->GetX6_DetN(i));
    fX6_Ch. push_back(2*(sData->GetX6_FStN(i) - 0)); // Front strip = 1 ~ 8 -> 0, 2, ..., 14
    fX6_E.  push_back(sData->GetX6_UpE(i));
    fX6_T.  push_back(sData->GetX6_T(i));

    fX6_Det.push_back(sData->GetX6_DetN(i));
    fX6_Ch. push_back(2*(sData->GetX6_FStN(i) - 0) + 1); // Front strip = 1 ~ 8 -> 1, 3, ..., 15
    fX6_E.  push_back(sData->GetX6_DwE(i));
    fX6_T.  push_back(sData->GetX6_T(i));

    fX6_Det.push_back(sData->GetX6_DetN(i));
    fX6_Ch. push_back(sData->GetX6_BStN(i) + 15); // Back strip = 1 ~ 4 -> 16 ~ 19
    fX6_E.  push_back(sData->GetX6_BkE(i));
    fX6_T.  push_back(sData->GetX6_T(i));}
  */}

//////////////////////////////////////////////////////////////////////
void TSTARKRaw::Dump() const { // to check the data
  std::cout << "========== Check X6 Raw Data ==============" << std::endl;
  std::cout << "  Total Size = " << GetX6_Mult() << std::endl;
  for (Int_t i = 0 ; i < GetX6_Mult() ; i++) {
    std::cout << " DetN = " << fX6_Det[i] << ", ";
    std::cout << " BStN = " << fX6_Ch[i] << ", ";
    std::cout << " E = " << fX6_E[i] << ", ";
    std::cout << " T = " << fX6_T[i] << std::endl;
  }
  std::cout << "=======================================" << std::endl;}
