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

#include "TSTARKData.h"
#include <iostream>

using namespace std;

ClassImp(TSTARKData)

    //////////////////////////////////////////////////////////////////////
    TSTARKData::TSTARKData() {}

//////////////////////////////////////////////////////////////////////
TSTARKData::~TSTARKData() {}

//////////////////////////////////////////////////////////////////////
void TSTARKData::Clear() {
  // X6
  fType.clear();
  fDetN.clear();
  fFStN.clear();
  fBStN.clear();
  fFrE.clear();
  fBkE.clear();
  fUpE.clear();
  fDwE.clear();
  fT.clear();
  fPos.clear();
  fCsIN.clear();
  fCsIE.clear();
}

//////////////////////////////////////////////////////////////////////
void TSTARKData::Print() const {
  std::cout << "========== Print STARK Data ==============" << std::endl;
  std::cout << "  Total Size = " << GetMult() << std::endl;
  for (Int_t i = 0; i < GetMult(); i++) {
    std::cout << " Type = " << fType[i] << ", ";
    std::cout << " DetN = " << fDetN[i] << ", ";
    std::cout << " HasCsI = " << (HasCsI(i) ? "Yes" : "No") << std::endl;
  }
  std::cout << "=======================================" << std::endl;
}

//////////////////////////////////////////////////////////////////////
void TSTARKData::Dump() const { // to check the data
  std::cout << "========== Check STARK Data ==============" << std::endl;
  std::cout << "  Total Size = " << GetMult() << std::endl;
  for (Int_t i = 0; i < GetMult(); i++) {
    std::cout << " Type = " << fType[i] << ", ";
    std::cout << " DetN = " << fDetN[i] << ", ";
    std::cout << " FStN = " << fFStN[i] << ", ";
    std::cout << " BStN = " << fBStN[i] << ", ";
    std::cout << " FrE = " << fFrE[i] << ", ";
    std::cout << " BkE = " << fBkE[i] << ", ";
    std::cout << " UpE = " << fUpE[i] << ", ";
    std::cout << " DwE = " << fDwE[i] << ", ";
    std::cout << " T = " << fT[i] << ", ";
    std::cout << " CsIN = ";
    for (Int_t j = 0; j < GetCsINArray(i).size(); j++) {
      std::cout << GetCsIN(i, j) << " ";
    }
    std::cout << ", ";
    std::cout << " CsIE = ";
    for (Int_t j = 0; j < GetCsIEArray(i).size(); j++) {
      std::cout << GetCsIE(i, j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "=======================================" << std::endl;
}

//////////////////////////////////////////////////////////////////////
void TSTARKData::SetCsI(const Int_t DetN, const Int_t numCsI, const Double_t E) {
  // fDetN에서 DetN이 같은 인덱스 찾기
  auto it = std::find(fDetN.begin(), fDetN.end(), DetN);
  if (it != fDetN.end()) {
    int idx = std::distance(fDetN.begin(), it);
    fCsIN[idx].push_back(numCsI);
    fCsIE[idx].push_back(E);
    return;
  }
  // 없을리가 없다. 에러 처리
  throw std::runtime_error("Detector not found in TSTARKData::SetCsI");
}