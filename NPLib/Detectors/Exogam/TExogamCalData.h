#ifndef __EXOGAMCALDATA__
#define __EXOGAMCALDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2024                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Exogam Cal data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include "TObject.h"
#include <map>
#include <vector>


class TExogamCalData : public TObject {

 public:
  TExogamCalData();
  ~TExogamCalData();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

 public:
  std::vector<uint16_t> cExo_Crystal;
  std::vector<float> cExo_E;
  std::vector<float> cExo_E_HG; // High gain x20
  std::vector<uint64_t> cExo_TS;
  std::vector<float> cExo_TDC;
  std::vector<float> cExo_BGO;
  std::vector<float> cExo_CsI;
  std::vector<float> cExo_Outer1;
  std::vector<float> cExo_Outer2;
  std::vector<float> cExo_Outer3;
  std::vector<float> cExo_Outer4;


  /////////////////////           SETTERS           ////////////////////////
  inline void SetExo(const uint16_t& Crystal,const float& Energy,
  const float& Energy_HG,const uint64_t& TS,const float& TDC,
  const float& BGO,const float& CsI,const float& Outer1,
  const float& Outer2,const float& Outer3,const float& Outer4) { 
  cExo_Crystal.push_back(Crystal);
  cExo_E.push_back(Energy);
  cExo_E_HG.push_back(Energy_HG);
  cExo_TS.push_back(TS);
  cExo_TDC.push_back(TDC);
  cExo_BGO.push_back(BGO);
  cExo_CsI.push_back(CsI);
  cExo_Outer1.push_back(Outer1);
  cExo_Outer2.push_back(Outer2);
  cExo_Outer3.push_back(Outer3);
  cExo_Outer4.push_back(Outer4);
  }
  /////////////////////           GETTERS           ////////////////////////
  inline uint16_t GetExoMult() { return cExo_Crystal.size(); }
  inline uint16_t GetExoCrystal(const uint16_t& i) const  { return cExo_Crystal[i]; }
  inline float GetExoE(const uint16_t& i) const  { return cExo_E[i]; }
  inline float GetExoEHG(const uint16_t& i) const  { return cExo_E_HG[i]; }
  inline uint64_t GetExoTS(const uint16_t& i) const  { return cExo_TS[i]; }
  inline float GetExoTDC(const uint16_t& i) const  { return cExo_TDC[i]; }
  inline float GetExoBGO(const uint16_t& i) const  { return cExo_BGO[i]; }
  inline float GetExoCsI(const uint16_t& i) const  { return cExo_CsI[i]; }
  inline float GetExoOuter1(const uint16_t& i) const { return cExo_Outer1[i]; }
  inline float GetExoOuter2(const uint16_t& i) const  { return cExo_Outer2[i]; }
  inline float GetExoOuter3(const uint16_t& i) const  { return cExo_Outer3[i]; }
  inline float GetExoOuter4(const uint16_t& i) const  { return cExo_Outer4[i]; }

  ClassDef(TExogamCalData, 1) // ExogamData structure
};

#endif
