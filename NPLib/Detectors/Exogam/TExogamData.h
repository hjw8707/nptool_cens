#ifndef __EXOGAMDATA__
#define __EXOGAMDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : march 2009                                               *
 * Last update    : july 2019                                                         *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Exogam Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: Added vectors for real energy/time values (double) (T.Goigoux CEA) *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include "TObject.h"
#include <map>
#include <vector>

class TExogamData : public TObject {

 public:
  TExogamData();
  ~TExogamData();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

 public:
  std::vector<uint16_t> fExo_Crystal;
  std::vector<uint16_t> fExo_E;
  std::vector<uint16_t> fExo_E_HG; // High gain x20
  std::vector<uint64_t> fExo_TS;
  std::vector<uint16_t> fExo_TDC;
  std::vector<uint16_t> fExo_BGO;
  std::vector<uint16_t> fExo_CsI;
  std::vector<uint16_t> fExo_Outer1;
  std::vector<uint16_t> fExo_Outer2;
  std::vector<uint16_t> fExo_Outer3;
  std::vector<uint16_t> fExo_Outer4;


  /////////////////////           SETTERS           ////////////////////////
  inline void SetExo(const uint16_t& Crystal,const uint16_t& Energy,
  const uint16_t& Energy_HG,const uint64_t& TS,const uint16_t& TDC,
  const uint16_t& BGO,const uint16_t& CsI,const uint16_t& Outer1,
  const uint16_t& Outer2,const uint16_t& Outer3,const uint16_t& Outer4) { 
  fExo_Crystal.push_back(Crystal);
  fExo_E.push_back(Energy);
  fExo_E_HG.push_back(Energy_HG);
  fExo_TS.push_back(TS);
  fExo_TDC.push_back(TDC);
  fExo_BGO.push_back(BGO);
  fExo_CsI.push_back(CsI);
  fExo_Outer1.push_back(Outer1);
  fExo_Outer2.push_back(Outer2);
  fExo_Outer3.push_back(Outer3);
  fExo_Outer4.push_back(Outer4);
  }
  /////////////////////           GETTERS           ////////////////////////
  inline uint16_t GetExoMult() { return fExo_Crystal.size(); }
  inline uint16_t GetExoCrystal(const uint16_t& i) const  { return fExo_Crystal[i]; }
  inline uint16_t GetExoE(const uint16_t& i) const  { return fExo_E[i]; }
  inline uint16_t GetExoEHG(const uint16_t& i) const  { return fExo_E_HG[i]; }
  inline uint64_t GetExoTS(const uint16_t& i) const  { return fExo_TS[i]; }
  inline uint16_t GetExoTDC(const uint16_t& i) const  { return fExo_TDC[i]; }
  inline uint16_t GetExoBGO(const uint16_t& i) const  { return fExo_BGO[i]; }
  inline uint16_t GetExoCsI(const uint16_t& i) const  { return fExo_CsI[i]; }
  inline uint16_t GetExoOuter1(const uint16_t& i) const { return fExo_Outer1[i]; }
  inline uint16_t GetExoOuter2(const uint16_t& i) const  { return fExo_Outer2[i]; }
  inline uint16_t GetExoOuter3(const uint16_t& i) const  { return fExo_Outer3[i]; }
  inline uint16_t GetExoOuter4(const uint16_t& i) const  { return fExo_Outer4[i]; }

  ClassDef(TExogamData, 1) // ExogamData structure
};
#endif
