#ifndef __ZDDCALDATA__
#define __ZDDCALDATA__
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <iostream>
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TZDDCalData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order
  // to allow multiplicity treatment
 private:
  // IC
  vector<float>   cZDD_IC_E;
  vector<uint16_t> cZDD_IC_N;
  vector<uint64_t> cZDD_IC_TS;

  // Plastic
  vector<float>   cZDD_PL_E;
  vector<uint16_t> cZDD_PL_N;
  vector<uint64_t> cZDD_PL_TS;

  // DC
  vector<float>   cZDD_DC_E;
  vector<uint16_t> cZDD_DC_N;
  vector<uint64_t> cZDD_DC_TS;

  // EXOZDD
  vector<float>   cZDD_EXO_E;
  vector<uint16_t> cZDD_EXO_N;
  vector<uint64_t> cZDD_EXO_TS;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
 public:
  TZDDCalData();
  ~TZDDCalData();

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
 public:
  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
 public:
  //////////////////////    SETTERS    ////////////////////////
  // Energy
  inline void SetZDDIC(const uint16_t& DetNumb, const float& Energy, const uint64_t& TimeStamp) {
    cZDD_IC_E.push_back(Energy);
    cZDD_IC_N.push_back(DetNumb);
    cZDD_IC_TS.push_back(TimeStamp);
  }; //!
  inline void SetZDDPL(const uint16_t& DetNumb, const float& Energy, const uint64_t& TimeStamp) {
    cZDD_PL_E.push_back(Energy);
    cZDD_PL_N.push_back(DetNumb);
    cZDD_PL_TS.push_back(TimeStamp);
  }; //!
  inline void SetZDDDC(const uint16_t& DetNumb, const float& Energy, const uint64_t& TimeStamp) {
    cZDD_DC_E.push_back(Energy);
    cZDD_DC_N.push_back(DetNumb);
    cZDD_DC_TS.push_back(TimeStamp);
  }; //!
  inline void SetZDDEXO(const uint16_t& DetNumb, const float& Energy, const uint64_t& TimeStamp) {
    cZDD_EXO_E.push_back(Energy);
    cZDD_EXO_N.push_back(DetNumb);
    cZDD_EXO_TS.push_back(TimeStamp);
  }; //!

  //////////////////////    GETTERS    ////////////////////////
  inline uint16_t	GetZDD_ICMult()		                const {return cZDD_IC_E.size();}
  inline float    GetZDD_ICE(const uint16_t& i) const { return cZDD_IC_E[i]; }
  inline uint16_t GetZDD_ICN(const uint16_t& i) const { return cZDD_IC_N[i]; }
  inline uint64_t GetZDD_ICTS(const uint16_t& i) const { return cZDD_IC_TS[i]; }
  
  inline uint16_t	GetZDD_PLMult()		                const {return cZDD_PL_E.size();}
  inline float    GetZDD_PLE(const uint16_t& i) const { return cZDD_PL_E[i]; }
  inline uint16_t GetZDD_PLN(const uint16_t& i) const { return cZDD_PL_N[i]; }
  inline uint64_t GetZDD_PLTS(const uint16_t& i) const { return cZDD_PL_TS[i]; }
  
  inline uint16_t	GetZDD_DCMult()		                const {return cZDD_DC_E.size();}
  inline float    GetZDD_DCE(const uint16_t& i) const { return cZDD_DC_E[i]; }
  inline uint16_t GetZDD_DCN(const uint16_t& i) const { return cZDD_DC_N[i]; }
  inline uint64_t GetZDD_DCTS(const uint16_t& i) const { return cZDD_DC_TS[i]; }
  
  inline uint16_t	GetZDD_EXOMult()		                const {return cZDD_EXO_E.size();}
  inline float    GetZDD_EXOE(const uint16_t& i) const { return cZDD_EXO_E[i]; }
  inline uint16_t GetZDD_EXON(const uint16_t& i) const { return cZDD_EXO_N[i]; }
  inline uint64_t GetZDD_EXOTS(const uint16_t& i) const { return cZDD_EXO_TS[i]; }

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TZDDCalData, 1) // ZDDData structure
};

#endif
