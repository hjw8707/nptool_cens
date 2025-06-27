#ifndef __STARKjrRAW__
#define __STARKjrRAW__
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
 *  This class hold STARKjr Raw data                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include <vector>
#include "TSTARKjrData.h"
#include "TObject.h"

using namespace std;

class TSTARKjrRaw : public TObject {
  // data members are held in vectors to allow multiplicity treatment
 private: 
  // X6
  vector<Int_t> fX6_Det; // Detector number
  vector<Int_t> fX6_Ch;  // front = 0 ~ 15, back = 16 ~ 19
  vector<Double_t> fX6_E; // Energy
  vector<Double_t> fX6_T; // Time

  // constructor and destructor
 public: 
  TSTARKjrRaw();
  ~TSTARKjrRaw();
  
  // inherited from TObject and overriden to avoid warnings
 public:
  void Clear();
  void Clear(const Option_t*) { Clear(); };
  void Dump() const;

  // setters & getters
  // prefer inline declaration to avoid unnecessary call of frequently used methods
  // add //! to avoid Root creating dictionary for the methods
 public:
  // X6
  void SetRaw(const TSTARKjrData* sData);

  // X6
  inline Int_t GetX6_Mult()        const {return fX6_Det.size();}
  inline Int_t GetX6_Det(Int_t i)  const {return fX6_Det[i]; }
  inline Int_t GetX6_Ch(Int_t i)   const {return fX6_Ch[i]; }
  inline Double_t GetX6_E(Int_t i) const { return fX6_E[i]; }
  inline Double_t GetX6_T(Int_t i) const { return fX6_T[i]; }

  // Required for Root dictionary
  ClassDef(TSTARKjrRaw,1) 
    };

#endif
