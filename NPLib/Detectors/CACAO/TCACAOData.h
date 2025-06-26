#ifndef __CACAODATA__
#define __CACAODATA__
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

// STL
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TCACAOData : public TObject {
 private: 
  vector<Int_t> fDetN;
  vector<Double_t> fE;
  vector<Double_t> fT;

 public:
   TCACAOData();
   virtual ~TCACAOData();

   void   Clear();
   void   Clear(const Option_t*) {};
   void   Dump() const;

   //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
 public:
  //////////////////////    SETTERS    ////////////////////////
  // Energy
  inline void Set(const Int_t DetN,
		  const Double_t E, const Double_t T) {
    fDetN.push_back(DetN);
    fE.push_back(E);
    fT.push_back(T);}
  
  // all
  inline Int_t GetMult() const {return fDetN.size();}
  inline Int_t GetDetN(Int_t i) const {return fDetN[i]; }
  inline Double_t GetE(Int_t i) const { return fE[i]; }
  inline Double_t GetT(Int_t i)   const { return fT[i]; }

   ClassDef(TCACAOData,1)  // CACAOData structure
};

#endif
