#ifndef __CSITLDATA__
#define __CSITLDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author:    contact address:                                      *
 *                                                                           *
 * Creation Date  :                                                          *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include <vector>

#include "TObject.h"
using namespace std ;


class TCsITlData : public TObject {
 private: 
  vector<Int_t> fDetN;
  vector<Double_t> fE;
  vector<Double_t> fT;

 public:
   TCsITlData();
   virtual ~TCsITlData();

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

   ClassDef(TCsITlData,1)  // CsITlData structure
};

#endif
