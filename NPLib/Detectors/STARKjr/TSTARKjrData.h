#ifndef __STARKjrDATA__
#define __STARKjrDATA__
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
#include <string>
#include "TObject.h"
#include "TVector3.h"

using namespace std;

class TSTARKjrData : public TObject {
  // data members are held in vectors to allow multiplicity treatment
 private: 
  // For all detector
  vector<Int_t>  fType; // 0 for X6, 1 for BB10, 2 for QQQ5
  vector<Int_t>  fDetN;
  vector<Int_t>  fFStN;
  vector<Int_t>  fBStN;
  vector<Double_t> fFrE;
  vector<Double_t> fBkE;
  vector<Double_t> fUpE; // for X6 only
  vector<Double_t> fDwE; // for X6 only
  vector<Double_t> fT;
  vector<TVector3> fPos; // relative to each detector

  
  // constructor and destructor
 public: 
  TSTARKjrData();
  ~TSTARKjrData();
  
  // inherited from TObject and overriden to avoid warnings
 public:
  void Clear();
  void Clear(const Option_t*) { Clear(); };
  void Dump() const;

  // setters & getters
  // prefer inline declaration to avoid unnecessary call of frequently used methods
  // add //! to avoid Root creating dictionary for the methods
 public:
  // for all
  inline void Set(const Int_t Type,
		  const Int_t DetN, const Int_t FStN, const Int_t BStN,
		  const Double_t FrE, const Double_t BkE,
		  const Double_t UpE, const Double_t DwE,
		  const Double_t T,
		  const Double_t posX, const Double_t posY,  const Double_t posZ) {
    fType.push_back(Type);
    fDetN.push_back(DetN);
    fFStN.push_back(FStN);
    fBStN.push_back(BStN);
    fFrE.push_back(FrE);
    fBkE.push_back(BkE);
    fUpE.push_back(UpE);
    fDwE.push_back(DwE);
    fT.push_back(T);
    fPos.push_back(TVector3(posX, posY, posZ));};

  // all
  inline Int_t GetMult() const {return fDetN.size();}
  inline Int_t GetType(Int_t i) const {return fType[i]; }
  inline Int_t GetDetN(Int_t i) const {return fDetN[i]; }
  inline Int_t GetFStN(Int_t i) const {return fFStN[i]; }
  inline Int_t GetBStN(Int_t i) const {return fBStN[i]; }
  inline Double_t GetFrE(Int_t i) const { return fFrE[i]; }
  inline Double_t GetBkE(Int_t i) const { return fBkE[i]; }
  inline Double_t GetUpE(Int_t i) const { return fUpE[i]; }
  inline Double_t GetDwE(Int_t i) const { return fDwE[i]; }
  inline Double_t GetT(Int_t i)   const { return fT[i]; }
  inline TVector3 GetPos(Int_t i)   const { return fPos[i]; }

  // Required for Root dictionary
  ClassDef(TSTARKjrData,1) 
    };

#endif
