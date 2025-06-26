#ifndef __TinaDATA__
#define __TinaDATA__
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
 *  This class hold Tina Raw data                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include <vector>
#include "TObject.h"

using namespace std;

class TTinaData : public TObject {
  // data members are held in vectors to allow multiplicity treatment
 private: 
  // initial particle KE
  double fKE;
  // TTT front
  vector<UShort_t>   fTTTfront_E_DetectorNbr;
  vector<UShort_t>   fTTTfront_E_StripNbr;
  vector<Double_t>   fTTTfront_Energy;
  vector<UShort_t>   fTTTfront_T_DetectorNbr;
  vector<UShort_t>   fTTTfront_T_StripNbr;
  vector<Double_t>   fTTTfront_Time;
  // TTT back
  vector<UShort_t>   fTTTback_E_DetectorNbr;
  vector<UShort_t>   fTTTback_E_StripNbr;
  vector<Double_t>   fTTTback_Energy;
  vector<UShort_t>   fTTTback_T_DetectorNbr;
  vector<UShort_t>   fTTTback_T_StripNbr;
  vector<Double_t>   fTTTback_Time;
  // Pad
  vector<UShort_t>   fPad_E_DetectorNbr;
  vector<UShort_t>   fPad_E_PadNbr;
  vector<Double_t>   fPad_Energy;
  vector<UShort_t>   fPad_T_DetectorNbr;
  vector<UShort_t>   fPad_T_PadNbr;
  vector<Double_t>   fPad_Time;
  // YY1 ring
  vector<UShort_t>   fYY1ring_E_DetectorNbr;
  vector<UShort_t>   fYY1ring_E_StripNbr;
  vector<Double_t>   fYY1ring_Energy;
  vector<UShort_t>   fYY1ring_T_DetectorNbr;
  vector<UShort_t>   fYY1ring_T_StripNbr;
  vector<Double_t>   fYY1ring_Time;
  // YY1 sector
  vector<UShort_t>   fYY1sector_E_DetectorNbr;
  vector<UShort_t>   fYY1sector_E_StripNbr;
  vector<Double_t>   fYY1sector_Energy;
  vector<UShort_t>   fYY1sector_T_DetectorNbr;
  vector<UShort_t>   fYY1sector_T_StripNbr;
  vector<Double_t>   fYY1sector_Time;
  // CsI
  vector<UShort_t>   fCsI_E_DetectorNbr;
  vector<Double_t>   fCsI_Energy;
  vector<UShort_t>   fCsI_T_DetectorNbr;
  vector<Double_t>   fCsI_Time;

  // constructor and destructor
 public: 
  TTinaData();
  ~TTinaData();
    
  // inherited from TObject and overriden to avoid warnings
 public:
  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;

  // setters & getters
  // prefer inline declaration to avoid unnecessary call of frequently used methods
  // add //! to avoid Root creating dictionary for the methods
 public:
  inline void SetKE(double ke) { fKE = ke; }
  // TTT front
  inline void SetTTTfrontEnergy(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy){
    fTTTfront_E_DetectorNbr.push_back(DetNbr);
    fTTTfront_E_StripNbr.push_back(StripNbr);
    fTTTfront_Energy.push_back(Energy);
  };//!
  inline void SetTTTfrontTime(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Time){
    fTTTfront_T_DetectorNbr.push_back(DetNbr);     
    fTTTfront_T_StripNbr.push_back(StripNbr);
    fTTTfront_Time.push_back(Time);
  };//!
    // TTT back
  inline void SetTTTbackEnergy(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy){
    fTTTback_E_DetectorNbr.push_back(DetNbr);
    fTTTback_E_StripNbr.push_back(StripNbr);
    fTTTback_Energy.push_back(Energy);
  };//!
  inline void SetTTTbackTime(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Time){
    fTTTback_T_DetectorNbr.push_back(DetNbr);     
    fTTTback_T_StripNbr.push_back(StripNbr);
    fTTTback_Time.push_back(Time);
  };//!
    // Pad
  inline void SetPadEnergy(const UShort_t& DetNbr,const UShort_t& PadNbr,const Double_t& Energy){
    fPad_E_DetectorNbr.push_back(DetNbr);
    fPad_E_PadNbr.push_back(PadNbr);
    fPad_Energy.push_back(Energy);
  };//!
  inline void SetPadTime(const UShort_t& DetNbr,const UShort_t& PadNbr,const Double_t& Time){
    fPad_T_DetectorNbr.push_back(DetNbr);     
    fPad_T_PadNbr.push_back(PadNbr);
    fPad_Time.push_back(Time);
  };//!
    // YY1 ring
  inline void SetYY1ringEnergy(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy){
    fYY1ring_E_DetectorNbr.push_back(DetNbr);
    fYY1ring_E_StripNbr.push_back(StripNbr);
    fYY1ring_Energy.push_back(Energy);
  };//!
  inline void SetYY1ringTime(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Time){
    fYY1ring_T_DetectorNbr.push_back(DetNbr);     
    fYY1ring_T_StripNbr.push_back(StripNbr);
    fYY1ring_Time.push_back(Time);
  };//!
    // YY1 sector
  inline void SetYY1sectorEnergy(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy){
    fYY1sector_E_DetectorNbr.push_back(DetNbr);
    fYY1sector_E_StripNbr.push_back(StripNbr);
    fYY1sector_Energy.push_back(Energy);
  };//!
  inline void SetYY1sectorTime(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Time){
    fYY1sector_T_DetectorNbr.push_back(DetNbr);     
    fYY1sector_T_StripNbr.push_back(StripNbr);
    fYY1sector_Time.push_back(Time);
  };//!
    // CsI
  inline void SetCsIEnergy(const UShort_t& DetNbr,const Double_t& Energy){
    fCsI_E_DetectorNbr.push_back(DetNbr);
    fCsI_Energy.push_back(Energy);
  };//!
  inline void SetCsITime(const UShort_t& DetNbr,const Double_t& Time){
    fCsI_T_DetectorNbr.push_back(DetNbr);     
    fCsI_Time.push_back(Time);
  };//!

  inline double GetKE() { return fKE; }
    // TTT front
  inline UShort_t GetTTTfront_E_Mult() const
  {return fTTTfront_E_DetectorNbr.size();}
  inline UShort_t GetTTTfront_E_DetNbr(const unsigned int &i) const 
  {return fTTTfront_E_DetectorNbr[i];}//!
  inline UShort_t GetTTTfront_E_StripNbr(const unsigned int &i) const 
  {return fTTTfront_E_StripNbr[i];}//!
  inline Double_t GetTTTfront_Energy(const unsigned int &i) const 
  {return fTTTfront_Energy[i];}//!
  inline UShort_t GetTTTfront_T_Mult() const
  {return fTTTfront_T_DetectorNbr.size();}
  inline UShort_t GetTTTfront_T_DetNbr(const unsigned int &i) const 
  {return fTTTfront_T_DetectorNbr[i];}//!
  inline UShort_t GetTTTfront_T_StripNbr(const unsigned int &i) const 
  {return fTTTfront_T_StripNbr[i];}//!
  inline Double_t GetTTTfront_Time(const unsigned int &i) const 
  {return fTTTfront_Time[i];}//!
  // TTT back
  inline UShort_t GetTTTback_E_Mult() const
  {return fTTTback_E_DetectorNbr.size();}
  inline UShort_t GetTTTback_E_DetNbr(const unsigned int &i) const 
  {return fTTTback_E_DetectorNbr[i];}//!
  inline UShort_t GetTTTback_E_StripNbr(const unsigned int &i) const 
  {return fTTTback_E_StripNbr[i];}//!
  inline Double_t GetTTTback_Energy(const unsigned int &i) const 
  {return fTTTback_Energy[i];}//!
  inline UShort_t GetTTTback_T_Mult() const
  {return fTTTback_T_DetectorNbr.size();}
  inline UShort_t GetTTTback_T_DetNbr(const unsigned int &i) const 
  {return fTTTback_T_DetectorNbr[i];}//!
  inline UShort_t GetTTTback_T_StripNbr(const unsigned int &i) const 
  {return fTTTback_T_StripNbr[i];}//!
  inline Double_t GetTTTback_Time(const unsigned int &i) const 
  {return fTTTback_Time[i];}//!
  // Pad
  inline UShort_t GetPad_E_Mult() const
  {return fPad_E_DetectorNbr.size();}
  inline UShort_t GetPad_E_DetNbr(const unsigned int &i) const 
  {return fPad_E_DetectorNbr[i];}//!
  inline UShort_t GetPad_E_PadNbr(const unsigned int &i) const 
  {return fPad_E_PadNbr[i];}//!
  inline Double_t GetPad_Energy(const unsigned int &i) const 
  {return fPad_Energy[i];}//!
  inline UShort_t GetPad_T_Mult() const
  {return fPad_T_DetectorNbr.size();}
  inline UShort_t GetPad_T_DetNbr(const unsigned int &i) const 
  {return fPad_T_DetectorNbr[i];}//!
  inline UShort_t GetPad_T_PadNbr(const unsigned int &i) const 
  {return fPad_T_PadNbr[i];}//!
  inline Double_t GetPad_Time(const unsigned int &i) const 
  {return fPad_Time[i];}//!
  // YY1 ring
  inline UShort_t GetYY1ring_E_Mult() const
  {return fYY1ring_E_DetectorNbr.size();}
  inline UShort_t GetYY1ring_E_DetNbr(const unsigned int &i) const 
  {return fYY1ring_E_DetectorNbr[i];}//!
  inline UShort_t GetYY1ring_E_StripNbr(const unsigned int &i) const 
  {return fYY1ring_E_StripNbr[i];}//!
  inline Double_t GetYY1ring_Energy(const unsigned int &i) const 
  {return fYY1ring_Energy[i];}//!
  inline UShort_t GetYY1ring_T_Mult() const
  {return fYY1ring_T_DetectorNbr.size();}
  inline UShort_t GetYY1ring_T_DetNbr(const unsigned int &i) const 
  {return fYY1ring_T_DetectorNbr[i];}//!
  inline UShort_t GetYY1ring_T_StripNbr(const unsigned int &i) const 
  {return fYY1ring_T_StripNbr[i];}//!
  inline Double_t GetYY1ring_Time(const unsigned int &i) const 
  {return fYY1ring_Time[i];}//!
  // YY1 sector
  inline UShort_t GetYY1sector_E_Mult() const
  {return fYY1sector_E_DetectorNbr.size();}
  inline UShort_t GetYY1sector_E_DetNbr(const unsigned int &i) const 
  {return fYY1sector_E_DetectorNbr[i];}//!
  inline UShort_t GetYY1sector_E_StripNbr(const unsigned int &i) const 
  {return fYY1sector_E_StripNbr[i];}//!
  inline Double_t GetYY1sector_Energy(const unsigned int &i) const 
  {return fYY1sector_Energy[i];}//!
  inline UShort_t GetYY1sector_T_Mult() const
  {return fYY1sector_T_DetectorNbr.size();}
  inline UShort_t GetYY1sector_T_DetNbr(const unsigned int &i) const 
  {return fYY1sector_T_DetectorNbr[i];}//!
  inline UShort_t GetYY1sector_T_StripNbr(const unsigned int &i) const 
  {return fYY1sector_T_StripNbr[i];}//!
  inline Double_t GetYY1sector_Time(const unsigned int &i) const 
  {return fYY1sector_Time[i];}//!
  // CsI
  inline UShort_t GetCsI_E_Mult() const
  {return fCsI_E_DetectorNbr.size();}
  inline UShort_t GetCsI_E_DetNbr(const unsigned int &i) const 
  {return fCsI_E_DetectorNbr[i];}//!
  inline Double_t GetCsI_Energy(const unsigned int &i) const 
  {return fCsI_Energy[i];}//!
  inline UShort_t GetCsI_T_Mult() const
  {return fCsI_T_DetectorNbr.size();}
  inline UShort_t GetCsI_T_DetNbr(const unsigned int &i) const 
  {return fCsI_T_DetectorNbr[i];}//!
  inline Double_t GetCsI_Time(const unsigned int &i) const 
  {return fCsI_Time[i];}//!

  // Required for Root dictionary
  ClassDef(TTinaData,1) 
    };

#endif
