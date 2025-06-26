#ifndef __SuperX3DATA__
#define __SuperX3DATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : january 2011                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *     This class holds the raw data storage for the SuperX3 detector from Micron *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ headers
#include <vector>
using namespace std;

// ROOT headers
#include "TObject.h"

class TSuperX3Data : public TObject {
 private:
  // Front
  // Up
  // Energy
  vector<UShort_t> fSuperX3_UpE_DetectorNbr;
  vector<UShort_t> fSuperX3_UpE_StripNbr;
  vector<Double_t> fSuperX3_UpE_Energy;
  // Time
  vector<UShort_t> fSuperX3_UpT_DetectorNbr;
  vector<UShort_t> fSuperX3_UpT_StripNbr;
  vector<Double_t> fSuperX3_UpT_Time;
  // Energy
  vector<UShort_t> fSuperX3_DownE_DetectorNbr;
  vector<UShort_t> fSuperX3_DownE_StripNbr;
  vector<Double_t> fSuperX3_DownE_Energy;
  // Time
  vector<UShort_t> fSuperX3_DownT_DetectorNbr;
  vector<UShort_t> fSuperX3_DownT_StripNbr;
  vector<Double_t> fSuperX3_DownT_Time;

  // Back
  // Energy
  vector<Double_t> fSuperX3_BackE_Energy;
  vector<UShort_t> fSuperX3_BackE_StripNbr;
  vector<UShort_t> fSuperX3_BackE_DetectorNbr;
  // Time
  vector<Double_t> fSuperX3_BackT_Time;
  vector<UShort_t> fSuperX3_BackT_StripNbr;
  vector<UShort_t> fSuperX3_BackT_DetectorNbr;

 public:
  TSuperX3Data();
  virtual ~TSuperX3Data();

  void Clear();
  void Clear(const Option_t*){};
  void Dump() const;

  /////////////////////           SETTERS           ////////////////////////
  // up
  void SetUpE(UShort_t DetNbr, UShort_t StripNbr, Double_t Energy) {
    fSuperX3_UpE_DetectorNbr.push_back(DetNbr);
    fSuperX3_UpE_StripNbr.push_back(StripNbr);
    fSuperX3_UpE_Energy.push_back(Energy);
  }
  void SetUpT(UShort_t DetNbr, UShort_t StripNbr, Double_t Time) {
    fSuperX3_UpT_DetectorNbr.push_back(DetNbr);
    fSuperX3_UpT_StripNbr.push_back(StripNbr);
    fSuperX3_UpT_Time.push_back(Time);
  }
  // down
  void SetDownE(UShort_t DetNbr, UShort_t StripNbr, Double_t Energy) {
    fSuperX3_DownE_DetectorNbr.push_back(DetNbr);
    fSuperX3_DownE_StripNbr.push_back(StripNbr);
    fSuperX3_DownE_Energy.push_back(Energy);
  }
  void SetDownT(UShort_t DetNbr, UShort_t StripNbr, Double_t Time) {
    fSuperX3_DownT_DetectorNbr.push_back(DetNbr);
    fSuperX3_DownT_StripNbr.push_back(StripNbr);
    fSuperX3_DownT_Time.push_back(Time);
  }

  // Back E
  void SetBackE(UShort_t DetNbr, UShort_t StripNbr, Double_t Energy) {
    fSuperX3_BackE_DetectorNbr.push_back(DetNbr);
    fSuperX3_BackE_StripNbr.push_back(StripNbr);
    fSuperX3_BackE_Energy.push_back(Energy);
  }
  // Back T
  void SetBackT(UShort_t DetNbr, UShort_t StripNbr, Double_t Time) {
    fSuperX3_BackT_DetectorNbr.push_back(DetNbr);
    fSuperX3_BackT_StripNbr.push_back(StripNbr);
    fSuperX3_BackT_Time.push_back(Time);
  }

  /////////////////////           GETTERS           ////////////////////////
  // DSSD
  // (Up, E)
  UShort_t GetUpEMult() const { return fSuperX3_UpE_DetectorNbr.size(); }
  UShort_t GetUpEDetectorNbr(Int_t i) const { return fSuperX3_UpE_DetectorNbr[i]; }
  UShort_t GetUpEStripNbr(Int_t i) const { return fSuperX3_UpE_StripNbr[i]; }
  Double_t GetUpEEnergy(Int_t i) const { return fSuperX3_UpE_Energy[i]; }
  // (Up, T)
  UShort_t GetUpTMult() const { return fSuperX3_UpT_DetectorNbr.size(); }
  UShort_t GetUpTDetectorNbr(Int_t i) const { return fSuperX3_UpT_DetectorNbr[i]; }
  UShort_t GetUpTStripNbr(Int_t i) const { return fSuperX3_UpT_StripNbr[i]; }
  Double_t GetUpTTime(Int_t i) const { return fSuperX3_UpT_Time[i]; }
  // (Down, E)
  UShort_t GetDownEMult() const { return fSuperX3_DownE_DetectorNbr.size(); }
  UShort_t GetDownEDetectorNbr(Int_t i) const { return fSuperX3_DownE_DetectorNbr[i]; }
  UShort_t GetDownEStripNbr(Int_t i) const { return fSuperX3_DownE_StripNbr[i]; }
  Double_t GetDownEEnergy(Int_t i) const { return fSuperX3_DownE_Energy[i]; }
  // (Down, T)
  UShort_t GetDownTMult() const { return fSuperX3_DownT_DetectorNbr.size(); }
  UShort_t GetDownTDetectorNbr(Int_t i) const { return fSuperX3_DownT_DetectorNbr[i]; }
  UShort_t GetDownTStripNbr(Int_t i) const { return fSuperX3_DownT_StripNbr[i]; }
  Double_t GetDownTTime(Int_t i) const { return fSuperX3_DownT_Time[i]; }

  // (Back, E)
  UShort_t GetBackEMult() const { return fSuperX3_BackE_DetectorNbr.size(); }
  UShort_t GetBackEDetectorNbr(Int_t i) const { return fSuperX3_BackE_DetectorNbr[i]; }
  UShort_t GetBackEStripNbr(Int_t i) const { return fSuperX3_BackE_StripNbr[i]; }
  Double_t GetBackEEnergy(Int_t i) const { return fSuperX3_BackE_Energy[i]; }
  // (Back, T)
  UShort_t GetBackTMult() const { return fSuperX3_BackT_DetectorNbr.size(); }
  UShort_t GetBackTDetectorNbr(Int_t i) const { return fSuperX3_BackT_DetectorNbr[i]; }
  UShort_t GetBackTStripNbr(Int_t i) const { return fSuperX3_BackT_StripNbr[i]; }
  Double_t GetBackTTime(Int_t i) const { return fSuperX3_BackT_Time[i]; }

  ClassDef(TSuperX3Data, 1) // TSuperX3Data raw data
};

#endif
