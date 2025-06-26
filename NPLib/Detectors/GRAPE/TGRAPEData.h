#ifndef __GRAPEDATA__
#define __GRAPEDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jw.hwang@ibs.re.kr       *
 *                                                                           *
 * Creation Date  : February 2021                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the GRAPE  raw data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include<stdlib.h>
#include <vector>
#include <map>
using namespace std ;

// ROOT
#include "TObject.h"

class TGRAPEData : public TObject {
private:
  // GRAPE
  // Energy
  vector<UShort_t> fGRAPE_Ge_GRAPENbr;
  vector<UShort_t> fGRAPE_Ge_CrystalNbr;
  vector<UShort_t> fGRAPE_Ge_SegmentNbr;
  vector<Double_t> fGRAPE_Ge_Energy;
  vector<Double_t> fGRAPE_Ge_TimeCFD;
  vector<Double_t> fGRAPE_Ge_TimeLED;

public:
  TGRAPEData();
  virtual ~TGRAPEData();
  
  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;
  
  /////////////////////           SETTERS           ////////////////////////
  inline void SetGeGRAPENbr(const UShort_t &GeGRAPENbr){fGRAPE_Ge_GRAPENbr.push_back(GeGRAPENbr); }
  inline void SetGeCrystalNbr(const UShort_t &GeCrystalNbr){fGRAPE_Ge_CrystalNbr.push_back(GeCrystalNbr);}
  inline void SetGeSegmentNbr(const UShort_t &GeSegmentNbr){fGRAPE_Ge_SegmentNbr.push_back(GeSegmentNbr);}
  inline void SetGeEnergy(const Double_t &GeEnergy){fGRAPE_Ge_Energy.push_back(GeEnergy);}
  inline void SetGeTimeCFD(const Double_t &GeTimeCFD){fGRAPE_Ge_TimeCFD.push_back(GeTimeCFD);}
  inline void SetGeTimeLED(const Double_t &GeTimeLED){fGRAPE_Ge_TimeLED.push_back(GeTimeLED);}


  /////////////////////           GETTERS           ////////////////////////
  inline UShort_t GetGeGRAPENbr(const unsigned int &i)   {return fGRAPE_Ge_GRAPENbr[i]; }
  inline UShort_t GetGeCrystalNbr(const unsigned int &i)  {return fGRAPE_Ge_CrystalNbr[i]; }
  inline UShort_t GetGeSegmentNbr(const unsigned int &i)  {return fGRAPE_Ge_SegmentNbr[i]; }

  inline Double_t GetGeEnergy(const unsigned int &i)      {return fGRAPE_Ge_Energy[i];}
  inline Double_t GetGeTimeCFD(const unsigned int &i)     {return fGRAPE_Ge_TimeCFD[i];}
  inline Double_t GetGeTimeLED(const unsigned int &i)     {return fGRAPE_Ge_TimeLED[i];}

  inline unsigned int GetMultiplicityGe()  {return fGRAPE_Ge_GRAPENbr.size();}
  
  ClassDef(TGRAPEData,1)  // GRAPEData structure
};

#endif

















