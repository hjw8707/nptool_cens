#ifndef __ASGARDDATA__
#define __ASGARDDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the ASGARD  raw data (Made for ASG10 card)              *
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
#include "TVector3.h"
class TASGARDData : public TObject {
private:
  // ASGARD
  // crystal-level data
  vector<UInt_t> fASGARD_CloverNbr;
  vector<UInt_t> fASGARD_CrystalNbr;
  vector<UInt_t> fASGARD_SegmentNbr; // highest-energy segment
  vector<Double_t> fASGARD_Energy; // summed for crystal
  vector<Double_t> fASGARD_MaxEnergySegment; // summed for crystal
  vector<Double_t> fASGARD_TimeLED;
  // from first interacton point, radians
  vector<Double_t> fASGARD_Theta_SemiTrue;  //info[5]
  vector<Double_t> fASGARD_Phi_SemiTrue; //info[6]
  vector<TVector3> fASGARD_position;
  

  // segment-level data
  vector<UInt_t> fASGARD_CloverNbr_sub;
  vector<UInt_t> fASGARD_CrystalNbr_sub;
  vector<UInt_t> fASGARD_SegmentNbr_sub; // individual segments
  vector<Double_t> fASGARD_Energy_sub; // individual energies
  vector<Double_t> fASGARD_TimeLED_sub;
  vector<Double_t> fASGARD_Theta_SemiTrue_sub;
  vector<Double_t> fASGARD_Phi_SemiTrue_sub;
  vector<TVector3> fASGARD_position_sub;
  

  // vector<UInt_t> fASG_BGO_CloverNbr;
  // vector<UInt_t> fASG_BGO_CrystalNbr;
  // vector<UInt_t> fASG_BGO_PmNbr;
  // vector<Double_t> fASG_BGO_Energy;
  // vector<Double_t> fASG_BGO_TimeCFD;
  // vector<Double_t> fASG_BGO_TimeLED;
  int last_clover = -1;
  int last_crystal = -1;
  

public:
  TASGARDData();
  virtual ~TASGARDData();
  
  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;

  inline void ResetLastCloverCrystal(){last_clover= -1; last_crystal = -1;}
  /////////////////////           SETTERS           ////////////////////////
  inline void AddGeCloverNbr(const UInt_t &GeCloverNbr){fASGARD_CloverNbr_sub.push_back(GeCloverNbr); }
  inline void AddGeCrystalNbr(const UInt_t &GeCrystalNbr){fASGARD_CrystalNbr_sub.push_back(GeCrystalNbr);}
  inline void AddGeSegmentNbr(const UInt_t &GeSegmentNbr){fASGARD_SegmentNbr_sub.push_back(GeSegmentNbr);}
  inline void AddGeEnergy(const Double_t &GeEnergy){fASGARD_Energy_sub.push_back(GeEnergy);}
  inline void AddGeTimeLED(const Double_t &GeTimeLED){fASGARD_TimeLED_sub.push_back(GeTimeLED);}
  inline void AddGePosition(const TVector3 &GePosition){fASGARD_position_sub.push_back(GePosition);}
  inline void AddGeThetaSemiTrue(const Double_t &GeThetaSemiTrue){fASGARD_Theta_SemiTrue_sub.push_back(GeThetaSemiTrue);}
  inline void AddGePhiSemiTrue(const Double_t &GePhiSemiTrue){fASGARD_Phi_SemiTrue_sub.push_back(GePhiSemiTrue);}
  
  // inline void SetGeCloverNbr(const UInt_t &GeCloverNbr){fASGARD_CloverNbr.push_back(GeCloverNbr); }
  // inline void SetGeCrystalNbr(const UInt_t &GeCrystalNbr){fASGARD_CrystalNbr.push_back(GeCrystalNbr);}
  // inline void SetGeEnergy(const Double_t &GeEnergy){fASGARD_Energy.push_back(GeEnergy);} // new energy
  //  inline void SetGeSegmentNbr(const UInt_t &GeSegmentNbr){fASGARD_SegmentNbr.push_back(GeSegmentNbr);} // maximum-energy core
  //  inline void SetGeTimeCFD(const Double_t &GeTimeCFD){fASGARD_TimeCFD.push_back(GeTimeCFD);}
  //  inline void SetGeTimeLED(const Double_t &GeTimeLED){fASGARD_TimeLED.push_back(GeTimeLED);}
  
  void SortSegmentData(); // add up segment energies and determine crystal energy with maximum-energy segment
  void SetNewGeData(UInt_t GeCloverNbr, UInt_t GeCrystalNbr, UInt_t GeSegmentNbr, Double_t GeEnergy, Double_t GeTimeLED, Double_t GeThetaSemiTrue, Double_t GePhiSemiTrue, TVector3 GePosition);
  void AddGeEnergy(const UInt_t &index, const Double_t &GeEnergy); // total energy
//  void SetGeTimeLED(const UInt_t &index, const Double_t &GeTimeLED)
  
  // inline void SetBGOCloverNbr(const UInt_t &BGOCloverNbr){fASG_BGO_CloverNbr.push_back(BGOCloverNbr); }
  // inline void SetBGOCrystalNbr(const UInt_t &BGOCrystalNbr){fASG_BGO_CrystalNbr.push_back(BGOCrystalNbr);}
  // inline void SetBGOPmNbr(const UInt_t &BGOPmNbr){fASG_BGO_PmNbr.push_back(BGOPmNbr);}
  // inline void SetBGOEnergy(const Double_t &BGOEnergy){fASG_BGO_Energy.push_back(BGOEnergy);}
  // inline void SetBGOTimeCFD(const Double_t &BGOTimeCFD){fASG_BGO_TimeCFD.push_back(BGOTimeCFD);}
  // inline void SetBGOTimeLED(const Double_t &BGOTimeLED){fASG_BGO_TimeLED.push_back(BGOTimeLED);}

  /////////////////////           GETTERS           ////////////////////////
  inline UInt_t GetGeCloverNbr(const unsigned int &i)   {return fASGARD_CloverNbr[i]; }
  inline UInt_t GetGeCrystalNbr(const unsigned int &i)  {return fASGARD_CrystalNbr[i]; }
  inline UInt_t GetGeSegmentNbr(const unsigned int &i)  {return fASGARD_SegmentNbr[i]; }
  
  inline Double_t GetGeEnergy(const unsigned int &i)      {return fASGARD_Energy[i];}
  //  inline Double_t GetGeTimeCFD(const unsigned int &i)     {return fASGARD_TimeCFD[i];}
  inline Double_t GetGeTimeLED(const unsigned int &i)     {return fASGARD_TimeLED[i];}
  inline Double_t GetGeThetaSemiTrue(const unsigned int &i)     {return fASGARD_Theta_SemiTrue[i];}
  inline Double_t GetGePhiSemiTrue(const unsigned int &i)     {return fASGARD_Phi_SemiTrue[i];}
  inline TVector3 GetGePosition(const unsigned int &i) {return fASGARD_position[i];}

  // inline UInt_t GetBGOCloverNbr(const unsigned int &i)   {return fASG_BGO_CloverNbr[i]; }
  // inline UInt_t GetBGOCrystalNbr(const unsigned int &i)  {return fASG_BGO_CrystalNbr[i]; }
  // inline UInt_t GetBGOPmNbr(const unsigned int &i)       {return fASG_BGO_PmNbr[i]; }
  // inline Double_t GetBGOEnergy(const unsigned int &i)      {return fASG_BGO_Energy[i];}
  // inline Double_t GetBGOTimeCFD(const unsigned int &i)     {return fASG_BGO_TimeCFD[i];}
  // inline Double_t GetBGOTimeLED(const unsigned int &i)     {return fASG_BGO_TimeLED[i];}

  inline unsigned int GetMultiplicityGe()  {return fASGARD_Energy.size();}
  //  inline unsigned int GetMultiplicityBGO()  {return fASG_BGO_CloverNbr.size();}
  
  ClassDef(TASGARDData,1)  // ASGARDData structure
};

#endif

















