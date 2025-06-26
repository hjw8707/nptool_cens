#ifndef __LEPSDATA__
#define __LEPSDATA__
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Leo Plagnol  contact address: plagnol@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : March 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold LEPS Raw data                                    *
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

class TLEPSData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // Energy
    vector<UShort_t>   fLEPS_E_DetectorNbr;
    vector<Double_t>   fLEPS_Energy;

    // Time
    vector<UShort_t>   fLEPS_T_DetectorNbr;
    vector<Double_t>   fLEPS_Time;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TLEPSData();
    ~TLEPSData();
    

  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public:
    void Clear();
    void Clear(const Option_t*) {};
    void Dump() const;


  //////////////////////////////////////////////////////////////
  // Getters and Setters
  // Prefer inline declaration to avoid unnecessary called of 
  // frequently used methods
  // add //! to avoid ROOT creating dictionnary for the methods
  public:
    //////////////////////    SETTERS    ////////////////////////
    // Energy
    inline void SetEnergy(const UShort_t& DetNbr,const Double_t& Energy){
      fLEPS_E_DetectorNbr.push_back(DetNbr);
      fLEPS_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetTime(const UShort_t& DetNbr,const Double_t& Time)	{
      fLEPS_T_DetectorNbr.push_back(DetNbr);     
      fLEPS_Time.push_back(Time);
    };//!


    //////////////////////    GETTERS    ////////////////////////
    // Energy
    inline UShort_t GetMultEnergy() const
      {return fLEPS_E_DetectorNbr.size();}
    inline UShort_t GetE_DetectorNbr(const unsigned int &i) const 
      {return fLEPS_E_DetectorNbr[i];}//!
    inline Double_t Get_Energy(const unsigned int &i) const 
      {return fLEPS_Energy[i];}//!

    // Time
    inline UShort_t GetMultTime() const
      {return fLEPS_T_DetectorNbr.size();}
    inline UShort_t GetT_DetectorNbr(const unsigned int &i) const 
      {return fLEPS_T_DetectorNbr[i];}//!
    inline Double_t Get_Time(const unsigned int &i) const 
      {return fLEPS_Time[i];}//!


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TLEPSData,1)  // LEPSData structure
};

#endif
