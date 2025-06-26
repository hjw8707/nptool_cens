#ifndef __TOGAXSI_GAGGDATA__
#define __TOGAXSI_GAGGDATA__
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_GAGG Raw data                                    *
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

class TTOGAXSI_GAGGData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  //private:
  public: 
    //Recoil detector
    // Energy
    vector<UShort_t>   fRecoil_E_DetectorNbr;
    vector<Double_t>   fRecoil_Energy;

    // Time
    vector<UShort_t>   fRecoil_T_DetectorNbr;
    vector<Double_t>   fRecoil_Time;

    //Cluster detector
    // Energy
    vector<UShort_t>   fCluster_E_DetectorNbr;
    vector<Double_t>   fCluster_Energy;

    // Time
    vector<UShort_t>   fCluster_T_DetectorNbr;
    vector<Double_t>   fCluster_Time;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TTOGAXSI_GAGGData();
    ~TTOGAXSI_GAGGData();
    

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
    // Recoil detector
    // Energy
    inline void SetRecoilEnergy(const unsigned& DetNbr,const Double_t& Energy) {
      fRecoil_E_DetectorNbr.push_back(DetNbr);
      fRecoil_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetRecoilTime(const unsigned& DetNbr,const Double_t& Time) {
      fRecoil_T_DetectorNbr.push_back(DetNbr);     
      fRecoil_Time.push_back(Time);
    };//!

    //Cluster detector
    // Energy
    inline void SetClusterEnergy(const unsigned& DetNbr,const Double_t& Energy) {
      fCluster_E_DetectorNbr.push_back(DetNbr);
      fCluster_Energy.push_back(Energy);
    };//!

    // Time
    inline void SetClusterTime(const unsigned& DetNbr,const Double_t& Time) {
      fCluster_T_DetectorNbr.push_back(DetNbr);     
      fCluster_Time.push_back(Time);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // Recoil detector
    // Energy
    inline UShort_t GetRecoilMultEnergy() const
      {return fRecoil_E_DetectorNbr.size();}
    inline UShort_t GetRecoil_E_DetectorNbr(const unsigned int &i) const 
      {return fRecoil_E_DetectorNbr[i];}//!
    inline Double_t GetRecoil_Energy(const unsigned int &i) const 
      {return fRecoil_Energy[i];}//!

    // Time
    inline UShort_t GetRecoilMultTime() const
      {return fRecoil_T_DetectorNbr.size();}
    inline UShort_t GetRecoil_T_DetectorNbr(const unsigned int &i) const 
      {return fRecoil_T_DetectorNbr[i];}//!
    inline Double_t GetRecoil_Time(const unsigned int &i) const 
      {return fRecoil_Time[i];}//!

    // Cluster detector
    // Energy
    inline UShort_t GetClusterMultEnergy() const
      {return fCluster_E_DetectorNbr.size();}
    inline UShort_t GetCluster_E_DetectorNbr(const unsigned int &i) const 
      {return fCluster_E_DetectorNbr[i];}//!
    inline Double_t GetCluster_Energy(const unsigned int &i) const 
      {return fCluster_Energy[i];}//!

    // Time
    inline UShort_t GetClusterMultTime() const
      {return fCluster_T_DetectorNbr.size();}
    inline UShort_t GetCluster_T_DetectorNbr(const unsigned int &i) const 
      {return fCluster_T_DetectorNbr[i];}//!
    inline Double_t GetCluster_Time(const unsigned int &i) const 
      {return fCluster_Time[i];}//!

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TTOGAXSI_GAGGData,1)  // TOGAXSI_GAGGData structure
};

#endif
