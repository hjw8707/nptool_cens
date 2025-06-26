#ifndef __VendetaDATA__
#define __VendetaDATA__
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Vendeta Raw data                                    *
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

class TVendetaData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    
    vector<UShort_t>   fVendeta_LG_DetectorNbr;
    vector<Double_t>   fVendeta_LG_Q1;
    vector<Double_t>   fVendeta_LG_Q2;
    vector<Double_t>   fVendeta_LG_Time;
    vector<Double_t>   fVendeta_LG_Qmax; 
    vector<bool>       fVendeta_LG_IsSat; 
    vector<bool>       fVendeta_LG_IsGamma; 

    vector<UShort_t>   fVendeta_HG_DetectorNbr;
    vector<Double_t>   fVendeta_HG_Q1;
    vector<Double_t>   fVendeta_HG_Q2;
    vector<Double_t>   fVendeta_HG_Time;
    vector<Double_t>   fVendeta_HG_Qmax; 
    vector<bool>       fVendeta_HG_IsSat; 
    vector<bool>       fVendeta_HG_IsGamma; 

    //////////////////////////////////////////////////////////////
    // Constructor and destructor
  public: 
    TVendetaData();
    ~TVendetaData();

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
    
    inline void SetLGDetectorNbr(const UShort_t& DetNbr) {fVendeta_LG_DetectorNbr.push_back(DetNbr);};//!
    inline void SetLGQ1(const Double_t& Q1) {fVendeta_LG_Q1.push_back(Q1);};//!
    inline void SetLGQ2(const Double_t& Q2) {fVendeta_LG_Q2.push_back(Q2);};//!
    inline void SetLGTime(const Double_t& Time) {fVendeta_LG_Time.push_back(Time);};//!
    inline void SetLGQmax(const Double_t& Qmax) {fVendeta_LG_Qmax.push_back(Qmax);};//
    inline void SetLGIsSat(const bool& IsSat) {fVendeta_LG_IsSat.push_back(IsSat);};//
    inline void SetLGIsGamma(const bool& IsGamma) {fVendeta_LG_IsGamma.push_back(IsGamma);};//

    inline void SetHGDetectorNbr(const UShort_t& DetNbr) {fVendeta_HG_DetectorNbr.push_back(DetNbr);};//!
    inline void SetHGQ1(const Double_t& Q1) {fVendeta_HG_Q1.push_back(Q1);};//!
    inline void SetHGQ2(const Double_t& Q2) {fVendeta_HG_Q2.push_back(Q2);};//!
    inline void SetHGTime(const Double_t& Time) {fVendeta_HG_Time.push_back(Time);};//!
    inline void SetHGQmax(const Double_t& Qmax) {fVendeta_HG_Qmax.push_back(Qmax);};//
    inline void SetHGIsSat(const bool& IsSat) {fVendeta_HG_IsSat.push_back(IsSat);};//
    inline void SetHGIsGamma(const bool& IsGamma) {fVendeta_HG_IsGamma.push_back(IsGamma);};//
  
    //////////////////////    GETTERS    ////////////////////////
    inline UShort_t GetLGMultEnergy() const {return fVendeta_LG_DetectorNbr.size();}
    inline UShort_t GetLGDetectorNbr(const unsigned int &i) const {return fVendeta_LG_DetectorNbr[i];}//!
    inline Double_t GetLGQ1(const unsigned int &i) const {return fVendeta_LG_Q1[i];}//!
    inline Double_t GetLGQ2(const unsigned int &i) const {return fVendeta_LG_Q2[i];}//!
    inline Double_t GetLGTime(const unsigned int &i) const {return fVendeta_LG_Time[i];}//!
    inline Double_t GetLGQmax(const unsigned int &i) const {return fVendeta_LG_Qmax[i];}//!
    inline bool     GetLGIsSat(const unsigned int &i) const {return fVendeta_LG_IsSat[i];}//!
    inline bool     GetLGIsGamma(const unsigned int &i) const {return fVendeta_LG_IsGamma[i];}//!

    inline UShort_t GetHGMultEnergy() const {return fVendeta_HG_DetectorNbr.size();}
    inline UShort_t GetHGDetectorNbr(const unsigned int &i) const {return fVendeta_HG_DetectorNbr[i];}//!
    inline Double_t GetHGQ1(const unsigned int &i) const {return fVendeta_HG_Q1[i];}//!
    inline Double_t GetHGQ2(const unsigned int &i) const {return fVendeta_HG_Q2[i];}//!
    inline Double_t GetHGTime(const unsigned int &i) const {return fVendeta_HG_Time[i];}//!
    inline Double_t GetHGQmax(const unsigned int &i) const {return fVendeta_HG_Qmax[i];}//!
    inline bool     GetHGIsSat(const unsigned int &i) const {return fVendeta_HG_IsSat[i];}//!
    inline bool     GetHGIsGamma(const unsigned int &i) const {return fVendeta_HG_IsGamma[i];}//!

    //////////////////////////////////////////////////////////////
    // Required for ROOT dictionnary
    ClassDef(TVendetaData,1)  // VendetaData structure
};

#endif
