#ifndef __PISTADATA__
#define __PISTADATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Raw data                                    *
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

class TPISTAData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // DE //
    vector<unsigned short> fPISTA_DE_DetectorNbr;
    vector<unsigned short> fPISTA_DE_StripNbr;
    vector<double> fPISTA_DE_StripTime;
    vector<double> fPISTA_DE_StripEnergy;
    vector<double> fPISTA_DE_BackEnergy;
    vector<double> fPISTA_DE_BackTime;
    vector<double> fPISTA_DE_BackDetector;

    // DE //
    vector<unsigned short> fPISTA_E_DetectorNbr;
    vector<unsigned short> fPISTA_E_StripNbr;
    vector<double> fPISTA_E_StripEnergy;
    vector<double> fPISTA_E_BackEnergy;
    vector<double> fPISTA_E_StripTime;
    vector<double> fPISTA_E_BackTime;
    vector<double> fPISTA_E_BackDetector;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TPISTAData();
    ~TPISTAData();
    

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
    // DE
    inline void SetPISTA_DE(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& StripE, const Double_t& BackE, const Double_t& StripT, const Double_t& BackT){
      SetPISTA_DE_DetectorNbr(DetNbr);
      SetPISTA_DE_StripNbr(StripNbr);
      SetPISTA_DE_StripEnergy(StripE);
      SetPISTA_DE_BackEnergy(BackE);
      SetPISTA_DE_StripTime(StripT);
      SetPISTA_DE_BackTime(BackT);
    };//!
    inline void SetPISTA_DE_DetectorNbr(const UShort_t& DetNbr){fPISTA_DE_DetectorNbr.push_back(DetNbr);};//!
    inline void SetPISTA_DE_StripNbr(const UShort_t& StripNbr){fPISTA_DE_StripNbr.push_back(StripNbr);};//!
    inline void SetPISTA_DE_StripEnergy(const Double_t& Energy){fPISTA_DE_StripEnergy.push_back(Energy);};//!
    inline void SetPISTA_DE_BackEnergy(const Double_t& Energy){fPISTA_DE_BackEnergy.push_back(Energy);};//!
    inline void SetPISTA_DE_StripTime(const Double_t& Time){fPISTA_DE_StripTime.push_back(Time);};//!
    inline void SetPISTA_DE_BackTime(const Double_t& Time){fPISTA_DE_BackTime.push_back(Time);};//!
    inline void SetPISTA_DE_BackDetector(const UShort_t& DetNbr){fPISTA_DE_BackDetector.push_back(DetNbr);};//!

    //////
    // E
    inline void SetPISTA_E(const UShort_t& DetNbr, const UShort_t& StripNbr, const Double_t& StripE, const Double_t& BackE, const Double_t& StripT, const Double_t& BackT){
      SetPISTA_E_DetectorNbr(DetNbr);
      SetPISTA_E_StripNbr(StripNbr);
      SetPISTA_E_StripEnergy(StripE);
      SetPISTA_E_BackEnergy(BackE);
      SetPISTA_E_StripTime(StripT);
      SetPISTA_E_BackTime(BackT);
    };//!
    inline void SetPISTA_E_DetectorNbr(const UShort_t& DetNbr){fPISTA_E_DetectorNbr.push_back(DetNbr);};//!
    inline void SetPISTA_E_StripNbr(const UShort_t& StripNbr){fPISTA_E_StripNbr.push_back(StripNbr);};//!
    inline void SetPISTA_E_StripEnergy(const Double_t& Energy){fPISTA_E_StripEnergy.push_back(Energy);};//!
    inline void SetPISTA_E_BackEnergy(const Double_t& Energy){fPISTA_E_BackEnergy.push_back(Energy);};//!
    inline void SetPISTA_E_StripTime(const Double_t& Time){fPISTA_E_StripTime.push_back(Time);};//!
    inline void SetPISTA_E_BackTime(const Double_t& Time){fPISTA_E_BackTime.push_back(Time);};//!
    inline void SetPISTA_E_BackDetector(const UShort_t& DetNbr){fPISTA_E_BackDetector.push_back(DetNbr);};//!

    //////////////////////    GETTERS    ////////////////////////
    // DE
    inline UShort_t GetPISTADEMult() const
      {return fPISTA_DE_StripNbr.size();}
    inline UShort_t GetPISTADEBackMult() const
      {return fPISTA_DE_BackTime.size();}
    inline UShort_t GetPISTA_DE_DetectorNbr(const unsigned int &i) const 
      {return fPISTA_DE_DetectorNbr[i];}//!
    inline UShort_t GetPISTA_DE_StripNbr(const unsigned int &i) const 
      {return fPISTA_DE_StripNbr[i];}//!
    inline Double_t GetPISTA_DE_StripEnergy(const unsigned int &i) const 
      {return fPISTA_DE_StripEnergy[i];}//!      
    inline Double_t GetPISTA_DE_BackEnergy(const unsigned int &i) const 
      {return fPISTA_DE_BackEnergy[i];}//!
    inline Double_t GetPISTA_DE_StripTime(const unsigned int &i) const 
      {return fPISTA_DE_StripTime[i];}//!      
    inline Double_t GetPISTA_DE_BackTime(const unsigned int &i) const 
      {return fPISTA_DE_BackTime[i];}//! 
    inline Double_t GetPISTA_DE_BackDetector(const unsigned int &i) const 
      {return fPISTA_DE_BackDetector[i];}//!
         
    //////
    // E
    inline UShort_t GetPISTAEMult() const
      {return fPISTA_E_StripNbr.size();}
    inline UShort_t GetPISTAEBackMult() const
      {return fPISTA_E_BackTime.size();}
    inline UShort_t GetPISTA_E_DetectorNbr(const unsigned int &i) const 
      {return fPISTA_E_DetectorNbr[i];}//!
    inline UShort_t GetPISTA_E_StripNbr(const unsigned int &i) const 
      {return fPISTA_E_StripNbr[i];}//!
    inline Double_t GetPISTA_E_StripEnergy(const unsigned int &i) const 
      {return fPISTA_E_StripEnergy[i];}//!      
    inline Double_t GetPISTA_E_BackEnergy(const unsigned int &i) const 
      {return fPISTA_E_BackEnergy[i];}//!
    inline Double_t GetPISTA_E_StripTime(const unsigned int &i) const 
      {return fPISTA_E_StripTime[i];}//!      
    inline Double_t GetPISTA_E_BackTime(const unsigned int &i) const 
      {return fPISTA_E_BackTime[i];}//!
    inline Double_t GetPISTA_E_BackDetector(const unsigned int &i) const 
      {return fPISTA_E_BackDetector[i];}//!
 
  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TPISTAData,1)  // PISTAData structure
};

#endif
