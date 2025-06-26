#ifndef __ICDATA__
#define __ICDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : Oct 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold IC Raw data                                    *
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

class TICData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<int> fIC_Section;
    vector<double> fIC_Charge;
    vector<long> fIC_TS;



  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TICData();
    ~TICData();
    

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
    // X setters
    inline void SetIC_Charge(double Charge){fIC_Charge.push_back(Charge);};//!
    inline void SetIC_Section(int sec){fIC_Section.push_back(sec);};//!
    inline void SetIC_TS(long TS){fIC_TS.push_back(TS);};//!


    //////////////////////    GETTERS    ////////////////////////
    inline UShort_t GetICMult() const
      {return fIC_Section.size();}
    inline UShort_t GetIC_Section(const unsigned int &i) const 
      {return fIC_Section[i];}//!
    inline Double_t GetIC_Charge(const unsigned int &i) const 
      {return fIC_Charge[i];}//!      
    inline long GetIC_TS(const unsigned int &i) const 
      {return fIC_TS[i];}//!      
     
    //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TICData,1)  // ICData structure
};

#endif
