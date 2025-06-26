#ifndef __FPMWDATA__
#define __FPMWDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr*
 *                                                                           *
 * Creation Date  : Oct 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FPMW Raw data                                    *
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

class TFPMWData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // X strips
    vector<int> fFPMW_DetX;
    vector<int> fFPMW_StripX;
    vector<double> fFPMW_ChargeX;
    // Y strips
    vector<int> fFPMW_DetY;
    vector<int> fFPMW_StripY;
    vector<double> fFPMW_ChargeY;



  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TFPMWData();
    ~TFPMWData();
    

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
    inline void SetFPMW_X(int DetNbr, int Wire, double Charge){
      SetFPMW_DetX(DetNbr);
      SetFPMW_StripX(Wire);
      SetFPMW_ChargeX(Charge);
    };//!
    inline void SetFPMW_DetX(int DetNbr){fFPMW_DetX.push_back(DetNbr);};//!
    inline void SetFPMW_StripX(int Wire){fFPMW_StripX.push_back(Wire);};//!
    inline void SetFPMW_ChargeX(double Charge){fFPMW_ChargeX.push_back(Charge);};//!

    // Y setters
    inline void SetFPMW_Y(int DetNbr, int Wire, double Charge){
      SetFPMW_DetY(DetNbr);
      SetFPMW_StripY(Wire);
      SetFPMW_ChargeY(Charge);
    };//!
    inline void SetFPMW_DetY(int DetNbr){fFPMW_DetY.push_back(DetNbr);};//!
    inline void SetFPMW_StripY(int Wire){fFPMW_StripY.push_back(Wire);};//!
    inline void SetFPMW_ChargeY(double Charge){fFPMW_ChargeY.push_back(Charge);};//!


    //////////////////////    GETTERS    ////////////////////////
    // X
    inline UShort_t GetFPMWMultX() const
      {return fFPMW_DetX.size();}
    inline UShort_t GetFPMW_DetX(const unsigned int &i) const 
      {return fFPMW_DetX[i];}//!
    inline UShort_t GetFPMW_StripX(const unsigned int &i) const 
      {return fFPMW_StripX[i];}//!
    inline Double_t GetFPMW_ChargeX(const unsigned int &i) const 
      {return fFPMW_ChargeX[i];}//!      
    // Y
    inline UShort_t GetFPMWMultY() const
      {return fFPMW_DetY.size();}
    inline UShort_t GetFPMW_DetY(const unsigned int &i) const 
      {return fFPMW_DetY[i];}//!
    inline UShort_t GetFPMW_StripY(const unsigned int &i) const 
      {return fFPMW_StripY[i];}//!
    inline Double_t GetFPMW_ChargeY(const unsigned int &i) const 
      {return fFPMW_ChargeY[i];}//!      
    
    //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TFPMWData,1)  // FPMWData structure
};

#endif
