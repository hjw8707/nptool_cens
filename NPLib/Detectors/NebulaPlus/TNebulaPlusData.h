#ifndef __NebulaPlusDATA__
#define __NebulaPlusDATA__
/*****************************************************************************
 * Copyright (C) 2009-2019   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: freddy flavigny  contact: flavigny@lpccaen.in2p3.fr      *
 *                                                                           *
 * Creation Date  : OCotber 2022                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold NebulaPlus Raw data                                      *
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

class TNebulaPlusData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    // UP // 
    // Charges and Time 
    vector<UShort_t>   fNebulaPlus_u_ID;
    vector<Double_t>   fNebulaPlus_u_Q;
    vector<Double_t>   fNebulaPlus_u_Q4;
    vector<Double_t>   fNebulaPlus_u_T;
    
    // DOWN // 
    // Charges and Time
    vector<UShort_t>   fNebulaPlus_d_ID;
    vector<Double_t>   fNebulaPlus_d_Q;
    vector<Double_t>   fNebulaPlus_d_Q4;
    vector<Double_t>   fNebulaPlus_d_T;
 

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TNebulaPlusData();
    ~TNebulaPlusData();
    

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
    // UP // 
    inline void SetUp(const Double_t& ID, const Double_t& Q, const Double_t& Q4, const Double_t& T){
    fNebulaPlus_u_ID.push_back(ID);
    fNebulaPlus_u_Q.push_back(Q);
    fNebulaPlus_u_Q4.push_back(Q4);
    fNebulaPlus_u_T.push_back(T);
    };//!

    // DOWN // 
    inline void SetDown(const Double_t& ID, const Double_t& Q, const Double_t& Q4, const Double_t& T){
    fNebulaPlus_d_ID.push_back(ID);
    fNebulaPlus_d_Q.push_back(Q);
    fNebulaPlus_d_Q4.push_back(Q4);
    fNebulaPlus_d_T.push_back(T);
    };//!

    //////////////////////    GETTERS    ////////////////////////
    // MULT //
    inline unsigned int GetMultUp() const
      {return fNebulaPlus_u_ID.size();};
    inline unsigned int GetMultDown() const
      {return fNebulaPlus_d_ID.size();};

    // Value // 
    // Up 
    inline unsigned short GetIDUp(const int& i) const
      {return fNebulaPlus_u_ID[i];};
    inline double GetQUp(const int& i) const
      {return fNebulaPlus_u_Q[i];};
    inline double GetQ4Up(const int& i) const
      {return fNebulaPlus_u_Q4[i];};
    inline double GetTUp(const int& i) const
      {return fNebulaPlus_u_T[i];};

    // Down 
    inline unsigned short GetIDDown(const int& i) const
      {return fNebulaPlus_d_ID[i];};
    inline double GetQDown(const int& i) const
      {return fNebulaPlus_d_Q[i];};
    inline double GetQ4Down(const int& i) const
      {return fNebulaPlus_d_Q4[i];};
    inline double GetTDown(const int& i) const
      {return fNebulaPlus_d_T[i];};

  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TNebulaPlusData,1)  // NebulaPlusData structure
};

#endif
