#ifndef __EpicDATA__
#define __EpicDATA__
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon  contact address: audrey.chatillon@cea.fr                        *
 *                                                                           *
 * Creation Date  : décembre 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Epic Raw data                                    *
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

class TEpicData : public TObject {

 public:
    // per TrackID

    struct EpicAnodeData {
      UShort_t    AnodeNbr;
      Double_t    Q1;
      Double_t    Qmax;
      Double_t    Xpos;
      Double_t    Ypos;
      Double_t    Zpos;

      // Getters
      UShort_t GetAnodeNbr() const { return AnodeNbr; }
      Double_t GetQ1() const { return Q1; }
      Double_t GetQmax() const { return Qmax; }
      Double_t GetXpos() const { return Xpos; }
      Double_t GetYpos() const { return Ypos; }
      Double_t GetZpos() const { return Zpos; }
   };

  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private: 
    vector<EpicAnodeData> fEpic_Data;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TEpicData();
    ~TEpicData();
    

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
    void Set(const UShort_t& n, const Double_t& q1, const Double_t& qmax, 
             const Double_t& x, const Double_t& y,  const Double_t& z) {
        fEpic_Data.push_back({n, q1, qmax, x, y, z});
    }
    const EpicAnodeData& operator[](const unsigned int& i) const {return fEpic_Data[i];}

    ////////////////////////    GETTERS    ////////////////////////
    inline UShort_t GetMultiplicity() const               {return fEpic_Data.size();}             // num of particles in the sensitive volume
    UShort_t GetAnodeNbr(const unsigned short& i) const   { return fEpic_Data[i].AnodeNbr; }      // sensitive volume
    Double_t GetQ1(const unsigned int& i) const           { return fEpic_Data[i].Q1; }            // Q1 on all the sensitive volume
    Double_t GetQmax(const unsigned int& i) const           { return fEpic_Data[i].Qmax; }            // Qmax on all the sensitive volume
    Double_t GetXpos(const unsigned int& i) const           { return fEpic_Data[i].Xpos; }            // x position at the first step in the sensitive volume
    Double_t GetYpos(const unsigned int& i) const           { return fEpic_Data[i].Ypos; }            // y position at the first step in the sensitive volume
    Double_t GetZpos(const unsigned int& i) const           { return fEpic_Data[i].Zpos; }            // z position at the first step in the  sensitive volume
    
  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TEpicData,1)  // EpicData structure
};

#endif
