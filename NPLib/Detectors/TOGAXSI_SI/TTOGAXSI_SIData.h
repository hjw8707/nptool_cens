#ifndef __TOGAXSI_SIDATA__
#define __TOGAXSI_SIDATA__
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
 *  This class hold TOGAXSI_SI Raw data                                    *
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

class TTOGAXSI_SIData : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment

  public:
    // InnerX Energy
    vector<unsigned short> fInnerX_E_DetectorNbr;
    vector<unsigned short> fInnerX_E_StripNbr;
    vector<double>	   fInnerX_E_Energy;
    // InnerZ Energy
    vector<unsigned short> fInnerZ_E_DetectorNbr;
    vector<unsigned short> fInnerZ_E_StripNbr;
    vector<double>	   fInnerZ_E_Energy;

    // OuterX Energy
    vector<unsigned short> fOuterX_E_DetectorNbr;
    vector<unsigned short> fOuterX_E_StripNbr;
    vector<double>	   fOuterX_E_Energy;
    // OuterZ Energy
    vector<unsigned short> fOuterZ_E_DetectorNbr;
    vector<unsigned short> fOuterZ_E_StripNbr;
    vector<double>	   fOuterZ_E_Energy;

    // ClusterInner Energy
    vector<unsigned short> fClusterInner_E_DetectorNbr;
    vector<unsigned short> fClusterInner_E_StripNbr;
    vector<double>	   fClusterInner_E_Energy;

    // ClusterX Energy
    vector<unsigned short> fClusterX1_E_DetectorNbr;
    vector<unsigned short> fClusterX1_E_StripNbr;
    vector<double>	   fClusterX1_E_Energy;

    // ClusteY Energy
    vector<unsigned short> fClusterY1_E_DetectorNbr;
    vector<unsigned short> fClusterY1_E_StripNbr;
    vector<double>	   fClusterY1_E_Energy;

    // ClusterX Energy
    vector<unsigned short> fClusterX2_E_DetectorNbr;
    vector<unsigned short> fClusterX2_E_StripNbr;
    vector<double>	   fClusterX2_E_Energy;

    // ClusteY Energy
    vector<unsigned short> fClusterY2_E_DetectorNbr;
    vector<unsigned short> fClusterY2_E_StripNbr;
    vector<double>	   fClusterY2_E_Energy;


  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TTOGAXSI_SIData();
    ~TTOGAXSI_SIData();
    

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
    // InnerX Energy
    inline void SetInnerXE(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fInnerX_E_DetectorNbr.push_back(DetNbr);
      fInnerX_E_StripNbr.push_back(StripNbr);
      fInnerX_E_Energy.push_back(Energy);
    };//!

    // InnerZ Energy
    inline void SetInnerZE(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fInnerZ_E_DetectorNbr.push_back(DetNbr);
      fInnerZ_E_StripNbr.push_back(StripNbr);
      fInnerZ_E_Energy.push_back(Energy);
    };//!

    // OuterX Energy
    inline void SetOuterXE(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fOuterX_E_DetectorNbr.push_back(DetNbr);
      fOuterX_E_StripNbr.push_back(StripNbr);
      fOuterX_E_Energy.push_back(Energy);
    };//!

    // OuterZ Energy
    inline void SetOuterZE(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fOuterZ_E_DetectorNbr.push_back(DetNbr);
      fOuterZ_E_StripNbr.push_back(StripNbr);
      fOuterZ_E_Energy.push_back(Energy);
    };//!

    // ClusterInner Energy Front
    inline void SetClusterInnerE(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fClusterInner_E_DetectorNbr.push_back(DetNbr);
      fClusterInner_E_StripNbr.push_back(StripNbr);
      fClusterInner_E_Energy.push_back(Energy);
    };//!

    // ClusterX Energy Front
    inline void SetClusterX1E(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fClusterX1_E_DetectorNbr.push_back(DetNbr);
      fClusterX1_E_StripNbr.push_back(StripNbr);
      fClusterX1_E_Energy.push_back(Energy);
    };//!

    // ClusterY Energy Front
    inline void SetClusterY1E(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fClusterY1_E_DetectorNbr.push_back(DetNbr);
      fClusterY1_E_StripNbr.push_back(StripNbr);
      fClusterY1_E_Energy.push_back(Energy);
    };//!

    // ClusterX Energy Front
    inline void SetClusterX2E(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fClusterX2_E_DetectorNbr.push_back(DetNbr);
      fClusterX2_E_StripNbr.push_back(StripNbr);
      fClusterX2_E_Energy.push_back(Energy);
    };//!

    // ClusterY Energy Front
    inline void SetClusterY2E(const unsigned& DetNbr, const unsigned short& StripNbr, const Double_t& Energy){
      fClusterY2_E_DetectorNbr.push_back(DetNbr);
      fClusterY2_E_StripNbr.push_back(StripNbr);
      fClusterY2_E_Energy.push_back(Energy);
    };//!


    //////////////////////    GETTERS    ////////////////////////
 // InnerX Energy
  inline unsigned short GetInnerXMultEnergy() const {return fInnerX_E_DetectorNbr.size();}
  inline unsigned short GetInnerX_E_DetectorNbr(const unsigned int &i) const {return fInnerX_E_DetectorNbr[i];}//!
  inline unsigned short GetInnerX_E_StripNbr(const unsigned int &i) const {return fInnerX_E_StripNbr[i];}//!
  inline Double_t GetInnerX_E_Energy(const unsigned int &i) const {return fInnerX_E_Energy[i];}//!

  // InnerZ Energy
  inline unsigned short GetInnerZMultEnergy() const {return fInnerZ_E_DetectorNbr.size();}
  inline unsigned short GetInnerZ_E_DetectorNbr(const unsigned int &i) const {return fInnerZ_E_DetectorNbr[i];}//!
  inline unsigned short GetInnerZ_E_StripNbr(const unsigned int &i) const {return fInnerZ_E_StripNbr[i];}//!
  inline Double_t GetInnerZ_E_Energy(const unsigned int &i) const {return fInnerZ_E_Energy[i];}//!

  // OuterX Energy
  inline unsigned short GetOuterXMultEnergy() const {return fOuterX_E_DetectorNbr.size();}
  inline unsigned short GetOuterX_E_DetectorNbr(const unsigned int &i) const {return fOuterX_E_DetectorNbr[i];}//!
  inline unsigned short GetOuterX_E_StripNbr(const unsigned int &i) const {return fOuterX_E_StripNbr[i];}//!
  inline Double_t GetOuterX_E_Energy(const unsigned int &i) const {return fOuterX_E_Energy[i];}//!

  // OuterZ Energy
  inline unsigned short GetOuterZMultEnergy() const {return fOuterZ_E_DetectorNbr.size();}
  inline unsigned short GetOuterZ_E_DetectorNbr(const unsigned int &i) const {return fOuterZ_E_DetectorNbr[i];}//!
  inline unsigned short GetOuterZ_E_StripNbr(const unsigned int &i) const {return fOuterZ_E_StripNbr[i];}//!
  inline Double_t GetOuterZ_E_Energy(const unsigned int &i) const {return fOuterZ_E_Energy[i];}//!

  // ClusterInner Energy
  inline unsigned short GetClusterInnerMultEnergy() const {return fClusterInner_E_DetectorNbr.size();}
  inline unsigned short GetClusterInner_E_DetectorNbr(const unsigned int &i) const {return fClusterInner_E_DetectorNbr[i];}//!
  inline unsigned short GetClusterInner_E_StripNbr(const unsigned int &i) const {return fClusterInner_E_StripNbr[i];}//!
  inline Double_t GetClusterInner_E_Energy(const unsigned int &i) const {return fClusterInner_E_Energy[i];}//!

  // ClusterX1 Energy
  inline unsigned short GetClusterX1MultEnergy() const {return fClusterX1_E_DetectorNbr.size();}
  inline unsigned short GetClusterX1_E_DetectorNbr(const unsigned int &i) const {return fClusterX1_E_DetectorNbr[i];}//!
  inline unsigned short GetClusterX1_E_StripNbr(const unsigned int &i) const {return fClusterX1_E_StripNbr[i];}//!
  inline Double_t GetClusterX1_E_Energy(const unsigned int &i) const {return fClusterX1_E_Energy[i];}//!

  // ClusterY1 Energy
  inline unsigned short GetClusterY1MultEnergy() const {return fClusterY1_E_DetectorNbr.size();}
  inline unsigned short GetClusterY1_E_DetectorNbr(const unsigned int &i) const {return fClusterY1_E_DetectorNbr[i];}//!
  inline unsigned short GetClusterY1_E_StripNbr(const unsigned int &i) const {return fClusterY1_E_StripNbr[i];}//!
  inline Double_t GetClusterY1_E_Energy(const unsigned int &i) const {return fClusterY1_E_Energy[i];}//!

  // ClusterX2 Energy
  inline unsigned short GetClusterX2MultEnergy() const {return fClusterX2_E_DetectorNbr.size();}
  inline unsigned short GetClusterX2_E_DetectorNbr(const unsigned int &i) const {return fClusterX2_E_DetectorNbr[i];}//!
  inline unsigned short GetClusterX2_E_StripNbr(const unsigned int &i) const {return fClusterX2_E_StripNbr[i];}//!
  inline Double_t GetClusterX2_E_Energy(const unsigned int &i) const {return fClusterX2_E_Energy[i];}//!

  // ClusterY2 Energy
  inline unsigned short GetClusterY2MultEnergy() const {return fClusterY2_E_DetectorNbr.size();}
  inline unsigned short GetClusterY2_E_DetectorNbr(const unsigned int &i) const {return fClusterY2_E_DetectorNbr[i];}//!
  inline unsigned short GetClusterY2_E_StripNbr(const unsigned int &i) const {return fClusterY2_E_StripNbr[i];}//!
  inline Double_t GetClusterY2_E_Energy(const unsigned int &i) const {return fClusterY2_E_Energy[i];}//!




  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TTOGAXSI_SIData,1)  // TOGAXSI_SIData structure
};

#endif
