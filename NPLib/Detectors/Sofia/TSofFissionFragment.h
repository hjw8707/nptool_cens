#ifndef __FissionFragmentDATA__
#define __FissionFragmentDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2021                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold FissionFragment Raw data                                    *
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

class TSofFissionFragment : public TObject {
  //////////////////////////////////////////////////////////////
  // data members are hold into vectors in order 
  // to allow multiplicity treatment
  private:
    vector<double> fFF_Z; 
    vector<int> fFF_iZ;
    vector<int> fFF_Plastic;
    vector<double> fFF_AoQ;
    vector<double> fFF_A;
    vector<double> fFF_Beta;
    vector<double> fFF_TOF;
    vector<double> fFF_Gamma;
    vector<double> fFF_Brho;
    vector<double> fFF_Brho_X;
    vector<double> fFF_Brho_Z;
    vector<double> fFF_Rho;
    vector<double> fFF_Omega;
    vector<double> fFF_DT;
    vector<int>    fFF_Section;
    vector<double> fFF_ThetaIn;
    vector<double> fFF_ThetaOut;
    vector<double> fFF_TofPosX;
    vector<double> fFF_TofPosY;
    vector<double> fFF_XB;
    vector<double> fFF_XC;
    vector<double> fFF_XD;
    vector<double> fFF_ZB;
    vector<double> fFF_ZC;
    vector<double> fFF_ZD;
    vector<double> fFF_X1;
    vector<double> fFF_Y1;
    vector<double> fFF_X2;
    vector<double> fFF_Y2;
    vector<double> fFF_X3lab;
    vector<double> fFF_Z3lab;
    vector<double> fFF_FlightPath;
    vector<double> fFF_Leff;
    vector<double> fFF_deff1;
    vector<double> fFF_deff2;
    double fFF_Zsum;
    int fFF_iZsum;
    int fFF_Cathode;

  //////////////////////////////////////////////////////////////
  // Constructor and destructor
  public: 
    TSofFissionFragment();
    ~TSofFissionFragment();
    

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
    inline void SetCathode(int val){fFF_Cathode = val;};//!
    inline void SetZsum(double val){fFF_Zsum = val;};//!
    inline void SetiZsum(int val){fFF_iZsum = val;};//!
    inline void SetPlastic(int val){fFF_Plastic.push_back(val);};//!
    inline void SetZ(double val){fFF_Z.push_back(val);};//!
    inline void SetiZ(int val){fFF_iZ.push_back(val);};//!
    inline void SetAoQ(double val){fFF_AoQ.push_back(val);};//!
    inline void SetA(double val){fFF_A.push_back(val);};//!
    inline void SetBeta(double val){fFF_Beta.push_back(val);};//!
    inline void SetTOF(double val){fFF_TOF.push_back(val);};//!
    inline void SetGamma(double val){fFF_Gamma.push_back(val);};//!
    inline void SetBrho(double val){fFF_Brho.push_back(val);};//!
    inline void SetBrhoX(double val){fFF_Brho_X.push_back(val);};//!
    inline void SetBrhoZ(double val){fFF_Brho_Z.push_back(val);};//!
    inline void SetRho(double val){fFF_Rho.push_back(val);};//!
    inline void SetOmega(double val){fFF_Omega.push_back(val);};//!
    inline void SetDT(double val){fFF_DT.push_back(val);};//!
    inline void SetSection(int val){fFF_Section.push_back(val);};//!
    inline void SetThetaIn(double val){fFF_ThetaIn.push_back(val);};//!
    inline void SetThetaOut(double val){fFF_ThetaOut.push_back(val);};//!
    inline void SetTofPosX(double val){fFF_TofPosX.push_back(val);};//!
    inline void SetTofPosY(double val){fFF_TofPosY.push_back(val);};//!
    inline void SetFlightPath(double val){fFF_FlightPath.push_back(val);};//!
    inline void SetLeff(double val){fFF_Leff.push_back(val);};//!
    inline void Setdeff1(double val){fFF_deff1.push_back(val);};//!
    inline void Setdeff2(double val){fFF_deff2.push_back(val);};//!
    inline void SetPosXB(double val){fFF_XB.push_back(val);};//!
    inline void SetPosXC(double val){fFF_XC.push_back(val);};//!
    inline void SetPosXD(double val){fFF_XD.push_back(val);};//!
    inline void SetPosZB(double val){fFF_ZB.push_back(val);};//!
    inline void SetPosZC(double val){fFF_ZC.push_back(val);};//!
    inline void SetPosZD(double val){fFF_ZD.push_back(val);};//!
    inline void SetPosX1(double val){fFF_X1.push_back(val);};//!
    inline void SetPosY1(double val){fFF_Y1.push_back(val);};//!
    inline void SetPosX2(double val){fFF_X2.push_back(val);};//!
    inline void SetPosY2(double val){fFF_Y2.push_back(val);};//!
    inline void SetPosX3lab(double val){fFF_X3lab.push_back(val);};//!
    inline void SetPosZ3lab(double val){fFF_Z3lab.push_back(val);};//!


    //////////////////////    GETTERS    ////////////////////////
    int GetMult() {return fFF_Z.size();}//!
    int GetMultTofPos() {return fFF_TofPosY.size();}//!
    /*int GetMultMwpc1() {return fFF_PosY1.size();}//!
    int GetMultMwpc2() {return fFF_PosY2.size();}//!
    int GetMultMwpc3() {return fFF_PosY3.size();}//!*/
    inline int GetCathode() const {return fFF_Cathode;}//! 
    inline double GetZsum() const {return fFF_Zsum;}//! 
    inline int GetiZsum() const {return fFF_iZsum;}//! 
    inline int GetPlastic(int i) const {return fFF_Plastic[i];}//! 
    inline double GetZ(int i) const {return fFF_Z[i];}//! 
    inline int GetiZ(int i) const {return fFF_iZ[i];}//! 
    inline double GetAoQ(int i) const {return fFF_AoQ[i];}//! 
    inline double GetA(int i) const {return fFF_A[i];}//! 
    inline double GetBeta(int i) const {return fFF_Beta[i];}//! 
    inline double GetTOF(int i) const {return fFF_TOF[i];}//! 
    inline double GetGamma(int i) const {return fFF_Gamma[i];}//! 
    inline double GetBrho(int i) const {return fFF_Brho[i];}//! 
    inline double GetBrhoX(int i) const {return fFF_Brho_X[i];}//! 
    inline double GetBrhoZ(int i) const {return fFF_Brho_Z[i];}//! 
    inline double GetRho(int i) const {return fFF_Rho[i];}//! 
    inline double GetOmega(int i) const {return fFF_Omega[i];}//! 
    inline double GetDT(int i) const {return fFF_DT[i];}//! 
    inline double GetSection(int i) const {return fFF_Section[i];}//! 
    inline double GetThetaIn(int i) const {return fFF_ThetaIn[i];}//! 
    inline double GetThetaOut(int i) const {return fFF_ThetaOut[i];}//! 
    inline double GetTofPosX(int i) const {return fFF_TofPosX[i];}//! 
    inline double GetTofPosY(int i) const {return fFF_TofPosY[i];}//! 
    inline double GetFlightPath(int i) const {return fFF_FlightPath[i];}//! 
    inline double GetLeff(int i) const {return fFF_Leff[i];}//! 
    inline double Getdeff1(int i) const {return fFF_deff1[i];}//! 
    inline double Getdeff2(int i) const {return fFF_deff2[i];}//! 
    inline double GetPosXB(int i) const {return fFF_XB[i];}//! 
    inline double GetPosXC(int i) const {return fFF_XC[i];}//! 
    inline double GetPosXD(int i) const {return fFF_XD[i];}//! 
    inline double GetPosZB(int i) const {return fFF_ZB[i];}//! 
    inline double GetPosZC(int i) const {return fFF_ZC[i];}//! 
    inline double GetPosZD(int i) const {return fFF_ZD[i];}//! 
    inline double GetPosX1(int i) const {return fFF_X1[i];}//! 
    inline double GetPosY1(int i) const {return fFF_Y1[i];}//! 
    inline double GetPosX2(int i) const {return fFF_X2[i];}//! 
    inline double GetPosY2(int i) const {return fFF_Y2[i];}//! 
    inline double GetPosX3lab(int i) const {return fFF_X3lab[i];}//! 
    inline double GetPosZ3lab(int i) const {return fFF_Z3lab[i];}//! 


  //////////////////////////////////////////////////////////////
  // Required for ROOT dictionnary
  ClassDef(TSofFissionFragment,1)  // FissionFragment structure
};

#endif
