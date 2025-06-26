#ifndef __ZDDDATA__
#define __ZDDDATA__
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Raw data                                    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <iostream>
#include <vector>
using namespace std;

// ROOT
#include "TObject.h"

class TZDDData : public TObject {
    //////////////////////////////////////////////////////////////
    // data members are hold into vectors in order
    // to allow multiplicity treatment
private:
    // IC
    vector<unsigned int> fZDD_IC_E;
    vector<unsigned int> fZDD_IC_N;
    vector<unsigned long long> fZDD_IC_TS;

    // Plastic
    vector<unsigned int> fZDD_PL_E;
    vector<unsigned int> fZDD_PL_N;
    vector<unsigned long long> fZDD_PL_TS;

    // DC
    vector<unsigned int> fZDD_DC_E;
    vector<unsigned int> fZDD_DC_N;
    vector<unsigned long long> fZDD_DC_TS;

    // DriftTime in DC
    vector<UShort_t> fZDD_Drift_DetectorNbr;
    vector<Double_t> fZDD_DriftTime;
    vector<Double_t> fZDD_X;

    // EXOZDD
    vector<unsigned int> fZDD_EXO_E;
    vector<unsigned int> fZDD_EXO_N;
    vector<unsigned long long> fZDD_EXO_TS;

    //////////////////////////////////////////////////////////////
    // Constructor and destructor
public:
    TZDDData();
    ~TZDDData();

    //////////////////////////////////////////////////////////////
    // Inherited from TObject and overriden to avoid warnings
public:
    void Clear();
    void Clear(const Option_t*){};
    void Dump() const;

    //////////////////////////////////////////////////////////////
    // Getters and Setters
    // Prefer inline declaration to avoid unnecessary called of
    // frequently used methods
    // add //! to avoid ROOT creating dictionnary for the methods
public:
    //////////////////////    SETTERS    ////////////////////////
    // Energy
    inline void SetZDDIC(const unsigned int& DetNumb, const unsigned int& Energy, const unsigned long long& TimeStamp) {
        fZDD_IC_E.push_back(Energy);
        fZDD_IC_N.push_back(DetNumb);
        fZDD_IC_TS.push_back(TimeStamp);
    }; //!
    inline void SetZDDPL(const unsigned int& DetNumb, const unsigned int& Energy, const unsigned long long& TimeStamp) {
        fZDD_PL_E.push_back(Energy);
        fZDD_PL_N.push_back(DetNumb);
        fZDD_PL_TS.push_back(TimeStamp);
    }; //!
    inline void SetZDDDC(const unsigned int& DetNumb, const unsigned int& Energy, const unsigned long long& TimeStamp) {
        fZDD_DC_E.push_back(Energy);
        fZDD_DC_N.push_back(DetNumb);
        fZDD_DC_TS.push_back(TimeStamp);
    }; //!
    inline void SetDrift(const UShort_t& DetNbr, const Double_t& DriftTime, const Double_t& X) {
        fZDD_Drift_DetectorNbr.push_back(DetNbr);
        fZDD_DriftTime.push_back(DriftTime);
        fZDD_X.push_back(X);
    }; //!
    inline void SetZDDEXO(const unsigned int& DetNumb, const unsigned int& Energy, const unsigned long long& TimeStamp) {
        fZDD_EXO_E.push_back(Energy);
        fZDD_EXO_N.push_back(DetNumb);
        fZDD_EXO_TS.push_back(TimeStamp);
    }; //!

    //////////////////////    GETTERS    ////////////////////////
    inline unsigned int	GetZDD_ICMult()		                const {return fZDD_IC_E.size();}
    inline unsigned int GetZDD_ICE(const unsigned int& i) const { return fZDD_IC_E[i]; }
    inline unsigned int GetZDD_ICN(const unsigned int& i) const { return fZDD_IC_N[i]; }
    inline unsigned long long GetZDD_ICTS(const unsigned int& i) const { return fZDD_IC_TS[i]; }

    inline unsigned int	GetZDD_PLMult()		                const {return fZDD_PL_E.size();}
    inline unsigned int GetZDD_PLE(const unsigned int& i) const { return fZDD_PL_E[i]; }
    inline unsigned int GetZDD_PLN(const unsigned int& i) const { return fZDD_PL_N[i]; }
    inline unsigned long long GetZDD_PLTS(const unsigned int& i) const { return fZDD_PL_TS[i]; }

    inline unsigned int	GetZDD_DCMult()		                const {return fZDD_DC_E.size();}
    inline unsigned int GetZDD_DCE(const unsigned int& i) const { return fZDD_DC_E[i]; }
    inline unsigned int GetZDD_DCN(const unsigned int& i) const { return fZDD_DC_N[i]; }
    inline unsigned long long GetZDD_DCTS(const unsigned int& i) const { return fZDD_DC_TS[i]; }

    inline UShort_t GetMultDrift() const { return fZDD_Drift_DetectorNbr.size(); }
    inline UShort_t GetDrift_DetectorNbr(const unsigned int& i) const { return fZDD_Drift_DetectorNbr[i];}
    inline Double_t Get_DriftTime(const unsigned int& i) const { return fZDD_DriftTime[i]; }

    inline unsigned int	GetZDD_EXOMult()		                const {return fZDD_EXO_E.size();}
    inline unsigned int GetZDD_EXOE(const unsigned int& i) const { return fZDD_EXO_E[i]; }
    inline unsigned int GetZDD_EXON(const unsigned int& i) const { return fZDD_EXO_N[i]; }
    inline unsigned long long GetZDD_EXOTS(const unsigned int& i) const { return fZDD_EXO_TS[i]; }

    //////////////////////////////////////////////////////////////
    // Required for ROOT dictionnary
    ClassDef(TZDDData, 1) // ZDDData structure
};

#endif
