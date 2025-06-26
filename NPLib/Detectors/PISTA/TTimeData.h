#ifndef __TIMEDATA__
#define __TIMEDATA__
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Théodore Efremov  contact address: theodore.efremov@cea.fr                        *
 *                                                                           *
 * Creation Date  : Oct 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Time Raw data                                    *
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

class TTimeData : public TObject {
    //////////////////////////////////////////////////////////////
    // data members are hold into vectors in order 
    // to allow multiplicity treatment
    private: 
        vector<long> fTS_MWPC13;
        vector<long> fTS_MWPC14;
        vector<long> fTS_MWPC23;
        vector<long> fTS_MWPC24;

        vector<float> fTime_MWPC13;
        vector<float> fTime_MWPC14;
        vector<float> fTime_MWPC23;
        vector<float> fTime_MWPC24;

        vector<double> fToff_DT13;
        vector<double> fToff_DT14;

        vector<short> fSection_MWPC3;
        vector<short> fSection_MWPC4;



        //////////////////////////////////////////////////////////////
        // Constructor and destructor
    public: 
        TTimeData();
        ~TTimeData();


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
        inline void SetTS_MWPC13(long time ){fTS_MWPC13.push_back(time);};//!
        inline void SetTS_MWPC14(long time ){fTS_MWPC14.push_back(time);};//!
        inline void SetTS_MWPC23(long time ){fTS_MWPC23.push_back(time);};//!
        inline void SetTS_MWPC24(long time ){fTS_MWPC24.push_back(time);};//!

        inline void SetTime_MWPC13(float time ){fTime_MWPC13.push_back(time);};//!
        inline void SetTime_MWPC14(float time ){fTime_MWPC14.push_back(time);};//!
        inline void SetTime_MWPC23(float time ){fTime_MWPC23.push_back(time);};//!
        inline void SetTime_MWPC24(float time ){fTime_MWPC24.push_back(time);};//!

        inline void SetToff_DT13(double offset ){fToff_DT13.push_back(offset);};//!
        inline void SetToff_DT14(double offset ){fToff_DT14.push_back(offset);};//!

        inline void SetSection_MWPC3(short section ){fSection_MWPC3.push_back(section);};//!
        inline void SetSection_MWPC4(short section ){fSection_MWPC4.push_back(section);};//!




        //////////////////////    GETTERS    ////////////////////////
        inline long GetTS_MWPC13(const unsigned int &i) const
        {return fTS_MWPC13.at(i) ;}//!
        inline long GetTS_MWPC14(const unsigned int &i) const
        {return fTS_MWPC14.at(i) ;}//!
        inline long GetTS_MWPC23(const unsigned int &i) const
        {return fTS_MWPC23.at(i) ;}//!
        inline long GetTS_MWPC24(const unsigned int &i) const
        {return fTS_MWPC24.at(i) ;}//!

        inline float GetTime_MWPC13(const unsigned int &i) const
        {return fTime_MWPC13.at(i) ;}//!
        inline float GetTime_MWPC14(const unsigned int &i) const
        {return fTime_MWPC14.at(i) ;}//!
        inline float GetTime_MWPC23(const unsigned int &i) const
        {return fTime_MWPC23.at(i) ;}//!
        inline float GetTime_MWPC24(const unsigned int &i) const
        {return fTime_MWPC24.at(i) ;}//!
        
        inline double GetToff_DT13(const unsigned int &i) const
        {return fToff_DT13.at(i) ;}//!
        inline double GetToff_DT14(const unsigned int &i) const
        {return fToff_DT14.at(i) ;}//!
          
        inline short GetSection_MWPC3(const unsigned int &i) const
        {return fSection_MWPC3.at(i) ;}//!
        inline short GetSection_MWPC4(const unsigned int &i) const
        {return fSection_MWPC4.at(i) ;}//!
          
        inline Int_t GetMWPC13Mult() const
        {return static_cast<int>(fTime_MWPC13.size());}         
        inline Int_t GetMWPC14Mult() const
        {return static_cast<int>(fTime_MWPC14.size());}        
        inline Int_t GetMWPC23Mult() const
        {return static_cast<int>(fTime_MWPC23.size());}        
        inline Int_t GetMWPC24Mult() const
        {return static_cast<int>(fTime_MWPC24.size());}

        //////////////////////////////////////////////////////////////
        // Required for ROOT dictionnary
        ClassDef(TTimeData,3)  // TimeData structure
};

#endif
