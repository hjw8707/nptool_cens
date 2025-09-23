#ifndef __ATOMXDATA__
#define __ATOMXDATA__
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: jungwoo  contact address: phyjics@gmail.com                        *
 *                                                                           *
 * Creation Date  : May 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ATOMX Raw data                                    *
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
#include "TVector3.h"

class TATOMXData : public TObject
{
    public: 
        std::vector<double> fEnergyLoss;
        std::vector<double> fTime;
        std::vector<TVector3> fPosition;
        int GetEnergyLossMultiplicity() const { return fEnergyLoss.size(); }

    public: 
        TATOMXData();
        ~TATOMXData();

    public:
        void Clear();
        void Clear(const Option_t*) {};
        void Dump(int max_points) const;
        void Dump() const;

        ClassDef(TATOMXData,1)
};

#endif
