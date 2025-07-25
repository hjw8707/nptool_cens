/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the ASGARD  raw data (Made for ASG10 card)              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <fstream>
#include <iostream>
#include <sstream>
using namespace std;

#include "TASGARDData.h"

ClassImp(TASGARDData)

    /////////////////////////
    TASGARDData::TASGARDData() {}

/////////////////////////
TASGARDData::~TASGARDData() {}

/////////////////////////
void TASGARDData::Clear() {
    fASGARD_CloverNbr.clear();
    fASGARD_CrystalNbr.clear();
    fASGARD_SegmentNbr.clear();
    fASGARD_Energy.clear();
    fASGARD_MaxEnergySegment.clear();
    fASGARD_TimeLED.clear();
    fASGARD_Theta_SemiTrue.clear();
    fASGARD_Phi_SemiTrue.clear();
    fASGARD_position.clear();

    fASGARD_CloverNbr_sub.clear();
    fASGARD_CrystalNbr_sub.clear();
    fASGARD_SegmentNbr_sub.clear();
    fASGARD_Energy_sub.clear();
    fASGARD_TimeLED_sub.clear();
    fASGARD_Theta_SemiTrue_sub.clear();
    fASGARD_Phi_SemiTrue_sub.clear();
    fASGARD_position_sub.clear();

    // fASG_BGO_CloverNbr.clear();
    // fASG_BGO_CrystalNbr.clear();
    // fASG_BGO_PmNbr.clear();
    // fASG_BGO_Energy.clear();
    // fASG_BGO_TimeCFD.clear();
    // fASG_BGO_TimeLED.clear();
}

/////////////////////////
void TASGARDData::SetNewGeData(UInt_t GeCloverNbr, UInt_t GeCrystalNbr, UInt_t GeSegmentNbr, Double_t GeEnergy,
                               Double_t GeTimeLED, Double_t GeThetaSemiTrue, Double_t GePhiSemiTrue,
                               TVector3 GePosition) {
    fASGARD_CloverNbr.push_back(GeCloverNbr);
    fASGARD_CrystalNbr.push_back(GeCrystalNbr);
    fASGARD_SegmentNbr.push_back(GeSegmentNbr);
    fASGARD_Energy.push_back(GeEnergy);
    fASGARD_MaxEnergySegment.push_back(GeEnergy);
    fASGARD_TimeLED.push_back(GeTimeLED);  // take time of maximum energy deposited (1st interaction assumed)
    fASGARD_Theta_SemiTrue.push_back(GeThetaSemiTrue);
    fASGARD_Phi_SemiTrue.push_back(GePhiSemiTrue);
    fASGARD_position.push_back(GePosition);
}
void TASGARDData::AddGeEnergy(const UInt_t &index, const Double_t &GeEnergy) {
    //  cout<<"index: "<<index<<endl;
    //  if (index > -1)
    fASGARD_Energy[index] += GeEnergy;
    // else // starting fresh
    //   fASGARD_Energy.push_back(GeEnergy);
}  // total energy

void TASGARDData::SortSegmentData() {                               // loop through
    for (unsigned int i = 0; i < fASGARD_Energy_sub.size(); i++) {  //
        // check for new clover, crystal number
        if (fASGARD_CloverNbr_sub[i] != last_clover ||
            fASGARD_CrystalNbr_sub[i] != last_crystal) {  // new clover/crystal, register
            SetNewGeData(fASGARD_CloverNbr_sub[i], fASGARD_CrystalNbr_sub[i], fASGARD_SegmentNbr_sub[i],
                         fASGARD_Energy_sub[i], fASGARD_TimeLED_sub[i], fASGARD_Theta_SemiTrue_sub[i],
                         fASGARD_Phi_SemiTrue_sub[i], fASGARD_position_sub[i]);
            last_clover = fASGARD_CloverNbr_sub[i];
            last_crystal = fASGARD_CrystalNbr_sub[i];
        } else {  // same clover & crystal: add energy, update time and check for higher segment energy
            int index_now = fASGARD_Energy.size() - 1;
            AddGeEnergy(index_now, fASGARD_Energy_sub[i]);
            // check maximum energy
            if (fASGARD_Energy_sub[i] > fASGARD_MaxEnergySegment[index_now]) {
                fASGARD_MaxEnergySegment[index_now] = fASGARD_Energy_sub[i];
                fASGARD_SegmentNbr[index_now] = fASGARD_SegmentNbr_sub[i];
                fASGARD_TimeLED[index_now] =
                    fASGARD_TimeLED_sub[i];  // take time of maximum energy deposited (1st interaction assumed)
                fASGARD_Theta_SemiTrue[index_now] =
                    fASGARD_Theta_SemiTrue_sub[i];  // take time of maximum energy deposited (1st interaction assumed)
                fASGARD_Phi_SemiTrue[index_now] =
                    fASGARD_Phi_SemiTrue_sub[i];  // take time of maximum energy deposited (1st interaction assumed)
                fASGARD_position[index_now] = fASGARD_position_sub[i];
            }
        }
    }
}
void TASGARDData::Dump() const {
    // Energy
    cout << "ASGARD_Mult = " << fASGARD_CloverNbr.size() << endl;

    // Front
    for (UInt_t i = 0; i < fASGARD_CloverNbr.size(); i++) {
        cout << "Clover: " << fASGARD_CloverNbr[i] << " Crystal: " << fASGARD_CrystalNbr[i]
             << " Energy: " << fASGARD_Energy[i] << " Time: " << fASGARD_TimeLED[i] << endl;
    }
}
