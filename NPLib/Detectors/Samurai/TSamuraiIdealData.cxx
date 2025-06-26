/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien MATTA  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   :                                                         *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
#include "TSamuraiIdealData.h"


ClassImp(TSamuraiIdealData)

TSamuraiIdealData::TSamuraiIdealData()
{
}



TSamuraiIdealData::~TSamuraiIdealData()
{
}



void TSamuraiIdealData::Clear()
{
   Detector_Number.clear();
   Dep_Energy.clear(); 
   Brho.clear();
   Pos_X.clear(); 
   Pos_Y.clear(); 
   Pos_Z.clear(); 
   Mom_Mag.clear(); 
   Mom_Theta.clear(); 
   Mom_Phi.clear();
   return;
}



void TSamuraiIdealData::Dump() const
{
   //This method is very useful for debugging and worth the dev.
   cout << "XXXXXXXXXXXXXXXXXXXXX TSamuraiIdealData: New Event XXXXXXXXXXXXXXXXX" << endl;

   // Brho
   unsigned int n = GetMult();
   for (unsigned short i = 0 ; i< n ; i++) {
      cout << "Detector Number " << Detector_Number[i] << " Deposited Energy: " << Dep_Energy[i]  << endl;
   }
   
}

void TSamuraiIdealData::SetData (short detector, double energy, double pos_x, double pos_y, 
   double pos_z, double mom_r, double mom_theta, double mom_phi, double brho){
   Detector_Number.push_back(detector);
   Dep_Energy.push_back(energy);
   Pos_X.push_back(pos_x);
   Pos_Y.push_back(pos_y);
   Pos_Z.push_back(pos_z);
   Mom_Mag.push_back(mom_r);
   Mom_Theta.push_back(mom_theta);
   Mom_Phi.push_back(mom_phi); 
   Brho.push_back(brho);
}
