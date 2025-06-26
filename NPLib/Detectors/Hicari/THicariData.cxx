/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F. Flavigny  contact : flavigny@lpccaen.in2p3.fr         *
 *                                                                           *
 * Creation Date  : april 2022                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold HiCARI Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <iostream>
using namespace std;

#include "THicariData.h"


ClassImp(THicariData)

THicariData::THicariData()
{
   // Default constructor
   Clear();
}



THicariData::~THicariData()
{
}

unsigned int  THicariData::Find(const unsigned int& cluindex, const unsigned int& cryindex, const unsigned int& segindex){
  for(unsigned int i = 0  ; i<fHi_Cluster.size(); i++){
    if(fHi_Cluster[i]==cluindex && fHi_Crystal[i]==cryindex && fHi_Segment[i]==segindex)
      return i;
  }
  return 99999;

}

void THicariData::Clear()
{
   fHi_Cluster.clear();
   fHi_Crystal.clear();
   fHi_Segment.clear();
   fHi_Energy.clear();
   fHi_Time.clear();
}



void THicariData::Dump() const
{
   cout << "XXXXXXXXXXXXXXXXXXXXXXXX New Event XXXXXXXXXXXXXXXXX" << endl;

   // RealValues (for simulation purposes)
   cout << "Gr_Mult = " << fHi_Cluster.size() << endl;
   for (unsigned int i = 0; i < fHi_Cluster.size(); i++) {
      cout << "CloverE: " << fHi_Cluster[i] << " CristalE: " << fHi_Crystal[i]; 
      cout << "Segment: " << fHi_Segment[i]; 
      cout << " Energy: " << fHi_Energy[i];
      cout << " Time: " <<  fHi_Time[i] << endl;
   }
}
