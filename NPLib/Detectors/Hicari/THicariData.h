#ifndef __HICARIDATA__
#define __HICARIDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: F.Flavigny  contact: flavigny@lpccaen.in2p3.fr           *
 *                                                                           *
 * Creation Date  : april 2022                                               *
 * Last update    : april 2022                                               *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Hicari Raw data                                          *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/


// STL
#include <vector>
using namespace std;
#include "TObject.h"


class THicariData : public TObject {
 private:
  // real energy value (for npsimulation)
  vector<unsigned int>	fHi_Cluster;
  vector<unsigned int>	fHi_Crystal;
  vector<unsigned int>	fHi_Segment;
  vector<double> 	fHi_Energy;
  vector<double> 	fHi_Time;
  //vector<vector<double>> 	fHi_SegEnergy;
  //vector<vector<unsigned int>> 	fHi_Seg;

 public:
  THicariData();
  virtual ~THicariData();

  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;

  unsigned int Find(const unsigned int&, const unsigned int&, const unsigned int&);
  void AddE(const unsigned int& index, const double en){fHi_Energy[index]+=en;}
  void ChangeE(const unsigned int& index, const double en){fHi_Energy[index]=en;}
  unsigned int GetMult(){return fHi_Cluster.size();}

  /////////////////////           SETTERS           ////////////////////////
    void	SetCluster(unsigned int clu)	{ fHi_Cluster.push_back(clu);}
    void	SetCrystal(unsigned int crys)	{ fHi_Crystal.push_back(crys);}
    void	SetSegment(unsigned int seg)	{ fHi_Segment.push_back(seg);}
    void	SetEnergy(double ener)	{ fHi_Energy.push_back(ener);}
    void	SetTime(double time)	{ fHi_Time.push_back(time);}

  /////////////////////           GETTERS           ////////////////////////
      unsigned int	GetCluster(Int_t i)	{return fHi_Cluster[i];}
      unsigned int	GetCrystal(Int_t i)	{return fHi_Crystal[i];}
      unsigned int	GetSegment(Int_t i)	{return fHi_Segment[i];}
      unsigned int	GetEnergy(Int_t i)	{return fHi_Energy[i];}
      unsigned int	GetTime(Int_t i)	{return fHi_Time[i];}

      ClassDef(THicariData,1)  // HicariData structure
	};

#endif
