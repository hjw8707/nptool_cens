/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Cyril Lenain   contact address: lenain@lpcaen.in2p3.fr *
 *                                                                           *
 * Creation Date   : November 2020                                           *
 * Last update     : November 2020                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class find multi-neutrons in an neutron detection array by          *
 *    rejecting "crosstalk"                                                  *
 *****************************************************************************/

#include"NPCrossTalk.h"
#include "Math/Factory.h"
#include "TError.h"
#include "TPad.h"
#include "TGraph.h"

#include <iostream>
using namespace std;
using namespace NPL;

////////////////////////////////////////////////////////////////////////////////
CrossTalk::CrossTalk(){
  // Radius of a cluster = Coef * dXdYdZ
  // this avoid error 
  gErrorIgnoreLevel = kError;

}
////////////////////////////////////////////////////////////////////////////////
CrossTalk::~CrossTalk(){
}

////////////////////////////////////////////////////////////////////////////////
void CrossTalk::AddHitVector(const vector<double>& X, const vector<double>& Y, const vector<double>& Z, const vector<double>& Q, const vector<double>& dX, const vector<double>& dY, const vector<double>& dZ, const vector<double> &T, const vector<int> &ID){


  HitX = &X;
  HitY = &Y;
  HitZ = &Z;

  HitdX = &dX;
  HitdY = &dY;
  HitdZ = &dZ;

  HitT = &T;
  HitQ = &Q;
  HitID = &ID;
  sizeHit = X.size();
}

#if __cplusplus > 201103L 
bool cmp(pair<int,double> & a, pair<int,double>&b){
    return a.second < b.second;
}
#else
bool cmp(pair<int,double>a, pair<int,double>b){
  return a.second < b.second;
}
#endif

vector<int> CrossTalk::ComputeCrossTalk(const double& Causality, const double& DistMin, const int& Option, TCutG* TCut1, TCutG* TCut2){

  const double C_light = 299.792458;
  static double x1,y1,z1,dx1,dy1,dz1,t1, q1;
  static double x2,y2,z2,dx2,dy2,dz2,t2;
  static double Dist, HyperSphere, Beta1;

  //Attribute a new index based on the time of arrival from the first to the last hit
  //Using vector of pairs (ID,Time)
  static vector<pair<int, double>> pair_SortedID;
  pair_SortedID.clear();
  for(int i = 0; i < sizeHit; i++){
    pair_SortedID.emplace_back(i, (*HitT)[i]);
  }
  //Sort pair vector (ID,Time) in Time
  sort(pair_SortedID.begin(), pair_SortedID.end(), cmp);

  m_SortedID.clear();
  //Put new ID sorted in a vector
  for(int i = 0; i < sizeHit; i++){
    m_SortedID.push_back(pair_SortedID[i].first);
  }

  // A different Cluster number (starting at 1) is assigned to each hit
  static vector<int> ID_ClustHit;
  ID_ClustHit.clear();
  for(int i =0 ; i < sizeHit; i++){
    ID_ClustHit.push_back(i+1);
  }

  //Test each Dist(n-n) to find clusters
  //When 2 neutrons are part of a the same cluster, change their ID_ClustHit for the lowest
  //Create a map to hold the Hit_ID for each cluster 
  mapOfClust.clear();
  for(int j = 0; j < sizeHit-1; j++){
    x1 = (*HitX)[m_SortedID[j]], y1 = (*HitY)[m_SortedID[j]], z1 = (*HitZ)[m_SortedID[j]], dx1 = (*HitdX)[m_SortedID[j]], dy1 = (*HitdY)[m_SortedID[j]], dz1 = (*HitdZ)[m_SortedID[j]], t1 = (*HitT)[m_SortedID[j]];   
    Beta1 = sqrt(x1*x1 + y1*y1 + z1*z1)/t1/C_light;
    for(int jj = j+1; jj < sizeHit ; jj++){
      x2 = (*HitX)[m_SortedID[jj]], y2 = (*HitY)[m_SortedID[jj]], z2 = (*HitZ)[m_SortedID[jj]], t2 = (*HitT)[m_SortedID[jj]], dx2 = (*HitdX)[m_SortedID[jj]];   
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      HyperSphere = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) + ((t2-t1)*Beta1*C_light)*((t2-t1)*Beta1*C_light) );
      if( (Dist < DistMin && Option == 1) || ( HyperSphere < coef && Option == 2) || ( (dx1+dx2 ==100. && TCut1->IsInside(t2-t1,Dist) && Option == 3  ) || ( (dx1+dx2 > 120) && TCut2->IsInside(t2-t1,Dist) && Option == 3) ) ){
        if(ID_ClustHit[jj] < ID_ClustHit[j]){
          ID_ClustHit[j] = ID_ClustHit[jj]; 
        }
        else if(ID_ClustHit[jj] > ID_ClustHit[j]){
          ID_ClustHit[jj] = ID_ClustHit[j];
        }
      }
    }
  }

  for(unsigned int i = 0; i < sizeHit; i++){
    mapOfClust[ID_ClustHit[i]].push_back(m_SortedID[i]);
  }

  //Put first hit in time of each cluster in a new map mapOfHead and remake 
  //the numbering of clusters (1,2,3...)
  static map< unsigned int, unsigned int > mapOfHead;
  static map< unsigned int, double > mapOfQtot;
  static unsigned int NbrOfClust;
  mapOfHead.clear();
  mapOfQtot.clear();
  NbrOfClust = mapOfClust.size();

  static unsigned int count, count2;
  count=1,count2=9;  m_HeadClust.clear(), m_QtotClust.clear();
  
  for(auto itr = mapOfClust.begin(); itr != mapOfClust.end(); itr++){
    mapOfHead[count] = itr->second[0];
    m_HeadClust.push_back(itr->second[0]);
    Clusters_3D->SetMarkerColor(count);
    Clusters_3D->SetLineColor(count);
    static double Qtot;
    Qtot = 0;
    for(int j =0; j < itr->second.size();j++){
      Qtot += (*HitQ)[itr->second[j]];
      /* cout << (*HitX)[itr->second[j]] << " " <<  (*HitdX)[itr->second[j]] << " " << (*HitdY)[itr->second[j]]<< ""<< (*HitdZ)[itr->second[j]]<<  endl; */
      Clusters_3D->SetPoint(count2, (*HitZ)[itr->second[j]], (*HitX)[itr->second[j]], (*HitY)[itr->second[j]]);
      Clusters_3D->SetPointError(count2, (*HitdZ)[itr->second[j]], (*HitdX)[itr->second[j]], (*HitdY)[itr->second[j]]);
      count2++;
    }
    m_QtotClust.push_back(Qtot);
    mapOfQtot[count] = Qtot;
    count++;   
  } 
  
  double v_n, Dmax, Qclust1, beta1, beta12;
  static vector<double> CrossTalk;
  CrossTalk.clear();
  //Test now the "Friend" clusters
  for(int i = 1; i < NbrOfClust; i++){ 
    x1 = (*HitX)[mapOfHead[i]], y1 = (*HitY)[mapOfHead[i]], z1 = (*HitZ)[mapOfHead[i]], t1 = (*HitT)[mapOfHead[i]]; 
    q1 = (*HitQ)[mapOfHead[i]];
    Qclust1 = mapOfQtot[i];
    v_n = sqrt(x1*x1+y1*y1+z1*z1)/t1;
    beta1 = v_n/C_light;
    for(int j = i+1; j <= NbrOfClust; j++){
      x2 = (*HitX)[mapOfHead[j]], y2 = (*HitY)[mapOfHead[j]], z2 = (*HitZ)[mapOfHead[j]], t2 = (*HitT)[mapOfHead[j]]; 
      Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
      Dmax = (t2-t1)*v_n;
      beta12 = Dist/(t2-t1)/C_light;
      if(beta1/beta12 > Causality){
        /* if( ( (*HitdX)[mapOfHead[i]] == 120. && Qclust1 < 96.9*beta1/beta12 - 70.5 )||(  (*HitdX)[mapOfHead[i]] == 50. && Qclust1 < 150*beta1/beta12 - 130 ) ){ */
        CrossTalk.push_back(j);
      }
      }
    }

    //elimate potential CrossTalk in mapOfHead
    static unsigned int sizeCT;
    sizeCT = CrossTalk.size();
    for(unsigned int i = 0; i < sizeCT; i++  ){
      mapOfHead.erase(CrossTalk[i]);
      mapOfQtot.erase(CrossTalk[i]);
    }

    //return vector of real neutrons IDs 
    m_Neutrons.clear();
    m_QtotNeut.clear();
    for(auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr++){
      m_Neutrons.push_back(itr->second);
    }
    for(auto itr = mapOfQtot.begin(); itr != mapOfQtot.end(); itr++){
      m_QtotNeut.push_back(itr->second);
    }
    return m_Neutrons;
  }

  /////////////////////////////////////////////////////
  ////////// KONDO CT /////////////////////////////////
  /////////////////////////////////////////////////////
  vector<int> CrossTalk::ComputeCrossTalk_Kondo(const double& Causality, const double& DistMin, const int& Option, TCutG* TCut1, TCutG* TCut2){

    const double C_light = 299.792458;
    static double x1,y1,z1,dx1,dy1,dz1,t1,q1,id1;
    static double x2,y2,z2,dx2,dy2,dz2,t2,q2,id2;
    static double Dt, Dist, HyperSphere, Beta1;

    //Attribute a new index based on the time arrive from the first to the last hit
    //Using vector of pairs (ID,Time)
    static vector<pair<int, double>> pair_SortedID;
    pair_SortedID.clear();
    for(int i = 0; i < sizeHit; i++){
      pair_SortedID.emplace_back(i, (*HitT)[i]);
    }
    //Sort pair vector (ID,Time) in Time
    sort(pair_SortedID.begin(), pair_SortedID.end(), cmp);

    m_SortedID.clear();
    //Put new ID sorted in a vector
    for(int i = 0; i < sizeHit; i++){
      m_SortedID.push_back(pair_SortedID[i].first);
    }

    // A different Cluster number (starting at 1) is assigned to each hit
    static vector<int> ID_ClustHit;
    ID_ClustHit.clear();
    for(int i =0 ; i < sizeHit; i++){
      ID_ClustHit.push_back(i+1);
    }

    //Test each Dist(n-n) to find clusters
    //When 2 neutrons are part of a the same cluster, change their ID_ClustHit for the lowest
    //Create a map to hold the Hit_ID for each cluster 
    mapOfClust.clear();

    ///////////////////////////////////////
    /// HyperSphere cleaning (same wall)
    //////////////////////////////////////

    for(int j = 0; j < sizeHit-1; j++){
      if(ID_ClustHit[j]==-1) continue;
      x1 = (*HitX)[m_SortedID[j]], y1 = (*HitY)[m_SortedID[j]], z1 = (*HitZ)[m_SortedID[j]];  
      dx1 = (*HitdX)[m_SortedID[j]], dy1 = (*HitdY)[m_SortedID[j]];
      dz1 = (*HitdZ)[m_SortedID[j]];
      t1 = (*HitT)[m_SortedID[j]];
      id1=(*HitID)[m_SortedID[j]];
      Beta1 = sqrt(x1*x1 + y1*y1 + z1*z1)/t1/C_light;
      for(int jj = j+1; jj < sizeHit ; jj++){
        if(ID_ClustHit[jj]==-1) continue;
        x2 = (*HitX)[m_SortedID[jj]], y2 = (*HitY)[m_SortedID[jj]], z2 = (*HitZ)[m_SortedID[jj]], t2 = (*HitT)[m_SortedID[jj]], dx2 = (*HitdX)[m_SortedID[jj]];   
        Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        Dt = t2-t1;
        id2=(*HitID)[m_SortedID[jj]];
        double b12 = Dist/Dt/C_light;
        if(id1 <= 400 && id2 <= 400){
          //Neuland same wall
          double HypeSph1 = sqrt(((Dt-0.2591)/1.165)*((Dt-0.2591)/1.165) + ((Dist-70.8)/105)*((Dist-70.8)/105));
          if(HypeSph1 < 1){
            ID_ClustHit[jj]=-1;
          }
        }
        else if( ((id1>400 && id1<=460)&&(id2>400&&id2<=460)) || (id1>460&&id2>460)  ){ 
          // Nebula same wall
          double HypeSph2 = sqrt(((Dt-0.6528)/1.608)*((Dt-0.6528)/1.608) + ((Dist-158.8)/151.5)*((Dist-158.8)/151.5));
          if(HypeSph2 < 1){
            ID_ClustHit[jj]=-1;
          }
        }
      }
    }
    
    /////////////////////////////////////////
    ///// Cleaning Gammas
    ///////////////////////////////////////////

    for(int j = 0; j < sizeHit-1; j++){
      if(ID_ClustHit[j]==-1) continue;
      x1 = (*HitX)[m_SortedID[j]], y1 = (*HitY)[m_SortedID[j]], z1 = (*HitZ)[m_SortedID[j]];  
      dx1 = (*HitdX)[m_SortedID[j]], dy1 = (*HitdY)[m_SortedID[j]];
      dz1 = (*HitdZ)[m_SortedID[j]];
      t1 = (*HitT)[m_SortedID[j]];
      id1= (*HitID)[m_SortedID[j]];
      Beta1 = sqrt(x1*x1 + y1*y1 + z1*z1)/t1/C_light;
      for(int jj = j+1; jj < sizeHit ; jj++){
        if(ID_ClustHit[jj]==-1) continue;
        x2 = (*HitX)[m_SortedID[jj]], y2 = (*HitY)[m_SortedID[jj]], z2 = (*HitZ)[m_SortedID[jj]], t2 = (*HitT)[m_SortedID[jj]], dx2 = (*HitdX)[m_SortedID[jj]];   
        Dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
        Dt = t2-t1;
        id2=(*HitID)[m_SortedID[jj]];
        double b12 = Dist/Dt/C_light;
        double q2 = (*HitQ)[m_SortedID[jj]];
        if(id1 <= 400 && id2 <= 400){
          //Neuland same wall
          if( ( abs(1/b12 -1) < 3*0.24 && q2 < 15 ) ||  1/b12 > 1 ){
            ID_ClustHit[jj] = -1;
          }
        }
        else if( ((id1>400 && id1<=460)&&(id2>400&&id2<=460)) || (id1>460&&id2>460)  ){ 
          // Nebula same wall
          if(( abs((1/b12)-1) < 3*0.49 && q2<15 ) || 1/b12 > 1 ){
            ID_ClustHit[jj]=-1;
          }
        }
        else{
          // Different Wall a
          if(Beta1/b12 > 1 || (abs(1/b12 - 1) < 3*0.17 && q2 < 15))
            ID_ClustHit[jj]=-1;
        }
        }
      }

    for(unsigned int i = 0; i < sizeHit; i++){
      if(ID_ClustHit[i]!=-1)
        mapOfClust[ID_ClustHit[i]].push_back(m_SortedID[i]);
    }

    //Put first hit in time of each cluster in a new map mapOfHead and remake 
    //the numbering of clusters (1,2,3...)
    static map< unsigned int, unsigned int > mapOfHead;
    static map< unsigned int, double > mapOfQtot;
    /* static map< unsigned int, unsigned int > mapOfClust; */
    static unsigned int NbrOfClust;
    mapOfHead.clear(), mapOfQtot.clear();
    NbrOfClust = mapOfClust.size();

    static unsigned int count, count2;
    count=1; m_HeadClust.clear();
    for(auto itr = mapOfClust.begin(); itr != mapOfClust.end(); itr++){
      mapOfHead[count] = itr->second[0];
      mapOfQtot[count] = 0;
      m_HeadClust.push_back(itr->second[0]);
      count++;   
    } 
    //return vector of real neutrons IDs 
    m_Neutrons.clear();
    m_QtotNeut.clear();
    m_QtotClust.clear();
    for(auto itr = mapOfHead.begin(); itr != mapOfHead.end(); itr++){
      m_Neutrons.push_back(itr->second);
      m_HeadClust.push_back(itr->second);
    }
    for(auto itr = mapOfQtot.begin(); itr != mapOfQtot.end(); itr++){
      m_QtotNeut.push_back(itr->second);
      m_QtotClust.push_back(itr->second);
    }
    return m_Neutrons;
  }

  vector<int> CrossTalk::GetSortedHits(){
    return m_SortedID;
  }
  vector<TGraph> CrossTalk::GetClusters(){
    return Clusters_2D;
  }
  TGraph2DErrors *CrossTalk::Get3DClusters(){
    return Clusters_3D;
  }
  vector<int> CrossTalk::GetHeadClust(){
    return m_HeadClust;
  }
  vector<double> CrossTalk::GetQtotClust(){
    return m_QtotClust;
  }
  vector<double> CrossTalk::GetQtotNeut(){
    return m_QtotNeut;
  }
