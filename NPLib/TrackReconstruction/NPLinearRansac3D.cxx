/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author :  Adrien Matta  contact address: matta@lpccaen.in2p3.fr  *
 *                                                                           *
 * Creation Date   : April 2021                                              *
 * Last update     : April 2021                                              *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Find linear tracks using the ransac method                               *
 *****************************************************************************/

#include "NPLinearRansac3D.h"
#include<iostream>
using namespace std;
using namespace NPL;


ClassImp(LinearCluster3D);

ClassImp(LinearRansac3D);
////////////////////////////////////////////////////////////////////////////////
void LinearCluster3D::LinearFit() {

  unsigned int sizeI = m_index.size();
  double min_distance = 1e20;
  unsigned int min_i = 0;
  unsigned int min_j = 0;

  double meanDX = 0;
  double meanDY = 0;
  double meanDZ = 0;
  double totalW = 0;
  double meanX = 0;
  double meanY = 0;
  double meanZ = 0;
  double w;


  for (unsigned int i = 0; i < sizeI; i++) {
    TVector3 Pi((*m_X)[m_index[i]], (*m_Y)[m_index[i]], (*m_Z)[m_index[i]]);
    for (unsigned int j = i + 1; j < sizeI; j++) {
      double dij = 0;
      TVector3 Pj((*m_X)[m_index[j]], (*m_Y)[m_index[j]], (*m_Z)[m_index[j]]);
      // compute the average distance of all point to this line
      for (unsigned int p = 0; p < sizeI; p++) {
        TVector3 Pp((*m_X)[m_index[p]], (*m_Y)[m_index[p]], (*m_Z)[m_index[p]]);
        dij += MinimumDistancePointLine(Pi, Pj, Pp);
      }

      w = 1. / dij;
      totalW += w;

      meanX += w * ((*m_X)[m_index[i]]);
      meanY += w * ((*m_Y)[m_index[i]]);
      meanZ += w * ((*m_Z)[m_index[i]]);

      meanDX += w * ((*m_X)[m_index[i]] - (*m_X)[m_index[j]]);
      meanDY += w * ((*m_Y)[m_index[i]] - (*m_Y)[m_index[j]]);
      meanDZ += w * ((*m_Z)[m_index[i]] - (*m_Z)[m_index[j]]);
    }
  }

  meanDX /= totalW;
  meanDY /= totalW;
  meanDZ /= totalW;
  meanX /= totalW;
  meanY /= totalW;
  meanZ /= totalW;

  m_P0 = TVector3(meanX, meanY, meanZ);
  m_Dir = (TVector3(meanDX, meanDY, meanDZ)).Unit();
  m_isFitted = true;
  // Cluster dir is arbitrary. we choose to always have positive Z dir
  if (m_Dir.Z() < 0)
    m_Dir = -m_Dir;
};

void LinearCluster3D::FastLinearFit() {

double xx=0, yy=0, zz=0, xx2=0, yy2=0, zz2=0, xy=0, xz=0, yz=0, qq=0;
unsigned int sizeI = m_index.size();
	for(unsigned int i = 0; i < sizeI; i++){
		xx+=(*m_X)[m_index[i]];
		yy+=(*m_Y)[m_index[i]];
		zz+=(*m_Z)[m_index[i]];
		xx2+=(*m_X)[m_index[i]]*(*m_X)[m_index[i]];
		yy2+=(*m_Y)[m_index[i]]*(*m_Y)[m_index[i]];
		zz2+=(*m_Z)[m_index[i]]*(*m_Z)[m_index[i]];
		xy+=(*m_X)[m_index[i]]*(*m_Y)[m_index[i]];
		xz+=(*m_X)[m_index[i]]*(*m_Z)[m_index[i]];
		yz+=(*m_Y)[m_index[i]]*(*m_Z)[m_index[i]];
		qq+=1;
	}

  double Xm = xx/qq, Ym = yy/qq, Zm = zz/qq;
  double Sxx = -Xm*Xm+xx2/qq, Syy = -Ym*Ym+yy2/qq, Szz = -Zm*Zm+zz2/qq;
  double Sxy = -Xm*Ym+xy/qq, Sxz = -Xm*Zm+xz/qq, Syz = -Ym*Zm+yz/qq;
  double theta = 0.5*atan(2*Sxy/(Sxx-Syy));

  double ct = cos(theta), st = sin(theta), ct2 = ct*ct, st2 = st*st, cst = cos(theta)*sin(theta);
  double K11 = (Syy+Szz)*ct2+(Sxx+Szz)*st2-2*Sxy*cst;
  double K22 = (Syy+Szz)*st2+(Sxx+Szz)*ct2+2*Sxy*cst;
  //double K12 = -Sxy*(ct2-st2)+(Sxx-Syy)*cst; //=0 by construction
  double K10 = Sxz*ct+Syz*st;
  double K01 = -Sxz*st+Syz*ct;
  double K00 = Sxx+Syy;

  double C2 = -K00-K11-K22;
  double C1 = K00*K11+K00*K22+K11*K22-K01*K01-K10*K10;
  double C0 = K01*K01*K11+K10*K10*K22-K00*K11*K22;

  double p = C1-C2*C2/3, q = 2*C2*C2*C2/27-C1*C2/3+C0;
  double R = q*q/4+p*p*p/27;
  double dm2 = 0;
  if(R>0)
  {
    dm2 = -C2/3+pow(-q/2+sqrt(R), 1./3.)+pow(-q/2-sqrt(R), 1./3.);
  }
  else
  {
    double rho = sqrt(-p*p*p/27), phi = acos(-q/(2*rho));  double R = q*q/4+p*p*p/27;

    double sqrt3rho = pow(rho, 1./3.);
    double d1 = -C2/3+2*sqrt3rho*cos(phi/3);
    double d2 = -C2/3+2*sqrt3rho*cos((phi+2*M_PI)/3);
    double d3 = -C2/3+2*sqrt3rho*cos((phi+4*M_PI)/3);
    dm2 = std::min(std::min(d1, d2), d3);
  }
  double a = -K10*ct/(K11-dm2)+K01*st/(K22-dm2);
  double b = -K10*st/(K11-dm2)-K01*ct/(K22-dm2);

  double denom = 1./(1+a*a+b*b);
  double u = ((1+b*b)*Xm-a*b*Ym+a*Zm)*denom;
  double v = (-a*b*Xm+(1+a*a)*Ym+b*Zm)*denom;
  double w = (a*Xm+b*Ym+(a*a+b*b)*Zm)*denom;

  m_P0 = TVector3(Xm, Ym, Zm);
  m_Dir = (TVector3(Xm-u, Ym-v, Zm-w)).Unit();
}

void LinearCluster3D::FastLinearFitShort() {

double xx=0, yy=0, zz=0, xx2=0, yy2=0, zz2=0, xy=0, xz=0, yz=0, qq=0;
unsigned int sizeI = m_index_short.size();
	for(unsigned int i = 0; i < sizeI; i++){
		xx+=(*m_X)[m_index_short[i]];
		yy+=(*m_Y)[m_index_short[i]];
		zz+=(*m_Z)[m_index_short[i]];
		xx2+=(*m_X)[m_index_short[i]]*(*m_X)[m_index_short[i]];
		yy2+=(*m_Y)[m_index_short[i]]*(*m_Y)[m_index_short[i]];
		zz2+=(*m_Z)[m_index_short[i]]*(*m_Z)[m_index_short[i]];
		xy+=(*m_X)[m_index_short[i]]*(*m_Y)[m_index_short[i]];
		xz+=(*m_X)[m_index_short[i]]*(*m_Z)[m_index_short[i]];
		yz+=(*m_Y)[m_index_short[i]]*(*m_Z)[m_index_short[i]];
		qq+=1;
	}

  double Xm = xx/qq, Ym = yy/qq, Zm = zz/qq;
  double Sxx = -Xm*Xm+xx2/qq, Syy = -Ym*Ym+yy2/qq, Szz = -Zm*Zm+zz2/qq;
  double Sxy = -Xm*Ym+xy/qq, Sxz = -Xm*Zm+xz/qq, Syz = -Ym*Zm+yz/qq;
  double theta = 0.5*atan(2*Sxy/(Sxx-Syy));

  double ct = cos(theta), st = sin(theta), ct2 = ct*ct, st2 = st*st, cst = cos(theta)*sin(theta);
  double K11 = (Syy+Szz)*ct2+(Sxx+Szz)*st2-2*Sxy*cst;
  double K22 = (Syy+Szz)*st2+(Sxx+Szz)*ct2+2*Sxy*cst;
  //double K12 = -Sxy*(ct2-st2)+(Sxx-Syy)*cst; //=0 by construction
  double K10 = Sxz*ct+Syz*st;
  double K01 = -Sxz*st+Syz*ct;
  double K00 = Sxx+Syy;

  double C2 = -K00-K11-K22;
  double C1 = K00*K11+K00*K22+K11*K22-K01*K01-K10*K10;
  double C0 = K01*K01*K11+K10*K10*K22-K00*K11*K22;

  double p = C1-C2*C2/3, q = 2*C2*C2*C2/27-C1*C2/3+C0;
  double R = q*q/4+p*p*p/27;
  double dm2 = 0;
  if(R>0)
  {
    dm2 = -C2/3+pow(-q/2+sqrt(R), 1./3.)+pow(-q/2-sqrt(R), 1./3.);
  }
  else
  {
    double rho = sqrt(-p*p*p/27), phi = acos(-q/(2*rho));  double R = q*q/4+p*p*p/27;

    double sqrt3rho = pow(rho, 1./3.);
    double d1 = -C2/3+2*sqrt3rho*cos(phi/3);
    double d2 = -C2/3+2*sqrt3rho*cos((phi+2*M_PI)/3);
    double d3 = -C2/3+2*sqrt3rho*cos((phi+4*M_PI)/3);
    dm2 = std::min(std::min(d1, d2), d3);
  }
  double a = -K10*ct/(K11-dm2)+K01*st/(K22-dm2);
  double b = -K10*st/(K11-dm2)-K01*ct/(K22-dm2);

  double denom = 1./(1+a*a+b*b);
  double u = ((1+b*b)*Xm-a*b*Ym+a*Zm)*denom;
  double v = (-a*b*Xm+(1+a*a)*Ym+b*Zm)*denom;
  double w = (a*Xm+b*Ym+(a*a+b*b)*Zm)*denom;

  m_P0 = TVector3(Xm, Ym, Zm);
  m_Dir = (TVector3(Xm-u, Ym-v, Zm-w)).Unit();
}

void LinearCluster3D::CheckDirection(const TVector3 & vertex){
    if(m_FurthestPoint == TVector3(0, 0, 0)){
      double max_distance = 0;
      for(int i = 0; i < m_index.size(); i++)
      {
        TVector3 aPoint((*m_X)[m_index[i]], (*m_Y)[m_index[i]], (*m_Z)[m_index[i]]);
        double distance = (aPoint-vertex).Mag();
        if(distance > max_distance){
          max_distance = distance;
          m_FurthestPoint = aPoint;
        }
      }
    }
    if((m_FurthestPoint-vertex).X()*m_Dir.X() < 0) this->Inverse();
}

void LinearCluster3D::CalculateTrackLength(const TVector3 & vertex, const double & percentRange){
  TH1D* Hrange = new TH1D("Range", "range", 256, 0, 512);
  TH1D* Hrangecumule = new TH1D("RangeCumulated", "rangecumule", 256, 0, 512);

  double total_charge = 0;
  double cumulated_charge = 0;
  double max_distance = 0;


  for(int i = 0; i < m_index.size(); i++)
  {
    TVector3 aPoint((*m_X)[m_index[i]], (*m_Y)[m_index[i]], (*m_Z)[m_index[i]]);
    double distance = (aPoint-vertex).Mag();
    total_charge+=(*m_Q)[m_index[i]];
    Hrange->Fill(distance, (*m_Q)[m_index[i]]);
    if(distance > max_distance){
      max_distance = distance;
      m_FurthestPoint = aPoint;
    }
    if(distance < m_ShortTrackMaxLength){
      m_index_short.push_back(m_index[i]);
    }
  }
  m_Charge = total_charge;
  //Inverse the track if it is not going in the good direction
  if((m_FurthestPoint-vertex).X()*m_Dir.X() < 0) this->Inverse();

  for(int bin = 1; bin < Hrange->GetNbinsX(); bin++)
  {
    cumulated_charge+=Hrange->GetBinContent(bin);
    Hrangecumule->SetBinContent(bin, cumulated_charge/total_charge);
  }
  Int_t nbin = 1;
  while(Hrangecumule->GetBinContent(nbin) < percentRange && nbin < Hrangecumule->GetNbinsX()-1)
  {
    nbin++;
  }
  Double_t delta, deltabinap, deltabinav;
  delta = Hrangecumule->GetBinContent(nbin)-Hrangecumule->GetBinContent(nbin-1);
  deltabinap = Hrangecumule->GetBinContent(nbin)-percentRange;
  deltabinav = percentRange-Hrangecumule->GetBinContent(nbin-1);
  m_TrackLength = deltabinap*Hrangecumule->GetBinCenter(nbin-1)/delta+deltabinav*Hrangecumule->GetBinCenter(nbin)/delta;
  delete Hrange;
  delete Hrangecumule;
}

////////////////////////////////////////////////////////////////////////////////
vector<LinearCluster3D> LinearRansac3D::TreatEvent(vector<double>& X, vector<double>& Y, vector<double>& Z) {
  cluster_id.clear();
  clusters.clear();
  unsigned int sizeX = X.size();
  cluster_id.resize(sizeX, 0);
  m_cluster.clear();
  m_assigned.clear();

  if (sizeX < m_min_count)
    return clusters;

  unsigned int p1, p2;
  double d;
  TVector3 Vp1, Vp2, D, P;
  m_iteration_d = 0;
  m_iteration = 0;

  while (m_iteration++ < m_max_iteration/* && m_assigned.size() < sizeX*/) {
    LinearCluster3D cluster(&X, &Y, &Z);
    // take 2 distant point randomly that has not been match before
    d = 0;
    m_iteration_d = 0;
    while (d < 3 * m_match_distance && m_iteration_d++ < m_max_iteration) {
      p1 = rand.Integer(sizeX);
      p2 = rand.Integer(sizeX);

      Vp1.SetXYZ(X[p1], Y[p1], Z[p1]);
      Vp2.SetXYZ(X[p2], Y[p2], Z[p2]);
      D = Vp1 - Vp2;
      d = D.Mag();
      if (d > m_max_distance)
        d = 0;
    }

    // loop over all points
    for (unsigned int i = 0; i < sizeX; i++) {
      P.SetXYZ(X[i], Y[i], Z[i]);
      if (MinimumDistancePointLine(Vp1, Vp2, P) < m_match_distance) {
        cluster.AddIndex(i);
        m_assigned.insert(i);
      }
    }
    // insert the newly formed cluster
    if (cluster.size() > m_min_count) {
      m_cluster.insert(cluster);
    }

  } // while
	
  // loop over the cluster starting with the biggest
  unsigned int current_cluster = 0;
  unsigned int index;
  for (auto it = m_cluster.begin(); it != m_cluster.end(); ++it) {
    current_cluster++;
    //  cout << current_cluster << endl;
    unsigned int sizeC = (*it).size();
    unsigned int cluster_size = 0;
    for (unsigned int i = 0; i < sizeC; i++) {
      // Assigned cluster id to identified points
      unsigned int index = (*it).GetIndex(i);
      if (!cluster_id[index]) {
        cluster_id[index] = current_cluster;
        cluster_size++;
      }
    }
    if (cluster_size < m_min_count) {
      for (unsigned int i = 0; i < sizeC; i++) {
        unsigned int index = (*it).GetIndex(i);
        // remove the assigned point
        if (cluster_id[index] == current_cluster) {
          cluster_id[index] = 0;
        }
      }
    }
    else {
      clusters.push_back(*it);
    }
  }
  return clusters;
}

