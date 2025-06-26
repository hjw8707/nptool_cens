/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien Matta  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : May  2021                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Samurai field map data                                   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "SamuraiFieldMap.h"
#include "NPPhysicalConstants.h"
#include "Math/Factory.h"
using namespace NPUNITS;

#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

ClassImp(SamuraiFieldMap);

SamuraiFieldMap::SamuraiFieldMap(){
  m_BrhoScan=NULL;
  m_min=ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad"); 
  m_func=ROOT::Math::Functor(this,&SamuraiFieldMap::Delta,1); 
  m_min->SetFunction(m_func);
  m_min->SetPrintLevel(-1);
  //default values
  m_StepTime = 1.*nanosecond;//propagation time interval size
  m_Limit = 1000;//maximum number of steps before giving up


}
////////////////////////////////////////////////////////////////////////////////
double SamuraiFieldMap::Delta(const double* parameter){
  static vector<TVector3>pos ;
  static TVector3 diff;
  //pos =Propagate(parameter[0],m_FitPosFDC0,m_FitDirFDC0,false); 
  pos =Propagate(parameter[0],m_FitPosFDC0,m_FitDirFDC0,true); 
  // Move the fdc2 pos from lab frame to fdc2 frame 
  pos.back().RotateY(-m_fdc2angle+m_angle);
  //double d = (pos.back().X()-m_FitPosFDC2.X())*(pos.back().X()-m_FitPosFDC2.X());
  // return d;
  diff = pos.back()-m_FitPosFDC2;
  return diff.Mag2();
}

////////////////////////////////////////////////////////////////////////////////
double SamuraiFieldMap::FindBrho(TVector3 p_fdc0,TVector3 d_fdc0,TVector3 p_fdc2,TVector3 d_fdc2){
  m_FitPosFDC0=p_fdc0;
  m_FitDirFDC0=d_fdc0;
  m_FitPosFDC2=p_fdc2;
  m_FitDirFDC2=d_fdc2;

  if(!m_BrhoScan)
    BrhoScan(1,10,0.1);
  // do a first guess based on fdc2 pos
  double b0[1] ={m_BrhoScan->Eval(p_fdc2.X())};
  //cout << "First guess Brho " << b0[0] << " "; //endl;


  m_min->Clear();
  m_min->SetPrecision(1e-6);
  m_min->SetMaxFunctionCalls(1000);
  m_min->SetLimitedVariable(0,"B",b0[0],0.1,1,10);
  m_min->Minimize();
  return m_min->X()[0];
}

////////////////////////////////////////////////////////////////////////////////
TGraph* SamuraiFieldMap::BrhoScan(double min, double max,double step){
  if(m_BrhoScan)
    delete m_BrhoScan;
  m_BrhoScan=new TGraph;
  unsigned int size = (max-min)/step;
  m_BrhoScan->Set(size);
  unsigned int i=0;
  TVector3 p(0,0,-3500);
  TVector3 d(0,0,1);
  p.RotateY(m_angle);
  d.RotateY(m_angle);
  for(double b = min ; b < max ; b+=step){
    vector<TVector3> pos= Propagate(b,p,d,false);
    pos.back().RotateY(-m_fdc2angle);
    m_BrhoScan->SetPoint(i++,pos.back().X(),b); 
  }
  m_BrhoScan->Sort();
  return m_BrhoScan;
}

////////////////////////////////////////////////////////////////////////////////
TVector3 SamuraiFieldMap::PropagateToFDC2(TVector3 pos, TVector3 dir){
  // go to FDC2 frame reference
  pos.RotateY(-m_fdc2angle);
  dir.RotateY(-m_fdc2angle);

  double deltaZ=m_fdc2R-pos.Z();
  dir*=deltaZ/dir.Z();
  pos+=dir;
  pos.SetX(pos.X());
  pos.RotateY(m_fdc2angle);
  return pos;
}

////////////////////////////////////////////////////////////////////////////////
std::vector< TVector3 > SamuraiFieldMap::Propagate(double Brho, TVector3 pos, TVector3 dir,bool store){
  pos.RotateY(m_angle);
  dir.RotateY(m_angle);
  dir=dir.Unit();
  // Property of a particle with the correct Brho:
  // We assume a 4He to compute v
  // The choice of the particle is of no importance
  static NPL::Particle N("4He");
  N.SetBrho(Brho);

  // track result
  static std::vector< TVector3 > track;
  track.clear();

  // starting point of the track
  if(store){
    pos.RotateY(-m_angle);
    track.push_back(pos);
    pos.RotateY(m_angle);
  }

  dir=dir.Unit();
  static double r;
  r = sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
  
  // number of step taken
  static unsigned int count;
  count = 0;

  // First propagate to r_max with one line
  while(r>m_Rmax && count<m_Limit){
    pos+=(r-m_Rmax)/cos(dir.Theta())*dir.Unit();
    r= 1.01*sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
  }

  if(r<=m_Rmax){ // success
    if(store){
      pos.RotateY(-m_angle);
      track.push_back(pos);
      pos.RotateY(m_angle);
    }
  }
  else {// failure
    return track;
  }

  static TVector3 xk1,xk2,xk3,xk4; // position
  static TVector3 pk1,pk2,pk3,pk4; // impulsion
  static TVector3 imp;
  static double K,m,P,px,py,pz;
  K = N.GetEnergy(); // kinetic energy
  m = N.Mass(); // mc2
  P = sqrt(K*K+2*K*m)/c_light; // P
  px = P*dir.X();//px
  py = P*dir.Y();//py
  pz = P*dir.Z();//pz
  imp = P*dir;

  while(r<=m_Rmax && count < m_Limit){
    func(N, pos           , imp            , xk1, pk1);
    func(N, pos+xk1*(m_StepTime/2.), imp+pk1*(m_StepTime/2.) , xk2, pk2);
    func(N, pos+xk2*(m_StepTime/2.), imp+pk2*(m_StepTime/2.) , xk3, pk3);
    func(N, pos+xk3*m_StepTime     , imp+pk3*m_StepTime      , xk4, pk4);
    pos +=(xk1+2*xk2+2*xk3+xk4)*(m_StepTime/6.);

    imp +=(pk1+2*pk2+2*pk3+pk4)*(m_StepTime/6.);
    if(store){
      pos.RotateY(-m_angle);
      track.push_back(pos);
      pos.RotateY(m_angle);
    }
    r = sqrt(pos.X()*pos.X()+pos.Z()*pos.Z());
    count++;
  }

  imp=imp.Unit();
  pos = PropagateToFDC2(pos, imp);
  pos.RotateY(-m_angle);
  track.push_back(pos);

  return track;

}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::func(NPL::Particle& N, TVector3 pos, TVector3 imp, TVector3& new_pos, TVector3& new_imp){
  static double px,py,pz,vx,vy,vz,Bx,By,Bz,q,P2,D,m2c4;
  static vector<double> B; 
  px=imp.X(); 
  py=imp.Y();
  pz=imp.Z();

  P2=imp.Mag2(); // P2
  m2c4 = N.Mass()*N.Mass();
  D=sqrt(m2c4+P2*c_squared); // sqrt(m2c4+P2c2)
  vx=px*c_squared/D;// pxc * c / D = pxc2/D
  vy=py*c_squared/D;
  vz=pz*c_squared/D;
  new_pos.SetX(vx);
  new_pos.SetY(vy);
  new_pos.SetZ(vz);
  B = InterpolateB(pos);
  Bx= B[0]; 
  By= B[1];
  Bz= B[2];
  q = N.GetZ()*eplus; // issue with the tesla/coulomb definition
  new_imp.SetX(q*(vy*Bz-vz*By));// q*pyc2*Bz/D -q*pzc2*By/D
  new_imp.SetY(q*(vz*Bx-vx*Bz));
  new_imp.SetZ(q*(vx*By-vy*Bx));
}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadMap(double angle,std::string file,unsigned int bin){
  m_bin=bin;
  m_angle=angle;
  if(file.find(".bin")!=std::string::npos)
    LoadBinary(file);
  else
    LoadAscii(file);
}
////////////////////////////////////////////////////////////////////////////////
std::vector<double> SamuraiFieldMap::GetB(std::vector<double>& pos){
  static vector<double> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;

  if(pos[0]<0)
    pos[0] = -pos[0];

  if(pos[2]<0)
    pos[2] = -pos[2];

  auto it=m_field.find(pos);
  if(it!=m_field.end()){
    return it->second;
  }
  else 
    return nullv;
}
////////////////////////////////////////////////////////////////////////////////
std::vector<double> SamuraiFieldMap::InterpolateB(const std::vector<double>& pos){
  static vector<double> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;

  if(pos[0]>0)
    x = pos[0];
  else
    x = -pos[0];

  y = pos[1];

  if(pos[2]>0)
    z = pos[2];
  else
    z = -pos[2];

  // out of bound 
  if(x<m_x_min || x>m_x_max)
    return nullv;
  if(y<m_y_min || y>m_y_max)
    return nullv;
  if(z<m_z_min || z>m_z_max)
    return nullv;

  double x0 = (double)((int)(x)/m_bin)*m_bin;
  if(x<=x0)
    x0=(double)((int)(x-m_bin)/m_bin)*m_bin; 

  double x1 = (double)((int)(x)/m_bin)*m_bin;
  if(x>=x1)
    x1=(double)((int)(x+m_bin)/m_bin)*m_bin; 

  double y0 = (double)((int)(y)/m_bin)*m_bin;
  if(y<=y0)
    y0=(double)((int)(y-m_bin)/m_bin)*m_bin; 

  double y1 = (double)((int)(y)/m_bin)*m_bin;
  if(y>=y1)
    y1=(double)((int)(y+m_bin)/m_bin)*m_bin; 

  double z0 = (double)((int)(z)/m_bin)*m_bin;
  if(z<=z0)
    z0=(double)((int)(z-m_bin)/m_bin)*m_bin; 

  double z1 = (double)((int)(z)/m_bin)*m_bin;
  if(z>=z1)
    z1=(double)((int)(z+m_bin)/m_bin)*m_bin; 

  // Apply symmetry
  // corrected coordinate for the symmetrised map
  double xx0,xx1,zz0,zz1;
  if(x0<0)
    xx0 = -x0;
  else 
    xx0 = x0;
  if(z0<0)
    zz0 = -z0;
  else
    zz0 = z0;
  if(x1<0)
    xx1 = -x1;
  else
    xx1 = x1;
  if(z1<0)
    zz1 = -z1;
  else
    zz1 = z1;

  //vector<double> X={xm,ym,zm};
  vector<double> X000={xx0,y0,zz0};
  vector<double> X111={xx1,y1,zz1};
  vector<double> X100={xx1,y0,zz0};
  vector<double> X010={xx0,y1,zz0};
  vector<double> X001={xx0,y0,zz1};
  vector<double> X101={xx1,y0,zz1};
  vector<double> X011={xx0,y1,zz1};
  vector<double> X110={xx1,y1,zz0};

  vector<double> C000;
  vector<double> C111;
  vector<double> C100;
  vector<double> C010;
  vector<double> C001;
  vector<double> C101;
  vector<double> C011;
  vector<double> C110;

  if(m_field.lower_bound(X000)!=m_field.end())
    C000=m_field.lower_bound(X000)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X111)!=m_field.end())
    C111=m_field.lower_bound(X111)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X100)!=m_field.end())
    C100=m_field.lower_bound(X100)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X010)!=m_field.end())
    C010=m_field.lower_bound(X010)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X001)!=m_field.end())
    C001=m_field.lower_bound(X001)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X101)!=m_field.end())
    C101=m_field.lower_bound(X101)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X011)!=m_field.end())
    C011=m_field.lower_bound(X011)->second;
  else 
    return nullv;

  if(m_field.lower_bound(X110)!=m_field.end())
    C110=m_field.lower_bound(X110)->second;
  else 
    return nullv;

  double xd = (x-x0)/(x1-x0);    
  double yd = (y-y0)/(y1-y0);    
  double zd = (z-z0)/(z1-z0);    
  double alphaX = 1-xd;
  double alphaY = 1-yd;
  double alphaZ = 1-zd;
  // X 
  vector<double> C00 = {C000[0]*alphaX+C100[0]*xd,C000[1]*alphaX+C100[1]*xd,C000[2]*alphaX+C100[2]*xd};
  vector<double> C01 = {C001[0]*alphaX+C101[0]*xd,C001[1]*alphaX+C101[1]*xd,C001[2]*alphaX+C101[2]*xd};
  vector<double> C10 = {C010[0]*alphaX+C110[0]*xd,C010[1]*alphaX+C110[1]*xd,C010[2]*alphaX+C110[2]*xd};
  vector<double> C11 = {C011[0]*alphaX+C111[0]*xd,C011[1]*alphaX+C111[1]*xd,C011[2]*alphaX+C111[2]*xd};
  // Y
  vector<double> C0  = {C00[0] *alphaY+C10[0] *yd,C00[1] *alphaY+C10[1] *yd,C00[2] *alphaY+C10[2] *yd};
  vector<double> C1  = {C01[0] *alphaY+C11[0] *yd,C01[1] *alphaY+C11[1] *yd,C01[2] *alphaY+C11[2] *yd};
  // Z
  vector<double> res = {C0[0]  *alphaZ+C1[0]  *zd,C0[1]  *alphaZ+C1[1]  *zd,C0[2]  *alphaZ+C1[2]  *zd};

  return res;
}
////////////////////////////////////////////////////////////////////////////////
std::vector<double> SamuraiFieldMap::InterpolateBalt(const std::vector<double>& pos){
  static vector<double> nullv ={0,0,0};
  // the map is only 1/4 of the detector so we apply symetrie:
  double x,y,z ;

  if(pos[0]>0)
    x = pos[0];
  else
    x = -pos[0];

  y = pos[1];

  if(pos[2]>0)
    z = pos[2];
  else
    z = -pos[2];

  // out of bound 
  if(x<m_x_min || x>m_x_max)
    return nullv;
  if(y<m_y_min || y>m_y_max)
    return nullv;
  if(z<m_z_min || z>m_z_max)
    return nullv;

  double xm = (double)((int)x/m_bin*m_bin);
  double ym = (double)((int)y/m_bin*m_bin);
  double zm = (double)((int)z/m_bin*m_bin);

  vector<double> p0={xm,ym,zm};
  vector<double> p1={xm+m_bin,ym,zm};
  vector<double> p2={xm,ym+m_bin,zm};
  vector<double> p3={xm,ym,zm+m_bin};
  vector<double> p4={xm-m_bin,ym,zm};
  vector<double> p5={xm,ym-m_bin,zm};
  vector<double> p6={xm,ym,zm-m_bin};

  vector<map<vector<double>,vector<double>>::iterator> it=
  { m_field.lower_bound(p0),
    m_field.lower_bound(p1),m_field.lower_bound(p2),m_field.lower_bound(p3),
    m_field.lower_bound(p4),m_field.lower_bound(p5),m_field.lower_bound(p6)};

  double Bx=0;
  double By=0;
  double Bz=0;
  double totalW=0;
  auto end=m_field.end();
  unsigned int size = it.size();
  for(unsigned int i = 0 ; i < size; i++){
    if(it[i]!=end){
      double d = 1e-6+sqrt( 
          (x-it[i]->first[0])*(x-it[i]->first[0])+
          (y-it[i]->first[1])*(y-it[i]->first[1])+
          (z-it[i]->first[2])*(z-it[i]->first[2]));

      Bx+=it[i]->second[0]/(d*d);
      By+=it[i]->second[1]/(d*d);
      Bz+=it[i]->second[2]/(d*d);
      totalW+=1./(d*d);
    }
  }
  vector<double> res = {Bx/totalW,By/totalW,Bz/totalW};
  return res;
}

////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadAscii(std::string file){
  ifstream in(file.c_str());
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Ascii Samurai field map " << file << endl; 
  double x,y,z,Bx,By,Bz;

  m_x_max=m_y_max=m_z_max=-1e32;
  m_x_min=m_y_min=m_z_min=1e32;
  unsigned int  count =0 ;

  // ignore 8 first line 
  string buffer;
  for(unsigned int i = 0 ; i < 8 ; i++){
    getline(in,buffer);
  }

  while(in >> x >> y >> z >> Bx >> By >> Bz){
    if(++count%50000==0)
      cout << "\r  - Loading " << count << " values " << flush; 
    vector<double> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<double> B = {Bx,By,Bz};
    m_field[p]=B;
    if(x<m_x_min)
      m_x_min=x;
    if(x>m_x_max)
      m_x_max=x;  
    if(y<m_y_min)
      m_y_min=y;
    if(y>m_y_max)
      m_y_max=y;  
    if(z<m_z_min)
      m_z_min=z;
    if(z>m_z_max)
      m_z_max=z;  
  }

  m_Rmax=m_x_max;
  cout << "\r  - " << count << " values loaded" << endl; 
  cout << "  - min(" << m_x_min <<";"<< m_y_min <<";" << m_z_min<< ") max(" << m_x_max <<";"<< m_y_max <<";" << m_z_max<< ")" << endl; 
  cout << "  - Rmax = " << m_Rmax << endl;
  in.close();
}
////////////////////////////////////////////////////////////////////////////////
void SamuraiFieldMap::LoadBinary(std::string file){
  ifstream in(file.c_str(),std::ifstream::binary);
  if(!in.is_open()){
    cout << "Error: failed to load samurai field map " << file << endl;
    exit(1);
  }

  cout << "//////// Loading Binary Samurai field map " << file << endl; 
  double x,y,z,Bx,By,Bz;

  m_x_max=m_y_max=m_z_max=-1e32;
  m_x_min=m_y_min=m_z_min=1e32;
  unsigned int  count =0 ;
  while(!in.eof()){

    if(++count%50000==0)
      cout << "\r  - Loading " << count << " values " << flush; 

    in.read((char*)&x,sizeof(x));
    in.read((char*)&y,sizeof(y));
    in.read((char*)&z,sizeof(z));
    in.read((char*)&Bx,sizeof(Bx));
    in.read((char*)&By,sizeof(By));
    in.read((char*)&Bz,sizeof(Bz));
    vector<double> p = {x,y,z};
    Bx*=tesla;
    By*=tesla;
    Bz*=tesla;
    vector<double> B = {Bx,By,Bz};
    m_field[p]=B;
    if(x<m_x_min)
      m_x_min=x;
    if(x>m_x_max)
      m_x_max=x;  
    if(y<m_y_min)
      m_y_min=y;
    if(y>m_y_max)
      m_y_max=y;  
    if(z<m_z_min)
      m_z_min=z;
    if(z>m_z_max)
      m_z_max=z;  
  }

  //default value for m_Rmax
  m_Rmax=m_x_max;

  cout << "\r  - " << count << " values loaded" << endl; 
  cout << "  - min(" << m_x_min <<";"<< m_y_min <<";" << m_z_min<< ") max(" << m_x_max <<";"<< m_y_max <<";" << m_z_max<< ")" << endl; 
  cout << "  - Rmax = " << m_Rmax << endl;
  in.close();
}
