#ifndef GladFieldMap_h
#define GladFieldMap_h
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2022                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold a basic field map and is used for Brho reconstruction    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>
#include <map>
using namespace std;

// ROOT
#include "TObject.h"
#include "TGraph.h"
#include "TVector3.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"

#include "NPParticle.h"

#define NX 80;
#define NY 40;
#define NZ 120;

class GladFieldMap{
  public: 
    GladFieldMap();
    ~GladFieldMap();
    
  private:
    // GLAD parameters
    map<vector<double>, vector<double>> m_field;
    vector<double> m_Bx;
    vector<double> m_By;
    vector<double> m_Bz;
    double m_Bmax;
    double m_B;
    double m_Current;
    TVector3 m_Glad_Entrance;
    TVector3 m_Glad_TurningPoint;
    double m_Leff[81][41];
    double m_Tilt;
    int m_Nx;
    int m_Ny;
    int m_Nz;
    double m_x_min;
    double m_y_min;
    double m_z_min;
    double m_x_max;
    double m_y_max;
    double m_z_max;
    double m_R_max;
    double m_Scale;
    int m_bin;
  private:
    // MWPC3 paramters
    double m_CentralTheta;
    TVector3 m_MWPC3_POS;
    double m_Angle_MWPC3;
  private:
    // Runge-Kunta 4 paramaters
    double m_dt;
    int m_Limit;

  private:
    TVector3 m_InitPos;
    TVector3 m_InitDir;
    TVector3 m_FinalPos;

  public:
    void SetGladEntrance(double x, double y, double z){m_Glad_Entrance = TVector3(x,y,z);}
    void SetGladTurningPoint(double x, double y, double z){m_Glad_TurningPoint = TVector3(x,y,z);}
    void SetLeff(int ix, int iy, double val) {m_Leff[ix][iy] = val;}
    void SetGladTiltAngle(double val) {m_Tilt = val;}
    void SetScale(double val) {m_Scale = val;}
    void SetCurrent(double val){
      m_Current = val;
      m_Scale = m_Current/3583.81;
      m_B = 2.2*m_Scale;
    }
    void SetBin(double val) {m_bin = val;}
    void SetTimeStep(double val) {m_dt = val;}

    void SetCentralTheta(double val) {m_CentralTheta = val;}
    void Set_MWPC3_Position(double x, double y, double z) {m_MWPC3_POS = TVector3(x,y,z);}
  
    void SetPropagationTimeInterval(double val) {m_dt = val;}
    void SetLimit(int val) {m_Limit = val;}

    void SetInitPos(TVector3 Pos) {m_InitPos = Pos;}
    void SetInitDir(TVector3 Dir) {m_InitPos = Dir;}

  public:
    double GetLeff(int ix, int iy) {return m_Leff[ix][iy];}
    double GetGladTiltAngle() {return m_Tilt;}
    TVector3 GetGladEntrance() {return m_Glad_Entrance;}
    TVector3 GetGladTurningPoint() {return m_Glad_TurningPoint;}
    double GetB() {return m_B;}
    double GetXmin() {return m_x_min;}
    double GetXmax() {return m_x_max;}
    double GetYmin() {return m_y_min;}
    double GetYmax() {return m_y_max;}
    double GetZmin() {return m_z_min;}
    double GetZmax() {return m_z_max;}
    double GetCentralTheta() {return m_CentralTheta;}
    double GetBin() {return m_bin;}
    double GetTimeStep() {return m_dt;}
    TVector3 Get_MWPC3_Position() {return m_MWPC3_POS;}
     
  public:
    void LoadMap(string filename);
    vector<double> InterpolateB(const vector<double>& pos);
    inline vector<double> InterpolateB(const TVector3& pos){
      vector<double> position={(double)pos.X(),(double)pos.Y(),(double)pos.Z()};
      return InterpolateB(position);
    };
    double GetB(TVector3 localpoint, string field_component);
    TGraph* BrhoScan(double Brho_min, double Brho_max, double Brho_step, TVector3 pos, TVector3 dir);
    TVector3 CalculateIntersectionPoint(vector<TVector3> vPos);
    vector<TVector3> Propagate(double Brho, TVector3 Pos, TVector3 Dir, bool store);
    double GetFlightPath(TVector3 vStart, double Brho, TVector3 Pos, TVector3 Dir);
    TVector3 PropagateToMWPC(TVector3 pos, TVector3 dir);
    void func(NPL::Particle& N, TVector3 Pos, TVector3 Imp, TVector3& xk, TVector3& pk);
    double FindBrho(TVector3 Pos_init, TVector3 Dir_init, TVector3 Pos_final);

  private:
    TGraph* m_BrhoScan;
    ROOT::Math::Minimizer* m_min;
    ROOT::Math::Functor m_func;
    double Delta(const double* parameter);
    bool IsInside(TVector3 localpoint, int& ix, int& iy, int& iz, double& dx, double& dy, double& dz);
    double Interpolate(double dx, double dy, double dz);

    double m_C0[2][2][2]; //!
    double m_C1[2][2]; //!
    double m_C2[2]; //!

  private:

    ClassDef(GladFieldMap,1)  
};

#endif
