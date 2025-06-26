/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr  *
 *                                                                           *
 * Creation Date  : February 2018                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Tina Treated  data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#include <sstream>
#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
#include "RootInput.h"
#include "RootOutput.h"
#include "TChain.h"
#include "TAsciiFile.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"
#include "NPInputParser.h"
#include "NPSystemOfUnits.h"
#include "TTinaPhysics.h"

using namespace NPUNITS;
using namespace TINA_LOCAL;
using namespace std;

ClassImp(TTinaPhysics)

///////////////////////////////////////////////////////////////////////////
TTinaPhysics::TTinaPhysics(){
    
    EventMultiplicity                  = 0;
    m_EventData                        = new TTinaData;
    m_PreTreatedData                   = new TTinaData;
    m_EventPhysics                     = this;
    m_Spectra                          = NULL;
    m_NumberOfTelescopes               = 0;
    m_NumberOfTTTPad                   = 0;
    m_MaximumStripMultiplicityAllowed  = 10;
    m_StripEnergyMatchingSigma         = 0.020;  // sf warning: check validity for Tina!
    m_StripEnergyMatchingNumberOfSigma = 3;      // sf warning: check validity for Tina!
    m_TTT_XE_RAW_Threshold = 0;                  // sf warning: used to be 8200
    m_TTT_YE_RAW_Threshold = 0;                  // sf warning: used to be 8200
    m_Pad_E_RAW_Threshold  = 0;                  // sf warning: used to be 8200
    m_TTT_XE_CAL_Threshold = 0;
    m_TTT_YE_CAL_Threshold = 0;
    m_Pad_E_CAL_Threshold  = 0;
    
    m_YY1_RE_RAW_Threshold = 0;                  // sf warning: used to be 8200
    m_YY1_SE_RAW_Threshold = 0;                  // sf warning: used to be 8200
    m_CsI_E_RAW_Threshold  = 0;                  // sf warning: used to be 8200
    m_YY1_RE_CAL_Threshold = 0;
    m_YY1_SE_CAL_Threshold = 0;
    m_CsI_E_CAL_Threshold  = 0;
    // m_Ignore_not_matching_CsI = false;
    m_Take_YE = true;                           // sf warning: check validity for Tina!
    m_Take_YT = true;                            // sf warning: check validity for Tina!
    m_Take_SE = false;
    m_Take_ST = true;
    m_Pad_Size       = 32;                       // sf warning: used for matching in Must, check validity for Tina!
    m_NumberOfStrips = 128;
    m_Face           = 97.22;
    // easm file gives 101.00 x 100.50 but micron data file 100.42 x 100.42 compare Tina.hh Taking the active area only.
    m_NumberOfRings  = 16;
    m_Face_Ring      = 80;
    
}

///////////////////////////////////////////////////////////////////////////
TTinaPhysics::~TTinaPhysics() {}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::AddDetector(TVector3 Pos, string Shape, double Translation){
    
  if (Shape == "TTTPad") {
    m_NumberOfTelescopes++;
    double StripPitch       = m_Face / m_NumberOfStrips;
    ///////////////////////////////////////////////////
    if(fabs(Pos.X())<1) Pos.SetX(0);
    if(fabs(Pos.Y())<1) Pos.SetY(0);
    if(fabs(Pos.Z())<1) Pos.SetZ(0);
    TVector3 Pos2D(fabs(Pos.X()),fabs(Pos.Y()),0);
    Pos2D=Pos2D.Unit();
    /*TVector3 C_X1_Y1(m_Face/2.,m_Face/2.-min(1.,fabs(Pos.X()))*m_Face,0);
      TVector3 C_X1_Y128(m_Face/2.-min(1.,fabs(Pos.X()))*m_Face,-m_Face/2.,0);
      TVector3 C_X128_Y1(-m_Face/2.+min(1.,fabs(Pos.X()))*m_Face,m_Face/2.,0);
      TVector3 C_X128_Y128(-m_Face/2.,-m_Face/2.+min(1.,fabs(Pos.X()))*m_Face,0);//*/
    TVector3 C_X1_Y1(m_Face/2.,m_Face/2.,0);
    TVector3 C_X1_Y128(-m_Face/2.+min(1.,fabs(Pos.X()))*m_Face,m_Face/2.-min(1.,fabs(Pos.X()))*m_Face,0);
    TVector3 C_X128_Y1(m_Face/2.-min(1.,fabs(Pos.X()))*m_Face,-m_Face/2.+min(1.,fabs(Pos.X()))*m_Face,0);
    TVector3 C_X128_Y128(-m_Face/2.,-m_Face/2.,0);
    TVector3 W(0,0,1);
        
    //    C_X1_Y1.Rotate(-TMath::Pi()/2.,W.Cross(Pos2D));
    //    C_X1_Y128.Rotate(-TMath::Pi()/2.,W.Cross(Pos2D));
    //    C_X128_Y1.Rotate(-TMath::Pi()/2.,W.Cross(Pos2D));
    //    C_X128_Y128.Rotate(-TMath::Pi()/2.,W.Cross(Pos2D));
    C_X1_Y1+=Pos;
    C_X1_Y128+=Pos;
    C_X128_Y1+=Pos;
    C_X128_Y128+=Pos;
        
    //vector U on detector face parallel to Y strips (NB: Y strips along X axis)
    TVector3 U = C_X128_Y1 - C_X1_Y1;
    double   Ushift = (U.Mag() - m_Face) / 2.;
    U = U.Unit();
    // vector V on detector face parallel to X strips
    TVector3 V      = C_X1_Y128 - C_X1_Y1;
    double   Vshift = (V.Mag() - m_Face) / 2.;
    V = V.Unit();//*/
        
        
    ///////////////////////////////////////////////////
    // define the corners in case of box geometry
    // TVector3 C_X1_Y1     = ...
    // TVector3 C_X1_Y128   = ...
    // TVector3 C_X128_Y1   = ...
    // TVector3 C_X128_Y128 = ...
    // vector U on detector face parallel to Y strips (NB: Y strips along X axis)
    // TVector3 U = C_X128_Y1 - C_X1_Y1;
    // double   Ushift = (U.Mag() -m_Face) / 2.;
    // U = U.Unit();
    // vector V on detector face parallel to X strips
    // TVector3 V      = C_X1_Y128 - C_X1_Y1;
    // double   Vshift = (V.Mag() - m_Face) / 2.;
    // V = V.Unit();
    // for wall geometry, simply use the following
    U = TVector3(1,0,0);
    V = TVector3(0,1,0);//
    // position vector of pixel defined by X Y strips
    TVector3 Pixel = TVector3(0,0,0);
    // position vector of X=1 Y=1 strip
    TVector3 Strip_1_1;
        
    vector<double> lineX;
    vector<double> lineY;
    vector<double> lineZ;
    vector<vector<double>> OneTelescopeStripPositionX;
    vector<vector<double>> OneTelescopeStripPositionY;
    vector<vector<double>> OneTelescopeStripPositionZ;
        
    // find position of (1,1) corner at bottom right
    // sf: png files define pixel (1,1) at top left
    // however (x,y) hit coordinates are transformed into index in NPImage.cxx
    // such that a hit in pixel (1,1) finds itself in Geant4 at bottom right
    // this is because (x,y) are defined in opposite way in 2D png vs 3D Geant4
    // png: o------>X     geant4:         Y
    //      |                             ^
    //      v                             |
    //      Y                     X<------o
        
    // TVector3 C_X1_Y1 = Pos - (U + V) * (m_Face/2.);
    //    Strip_1_1 = C_X1_Y1 + (U + V) * (StripPitch/2.);
    Strip_1_1 = C_X128_Y128 + (U + V) * (StripPitch/2.);
    TVector3 Strip_128_128 = C_X1_Y1 + (U + V) * (StripPitch/2.+StripPitch*m_NumberOfStrips);
        
    // Strip_1_1 += U * Ushift + V * Vshift;
    for (int i = 0; i < m_NumberOfStrips; ++i) {
      lineX.clear();
      lineY.clear();
      lineZ.clear();
      for (int j = 0; j < m_NumberOfStrips; ++j) {
	Pixel = Strip_1_1 + StripPitch * (i * U + j * V);
	lineX.push_back(Pixel.X());
	lineY.push_back(Pixel.Y());
	lineZ.push_back(Pixel.Z());
	if ( (i == 0 || i == (m_NumberOfStrips - 1))
	     && (j == 0 || j == (m_NumberOfStrips - 1))) {
	  std::cout << m_NumberOfTelescopes << " ";
	  std::cout << i << " ";
	  std::cout << j << ": (";
	  std::cout << Pixel.X() << ", ";
	  std::cout << Pixel.Y() << ", ";
	  std::cout << Pixel.Z() << ") " << std::endl;}
      }
      OneTelescopeStripPositionX.push_back(lineX);
      OneTelescopeStripPositionY.push_back(lineY);
      OneTelescopeStripPositionZ.push_back(lineZ);
    }
    m_StripPositionX.push_back(OneTelescopeStripPositionX);
    m_StripPositionY.push_back(OneTelescopeStripPositionY);
    m_StripPositionZ.push_back(OneTelescopeStripPositionZ);
        
    m_NumberOfTTTPad = m_NumberOfTelescopes;
  }
    
  else if (Shape == "YY1CsI") {
    m_NumberOfTelescopes++;

    double StripPitch       = m_Face_Ring / m_NumberOfRings;
    ///////////////////////////////////////////////////
    if(fabs(Pos.X())<1e-10) Pos.SetX(0);
    if(fabs(Pos.Y())<1e-10) Pos.SetY(0);
    if(fabs(Pos.Z())<1e-10) Pos.SetZ(0);
    TVector3 Pos2D(fabs(Pos.X()),fabs(Pos.Y()),0);
    Pos2D=Pos2D.Unit();
    int sign=1;
    if(Translation!=0) sign=Translation/fabs(Translation);
    TVector3 FirstRing(0,-sign*m_Face_Ring/2.,0) ;
    TVector3 LastRing (0,sign*m_Face_Ring/2.,0) ;
    TVector3 V(0,1,0);
    TVector3 W(0,0,1);
        
    //extract angle parameter
    TVector3 Trans (0,0,Translation);
    Pos-=Trans;
    double Theta = Pos.Theta();
    double Phi = Pos.Phi();
    //cout << "Theta=" << Theta << " Phi=" << Phi << endl;
    /*note
      Vector u=(sin(Phi),-cos(Phi),0)
      Vector v=(sin(Theta+90)cos(Phi),sin(Theta+90)sin(Phi),cos(Theta+90))=(cos(Theta)cos(Phi),cos(Theta)sin(Phi),-sin(Theta))
      Vector w=(sin(Theta)cos(Phi),sin(Theta)sin(Phi),cos(Theta))
      RotationMatrix Rot=(u,v,w)  */
        
    FirstRing.Rotate(Phi,W);
    LastRing.Rotate(Phi,W);
    V.Rotate(Phi,W);
    FirstRing.Rotate(Theta,V);
    LastRing.Rotate(Theta,V);
    W.Rotate(Theta,V);
    FirstRing.Rotate(-TMath::Pi()/2.,W);
    LastRing.Rotate(-TMath::Pi()/2.,W);
    /*
      TVector3 FirstRing_Rotated(cos(Theta)*cos(Phi)*(m_Face_Ring/2.),cos(Theta)*sin(Phi)*(m_Face_Ring/2.),-sin(Theta)*(m_Face_Ring/2.));
      TVector3 LastRing_Rotated(-cos(Theta)*cos(Phi)*(m_Face_Ring/2.),-cos(Theta)*sin(Phi)*(m_Face_Ring/2.),sin(Theta)*(m_Face_Ring/2.));
    */
        
        
    Pos+=Trans;
    //cout << "Pos.X()=" << Pos.X() << " Pos.Y()=" << Pos.Y() << " Pos.Z()=" << Pos.Z() << endl;
        
    FirstRing+=Pos;
    LastRing+=Pos;
        
    TVector3 FirstRing_Rotated_Translated = FirstRing;
    TVector3 LastRing_Rotated_Translated = LastRing;
        
        
    //vector U on detector face parallel to Y strips (NB: Y strips along X axis)
    TVector3 P = LastRing_Rotated_Translated - FirstRing_Rotated_Translated;
    //cout << "P.X()=" << P.X() << " P.Y()=" << P.Y() << " P.Z()=" << P.Z() << endl;
    double   Pshift = (P.Mag() - m_Face_Ring) / 2.;
    P = P.Unit();
        
    // position vector of pixel defined by X Y strips
    TVector3 Pixel = TVector3(0,0,0);
    // position vector of X=1 Y=1 strip
    TVector3 Strip_1;
        
    vector<double> lineX;
    vector<double> lineY;
    vector<double> lineZ;
        
        
    // TVector3 C_X1_Y1 = Pos - (U + V) * (m_Face/2.);
    Strip_1 = FirstRing_Rotated_Translated + P * (StripPitch/2.);
    TVector3 Strip_16 = FirstRing_Rotated_Translated + P * (StripPitch/2.+StripPitch*m_NumberOfRings);
    //cout << "Strip_1.Z()=" << Strip_1.Z() << " Strip_16.Z()=" << Strip_16.Z() << endl;
    // Strip_1_1 += U * Ushift + V * Vshift;
        
    lineX.clear();
    lineY.clear();
    lineZ.clear();
        
    for (int i = 0; i < m_NumberOfRings; ++i) {
      Pixel = Strip_1 + StripPitch * (i * P);
      lineX.push_back(Pixel.X());
      lineY.push_back(Pixel.Y());
      lineZ.push_back(Pixel.Z());
    }
    m_RingPositionX.push_back(lineX);
    m_RingPositionY.push_back(lineY);
    m_RingPositionZ.push_back(lineZ);
  }
    
  else {
    // sf warning: other detectors to be done
    cout << "Warning: no geometry implemented for " << Shape << endl;
  }
    
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::AddDetector(double R, double Theta, double Phi, double Translation, string Shape){
    TVector3 Pos(R*sin(Theta)*cos(Phi),R*sin(Theta)*sin(Phi),R*cos(Theta)+Translation);
    AddDetector(Pos,Shape,Translation);
} 

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::BuildSimplePhysicalEvent() {
    BuildPhysicalEvent();
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::BuildPhysicalEvent() {
    
    /* following is sample code from detector wizard to match energy and time
     unsigned int mysizeE = m_PreTreatedData->GetMultEnergy();
     unsigned int mysizeT = m_PreTreatedData->GetMultTime();
     for (UShort_t e = 0; e < mysizeE ; e++) {
     for (UShort_t t = 0; t < mysizeT ; t++) {
     if (m_PreTreatedData->GetE_DetectorNbr(e) == m_PreTreatedData->GetT_DetectorNbr(t)) {
     DetectorNumber.push_back(m_PreTreatedData->GetE_DetectorNbr(e));
     Energy.push_back(m_PreTreatedData->Get_Energy(e));
     Time.push_back(m_PreTreatedData->Get_Time(t));
     }
     }
     }*/
    
    // apply thresholds and calibration

    PreTreat();

    // m_multimatch = false;
    bool check_Pad = false;
    bool check_CsI = false;
    
    m_StripXEMult = m_PreTreatedData->GetTTTback_E_Mult();
    m_StripYEMult = m_PreTreatedData->GetTTTfront_E_Mult();
    m_StripXTMult = m_PreTreatedData->GetTTTback_T_Mult();
    m_StripYTMult = m_PreTreatedData->GetTTTfront_T_Mult();
    m_PadEMult    = m_PreTreatedData->GetPad_E_Mult();
    m_PadTMult    = m_PreTreatedData->GetPad_T_Mult();
    
    m_YY1RingEMult = m_PreTreatedData->GetYY1ring_E_Mult();
    m_YY1SectorEMult =m_PreTreatedData->GetYY1sector_E_Mult();
    m_YY1RingTMult =m_PreTreatedData->GetYY1ring_T_Mult();
    m_YY1SectorTMult =m_PreTreatedData->GetYY1sector_T_Mult();
    m_CsIEMult = m_PreTreatedData->GetCsI_E_Mult();
    m_CsITMult = m_PreTreatedData->GetCsI_T_Mult();
    
    // return the matching couples and set the type of matching. In case of a multiple match in one
    // detector, m_match_type is set to 1, else it is set to 2
    vector<TVector2> couple = Match_X_Y();
    unsigned int couple_size = couple.size();
    EventMultiplicity = std::min(m_StripXEMult, m_StripYEMult);
    
    //if(EventMultiplicity>1)

    for (unsigned int i = 0; i < couple_size; ++i) {
        
        if (m_match_type[i] == 1) {
            
            check_Pad = false;
            int N = m_PreTreatedData->GetTTTback_E_DetNbr(couple[i].X());
            int X = m_PreTreatedData->GetTTTback_E_StripNbr(couple[i].X());
            int Y = m_PreTreatedData->GetTTTfront_E_StripNbr(couple[i].Y());
            double TTT_X_E = m_PreTreatedData->GetTTTback_Energy(couple[i].X());
            double TTT_Y_E = m_PreTreatedData->GetTTTfront_Energy(couple[i].Y());
            
            // search for associated time
            double TTT_X_T = -1000;
            for (unsigned int t = 0; t < m_StripXTMult; ++t) {
                if (m_PreTreatedData->GetTTTback_T_StripNbr(couple[i].X())
                    == m_PreTreatedData->GetTTTback_T_StripNbr(t)) {
                    TTT_X_T = m_PreTreatedData->GetTTTback_Time(t);
                    break;
                }
            }
            double TTT_Y_T = -1000;
            for (unsigned int t = 0; t < m_StripYTMult; ++t) {
                if (m_PreTreatedData->GetTTTfront_T_StripNbr(couple[i].Y())
                    == m_PreTreatedData->GetTTTfront_T_StripNbr(t)) {
                    TTT_Y_T = m_PreTreatedData->GetTTTfront_Time(t);
                    break;
                }
            }
            
            for (unsigned int j = 0; j < m_PadEMult; ++j) {
                //     cout << "m_PreTreatedData->GetPad_E_DetNbr(j)=" << m_PreTreatedData->GetPad_E_DetNbr(j) << " N=" << N << endl;
                if (m_PreTreatedData->GetPad_E_DetNbr(j) == N) {
                    // if (Match_Si_CsI(X, Y, m_PreTreatedData->GetMMCsIECristalNbr(j))) {
                    Pad_N.push_back(m_PreTreatedData->GetPad_E_PadNbr(j));
                    Pad_E.push_back(m_PreTreatedData->GetPad_Energy(j));
                    Pad_T.push_back(-1000);
                    // look for associated time
                    for (unsigned int k = 0; k < m_PadTMult; ++k) {
                        // same crystal, same detector
                        // if (m_PreTreatedData->GetMMCsIECristalNbr(j)== m_PreTreatedData->GetMMCsITCristalNbr(k)) {
                        Pad_T[Pad_T.size()-1] = m_PreTreatedData->GetPad_Time(j);
                        break;
                        //}
                    }
                    check_Pad = true;
                    //}
                }
            }
            
            SquareTelescopeNumber.push_back(N);
            TTT_X.push_back(X);
            TTT_Y.push_back(Y);
            // store both values for checking purpose
            // TTT_EX.push_back(TTT_X_E);
            // TTT_TX.push_back(TTT_X_T);
            // TTT_EY.push_back(TTT_Y_E);
            // TTT_TY.push_back(TTT_Y_T);
            if (m_Take_YE)
                TTT_E.push_back(TTT_Y_E);
            else
                TTT_E.push_back(TTT_X_E);
            if (m_Take_YT)
                TTT_T.push_back(TTT_Y_T);
            else
                TTT_T.push_back(TTT_X_T);
            if (!check_Pad) {
                Pad_N.push_back(0);
                Pad_E.push_back(-1000);
                Pad_T.push_back(-1000);
            }
            EventType.push_back(1);
        }
    }
    
    
    for (unsigned int i = 0; i < m_PreTreatedData->GetYY1ring_E_Mult(); ++i) {
        
        check_CsI = false;
        int N = m_PreTreatedData->GetYY1ring_E_DetNbr(i);
        int R = m_PreTreatedData->GetYY1ring_E_StripNbr(i);
        int S = m_PreTreatedData->GetYY1sector_E_StripNbr(i);
        double YY1_R_E = m_PreTreatedData->GetYY1ring_Energy(i);
        double YY1_S_E = m_PreTreatedData->GetYY1sector_Energy(i);
        
        // search for associated time
        double YY1_R_T = -1000;
        for (unsigned int t = 0; t < m_YY1RingTMult; ++t) {
            if (m_PreTreatedData->GetYY1ring_T_StripNbr(i)
                == m_PreTreatedData->GetYY1ring_T_StripNbr(t)) {
                YY1_R_T = m_PreTreatedData->GetYY1ring_Time(t);
                break;
            }
        }
        double YY1_S_T = -1000;
        for (unsigned int t = 0; t < m_YY1SectorTMult; ++t) {
            if (m_PreTreatedData->GetYY1sector_T_StripNbr(i)
                == m_PreTreatedData->GetYY1sector_T_StripNbr(t)) {
                YY1_S_T = m_PreTreatedData->GetYY1sector_Time(t);
                break;
            }
        }
        
        for (unsigned int j = 0; j < m_CsIEMult; ++j) {
            //     cout << "m_PreTreatedData->GetPad_E_DetNbr(j)=" << m_PreTreatedData->GetPad_E_DetNbr(j) << " N=" << N << endl;
            if (m_PreTreatedData->GetCsI_E_DetNbr(j) == N) {
                // if (Match_Si_CsI(X, Y, m_PreTreatedData->GetMMCsIECristalNbr(j))) {
                CsI_N.push_back(m_PreTreatedData->GetCsI_E_DetNbr(j));
                CsI_E.push_back(m_PreTreatedData->GetCsI_Energy(j));
                CsI_T.push_back(-1000);
                // look for associated time
                for (unsigned int k = 0; k < m_CsITMult; ++k) {
                    // same crystal, same detector
                    // if (m_PreTreatedData->GetMMCsIECristalNbr(j)== m_PreTreatedData->GetMMCsITCristalNbr(k)) {
                    CsI_T[CsI_T.size()-1] = m_PreTreatedData->GetCsI_Time(j);
                    break;
                    //}
                }
                check_CsI = true;
                //}
            }
        }
        TrapezTelescopeNumber.push_back(N);
        YY1_R.push_back(R);
        YY1_S.push_back(S);
        // store both values for checking purpose
        // TTT_EX.push_back(TTT_X_E);
        // TTT_TX.push_back(TTT_X_T);
        // TTT_EY.push_back(TTT_Y_E);
        // TTT_TY.push_back(TTT_Y_T);
        if (m_Take_SE)
            YY1_E.push_back(YY1_S_E);
        else
            YY1_E.push_back(YY1_R_E);
        if (m_Take_ST)
            YY1_T.push_back(YY1_S_T);
        else
            YY1_T.push_back(YY1_R_T);
        if (!check_CsI) {
            CsI_N.push_back(0);
            CsI_E.push_back(-1000);
            CsI_T.push_back(-1000);
        }
        EventType.push_back(1);
        
    }//
    
    return;
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::PreTreat() {
    // apply thresholds and calibrations and test for disabled channels
    
    /* following is sample code from detector wizard to calibrate the data
     static CalibrationManager* Cal = CalibrationManager::getInstance();
     unsigned int mysize = m_EventData->GetMult();
     for (UShort_t i = 0; i < mysize ; ++i) {
     if (m_EventData->Get_Energy(i) > m_E_RAW_Threshold) {
     Double_t Energy = Cal->ApplyCalibration("Tina/ENERGY"+NPL::itoa(m_EventData->GetE_DetectorNbr(i)),m_EventData->Get_Energy(i));
     if (Energy > m_E_Threshold) {
     m_PreTreatedData->SetEnergy(m_EventData->GetE_DetectorNbr(i), Energy);
     }
     }
     }
     mysize = m_EventData->GetMultTime();
     for (UShort_t i = 0; i < mysize; ++i) {
     Double_t Time= Cal->ApplyCalibration("Tina/TIME"+NPL::itoa(m_EventData->GetT_DetectorNbr(i)),m_EventData->Get_Time(i));
     m_PreTreatedData->SetTime(m_EventData->GetT_DetectorNbr(i), Time);
     }*/
    
    ClearPreTreatedData();

    m_StripXEMult = m_EventData->GetTTTback_E_Mult();
    m_StripYEMult = m_EventData->GetTTTfront_E_Mult();
    m_StripXTMult = m_EventData->GetTTTback_T_Mult();
    m_StripYTMult = m_EventData->GetTTTfront_T_Mult();
    m_PadEMult    = m_EventData->GetPad_E_Mult();
    m_PadTMult    = m_EventData->GetPad_T_Mult();
    
    m_YY1RingEMult = m_EventData->GetYY1ring_E_Mult();
    m_YY1SectorEMult = m_EventData->GetYY1sector_E_Mult();
    m_YY1RingTMult = m_EventData->GetYY1ring_T_Mult();
    m_YY1SectorTMult = m_EventData->GetYY1sector_T_Mult();
    m_CsIEMult    = m_EventData->GetCsI_E_Mult();
    m_CsITMult    = m_EventData->GetCsI_T_Mult();
    // map<int, int> hitX;
    // map<int, int> hitY;

    // TTT YE
    for (unsigned int i = 0; i < m_StripYEMult; ++i) {
        if (m_EventData->GetTTTfront_Energy(i) > m_TTT_YE_RAW_Threshold
            && IsValidChannel(0, m_EventData->GetTTTfront_E_DetNbr(i),m_EventData->GetTTTfront_E_StripNbr(i))){
            double EY = fTTT_YE(m_EventData,i);
            if (EY > m_TTT_YE_CAL_Threshold)
                m_PreTreatedData->SetTTTfrontEnergy(m_EventData->GetTTTfront_E_DetNbr(i),m_EventData->GetTTTfront_E_StripNbr(i),EY);
        }
    }
    
    // TTT YT
    for (unsigned int i = 0; i < m_StripYTMult; ++i) {
        if (IsValidChannel(0, m_EventData->GetTTTfront_T_DetNbr(i),m_EventData->GetTTTfront_T_StripNbr(i)))
            m_PreTreatedData->SetTTTfrontTime(m_EventData->GetTTTfront_T_DetNbr(i),m_EventData->GetTTTfront_T_StripNbr(i),
                                              fTTT_YT(m_EventData,i));
    }

    // TTT XE
    for (unsigned int i = 0; i < m_StripXEMult; ++i) {
        // sf warning: Must2 uses here "<" but replaced by ">" for Tina
        if (m_EventData->GetTTTback_Energy(i) > m_TTT_XE_RAW_Threshold
            && IsValidChannel(1, m_EventData->GetTTTback_E_DetNbr(i),m_EventData->GetTTTback_E_StripNbr(i))){
            double EX = fTTT_XE(m_EventData,i);
            if (EX > m_TTT_XE_CAL_Threshold)
                m_PreTreatedData->SetTTTbackEnergy(m_EventData->GetTTTback_E_DetNbr(i),m_EventData->GetTTTback_E_StripNbr(i),EX);
        }
    }
    
    // TTT XT
    for (unsigned int i = 0; i < m_StripXTMult; ++i) {
        if (IsValidChannel(1, m_EventData->GetTTTback_T_DetNbr(i),m_EventData->GetTTTback_T_StripNbr(i)))
            m_PreTreatedData->SetTTTbackTime(m_EventData->GetTTTback_T_DetNbr(i),m_EventData->GetTTTback_T_StripNbr(i),
                                             fTTT_XT(m_EventData,i));
    }

    // Pad E
    for (unsigned int i = 0; i < m_PadEMult; ++i) {
        if (m_EventData->GetPad_Energy(i) > m_Pad_E_RAW_Threshold
            // sf warning: pad number to be implemented for Tina
            // && IsValidChannel(3, m_EventData->GetPad_E_DetNbr(i),m_EventData->GetMMCsIECristalNbr(i))) {
            && IsValidChannel(2, m_EventData->GetPad_E_DetNbr(i),0)){
            double EPad = fPad_E(m_EventData,i);
            if (EPad > m_Pad_E_CAL_Threshold)
                m_PreTreatedData->SetPadEnergy(m_EventData->GetPad_E_DetNbr(i),m_EventData->GetPad_E_PadNbr(i),EPad);
        }
    }
    
    // Pad T
    for (unsigned int i = 0; i < m_PadTMult; ++i) {
        // sf warning: pad number to be implemented for Tina
        // if (IsValidChannel(3, m_EventData->GetMMCsITDetectorNbr(i),m_EventData->GetMMCsITCristalNbr(i)))
        if (IsValidChannel(2, m_EventData->GetPad_T_DetNbr(i),0))
            m_PreTreatedData->SetPadTime(m_EventData->GetPad_T_DetNbr(i),m_EventData->GetPad_T_PadNbr(i),fPad_T(m_EventData,i));
    }

    // YY1 ring E
    for (unsigned int i = 0; i < m_YY1RingEMult; ++i) {
        // sf warning: Must2 uses here "<" but replaced by ">" for Tina
        if (m_EventData->GetYY1ring_Energy(i) > m_YY1_RE_RAW_Threshold
            && IsValidChannel(0, m_EventData->GetYY1ring_E_DetNbr(i),m_EventData->GetYY1ring_E_StripNbr(i))){
            double ER = fYY1_RE(m_EventData,i);
            if (ER > m_YY1_RE_CAL_Threshold)
                m_PreTreatedData->SetYY1ringEnergy(m_EventData->GetYY1ring_E_DetNbr(i),m_EventData->GetYY1ring_E_StripNbr(i),ER);
        }
    }
    
    // YY1 ring T
    for (unsigned int i = 0; i < m_YY1RingTMult; ++i) {
        if (IsValidChannel(0, m_EventData->GetYY1ring_T_DetNbr(i),m_EventData->GetYY1ring_T_StripNbr(i)))
            m_PreTreatedData->SetYY1ringTime(m_EventData->GetYY1ring_T_DetNbr(i),m_EventData->GetYY1ring_T_StripNbr(i),
                                             fYY1_RT(m_EventData,i));
    }

    // YY1 sector E
    for (unsigned int i = 0; i < m_YY1SectorEMult; ++i) {
        // sf warning: Must2 uses here "<" but replaced by ">" for Tina
        if (m_EventData->GetYY1sector_Energy(i) > m_YY1_SE_RAW_Threshold
            && IsValidChannel(1, m_EventData->GetYY1sector_E_DetNbr(i),m_EventData->GetYY1sector_E_StripNbr(i))){
            double ES = fYY1_SE(m_EventData,i);
            if (ES > m_YY1_SE_CAL_Threshold)
                m_PreTreatedData->SetYY1sectorEnergy(m_EventData->GetYY1sector_E_DetNbr(i),m_EventData->GetYY1sector_E_StripNbr(i),ES);
        }
    }
    
    // YY1 sector T
    for (unsigned int i = 0; i < m_YY1SectorTMult; ++i) {
        if (IsValidChannel(1, m_EventData->GetYY1sector_T_DetNbr(i),m_EventData->GetYY1sector_T_StripNbr(i)))
            m_PreTreatedData->SetYY1sectorTime(m_EventData->GetYY1sector_T_DetNbr(i),m_EventData->GetYY1sector_T_StripNbr(i),
                                               fYY1_ST(m_EventData,i));
    }

    // CsI E
    for (unsigned int i = 0; i < m_CsIEMult; ++i) {
        if (m_EventData->GetCsI_Energy(i) > m_CsI_E_RAW_Threshold
            // sf warning: pad number to be implemented for Tina
            // && IsValidChannel(3, m_EventData->GetPad_E_DetNbr(i),m_EventData->GetMMCsIECristalNbr(i))) {
            && IsValidChannel(2, m_EventData->GetCsI_E_DetNbr(i),0)){
            double ECsI = fCsI_E(m_EventData,i);
            if (ECsI > m_CsI_E_CAL_Threshold)
                m_PreTreatedData->SetCsIEnergy(m_EventData->GetCsI_E_DetNbr(i),ECsI);
        }
    }
    
    // CsI T
    for (unsigned int i = 0; i < m_CsITMult; ++i) {
        // sf warning: pad number to be implemented for Tina
        // if (IsValidChannel(3, m_EventData->GetMMCsITDetectorNbr(i),m_EventData->GetMMCsITCristalNbr(i)))
        if (IsValidChannel(2, m_EventData->GetCsI_T_DetNbr(i),0))
            m_PreTreatedData->SetCsITime(m_EventData->GetCsI_T_DetNbr(i),fCsI_T(m_EventData,i));
    }//
    
    return;
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TTinaPhysics::Match_X_Y() {
    vector<TVector2> ArrayOfGoodCouple;
    ArrayOfGoodCouple.clear();
    
    // sf warning: check that for Tina we define indeed front=X and back=Y
    // bm: front=Y, back=X.
    m_StripXEMult = m_PreTreatedData->GetTTTback_E_Mult();
    m_StripYEMult = m_PreTreatedData->GetTTTfront_E_Mult();
    double matchSigma  = m_StripEnergyMatchingSigma;
    double NmatchSigma = m_StripEnergyMatchingNumberOfSigma;
    
    // remove very high multiplicity events that are not physical anyway
    if (m_StripXEMult > m_MaximumStripMultiplicityAllowed
        || m_StripYEMult > m_MaximumStripMultiplicityAllowed) {
        return ArrayOfGoodCouple;
    }
    
    // get detector multiplicity
    for (unsigned short i = 0; i < m_StripXEMult; i++) {
        unsigned short N = m_PreTreatedData->GetTTTback_E_DetNbr(i);
        m_StripXMultDet[N] += 1;
    }
    for (unsigned short j = 0; j < m_StripYEMult; j++) {
        unsigned short N = m_PreTreatedData->GetTTTfront_E_DetNbr(j);
        m_StripYMultDet[N] += 1;
    }
    
    for (unsigned short i = 0; i < m_StripXEMult; i++) {
        for (unsigned short j = 0; j < m_StripYEMult; j++) {
            
            unsigned short StripXDetNbr = m_PreTreatedData->GetTTTback_E_DetNbr(i);
            unsigned short StripYDetNbr = m_PreTreatedData->GetTTTfront_E_DetNbr(j);
            // if same detector then check the energy
            if (StripXDetNbr == StripYDetNbr) {
                unsigned short DetNbr = StripXDetNbr;
                double StripXEnergy = m_PreTreatedData->GetTTTback_Energy(i);
                double StripXNbr    = m_PreTreatedData->GetTTTback_E_StripNbr(i);
                double StripYEnergy = m_PreTreatedData->GetTTTfront_Energy(j);
                double StripYNbr    = m_PreTreatedData->GetTTTfront_E_StripNbr(j);
                // look if energy matches
                // FIXME should be proportional to the energy loss in the DSSDs
                // if (abs(StripXEnergy - StripYEnergy) < 0.09 * (std::max(StripXEnergy, StripYEnergy))) {
                // negligible correction according to Adrien
                // sf warning: check that this also holds for Tina?!
                if (abs((StripXEnergy - StripYEnergy) / 2.) < NmatchSigma * matchSigma) {
                    // if the event is between two CsI crystals, it is rejected
                    // if (m_Ignore_not_matching_CsI) {
                    //   bool check_validity = false;
                    //   for (unsigned int hh = 0; hh < 16; ++hh) {
                    //     if (Match_Si_CsI(StripXNbr, StripYNbr, hh + 1)) {
                    //       check_validity = true;
                    //     }
                    //   }
                    //   if (check_validity) {
                    //     ArrayOfGoodCouple.push_back(TVector2(i, j));
                    //   }
                    // }
                    ArrayOfGoodCouple.push_back(TVector2(i,j));
                    m_NMatchDet[DetNbr] += 1;
                }
            }
        }
    }
    
    unsigned int couple_size = ArrayOfGoodCouple.size();
    for (unsigned int i = 0; i < couple_size; ++i) {
        int N = m_PreTreatedData->GetTTTback_E_DetNbr(ArrayOfGoodCouple[i].X());
        CheckEvent(N);
        if (m_OrderMatch == 2)
            m_match_type.push_back(2);
        else
            m_match_type.push_back(1);
    }
    
    return ArrayOfGoodCouple;
}

////////////////////////////////////////////////////////////////////////////
void TTinaPhysics::CheckEvent(int N) {
    if (m_NMatchDet[N] > m_StripXMultDet[N] || m_NMatchDet[N] > m_StripYMultDet[N]) {
        m_OrderMatch = 2;
    }
    else {
        m_OrderMatch = 1;
    }
}

////////////////////////////////////////////////////////////////////////////
bool TTinaPhysics::IsValidChannel(const int& DetectorType,
                                  const int& telescope, const int& channel) {
    // sf warning: this routine is currently disabled
    if (DetectorType == 0)
        //return *(m_XChannelStatus[telescope - 1].begin() + channel - 1);
        return true;
    else if (DetectorType == 1)
        //return *(m_YChannelStatus[telescope - 1].begin() + channel - 1);
        return true;
    if (DetectorType == 2)
        //return *(m_PadChannelStatus[telescope - 1].begin() + channel - 1);
        return true;
    else
        return false;
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::ReadAnalysisConfig() {
    // sf warning: routine remains to be done
    bool ReadingStatus = false;
    
    // open analysis configuration file
    string FileName = "./configs/ConfigTina.dat";
    ifstream AnalysisConfigFile;
    AnalysisConfigFile.open(FileName.c_str());
    
    if (!AnalysisConfigFile.is_open()) {
        cout << "No ConfigTina.dat file found: default parameters loaded for analysis " << FileName << endl;
        return;
    }
    cout << "Loading user parameters for analysis from file ConfigTina.dat" << endl;
    
    // save it in a TAsciiFile
    TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
    asciiConfig->AppendLine("%%% ConfigTina.dat %%%");
    asciiConfig->Append(FileName.c_str());
    asciiConfig->AppendLine("");
    // read analysis config file
    string LineBuffer,DataBuffer,whatToDo;
    while (!AnalysisConfigFile.eof()) {
        // pick-up next line
        getline(AnalysisConfigFile, LineBuffer);
        // search for "header"
        string name = "ConfigTina";
        if (LineBuffer.compare(0, name.length(), name) == 0)
            ReadingStatus = true;
        // loop on tokens and data
        while (ReadingStatus ) {
            whatToDo="";
            AnalysisConfigFile >> whatToDo;
            // search for comment symbol (%)
            if (whatToDo.compare(0, 1, "%") == 0) {
                AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
            }
            else if (whatToDo=="E_RAW_THRESHOLD") {
                AnalysisConfigFile >> DataBuffer;
                m_TTT_XE_RAW_Threshold = atof(DataBuffer.c_str());
                cout << whatToDo << " " << m_TTT_XE_RAW_Threshold << endl;
            }
            else if (whatToDo=="E_THRESHOLD") {
                AnalysisConfigFile >> DataBuffer;
                m_TTT_XE_CAL_Threshold = atof(DataBuffer.c_str());
                cout << whatToDo << " " << m_TTT_XE_CAL_Threshold << endl;
            }
            else {
                ReadingStatus = false;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::Clear() {
    
    EventMultiplicity = 0;
    EventType.clear();
    SquareTelescopeNumber.clear();
    TrapezTelescopeNumber.clear();
    TTT_E.clear();
    TTT_T.clear();
    TTT_X.clear();
    TTT_Y.clear();
    //TTT_EX.clear();
    //TTT_TX.clear();
    //TTT_EY.clear();
    //TTT_TY.clear();
    //TelescopeNumber_X.clear();
    //TelescopeNumber_Y.clear();
    Pad_E.clear();
    Pad_T.clear();
    Pad_N.clear();
    YY1_E.clear();
    YY1_T.clear();
    YY1_R.clear();
    YY1_S.clear();
    //TTT_EX.clear();
    //TTT_TX.clear();
    //TTT_EY.clear();
    //TTT_TY.clear();
    //TelescopeNumber_X.clear();
    //TelescopeNumber_Y.clear();
    CsI_E.clear();
    CsI_T.clear();
    CsI_N.clear();
    Energy.clear();
    Time.clear();
    
    m_OrderMatch = 0;
    m_match_type.clear();
    m_NMatchDet.clear();
    m_StripXMultDet.clear();
    m_StripYMultDet.clear();
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::ReadConfiguration(NPL::InputParser parser) {
    // sf warning: routine remains to be done
    vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("Tina");
    if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << "//// " << blocks.size() << " detectors found " << endl;
    
    vector<string> cart = {"POS","Shape"};
    vector<string> sphe = {"R","Theta","Phi","T","Shape"};
    
    for(unsigned int i = 0 ; i < blocks.size() ; i++){
        if(blocks[i]->HasTokenList(cart)){
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  Tina " << i+1 <<  endl;
            TVector3 Pos = blocks[i]->GetTVector3("POS","mm");
            string Shape = blocks[i]->GetString("Shape");
            double Translation = 0;
            AddDetector(Pos,Shape,Translation);
        }
        else if(blocks[i]->HasTokenList(sphe)){
            if(NPOptionManager::getInstance()->GetVerboseLevel())
                cout << endl << "////  Tina " << i+1 <<  endl;
            double R = blocks[i]->GetDouble("R","mm");
            double Theta = blocks[i]->GetDouble("Theta","deg");
            double Phi = blocks[i]->GetDouble("Phi","deg");
            double T = blocks[i]->GetDouble("T","mm");
            string Shape = blocks[i]->GetString("Shape");
            // sf warning: need to implement AddDetector(...T...) for YY1
            AddDetector(R,Theta,Phi,T,Shape);
        }
        else{
            cout << "Error: check your input file formatting" << endl;
            exit(1);
        }
    }
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::InitSpectra() {
    m_Spectra = new TTinaSpectra(m_NumberOfTelescopes);
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::FillSpectra() {
    m_Spectra -> FillRawSpectra(m_EventData);
    m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
    m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::CheckSpectra() {
    m_Spectra->CheckSpectra();
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::ClearSpectra() {
}

///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TTinaPhysics::GetSpectra() {
    if(m_Spectra)
        return m_Spectra->GetMapHisto();
    else{
        map< string,TH1*> empty;
        return empty;
    }
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::WriteSpectra() {
    m_Spectra->WriteSpectra();
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::AddParameterToCalibrationManager() {
    // sf warning: routine remains to be done
    CalibrationManager* Cal = CalibrationManager::getInstance();
    for (int i = 0; i < m_NumberOfTelescopes; ++i) {
        Cal->AddParameter("Tina", "D"+ NPL::itoa(i+1)+"_ENERGY","Tina_D"+ NPL::itoa(i+1)+"_ENERGY");
        Cal->AddParameter("Tina", "D"+ NPL::itoa(i+1)+"_TIME","Tina_D"+ NPL::itoa(i+1)+"_TIME");
    }
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::InitializeRootInputRaw() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("Tina",true);
    inputChain->SetBranchAddress("Tina",&m_EventData);
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::InitializeRootInputPhysics() {
    TChain* inputChain = RootInput::getInstance()->GetChain();
    inputChain->SetBranchStatus("Tina", true);
    inputChain->SetBranchStatus("EventMultiplicity", true);
    inputChain->SetBranchStatus("EventType", true);
    inputChain->SetBranchStatus("SquareTelescopeNumber", true);
    inputChain->SetBranchStatus("TTT_E", true);
    inputChain->SetBranchStatus("TTT_T", true);
    inputChain->SetBranchStatus("TTT_X", true);
    inputChain->SetBranchStatus("TTT_Y", true);
    // sf warning: following commented since these vectors are not filled in BuildPhysicalEvent
    //inputChain->SetBranchStatus("TTT_EX", true);
    //inputChain->SetBranchStatus("TTT_TX", true);
    //inputChain->SetBranchStatus("TTT_EY", true);
    //inputChain->SetBranchStatus("TTT_TY", true);
    //inputChain->SetBranchStatus("TelescopeNumber_X", true);
    //inputChain->SetBranchStatus("TelescopeNumber_Y", true);
    inputChain->SetBranchStatus("Pad_E", true);
    inputChain->SetBranchStatus("Pad_T", true);
    inputChain->SetBranchStatus("Pad_N", true);
    
    inputChain->SetBranchStatus("TrapezTelescopeNumber", true);
    inputChain->SetBranchStatus("YY1_E", true);
    inputChain->SetBranchStatus("YY1_T", true);
    inputChain->SetBranchStatus("YY1_R", true);
    inputChain->SetBranchStatus("YY1_S", true);
    inputChain->SetBranchStatus("CsI_E", true);
    inputChain->SetBranchStatus("CsI_T", true);
    inputChain->SetBranchStatus("CsI_N", true);
    inputChain->SetBranchStatus("TotalEnergy", true);
    inputChain->SetBranchAddress("Tina", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
void TTinaPhysics::InitializeRootOutput() {
    TTree* outputTree = RootOutput::getInstance()->GetTree();
    outputTree->Branch("Tina", "TTinaPhysics", &m_EventPhysics);
}

///////////////////////////////////////////////////////////////////////////
// specific to Tina         ///////////////////////////////////////////////
TVector3 TTinaPhysics::GetPositionOfInteraction(const int i) const {
    // cout << "SquareTelescopeNumber[i]=" << SquareTelescopeNumber[i] << " TTT_X[i]=" << TTT_X[i] << " TTT_Y[i]=" << TTT_Y[i] << endl;
    TVector3 Position
    = TVector3(GetStripPositionX(SquareTelescopeNumber[i], TTT_X[i], TTT_Y[i]),
               GetStripPositionY(SquareTelescopeNumber[i], TTT_X[i], TTT_Y[i]),
               GetStripPositionZ(SquareTelescopeNumber[i], TTT_X[i], TTT_Y[i]));
    
    return (Position);
}

///////////////////////////////////////////////////////////////////////////
// specific to Tina         ///////////////////////////////////////////////
TVector3 TTinaPhysics::GetRingPositionOfInteraction(const int i) const {
    TVector3 RingPosition
    = TVector3(GetRingPositionX(TrapezTelescopeNumber[i], YY1_R[i]),
               GetRingPositionY(TrapezTelescopeNumber[i], YY1_R[i]),
               GetRingPositionZ(TrapezTelescopeNumber[i], YY1_R[i]));
    return (RingPosition);
}

///////////////////////////////////////////////////////////////////////////
// specific to Tina         ///////////////////////////////////////////////
TVector3 TTinaPhysics::GetTelescopeNormal(const int i) const {
    
    TVector3 U = TVector3(GetStripPositionX(SquareTelescopeNumber[i], 128, 1),
                          GetStripPositionY(SquareTelescopeNumber[i], 128, 1),
                          GetStripPositionZ(SquareTelescopeNumber[i], 128, 1))
    - TVector3(GetStripPositionX(SquareTelescopeNumber[i], 1, 1),
               GetStripPositionY(SquareTelescopeNumber[i], 1, 1),
               GetStripPositionZ(SquareTelescopeNumber[i], 1, 1));
    TVector3 V = TVector3(GetStripPositionX(SquareTelescopeNumber[i], 128, 128),
                          GetStripPositionY(SquareTelescopeNumber[i], 128, 128),
                          GetStripPositionZ(SquareTelescopeNumber[i], 128, 128))
    - TVector3(GetStripPositionX(SquareTelescopeNumber[i], 128, 1),
               GetStripPositionY(SquareTelescopeNumber[i], 128, 1),
               GetStripPositionZ(SquareTelescopeNumber[i], 128, 1));
    TVector3 Normal = U.Cross(V);
    
    return (Normal.Unit());
}

///////////////////////////////////////////////////////////////////////////
// specific to Tina         ///////////////////////////////////////////////
TVector3 TTinaPhysics::GetRingNormal(const int i) const {
    
    TVector3 P = TVector3(GetRingPositionX(TrapezTelescopeNumber[i], 16),
                          GetRingPositionY(TrapezTelescopeNumber[i], 16),
                          GetRingPositionZ(TrapezTelescopeNumber[i], 16))
    - TVector3(GetRingPositionX(TrapezTelescopeNumber[i], 1),
               GetRingPositionY(TrapezTelescopeNumber[i], 1),
               GetRingPositionZ(TrapezTelescopeNumber[i], 1));
    
    TVector3 Q = TVector3(0,0,1).Cross(P);
    TVector3 Normal = P.Cross(Q);
    
    return (Normal.Unit());
}

///////////////////////////////////////////////////////////////////////////
namespace TINA_LOCAL {

double fTTT_XE(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetTTTback_E_DetNbr(i));
    name += "_TTT_X";
    name += NPL::itoa(m_EventData->GetTTTback_E_StripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetTTTback_Energy(i));
}

double fTTT_XT(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetTTTback_T_DetNbr(i));
    name += "_TTT_X";
    name += NPL::itoa(m_EventData->GetTTTback_T_StripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetTTTback_Time(i));
}

double fTTT_YE(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetTTTfront_E_DetNbr(i));
    name += "_TTT_Y";
    name += NPL::itoa(m_EventData->GetTTTfront_E_StripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetTTTfront_Energy(i));
}

double fTTT_YT(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetTTTfront_T_DetNbr(i));
    name += "_TTT_Y";
    name += NPL::itoa(m_EventData->GetTTTfront_T_StripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetTTTfront_Time(i));
}

double fPad_E(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetPad_E_DetNbr(i));
    name += "_Pad";
    //name += NPL::itoa(m_EventData->GetMMCsIECristalNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetPad_Energy(i));
}

double fPad_T(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetPad_T_DetNbr(i));
    name += "_Pad";
    //name += NPL::itoa(m_EventData->GetMMCsITCristalNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetPad_Time(i));
}
double fYY1_RE(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetYY1ring_E_DetNbr(i));
    name += "_YY1_R";
    name += NPL::itoa(m_EventData->GetYY1ring_E_StripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetYY1ring_Energy(i));
}

double fYY1_RT(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetYY1ring_T_DetNbr(i));
    name += "_YY1_R";
    name += NPL::itoa(m_EventData->GetYY1ring_T_StripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetYY1ring_Time(i));
}

double fYY1_SE(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetYY1sector_E_DetNbr(i));
    name += "_YY1_S";
    name += NPL::itoa(m_EventData->GetYY1sector_E_StripNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetYY1sector_Energy(i));
}

double fYY1_ST(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetYY1sector_T_DetNbr(i));
    name += "_YY1_S";
    name += NPL::itoa(m_EventData->GetYY1sector_T_StripNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetYY1sector_Time(i));
}

double fCsI_E(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetCsI_E_DetNbr(i));
    name += "_CsI";
    //name += NPL::itoa(m_EventData->GetMMCsIECristalNbr(i));
    name += "_E";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetCsI_Energy(i));
}

double fCsI_T(const TTinaData* m_EventData, const int& i) {
    static string name;
    name = "TINA/T";
    name += NPL::itoa(m_EventData->GetCsI_T_DetNbr(i));
    name += "_CsI";
    //name += NPL::itoa(m_EventData->GetMMCsITCristalNbr(i));
    name += "_T";
    return CalibrationManager::getInstance()->ApplyCalibration(name, m_EventData->GetCsI_Time(i));
}
} // namespace MUST2_LOCAL

////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TTinaPhysics::Construct() {
    // construct method to be passed to the DetectorFactory
    return (NPL::VDetector*) new TTinaPhysics();
}

////////////////////////////////////////////////////////////////////////////////
extern "C"{
class proxy_Tina{
    // register the construct method to the factory
public:
    proxy_Tina(){
        NPL::DetectorFactory::getInstance()->AddToken("Tina","Tina");
        NPL::DetectorFactory::getInstance()->AddDetector("Tina",TTinaPhysics::Construct);
    }
};

proxy_Tina p_Tina;
}
