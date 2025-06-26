/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_SI Treated  data                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include "TTOGAXSI_SIPhysics.h"

//   STL
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <limits>
using namespace std;

//   NPL
#include "RootInput.h"
#include "RootOutput.h"
#include "NPDetectorFactory.h"
#include "NPOptionManager.h"

#include "Math/Transform3D.h"
#include "Math/RotationZYX.h"
#include "NPSystemOfUnits.h"
using namespace NPUNITS;

//   ROOT
#include "TChain.h"

ClassImp(TTOGAXSI_SIPhysics)


///////////////////////////////////////////////////////////////////////////
/*
TTOGAXSI_SIPhysics::TTOGAXSI_SIPhysics()
   : m_EventData(new TTOGAXSI_SIData),
     m_PreTreatedData(new TTOGAXSI_SIData),
     m_EventPhysics(this),
     m_Spectra(0),
     m_E_RAW_Threshold(0), // adc channels
     m_E_Threshold(0),     // MeV
     m_NumberOfDetectors(0) {
}
*/
TTOGAXSI_SIPhysics::TTOGAXSI_SIPhysics(){
  EventMultiplicity = 0;
  EventMultiplicity_cluster = 0;
  m_EventData = new TTOGAXSI_SIData;
  m_PreTreatedData = new TTOGAXSI_SIData;
  m_EventPhysics = this;
  m_Spectra = NULL;
  m_E_RAW_Threshold = 0; // adc channels
  m_E_Threshold = 0;

  m_NumberOfInnerXDetectors = 0;
  m_NumberOfInnerZDetectors = 0;
  m_NumberOfOuterXDetectors = 0;
  m_NumberOfOuterZDetectors = 0;
  m_NumberOfClusterInnerDetectors = 0;
  m_NumberOfClusterX1Detectors = 0;
  m_NumberOfClusterY1Detectors = 0;
  m_NumberOfClusterX2Detectors = 0;
  m_NumberOfClusterY2Detectors = 0;

  m_MaximumStripMultiplicityAllowed = 10;
  m_StripEnergyMatching = 0.100;
  ///////////////////
  //Inner Detector //
  ///////////////////
  //WaferX parameter
  InnerX_Wafer_Length=78.4*mm;
  InnerX_Wafer_Width=51*mm;
  InnerX_Wafer_Thickness=0.1*mm;

  // PCB parameter
  InnerX_PCB_Thickness=3*mm;
  InnerX_Wafer_LongitudinalStrips = 255;
//  InnerX_Wafer_LongitudinalStrips = 64;
  InnerX_Wafer_TransverseStrips = 1;

  InnerX_PCB_Length = 84.9*mm;
  InnerX_PCB_Width = 55.83*mm;

  //WaferZ parameter
  InnerZ_Wafer_Length=78.4*mm;
  InnerZ_Wafer_Width=51*mm;
  InnerZ_Wafer_Thickness=0.1*mm;

  // PCB parameter
  InnerZ_PCB_Thickness=3*mm;
  InnerZ_Wafer_LongitudinalStrips = 1;
  InnerZ_Wafer_TransverseStrips = 392;
//  InnerZ_Wafer_TransverseStrips = 98;

  InnerZ_PCB_Length = 84.9*mm;
  InnerZ_PCB_Width = 55.83*mm;

  ///////////////////
  //Outer Detector //
  ///////////////////
  //WaferX parameter
  OuterX_Wafer_Length=78.4*mm;
  OuterX_Wafer_Width=51*mm;
  OuterX_Wafer_Thickness=0.1*mm;

  // PCB parameter
  OuterX_PCB_Thickness=3*mm;
  OuterX_Wafer_LongitudinalStrips = 392;
//  OuterX_Wafer_LongitudinalStrips = 98;
  OuterX_Wafer_TransverseStrips = 1;

  OuterX_PCB_Length = 107.16*mm;
  OuterX_PCB_Width = 84.9*mm;
  OuterX_PCB_gap = 3.5*mm;

  //WaferZ parameter
  OuterZ_Wafer_Length=78.4*mm;
  OuterZ_Wafer_Width=51*mm;
  OuterZ_Wafer_Thickness=0.1*mm;

  // PCB parameter
  OuterZ_PCB_Thickness=3*mm;
  OuterZ_Wafer_LongitudinalStrips = 1;
  OuterZ_Wafer_TransverseStrips = 255;
//  OuterZ_Wafer_TransverseStrips = 64;

  OuterZ_PCB_Length = 108.6*mm;
  OuterZ_PCB_Width = 83.3*mm;
  OuterZ_PCB_gap = 3.5*mm;

  /////////////////////
  //Cluster Detector //
  /////////////////////
  //Wafer parameter
  ClusterInner_Wafer_Base = 56*mm;
  ClusterInner_Wafer_Top = 7*mm;
  ClusterInner_Wafer_Height = 80*mm;
  ClusterInner_Wafer_Thickness = 0.1*mm;
 
  ClusterInner_ActiveWafer_Base = 88.*mm;
  ClusterInner_ActiveWafer_Top = 12.*mm;
  ClusterInner_ActiveWafer_Height = 66.*mm;
  ClusterInner_ActiveWafer_Thickness = 100*mm;
  
  ClusterInner_Wafer_TransverseStrips = 128;
  ClusterInner_Wafer_LongitudinalStrips = 128;


  //Cluster X Parameter
  ClusterX_Wafer_Length=78.4*mm;
  ClusterX_Wafer_Width=51*mm;
  ClusterX_Wafer_Thickness=0.1*mm;

  // PCB parameter
  ClusterX_PCB_Thickness=3*mm;
  ClusterX_Wafer_LongitudinalStrips = 392;
//  ClusterX_Wafer_LongitudinalStrips = 98;
  ClusterX_Wafer_TransverseStrips = 1;

  ClusterX_PCB_Length = 107.16*mm;
  ClusterX_PCB_Width = 84.9*mm;
  ClusterX_PCB_gap = 3.5*mm;


  //Cluster Y Parameter
  ClusterY_Wafer_Length=78.4*mm;
  ClusterY_Wafer_Width=51*mm;
  ClusterY_Wafer_Thickness=0.1*mm;

  // PCB parameter
  ClusterY_PCB_Thickness=3*mm;
  ClusterY_Wafer_LongitudinalStrips = 1;
  ClusterY_Wafer_TransverseStrips = 255;
//  ClusterY_Wafer_TransverseStrips = 64;

  ClusterY_PCB_Length = 108.6*mm;
  ClusterY_PCB_Width = 83.3*mm;
  ClusterY_PCB_gap = 3.5*mm;
}

///////////////////////////////////////////////////////////////////////////
/// A usefull method to bundle all operation to add a detector
/* void TTOGAXSI_SIPhysics::AddDetector(TVector3 , string ){
  // In That simple case nothing is done
  // Typically for more complex detector one would calculate the relevant 
  // positions (stripped silicon) or angles (gamma array)
  m_NumberOfDetectors++;
} 
*/

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddInnerXDetector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfInnerXDetectors++;
  cout << m_NumberOfInnerXDetectors << endl;

  double LongitudinalPitch = InnerX_Wafer_Width/InnerX_Wafer_LongitudinalStrips;

  // Vector C position of detector face center
  double Recess = 0.5*InnerX_Wafer_Thickness;
  TVector3 C(0,R-Recess,Z);
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position inner silicon detector
  //Moving to corner of the silicon
     TVector3 P_1_1 = C
    -V*0.5*InnerX_Wafer_Width
    +V*0.5*LongitudinalPitch; // middle of strip
  
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<InnerX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    
    P = P_1_1 + Ref + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());
  
//    cout << i << ": " << lineX[0] << endl;

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }
  
  m_InnerXStripPositionX.push_back(OneDetectorStripPositionX);  
  m_InnerXStripPositionY.push_back(OneDetectorStripPositionY);  
  m_InnerXStripPositionZ.push_back(OneDetectorStripPositionZ);  

} 


///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddInnerZDetector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfInnerZDetectors++;

  double TransversePitch = InnerZ_Wafer_Length/InnerZ_Wafer_TransverseStrips;

  // Vector C position of detector face center
  double Recess = 0.5*InnerZ_Wafer_Thickness;
  TVector3 C(0,R-Recess,Z);
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position inner silicon detector
  //Moving to corner of the silicon
/*
      TVector3 P_1_1 = C
    -U*0.5*InnerZ_Wafer_Length
    +U*0.5*TransversePitch // middle of strip
    -V*0.5*InnerZ_Wafer_Width
    +V*0.5*LongitudinalPitch; // middle of strip
*/
      TVector3 P_1_1 = C
    -U*0.5*InnerZ_Wafer_Length
    +U*0.5*TransversePitch; // middle of strip


  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<InnerZ_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
  
    P = P_1_1 + Ref + i*U*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }
  
  m_InnerZStripPositionX.push_back(OneDetectorStripPositionX);  
  m_InnerZStripPositionY.push_back(OneDetectorStripPositionY);  
  m_InnerZStripPositionZ.push_back(OneDetectorStripPositionZ);  
  
} 


///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddOuterXDetector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfOuterXDetectors++;

  double LongitudinalPitch = OuterX_Wafer_Length/OuterX_Wafer_LongitudinalStrips;

  // Vector C position of detector face center
  double Recess = 0.5*OuterX_Wafer_Thickness;
  TVector3 C(0,R-Recess,Z);
  C.RotateZ(-Phi);


  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position 1st outer silicon detector:
  //Moving to corner of the silicon
  TVector3 P_1_1 = C
    -V*0.5*OuterX_Wafer_Length
    +V*0.5*LongitudinalPitch; // middle of strip

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<OuterX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    P = P_1_1 + Ref + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

    cout << i << ": " << lineX[0] << endl;

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_OuterXStripPositionX.push_back(OneDetectorStripPositionX);
  m_OuterXStripPositionY.push_back(OneDetectorStripPositionY);
  m_OuterXStripPositionZ.push_back(OneDetectorStripPositionZ);

  // Position 2nd outer silicon outer detector
  // Moving to corner of the silicon
  P_1_1 = C
    -V*0.5*OuterX_Wafer_Length
    +V*0.5*LongitudinalPitch; // middle of strip

  for(int i=0; i<OuterX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    P = P_1_1 + Ref + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

    cout << i << ": " << lineX[0] << endl;

    m_OuterXStripPositionX[m_NumberOfOuterXDetectors-1].push_back(lineX);
    m_OuterXStripPositionY[m_NumberOfOuterXDetectors-1].push_back(lineY);
    m_OuterXStripPositionZ[m_NumberOfOuterXDetectors-1].push_back(lineZ);

  }

} 

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddOuterZDetector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfOuterZDetectors++;

  double TransversePitch = OuterZ_Wafer_Width/OuterZ_Wafer_TransverseStrips;

  // Vector C position of detector face center
  double Recess = 0.5*OuterZ_Wafer_Thickness;
  TVector3 C(0,R-Recess,Z);
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position 1st outer silicon detector:
  //Moving to corner of the silicon
  TVector3 P_1_1 = C
    -U*0.5*OuterZ_PCB_gap //In between wafer
    -U*OuterZ_Wafer_Width // External wafer edge
    +U*0.5*TransversePitch; // middle of strip

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<OuterZ_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    P = P_1_1 + Ref + i*U*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());


    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_OuterZStripPositionX.push_back(OneDetectorStripPositionX);
  m_OuterZStripPositionY.push_back(OneDetectorStripPositionY);
  m_OuterZStripPositionZ.push_back(OneDetectorStripPositionZ);

  // Position 2nd outer silicon outer detector
  // Moving to corner of the silicon
  P_1_1 = C
    +U*0.5*OuterZ_PCB_gap //In between wafer
    +U*0.5*TransversePitch; // middle of strip

  for(int i=0; i<OuterZ_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    P = P_1_1 + Ref + i*U*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());


    m_OuterZStripPositionX[m_NumberOfOuterZDetectors-1].push_back(lineX);
    m_OuterZStripPositionY[m_NumberOfOuterZDetectors-1].push_back(lineY);
    m_OuterZStripPositionZ[m_NumberOfOuterZDetectors-1].push_back(lineZ);

  }

} 


///////////////////////////////////////////////////////////////////////////
/*
void TTOGAXSI_SIPhysics::AddClusterDetector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfClusterDetectors++;

  double LongitudinalPitch = Cluster_ActiveWafer_Height/Cluster_ActiveWafer_LongitudinalStrips;
  double TransversePitch = Cluster_ActiveWafer_Base/Cluster_ActiveWafer_TransverseStrips;

  // Vector C position of detector face center
  //double Recess = (Inner_PCB_Thickness-Inner_PCB_Step-Inner_Wafer_Thickness);
  double Recess = (Cluster_ActiveWafer_Thickness);
  TVector3 C(0,R+Recess,Z);
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position inner silicon detector
  //Moving to corner of the silicon
     TVector3 P_1_1 = C
    -U*0.5*Cluster_ActiveWafer_Height
    +U*0.5*TransversePitch // middle of strip
    -V*0.5*Cluster_ActiveWafer_Base
    +V*0.5*LongitudinalPitch; // middle of strip
  
  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<Cluster_ActiveWafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    for(int j=0; j<Cluster_ActiveWafer_LongitudinalStrips; j++){
      P = P_1_1 + Ref + i*U*TransversePitch + j*V*LongitudinalPitch;
      lineX.push_back(P.X());
      lineY.push_back(P.Y());
      lineZ.push_back(P.Z());
    }
    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }
  
  m_ClusterStripPositionX.push_back(OneDetectorStripPositionX);  
  m_ClusterStripPositionY.push_back(OneDetectorStripPositionY);  
  m_ClusterStripPositionZ.push_back(OneDetectorStripPositionZ);  
  
} 
*/

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddClusterX1Detector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfClusterX1Detectors++;

  double LongitudinalPitch = ClusterX_Wafer_Length/ClusterX_Wafer_LongitudinalStrips;


  // Vector C position of detector face center
  double Recess = 0.5*ClusterX_Wafer_Thickness;
  TVector3 C(R,0,Z - Recess);
//  TVector3 C(0,R,Z);

  C.RotateZ(-Phi);
  cout << "C vector: " << C.x() << endl;

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  cout << "W vector: " << W.y() << endl;


  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);
  cout << "V vector: " << V.x() << endl;

  //Position 1st outer silicon detector:
  //Moving to corner of the silicon
  TVector3 P_1_1 = C
    -V*0.5*ClusterX_Wafer_Length
    +V*0.5*LongitudinalPitch; // middle of strip


  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<ClusterX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    P = P_1_1 + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

    cout << i << ": " << lineX[0] << endl;

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_ClusterX1StripPositionX.push_back(OneDetectorStripPositionX);
  m_ClusterX1StripPositionY.push_back(OneDetectorStripPositionY);
  m_ClusterX1StripPositionZ.push_back(OneDetectorStripPositionZ);

  // Position 2nd outer silicon outer detector
  // Moving to corner of the silicon
  P_1_1 = C
    -V*0.5*ClusterX_Wafer_Length
    +V*0.5*LongitudinalPitch; // middle of strip

  for(int i=0; i<ClusterX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    P = P_1_1 + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

//    cout << i << ": " << lineX[0] << endl;

    m_ClusterX1StripPositionX[m_NumberOfClusterX1Detectors-1].push_back(lineX);
    m_ClusterX1StripPositionY[m_NumberOfClusterX1Detectors-1].push_back(lineY);
    m_ClusterX1StripPositionZ[m_NumberOfClusterX1Detectors-1].push_back(lineZ);

  }

} 

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddClusterY1Detector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfClusterY1Detectors++;

  double TransversePitch = ClusterY_Wafer_Width/ClusterY_Wafer_TransverseStrips;

  // Vector C position of detector face center
  double Recess = 0.5*ClusterY_Wafer_Thickness;
  TVector3 C(R,0,Z-Recess);

  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position 1st outer silicon detector:
  //Moving to corner of the silicon
  TVector3 P_1_1 = C
    -W*0.5*ClusterY_PCB_gap //In between wafer
    -W*ClusterY_Wafer_Width
    +W*0.5*TransversePitch; // middle of strip

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<ClusterY_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    P = P_1_1 + i*W*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

//    cout << i << ": " << lineY[0] << endl;

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_ClusterY1StripPositionX.push_back(OneDetectorStripPositionX);
  m_ClusterY1StripPositionY.push_back(OneDetectorStripPositionY);
  m_ClusterY1StripPositionZ.push_back(OneDetectorStripPositionZ);

  // Position 2nd outer silicon outer detector
  // Moving to corner of the silicon
  P_1_1 = C
    +W*0.5*ClusterY_PCB_gap //In between wafer
    +W*0.5*TransversePitch; // middle of strip

  for(int i=0; i<ClusterY_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    P = P_1_1 + i*W*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

//    cout << i << ": " << lineY[0] << endl;

    m_ClusterY1StripPositionX[m_NumberOfClusterY1Detectors-1].push_back(lineX);
    m_ClusterY1StripPositionY[m_NumberOfClusterY1Detectors-1].push_back(lineY);
    m_ClusterY1StripPositionZ[m_NumberOfClusterY1Detectors-1].push_back(lineZ);

  }
} 



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddClusterX2Detector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfClusterX2Detectors++;

  double LongitudinalPitch = ClusterX_Wafer_Length/ClusterX_Wafer_LongitudinalStrips;

  // Vector C position of detector face center
  double Recess = 0.5*ClusterX_Wafer_Thickness;
  TVector3 C(R,0,Z-Recess);
//  TVector3 C(0,R,Z);
  C.RotateZ(-Phi);

  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);

  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position 1st outer silicon detector:
  //Moving to corner of the silicon
  TVector3 P_1_1 = C
    -V*0.5*ClusterX_Wafer_Length
    +V*0.5*LongitudinalPitch; // middle of strip


  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<ClusterX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    P = P_1_1 + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

//    cout << i << ": " << lineX[0] << endl;

    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_ClusterX2StripPositionX.push_back(OneDetectorStripPositionX);
  m_ClusterX2StripPositionY.push_back(OneDetectorStripPositionY);
  m_ClusterX2StripPositionZ.push_back(OneDetectorStripPositionZ);

  // Position 2nd outer silicon outer detector
  // Moving to corner of the silicon
  P_1_1 = C
    -V*0.5*ClusterX_Wafer_Length
    +V*0.5*LongitudinalPitch; // middle of strip

  for(int i=0; i<ClusterX_Wafer_LongitudinalStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    P = P_1_1 + i*V*LongitudinalPitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());

//    cout << i << ": " << lineX[0] << endl;

    m_ClusterX2StripPositionX[m_NumberOfClusterX2Detectors-1].push_back(lineX);
    m_ClusterX2StripPositionY[m_NumberOfClusterX2Detectors-1].push_back(lineY);
    m_ClusterX2StripPositionZ[m_NumberOfClusterX2Detectors-1].push_back(lineZ);

  }

} 

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddClusterY2Detector(double R, double Z, double Phi, TVector3 Ref){
  m_NumberOfClusterY2Detectors++;

  double TransversePitch = ClusterY_Wafer_Width/ClusterY_Wafer_TransverseStrips;

  // Vector C position of detector face center
  double Recess = 0.5*ClusterY_Wafer_Thickness;
  TVector3 C(R,0,Z-Recess);
  C.RotateZ(-Phi);


  // Vector W normal to detector face (pointing to the back)
  TVector3 W(0,1,0);
  W.RotateZ(-Phi);
  // Vector U on detector face (parallel to Z axis/longitudinal strips)
  TVector3 U = TVector3(0,0,1);
  // Vector V on detector face (parallel to transverse strips)
  TVector3 V = W.Cross(U);

  //Position 1st outer silicon detector:
  //Moving to corner of the silicon
  TVector3 P_1_1 = C
    -W*0.5*ClusterY_PCB_gap //In between wafer
    -W*ClusterY_Wafer_Width
    +W*0.5*TransversePitch; // middle of strip

  vector<double> lineX;
  vector<double> lineY;
  vector<double> lineZ;

  vector<vector<double>> OneDetectorStripPositionX;
  vector<vector<double>> OneDetectorStripPositionY;
  vector<vector<double>> OneDetectorStripPositionZ;

  TVector3 P;
  for(int i=0; i<ClusterY_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();
    P = P_1_1 + i*W*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());


    OneDetectorStripPositionX.push_back(lineX);
    OneDetectorStripPositionY.push_back(lineY);
    OneDetectorStripPositionZ.push_back(lineZ);
  }

  m_ClusterY2StripPositionX.push_back(OneDetectorStripPositionX);
  m_ClusterY2StripPositionY.push_back(OneDetectorStripPositionY);
  m_ClusterY2StripPositionZ.push_back(OneDetectorStripPositionZ);

  // Position 2nd outer silicon outer detector
  // Moving to corner of the silicon
  P_1_1 = C
    +W*0.5*ClusterY_PCB_gap //In between wafer
    +W*0.5*TransversePitch; // middle of strip

  for(int i=0; i<ClusterY_Wafer_TransverseStrips; i++){
    lineX.clear();
    lineY.clear();
    lineZ.clear();

    P = P_1_1 + i*W*TransversePitch;
    lineX.push_back(P.X());
    lineY.push_back(P.Y());
    lineZ.push_back(P.Z());


    m_ClusterY2StripPositionX[m_NumberOfClusterY2Detectors-1].push_back(lineX);
    m_ClusterY2StripPositionY[m_NumberOfClusterY2Detectors-1].push_back(lineY);
    m_ClusterY2StripPositionZ[m_NumberOfClusterY2Detectors-1].push_back(lineZ);

  }

} 




///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetInnerXPositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetInnerXStripPositionX(DetectorNumberInnerX[i], InnerXStrip[i], 1),
	GetInnerXStripPositionY(DetectorNumberInnerX[i], InnerXStrip[i], 1),
	GetInnerXStripPositionZ(DetectorNumberInnerX[i], InnerXStrip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetInnerZPositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetInnerZStripPositionX(DetectorNumberInnerZ[i], InnerZStrip[i], 1),
	GetInnerZStripPositionY(DetectorNumberInnerZ[i], InnerZStrip[i], 1),
	GetInnerZStripPositionZ(DetectorNumberInnerZ[i], InnerZStrip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetOuterXPositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetOuterXStripPositionX(DetectorNumberOuterX[i], OuterXStrip[i], 1),
	GetOuterXStripPositionY(DetectorNumberOuterX[i], OuterXStrip[i], 1),
	GetOuterXStripPositionZ(DetectorNumberOuterX[i], OuterXStrip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetOuterZPositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetOuterZStripPositionX(DetectorNumberOuterZ[i], OuterZStrip[i], 1),
	GetOuterZStripPositionY(DetectorNumberOuterZ[i], OuterZStrip[i], 1),
	GetOuterZStripPositionZ(DetectorNumberOuterZ[i], OuterZStrip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
/*
TVector3 TTOGAXSI_SIPhysics::GetClusterPositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetClusterStripPositionX(DetectorNumber[i], ClusterStripT[i], ClusterStripL[i]),
	GetClusterStripPositionY(DetectorNumber[i], ClusterStripT[i], ClusterStripL[i]),
	GetClusterStripPositionZ(DetectorNumber[i], ClusterStripT[i], ClusterStripL[i]));

  return Position;
}
*/

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetClusterX1PositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetClusterX1StripPositionX(DetectorNumberClusterX1[i], ClusterX1Strip[i], 1),
	GetClusterX1StripPositionY(DetectorNumberClusterX1[i], ClusterX1Strip[i], 1),
	GetClusterX1StripPositionZ(DetectorNumberClusterX1[i], ClusterX1Strip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetClusterY1PositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetClusterY1StripPositionX(DetectorNumberClusterY1[i], ClusterY1Strip[i], 1),
	GetClusterY1StripPositionY(DetectorNumberClusterY1[i], ClusterY1Strip[i], 1),
	GetClusterY1StripPositionZ(DetectorNumberClusterY1[i], ClusterY1Strip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetClusterX2PositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetClusterX2StripPositionX(DetectorNumberClusterX2[i], ClusterX2Strip[i], 1),
	GetClusterX2StripPositionY(DetectorNumberClusterX2[i], ClusterX2Strip[i], 1),
	GetClusterX2StripPositionZ(DetectorNumberClusterX2[i], ClusterX2Strip[i], 1));

  return Position;
}

///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetClusterY2PositionOfInteraction(const int i){
  TVector3 Position = TVector3(
	GetClusterY2StripPositionX(DetectorNumberClusterY2[i], ClusterY2Strip[i], 1),
	GetClusterY2StripPositionY(DetectorNumberClusterY2[i], ClusterY2Strip[i], 1),
	GetClusterY2StripPositionZ(DetectorNumberClusterY2[i], ClusterY2Strip[i], 1));

  return Position;
}


///////////////////////////////////////////////////////////////////////////
TVector3 TTOGAXSI_SIPhysics::GetDetectorNormal(const int i){
  return TVector3(0,0,0);
}


///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::BuildSimplePhysicalEvent() {
  BuildPhysicalEvent();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::BuildPhysicalEvent() {
  // apply thresholds and calibration
  PreTreat();


  vector<TVector2> inner = MatchInner();
  vector<TVector2> outer = MatchOuter();
  vector<TVector2> cluster1 = MatchCluster1();
  vector<TVector2> cluster2 = MatchCluster2();
 

//  cout << inner.size() << endl;
//  cout << outer.size() << endl; 

  for(unsigned int i=0; i<inner.size(); i++){
      int N = m_PreTreatedData->GetInnerX_E_DetectorNbr(inner[i].X());
      int innerX = m_PreTreatedData->GetInnerX_E_StripNbr(inner[i].X());
      int innerZ = m_PreTreatedData->GetInnerZ_E_StripNbr(inner[i].Y()); 

      double innerXE = m_PreTreatedData->GetInnerX_E_Energy(inner[i].X());
      double innerZE = m_PreTreatedData->GetInnerZ_E_Energy(inner[i].Y());
      // look for outer
      double outerXE = 0;
      double outerZE = 0;
      int outerX=0;
      int outerZ=0;
//	  cout << outer.size() << endl;     
      for(unsigned int j=0; j<outer.size(); j++){
        if(m_PreTreatedData->GetInnerX_E_DetectorNbr(outer[j].X())==N){
	  outerXE = m_PreTreatedData->GetOuterX_E_Energy(outer[j].X());
          outerZE = m_PreTreatedData->GetOuterZ_E_Energy(outer[j].Y());
          outerX = m_PreTreatedData->GetOuterX_E_StripNbr(outer[j].X());
          outerZ = m_PreTreatedData->GetOuterZ_E_StripNbr(outer[j].Y());
	}
      }	

      if(outerXE) {
	
	EventMultiplicity++;
	DetectorNumberInnerX.push_back(N);
	DetectorNumberInnerZ.push_back(N);
	DetectorNumberOuterX.push_back(N);
	DetectorNumberOuterZ.push_back(N);

	InnerXStrip.push_back(innerX);
	InnerXE.push_back(innerXE);
	InnerXPosX.push_back(GetInnerXPositionOfInteraction(EventMultiplicity-1).x());
	InnerXPosY.push_back(GetInnerXPositionOfInteraction(EventMultiplicity-1).y());
	InnerXPosZ.push_back(GetInnerXPositionOfInteraction(EventMultiplicity-1).z());

	InnerZStrip.push_back(innerZ);
	InnerZE.push_back(innerZE);
	InnerZPosX.push_back(GetInnerZPositionOfInteraction(EventMultiplicity-1).x());
	InnerZPosY.push_back(GetInnerZPositionOfInteraction(EventMultiplicity-1).y());
	InnerZPosZ.push_back(GetInnerZPositionOfInteraction(EventMultiplicity-1).z());


	OuterXStrip.push_back(outerX);
	OuterXE.push_back(outerXE);
	OuterXPosX.push_back(GetOuterXPositionOfInteraction(EventMultiplicity-1).x());
	OuterXPosY.push_back(GetOuterXPositionOfInteraction(EventMultiplicity-1).y());
	OuterXPosZ.push_back(GetOuterXPositionOfInteraction(EventMultiplicity-1).z());


	OuterZStrip.push_back(outerZ);
	OuterZE.push_back(outerZE);
	OuterZPosX.push_back(GetOuterZPositionOfInteraction(EventMultiplicity-1).x());
	OuterZPosY.push_back(GetOuterZPositionOfInteraction(EventMultiplicity-1).y());
	OuterZPosZ.push_back(GetOuterZPositionOfInteraction(EventMultiplicity-1).z());

      }

  }


  for(unsigned int i=0; i<cluster1.size(); i++){
      int N_cluster = m_PreTreatedData->GetClusterX1_E_DetectorNbr(cluster1[i].X());
      int clusterX1 = m_PreTreatedData->GetClusterX1_E_StripNbr(cluster1[i].X());
      int clusterY1 = m_PreTreatedData->GetClusterY1_E_StripNbr(cluster1[i].Y()); 

      double clusterX1E = m_PreTreatedData->GetClusterX1_E_Energy(cluster1[i].X());
      double clusterY1E = m_PreTreatedData->GetClusterY1_E_Energy(cluster1[i].Y());


      // look for outer
      double clusterX2E = 0;
      double clusterY2E = 0;
      int clusterX2=0;
      int clusterY2=0;

      for(unsigned int j=0; j<cluster2.size(); j++){
        if(m_PreTreatedData->GetClusterX1_E_DetectorNbr(cluster2[j].X())==N_cluster){
	  clusterX2E = m_PreTreatedData->GetClusterX2_E_Energy(cluster2[j].X());
          clusterY2E = m_PreTreatedData->GetClusterY2_E_Energy(cluster2[j].Y());
          clusterX2 = m_PreTreatedData->GetClusterX2_E_StripNbr(cluster2[j].X());
          clusterY2 = m_PreTreatedData->GetClusterY2_E_StripNbr(cluster2[j].Y());
	}
      }	

      if(clusterX2E) {
	EventMultiplicity_cluster++;
	DetectorNumberClusterX1.push_back(N_cluster);
	DetectorNumberClusterY1.push_back(N_cluster);
	DetectorNumberClusterX2.push_back(N_cluster);
	DetectorNumberClusterY2.push_back(N_cluster);

	
	ClusterX1Strip.push_back(clusterX1);
	ClusterX1E.push_back(clusterX1E);
	ClusterX1PosX.push_back(GetClusterX1PositionOfInteraction(EventMultiplicity_cluster-1).x());
	ClusterX1PosY.push_back(GetClusterX1PositionOfInteraction(EventMultiplicity_cluster-1).y());
	ClusterX1PosZ.push_back(GetClusterX1PositionOfInteraction(EventMultiplicity_cluster-1).z());

	ClusterY1Strip.push_back(clusterY1);
	ClusterY1E.push_back(clusterY1E);
	ClusterY1PosX.push_back(GetClusterY1PositionOfInteraction(EventMultiplicity_cluster-1).x());
	ClusterY1PosY.push_back(GetClusterY1PositionOfInteraction(EventMultiplicity_cluster-1).y());
	ClusterY1PosZ.push_back(GetClusterY1PositionOfInteraction(EventMultiplicity_cluster-1).z());

	ClusterX2Strip.push_back(clusterX2);
	ClusterX2E.push_back(clusterX2E);
	ClusterX2PosX.push_back(GetClusterX2PositionOfInteraction(EventMultiplicity_cluster-1).x());
	ClusterX2PosY.push_back(GetClusterX2PositionOfInteraction(EventMultiplicity_cluster-1).y());
	ClusterX2PosZ.push_back(GetClusterX2PositionOfInteraction(EventMultiplicity_cluster-1).z());

	ClusterY2Strip.push_back(clusterY2);
	ClusterY2E.push_back(clusterY2E);
	ClusterY2PosX.push_back(GetClusterY2PositionOfInteraction(EventMultiplicity_cluster-1).x());
	ClusterY2PosY.push_back(GetClusterY2PositionOfInteraction(EventMultiplicity_cluster-1).y());
	ClusterY2PosZ.push_back(GetClusterY2PositionOfInteraction(EventMultiplicity_cluster-1).z());

      }

  }


}


///////////////////////////////////////////////////////////////////////////
vector<TVector2> TTOGAXSI_SIPhysics::MatchInner(){
  vector<TVector2> ArrayOfGoodCouple;


  static unsigned int m_XEMult, m_ZEMult;
  m_XEMult = m_PreTreatedData->GetInnerXMultEnergy();
  m_ZEMult = m_PreTreatedData->GetInnerZMultEnergy();


  if(m_XEMult>m_MaximumStripMultiplicityAllowed || m_ZEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_XEMult; i++){
// 	cout << "_____X: " << i << endl;
    for(unsigned int j=0; j<m_ZEMult; j++){
      
      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetInnerX_E_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetInnerZ_E_DetectorNbr(j);

//      cout << "X Det Inner: " << XDetNbr << endl;
//      cout << "Z Det Inner: " << YDetNbr << endl;

//      cout << "X Energy Inner: " << m_PreTreatedData->GetInnerX_E_Energy(i) << endl;
//      cout << "Z Energy Inner : " << m_PreTreatedData->GetInnerZ_E_Energy(i) << endl;
//      cout << "Diff Energy Inner: " << abs(m_PreTreatedData->GetInnerX_E_Energy(i)-m_PreTreatedData->GetInnerZ_E_Energy(i))/2 << endl;

      // if same detector check energy
      if(XDetNbr == YDetNbr){
	// Declaration of variable for clarity
	double XE = m_PreTreatedData->GetInnerX_E_Energy(i);
	double ZE = m_PreTreatedData->GetInnerZ_E_Energy(j);

	// look if energy matches
	if(abs(XE-ZE)/2.<m_StripEnergyMatching){
	  ArrayOfGoodCouple.push_back(TVector2(i,j));
	}
      }
    }
  }

  return ArrayOfGoodCouple;
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TTOGAXSI_SIPhysics::MatchOuter(){
  vector<TVector2> ArrayOfGoodCouple;

  static unsigned int m_XEMult, m_ZEMult;
  m_XEMult = m_PreTreatedData->GetOuterXMultEnergy();
  m_ZEMult = m_PreTreatedData->GetOuterZMultEnergy();

  if(m_XEMult>m_MaximumStripMultiplicityAllowed || m_ZEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_XEMult; i++){
    for(unsigned int j=0; j<m_ZEMult; j++){
      
      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetOuterX_E_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetOuterZ_E_DetectorNbr(j);

      // if same detector check energy
      if(XDetNbr == YDetNbr){
	// Declaration of variable for clarity
	double XE = m_PreTreatedData->GetOuterX_E_Energy(i);
	double ZE = m_PreTreatedData->GetOuterZ_E_Energy(j);

	// look if energy matches
	if(abs(XE-ZE)/2.<m_StripEnergyMatching){
	  ArrayOfGoodCouple.push_back(TVector2(i,j));
	}
      }
    }
  }

  return ArrayOfGoodCouple;
}

///////////////////////////////////////////////////////////////////////////
vector<TVector2> TTOGAXSI_SIPhysics::MatchCluster1(){
  vector<TVector2> ArrayOfGoodCouple;

//  cout << "______Event start________" << endl;

  static unsigned int m_XEMult, m_YEMult;
  m_XEMult = m_PreTreatedData->GetClusterX1MultEnergy();
  m_YEMult = m_PreTreatedData->GetClusterY1MultEnergy();

//  cout << m_XEMult << endl;
//  cout << m_ZEMult << endl;

  if(m_XEMult>m_MaximumStripMultiplicityAllowed || m_YEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_XEMult; i++){
// 	cout << "_____X: " << i << endl;
    for(unsigned int j=0; j<m_YEMult; j++){
      
      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetClusterX1_E_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetClusterY1_E_DetectorNbr(j);

//      cout << "X Det: " << XDetNbr << endl;
//      cout << "Y Det: " << YDetNbr << endl;

//      cout << "X Energy: " << m_PreTreatedData->GetClusterX1_E_Energy(i) << endl;
//      cout << "Y Energy: " << m_PreTreatedData->GetClusterY1_E_Energy(i) << endl;
//      cout << "Diff Energy: " << abs(m_PreTreatedData->GetClusterX1_E_Energy(i)-m_PreTreatedData->GetClusterY1_E_Energy(i))/2 << endl;

      // if same detector check energy
      if(XDetNbr == YDetNbr){
	// Declaration of variable for clarity
	double XE = m_PreTreatedData->GetClusterX1_E_Energy(i);
	double YE = m_PreTreatedData->GetClusterY1_E_Energy(j);

	// look if energy matches
	if(abs(XE-YE)/2.<m_StripEnergyMatching){
	  ArrayOfGoodCouple.push_back(TVector2(i,j));
	}
      }
    }
  }

  return ArrayOfGoodCouple;
}


///////////////////////////////////////////////////////////////////////////
vector<TVector2> TTOGAXSI_SIPhysics::MatchCluster2(){
  vector<TVector2> ArrayOfGoodCouple;

//  cout << "______Event start________" << endl;

  static unsigned int m_XEMult, m_YEMult;
  m_XEMult = m_PreTreatedData->GetClusterX2MultEnergy();
  m_YEMult = m_PreTreatedData->GetClusterY2MultEnergy();

//  cout << m_XEMult << endl;
//  cout << m_ZEMult << endl;

  if(m_XEMult>m_MaximumStripMultiplicityAllowed || m_YEMult>m_MaximumStripMultiplicityAllowed){
    return ArrayOfGoodCouple;
  }

  for(unsigned int i=0; i<m_XEMult; i++){
// 	cout << "_____X: " << i << endl;
    for(unsigned int j=0; j<m_YEMult; j++){
      
      // Declaration of variable for clarity
      int XDetNbr = m_PreTreatedData->GetClusterX2_E_DetectorNbr(i);
      int YDetNbr = m_PreTreatedData->GetClusterY2_E_DetectorNbr(j);

//      cout << "X Det Cluster 2: " << XDetNbr << endl;
//      cout << "Y Det Cluster 2: " << YDetNbr << endl;

//      cout << "X Energy Cluster 2: " << m_PreTreatedData->GetClusterX2_E_Energy(i) << endl;
//      cout << "Y Energy Cluster 2: " << m_PreTreatedData->GetClusterY2_E_Energy(i) << endl;
//      cout << "Diff Energy Cluster 2: " << abs(m_PreTreatedData->GetClusterX2_E_Energy(i)-m_PreTreatedData->GetClusterY2_E_Energy(i))/2 << endl;

      // if same detector check energy
      if(XDetNbr == YDetNbr){
	// Declaration of variable for clarity
	double XE = m_PreTreatedData->GetClusterX2_E_Energy(i);
	double YE = m_PreTreatedData->GetClusterY2_E_Energy(j);

	// look if energy matches
	if(abs(XE-YE)/2.<m_StripEnergyMatching){
	  ArrayOfGoodCouple.push_back(TVector2(i,j));
	}
      }
    }
  }

  return ArrayOfGoodCouple;
}




///////////////////////////////////////////////////////////////////////////
int TTOGAXSI_SIPhysics::CheckEvent(){
  //Check the size of the different elements
  if(m_PreTreatedData->GetInnerXMultEnergy() == m_PreTreatedData->GetInnerZMultEnergy() )
    return 1;
  else
   return -1;
}


///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::PreTreat() {
  // This method typically applies thresholds and calibrations
  // Migh

  // clear pre-treated object
  ClearPreTreatedData();

  // instantiate CalibrationManager
  static CalibrationManager* Cal = CalibrationManager::getInstance();

  //////
  // InnerX Energy
  unsigned int sizeX = m_EventData->GetInnerXMultEnergy();
  for (UShort_t i = 0; i < sizeX ; ++i) {
    if (m_EventData->GetInnerX_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetInnerX_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetInnerXE(m_EventData->GetInnerX_E_DetectorNbr(i), m_EventData->GetInnerX_E_StripNbr(i), Energy);
      }
    }
  }

  // InnerZ Energy
  unsigned int sizeZ = m_EventData->GetInnerZMultEnergy();
  for (UShort_t i = 0; i < sizeZ ; ++i) {
    if (m_EventData->GetInnerZ_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetInnerZ_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetInnerZE(m_EventData->GetInnerZ_E_DetectorNbr(i), m_EventData->GetInnerZ_E_StripNbr(i), Energy);
      }
    }
  }

  //////
  // OuterX Energy
  sizeX = m_EventData->GetOuterXMultEnergy();
  for (UShort_t i = 0; i < sizeX ; ++i) {
    if (m_EventData->GetOuterX_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetOuterX_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetOuterXE(m_EventData->GetOuterX_E_DetectorNbr(i), m_EventData->GetOuterX_E_StripNbr(i), Energy);
      }
    }
  }

  // OuterZ Energy
  sizeZ = m_EventData->GetOuterZMultEnergy();
  for (UShort_t i = 0; i < sizeZ ; ++i) {
    if (m_EventData->GetOuterZ_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetOuterZ_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetOuterZE(m_EventData->GetOuterZ_E_DetectorNbr(i), m_EventData->GetOuterZ_E_StripNbr(i), Energy);
      }
    }
  }

  //////
  // Cluster Stage Energy
/*
  sizeFront = m_EventData->GetClusterMultTEnergy();
  for (UShort_t i = 0; i < sizeFront ; ++i) {
    if (m_EventData->GetCluster_TE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetCluster_TE_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterTE(m_EventData->GetCluster_TE_DetectorNbr(i), m_EventData->GetCluster_TE_StripNbr(i), Energy);
      }
    }
  }

  sizeBack = m_EventData->GetClusterMultLEnergy();
  for (UShort_t i = 0; i < sizeBack ; ++i) {
    if (m_EventData->GetCluster_LE_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetCluster_LE_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterLE(m_EventData->GetCluster_LE_DetectorNbr(i), m_EventData->GetCluster_LE_StripNbr(i), Energy);
      }
    }
  }
*/


  // ClusterX1 Energy
  unsigned int sizeClusterX1 = m_EventData->GetClusterX1MultEnergy();
  for (UShort_t i = 0; i < sizeClusterX1 ; ++i) {
    if (m_EventData->GetClusterX1_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetClusterX1_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterX1E(m_EventData->GetClusterX1_E_DetectorNbr(i), m_EventData->GetClusterX1_E_StripNbr(i), Energy);
      }
    }
  }

  // ClusterY1 Energy
  unsigned int sizeClusterY1 = m_EventData->GetClusterY1MultEnergy();
  for (UShort_t i = 0; i < sizeClusterY1 ; ++i) {
    if (m_EventData->GetClusterY1_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetClusterY1_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterY1E(m_EventData->GetClusterY1_E_DetectorNbr(i), m_EventData->GetClusterY1_E_StripNbr(i), Energy);
      }
    }
  }

  // ClusterX2 Energy
  unsigned int sizeClusterX2 = m_EventData->GetClusterX2MultEnergy();
  for (UShort_t i = 0; i < sizeClusterX2 ; ++i) {
    if (m_EventData->GetClusterX2_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetClusterX2_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterX2E(m_EventData->GetClusterX2_E_DetectorNbr(i), m_EventData->GetClusterX2_E_StripNbr(i), Energy);
      }
    }
  }

  // ClusterY2 Energy
  unsigned int sizeClusterY2 = m_EventData->GetClusterY2MultEnergy();
  for (UShort_t i = 0; i < sizeClusterY2 ; ++i) {
    if (m_EventData->GetClusterY2_E_Energy(i) > m_E_RAW_Threshold) {
      Double_t Energy = m_EventData->GetClusterY2_E_Energy(i);
      if (Energy > m_E_Threshold) {
        m_PreTreatedData->SetClusterY2E(m_EventData->GetClusterY2_E_DetectorNbr(i), m_EventData->GetClusterY2_E_StripNbr(i), Energy);
      }
    }
  }


}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::ReadAnalysisConfig() {
  bool ReadingStatus = false;

  // path to file
  string FileName = "./configs/ConfigTOGAXSI_SI.dat";

  // open analysis config file
  ifstream AnalysisConfigFile;
  AnalysisConfigFile.open(FileName.c_str());

  if (!AnalysisConfigFile.is_open()) {
    cout << " No ConfigTOGAXSI_SI.dat found: Default parameter loaded for Analayis " << FileName << endl;
    return;
  }
  cout << " Loading user parameter for Analysis from ConfigTOGAXSI_SI.dat " << endl;

  // Save it in a TAsciiFile
  TAsciiFile* asciiConfig = RootOutput::getInstance()->GetAsciiFileAnalysisConfig();
  asciiConfig->AppendLine("%%% ConfigTOGAXSI_SI.dat %%%");
  asciiConfig->Append(FileName.c_str());
  asciiConfig->AppendLine("");
  // read analysis config file
  string LineBuffer,DataBuffer,whatToDo;
  while (!AnalysisConfigFile.eof()) {
    // Pick-up next line
    getline(AnalysisConfigFile, LineBuffer);

    // search for "header"
    string name = "ConfigTOGAXSI_SI";
    if (LineBuffer.compare(0, name.length(), name) == 0) 
      ReadingStatus = true;

    // loop on tokens and data
    while (ReadingStatus ) {
      whatToDo="";
      AnalysisConfigFile >> whatToDo;

      // Search for comment symbol (%)
      if (whatToDo.compare(0, 1, "%") == 0) {
        AnalysisConfigFile.ignore(numeric_limits<streamsize>::max(), '\n' );
      }

      else if (whatToDo=="E_RAW_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_RAW_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_RAW_Threshold << endl;
      }

      else if (whatToDo=="E_THRESHOLD") {
        AnalysisConfigFile >> DataBuffer;
        m_E_Threshold = atof(DataBuffer.c_str());
        cout << whatToDo << " " << m_E_Threshold << endl;
      }

      else if (whatToDo=="E_FRONTBACK_MATCHING") {
        AnalysisConfigFile >> DataBuffer;
	m_StripEnergyMatching = atof(DataBuffer.c_str());
	cout << whatToDo << " " << m_StripEnergyMatching << endl;
      }

      else {
        ReadingStatus = false;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::Clear() {

  EventMultiplicity = 0;
  EventMultiplicity_cluster = 0;

  //DSSD
  //DetectorNumber.clear();
  DetectorNumberInnerX.clear();
  DetectorNumberInnerZ.clear();
  DetectorNumberOuterX.clear();
  DetectorNumberOuterZ.clear();
  DetectorNumberClusterX1.clear();
  DetectorNumberClusterY1.clear();
  DetectorNumberClusterX2.clear();
  DetectorNumberClusterY2.clear();

  InnerXStrip.clear();
  InnerZStrip.clear();
  OuterXStrip.clear();
  OuterZStrip.clear();
  ClusterX1Strip.clear();
  ClusterY1Strip.clear();
  ClusterX2Strip.clear();
  ClusterY2Strip.clear();  
//  ClusterStripT.clear();
//  ClusterStripL.clear();
  InnerXE.clear();
  InnerZE.clear();
  OuterXE.clear();
  OuterZE.clear();
  ClusterX1E.clear();
  ClusterY1E.clear();
  ClusterX2E.clear();
  ClusterY2E.clear();
//  ClusterTE.clear();
//  ClusterLE.clear();

  //Position Information
  InnerXPosX.clear();
  InnerXPosY.clear();
  InnerXPosZ.clear();
  InnerZPosX.clear();
  InnerZPosY.clear();
  InnerZPosZ.clear();
  OuterXPosX.clear();
  OuterXPosY.clear();
  OuterXPosZ.clear();
  OuterZPosX.clear();
  OuterZPosY.clear();
  OuterZPosZ.clear();
/*
  ClusterPosX.clear();
  ClusterPosY.clear();
  ClusterPosZ.clear();
*/
  ClusterX1PosX.clear();
  ClusterX1PosY.clear();
  ClusterX1PosZ.clear();
  ClusterY1PosX.clear();
  ClusterY1PosY.clear();
  ClusterY1PosZ.clear();
  ClusterX2PosX.clear();
  ClusterX2PosY.clear();
  ClusterX2PosZ.clear();
  ClusterY2PosX.clear();
  ClusterY2PosY.clear();
  ClusterY2PosZ.clear();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::ReadConfiguration(NPL::InputParser parser) {


  //InnerX Tracker     
  vector<NPL::InputBlock*> blocks_innerX = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","InnerX");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_innerX.size() << "inner X detectors found " << endl; 

  vector<string> coord = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_innerX.size() ; i++){
    if(blocks_innerX[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI inner X detector " << i+1 <<  endl;
    
      double R = blocks_innerX[i]->GetDouble("Radius","mm");
      double Z = blocks_innerX[i]->GetDouble("Z","mm");
      double Phi = blocks_innerX[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_innerX[i]->GetTVector3("Ref","mm");
      AddInnerXDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting on " << i+1 << " inner block " <<endl;
      exit(1);
    }
  }

  //InnerZ Tracker     
  vector<NPL::InputBlock*> blocks_innerZ = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","InnerZ");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_innerZ.size() << "inner Z detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_innerZ.size() ; i++){
    if(blocks_innerZ[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI inner Z detector " << i+1 <<  endl;
    
      double R = blocks_innerZ[i]->GetDouble("Radius","mm");
      double Z = blocks_innerZ[i]->GetDouble("Z","mm");
      double Phi = blocks_innerZ[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_innerZ[i]->GetTVector3("Ref","mm");
      AddInnerZDetector(R,Z,Phi,Ref);
    }
    else{
      cout << "ERROR: check your input file formatting on " << i+1 << " inner block " <<endl;
      exit(1);
    }
  }


  //OuterX Tracker     
  vector<NPL::InputBlock*> blocks_outerX = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","OuterX");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_outerX.size() << "outer X detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_outerX.size() ; i++){
    if(blocks_outerX[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI outer X detector " << i+1 <<  endl;
    
      double R = blocks_outerX[i]->GetDouble("Radius","mm");
      double Z = blocks_outerX[i]->GetDouble("Z","mm");
      double Phi = blocks_outerX[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_outerX[i]->GetTVector3("Ref","mm");
      AddOuterXDetector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " outer block " <<endl;
      exit(1);
    }
  }


  //OuterZ Tracker     
  vector<NPL::InputBlock*> blocks_outerZ = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","OuterZ");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_outerZ.size() << "outer Z detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_outerZ.size() ; i++){
    if(blocks_outerZ[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI outer Z detector " << i+1 <<  endl;
    
      double R = blocks_outerZ[i]->GetDouble("Radius","mm");
      double Z = blocks_outerZ[i]->GetDouble("Z","mm");
      double Phi = blocks_outerZ[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_outerZ[i]->GetTVector3("Ref","mm");
      AddOuterZDetector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " outer block " <<endl;
      exit(1);
    }
  }

  //Cluster Tracker     
/*
  vector<NPL::InputBlock*> blocks_cluster = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","Cluster");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_cluster.size() << "cluster detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_cluster.size() ; i++){
    if(blocks_cluster[i]->HasTokenList(coord)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI cluster detector " << i+1 <<  endl;
    
      double R = blocks_cluster[i]->GetDouble("Radius","mm");
      double Z = blocks_cluster[i]->GetDouble("Z","mm");
      double Phi = blocks_cluster[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_cluster[i]->GetTVector3("Ref","mm");
      AddClusterDetector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " cluster block " <<endl;
      exit(1);
    }
  }
*/


  //ClusterX Tracker     
  vector<NPL::InputBlock*> blocks_clusterX1 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterX1");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterX1.size() << "cluster X1 detectors found " << endl; 

  vector<string> coord_cluster = {"Radius","Z","Phi","Ref"};

  for(unsigned int i = 0 ; i < blocks_clusterX1.size() ; i++){
    if(blocks_clusterX1[i]->HasTokenList(coord_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI cluster X1 detector " << i+1 <<  endl;
    
      double R = blocks_clusterX1[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterX1[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterX1[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_clusterX1[i]->GetTVector3("Ref","mm");
      AddClusterX1Detector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " cluster X1 block " <<endl;
      exit(1);
    }
  }

  //ClusterY Tracker     
  vector<NPL::InputBlock*> blocks_clusterY1 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterY1");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterY1.size() << "cluster Y1 detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_clusterY1.size() ; i++){
    if(blocks_clusterY1[i]->HasTokenList(coord_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI cluster Y1 detector " << i+1 <<  endl;
    
      double R = blocks_clusterY1[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterY1[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterY1[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_clusterY1[i]->GetTVector3("Ref","mm");
      AddClusterY1Detector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " cluster Y1 block " <<endl;
      exit(1);
    }
  }


  //ClusterX2 Tracker     
  vector<NPL::InputBlock*> blocks_clusterX2 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterX2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterX2.size() << "cluster X2 detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_clusterX2.size() ; i++){
    if(blocks_clusterX2[i]->HasTokenList(coord_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI cluster X2 detector " << i+1 <<  endl;
    
      double R = blocks_clusterX2[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterX2[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterX2[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_clusterX2[i]->GetTVector3("Ref","mm");
      AddClusterX2Detector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " cluster X2 block " <<endl;
      exit(1);
    }
  }

  //ClusterY Tracker     
  vector<NPL::InputBlock*> blocks_clusterY2 = parser.GetAllBlocksWithTokenAndValue("TOGAXSI_SI","ClusterY2");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << "//// " << blocks_clusterY2.size() << "cluster Y2 detectors found " << endl; 

  for(unsigned int i = 0 ; i < blocks_clusterY2.size() ; i++){
    if(blocks_clusterY2[i]->HasTokenList(coord_cluster)){
      if(NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "////  TOGAXSI_SI cluster Y2 detector " << i+1 <<  endl;
    
      double R = blocks_clusterY2[i]->GetDouble("Radius","mm");
      double Z = blocks_clusterY2[i]->GetDouble("Z","mm");
      double Phi = blocks_clusterY2[i]->GetDouble("Phi","deg");
      TVector3 Ref = blocks_clusterY2[i]->GetTVector3("Ref","mm");
      AddClusterY2Detector(R,Z,Phi,Ref);
    }
    else{

      cout << "ERROR: check your input file formatting on " << i+1 << " cluster Y2 block " <<endl;
      exit(1);
    }
  }



  ReadAnalysisConfig();

}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::InitSpectra() {
//  m_Spectra = new TTOGAXSI_SISpectra(m_NumberOfDetectors);
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::FillSpectra() {
  m_Spectra -> FillRawSpectra(m_EventData);
  m_Spectra -> FillPreTreatedSpectra(m_PreTreatedData);
  m_Spectra -> FillPhysicsSpectra(m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::CheckSpectra() {
  m_Spectra->CheckSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::ClearSpectra() {
  // To be done
}



///////////////////////////////////////////////////////////////////////////
map< string , TH1*> TTOGAXSI_SIPhysics::GetSpectra() {
  if(m_Spectra)
    return m_Spectra->GetMapHisto();
  else{
    map< string , TH1*> empty;
    return empty;
  }
}

///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::WriteSpectra() {
  m_Spectra->WriteSpectra();
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::AddParameterToCalibrationManager() {
  CalibrationManager* Cal = CalibrationManager::getInstance();
/*  for (int i = 0; i < m_NumberOfDetectors; ++i) {
    Cal->AddParameter("TOGAXSI_SI", "D"+ NPL::itoa(i+1)+"_ENERGY","TOGAXSI_SI_D"+ NPL::itoa(i+1)+"_ENERGY");
    Cal->AddParameter("TOGAXSI_SI", "D"+ NPL::itoa(i+1)+"_TIME","TOGAXSI_SI_D"+ NPL::itoa(i+1)+"_TIME");
  }
*/
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::InitializeRootInputRaw() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchStatus("TOGAXSI_SI",  true );
  inputChain->SetBranchAddress("TOGAXSI_SI", &m_EventData );
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::InitializeRootInputPhysics() {
  TChain* inputChain = RootInput::getInstance()->GetChain();
  inputChain->SetBranchAddress("TOGAXSI_SI", &m_EventPhysics);
}



///////////////////////////////////////////////////////////////////////////
void TTOGAXSI_SIPhysics::InitializeRootOutput() {
  TTree* outputTree = RootOutput::getInstance()->GetTree();
  outputTree->Branch("TOGAXSI_SI", "TTOGAXSI_SIPhysics", &m_EventPhysics);
}



////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VDetector* TTOGAXSI_SIPhysics::Construct() {
  return (NPL::VDetector*) new TTOGAXSI_SIPhysics();
}



////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C"{
  class proxy_TOGAXSI_SI{
    public:
      proxy_TOGAXSI_SI(){
        NPL::DetectorFactory::getInstance()->AddToken("TOGAXSI_SI","TOGAXSI_SI");
        NPL::DetectorFactory::getInstance()->AddDetector("TOGAXSI_SI",TTOGAXSI_SIPhysics::Construct);
      }
  };

  proxy_TOGAXSI_SI p_TOGAXSI_SI;
}

