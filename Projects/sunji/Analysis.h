#ifndef Analysis_h
#define Analysis_h
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: XAUTHORX  contact address: XMAILX                        *
 *                                                                           *
 * Creation Date  : XMONTHX XYEARX                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class describe sunji analysis project using STARK X6 and CsI       *
 *  for missing mass calculation                                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "NPBeam.h"
#include "NPEnergyLoss.h"
#include "NPReaction.h"
#include "NPVAnalysis.h"
#include "RootInput.h"
#include "RootOutput.h"
#include "TCsIPhysics.h"
#include "TInitialConditions.h"
#include "TReactionConditions.h"
#include "TSTARKPhysics.h"
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>
#include <vector>

class Analysis : public NPL::VAnalysis {
 public:
  Analysis();
  ~Analysis();

 public:
  void Init();
  void TreatEvent();
  void End();

  void InitOutputBranch();
  void InitInputBranch();
  void ReInitValue();
  static NPL::VAnalysis* Construct();

 private:
  ////////////////////////////////////////////////////////////
  // vertex and beam information
  TVector3 vert;
  double BeamReacE;
  TVector3 BeamDirection;
  TVector3 BeamImpact;
  TRandom3 Rand;
  int Run;

  ////////////////////////////////////////////////////////////
  // Reaction information
  NPL::Reaction myReaction;
  double OriginalBeamEnergy; // AMeV
  double TargetThickness;
  double WindowsThickness;

  // Energy loss tables
  NPL::EnergyLoss* BeamTarget;
  NPL::EnergyLoss* BeamWindow;
  NPL::EnergyLoss* OutgoingTarget;
  NPL::EnergyLoss* OutgoingWindow;

  ////////////////////////////////////////////////////////////
  // for STARK X6
  int starkM;
  int starkType[20]; // type = 0 (X6), 1 (BB10), 2 (QQQ5)
  int starkDetN[20];
  int starkFStrN[20], starkBStrN[20];
  double starkUppE[20], starkDwnE[20], starkSumE[20];
  TVector3 starkHitPos[20];
  double starkThetaLab[20];
  double starkELab[20];

  // Branches and detectors
  TSTARKPhysics* STARK;
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Missing Mass calculation
  double MissingMass;        // Missing mass in MeV/c^2
  double MissingMassSq;      // Missing mass squared
  TVector3 OutgoingMomentum; // Outgoing particle momentum
  double OutgoingEnergy;     // Outgoing particle energy
  double OutgoingThetaLab;   // Outgoing particle lab angle
  TLorentzVector Beam4Vector;
  TLorentzVector Target4Vector;
  TLorentzVector Outgoing4Vector;
  TLorentzVector Recoil4Vector;
  TLorentzVector Missing4Vector;
  ////////////////////////////////////////////////////////////

  TInitialConditions* myInit;
  TReactionConditions* myReac;
};
#endif
