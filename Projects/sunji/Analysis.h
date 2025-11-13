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

#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>

#include <vector>

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
  double OriginalBeamEnergy;  // AMeV
  double TargetThickness;
  double WindowsThickness;

  // Energy loss tables
  NPL::EnergyLoss* BeamTarget;
  NPL::EnergyLoss* BeamWindow;
  NPL::EnergyLoss* OutgoingTarget;
  NPL::EnergyLoss* OutgoingWindow;

  ////////////////////////////////////////////////////////////
  // Branches and detectors
  TSTARKPhysics* STARK;
  ////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////
  // Missing Mass calculation
  Double_t dE;
  Double_t E;
  TVector3 hitPos;
  double MissingMass;             // Missing mass in MeV/c^2
  double RecoilExcitationEnergy;  // Recoil excitation energy in MeV
  TVector3 OutgoingMomentum;      // Outgoing particle momentum
  double OutgoingEnergy;          // Outgoing particle energy
  double OutgoingThetaLab;        // Outgoing particle lab angle
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
