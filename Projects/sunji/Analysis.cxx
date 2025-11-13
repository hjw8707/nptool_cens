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

#include <iostream>

#include "Analysis.h"

#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPFunction.h"
#include "NPOptionManager.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis()
    : BeamReacE(0), OriginalBeamEnergy(0), TargetThickness(0), WindowsThickness(0), BeamTarget(NULL), BeamWindow(NULL),
      OutgoingTarget(NULL), OutgoingWindow(NULL), starkM(0), STARK(NULL), MissingMass(0), MissingMassSq(0),
      OutgoingEnergy(0), OutgoingThetaLab(0), myInit(NULL), myReac(NULL) {}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis() {
  if (BeamTarget)
    delete BeamTarget;
  if (BeamWindow)
    delete BeamWindow;
  if (OutgoingTarget)
    delete OutgoingTarget;
  if (OutgoingWindow)
    delete OutgoingWindow;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
  // initialize input and output branches
  cout << "!!!!!!!!!!!!!!! Initializing sunji analysis !!!!!!!!!!!!!!!!!!" << endl;

  // get STARK and CsI detector objects
  std::vector<std::string> detList = m_DetectorManager->GetDetectorList();
  for (auto it = detList.begin(); it != detList.end(); it++) {
    if ((*it) == "STARK")
      STARK = static_cast<TSTARKPhysics*>(m_DetectorManager->GetDetector("STARK"));
  }

  InitInputBranch();
  InitOutputBranch();

  // get reaction information
  myReaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = myReaction.GetBeamEnergy();
  std::cout << "Original beam energy: " << OriginalBeamEnergy << " AMeV" << std::endl;

  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  WindowsThickness = m_DetectorManager->GetFrontThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  string WindowsMaterial = m_DetectorManager->GetFrontMaterial();

  std::cout << "Target material: " << TargetMaterial << std::endl;
  std::cout << "Target thickness: " << TargetThickness << " mm" << std::endl;
  std::cout << "Windows material: " << WindowsMaterial << std::endl;
  std::cout << "Windows thickness: " << WindowsThickness << " mm" << std::endl;

  // energy losses
  string outgoing = NPL::ChangeNameToG4Standard(myReaction.GetParticle3()->GetName());
  string beam = NPL::ChangeNameToG4Standard(myReaction.GetParticle1()->GetName());

  BeamTarget = NULL;
  BeamWindow = NULL;
  OutgoingTarget = NULL;
  OutgoingWindow = NULL;

  if (OriginalBeamEnergy > 0 && !TargetMaterial.empty()) {
    BeamTarget = new NPL::EnergyLoss(beam + "_" + TargetMaterial + ".G4table", "G4Table", 100);
    OutgoingTarget = new NPL::EnergyLoss(outgoing + "_" + TargetMaterial + ".G4table", "G4Table", 100);
  }

  if (OriginalBeamEnergy > 0 && WindowsThickness > 0 && !WindowsMaterial.empty()) {
    BeamWindow = new NPL::EnergyLoss(beam + "_" + WindowsMaterial + ".G4table", "G4Table", 100);
    OutgoingWindow = new NPL::EnergyLoss(outgoing + "_" + WindowsMaterial + ".G4table", "G4Table", 100);
  }

  // initialize random number generator
  Rand = TRandom3();

  std::cout << "Analysis initialization complete!" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
  // Reinitiate calculated variable
  ReInitValue();

  // Get reaction vertex and beam information
  if (!myReac || !myInit)
    return;

  vert = myReac->GetVertexPosition();
  if (TMath::IsNaN(vert.X())) // no reaction at the target
    return;

  BeamImpact = vert;
  BeamDirection = myReac->GetBeamDirection();
  BeamReacE = myReac->GetBeamEnergy();
  myReaction.SetBeamEnergy(BeamReacE);

  // Get beam and target 4-vectors
  NPL::Particle* beam = myReaction.GetParticle1();
  NPL::Particle* target = myReaction.GetParticle2();
  NPL::Particle* outgoing = myReaction.GetParticle3();
  NPL::Particle* recoil = myReaction.GetParticle4();

  double beamMass = beam->Mass();
  double targetMass = target->Mass();
  double outgoingMass = outgoing->Mass();
  double recoilMass = recoil->Mass();

  // Beam 4-vector (in lab frame)
  TVector3 beamMom(0, 0, sqrt(BeamReacE * BeamReacE + 2 * BeamReacE * beamMass));
  double beamE = BeamReacE + beamMass;
  Beam4Vector = TLorentzVector(beamMom, beamE);

  // Target 4-vector (at rest in lab frame)
  Target4Vector = TLorentzVector(0, 0, 0, targetMass);

  // Process STARK X6 hits (only X6 type = 0)
  starkM = 0;
  if (STARK && STARK->nhit > 0) {
    for (int i = 0; i < STARK->nhit; i++) {
      // Only process X6 detectors (type = 0)
      if (STARK->type[i] != 0)
        continue;

      if (STARK->sumE[i] > 0 && STARK->hPosArr.size() > i) {
        starkType[starkM] = STARK->type[i];
        starkDetN[starkM] = STARK->detN[i];
        starkFStrN[starkM] = STARK->fStrN[i];
        starkBStrN[starkM] = STARK->bStrN[i];
        starkUppE[starkM] = STARK->uppE[i];
        starkDwnE[starkM] = STARK->dwnE[i];
        starkSumE[starkM] = STARK->sumE[i];
        starkHitPos[starkM] = STARK->hPosArr[i];

        // Calculate lab angle
        TVector3 hitDir = starkHitPos[starkM] - BeamImpact;
        starkThetaLab[starkM] = hitDir.Angle(BeamDirection);

        // Energy loss correction for outgoing particle
        double thetaNormal = hitDir.Angle(TVector3(0, 0, 1));
        double measuredE = starkSumE[starkM]; // MeV

        // Correct for energy loss in target and window
        if (OutgoingTarget) {
          measuredE = OutgoingTarget->EvaluateInitialEnergy(measuredE, TargetThickness * 0.5, thetaNormal);
        }
        if (OutgoingWindow) {
          measuredE = OutgoingWindow->EvaluateInitialEnergy(measuredE, WindowsThickness, thetaNormal);
        }

        starkELab[starkM] = measuredE;
        starkM++;
      }
    }
  }

  // // Calculate Missing Mass using STARK X6 data
  // if (starkM > 0) {
  //   // Use the highest energy X6 hit
  //   int bestHit = 0;
  //   double maxE = 0;
  //   for (int i = 0; i < starkM; i++) {
  //     if (starkELab[i] > maxE) {
  //       maxE = starkELab[i];
  //       bestHit = i;
  //     }
  //   }

  //   OutgoingEnergy = starkELab[bestHit];
  //   OutgoingThetaLab = starkThetaLab[bestHit];
  //   TVector3 hitDir = starkHitPos[bestHit] - BeamImpact;
  //   hitDir = hitDir.Unit();

  //   // Calculate outgoing particle momentum
  //   double outgoingKE = OutgoingEnergy; // MeV
  //   double outgoingP = sqrt(outgoingKE * outgoingKE + 2 * outgoingKE * outgoingMass);
  //   OutgoingMomentum = hitDir * outgoingP;

  //   // Outgoing particle 4-vector
  //   double outgoingE = outgoingKE + outgoingMass;
  //   Outgoing4Vector = TLorentzVector(OutgoingMomentum, outgoingE);

  //   // Missing 4-vector = Beam + Target - Outgoing
  //   Missing4Vector = Beam4Vector + Target4Vector - Outgoing4Vector;

  //   // Missing mass squared
  //   MissingMassSq = Missing4Vector.M2();

  //   // Missing mass (take square root if positive)
  //   if (MissingMassSq >= 0) {
  //     MissingMass = sqrt(MissingMassSq);
  //   }
  //   else {
  //     MissingMass = -sqrt(-MissingMassSq); // imaginary mass (shouldn't happen physically)
  //   }
  // }
  // else {
  //   // No STARK hit, set missing mass to NaN
  //   MissingMass = TMath::QuietNaN();
  //   MissingMassSq = TMath::QuietNaN();
  // }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("vert", "TVector3", &vert);
  RootOutput::getInstance()->GetTree()->Branch("BeamReacE", &BeamReacE, "BeamReacE/D");
  RootOutput::getInstance()->GetTree()->Branch("BeamDirection", "TVector3", &BeamDirection);
  RootOutput::getInstance()->GetTree()->Branch("BeamImpact", "TVector3", &BeamImpact);
  RootOutput::getInstance()->GetTree()->Branch("Run", &Run, "Run/I");

  // Missing Mass branches
  RootOutput::getInstance()->GetTree()->Branch("MissingMass", &MissingMass, "MissingMass/D");
  RootOutput::getInstance()->GetTree()->Branch("MissingMassSq", &MissingMassSq, "MissingMassSq/D");
  RootOutput::getInstance()->GetTree()->Branch("OutgoingEnergy", &OutgoingEnergy, "OutgoingEnergy/D");
  RootOutput::getInstance()->GetTree()->Branch("OutgoingThetaLab", &OutgoingThetaLab, "OutgoingThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("OutgoingMomentum", "TVector3", &OutgoingMomentum);
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::InitInputBranch() {
  RootInput::getInstance()->GetChain()->SetBranchStatus("Run", true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("Run", &Run);
  RootInput::getInstance()->GetChain()->SetBranchStatus("InitialConditions", true);
  RootInput::getInstance()->GetChain()->SetBranchAddress("InitialConditions", &myInit);
  myInit = new TInitialConditions();

  if (RootInput::getInstance()->GetChain()->GetBranch("ReactionConditions")) {
    RootInput::getInstance()->GetChain()->SetBranchStatus("ReactionConditions", true);
    RootInput::getInstance()->GetChain()->SetBranchAddress("ReactionConditions", &myReac);
    myReac = new TReactionConditions();
  }
  else
    myReac = NULL;
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue() {
  vert.SetXYZ(0, 0, 0);
  BeamReacE = 0;
  BeamDirection.SetXYZ(0, 0, 1);
  BeamImpact.SetXYZ(0, 0, 0);

  MissingMass = TMath::QuietNaN();
  MissingMassSq = TMath::QuietNaN();
  OutgoingEnergy = TMath::QuietNaN();
  OutgoingThetaLab = TMath::QuietNaN();
  OutgoingMomentum.SetXYZ(0, 0, 0);

  starkM = 0;
  for (int i = 0; i < 20; i++) {
    starkType[i] = 0;
    starkDetN[i] = 0;
    starkFStrN[i] = 0;
    starkBStrN[i] = 0;
    starkUppE[i] = TMath::QuietNaN();
    starkDwnE[i] = TMath::QuietNaN();
    starkSumE[i] = TMath::QuietNaN();
    starkHitPos[i].SetXYZ(0, 0, 0);
    starkThetaLab[i] = TMath::QuietNaN();
    starkELab[i] = TMath::QuietNaN();
  }
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct() { return (NPL::VAnalysis*)new Analysis(); }

////////////////////////////////////////////////////////////////////////////////
//            Registering the construct method to the factory                 //
////////////////////////////////////////////////////////////////////////////////
extern "C" {
class proxy_analysis {
 public:
  proxy_analysis() { NPL::AnalysisFactory::getInstance()->SetConstructor(Analysis::Construct); }
};

proxy_analysis p_analysis;
}
