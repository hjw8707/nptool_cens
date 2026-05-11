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

#include "Analysis.h"

#include <iostream>

#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPFunction.h"
#include "NPOptionManager.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis()
    : BeamReacE(0),
      OriginalBeamEnergy(0),
      TargetThickness(0),
      WindowsThickness(0),
      BeamTarget(NULL),
      BeamWindow(NULL),
      OutgoingTarget(NULL),
      OutgoingWindow(NULL),
      STARK(NULL),
      dE(TMath::QuietNaN()),
      E(TMath::QuietNaN()),
      hitPos(0, 0, 0),
      MissingMass(TMath::QuietNaN()),
      RecoilExcitationEnergy(TMath::QuietNaN()),
      OutgoingEnergy(TMath::QuietNaN()),
      OutgoingThetaLab(TMath::QuietNaN()),
      myInit(NULL),
      myReac(NULL) {}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis() {
  if (BeamTarget) delete BeamTarget;
  if (BeamWindow) delete BeamWindow;
  if (OutgoingTarget) delete OutgoingTarget;
  if (OutgoingWindow) delete OutgoingWindow;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
  // initialize input and output branches
  cout << "!!!!!!!!!!!!!!! Initializing sunji analysis !!!!!!!!!!!!!!!!!!" << endl;

  // get STARK and CsI detector objects
  std::vector<std::string> detList = m_DetectorManager->GetDetectorList();
  for (auto it = detList.begin(); it != detList.end(); it++) {
    if ((*it) == "STARK") STARK = static_cast<TSTARKPhysics*>(m_DetectorManager->GetDetector("STARK"));
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
  if (!myReac || !myInit) return;

  vert = myReac->GetVertexPosition();
  if (TMath::IsNaN(vert.X()))  // no reaction at the target
    return;

  BeamImpact = vert;
  // BeamImpact.SetXYZ(0, 0, 0);  // assume the beam impact is at the origin
  // BeamDirection = myReac->GetBeamDirection();
  BeamDirection.SetXYZ(0, 0, 1);  // assume the beam direction is along the z-axis
  BeamReacE = myReac->GetBeamEnergy();
  myReaction.SetBeamEnergy(BeamReacE);

  // Get beam and target 4-vectors
  NPL::Particle* beam = myReaction.GetParticle1();
  NPL::Particle* target = myReaction.GetParticle2();
  NPL::Particle* outgoing = myReaction.GetParticle3();  // proton (light ejectile)
  NPL::Particle* recoil = myReaction.GetParticle4();    // nucleus (heavy ejectile)

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
  if (STARK && STARK->nGroup > 0) {
    dE = STARK->groupdE[0][0];
    E = STARK->groupEwithCsI[0];
    hitPos = STARK->sPosArr[STARK->groupFirstHitIdx[0]];
  } else {
    return;
  }
  // Calculate lab angle
  TVector3 hitDir = hitPos - BeamImpact;
  // TVector3 hitDir = myReac->GetParticleDirection(0);
  OutgoingThetaLab = hitDir.Angle(BeamDirection);
  TVector3 targetNormal = TVector3(-1, 0, 1);
  double thetaNormal = hitDir.Angle(targetNormal);
  // Correct for energy loss in target and window -> not implemented yet
  double measuredE = E;
  if (OutgoingWindow) {
    measuredE = OutgoingWindow->EvaluateInitialEnergy(measuredE, WindowsThickness, thetaNormal);
  }
  if (OutgoingTarget) {
    measuredE = OutgoingTarget->EvaluateInitialEnergy(measuredE, TargetThickness * 0.5, thetaNormal);
  }

  // Calculate outgoing particle momentum
  // OutgoingEnergy = myReac->GetKineticEnergy(0);  // OutgoingEnergy;  // MeV
  OutgoingEnergy = measuredE;
  double outgoingP = sqrt(OutgoingEnergy * OutgoingEnergy + 2 * OutgoingEnergy * outgoingMass);
  OutgoingMomentum = hitDir.Unit() * outgoingP;

  // Outgoing particle 4-vector
  double outgoingE = OutgoingEnergy + outgoingMass;
  Outgoing4Vector = TLorentzVector(OutgoingMomentum, outgoingE);

  // Missing 4-vector = Beam + Target - Outgoing
  Missing4Vector = Beam4Vector + Target4Vector - Outgoing4Vector;

  // Missing mass squared
  // Calculate recoil excitation energy from missing4vector and recoil mass
  MissingMass = Missing4Vector.M();
  RecoilExcitationEnergy = MissingMass - recoilMass;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("myReac", "RC", &myReac);
  RootOutput::getInstance()->GetTree()->Branch("vert", "TVector3", &vert);
  RootOutput::getInstance()->GetTree()->Branch("BeamReacE", &BeamReacE, "BeamReacE/D");
  RootOutput::getInstance()->GetTree()->Branch("BeamDirection", "TVector3", &BeamDirection);
  RootOutput::getInstance()->GetTree()->Branch("BeamImpact", "TVector3", &BeamImpact);
  RootOutput::getInstance()->GetTree()->Branch("Run", &Run, "Run/I");

  // Missing Mass branches
  RootOutput::getInstance()->GetTree()->Branch("dE", &dE, "dE/D");
  RootOutput::getInstance()->GetTree()->Branch("E", &E, "E/D");
  RootOutput::getInstance()->GetTree()->Branch("hitPos", "TVector3", &hitPos);
  RootOutput::getInstance()->GetTree()->Branch("MissingMass", &MissingMass, "MissingMass/D");
  RootOutput::getInstance()->GetTree()->Branch("RecoilExcitationEnergy", &RecoilExcitationEnergy,
                                               "RecoilExcitationEnergy/D");
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
  } else
    myReac = NULL;
}
////////////////////////////////////////////////////////////////////////////////
void Analysis::ReInitValue() {
  vert.SetXYZ(0, 0, 0);
  BeamReacE = 0;
  BeamDirection.SetXYZ(0, 0, 1);
  BeamImpact.SetXYZ(0, 0, 0);

  MissingMass = TMath::QuietNaN();
  RecoilExcitationEnergy = TMath::QuietNaN();
  OutgoingEnergy = TMath::QuietNaN();
  OutgoingThetaLab = TMath::QuietNaN();
  OutgoingMomentum.SetXYZ(0, 0, 0);
  dE = TMath::QuietNaN();
  E = TMath::QuietNaN();
  hitPos.SetXYZ(0, 0, 0);
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
