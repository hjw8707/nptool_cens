/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author : Adrien MAtta contact: matta@lpccaen.in2p3.fr            *
 *                                                                           *
 * Creation Date   : April 2024                                              *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "NPCore.h"
#include "NPFunction.h"
#include "NPInelasticBreakup.h"
#include "NPOptionManager.h"

// Use CLHEP System of unit and Physical Constant
#include "NPGlobalSystemOfUnits.h"
#include "NPPhysicalConstants.h"
#ifdef NP_SYSTEM_OF_UNITS_H
using namespace NPUNITS;
#endif

// ROOT
#include "TF1.h"

ClassImp(InelasticBreakup)

    //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
    InelasticBreakup::InelasticBreakup() {
  //------------- Default Constructor -------------

  //
  fBeamEnergy = 0;
  fThetaCM = 0;
  fQValue = 0;
  //  fVerboseLevel = NPOptionManager::getInstance()->GetVerboseLevel();
  fVerboseLevel = 1; // NPOptionManager::getInstance()->GetVerboseLevel();

  fCrossSectionHist = NULL;
  fDoubleDifferentialCrossSectionHist = NULL;

  fLabCrossSection = false; // flag if the provided cross-section is in the lab or not
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
InelasticBreakup::~InelasticBreakup() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void InelasticBreakup::GenerateEvent(double& ThetaLight, double& PhiLight, double& KineticEnergyLight,
                                     double& ThetaLabTarget, double& KineticEnergyLabTarget, double& ThetaLabHeavy,
                                     double& KineticEnergyLabHeavy) {
  // 2-body relativistic kinematics: direct + inverse
  // EnergieLabTarget,4 : lab energy in MeV of the 2 ejectiles
  // ThetaLabTarget,4   : angles in rad
  // case of inverse kinematics

  double theta = fThetaCM;
  if (fParticle1.Mass() > fParticle2.Mass()) {
    theta = M_PI - fThetaCM;
    // fThetaCM = M_PI - fThetaCM;
  }

  // shoot light first
  ThetaLight = -1;
  KineticEnergyLight = -1;
  while (ThetaLight < 0 || KineticEnergyLight < 0) {
    ThetaLight = gRandom->Gaus(fMeanAngleLight, fSigmaAngleLight);
    KineticEnergyLight = gRandom->Gaus(fMeanEnergyLight, fSigmaEnergyLight);
  }

  PhiLight = gRandom->Uniform() * 2. * M_PI;
  double totalELight = fParticle3.Mass() + KineticEnergyLight;
  double momentumLight = sqrt(KineticEnergyLight * KineticEnergyLight + 2 * fParticle3.Mass() * KineticEnergyLight);
  TLorentzVector lightLV(momentumLight * sin(ThetaLight) * cos(PhiLight),
                         momentumLight * sin(ThetaLight) * sin(PhiLight), momentumLight * cos(ThetaLight), totalELight);

  // Substract light LV to Beam
  double momentumBeam = sqrt(fBeamEnergy * fBeamEnergy + 2 * fBeam.Mass() * fBeamEnergy);
  TLorentzVector beamLV(0, 0, momentumBeam, fBeamEnergy + fBeam.Mass());

  auto newLV = beamLV - lightLV;

  // Two body kine elastic scattering with heavy
  fReaction.SetBeamEnergy(newLV.E() - newLV.Mag());
  fReaction.ShootRandomThetaCM();
  fReaction.KineRelativistic(ThetaLabTarget, KineticEnergyLabTarget, ThetaLabHeavy, KineticEnergyLabHeavy);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void InelasticBreakup::ReadConfigurationFile(string Path) {
  ifstream InelasticBreakupFile;
  string GlobalPath = getenv("NPTOOL");
  string StandardPath = GlobalPath + "/Inputs/EventGenerator/" + Path;
  InelasticBreakupFile.open(Path.c_str());
  if (!InelasticBreakupFile.is_open()) {
    InelasticBreakupFile.open(StandardPath.c_str());
    if (InelasticBreakupFile.is_open()) {
      Path = StandardPath;
    }
    else {
      cout << "InelasticBreakup File " << Path << " not found" << endl;
      exit(1);
    }
  }
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}
////////////////////////////////////////////////////////////////////////////////
Particle InelasticBreakup::GetParticle(string name, NPL::InputParser parser) {
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("DefineParticle", name);
  unsigned int size = blocks.size();
  if (size == 0)
    return NPL::Particle(name);
  else if (size == 1) {
    cout << " -- User defined nucleus " << name << " -- " << endl;
    vector<string> token = {"SubPart", "BindingEnergy"};
    if (blocks[0]->HasTokenList(token)) {
      NPL::Particle N(name, blocks[0]->GetVectorString("SubPart"), blocks[0]->GetDouble("BindingEnergy", "MeV"));
      if (blocks[0]->HasToken("ExcitationEnergy"))
        N.SetExcitationEnergy(blocks[0]->GetDouble("ExcitationEnergy", "MeV"));
      if (blocks[0]->HasToken("SpinParity"))
        N.SetSpinParity(blocks[0]->GetString("SpinParity").c_str());
      if (blocks[0]->HasToken("Spin"))
        N.SetSpin(blocks[0]->GetDouble("Spin", ""));
      if (blocks[0]->HasToken("Parity"))
        N.SetParity(blocks[0]->GetString("Parity").c_str());
      if (blocks[0]->HasToken("LifeTime"))
        N.SetLifeTime(blocks[0]->GetDouble("LifeTime", "s"));

      cout << " -- -- -- -- -- -- -- -- -- -- --" << endl;
      return N;
    }
  }
  else {
    NPL::SendErrorAndExit("NPL::InelasticBreakup", "Too many nuclei define with the same name");
  }

  return (NPL::Particle());
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
void InelasticBreakup::ReadConfigurationFile(NPL::InputParser parser) {

  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("InelasticBreakup");
  if (blocks.size() > 0 && NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// Inelastic Breakup reaction found " << endl;

  vector<string> token1 = {"Beam",           "Target",          "Light",
                           "Heavy",          "MeanEnergyLight", "SigmaEnergyLight",
                           "MeanAngleLight", "SigmaAngleLight", "CrossSectionPath"};
  double CSHalfOpenAngleMin = 0 * deg;
  double CSHalfOpenAngleMax = 180 * deg;
  for (unsigned int i = 0; i < blocks.size(); i++) {
    if (blocks[i]->HasTokenList(token1)) {
      int v = NPOptionManager::getInstance()->GetVerboseLevel();
      NPOptionManager::getInstance()->SetVerboseLevel(0);
      // Read the beam block
      fBeam.ReadConfigurationFile(parser);
      NPOptionManager::getInstance()->SetVerboseLevel(v);

      fBeamEnergy = fBeam.GetEnergy();
      // set the particle
      fReaction.SetParticle1(fBeam);
      // fParticle1 = GetParticle(blocks[i]->GetString("Beam"), parser);
      fParticle2 = GetParticle(blocks[i]->GetString("Target"), parser);
      fParticle3 = GetParticle(blocks[i]->GetString("Light"), parser);
      fParticle4 = GetParticle(blocks[i]->GetString("Heavy"), parser);

      fReaction.GetParticle1()->SetUp(blocks[i]->GetString("Heavy"));
      fReaction.SetParticle2(GetParticle(blocks[i]->GetString("Target"), parser));
      fReaction.SetParticle3(GetParticle(blocks[i]->GetString("Target"), parser));
      fReaction.SetParticle4(GetParticle(blocks[i]->GetString("Heavy"), parser));
      fReaction.SetBeamEnergy(fBeamEnergy);
      fReaction.PreventDirect();
      fReaction.initializePrecomputeVariable();

      fParticle3 = GetParticle(blocks[i]->GetString("Light"), parser);
      fMeanEnergyLight = blocks[i]->GetDouble("MeanEnergyLight", "MeV");
      fSigmaEnergyLight = blocks[i]->GetDouble("SigmaEnergyLight", "MeV");
      fMeanAngleLight = blocks[i]->GetDouble("MeanAngleLight", "deg");
      fSigmaAngleLight = blocks[i]->GetDouble("SigmaAngleLight", "deg");
    }
    else {
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }

    if (blocks[i]->HasToken("ExcitationEnergyLight"))
      fExcitationLight = blocks[i]->GetDouble("ExcitationEnergyLight", "MeV");

    if (blocks[i]->HasToken("ExcitationEnergyHeavy"))
      fExcitationHeavy = blocks[i]->GetDouble("ExcitationEnergyHeavy", "MeV");

    if (blocks[i]->HasToken("CrossSectionPath")) {
      vector<string> file = blocks[i]->GetVectorString("CrossSectionPath");
      TH1D* CStemp = Read1DProfile(file[0], file[1]);

      // multiply CStemp by sin(theta)
      TF1* fsin = new TF1("sin", Form("1/(sin(x*%f/180.))", M_PI), 0, 180);
      CStemp->Divide(fsin, 1);
      fReaction.SetCrossSectionHist(CStemp);
      delete fsin;
    }

    if (blocks[i]->HasToken("LabCrossSectionPath")) {
      fLabCrossSection = true;

      vector<string> file = blocks[i]->GetVectorString("LabCrossSectionPath");
      TH1D* CStemp = Read1DProfile(file[0], file[1]);

      // multiply CStemp by sin(theta)
      TF1* fsin = new TF1("sin", Form("1/(sin(x*%f/180.))", M_PI), 0, 180);
      CStemp->Divide(fsin, 1);
      SetCrossSectionHist(CStemp);
      delete fsin;
    }

    if (blocks[i]->HasToken("DoubleDifferentialCrossSectionPath")) {
      vector<string> file = blocks[i]->GetVectorString("DoubleDifferentialCrossSectionPath");
      TH2F* CStemp = Read2DProfile(file[0], file[1]);

      // multiply CStemp by sin(theta)
      // X axis is theta CM
      // Y axis is beam energy
      // Division affect only X axis
      TF1* fsin = new TF1("sin", Form("1/(sin(x*%f/180.))", M_PI), 0, 180);
      CStemp->Divide(fsin, 1);

      SetDoubleDifferentialCrossSectionHist(CStemp);
      delete fsin;
    }

    if (blocks[i]->HasToken("HalfOpenAngleMin")) {
      CSHalfOpenAngleMin = blocks[i]->GetDouble("HalfOpenAngleMin", "deg");
    }
    if (blocks[i]->HasToken("HalfOpenAngleMax")) {
      CSHalfOpenAngleMax = blocks[i]->GetDouble("HalfOpenAngleMax", "deg");
    }
  }
  SetCSAngle(CSHalfOpenAngleMin / deg, CSHalfOpenAngleMax / deg);
  cout << "\033[0m";
}

////////////////////////////////////////////////////////////////////////////////////////////
void InelasticBreakup::SetCSAngle(double CSHalfOpenAngleMin, double CSHalfOpenAngleMax) {
  if (fCrossSectionHist) {
    for (int i = 0; i < fCrossSectionHist->GetNbinsX(); i++) {
      if (fCrossSectionHist->GetBinCenter(i) > CSHalfOpenAngleMax ||
          fCrossSectionHist->GetBinCenter(i) < CSHalfOpenAngleMin) {
        fCrossSectionHist->SetBinContent(i, 0);
      }
    }
  }
}

