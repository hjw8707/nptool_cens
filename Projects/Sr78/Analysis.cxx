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
 *  This class describe  TiNA analysis project                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

#include "Analysis.h"

#include <iostream>
#include <numeric>

#include "NPAnalysisFactory.h"
#include "NPDetectorManager.h"
#include "NPFunction.h"
#include "NPOptionManager.h"
#include "TRegexp.h"
#include "TString.h"
#include "TTrackInfo.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
Analysis::Analysis()
    : BeamReacE(0),
      Ex(0),
      Qval(0),
      ELab(0),
      ThetaLab(0),
      ThetaCM(0),
      beamElossT(0),
      lightElossT(0),
      beamElossW(0),
      lightElossW(0),
      X(0),
      Y(0),
      Z(0),
      dE(0),
      dTheta(0),
      grapeM(0),
      grapeMadd(0),
      daliM(0),
      daliMadd(0),
      cacaoM(0),
      cacaoMadd(0),
      TiNA(NULL),
      GRAPE(NULL),
      DALI(NULL),
      CACAO(NULL),
      qState(-1),
      myInit(NULL),
      myReac(NULL) {}
////////////////////////////////////////////////////////////////////////////////
Analysis::~Analysis() {
  if (LightTarget) delete LightTarget;
  if (BeamTarget) delete BeamTarget;
  if (LightWindow) delete LightWindow;
  if (BeamWindow) delete BeamWindow;
  if (LightAl) delete LightAl;
  if (LightSi) delete LightSi;
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::Init() {
  // initialize input and output branches
  cout << "!!!!!!!!!!!!!!! Initializing analysis !!!!!!!!!!!!!!!!!!" << endl;
  // TTrackInfo* trackInfo = new TTrackInfo();
  // trackInfo->Dump();

  TReactionConditions* reactionConditions = new TReactionConditions();
  reactionConditions->Dump();

  TInitialConditions* initialConditions = new TInitialConditions();
  initialConditions->Dump();

  // get TiNA objects
  std::vector<std::string> detList = m_DetectorManager->GetDetectorList();
  for (auto it = detList.begin(); it != detList.end(); it++) {
    if ((*it) == "Tina")
      TiNA = static_cast<TTinaPhysics*>(m_DetectorManager->GetDetector("Tina"));
    else if ((*it) == "GRAPE")
      GRAPE = static_cast<TGRAPEPhysics*>(m_DetectorManager->GetDetector("GRAPE"));
    else if ((*it) == "Dali2")
      DALI = static_cast<TDali2Physics*>(m_DetectorManager->GetDetector("Dali2"));
    else if ((*it) == "CACAO")
      CACAO = static_cast<TCACAOPhysics*>(m_DetectorManager->GetDetector("CACAO"));
  }

  InitInputBranch();
  InitOutputBranch();

  // get reaction information
  myReaction.ReadConfigurationFile(NPOptionManager::getInstance()->GetReactionFile());
  OriginalBeamEnergy = myReaction.GetBeamEnergy();
  OriginalBeamImpact = TVector3(static_cast<Beam*>(myReaction.GetParticle1())->GetMeanX(),
                                static_cast<Beam*>(myReaction.GetParticle1())->GetMeanY(), 0);
  OriginalBeamDirection = TVector3(tan(static_cast<Beam*>(myReaction.GetParticle1())->GetMeanThetaX()),
                                   tan(static_cast<Beam*>(myReaction.GetParticle1())->GetMeanPhiY()), 1);

  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OriginalBeamEnergy=" << OriginalBeamEnergy << endl;
  cout << "OriginalBeamImpact.X()=" << OriginalBeamImpact.X()
       << " OriginalBeamDirection.X()=" << OriginalBeamDirection.X() << endl;

  // target thickness
  TargetThickness = m_DetectorManager->GetTargetThickness();
  // TargetThickness = 1.5;
  WindowsThickness = m_DetectorManager->GetFrontThickness();
  string TargetMaterial = m_DetectorManager->GetTargetMaterial();
  string WindowsMaterial = m_DetectorManager->GetFrontMaterial();

  std::cout << "Target material: " << TargetMaterial << std::endl;
  std::cout << "Target thickness: " << TargetThickness << std::endl;

  std::cout << "Windows material: " << WindowsMaterial << std::endl;
  std::cout << "Windows thickness: " << WindowsThickness << std::endl;

  // energy losses
  string light = NPL::ChangeNameToG4Standard(myReaction.GetParticle3()->GetName());
  string beam = NPL::ChangeNameToG4Standard(myReaction.GetParticle1()->GetName());

  LightTarget = NULL;
  BeamTarget = NULL;
  LightAl = NULL;
  LightSi = NULL;
  LightWindow = NULL;
  BeamWindow = NULL;

  if (OriginalBeamEnergy > 0 && !TargetMaterial.empty()) {
    LightTarget = new NPL::EnergyLoss(light + "_" + TargetMaterial + ".G4table", "G4Table", 100);
    BeamTarget = new NPL::EnergyLoss(beam + "_" + TargetMaterial + ".G4table", "G4Table", 100);
  }
  if (OriginalBeamEnergy > 0) {
    // LightAl = NPL::EnergyLoss(light+"_Al.G4table","G4Table",100);
    LightSi = new NPL::EnergyLoss(light + "_Si.G4table", "G4Table", 100);
    // LightSi = NPL::EnergyLoss(light+"_Si.SRIM","SRIM",100);
  }

  if (OriginalBeamEnergy > 0 && WindowsThickness) {
    BeamWindow = new NPL::EnergyLoss(beam + "_" + WindowsMaterial + ".G4table", "G4Table", 100);
    LightWindow = new NPL::EnergyLoss(light + "_" + WindowsMaterial + ".G4table", "G4Table", 100);
  } else {
    BeamWindow = NULL;
    LightWindow = NULL;
  }

  // initialize various parameters
  Rand = TRandom3();
  ThetaNormalTarget = 0;
  ThetaTiNASurface = 0;
  Energy = 0;
  ThetaTiNASurface = 0;

  ////////////////////////////////////////////////////////////
  // make addback table for DALI & GRAPE
  Double_t distLimit = 10. * cm;

  if (DALI) {
    ofstream fout;
    fout.open("DALIAddbackTable.out");
    for (int i = 0; i < DALI->GetNumberOfDetectors(); i++) {
      fout << (i + 1) << " ";
      TVector3 aPos = DALI->GetDALIPosition(i);
      std::vector<int> tempTable;
      int counter = 0;
      for (int j = 0; j < DALI->GetNumberOfDetectors(); j++) {
        if (j == i) continue;
        TVector3 bPos = DALI->GetDALIPosition(j);
        double dist = (aPos - bPos).Mag();
        if (dist < distLimit) {
          fout << (j + 1) << " ";
          tempTable.push_back(j + 1);
        }
      }
      daliabTable.push_back(tempTable);
      fout << "\n";
    }
  }

  if (CACAO) {
    ofstream fout;
    fout.open("CACAOAddbackTable.out");
    for (int i = 0; i < CACAO->GetNumberOfDetectors(); i++) {
      for (int j = 0; j < 4; j++) {
        int csiNumber = (i + 1) * 100 + (j + 1);
        fout << csiNumber << " ";
        TVector3 aPos = CACAO->GetCsIPosition(i, j);
        std::vector<int> tempTable;
        int counter = 0;
        for (int k = 0; k < CACAO->GetNumberOfDetectors(); k++) {
          for (int l = 0; l < 4; l++) {
            int csiNumber2 = (k + 1) * 100 + (l + 1);
            if (csiNumber2 == csiNumber) continue;
            TVector3 bPos = CACAO->GetCsIPosition(k, l);
            double dist = (aPos - bPos).Mag();
            if (dist < distLimit) {
              fout << csiNumber2 << " ";
              tempTable.push_back(csiNumber2);
            }
          }
        }
        cacaoabTable[csiNumber] = tempTable;
        fout << "\n";
      }
    }
  }

  ////////////////////////////////////////////////////////////
  // Load charge state probability
  ifstream fin("charge.dat");
  double sum = 0;
  while (fin.good()) {
    double temp;
    fin >> temp;
    sum += temp;
    csProb.push_back(sum);
  }
  fin.close();

  std::cout << " Charge state distribution " << std::endl;
  for (int i = 0; i < csProb.size(); i++) {
    std::cout << "Q = Z - " << i << ": " << csProb[i] << std::endl;
  }
  ////////////////////////////////////////////////////////////
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::TreatEvent() {
  // Reinitiate calculated variable
  ReInitValue();

  // double zImpact = Rand.Uniform(-TargetThickness*0.5,TargetThickness*0.5);
  double zImpact = 0;
  TVector3 BeamPos(myInit->GetIncidentPositionX(), myInit->GetIncidentPositionY(), myInit->GetIncidentPositionZ());
  TVector3 BeamDir = myInit->GetBeamDirection();
  BeamPos += (BeamDir * ((zImpact - BeamPos.z()) / BeamDir.z()));  // at target center (z = 0)
  BeamPos[0] += Rand.Gaus(0, 0.1 * mm);                            // X resolution at the focal plane
  BeamPos[1] += Rand.Gaus(0, 0.1 * mm);                            // Y resolution at the focal plane

  double BeamNormalTarget, E_beam;
  bool flagPhys = false;
  // std::cout << "before myreac" << std::endl;
  if (myReac) {
    vert = myReac->GetVertexPosition();
    if (TMath::IsNaN(vert.X()))  // no reaction at the target
      return;
    if (flagPhys) {  // in reality
      BeamImpact = vert;
      BeamDirection = myReac->GetBeamDirection();
      BeamReacE = myReac->GetBeamEnergy();
      myReaction.SetBeamEnergy(BeamReacE);
    } else {  // in measurement
      BeamImpact = BeamPos;
      BeamDirection = BeamDir;

      BeamNormalTarget = BeamDirection.Angle(TVector3(0, 0, 1));
      E_beam = myInit->GetIncidentInitialKineticEnergy();  // OriginalBeamEnergy;
      E_beam *= Rand.Gaus(1, 0.01);                        // Energy resolution (0.01 = 1%)
      if (BeamWindow) {
        double Eafter = BeamWindow->Slow(E_beam, WindowsThickness, BeamNormalTarget);
        beamElossW = E_beam - Eafter;
        E_beam = Eafter;
      }
      if (BeamTarget) {
        double Eafter = BeamTarget->Slow(E_beam, TargetThickness * 0.5 - zImpact,
                                         // TargetThickness*0.5-BeamImpact.z(),
                                         BeamNormalTarget);
        beamElossT = E_beam - Eafter;
        E_beam = Eafter;
      }

      BeamReacE = E_beam;
      myReaction.SetBeamEnergy(E_beam);
    }
  }
  // std::cout << "after myreac" << std::endl;
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  //////////////////////////// LOOP on TTT //////////////////
  vector<double> HitEvent[9];
  /////////////////  Add multiple strips hit per detector together
  if (flagPhys) {
    ThetaNormalTarget = 0;
    HitDirection = myReac->GetParticleDirection(0);
    Double_t zdiff = 200. - vert.Z();
    TDir = HitDirection;
    TDir.SetMag(zdiff / TDir.Z());
    TDir += vert;
    ThetaLab = HitDirection.Angle(BeamDirection);
    ThetaNormalTarget = HitDirection.Angle(TVector3(0, 0, 1));

    ELab = myReac->GetKineticEnergy(0);
    Ex = myReaction.ReconstructRelativistic(ELab, ThetaLab);

    ThetaCM = myReaction.EnergyLabToThetaCM(ELab, ThetaLab) / deg;
    ThetaLab = ThetaLab / deg;

    QValue = myReaction.GetQValue();
    QValue = QValue - Ex;
  }
  // std::cout << "after flagPhys" << std::endl;
  if (TiNA) {
    for (int i = 0; i < TiNA->TTT_E.size(); i++)
      if (TiNA->TTT_E[i] > 0 && TiNA->TTT_X[i] > 0 && TiNA->TTT_X[i] < 129 && TiNA->TTT_Y[i] > 0 &&
          TiNA->TTT_Y[i] < 129) {
        if (HitEvent[0].size() == 0) {
          HitEvent[0].push_back(TiNA->SquareTelescopeNumber.at(i));
          HitEvent[1].push_back(TiNA->TTT_E[i]);
          HitEvent[2].push_back(TiNA->TTT_X[i]);
          HitEvent[3].push_back(TiNA->TTT_Y[i]);
          HitEvent[4].push_back(max(0., TiNA->Pad_E[i]));
          HitEvent[5].push_back(1);  // record multiplicity on detector
          HitEvent[6].push_back(TiNA->GetPositionOfInteraction(i).X());
          HitEvent[7].push_back(TiNA->GetPositionOfInteraction(i).Y());
          HitEvent[8].push_back(TiNA->GetPositionOfInteraction(i).Z());
        } else
          for (unsigned short j = 0; j < HitEvent[0].size(); j++)
            if (TiNA->SquareTelescopeNumber.at(i) == HitEvent[0].at(j)) {  // summing for the same telescope data
              HitEvent[1].at(j) += TiNA->TTT_E[i];
              HitEvent[2].at(j) += TiNA->TTT_X[i];
              HitEvent[3].at(j) += TiNA->TTT_Y[i];
              if (HitEvent[4].at(j) <= 0) HitEvent[4].at(j) += max(0., TiNA->Pad_E[i]);
              HitEvent[5].at(j) += 1;
              HitEvent[6].at(j) += TiNA->GetPositionOfInteraction(i).X();
              HitEvent[7].at(j) += TiNA->GetPositionOfInteraction(i).Y();
              HitEvent[8].at(j) += TiNA->GetPositionOfInteraction(i).Z();
              break;
            } else if (j == HitEvent[0].size() - 1) {  // new element for the diff. telescope data
              HitEvent[0].push_back(TiNA->SquareTelescopeNumber.at(i));
              HitEvent[1].push_back(TiNA->TTT_E[i]);
              HitEvent[2].push_back(TiNA->TTT_X[i]);
              HitEvent[3].push_back(TiNA->TTT_Y[i]);
              HitEvent[4].push_back(max(0., TiNA->Pad_E[i]));
              HitEvent[5].push_back(1);  // record multiplicity on detector
              HitEvent[6].push_back(TiNA->GetPositionOfInteraction(i).X());
              HitEvent[7].push_back(TiNA->GetPositionOfInteraction(i).Y());
              HitEvent[8].push_back(TiNA->GetPositionOfInteraction(i).Z());
            }
      }

    /////////////////  Add multiple strips together
    for (int i = 0; i < HitEvent[0].size(); i++) {
      SquareMultiplicity = HitEvent[5].at(i);  // multiplicity
      double Pad_E = HitEvent[4].at(i);
      if (Pad_E < 0) Pad_E = 0;
      //      X = HitEvent[6].at(i)/SquareMultiplicity;
      //      Y = HitEvent[7].at(i)/SquareMultiplicity;
      X = HitEvent[6].at(i) / SquareMultiplicity - 0.796;
      Y = HitEvent[7].at(i) / SquareMultiplicity - 0.791;
      Z = HitEvent[8].at(i) / SquareMultiplicity;
      HitPosition.SetXYZ(X, Y, Z);

      //==================================================//
      // Part 1 : Impact Angle
      ThetaTiNASurface = 0;
      ThetaNormalTarget = 0;

      HitDirection = HitPosition - BeamImpact;
      ThetaLab = HitDirection.Angle(BeamDirection);
      ThetaNormalTarget = HitDirection.Angle(TVector3(0, 0, 1));
      //==================================================//

      //==================================================//
      // Part 2 : Impact Energy
      Energy = ELab = 0;
      Energy = HitEvent[1].at(i) + Pad_E;

      // Target Correction
      if (LightWindow) {
        double Eafter = LightWindow->EvaluateInitialEnergy(Energy, WindowsThickness, ThetaNormalTarget);
        lightElossW = Eafter - Energy;
        ELab = Eafter;
      } else
        ELab = Energy;

      if (LightTarget) {
        double Eafter = LightTarget->EvaluateInitialEnergy(ELab, TargetThickness * 0.5 - zImpact,
                                                           // TargetThickness*0.5 - BeamImpact.z(),
                                                           ThetaNormalTarget);
        lightElossT = Eafter - ELab;
        ELab = Eafter;
      }

      ///////////////////////////////////////////////////////////////////
      //      ELab = myReac->GetKineticEnergy(0); // !!!! real triton energy !!!!
      ///////////////////////////////////////////////////////////////////

      dE = myReac->GetKineticEnergy(0) - ELab;
      dTheta = (myInit->GetThetaLab_WorldFrame(0) * deg - ThetaLab) / deg;
      //==================================================//

      //==================================================//
      // Part 3 : Excitation Energy Calculation
      Ex = myReaction.ReconstructRelativistic(ELab, ThetaLab);
      //==================================================//

      //==================================================//
      // Part 4 : Theta CM Calculation
      ThetaCM = myReaction.EnergyLabToThetaCM(ELab, ThetaLab) / deg;
      ThetaLab = ThetaLab / deg;
      //==================================================//

      //==================================================//
      // Part 5 : QValue Calculation
      QValue = myReaction.GetQValue();
      QValue = QValue - Ex;
      //==================================================//
    }  //

    for (int i = 0; i < 9; i++) HitEvent[i].clear();

    //////////////////////////////////////////////////////////
    ///////////     LOOP ON YY1 ///////////////
    //////////////////////////////////////////////////////
    for (int i = 0; i < TiNA->YY1_E.size(); i++)
      if (TiNA->YY1_E[i] > 0 && TiNA->YY1_R[i] > 0) {
        if (HitEvent[0].size() == 0) {
          HitEvent[0].push_back(TiNA->TrapezTelescopeNumber.at(i));
          HitEvent[1].push_back(TiNA->YY1_E[i]);
          HitEvent[2].push_back(TiNA->YY1_R[i]);
          HitEvent[3].push_back(TiNA->YY1_S[i]);
          HitEvent[4].push_back(max(0., TiNA->CsI_E[i]));
          HitEvent[5].push_back(1);  // record multiplicity on detector
          HitEvent[6].push_back(TiNA->GetRingPositionOfInteraction(i).X());
          HitEvent[7].push_back(TiNA->GetRingPositionOfInteraction(i).Y());
          HitEvent[8].push_back(TiNA->GetRingPositionOfInteraction(i).Z());
        } else
          for (unsigned short j = 0; j < HitEvent[0].size(); j++)
            if (TiNA->TrapezTelescopeNumber.at(i) == HitEvent[0].at(j)) {
              HitEvent[1].at(j) += TiNA->YY1_E[i];
              HitEvent[2].at(j) += TiNA->YY1_R[i];
              HitEvent[3].at(j) += TiNA->YY1_S[i];
              if (HitEvent[4].at(j) <= 0) HitEvent[4].at(j) += max(0., TiNA->CsI_E[i]);
              HitEvent[5].at(j) += 1;
              HitEvent[6].at(j) += TiNA->GetRingPositionOfInteraction(i).X();
              HitEvent[7].at(j) += TiNA->GetRingPositionOfInteraction(i).Y();
              HitEvent[8].at(j) += TiNA->GetRingPositionOfInteraction(i).Z();
              break;
            } else if (j == HitEvent[0].size() - 1) {
              HitEvent[0].push_back(TiNA->TrapezTelescopeNumber.at(i));
              HitEvent[1].push_back(TiNA->YY1_E[i]);
              HitEvent[2].push_back(TiNA->YY1_R[i]);
              HitEvent[3].push_back(TiNA->YY1_S[i]);
              HitEvent[4].push_back(max(0., TiNA->CsI_E[i]));
              HitEvent[5].push_back(1);  // record multiplicity on detector
              HitEvent[6].push_back(TiNA->GetRingPositionOfInteraction(i).X());
              HitEvent[7].push_back(TiNA->GetRingPositionOfInteraction(i).Y());
              HitEvent[8].push_back(TiNA->GetRingPositionOfInteraction(i).Z());
            }
      }

    for (int i = 0; i < HitEvent[0].size(); i++) {
      TrapMultiplicity = HitEvent[5].at(i);
      double CsI_E = HitEvent[4].at(i);
      if (CsI_E < 0) CsI_E = 0;
      X = HitEvent[6].at(i) / TrapMultiplicity;
      Y = HitEvent[7].at(i) / TrapMultiplicity;
      Z = HitEvent[8].at(i) / TrapMultiplicity;
      HitPosition.SetXYZ(X, Y, Z);

      //==================================================//
      // Part 1 : Impact Angle
      ThetaTiNASurface = 0;
      ThetaNormalTarget = 0;

      HitDirection = HitPosition - BeamImpact;
      ThetaLab = HitDirection.Angle(BeamDirection);
      ThetaNormalTarget = HitDirection.Angle(TVector3(0, 0, 1));
      //==================================================//

      //==================================================//
      // Part 2 : Impact Energy
      Energy = ELab = 0;
      Energy = HitEvent[1].at(i) + CsI_E;

      if (LightTarget)
        ELab = LightTarget->EvaluateInitialEnergy(Energy, TargetThickness * 0.5 - zImpact,
                                                  // TargetThickness*0.5-BeamImpact.z(),
                                                  ThetaNormalTarget);
      else
        ELab = Energy;

      if (LightWindow) ELab = LightWindow->EvaluateInitialEnergy(ELab, WindowsThickness, ThetaNormalTarget);

      dE = 0;  // myInit->GetKineticEnergy(3) - ELab ;
      dTheta = 0;
      //==================================================//

      //==================================================//
      // Part 3 : Excitation Energy Calculation
      Ex = myReaction.ReconstructRelativistic(ELab, ThetaLab);
      //==================================================//

      //==================================================//
      // Part 4 : Theta CM Calculation
      ThetaCM = myReaction.EnergyLabToThetaCM(ELab, ThetaLab) / deg;
      ThetaLab = ThetaLab / deg;
      //==================================================//

      //==================================================//
      // Part 5 : QValue Calculation
      QValue = myReaction.GetQValue();
      QValue = QValue - Ex;
      //==================================================//
    }  // end loop TiNA//*/

    for (int i = 0; i < 9; i++) HitEvent[i].clear();
  }
  // std::cout << "start gamma" << std::endl;
  ////////////////////////////////////////////////////////////
  // Gamma detector
  ////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////
  // Fragment momentum vector?
  TRegexp regexNu("^[a-zA-Z]+[0-9]+");
  TRegexp regexEl("^[a-zA-Z]+");
  NPL::Particle frag;  //(38, 78);
  TVector3 fragBetaVec;
  if (myReac) {
    for (int i = 0; i < myReac->GetParticleMultiplicity(); i++) {
      TString pName(myReac->GetParticleName(i));
      Ssiz_t len = 0;
      if (regexNu.Index(pName, &len) < 0)  // not matched to the regex (not nuclei)
        continue;

      std::string nuName(pName(0, len));
      regexEl.Index(pName, &len);
      std::string nuEl(pName(0, len));
      std::string nuMa(pName(len, nuName.size() - len));
      nuName = nuMa;
      nuName += nuEl;
      frag.SetUp(nuName);

      // double fragKE = myReac->GetKineticEnergy(i);
      frag.SetKineticEnergy(myReac->GetKineticEnergy(i));
      fragBetaVec = myReac->GetParticleDirection(i);  // only direction
      fragMom = fragBetaVec;                          // only direction
      break;
    }

    if (fragBetaVec.Mag() > 0) {  // nuclei found
      ////////////////////////////////////////////////////////////
      // charge state selection
      Double_t prob = Rand.Uniform(0, 1);
      for (qState = 0; qState < csProb.size(); qState++) {
        if (prob < csProb[qState]) break;
      }
      frag.EnergyToBrho(frag.GetZ() - qState);                          // frag.GetZ() - 1, 2, ...
      frag.SetBrho(Rand.Gaus(frag.GetBrho(), frag.GetBrho() * 0.001));  // 0.1% resolution
      frag.BrhoToEnergy(frag.GetZ() - qState);                          // frag.GetZ() - 1, 2, ...
      frag.EnergyToBeta();
      ////////////////////////////////////////////////////////////

      fragMom.SetMag(frag.GetBrho());
      fragBetaVec.SetMag(frag.GetBeta());
    }
  }
  fragBeta = fragBetaVec;

  // std::cout << "start dali" << std::endl;
  //////////////////////////////////////////////////
  // Dali
  //////////////////////////////////////////////////
  if (DALI && DALI->DetectorNumber.size() > 0) {
    for (int i = 0; i < DALI->DetectorNumber.size(); i++) {
      if (DALI->Energy[i] > 0 && DALI->Time[i] > 0) {
        daliCN[daliM] = DALI->DetectorNumber[i];
        daliE[daliM] = DALI->Energy[i];
        daliEdc[daliM] = GetDopplerCorrectedEnergy(daliE[daliM], DALI->GetDALIPosition(daliCN[daliM] - 1), fragBetaVec);
        daliT[daliM++] = DALI->Time[i];
      }
    }

    //////////////////////////////////////////////////
    // addback
    // add hit modules
    for (int i = 0; i < DALI->DetectorNumber.size(); i++) {
      if (DALI->Energy[i] > 0 && DALI->Time[i] > 0) {
        daliabN.push_back(DALI->DetectorNumber[i]);
        daliabE.push_back(DALI->Energy[i]);
      }
    }
    // sort with E (using index array)
    std::vector<double> daliabETemp = daliabE;
    std::vector<int> sorted;
    for (size_t i = 0; i < daliabETemp.size(); i++) {
      sorted.push_back(std::max_element(daliabETemp.begin(), daliabETemp.end()) - daliabETemp.begin());
      daliabETemp[sorted[i]] = -1;
    }

    Double_t tempE[100], tempEdc[100], tempN[100];
    while (!sorted.empty()) {
      Double_t en = daliabE[sorted[0]];
      Int_t num = daliabN[sorted[0]];
      sorted.erase(sorted.begin());

      for (size_t j = 0; j < sorted.size(); j++) {
        Int_t ffindex = std::find(daliabTable[num - 1].begin(), daliabTable[num - 1].end(), daliabN[sorted[j]]) -
                        daliabTable[num - 1].begin();
        if (ffindex < daliabTable[num - 1].size()) {
          en += daliabE[sorted[j]];
          sorted.erase(sorted.begin() + j);
          j--;
          if (sorted.empty()) break;
        }
      }

      tempN[daliMadd] = num;
      tempE[daliMadd] = en;
      tempEdc[daliMadd++] = GetDopplerCorrectedEnergy(en, DALI->GetDALIPosition(num - 1), fragBetaVec);
    }

    //////////////////////////////////////////////////
    // re-sort
    Int_t ind[100];
    TMath::Sort(daliMadd, tempEdc, ind);
    for (Int_t i = 0; i < daliMadd; i++) {
      daliCNadd[i] = tempN[ind[i]];
      daliEadd[i] = tempE[ind[i]];
      daliEdcadd[i] = tempEdc[ind[i]];
    }
    //////////////////////////////////////////////////
  }
  ////////////////////////////////////////////////////////////
  // std::cout << "start cacao" << std::endl;
  //////////////////////////////////////////////////
  // cacao
  //////////////////////////////////////////////////
  if (CACAO && CACAO->nhit > 0) {
    for (int i = 0; i < CACAO->nhit; i++) {
      if (CACAO->E[i] > 0) {  //&& CACAO->T[i] > 0) { What is the negative T?
        int detNumber = CACAO->detN[i];
        int csINumber = CACAO->CsIN[i];
        cacaoCN[cacaoM] = detNumber * 100 + csINumber;
        ////////////////////////////////////////
        // position
        cacaoR[cacaoM] = CACAO->GetCsIPosition(detNumber - 1, csINumber - 1).Mag();
        cacaoTheta[cacaoM] = CACAO->GetCsIPosition(detNumber - 1, csINumber - 1).Theta();
        ////////////////////////////////////////
        cacaoE[cacaoM] = CACAO->E[i];
        cacaoEdc[cacaoM] =
            GetDopplerCorrectedEnergy(cacaoE[cacaoM], CACAO->GetCsIPosition(detNumber - 1, csINumber - 1), fragBetaVec);
        cacaoT[cacaoM++] = CACAO->T[i];

        // add hit modules for addback
        cacaoabN.push_back(detNumber * 100 + csINumber);
        cacaoabE.push_back(CACAO->E[i]);
      }
    }
    //////////////////////////////////////////////////
    // addback
    // sort with E (using index array)
    std::vector<double> cacaoabETemp = cacaoabE;
    std::vector<int> sorted;
    for (size_t i = 0; i < cacaoabETemp.size(); i++) {
      sorted.push_back(std::max_element(cacaoabETemp.begin(), cacaoabETemp.end()) - cacaoabETemp.begin());
      cacaoabETemp[sorted[i]] = -1;
    }

    Double_t tempE[100], tempEdc[100], tempN[100];
    while (!sorted.empty()) {
      Double_t en = cacaoabE[sorted[0]];
      Int_t num = cacaoabN[sorted[0]];
      sorted.erase(sorted.begin());

      for (size_t j = 0; j < sorted.size(); j++) {
        Int_t ffindex = std::find(cacaoabTable[num].begin(), cacaoabTable[num].end(), cacaoabN[sorted[j]]) -
                        cacaoabTable[num].begin();
        if (ffindex < cacaoabTable[num].size()) {
          en += cacaoabE[sorted[j]];
          sorted.erase(sorted.begin() + j);
          j--;
          if (sorted.empty()) break;
        }
      }
      tempN[cacaoMadd] = num;
      tempE[cacaoMadd] = en;
      Int_t detNumber = num / 100;
      Int_t csINumber = num % 100;
      tempEdc[cacaoMadd++] =
          GetDopplerCorrectedEnergy(en, CACAO->GetCsIPosition(detNumber - 1, csINumber - 1), fragBetaVec);
    }

    //////////////////////////////////////////////////
    // re-sort
    Int_t ind[100];
    TMath::Sort(cacaoMadd, tempEdc, ind);
    for (Int_t i = 0; i < cacaoMadd; i++) {
      cacaoCNadd[i] = tempN[ind[i]];
      cacaoEadd[i] = tempE[ind[i]];
      cacaoEdcadd[i] = tempEdc[ind[i]];
      ////////////////////////////////////////
      // position
      cacaoRadd[i] = CACAO->GetDetPosition(cacaoCNadd[i] - 1).Mag();
      cacaoThetaadd[i] = CACAO->GetDetPosition(cacaoCNadd[i] - 1).Theta();
      ////////////////////////////////////////
    }
    //////////////////////////////////////////////////
  }
  ////////////////////////////////////////////////////////////
  // std::cout << "start grape" << std::endl;
  //////////////////////////////////////////////////
  // GRAPE
  //////////////////////////////////////////////////
  if (GRAPE) {
    for (int i = 0; i < GRAPE->Gamma_Energy.size(); i++) {
      if (GRAPE->Gamma_Energy[i] > 0 && GRAPE->Gamma_Time[i] > 0) {
        grapeCN[grapeM] = GRAPE->GRAPE_Number[i];
        grapeSN[grapeM] = GRAPE->Segment_Number[i];
        grapeE[grapeM] = GRAPE->Gamma_Energy[i];
        // grapeEdc[grapeM] = GetDopplerCorrectedEnergy(
        //     grapeE[grapeM], GRAPE->GetSegmentPosition(grapeCN[grapeM], 0, grapeSN[grapeM]), fragBetaVec);
        grapeT[grapeM++] = GRAPE->Gamma_Time[i];
      }
    }

    Double_t tempE[100], tempEdc[100], tempCN[100], tempSN[100];
    ////////////////////////////////////////////////////////////
    // addback for each detector (module)
    for (int i = 0; i < GRAPE->GetNumberOfDetectors(); i++) {
      std::vector<double> grapeabE;
      std::vector<int> grapeabSN;
      for (int j = 0; j < GRAPE->Gamma_Energy.size(); j++) {
        if (GRAPE->Gamma_Energy[j] > 0 && GRAPE->Gamma_Time[j] > 0) {
          if (GRAPE->GRAPE_Number[j] == (i + 1)) {
            grapeabE.push_back(GRAPE->Gamma_Energy[j]);
            grapeabSN.push_back(GRAPE->Segment_Number[j]);
          }
        }
      }
      if (grapeabE.empty()) continue;

      // find maximum segment
      size_t maxJ = std::max_element(grapeabE.begin(), grapeabE.end()) - grapeabE.begin();
      // TVector3 Pos = GRAPE->GetSegmentPosition(i + 1, 0, grapeabSN[maxJ]);
      TVector3 Pos;

      tempCN[grapeMadd] = i + 1;
      tempSN[grapeMadd] = grapeabSN[maxJ];
      tempE[grapeMadd] = std::accumulate(grapeabE.begin(), grapeabE.end(), 0.);
      tempEdc[grapeMadd] = GetDopplerCorrectedEnergy(tempE[grapeMadd], Pos, fragBetaVec);
      grapeMadd++;
    }

    //////////////////////////////////////////////////
    // re-sort
    Int_t ind[100];
    TMath::Sort(grapeMadd, tempEdc, ind);
    for (Int_t i = 0; i < grapeMadd; i++) {
      grapeCNadd[i] = tempCN[ind[i]];
      grapeCNadd[i] = tempSN[ind[i]];
      grapeEadd[i] = tempE[ind[i]];
      grapeEdcadd[i] = tempEdc[ind[i]];
    }
    //////////////////////////////////////////////////
  }
}

////////////////////////////////////////////////////////////////////////////////
void Analysis::End() {}
////////////////////////////////////////////////////////////////////////////////
void Analysis::InitOutputBranch() {
  RootOutput::getInstance()->GetTree()->Branch("vert", "TVector3", &vert);
  RootOutput::getInstance()->GetTree()->Branch("fragBeta", "TVector3", &fragBeta);
  RootOutput::getInstance()->GetTree()->Branch("fragMom", "TVector3", &fragMom);
  RootOutput::getInstance()->GetTree()->Branch("BeamReacE", &BeamReacE, "BeamReacE/D");
  RootOutput::getInstance()->GetTree()->Branch("TDir", "TVector3", &TDir);

  RootOutput::getInstance()->GetTree()->Branch("Ex", &Ex, "Ex/D");
  RootOutput::getInstance()->GetTree()->Branch("ELab", &ELab, "ELab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaLab", &ThetaLab, "ThetaLab/D");
  RootOutput::getInstance()->GetTree()->Branch("ThetaCM", &ThetaCM, "ThetaCM/D");

  RootOutput::getInstance()->GetTree()->Branch("beamElossT", &beamElossT, "beamElossT/D");
  RootOutput::getInstance()->GetTree()->Branch("lightElossT", &lightElossT, "lightElossT/D");
  RootOutput::getInstance()->GetTree()->Branch("beamElossW", &beamElossW, "beamElossW/D");
  RootOutput::getInstance()->GetTree()->Branch("lightElossW", &lightElossW, "lightElossW/D");

  RootOutput::getInstance()->GetTree()->Branch("Run", &Run, "Run/I");
  RootOutput::getInstance()->GetTree()->Branch("HitPos", "TVector3", &HitPosition);
  RootOutput::getInstance()->GetTree()->Branch("HitDir", "TVector3", &HitDirection);
  RootOutput::getInstance()->GetTree()->Branch("TrapMultiplicity", &TrapMultiplicity, "TrapMultiplicity/I");
  RootOutput::getInstance()->GetTree()->Branch("SquareMultiplicity", &SquareMultiplicity, "SquareMultiplicity/I");
  RootOutput::getInstance()->GetTree()->Branch("dE", &dE, "dE/D");
  RootOutput::getInstance()->GetTree()->Branch("dTheta", &dTheta, "dTheta/D");
  RootOutput::getInstance()->GetTree()->Branch("QValue", &QValue, "QValue/D");

  RootOutput::getInstance()->GetTree()->Branch("qState", &qState, "qState/I");

  if (GRAPE) {
    RootOutput::getInstance()->GetTree()->Branch("grapeM", &grapeM, "grapeM/I");
    RootOutput::getInstance()->GetTree()->Branch("grapeCN", grapeCN, "grapeCN[grapeM]/I");
    RootOutput::getInstance()->GetTree()->Branch("grapeSN", grapeSN, "grapeSN[grapeM]/I");
    RootOutput::getInstance()->GetTree()->Branch("grapeE", grapeE, "grapeE[grapeM]/D");
    RootOutput::getInstance()->GetTree()->Branch("grapeEdc", grapeEdc, "grapeEdc[grapeM]/D");
    RootOutput::getInstance()->GetTree()->Branch("grapeT", grapeT, "grapeT[grapeM]/D");

    RootOutput::getInstance()->GetTree()->Branch("grapeMadd", &grapeMadd, "grapeMadd/I");
    RootOutput::getInstance()->GetTree()->Branch("grapeCNadd", grapeCNadd, "grapeCNadd[grapeMadd]/I");
    RootOutput::getInstance()->GetTree()->Branch("grapeSNadd", grapeSNadd, "grapeSNadd[grapeMadd]/I");
    RootOutput::getInstance()->GetTree()->Branch("grapeEadd", grapeEadd, "grapeEadd[grapeMadd]/D");
    RootOutput::getInstance()->GetTree()->Branch("grapeEdcadd", grapeEdcadd, "grapeEdcadd[grapeMadd]/D");
    RootOutput::getInstance()->GetTree()->Branch("grapeTadd", grapeTadd, "grapeTadd[grapeMadd]/D");
  }

  if (DALI) {
    RootOutput::getInstance()->GetTree()->Branch("daliM", &daliM, "daliM/I");
    RootOutput::getInstance()->GetTree()->Branch("daliCN", daliCN, "daliCN[daliM]/I");
    RootOutput::getInstance()->GetTree()->Branch("daliE", daliE, "daliE[daliM]/D");
    RootOutput::getInstance()->GetTree()->Branch("daliEdc", daliEdc, "daliEdc[daliM]/D");
    RootOutput::getInstance()->GetTree()->Branch("daliT", daliT, "daliT[daliM]/D");

    RootOutput::getInstance()->GetTree()->Branch("daliMadd", &daliMadd, "daliMadd/I");
    RootOutput::getInstance()->GetTree()->Branch("daliCNadd", daliCNadd, "daliCNadd[daliMadd]/I");
    RootOutput::getInstance()->GetTree()->Branch("daliEadd", daliEadd, "daliEadd[daliMadd]/D");
    RootOutput::getInstance()->GetTree()->Branch("daliEdcadd", daliEdcadd, "daliEdcadd[daliMadd]/D");
    RootOutput::getInstance()->GetTree()->Branch("daliTadd", daliTadd, "daliTadd[daliMadd]/D");
  }

  if (CACAO) {
    RootOutput::getInstance()->GetTree()->Branch("cacaoM", &cacaoM, "cacaoM/I");
    RootOutput::getInstance()->GetTree()->Branch("cacaoCN", cacaoCN, "cacaoCN[cacaoM]/I");
    RootOutput::getInstance()->GetTree()->Branch("cacaoR", cacaoR, "cacaoR[cacaoM]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoTheta", cacaoTheta, "cacaoTheta[cacaoM]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoE", cacaoE, "cacaoE[cacaoM]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoEdc", cacaoEdc, "cacaoEdc[cacaoM]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoT", cacaoT, "cacaoT[cacaoM]/D");

    RootOutput::getInstance()->GetTree()->Branch("cacaoMadd", &cacaoMadd, "cacaoMadd/I");
    RootOutput::getInstance()->GetTree()->Branch("cacaoCNadd", cacaoCNadd, "cacaoCNadd[cacaoMadd]/I");
    RootOutput::getInstance()->GetTree()->Branch("cacaoRadd", cacaoRadd, "cacaoRadd[cacaoM]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoThetaadd", cacaoThetaadd, "cacaoThetaadd[cacaoM]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoEadd", cacaoEadd, "cacaoEadd[cacaoMadd]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoEdcadd", cacaoEdcadd, "cacaoEdcadd[cacaoMadd]/D");
    RootOutput::getInstance()->GetTree()->Branch("cacaoTadd", cacaoTadd, "cacaoTadd[cacaoMadd]/D");
  }
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
  fragBeta.SetXYZ(0, 0, 0);
  fragMom.SetXYZ(0, 0, 0);
  TDir.SetXYZ(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  BeamReacE = 0;

  Ex = TMath::QuietNaN();
  ELab = TMath::QuietNaN();
  ThetaLab = TMath::QuietNaN();
  ThetaCM = TMath::QuietNaN();
  X = TMath::QuietNaN();
  Y = TMath::QuietNaN();
  Z = TMath::QuietNaN();
  HitPosition.SetXYZ(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  HitDirection.SetXYZ(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  TrapMultiplicity = 0;
  SquareMultiplicity = 0;
  dE = TMath::QuietNaN();
  dTheta = TMath::QuietNaN();

  beamElossT = TMath::QuietNaN();
  lightElossT = TMath::QuietNaN();
  beamElossW = TMath::QuietNaN();
  lightElossW = TMath::QuietNaN();

  qState = -1;

  daliabN.clear();
  daliabE.clear();

  cacaoabN.clear();
  cacaoabE.clear();

  grapeM = grapeMadd = daliM = daliMadd = cacaoM = cacaoMadd = 0;
  for (int i = 0; i < 400; i++) {
    grapeCN[i] = grapeSN[i] = grapeCNadd[i] = grapeSNadd[i] = 0;
    daliCN[i] = daliCNadd[i] = 0;
    cacaoCN[i] = cacaoCNadd[i] = 0;
    grapeE[i] = grapeEdc[i] = grapeT[i] = TMath::QuietNaN();
    grapeEadd[i] = grapeEdcadd[i] = grapeTadd[i] = TMath::QuietNaN();
    daliE[i] = daliEdc[i] = daliT[i] = TMath::QuietNaN();
    daliEadd[i] = daliEdcadd[i] = daliTadd[i] = TMath::QuietNaN();
    cacaoE[i] = cacaoEdc[i] = cacaoT[i] = TMath::QuietNaN();
    cacaoEadd[i] = cacaoEdcadd[i] = cacaoTadd[i] = TMath::QuietNaN();
    cacaoR[i] = cacaoTheta[i] = cacaoRadd[i] = cacaoThetaadd[i] = TMath::QuietNaN();
  }
}

////////////////////////////////////////////////////////////////////////////////
//            Construct Method to be pass to the DetectorFactory              //
////////////////////////////////////////////////////////////////////////////////
NPL::VAnalysis* Analysis::Construct() { return (NPL::VAnalysis*)new Analysis(); }

double Analysis::GetDopplerCorrectedEnergy(double energy, TVector3 position, TVector3 beta) {
  // renorm pos vector
  TLorentzVector m_GammaLV;
  position.SetMag(1);
  m_GammaLV.SetPx(energy * position.X());
  m_GammaLV.SetPy(energy * position.Y());
  m_GammaLV.SetPz(energy * position.Z());
  m_GammaLV.SetE(energy);
  m_GammaLV.Boost(-beta);
  return m_GammaLV.Energy();
}

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
