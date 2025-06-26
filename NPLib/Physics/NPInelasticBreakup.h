#ifndef NPInelasticBreakUp_h
#define NPInelasticBreakUp_h
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
// C++ header
#include <string>

// NPL
#include "NPBeam.h"
#include "NPInputParser.h"
#include "NPParticle.h"
#include "NPReaction.h"
using namespace NPL;

// ROOT header
#include "NPReaction.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TRandom.h"
#include "TVector3.h"

using namespace std;

namespace NPL {
  class InelasticBreakup {

   public: // Constructors and Destructors
    InelasticBreakup();
    ~InelasticBreakup();

   public: // Various Method
    Particle GetParticle(string name, NPL::InputParser parser);
    void ReadConfigurationFile(string Path);
    void ReadConfigurationFile(NPL::InputParser);

   private:
    TRandom* RandGen;

   private:
    int fVerboseLevel;

   private: // use for Monte Carlo simulation
    bool fLabCrossSection;

   private:
    NPL::Reaction fReaction; // two body part of the kinematic
    Beam fBeam;              // Beam
    Beam fParticle1;         // Beam
    Particle fParticle2;     // Target
    Particle fParticle3;     // Light ejectile
    Particle fParticle4;     // Heavy ejectile
    double fMeanEnergyLight;
    double fSigmaEnergyLight;
    double fMeanAngleLight;
    double fSigmaAngleLight;
    double fExcitationHeavy;
    double fExcitationLight;

    double fQValue;                            // Q-value in MeV
    double fBeamEnergy;                        // Beam energy in MeV
    double fThetaCM;                           // Center-of-mass angle in radian
    TH1D* fCrossSectionHist;                   // Differential cross section in CM frame
    TH2F* fDoubleDifferentialCrossSectionHist; // Diff. CS CM frame vs Beam E
    TH1D* fExcitationEnergyHist;               // Distribution of Excitation energy
   public:
    // Getters and Setters
    void SetBeamEnergy(const double& eBeam) { fBeamEnergy = eBeam; }
    void SetThetaCM(const double& angle) { fThetaCM = angle; }

    void SetExcitationLight(const double& exci) { fExcitationLight = exci; }
    void SetExcitationHeavy(const double& exci) { fExcitationHeavy = exci; }
    void SetVerboseLevel(const int& verbose) { fVerboseLevel = verbose; }
    void SetCrossSectionHist(TH1D* CrossSectionHist) {
      delete fCrossSectionHist;
      fCrossSectionHist = CrossSectionHist;
    }

    void SetDoubleDifferentialCrossSectionHist(TH2F* CrossSectionHist) {
      fDoubleDifferentialCrossSectionHist = CrossSectionHist;
    }
    double GetBeamEnergy() const { return fBeamEnergy; }
    double GetThetaCM() const { return fThetaCM; }
    double GetQValue() const { return fQValue; }
    Particle* GetParticleBeam() { return &fBeam; }
    Particle* GetParticle1() { return &fParticle1; }
    Particle* GetParticleTarget() { return &fParticle2; }
    Particle* GetParticleLight() { return &fParticle3; }
    Particle* GetParticleHeavy() { return &fParticle4; }
    /* Particle* GetNucleus1() { return GetParticle1(); }
     Particle* GetNucleus2() { return GetParticle2(); }
     Particle* GetNucleus3() { return GetParticle3(); }
     Particle* GetNucleus4() { return GetParticle4(); }
 */
    TH1D* GetCrossSectionHist() const { return fCrossSectionHist; }
    int GetVerboseLevel() const { return fVerboseLevel; }

   public:
    // Modify the CS histo so cross section shoot is within ]HalfOpenAngleMin,HalfOpenAngleMax[
    void SetCSAngle(double CSHalfOpenAngleMin, double CSHalfOpenAngleMax);

   public: // Kinematics
    // Check that the reaction is alowed
    bool IsLabCrossSection() { return fLabCrossSection; };

    // Use fCrossSectionHist to shoot a Random ThetaCM and set fThetaCM to this value
    double ShootRandomThetaCM();

    // Compute ThetaLab and EnergyLab for product of reaction
    void GenerateEvent(double& ThetaLight, double& PhiLight, double& KineticEnergyLight, double& ThetaLab3,
                       double& KineticEnergyLab3, double& ThetaLab4, double& KineticEnergyLab4);

    void SetParticle3(double EnergyLab, double ThetaLab);

    void SetParticle3(NPL::Particle p) { fParticle3 = p; };

    ClassDef(InelasticBreakup, 0)
  };
} // namespace NPL
#endif
