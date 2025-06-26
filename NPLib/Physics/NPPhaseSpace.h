#ifndef NPPhaseSpace_h
#define NPPhaseSpace_h
/*****************************************************************************
 * Copyright (C) 2009-2022    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 *                                                                           *
 * Original Author : Adrien Matta contact: matta@lpccaen.in2p3.fr            *
 *                                                                           *
 * Creation Date   : February 2022                                           *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class simulate a multibody phase space base on a given reaction     *
 *  and exit channel                                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
// C++ header
#include <string>

// NPL
#include "NPParticle.h"
#include "NPBeam.h"
#include "NPInputParser.h"
using namespace NPL;

// ROOT header
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TGenPhaseSpace.h"
#include "TVector3.h"
#include "TRandom.h"

using namespace std;

namespace NPL{

  class PhaseSpace{

    public:  // Constructors and Destructors
      PhaseSpace(){};
      ~PhaseSpace(){};

    public:  // Various Method
      Particle GetParticle(string name, NPL::InputParser parser);
      void ReadConfigurationFile(string Path);
      void ReadConfigurationFile(NPL::InputParser);

    private:
      Beam     fBeam;                      //! Beam
      Particle  fTarget;                   //! Target
      std::vector<Particle>  fDaughters;   //! Daughters particle
      std::vector<double>    fExcitations; //! Excitation energy in MeV
      double   fBeamEnergy;                //! Beam energy in MeV
      bool fFermi;                         //! If set to true use the fermi surface weighting of the phase space
      TGenPhaseSpace fPhaseSpace;          //! The phase space object

    public:
      // Getters and Setters
      bool SetBeamLV(const TLorentzVector& LV);
      bool SetBeamLV(double dirX,double dirY,double dirZ,double Kinetic);
      Particle* GetBeam(){return &fBeam;};
      double Generate(){return fPhaseSpace.Generate();}
      TLorentzVector* GetDecayLV(unsigned int& n){return fPhaseSpace.GetDecay(n);}
      size_t GetDecaySize(){return fDaughters.size();}
      Particle* GetParticle(unsigned int& n){return &fDaughters[n];}
      double  GetExcitation(unsigned int& n){return fExcitations[n];}
      
    private: // intern precompute variable
      void Init();
      TLorentzVector fBeamLV;
      TLorentzVector fTargetLV;
      TLorentzVector fTotalLV;

      std::vector<double> fMasses;

      ClassDef(PhaseSpace,0)

  };
}
#endif
