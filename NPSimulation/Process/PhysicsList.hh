/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    : October 2023                                             *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Modular Physics list calling Geant4 reference list                       *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <map>
#include <string>

#include "G4EmConfigurator.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserPhysicsList.hh"
#include "globals.hh"

// Particle definition
#include "G4IonConstructor.hh"

// Decay
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

// EM
#include "G4EmExtraPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmParameters.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4StoppingPhysics.hh"

// Optical
#include "G4OpticalPhysics.hh"

// Drift Electron
#include "G4DriftElectronPhysics.hh"

// Hadronique
#include "G4HadronElasticPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonElasticPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#if NPS_GEANT4_VERSION_MAJOR > 9
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "INCLXXPhysicsListHelper.hh"
#endif
#include "G4HadronPhysicsQGSP_BERT.hh"

class G4VPhysicsConstructor;

class PhysicsList : public G4VModularPhysicsList {
   public:
    PhysicsList();
    virtual ~PhysicsList();

    void ReadConfiguration(std::string filename);
    void ConstructParticle();
    void SetCuts();
    void ConstructProcess();
    void AddStepMax();
    void AddParametrisation();
    void BiasCrossSectionByFactor(double factor);
    void AddIonGasModels();
    void AddPAIModel(const G4String& modname);
    void NewPAIModel(const G4ParticleDefinition* part, const G4String& modname, const G4String& procname);
    void AddPackage(const G4String& name);
    void AddLevelData();

   private:
    std::map<std::string, G4VPhysicsConstructor*> m_PhysList;

   private:  // Physics List
    G4OpticalPhysics* opticalPhysicsList;
    G4DriftElectronPhysics* driftElectronPhysicsList;
    G4VPhysicsConstructor* emPhysicsList;
    G4VPhysicsConstructor* decay_List;
    G4VPhysicsConstructor* radioactiveDecay_List;
    G4EmConfigurator* emConfig;
    G4EmParameters* em_parameters;

   private:  // Physics option
    std::string m_EmList;
    double m_IonBinaryCascadePhysics;
    double m_NPIonInelasticPhysics;
    double m_EmExtraPhysics;
    double m_HadronElasticPhysics;
    double m_HadronInelasticPhysics;
    double m_StoppingPhysics;
    double m_OpticalPhysics;
    double m_DriftElectronPhysics;
    double m_HadronPhysicsQGSP_BIC_HP;
    double m_HadronPhysicsQGSP_BERT_HP;
    double m_HadronPhysicsINCLXX;
    bool m_INCLXXPhysics;
    double m_HadronPhysicsQGSP_INCLXX_HP;
    double m_HadronPhysicsQGSP_INCLXX;
    double m_HadronPhysicsFTFP_INCLXX_HP;
    double m_HadronPhysicsFTFP_INCLXX;
    double m_Decay;
    double m_IonGasModels;
    double m_pai;
    double m_pai_photon;
    double m_Menate_R;
    double m_NeutronHP;
    double m_Shielding;
    double m_ShieldingLEND;
    double m_Shielding_HPT;

    std::map<int, std::string> m_LevelDataFiles;
};

#endif
