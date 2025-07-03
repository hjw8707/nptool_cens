/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    : October 2015                                             *
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

#include <fstream>
#include <iostream>
// NPS
#include "PhysicsList.hh"
// NPL
#include "NPInputParser.h"
#include "NPNucleus.h"
#include "NPOptionManager.h"

// Local physic directly implemented from the Hadronthrapy example
// Physic dedicated to the ion-ion inelastic processes
#include "Decay.hh"
#include "menate_R.hh"
#ifndef GEANT4_11
#include "NPIonIonInelasticPhysic.hh"
#endif
// G4
#include "G4BetheBlochIonGasModel.hh"
#include "G4BraggIonGasModel.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4IonFluctuations.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4IonPhysics.hh"
#include "G4LevelManager.hh"
#include "G4LossTableManager.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4NuclearLevelData.hh"
#include "G4NuclideTable.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4PhysListFactory.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4RadioactiveDecay.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4RunManager.hh"
#include "G4StepLimiter.hh"
#include "G4StoppingPhysics.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4UniversalFluctuation.hh"
#include "G4VPhysicsConstructor.hh"
#ifdef USE_NEUTRONHP
#include "NeutronHPphysics.hh"
#endif

/////////////////////////////////////////////////////////////////////////////
PhysicsList::PhysicsList() : G4VUserPhysicsList() {
    m_EmList = "Option4";
    defaultCutValue = 1 * mm;  // 0.2*mm;
    opticalPhysicsList = NULL;
    driftElectronPhysicsList = NULL;
    ReadConfiguration("PhysicsListOption.txt");
    G4LossTableManager::Instance();
    SetVerboseLevel(0);

    // ******     Definition of defaults for the physics processes *****
    // ******     in case no physics is called by the macro file   *****
    //
    // The default physics corresponds to the actual QGSP_BIC_HP list
    // but with the following differences:
    // --> G4EmStandardPhysics_option4 for the electromagnetic processes
    //     is used n place of the less accurate G4EmStandardPhysics
    // --> The G4RadioactiveDecayPhysics is add
    // --> G4HadronPhysicsQGSP_BIC is used in place of G4HadronPhysicsQGSP_BIC_HP
    // --> G4HadronElasticPhysics is used in place of G4HadronElasticPhysics_HP

    // Elecromagnetic physics
    // Using the more accurate option4

    emPhysicsList = NULL;
    if (m_EmList == "INCLXX_EM") {
        cout << "Choosing to use INCLXX included EMList (probably standard)" << endl;
    } else if (m_EmList == "Option1") {
        emPhysicsList = new G4EmStandardPhysics_option1();
        cout << "//// Using G4EmStandardPhysics_option1 Physics List ////" << endl;
    } else if (m_EmList == "Option2") {
        emPhysicsList = new G4EmStandardPhysics_option2();
        cout << "//// Using G4EmStandardPhysics_option2 Physics List ////" << endl;
    } else if (m_EmList == "Option3") {
        emPhysicsList = new G4EmStandardPhysics_option3();
        cout << "//// Using G4EmStandardPhysics_option3 Physics List ////" << endl;
    } else if (m_EmList == "Option4") {
        emPhysicsList = new G4EmStandardPhysics_option4();
        cout << "//// Using G4EmStandardPhysics_option4 Physics List ////" << endl;
    } else if (m_EmList == "Standard") {
        emPhysicsList = new G4EmStandardPhysics();
        cout << "//// Using G4EmStandardPhysics default EM constructor Physics List ////" << endl;
    } else if (m_EmList == "Livermore") {
        emPhysicsList = new G4EmLivermorePhysics();
        cout << "//// Using G4EmLivermorePhysics Physics List ////" << endl;
    } else if (m_EmList == "Penelope") {
        emPhysicsList = new G4EmPenelopePhysics();
        cout << "//// Using G4EmPenelopePhysics Physics List ////" << endl;
    } else {
        std::cout << "\r\033[1;31mERROR: User given physics list " << m_EmList
                  << " is not supported, option are Option4 Livermore Penelope\033[0m" << std::endl;
        exit(1);
    }

    RegisterPhysics(emPhysicsList);

    // Hadronic physics
    if (m_IonBinaryCascadePhysics) {
        RegisterPhysics(new G4IonBinaryCascadePhysics());
        cout << "//// Using G4IonBinaryCascadePhysics Physics List ////" << endl;
    }
    if (m_EmExtraPhysics) {
        cout << "//// Using G4IonPhysics and G4EmExtraPhysics Physics List ////" << endl;
        // G4ChargeExchangePhysics* ceph = new G4ChargeExchangePhysics();
        // ceph->SetCrossSectionFactor(10);
        RegisterPhysics(new G4IonPhysics());
        RegisterPhysics(new G4EmExtraPhysics());
    }
    if (m_HadronElasticPhysics) {
        RegisterPhysics(new G4HadronElasticPhysics());
        RegisterPhysics(new G4IonElasticPhysics());
        cout << "//// Using G4HadronElasticPhysics Physics List ////" << endl;
        cout << "//// Using G4IonElasticPhysics Physics List ////" << endl;
    }

#ifndef GEANT4_11
    if (m_NPIonInelasticPhysics) {
        RegisterPhysics(new NPIonIonInelasticPhysic());
        cout << "//// Using NPIonIonInelasticPhysic Physics List ////" << endl;
    }
#endif

#ifdef USE_NEUTRONHP
    if (m_NeutronHP) {
        // RegisterPhysics(new NeutronHPphysics("neutronHP"));
        m_PhysList["NeutronHPphysics"] = new NeutronHPphysics("neutronHP");
        cout << "//// Using NeutronHPPhysics ////" << endl;
    }
#endif
    if (m_StoppingPhysics) {
        RegisterPhysics(new G4StoppingPhysics());
        cout << "//// Using G4StoppingPhysics Physics List ////" << endl;
    }

    if (m_HadronPhysicsINCLXX) {
        RegisterPhysics(new G4HadronPhysicsINCLXX());
        RegisterPhysics(new G4IonINCLXXPhysics());
        cout << "//// Using INCLXX Physics List ////" << endl;
    }

    if (m_HadronPhysicsQGSP_BERT_HP) {
        RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP());
        cout << "//// Using QGSP_BERT_HP Physics List ////" << endl;
    }

    if (m_HadronPhysicsQGSP_BIC_HP) {
#if NPS_GEANT4_VERSION_MAJOR > 9
        RegisterPhysics(new G4HadronPhysicsQGSP_BIC_HP());
        cout << "//// Using QGSP_BIC_HP Physics List ////" << endl;
#else
        std::cout << "\r\032[1;31m Warning: physics list HadronPhysicsQGSP_BIC_HP requires Geant4 10, process not "
                     "activated \033[0m"
                  << std::endl;
#endif
    }

    if (m_HadronPhysicsQGSP_INCLXX_HP) {
#if NPS_GEANT4_VERSION_MAJOR > 9
        AddPackage("QGSP_INCLXX_HP");
#else
        std::cout << "\r\032[1;31m Warning: physics list HadronPhysicsQGSP_INCLXX_HP requires Geant4 10, process not "
                     "activated \033[0m"
                  << std::endl;
#endif
    }

    // Optical Photon for scintillator simulation
    if (m_OpticalPhysics) {
        cout << "//// Using Optical Photon Physics List ////" << endl;
        RegisterPhysics(new G4OpticalPhysics(0));
        /*opticalPhysicsList->SetMaxNumPhotonsPerStep(100);
        opticalPhysicsList->SetScintillationYieldFactor(0.1);
        opticalPhysicsList->SetTrackSecondariesFirst(kScintillation,true);
        opticalPhysicsList->SetTrackSecondariesFirst(kCerenkov,true);
    */
    }

    // Drift electron for gazeous detector simulation
    if (m_DriftElectronPhysics) {
        cout << "//// Using Drift Electron Physics List ////" << endl;
        driftElectronPhysicsList = new G4DriftElectronPhysics(0);
        driftElectronPhysicsList->SetMaxNumDriftElectronPerStep(1e6);
        RegisterPhysics(driftElectronPhysicsList);
    }
    if (m_IonGasModels) {
        cout << "//// Using Ion Gas Model Physics List ////" << endl;
        AddIonGasModels();
    }

    if (m_pai) {
        cout << "//// Using PAI Model Physics List ////" << endl;
        AddPAIModel("pai");
    }
    if (m_pai_photon) {
        cout << "//// Using PAI PHOTON Model Physics List ////" << endl;
        AddPAIModel("pai_photon");
    }

    // Decay physics
    // Add Radioactive decay
    // Gamma decay of known states
    if (m_Decay) {
        cout << "//// Radioactive decay activated ////" << endl;
        decay_List = new G4DecayPhysics();
        radioactiveDecay_List = new G4RadioactiveDecayPhysics();
        RegisterPhysics(decay_List);
        RegisterPhysics(radioactiveDecay_List);

        //////////////////////////////////////////////////////////
        //// Level data
        //////////////////////////////////////////////////////////
        for (auto& [ZA, file] : m_LevelDataFiles) {
            G4NuclearLevelData::GetInstance()->AddPrivateData(ZA % 1000, int(ZA / 1000), file);
            G4cout << " Level data for " << ZA % 1000 << " " << int(ZA / 1000) << " " << file << G4endl;
        }
        //////////////////////////////////////////////////////////

    } else {
        decay_List = 0;
        radioactiveDecay_List = 0;
    }

    if (m_Shielding) {
#if NPS_GEANT4_VERSION_MAJOR > 9
        AddPackage("Shielding");
#else
        std::cout << "\r\032[1;31m Warning: physics list m_Shielding requires Geant4 10, process not "
                     "activated \033[0m"
                  << std::endl;
#endif
    }
    if (m_ShieldingLEND) {
#if NPS_GEANT4_VERSION_MAJOR > 10
        AddPackage("ShieldingLEND");
#else
        std::cout << "\r\032[1;31m Warning: physics list m_ShieldingLEND requires Geant4 11, process not "
                     "activated \033[0m"
                  << std::endl;
#endif
    }
    if (m_Shielding_HPT) {
#if NPS_GEANT4_VERSION_MAJOR > 10
        AddPackage("Shielding_HPT");
#else
        std::cout << "\r\032[1;31m Warning: physics list m_Shielding_HPT requires Geant4 11, process not "
                     "activated \033[0m"
                  << std::endl;
#endif
    }

    emConfig = G4LossTableManager::Instance()->EmConfigurator();
}

////////////////////////////////////////////////////////////////////////////////
void PhysicsList::ReadConfiguration(std::string filename) {
    m_IonBinaryCascadePhysics = 0;
    m_EmExtraPhysics = 0;
    m_HadronElasticPhysics = 0;
    m_NPIonInelasticPhysics = 0;
    m_StoppingPhysics = 0;
    m_OpticalPhysics = 0;
    m_DriftElectronPhysics = 0;
    m_HadronPhysicsQGSP_BIC_HP = 0;
    m_HadronPhysicsQGSP_BERT_HP = 0;
    m_HadronPhysicsINCLXX = 0;
    m_INCLXXPhysics = 0;
    m_HadronPhysicsQGSP_INCLXX_HP = 0;
    m_HadronPhysicsQGSP_INCLXX = 0;
    m_HadronPhysicsFTFP_INCLXX_HP = 0;
    m_HadronPhysicsFTFP_INCLXX = 0;
    m_Decay = 0;
    m_IonGasModels = 0;
    m_pai = 0;
    m_pai_photon = 0;
    m_Menate_R = 0;
    m_NeutronHP = 0;
    m_Shielding = 0;
    m_ShieldingLEND = 0;
    m_Shielding_HPT = 0;

    std::ifstream file(filename.c_str());
    if (!file.is_open()) {
        std::cout << "WARNING: Could not find Physics List option file " << filename << " | Using default Physics List"
                  << std::endl;
        return;
    }

    std::cout << "Reading Physics list option file " << filename << std::endl;

    std::string line;
    while (std::getline(file, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));  // 왼쪽 공백 제거
        line.erase(line.find_last_not_of(" \t\r\n") + 1);  // 오른쪽 공백 제거
        if (line.empty() || line[0] == '%') continue;      // 빈 줄이거나 주석 줄인 경우 무시
        std::istringstream iss(line);
        std::string word;
        std::vector<std::string> str_values;
        while (iss >> word) str_values.push_back(word);
        if (str_values.size() < 2) continue;  // 이름과 값이 없는 줄 무시

        std::string name = str_values[0];  // 첫 번째 단어를 이름으로 사용
        double value;
        try {
            value = std::stod(str_values[1]);
        } catch (const std::invalid_argument& e) {
            value = 0;
        }

        if (name == "EmPhysicsList")
            m_EmList = str_values[1];
        else if (name == "DefaultCutOff")
            defaultCutValue = value * mm;
        else if (name == "IonBinaryCascadePhysics")
            m_IonBinaryCascadePhysics = value;
        else if (name == "NPIonInelasticPhysics")
            m_NPIonInelasticPhysics = value;
        else if (name == "EmExtraPhysics")
            m_EmExtraPhysics = value;
        else if (name == "HadronElasticPhysics")
            m_HadronElasticPhysics = value;
        else if (name == "StoppingPhysics")
            m_StoppingPhysics = value;
        else if (name == "OpticalPhysics")
            m_OpticalPhysics = value;
        else if (name == "DriftElectronPhysics")
            m_DriftElectronPhysics = value;
        else if (name == "HadronPhysicsQGSP_BIC_HP")
            m_HadronPhysicsQGSP_BIC_HP = value;
        else if (name == "HadronPhysicsQGSP_BERT_HP")
            m_HadronPhysicsQGSP_BERT_HP = value;
        else if (name == "HadronPhysicsINCLXX")
            m_HadronPhysicsINCLXX = value;
        else if (name == "HadronPhysicsQGSP_INCLXX_HP") {
            m_HadronPhysicsQGSP_INCLXX_HP = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "HadronPhysicsQGSP_INCLXX") {
            m_HadronPhysicsQGSP_INCLXX = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "HadronPhysicsFTFP_INCLXX_HP") {
            m_HadronPhysicsFTFP_INCLXX_HP = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "HadronPhysicsFTFP_INCLXX") {
            m_HadronPhysicsFTFP_INCLXX = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "Decay")
            m_Decay = value;
        else if (name == "IonGasModels")
            m_IonGasModels = value;
        else if (name == "pai")
            m_pai = value;
        else if (name == "pai_photon")
            m_pai_photon = value;
        else if (name == "Menate_R")
            m_Menate_R = value;
        else if (name == "NeutronHP")
            m_NeutronHP = value;
        else if (name == "Shielding") {
            m_Shielding = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "Shielding_HPT") {
            m_Shielding_HPT = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "ShieldingLEND") {
            m_ShieldingLEND = value;
            if (value) m_INCLXXPhysics = true;
        } else if (name == "LevelData") {
            int Z = std::stoi(str_values[1]);
            int A = std::stoi(str_values[2]);
            m_LevelDataFiles[Z + 1000 * A] = str_values[3];
        } else
            std::cout << "WARNING: Physics List Token '" << name << "' unknown. Token is ignored." << std::endl;
    }

    // Most special process need decay to be activated
    if ((m_IonBinaryCascadePhysics || m_EmExtraPhysics || m_HadronElasticPhysics || m_NPIonInelasticPhysics ||
         m_StoppingPhysics || m_HadronPhysicsQGSP_BIC_HP || m_HadronPhysicsQGSP_BERT_HP || m_HadronPhysicsINCLXX) &&
        !m_Decay) {
        m_Decay = true;
        std::cout << "Information: Selected process require Decay to be activated." << std::endl;
    }
}
/////////////////////////////////////////////////////////////////////////////
PhysicsList::~PhysicsList() {}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddPackage(const G4String& name) {
    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = factory.GetReferencePhysList(name);
    G4int i = 0;
    const G4VPhysicsConstructor* elem = phys->GetPhysics(i);
    G4VPhysicsConstructor* tmp = const_cast<G4VPhysicsConstructor*>(elem);
    while (elem != 0) {
        RegisterPhysics(tmp);
        elem = phys->GetPhysics(++i);
        tmp = const_cast<G4VPhysicsConstructor*>(elem);
    }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddIonGasModels() {
    G4ParticleTable::G4PTblDicIterator* particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
    particleIterator->reset();
    while ((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4String partname = particle->GetParticleName();
        if (partname == "p" || partname == "proton" || partname == "alpha" || partname == "He3" ||
            partname == "GenericIon") {
            G4BraggIonGasModel* mod1 = new G4BraggIonGasModel();
            G4BetheBlochIonGasModel* mod2 = new G4BetheBlochIonGasModel();
            G4double eth = 2. * MeV * particle->GetPDGMass() / proton_mass_c2;
            emConfig->SetExtraEmModel(partname, "ionIoni", mod1, "", 0.0, eth, new G4IonFluctuations());
            emConfig->SetExtraEmModel(partname, "ionIoni", mod2, "", eth, 100 * TeV, new G4UniversalFluctuation());
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddPAIModel(const G4String& modname) {
    G4ParticleTable::G4PTblDicIterator* particleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
    particleIterator->reset();
    while ((*particleIterator)()) {
        G4ParticleDefinition* particle = particleIterator->value();
        G4String partname = particle->GetParticleName();
        if (partname == "e-" || partname == "e+") {
            NewPAIModel(particle, modname, "eIoni");
        } else if (partname == "mu-" || partname == "mu+") {
            NewPAIModel(particle, modname, "muIoni");
        } else if (partname == "proton" || partname == "pi+" || partname == "pi-") {
            NewPAIModel(particle, modname, "hIoni");
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::NewPAIModel(const G4ParticleDefinition* part, const G4String& modname, const G4String& procname) {
    G4String partname = part->GetParticleName();
    if (modname == "pai") {
        G4PAIModel* pai = new G4PAIModel(part, "PAIModel");
        emConfig->SetExtraEmModel(partname, procname, pai, "GasDetector", 0.0, 100. * TeV, pai);
    }

    else if (modname == "pai_photon") {
        G4PAIPhotModel* pai = new G4PAIPhotModel(part, "PAIPhotModel");
        emConfig->SetExtraEmModel(partname, procname, pai, "GasDetector", 0.0, 100. * TeV, pai);
    }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::ConstructParticle() {
    if (m_INCLXXPhysics) {
        this->G4VModularPhysicsList::ConstructParticle();
    }
    if (m_OpticalPhysics) {
        ((G4VPhysicsConstructor*)opticalPhysicsList)->ConstructParticle();
    }
    if (m_DriftElectronPhysics) {
        ((G4VPhysicsConstructor*)driftElectronPhysicsList)->ConstructParticle();
    }
    if (decay_List) {
        decay_List->ConstructParticle();
        radioactiveDecay_List->ConstructParticle();
    } else {
        // If decay is not activated we have to declare the particle ourself
        G4He3::He3Definition();
        G4Deuteron::DeuteronDefinition();
        G4Triton::TritonDefinition();
        G4Alpha::AlphaDefinition();
        G4Proton::ProtonDefinition();
        G4AntiProton::AntiProtonDefinition();
        G4Neutron::NeutronDefinition();
        G4Electron::ElectronDefinition();
        G4Positron::PositronDefinition();
        G4NeutrinoE::NeutrinoEDefinition();
        G4MuonPlus::MuonPlusDefinition();
        G4MuonMinus::MuonMinusDefinition();
        G4AntiNeutrinoE::AntiNeutrinoEDefinition();
        G4Geantino::GeantinoDefinition();
        G4ChargedGeantino::ChargedGeantinoDefinition();
        G4Gamma::GammaDefinition();
        //  mesons
        G4PionPlus::PionPlusDefinition();
        G4PionMinus::PionMinusDefinition();
        G4PionZero::PionZeroDefinition();
        G4Eta::EtaDefinition();
        G4EtaPrime::EtaPrimeDefinition();
        //  G4RhoZero::RhoZeroDefinition();
        G4KaonPlus::KaonPlusDefinition();
        G4KaonMinus::KaonMinusDefinition();
        G4KaonZero::KaonZeroDefinition();
        G4AntiKaonZero::AntiKaonZeroDefinition();
        G4KaonZeroLong::KaonZeroLongDefinition();
        G4KaonZeroShort::KaonZeroShortDefinition();
        // Ion
        G4IonConstructor ionConstructor;
        ionConstructor.ConstructParticle();
    }
}
/////////////////////////////////////////////////////////////////////////////
void PhysicsList::ConstructProcess() {
    if (m_INCLXXPhysics) {
        // G4VModularPhysicsList::ReplacePhysics(new G4EmStandardPhysics_option4());
        this->G4VModularPhysicsList::ConstructProcess();
    } else {
        // Transportation
        AddTransportation();
        // Electromagnetic physics
        emPhysicsList->ConstructProcess();
        if (opticalPhysicsList) {
            ((G4VPhysicsConstructor*)opticalPhysicsList)->ConstructProcess();
        }
        if (driftElectronPhysicsList) {
            ((G4VPhysicsConstructor*)driftElectronPhysicsList)->ConstructProcess();
        }
        if (decay_List) ((G4DecayPhysics*)decay_List)->ConstructProcess();
        if (radioactiveDecay_List) ((G4RadioactiveDecayPhysics*)radioactiveDecay_List)->ConstructProcess();

        // Hadronic physics
        for (auto it = m_PhysList.begin(), end = m_PhysList.end(); it != end; it++) {
            it->second->ConstructProcess();
        }

        BiasCrossSectionByFactor(m_IonBinaryCascadePhysics);
        em_parameters = G4EmParameters::Instance();
        em_parameters->SetFluo(true);
        em_parameters->SetAuger(true);
#if NPS_GEANT4_VERSION_MAJOR > 10
#if NPS_GEANT4_VERSION_MINOR > 2
        em_parameters->SetDeexActiveRegion("DefaultRegionForTheWorld", true, true, true);
        em_parameters->SetDeexcitationIgnoreCut(true);
#endif
#endif
    }
    AddParametrisation();
    AddLevelData();
    return;
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddStepMax() {}
/////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddParametrisation() {
    G4FastSimulationManagerProcess* BeamReaction = new G4FastSimulationManagerProcess("NPSimulationProcess");

    // For 10.3 and higher
#ifndef theParticleIterator
    G4ParticleTable::G4PTblDicIterator* theParticleIterator = GetParticleIterator();
#endif

    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        std::string name = particle->GetParticleName();
        pmanager->AddDiscreteProcess(BeamReaction);
        // Add a Step limiter to the beam particle.
        // This will be used to limit the step of the beam in the target
        pmanager->AddProcess(new G4StepLimiter, -1, -1, 5);

        // Add menate R to the process manager for neutrons
        // (if selected as a option)
        if (m_Menate_R > 0 && name == "neutron") {
            menate_R* theMenate = new menate_R("menateR_neutron");
            theMenate->SetMeanFreePathCalcMethod("ORIGINAL");
            pmanager->AddDiscreteProcess(theMenate);
            std::cout << "||--------------------------------------------------||" << std::endl;
            std::cout << "         MENATE_R Added to Process Manager            " << std::endl;
            std::cout << "||--------------------------------------------------||" << std::endl;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////
void PhysicsList::AddLevelData() {
    NPL::InputParser parser(NPOptionManager::getInstance()->GetReactionFile());
    vector<NPL::InputBlock*> blocks;
    blocks = parser.GetAllBlocksWithToken("LevelData");
    // Loop over all blocks
    if (blocks.size() > 0 && NPOptionManager::getInstance()->GetVerboseLevel())
        cout << endl << "//// Reading user level data " << endl;

    for (unsigned int i = 0; i < blocks.size(); i++) {
        NPL::Nucleus N(blocks[i]->GetMainValue());
        int Z = N.GetZ();
        int A = N.GetA();
        cout << " Nuclei: " << N.GetName() << endl;
        string level_path = blocks[i]->GetString("Path");
        G4IonTable* ionTable = G4IonTable::GetIonTable();
        G4NuclearLevelData* levelData = G4NuclearLevelData::GetInstance();
        levelData->AddPrivateData(Z, A, level_path);
        const G4LevelManager* levelManager = levelData->GetLevelManager(Z, A);
        G4int Nentries = levelManager->NumberOfTransitions() + 1;
        cout << " Found " << Nentries << " states including g.s." << endl;
        for (G4int j = 1; j < Nentries; j++) {  // Excited states
            cout << " - Level " << j << " energy = " << levelManager->LevelEnergy(j)
                 << " MeV \t lifetime = " << levelManager->LifeTime(j)
                 << " ns \t half-life = " << levelManager->LifeTime(j) * log(2) << " ns" << endl;
            G4ParticleDefinition* excitedState = ionTable->GetIon(Z, A, levelManager->LevelEnergy(j));
            excitedState->SetPDGStable(false);
            excitedState->SetPDGLifeTime(levelManager->LifeTime(j));
        }
    }
}

/////////////////////////////////////////////////////////////////////////////
void PhysicsList::SetCuts() {
    if (verboseLevel > 0) {
        G4cout << "PhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue, "Length") << G4endl;
    }

    // Special Cut for optical photon to be emmitted
    SetCutsWithDefault();

    if (m_NeutronHP) SetCutValue(0 * mm, "proton");
}

////////////////////////////////////////////////////////////////////////////////
/////// Friend Method for CS biasing
////////////////////////////////////////////////////////////////////////////////
void PhysicsList::BiasCrossSectionByFactor(double factor) { factor++; }
