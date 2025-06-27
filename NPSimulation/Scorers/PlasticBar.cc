/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : February 2013                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  File old the scorer specific to the Sharc Detector                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new type of scorer is aim to become the standard for DSSD,SSSD and   *
 * PAD detector (any Silicon Detector)                                       *
 *****************************************************************************/
#include "PlasticBar.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "RootOutput.h"
using namespace PlasticBar;
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
unsigned int PlasticBarData::CalculateIndex(const vector<unsigned int>& level) {

  unsigned int size = level.size();
  unsigned int result = 0;
  unsigned int multiplier = 1;
  for (unsigned int i = 0; i < size; i++) {
    result += level[i] * multiplier;
    multiplier *= 1000;
  }
  return result;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double PS_PlasticBar::EnergyToLight(double Energy, int Z) {
  // Get Energy to Light conversion factor depending on the particle doing the ionisation
  double a0 = t_Convertor[Z][0];
  double a1 = t_Convertor[Z][1];
  double a2 = t_Convertor[Z][2];
  double a3 = t_Convertor[Z][3];
  double Light = a0 * Energy - a1 * (1. - std::exp(-a2 * std::pow(Energy, a3)));
  return (std::max(Light, 0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::Compute_Light() {

  std::vector<PlasticBarData>::iterator it = m_Data.begin();
  std::vector<std::array<double, 7>>::iterator it_Energy = t_Energy_by_ChargeNumber.begin();
  std::vector<std::array<double, 7>>::iterator it_Light = t_Light_by_ChargeNumber.begin();

  for (; it_Energy < t_Energy_by_ChargeNumber.end(); it_Energy++, it_Light++, it++) {
    // cout << "Barre ID : " << (*it).GetIndex() << endl;
    double Light = 0, Energy = 0;
    double ConvertedLight = 0;
    // cout << "   Energy   Light   Z   " << endl;
    for (int i = 0; i < 7; i++) {
      Energy += (*it_Energy)[i];
      if (t_DoConversion) {
        ConvertedLight = EnergyToLight((*it_Energy)[i], i);
        (*it_Light)[i] = ConvertedLight;
        Light += ConvertedLight;
        t_TotalLight_by_Z[i] += ConvertedLight;
      }
      // cout << "   " << Energy;
      // cout << "   " << Light;
      // cout << "   " << i << endl;
    }
    // cout << endl;
    (*it).SetEnergy(Energy);
    if (t_DoConversion)
      (*it).SetLight(Light);
    else
      (*it).SetLight(Energy);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
vector<PlasticBarData>::iterator PlasticBarDataVector::find(const unsigned int& index) {
  for (vector<PlasticBarData>::iterator it = m_Data.begin(); it != m_Data.end(); it++) {
    if ((*it).GetIndex() == index)
      return it;
  }
  return m_Data.end();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_PlasticBar::PS_PlasticBar(G4String name, vector<G4int> NestingLevel, G4int depth) : G4VPrimitiveScorer(name, depth) {
  m_NestingLevel = NestingLevel;
  auto tree = RootOutput::getInstance()->GetTree();
  tree->Branch("PlasticBar_Time", &t_Time);
  tree->Branch("PlasticBar_PositionX", &t_PositionX);
  tree->Branch("PlasticBar_PositionY", &t_PositionY);
  tree->Branch("PlasticBar_TotalEnergy_by_Z", &t_TotalEnergy_by_Z);
  tree->Branch("PlasticBar_TotalLight_by_Z", &t_TotalLight_by_Z);
  t_TotalEnergy_by_Z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  t_TotalLight_by_Z = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
PS_PlasticBar::~PS_PlasticBar() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool PS_PlasticBar::ProcessHits(G4Step* aStep, G4TouchableHistory*) {

  // Contain Energy, Time + as many copy number as nested volume
  unsigned int mysize = m_NestingLevel.size();

  // first lets consider a step going inside the bar, the interaction point can be associated to the position of the end
  // of the step (which helps avoiding having all the steps at the entrance of the bar)
  G4String processname = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  if (processname != "Transportation") {
    G4String particlename = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
    int Z = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();

    // under the energy cut, neutrons can deposit energy via hadIon. But the conversion cant be done since we dont know
    // the recoil particle. We chose to neglect their effect. (validated because deposited energy is less than 1MeV,
    // which is the PM threshold)
    if (particlename != "neutron" && Z < 7) {
      t_Energy = aStep->GetTotalEnergyDeposit();
      t_Time = aStep->GetPostStepPoint()->GetGlobalTime();
      t_PositionX = aStep->GetPostStepPoint()->GetPosition()[0];
      t_PositionY = aStep->GetPostStepPoint()->GetPosition()[1];
      t_Level.clear();

      for (unsigned int i = 0; i < mysize; i++) {
        t_Level.push_back(aStep->GetPostStepPoint()->GetTouchableHandle()->GetCopyNumber(m_NestingLevel[i]));
      }

      if (t_Time < 150) {
        // Check if the bar has been hit before, if not, set the first interaction point.
        std::vector<PlasticBarData>::iterator it;
        it = m_Data.find(PlasticBarData::CalculateIndex(t_Level));

        if (it == m_Data.end()) {
          m_Data.Set(0.0, 0.0, t_PositionY, t_Time, t_Level, t_PositionX);
	  
          AddEntry();                                                // Adds null entry to t_Energy_by_ChargeNumber
          auto lastEnergyEntry = t_Energy_by_ChargeNumber.end() - 1; // iterator
          (*lastEnergyEntry)[Z] += t_Energy;
          t_TotalEnergy_by_Z[Z] += t_Energy;
        }
        else {
          t_Energy_by_ChargeNumber[it - m_Data.begin()][Z] += t_Energy;
          t_TotalEnergy_by_Z[Z] += t_Energy;
        }
      }
    }
    // if(it->GetLight()>160){
    //   cout << "{!!!}  Deposited Light : " << it->GetLight() << endl;
    // }
    // it->Add_Light(t_Light);
    // it->Set_Light(EnergyToLight(it->GetEnergy(),Z));
    // cout << particlename << " of Z = " << Z << " has left :" << endl;
    // cout << t_Energy << " MeV" << endl;
    // cout << EnergyToLight(it->GetEnergy(),Z)*MeV << " MeVEE" << endl;
    // now, if the next interaction is a Transportation, and there is an energy deposition, it means that we are
    // observing the last step. To avoid having all the positions on the exit of the bar, lets consider the
    // PreStepPosition instead.
  } //end "if process transportation"
  return TRUE;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::Initialize(G4HCofThisEvent*) { clear(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::EndOfEvent(G4HCofThisEvent*) { Compute_Light(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::clear() {
  t_Energy_by_ChargeNumber.clear();
  t_Light_by_ChargeNumber.clear();
  m_Data.clear();
  t_Level.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PS_PlasticBar::PrintAll() {}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
