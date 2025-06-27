/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : January 2009                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This event Generator is used to simulated AlphaDecay ion Source           *
 *  Very usefull to figure out Geometric Efficacity of experimental Set-Up   *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++
#include<limits>

// G4 headers
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
// G4 headers including CLHEP headers
// for generating random numbers
#include "Randomize.hh"

// NPS headers
#include "EventGeneratorAlphaDecay.hh"

// NPL headers
#include "RootOutput.h"
#include "NPNucleus.h"
#include "NPOptionManager.h"
#include "NPFunction.h"
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorAlphaDecay::EventGeneratorAlphaDecay(){
  // NPTool path
  string GlobalPath = getenv("NPTOOL");
  string StandardPath = GlobalPath + "/Inputs/EventGenerator/";

  m_EnergyLow    =  0  ;
  m_EnergyHigh   =  0  ;
  m_x0           =  0  ;
  m_y0           =  0  ;
  m_z0           =  0  ;
  m_SigmaX       =  0  ;
  m_SigmaY       =  0  ;
  m_SigmaZ       =  0  ;
  m_particle     = NULL;
  m_ExcitationEnergy = 0 ;
  m_direction    = 'z' ;
  m_SourceProfile = "Gauss" ;
  m_ActivityBq   = 1.  ;
  m_TimeWindow   = 50.e-9 ;

  m_ParticleStack = ParticleStack::getInstance();

}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorAlphaDecay::~EventGeneratorAlphaDecay(){
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorAlphaDecay::ReadConfiguration(NPL::InputParser parser){
  
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("AlphaDecay");
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// AlphaDecay reaction found " << endl;

  vector<string> token = {"EnergyLow","EnergyHigh","HalfOpenAngleMin","HalfOpenAngleMax","x0","y0","z0","ActivityBq"};
  cout << "block.size() = " << blocks.size() << endl;
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    cout << "[" << i << "]" << endl;
    if(blocks[i]->HasTokenList(token)){

      m_EnergyLow         = blocks[i]->GetDouble("EnergyLow","MeV");        
      m_EnergyHigh        = blocks[i]->GetDouble("EnergyHigh","MeV");
      m_HalfOpenAngleMin  = blocks[i]->GetDouble("HalfOpenAngleMin","deg"); 
      m_HalfOpenAngleMax  = blocks[i]->GetDouble("HalfOpenAngleMax","deg");
      m_x0                = blocks[i]->GetDouble("x0","mm");                
      m_y0                = blocks[i]->GetDouble("y0","mm");
      m_z0                = blocks[i]->GetDouble("z0","mm");
      m_SourceProfile     = blocks[i]->GetString("SourceProfile");          
      m_SigmaX            = blocks[i]->GetDouble("SigmaX","mm");            
      m_SigmaY            = blocks[i]->GetDouble("SigmaY","mm");
      m_SigmaZ            = blocks[i]->GetDouble("SigmaZ","mm");
      m_direction         = blocks[i]->GetString("Direction");              
      m_ExcitationEnergy  = blocks[i]->GetDouble("ExcitationEnergy","MeV"); 
      m_ActivityBq        = blocks[i]->GetDouble("ActivityBq","Bq");
      //m_ActivityBq        = blocks[i]->GetDouble("ActivityBq","becquerel"); //FIXME:  does not work
      m_TimeWindow        = blocks[i]->GetDouble("TimeWindow","s");         
      cout << "m_parameters done " << endl;
    }
    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorAlphaDecay::GenerateEvent(G4Event* evt){

  double alpha_t = 0.;
  double step_t  = 1.e-2 / m_ActivityBq; // step size in [s] to have a maximum proba of 1.e-2
  int    step_n  = ((int)(m_TimeWindow / step_t) == 0) ? 1 : (int)(m_TimeWindow / step_t);
  bool   DecayOn = true;

  for(int i=0; i<step_n; i++){
    
    m_particle=NULL;
    if(m_particle==NULL){
      m_particle = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(2, 4,m_ExcitationEnergy);
      if(i>0){
        alpha_t = i*step_t;
        DecayOn = !(RandFlat::shoot() >= 1.e-2);
      }
    }
    if(!DecayOn) continue;

    G4double cos_theta_min   = cos(m_HalfOpenAngleMin);
    G4double cos_theta_max   = cos(m_HalfOpenAngleMax);
    G4double cos_theta       = cos_theta_min + (cos_theta_max - cos_theta_min) * RandFlat::shoot();
    G4double theta           = acos(cos_theta)                                                   ;
    G4double phi             = RandFlat::shoot() * 2 * pi                                        ;
    G4double particle_energy = m_EnergyLow + RandFlat::shoot() * (m_EnergyHigh - m_EnergyLow)    ;

    G4double x0, y0, z0;
    G4double momentum_x, momentum_y, momentum_z;

    if(m_SourceProfile=="Flat"){
      double rand_r     = m_SigmaX*sqrt(RandFlat::shoot());
      double rand_theta = RandFlat::shoot() * 2. * pi;
      x0 = m_x0 + rand_r * cos(rand_theta);
      y0 = m_y0 + rand_r * sin(rand_theta);
      z0 = RandFlat::shoot(m_z0*mm-m_SigmaZ*mm, m_z0*mm+m_SigmaZ*mm);
    }
    else{ // by default gaussian source profile is considered 
      x0 = RandGauss::shoot(m_x0,m_SigmaX);
      y0 = RandGauss::shoot(m_y0,m_SigmaY);
      z0 = RandGauss::shoot(m_z0,m_SigmaZ);
    }
    if(m_direction == 'z')
    {
      // Direction of particle, energy and laboratory angle
      momentum_x = sin(theta) * cos(phi)  ;
      momentum_y = sin(theta) * sin(phi)  ;
      momentum_z = cos(theta)             ;

    }
    else if(m_direction == 'y')
    {
      // Direction of particle, energy and laboratory angle
      momentum_z = sin(theta) * cos(phi)  ;
      momentum_x = sin(theta) * sin(phi)  ;
      momentum_y = cos(theta)             ;
    }
    else // = 'x'
    {
      // Direction of particle, energy and laboratory angle
      momentum_y = sin(theta) * cos(phi)  ;
      momentum_z = sin(theta) * sin(phi)  ;
      momentum_x = cos(theta)             ;
    }
    //cout << "add an alpha particle to stack with global time = " << alpha_t << " s ==> " << alpha_t*1.e9 << " ns" << endl;
    NPS::Particle particle(m_particle, theta, particle_energy, G4ThreeVector(momentum_x, momentum_y, momentum_z), G4ThreeVector(x0, y0, z0),alpha_t*1.e9,true);
    m_ParticleStack->AddParticleToStack(particle);
  }

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorAlphaDecay::InitializeRootOutput(){

}
