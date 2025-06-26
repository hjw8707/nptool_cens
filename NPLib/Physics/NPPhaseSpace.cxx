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

#include"NPPhaseSpace.h"
#include"NPOptionManager.h"
#include "NPCore.h"

#include<fstream>
#include<iostream>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void NPL::PhaseSpace::ReadConfigurationFile(string Path){
  ifstream ReactionFile;
  string GlobalPath = getenv("NPTOOL");
  string StandardPath = GlobalPath + "/Inputs/EventGenerator/" + Path;
  ReactionFile.open(Path.c_str());
  if (!ReactionFile.is_open()) {
    ReactionFile.open(StandardPath.c_str());
    if(ReactionFile.is_open()) {
      Path = StandardPath;
    }
    else {cout << "Reaction File " << Path << " not found" << endl;exit(1);}
  }
  NPL::InputParser parser(Path);
  ReadConfigurationFile(parser);
}

////////////////////////////////////////////////////////////////////////////////
void NPL::PhaseSpace::ReadConfigurationFile(NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("PhaseSpace");
  if(blocks.size()>0 && NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// Phase space reaction found " << endl;

  vector<string> token = {"Beam","Target","Daughters","ExcitationEnergies","Fermi"};
  fDaughters.clear(); 
  fMasses.clear();
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      int v = NPOptionManager::getInstance()->GetVerboseLevel();
      NPOptionManager::getInstance()->SetVerboseLevel(0);
      fBeam.ReadConfigurationFile(parser);
      NPOptionManager::getInstance()->SetVerboseLevel(v);
      fBeamEnergy= fBeam.GetEnergy();
      fTarget= GetParticle(blocks[i]->GetString("Target"),parser);
      fTargetLV=TLorentzVector(0,0,0,fTarget.Mass());
      vector<string> vDaughters= blocks[i]->GetVectorString("Daughters");
      fExcitations= blocks[i]->GetVectorDouble("ExcitationEnergies","MeV");
      fFermi=blocks[i]->GetInt("Fermi"); 
      for(auto d:vDaughters){
        fDaughters.push_back(GetParticle(d,parser));
        fMasses.push_back(fDaughters.rbegin()->Mass());
      }
     SetBeamLV(0,0,1,fBeamEnergy);
    }

    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }

  }
  }
////////////////////////////////////////////////////////////////////////////////
Particle NPL::PhaseSpace::GetParticle(string name, NPL::InputParser parser){
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithTokenAndValue("DefineParticle",name);
  unsigned int size = blocks.size();
  if(size==0)
    return NPL::Particle(name);
  else if(size==1){
    cout << " -- User defined nucleus " << name << " -- " << endl;
    vector<string> token = {"SubPart","BindingEnergy"};
    if(blocks[0]->HasTokenList(token)){
      NPL::Particle N(name,blocks[0]->GetVectorString("SubPart"),blocks[0]->GetDouble("BindingEnergy","MeV"));
      if(blocks[0]->HasToken("ExcitationEnergy"))
        N.SetExcitationEnergy(blocks[0]->GetDouble("ExcitationEnergy","MeV"));
      if(blocks[0]->HasToken("SpinParity"))
        N.SetSpinParity(blocks[0]->GetString("SpinParity").c_str());
      if(blocks[0]->HasToken("Spin"))
        N.SetSpin(blocks[0]->GetDouble("Spin",""));
      if(blocks[0]->HasToken("Parity"))
        N.SetParity(blocks[0]->GetString("Parity").c_str());
      if(blocks[0]->HasToken("LifeTime"))
        N.SetLifeTime(blocks[0]->GetDouble("LifeTime","s"));

      cout << " -- -- -- -- -- -- -- -- -- -- --" << endl;
      return N;
    }
  }
  else{
    NPL::SendErrorAndExit("NPL::Reaction","Too many nuclei define with the same name");
  }

  return (NPL::Particle());
}
////////////////////////////////////////////////////////////////////////////////
bool NPL::PhaseSpace::SetBeamLV(const TLorentzVector& LV){
    fBeamLV=LV;
    fTotalLV = fBeamLV+fTargetLV;
    static std::string option;
    option="";
    if(fFermi)
      option="Fermi";

    return fPhaseSpace.SetDecay(fTotalLV,fMasses.size(),fMasses.data(),option.c_str());

  }
////////////////////////////////////////////////////////////////////////////////
bool NPL::PhaseSpace::SetBeamLV(double dirX,double dirY,double dirZ,double Kinetic){
     fBeamEnergy=Kinetic;
     fBeam.SetKineticEnergy(Kinetic);
     double E = fBeam.GetGamma()*fBeam.Mass();
     double pc=sqrt(E*E-fBeam.Mass()*fBeam.Mass());
     TVector3 p(dirX,dirY,dirZ);
     p=pc*p.Unit();
     TLorentzVector beamLV(p,E);
     return SetBeamLV(beamLV);
  }
