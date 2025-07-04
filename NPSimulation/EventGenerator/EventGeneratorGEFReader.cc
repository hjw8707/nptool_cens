/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Benoit MAUSS  contact address: benoit.mauss@cea.fr       *
 *                                                                           *
 * Creation Date  : October 2024                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This event Generator is used to simulate fission events based on         *
 *  data files generated by the GEF model code (General Description of       *
 *  Fission Observables)                                                     *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++
#include<limits>
#include<iostream>
#include<sstream>

// G4 headers
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
// G4 headers including CLHEP headers
// for generating random numbers
#include "Randomize.hh"

// NPS headers
#include "EventGeneratorGEFReader.hh"

// NPL headers
#include "RootOutput.h"
#include "NPNucleus.h"
#include "NPOptionManager.h"
#include "NPFunction.h"
using namespace CLHEP;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorGEFReader::EventGeneratorGEFReader(){
  m_ParticleStack = ParticleStack::getInstance();
  event_ID=0;

  m_isTwoBody = false;
  HasInputDataFile = false;
  m_FissionConditions = new TFissionConditions();
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorGEFReader::~EventGeneratorGEFReader(){
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
EventGeneratorGEFReader::SourceParameters::SourceParameters(){
  m_x0           =  0  ;
  m_y0           =  0  ;
  m_z0           =  0  ;
  m_SigmaX       =  0  ;
  m_SigmaY       =  0  ;
  m_SigmaZ       =  0  ;
  m_Boost        =  0  ;
  m_direction    = 'z' ;
  m_BeamProfile  = "Gauss" ;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorGEFReader::AttachFissionConditions(){
  if(RootOutput::getInstance()->GetTree()->FindBranch("FissionConditions"))
    RootOutput::getInstance()->GetTree()->SetBranchAddress("FissionConditions", &m_FissionConditions);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorGEFReader::ReadConfiguration(NPL::InputParser parser){

  AttachFissionConditions();
  if(!RootOutput::getInstance()->GetTree()->FindBranch("FissionConditions")){
    RootOutput::getInstance()->GetTree()->Branch("FissionConditions", "TFissionConditions", &m_FissionConditions);
  }
  vector<NPL::InputBlock*> blocks = parser.GetAllBlocksWithToken("GEFReader");
  m_Parameters.reserve(blocks.size());
  if(NPOptionManager::getInstance()->GetVerboseLevel())
    cout << endl << "\033[1;35m//// GEFReader source found " << endl;

  vector<string> token = {"x0","y0","z0","Particle"};
  for(unsigned int i = 0 ; i < blocks.size() ; i++){
    if(blocks[i]->HasTokenList(token)){
      m_Parameters.push_back(SourceParameters());
      const vector<SourceParameters>::reverse_iterator it = m_Parameters.rbegin();

      if(blocks[i]->HasToken("GEFversion"))
      {
        it->m_GEFversion=blocks[i]->GetDouble("GEFversion","");
        if     (2023.11 <= it->m_GEFversion && it->m_GEFversion <= 2023.22) {
          version_shift=0;
          version_ff_shift=0;
        }
        else if(2023.31 <= it->m_GEFversion && it->m_GEFversion <= 2023.33) {
          version_shift=1;
          version_ff_shift=0;
        }
        else if(2024.11 <= it->m_GEFversion && it->m_GEFversion <= 2024.11){
          version_shift=1;
          version_ff_shift=1;
        }
        else {
          cout << "!!!!!!!!!!!!!!!!! ERROR: GEF version not verified !!!!!!!!!!!!!!!!!!" << endl;
          cout << "Check the column positions of your GEF file and" << endl;
          cout << "add the version to the correct shift in EventGeneratorGEFReader.cc" << endl;
          exit(1);
        }
      }
      if(blocks[i]->HasToken("InputDataFile")){
        vector<string> file = blocks[i]->GetVectorString("InputDataFile");
        fInputDataFile.open(file.at(0).c_str());
        HasInputDataFile = true;

        string line;
        int    count=0;

        while(count==0)
        {
          getline(fInputDataFile,line);
          //cout << "line : " << line << endl;
          if (line.find("*")==0) continue;
          
          istringstream iss(line);
          count = std::distance(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>());
          //cout << "line : " << line << endl;
          
          if(count>0) 
          {// Storing the first line to be read by GenerateEvent
            iss.clear();              
            iss.seekg(0, iss.beg);    
            double data;
            while(iss >> data) {
              //cout << "data = " << data << endl;
              LastLine.push_back(data);
            }
          }
        }
      }

      if(blocks[i]->HasToken("Direction")){
        it->m_direction=blocks[i]->GetString("Direction");
      }
      it->m_x0                    = blocks[i]->GetDouble("x0","mm");
      it->m_y0                    = blocks[i]->GetDouble("y0","mm");
      it->m_z0                    = blocks[i]->GetDouble("z0","mm");
      vector<string> particleName = blocks[i]->GetVectorString("Particle");
      for(unsigned int j = 0 ; j < particleName.size() ; j++){
        if(particleName[j]=="proton"){ it->m_particleName.push_back("proton")  ;}
        else if(particleName[j]=="gamma") { it->m_particleName.push_back("gamma") ;}
        else if(particleName[j]=="neutron") {it->m_particleName.push_back("neutron") ;}
        else it->m_particleName.push_back(particleName[j]);
      }
      if(blocks[i]->HasToken("TwoBodyReaction")){
        string my_reaction = blocks[i]->GetString("TwoBodyReaction");
        m_TwoBodyReaction = new NPL::Reaction(my_reaction);
        m_isTwoBody = true;
      }

      //if(blocks[i]->HasToken("KineticEnergy_FS"))
      //it->m_Boost=blocks[i]->GetDouble("KineticEnergy_FS","MeV");
      if(blocks[i]->HasToken("FissioningSystem"))
      {
        it->m_FissioningSystemName=blocks[i]->GetString("FissioningSystem");
        m_FissioningSystem=new NPL::Particle(it->m_FissioningSystemName);
        m_FissioningSystem->SetKineticEnergy(0,0,0);
      }
      if(blocks[i]->HasToken("BeamProfile"))
        it->m_BeamProfile=blocks[i]->GetString("BeamProfile");
      if(blocks[i]->HasToken("SigmaX"))
        it->m_SigmaX=blocks[i]->GetDouble("SigmaX","mm");
      if(blocks[i]->HasToken("SigmaY"))
        it->m_SigmaY=blocks[i]->GetDouble("SigmaY","mm");
      if(blocks[i]->HasToken("SigmaZ"))
        it->m_SigmaZ=blocks[i]->GetDouble("SigmaZ","um");  // attention in um

    }
    else{
      cout << "ERROR: check your input file formatting \033[0m" << endl;
      exit(1);
    }
  }

  for(auto& par : m_Parameters) {
    for(auto partName: par.m_particleName){
      AllowedParticles.push_back(partName);
    }
    cout << "///////// Warning: Only ";
    for(auto particle: AllowedParticles) cout << particle << ", " ;
    cout << "will be simulated" << endl;

    if(AllowedParticles.size()==0) cout << "//////// All particles of the file will be simulated" << endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorGEFReader::GetBoostFromTwoBodyReaction(double Ex){
  double Theta3, E3;
  double Theta4, E4;
  double ThetaCM = RandFlat::shoot() * 2 * pi;
  m_TwoBodyReaction->SetExcitation4(Ex);
  m_TwoBodyReaction->SetThetaCM(ThetaCM);
  m_TwoBodyReaction->KineRelativistic(Theta3,E3,Theta4,E4);

  double Phi4 = RandFlat::shoot() * 2 * pi;
  double Phi3 = Phi4 - pi;
  m_FissioningSystem->SetKineticEnergy(E4,Theta4,Phi4);

  m_TwoBodyReaction->SetParticle3(E3,Theta3,Phi3);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorGEFReader::GenerateEvent(G4Event*){
  m_FissionConditions->Clear();
  for(auto& par : m_Parameters) {
    int fileParticleMultiplicity=-1;

    vector<double>                ELab;
    vector<double>                CosThetaLab;
    vector<double>                PhiLab;
    vector<G4ParticleDefinition*> vParticle;

    vector<int>    ZFrag;
    vector<int>    AFrag;
    vector<double> ELabFrag;
    vector<double> cos_thetaFrag;
    vector<double> PhiFrag;
    vector<int>    NneutronsEmission; // may be used to correct the energy of the fragments
    double         CN_ExcitationEnergy;

    TVector3 LightFragmentDirection(0,0,1.);
    TVector3 HeavyFragmentDirection(0,0,1.);
    TLorentzVector ImpulsionFrag[2];

    NPL::Particle *Neutron=new NPL::Particle("1n");
    NPL::Particle *Proton=new NPL::Particle("1H");
    NPL::Particle *Gamma=new NPL::Particle("gamma");
    NPL::Particle *Fragment;
    NPL::Particle *Ejectile;
    TVector3 FissioningSystemBoost;
    int it=0; 

    int index_Z_CN  =  0 + version_shift;                    // Zsad
    int index_A_CN  =  1 + version_shift;                    // Asad
    int index_Ex    = 24 + version_shift + version_ff_shift; // E@fission
    int data_length = 25 + version_shift + version_ff_shift; // data length of CN and FF, without particle-list info prior to fission
    
    int index_Zf[2];
    int index_Af[2];
    int index_ELabf[2];
    int index_nf[2];
    for(int frag=0; frag<2; frag++){
      index_Zf[frag] = 2 + frag + version_shift + version_ff_shift;
      index_Af[frag] = 6 + frag + version_shift + version_ff_shift;
      index_ELabf[frag] = 12 + 3*frag + version_shift + version_ff_shift;
      index_nf[frag] = 20 + frag + version_shift + version_ff_shift;
    }

    if(HasInputDataFile && !fInputDataFile.eof())
    {
      fileParticleMultiplicity=0;
      while(it==0 || LastLine.size()==1 || (LastLine.size()>0 && LastLine.at(1)<50))
      {
        it++;
	// === Read current line
        vector<double> DataLine;
        DataLine = LastLine;
        LastLine.clear();
        //cout << "line starts with " << DataLine.at(0) << " length : " << DataLine.size() << endl;
        //for (int d=0; d<DataLine.size(); d++){
        //	cout << d << " " << DataLine.at(d) << endl;
        //}

	// === Read data of the next line
        string         line;
        getline(fInputDataFile,line);
        istringstream  iss(line);
        double         data;
        //cout << "    it=" << it << " line=" << line << endl;
        while(iss >> data) LastLine.push_back(data);

	
	// === Extract data of the current line
        // FIRST LINE OF THE FISSION EVENT : 
        // data on CN, FFl, FFh and possibly pre-fission emission
        if(DataLine.size()>=data_length && DataLine.at(index_Z_CN)>50)
        {
          double Ex  = DataLine.at(index_Ex); // Excitation energy at fission
          if(m_isTwoBody){
            GetBoostFromTwoBodyReaction(Ex);
            Ejectile = m_TwoBodyReaction->GetParticle3();
          }
          FissioningSystemBoost = m_FissioningSystem->GetEnergyImpulsion().BoostVector();

          // Compound nucleus (after pre-fission emission)
          int Z_CN = DataLine.at(index_Z_CN);
          int A_CN = DataLine.at(index_A_CN);

          m_FissionConditions->SetZ_CN(Z_CN);
          m_FissionConditions->SetA_CN(A_CN);
          m_FissionConditions->SetEx_CN(Ex);
          m_FissionConditions->SetELab_CN(m_FissioningSystem->GetEnergy());
          m_FissionConditions->SetThetaLab_CN(m_FissioningSystem->GetEnergyImpulsion().Theta());

          // Fragment emission
          for(int frag=0;frag<2;frag++)
          {
            int            Zf = DataLine.at(index_Zf[frag]);
            int            Af = DataLine.at(index_Af[frag]); // Post-scission
            double      ELabf = DataLine.at(index_ELabf[frag]);
            double cos_thetaf = DataLine.at(index_ELabf[frag]+1);
            double       Phif = DataLine.at(index_ELabf[frag]+2) *M_PI/180;
            Fragment=new NPL::Particle(Zf,Af);
            Fragment->SetKineticEnergy(ELabf);
            Fragment->EnergyToEnergyImpulsion(ELabf,acos(cos_thetaf),Phif);

            ImpulsionFrag[frag] = Fragment->GetEnergyImpulsion();
            ImpulsionFrag[frag].Boost(FissioningSystemBoost);

            if(frag==0) LightFragmentDirection.SetMagThetaPhi(1.,acos(cos_thetaf),Phif);

            ELabf      = ImpulsionFrag[frag].E() - Fragment->Mass();
            cos_thetaf = ImpulsionFrag[frag].CosTheta();

            ZFrag            .push_back( Zf );
            AFrag            .push_back( Af );
            ELabFrag         .push_back( ELabf );
            cos_thetaFrag    .push_back( cos_thetaf);
            PhiFrag          .push_back( Phif);
            NneutronsEmission.push_back( DataLine.at(index_nf[frag]));
            //delete Fragment;

            // Filling fission conditions
            m_FissionConditions->SetFragmentZ(Zf);
            m_FissionConditions->SetFragmentA(Af);
            m_FissionConditions->SetFragmentKineticEnergy(ELabf);
            m_FissionConditions->SetFragmentBrho(Fragment->GetBrho());
            m_FissionConditions->SetFragmentTheta(acos(cos_thetaf));
            m_FissionConditions->SetFragmentPhi(Phif);
          } // end of for(frag)
          CN_ExcitationEnergy = DataLine.at(index_Ex);

          // Particle emission before fission
          if(DataLine.size()>data_length)
          {
            string ParticleList = to_string(DataLine.at(data_length));
            int idx=0;
            for(char type: ParticleList) // read one digit of the particle list 
              if((type=='1' || type=='3' || type=='5')
                  && (AllowedParticles.size()==0 || std::find(AllowedParticles.begin(), AllowedParticles.end(), "neutron") != AllowedParticles.end()))
              {// Pre-fission neutron treatment:

                double      ELabn = DataLine.at(data_length + 1 + idx*3);
                double cos_thetan = DataLine.at(data_length + 2 + idx*3);
                double       Phin = RandFlat::shoot() * 2 * pi;
                Neutron->SetKineticEnergy(ELabn,acos(cos_thetan),Phin);

                TLorentzVector ImpulsionNeutron = Neutron->GetEnergyImpulsion();
                ImpulsionNeutron.Boost(FissioningSystemBoost);
                ELabn      = ImpulsionNeutron.E() - Neutron->Mass();
                cos_thetan = ImpulsionNeutron.CosTheta();			  

                ELab       .push_back(ELabn);
                CosThetaLab.push_back(cos_thetan);
                PhiLab     .push_back(Phin);
                vParticle  .push_back(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));

                fileParticleMultiplicity++;
                idx++;
              }
              else if((type=='2' || type=='4')
                  && (AllowedParticles.size()==0 || std::find(AllowedParticles.begin(), AllowedParticles.end(), "proton") != AllowedParticles.end()))
              {// Pre-fission proton treatment:
                double      ELabp = DataLine.at(data_length + idx*3);
                double cos_thetap = DataLine.at(data_length + idx*3);
                double       Phip = RandFlat::shoot() * 2 * pi;
                Proton->SetKineticEnergy(ELabp,acos(cos_thetap),Phip);

                TLorentzVector ImpulsionProton = Proton->GetEnergyImpulsion();
                ImpulsionProton.Boost(FissioningSystemBoost);
                ELabp      = ImpulsionProton.E() - Proton->Mass();
                cos_thetap = ImpulsionProton.CosTheta();			  

                ELab       .push_back(ELabp);
                CosThetaLab.push_back(cos_thetap);
                PhiLab     .push_back(RandFlat::shoot() * 2 * pi);
                vParticle  .push_back(G4ParticleTable::GetParticleTable()->FindParticle("proton"));

                fileParticleMultiplicity++;
                idx++;
              }
              else if(type=='0') break; // if digit is 0, end of particle emission
          }// end of if(DataLine.size()>data_length)
        }// end of if(DataLine.size()>=data_length && DataLine.at(index_Z_CN)>50)

        // ELSE IF THE OPTIONS OF PROMPT-NEUTRON EMISSION BY FF IS SWITCH ON IN GEF INPUT
        // NEUTRONS: first line has 0 at DataLine.(0) and gives the list of the post-scission neutrons  
        else if(DataLine.size()>0 && DataLine.at(0)==0 
            && (AllowedParticles.size()==0 || std::find(AllowedParticles.begin(), AllowedParticles.end(), "neutron") != AllowedParticles.end())){
          for(int it=0;it<(DataLine.size()-1)/3;it++)
          {// Promt fission neutron treatment

            double ELabn      = DataLine.at(1+3*it);
            double cos_thetan = DataLine.at(2+3*it);
            double Phin       = DataLine.at(3+3*it) *M_PI/180.;

            TVector3 NeutronLabAngle(0,0,1.);
            NeutronLabAngle.SetMagThetaPhi(1.,acos(cos_thetan),Phin);
            NeutronLabAngle.RotateUz(LightFragmentDirection);
            Neutron->SetKineticEnergy(ELabn,NeutronLabAngle.Theta(),NeutronLabAngle.Phi());

            TLorentzVector ImpulsionNeutron = Neutron->GetEnergyImpulsion();
            ImpulsionNeutron.Boost(FissioningSystemBoost);
            ELabn      = ImpulsionNeutron.E()-Neutron->Mass();
            cos_thetan = ImpulsionNeutron.CosTheta();
            Phin       = NeutronLabAngle.Phi();

            ELab       .push_back(ELabn);
            CosThetaLab.push_back(cos_thetan);
            PhiLab     .push_back(Phin);
            vParticle  .push_back(G4ParticleTable::GetParticleTable()->FindParticle("neutron"));

            fileParticleMultiplicity++;
          }// end of for(it)
	}// end of else if prompt-neutrons

        // ELSE IF THE OPTIONS OF PROMPT GAMMA EMISSION BY FF IS SWITCH ON IN GEF INPUT
        // gamma from light (DataLine.at(0) in [3..5]) and heavy (DataLine.at(0) in[6..8])
        else if(DataLine.size()>0 && (DataLine.at(0)>=3 && DataLine.at(0)<=8)
            && (AllowedParticles.size()==0 || std::find(AllowedParticles.begin(), AllowedParticles.end(), "gamma") != AllowedParticles.end())){
          for(int it=0;it<(DataLine.size()-1);it++)
          {// Prompt fission gamma treatment

            if((DataLine.at(0)==5 || DataLine.at(0)==8) && (to_string(DataLine.at(it+1))=="GS" || to_string(DataLine.at(it+1)).at(0)=='M'))
            {
              cout << "DataLine.at(it+1)=" << DataLine.at(it+1) << endl;
              break;
            }

            double ELabg      = DataLine.at(it+1);
            double cos_thetag = RandFlat::shoot(-1.,1.);
            double Phig       = RandFlat::shoot() * 2 * pi;

            Gamma->SetKineticEnergy(ELabg,acos(cos_thetag),Phig);
            TLorentzVector ImpulsionGamma = Gamma->GetEnergyImpulsion();
            if(DataLine.at(0)<=5)
              ImpulsionGamma.Boost(ImpulsionFrag[0].BoostVector());
            else if(DataLine.at(0)<=8)
              ImpulsionGamma.Boost(ImpulsionFrag[1].BoostVector());

            ImpulsionGamma.Boost(FissioningSystemBoost);
            ELabg      = ImpulsionGamma.E();
            cos_thetag = ImpulsionGamma.CosTheta();		  

            ELab       .push_back( ELabg);
            CosThetaLab.push_back( cos_thetag); // Need to add the fragment boost to the gamma emission. 
            PhiLab     .push_back( Phig);
            vParticle  .push_back(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));

            fileParticleMultiplicity++;
          }// end of for (it)
	}// end of else if gamma-ray
      } // end of while(it==0 || LastLine.size()==1 || (LastLine.size()>0 && LastLine.at(1)<50))
      
      if(fInputDataFile.eof())
      {
        cout << "Problem with Input data file ? no more data to put in" << endl;
        return;
      }
    }// end of if(HasInputDataFile && !fInputDataFile.eof())
    else
    {
      cout << "Problem with Input data file ? no data to put in" << endl;
      return;
    }



    // === add particles into the stack of the data extracted from the current line 
    G4double x0, y0, z0;
    TVector3 Momentum;
    
    if(par.m_BeamProfile=="Flat"){
      double rand_r     = par.m_SigmaX*sqrt(RandFlat::shoot());
      double rand_theta = RandFlat::shoot() * 2. * pi;
      x0 = par.m_x0 + rand_r * cos(rand_theta);
      y0 = par.m_y0 + rand_r * sin(rand_theta);
      z0 = RandFlat::shoot(par.m_z0*mm-par.m_SigmaZ*mm, par.m_z0*mm+par.m_SigmaZ*mm);
      //cout << "Circular Flat beam profile: x0 = " << x0 << ", y0 = " << y0 << ", z0 = " << z0 << endl;
    }
    else{ // by default gaussian beam profile is considered 
      x0 = RandGauss::shoot(par.m_x0,par.m_SigmaX);
      y0 = RandGauss::shoot(par.m_y0,par.m_SigmaY);
      z0 = RandGauss::shoot(par.m_z0,par.m_SigmaZ);
    }
    
    if(AllowedParticles.size()==0 || std::find(AllowedParticles.begin(), AllowedParticles.end(), "fragments") != AllowedParticles.end())
      for(int i_m=0;i_m<2;i_m++)
	{
	  G4double theta                          = acos(cos_thetaFrag.at(i_m)) ;
	  G4double phi                            = PhiFrag.at(i_m)         ;
	  G4double particle_energy                = ELabFrag.at(i_m)           ;
	  G4ParticleDefinition* GeneratedParticle = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(ZFrag.at(i_m),AFrag.at(i_m),0.);

	  Momentum = ShootParticle(theta,phi,par.m_direction);
	  //cout << "particle_energy=" << particle_energy << endl;
	  NPS::Particle particle(GeneratedParticle, theta,particle_energy,G4ThreeVector(Momentum.x(), Momentum.y(), Momentum.z()),G4ThreeVector(x0, y0, z0));
	  m_ParticleStack->AddParticleToStack(particle);
       }

    for(int i_m=0;i_m<fileParticleMultiplicity;i_m++)
      {
	G4double theta                          = acos(CosThetaLab.at(i_m)) ;
	G4double phi                            = PhiLab.at(i_m)            ;
	G4double particle_energy                = ELab.at(i_m)              ;
	G4ParticleDefinition* GeneratedParticle = vParticle.at(i_m)         ;

	Momentum = ShootParticle(theta,phi,par.m_direction);

	NPS::Particle particle(GeneratedParticle, theta,particle_energy,G4ThreeVector(Momentum.x(), Momentum.y(), Momentum.z()),G4ThreeVector(x0, y0, z0));
	if(particle_energy<=20.0) // NeutronHP crashes for neutrons of higher energies. Need to use a different physics list.
	  m_ParticleStack->AddParticleToStack(particle);
      }
   
    if(m_isTwoBody){
      G4double theta                          = Ejectile->GetImpulsion().Theta();
      G4double phi                            = Ejectile->GetImpulsion().Phi();
      G4double particle_energy                = Ejectile->GetEnergy();
      G4ParticleDefinition* GeneratedParticle = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIon(Ejectile->GetZ(),Ejectile->GetA(),0.);

      Momentum = ShootParticle(theta,phi,par.m_direction);
      //cout << "particle_energy=" << particle_energy << endl;
      NPS::Particle particle(GeneratedParticle, theta,particle_energy,G4ThreeVector(Momentum.x(), Momentum.y(), Momentum.z()),G4ThreeVector(x0, y0, z0));
      m_ParticleStack->AddParticleToStack(particle);
    }
  }// end of for(par)
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void EventGeneratorGEFReader::InitializeRootOutput(){

}

TVector3 EventGeneratorGEFReader::ShootParticle(double theta,double phi,TString direction){

  G4double momentum_x, momentum_y, momentum_z;
  if(direction == 'z')
  {
    // Direction of particle, energy and laboratory angle
    momentum_x = sin(theta) * cos(phi)  ;
    momentum_y = sin(theta) * sin(phi)  ;
    momentum_z = cos(theta)             ;

  }
  else if(direction == 'y')
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

  TVector3 Momentum(momentum_x,momentum_y,momentum_z);

  return Momentum;
}
