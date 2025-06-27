#ifndef GaseousDetectorScorers_h
#define GaseousDetectorScorers_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey Chatillon, contact address: audrey.chatillon@cea.fr *
 *                                                                           *
 * Creation Date  : January 2025                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment: *** COMPLETELY BASED ON AND COPIED FROM  CalorimeterScorers ***  *
 *          *** from Adrien MATTA, contact adress matta@lpccaen.in2p3.f ***  *
 *****************************************************************************/
#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"

#include <map>
using namespace std;
using namespace CLHEP;

namespace GaseousDetectorScorers {
  
  class GaseousDetectorData{
  // --- ------------------------------------------- --- //
  // --- Hold info of one hit ---------------------- --- //
  // ---  one hit = all the interaction steps of the --- //
  // ---            particle in the  SensitiveVolume --- //
  // --- ------------------------------------------- --- //
      public:
        GaseousDetectorData(const double& Energy,
                            const double& Time,
                            const double& Xpos,
                            const double& Ypos,
                            const double& Zpos,
                            const vector<unsigned int>& Nesting,
                            const vector<string>& name,
                            const vector<int>&    tID,
                            const vector<int>&    A,
                            const vector<int>&    Z,
                            const vector<double>& posZ, 
                            const vector<double>& DE,
                            const vector<double>& Tstep){
          m_Index=CalculateIndex(Nesting);
          m_Level=Nesting;
          m_Energy=Energy;
          m_Time=Time;
          m_Xpos=Xpos;
          m_Ypos=Ypos;
          m_Zpos=Zpos;
          m_ParticleName=name;
          m_ParticleA=A;
          m_ParticleZ=Z;
          m_TrackID=tID;
          m_StepPosZ=posZ;
          m_EnergyLossPerStep=DE;
          m_StepTime=Tstep;
          };
        ~GaseousDetectorData(){};

      private:
        unsigned int         m_Index;
        vector<unsigned int> m_Level;
        double         m_Energy;
        double         m_Time;
        double         m_Xpos;
        double         m_Ypos;
        double         m_Zpos;
        vector<string> m_ParticleName;
        vector<int>    m_ParticleA;
        vector<int>    m_ParticleZ;
        vector<int>    m_TrackID;
        vector<double> m_StepPosZ;
        vector<double> m_EnergyLossPerStep;
        vector<double> m_StepTime;

      public:
        static unsigned int CalculateIndex(const vector<unsigned int>& Nesting);

      public:
        inline unsigned int GetIndex() const {return m_Index;}
        inline vector<unsigned int> GetLevel() const {return m_Level;} 
        inline double GetEnergy() const {return m_Energy;}
        inline double GetTime() const {return m_Time;}
        inline double GetXpos() const {return m_Xpos;}
        inline double GetYpos() const {return m_Ypos;}
        inline double GetZpos() const {return m_Zpos;}
        inline vector<string> GetParticleName() const {return m_ParticleName;}
        inline vector<int>    GetParticleA() const {return m_ParticleA;}
        inline vector<int>    GetParticleZ() const {return m_ParticleZ;}
        inline vector<int>    GetTrackID() const {return m_TrackID;}
        inline vector<double> GetStepPosZ() const {return m_StepPosZ;}
        inline vector<double> GetEnergyLossPerStep() const {return m_EnergyLossPerStep;}
        inline vector<double> GetStepTime() const {return m_StepTime;}


      public:
        // to be used in the PS_GaseousDetector::ProcessHits(G4Step*, G4TouchableHistory*)
        // add the current step energy loss to the total energy of the particle
        void Add(const double& DE) {m_Energy+=DE;}; 
        // add specific information at each step
        inline void SetParticleName(const string& name){m_ParticleName.push_back(name);}
        inline void SetParticleA(const int& mass){m_ParticleA.push_back(mass);}
        inline void SetParticleZ(const int& charge){m_ParticleZ.push_back(charge);}
        inline void SetTrackID(const int& id){m_TrackID.push_back(id);}
        inline void SetStepPosZ(const double& z){m_StepPosZ.push_back(z);}
        inline void SetEnergyLossPerStep(const double& DE){m_EnergyLossPerStep.push_back(DE);}
        inline void SetStepTime(const double& z){m_StepTime.push_back(z);}
  };

  class GaseousDetectorDataVector{
  // --- ----------------------- --- //
  // --- Manage a vector of hits --- //
  // --- ----------------------- --- //
    public:
      GaseousDetectorDataVector(){};
      ~GaseousDetectorDataVector(){};

    private:
      vector<GaseousDetectorData> m_Data;

    public:
      vector<GaseousDetectorData>::iterator find(const unsigned int& index) ;
      inline void clear(){m_Data.clear();} ;
      inline vector<GaseousDetectorData>::iterator end() {return m_Data.end();};
      inline vector<GaseousDetectorData>::iterator begin() {return m_Data.begin();};
      inline unsigned int size() {return m_Data.size();};
      inline void Add(const unsigned int& index,const double& Energy) {find(index)->Add(Energy);};
      inline void Set(const double& Energy, 
                      const double& Time, 
                      const double& Xpos, 
                      const double& Ypos, 
                      const double& Zpos, 
                      const vector<unsigned int>& Nesting, 
                      const vector<string>& name, 
                      const vector<int>& tMass, 
                      const vector<int>& tCharge, 
                      const vector<int>& tID, 
                      const vector<double>& z, 
                      const vector<double>& DE,
                      const vector<double>& t) {
        m_Data.push_back(GaseousDetectorData(Energy,Time,Xpos,Ypos,Zpos,Nesting,name,tMass,tCharge,tID,z,DE,t));
      };
      const GaseousDetectorData* operator[](const unsigned int& i) const {return &m_Data[i];};
  };


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  class PS_GaseousDetector : public G4VPrimitiveScorer{
        
    public: // with description
        PS_GaseousDetector(G4String name, vector<G4int> NestingLevel,G4int depth=0);
        ~PS_GaseousDetector();

        
    protected: // with description
        // ProcessHists is calles at each step of the simulation
        G4bool ProcessHits(G4Step*, G4TouchableHistory*);
        
    public:
        void Initialize(G4HCofThisEvent*);
        void EndOfEvent(G4HCofThisEvent*);
        void clear();
        void DrawAll();
        void PrintAll();
        
    private: // How much level of volume nesting should be considered
        // Give the list of the nesting level at which the copy number should be return.
        // 0 is the lowest level possible (the actual volume copy number in which the interaction happen)
        vector<G4int> m_NestingLevel;
        
    private: 
        GaseousDetectorDataVector m_Data;
        double t_Energy; // t_Energy is incremented at each step
        double t_Time;   // time at first step
        double t_Xpos;   // position of the vertex
        double t_Ypos;   // position of the vertex
        double t_Zpos;   // position of the vertex
        // all these vectors are used to accumulate info at each step
        vector<string> t_ParticleName;
        vector<int>    t_ParticleA;
        vector<int>    t_ParticleZ;
        vector<int>    t_TrackID;
        vector<double> t_StepPosZ;
        vector<double> t_EnergyLossPerStep;
        vector<double> t_StepTime;
        // unique vector which is cleared at each step
        vector<unsigned int> t_Level; 
    public:
      inline unsigned int  GetMult() {return m_Data.size();};
      inline double GetEnergy(const unsigned int& i) {return m_Data[i]->GetEnergy();};
      inline double GetTime(const unsigned int& i) {return m_Data[i]->GetTime();};
      inline double GetXpos(const unsigned int& i) {return m_Data[i]->GetXpos();};
      inline double GetYpos(const unsigned int& i) {return m_Data[i]->GetYpos();};
      inline double GetZpos(const unsigned int& i) {return m_Data[i]->GetZpos();};
      inline vector<string> GetParticleName(const unsigned int& i) const {return m_Data[i]->GetParticleName();}
      inline vector<int>    GetParticleA(const unsigned int& i) const {return m_Data[i]->GetParticleA();}
      inline vector<int>    GetParticleZ(const unsigned int& i) const {return m_Data[i]->GetParticleZ();}
      inline vector<int>    GetTrackID(const unsigned int& i) const {return m_Data[i]->GetTrackID();}
      inline vector<double> GetStepPosZ(const unsigned int& i) const {return m_Data[i]->GetStepPosZ();}
      inline vector<double> GetEnergyLossPerStep(const unsigned int& i) const {return m_Data[i]->GetEnergyLossPerStep();}
      inline vector<double> GetStepTime(const unsigned int& i) const {return m_Data[i]->GetStepTime();}
      inline vector<unsigned int> GetLevel(const unsigned int& i) {return m_Data[i]->GetLevel();};
    };
}


#endif
