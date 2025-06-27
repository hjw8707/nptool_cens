#ifndef ProcessScorers_h
#define ProcessScorers_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. MORFOUACE contact address: pierre.morfouace2@cea.fr   *
 *                                                                           *
 * Creation Date  : July 2020                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  Scorer specific to the processes occuring in the simulation              *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * The scorer hold the processes name                                        *
 *****************************************************************************/
#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"
//#include "NPSecondaries.hh"

#include <map>
#include <set>
using namespace std;
using namespace CLHEP;

namespace ProcessScorers {
  // Hold One hit info

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  class PS_Process : public G4VPrimitiveScorer{

    public: // with description
      PS_Process(G4String name, G4int depth=0);
      ~PS_Process();


    protected: // with description
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);

    public:
      void Initialize(G4HCofThisEvent*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

    private: 
      // How much level of volume nesting should be considered
      // Give the list of the nesting level at which the copy number should be return.
      // 0 is the lowest level possible (the actual volume copy number in which the interaction happen)

    private: 
      //double t_Time;
      //vector<unsigned int> t_Level;
      vector<int> t_trackID;
      vector<int> t_parentID;
      vector<string> t_particlename;
      vector<string> t_volumename;
      vector<string> t_processname;
      vector<double> t_processtime;
      vector<double> t_gamma_energy;
      vector<double> t_proton_energy;
      vector<double> t_proton_time;
      vector<int> t_FC_process;
      vector<int> t_efficiency;
      vector<double> t_pos_Z;
      vector<double> t_kinetic_energy;
      vector<double> t_energy_deposit;

      std::set<int> HasBeenTracked; //map<ID,HasBeenTrackedStatus>
    public:
      inline unsigned int GetMult() {return t_processname.size();};
      inline int GetTrackID(const unsigned int& i) {return t_trackID[i];};
      inline int GetParentID(const unsigned int& i) {return t_parentID[i];};
      inline string GetProcessName(const unsigned int& i) {return t_processname[i];};
      inline string GetParticleName(const unsigned int& i) {return t_particlename[i];};
      inline string GetVolumeName(const unsigned int& i) {return t_volumename[i];};
      inline double GetProcessTime(const unsigned int& i) {return t_processtime[i];};
      inline double GetEnergyDeposit(const unsigned int& i) {return t_energy_deposit[i];};
      inline vector<double> GetGammaEnergy() {return t_gamma_energy;};
      inline vector<double> GetKinEnergy() {return t_kinetic_energy;};
      inline vector<double> GetProtonEnergy() {return t_proton_energy;};
      inline vector<double> GetProtonTime() {return t_proton_time;};
      inline vector<int> GetFCProcess() {return t_FC_process;};
   //inline double GetTime(const unsigned int& i) {return m_Data[i]->GetTime();};
      //inline vector<unsigned int> GetLevel(const unsigned int& i) {return m_Data[i]->GetLevel();};
  };
}


#endif
