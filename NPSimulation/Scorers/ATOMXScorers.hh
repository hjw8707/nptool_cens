#ifndef ATOMXScorers_h
#define ATOMXScorers_h 1
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
 *  File old the scorer to record Hit energy,time and position               *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include <map>

#include "G4RunManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "NPImage.h"
using namespace std;
using namespace CLHEP;

namespace ATOMXScorers {

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class ATOMXData {
   public:
    ATOMXData() { m_Index = 0; };
    ATOMXData(const unsigned int& Index, const double& Energy, const double& Time, const double& PositionX,
                    const double& PositionY, const double& PositionZ) {
        m_Index = Index;
        m_Energy = Energy;
        m_Time = Time;
        m_PositionX = PositionX;
        m_PositionY = PositionY;
        m_PositionZ = PositionZ;
    }

    ATOMXData(const unsigned int& Index, const double& Energy, const double& Time, const double& PositionX,
                    const double& PositionY, const double& PositionZ, const int& PDG, const double& KineticEnergy) {
        m_Index = Index;
        m_Energy = Energy;
        m_Time = Time;
        m_PositionX = PositionX;
        m_PositionY = PositionY;
        m_PositionZ = PositionZ;
        m_PDG = PDG;
        m_KineticEnergy = KineticEnergy;
    }

    ~ATOMXData() {};

   private:
    unsigned int m_Index;
    double m_Energy;
    double m_KineticEnergy;
    double m_Time;
    double m_PositionX;
    double m_PositionY;
    double m_PositionZ;
    int m_PDG;

   public:
    unsigned int GetIndex() const { return m_Index; };
    double GetEnergy() const { return m_Energy; };
    double GetKineticEnergy() const { return m_KineticEnergy; };
    double GetTime() const { return m_Time; };
    double GetPositionX() const { return m_PositionX; };
    double GetPositionY() const { return m_PositionY; };
    double GetPositionZ() const { return m_PositionZ; };
    int GetPDG() const { return m_PDG; };

   public:
    void Set(const unsigned int& Index, const double& Energy, const double& Time, const double& PositionX,
             const double& PositionY, const double& PositionZ) {
        m_Index = Index;
        m_Energy = Energy;
        m_Time = Time;
        m_PositionX = PositionX;
        m_PositionY = PositionY;
        m_PositionZ = PositionZ;
    }

    void Set(const unsigned int& Index, const double& Energy, const double& Time, const double& PositionX,
             const double& PositionY, const double& PositionZ, const int& PDG, const double& KineticEnergy) {
        m_Index = Index;
        m_Energy = Energy;
        m_KineticEnergy = KineticEnergy;
        m_Time = Time;
        m_PositionX = PositionX;
        m_PositionY = PositionY;
        m_PositionZ = PositionZ;
        m_PDG = PDG;
    }
    void Add(const double& Energy) { m_Energy += Energy; };
    unsigned int GetIndex() { return m_Index; };
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Manage a vector of ATOMX hit
class ATOMXDataVector {
   public:
    ATOMXDataVector() {};
    ~ATOMXDataVector() {};

   private:
    vector<ATOMXData> m_Data;

   public:
    vector<ATOMXData>::iterator find(const unsigned int& index);
    void clear() { m_Data.clear(); };
    vector<ATOMXData>::iterator end() { return m_Data.end(); };
    vector<ATOMXData>::iterator begin() { return m_Data.begin(); };
    unsigned int size() { return m_Data.size(); };
    void Add(const unsigned int& index, const double& Energy) { find(index)->Add(Energy); };

    void Set(const unsigned int& index, const double& Energy, const double& Time, const double& PositionX,
             const double& PositionY, const double& PositionZ) {
        m_Data.push_back(ATOMXData(index, Energy, Time, PositionX, PositionY, PositionZ));
    };

    void Set(const unsigned int& index, const double& Energy, const double& Time, const double& PositionX,
             const double& PositionY, const double& PositionZ, const int& PDG, const double& KineticEnergy) {
        m_Data.push_back(ATOMXData(index, Energy, Time, PositionX, PositionY, PositionZ, PDG, KineticEnergy));
    };

    ATOMXData* operator[](const unsigned int& i) { return &m_Data[i]; };
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class PS_ATOMX : public G4VPrimitiveScorer {
   public:  // with description
    PS_ATOMX(G4String name, G4int depth = 0);
    ~PS_ATOMX() {};

   protected:  // with description
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);

   public:
    void Initialize(G4HCofThisEvent*);
    void EndOfEvent(G4HCofThisEvent*);
    void clear();
    void DrawAll() {};
    void PrintAll() {};

    // Level at which to find the copy number linked to the detector number
    G4int m_Level;

   private:
    ATOMXDataVector m_DataVector;
    G4ThreeVector t_Position;
   public:
    inline unsigned int GetMult() { return m_DataVector.size(); };
    inline double GetEnergy(const unsigned int& i) { return m_DataVector[i]->GetEnergy(); };
    inline double GetTime(const unsigned int& i) { return m_DataVector[i]->GetTime(); };
    inline double GetKineticEnergy(const unsigned int& i) { return m_DataVector[i]->GetKineticEnergy(); };
    inline double GetPositionX(const unsigned int& i) { return m_DataVector[i]->GetPositionX(); };
    inline double GetPositionY(const unsigned int& i) { return m_DataVector[i]->GetPositionY(); };
    inline double GetPositionZ(const unsigned int& i) { return m_DataVector[i]->GetPositionZ(); };
    inline double GetPDG(const unsigned int& i) { return m_DataVector[i]->GetPDG(); };
};

}  // namespace ATOMXScorers

#endif
