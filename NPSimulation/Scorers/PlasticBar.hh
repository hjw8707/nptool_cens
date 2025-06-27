#ifndef PlasticBar_h
#define PlasticBar_h 1
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
 *  File old the scorer specific to the Silicon Detector                     *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 * This new style of scorer is aim to become the standard way of doing scorer*
 * in NPTool.                                                                *
 *The index is build using the TrackID, Detector Number and Strip Number.    *
 *The scorer Hold Energy and time together                                   *
 *Only one scorer is needed for a detector                                   *
 *****************************************************************************/
#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"
// #include "NPSecondaries.hh"

#include <array>
#include <map>
using namespace CLHEP;

namespace PlasticBar {

  // Hold One hit info
  class PlasticBarData {
   public:  
    PlasticBarData(const double& Energy, const double& Light, const double& PositionY, const double& Time,
                   const std::vector<unsigned int>& Nesting, const double& PositionX)
    {
      m_Index = CalculateIndex(Nesting);
      m_Level = Nesting;
      m_Energy = Energy;
      m_Light = Light;
      m_PositionX = PositionX;
      m_PositionY = PositionY;
      m_Time = Time;
    };
    ~PlasticBarData(){};

   private:
    unsigned int m_Index;
    std::vector<unsigned int> m_Level;
    double m_Energy;
    double m_Light;
    double m_PositionX;
    double m_PositionY;
    double m_Time;

   public:
    static unsigned int CalculateIndex(const std::vector<unsigned int>& Nesting);

   public:
    inline unsigned int GetIndex() const { return m_Index; };
    inline std::vector<unsigned int> GetLevel() const { return m_Level; };
    inline double GetEnergy() const { return m_Energy; };
    inline double GetLight() const { return m_Light; };
    inline double GetPositionX() const {return m_PositionX; };
    inline double GetPositionY() const { return m_PositionY; };
    inline double GetTime() const { return m_Time; };

   public:
    inline void SetEnergy(const double& Energy) { m_Energy = Energy; };
    inline void SetLight(const double& Light) { m_Light = Light; };
  };

  // Manage a vector of DSSD hit
  class PlasticBarDataVector {
   public:
    PlasticBarDataVector(){};
    ~PlasticBarDataVector(){};

   private:
    std::vector<PlasticBarData> m_Data;

   public:
    std::vector<PlasticBarData>::iterator find(const unsigned int& index);
    inline void clear() { m_Data.clear(); };
    inline std::vector<PlasticBarData>::iterator end() { return m_Data.end(); };
    inline std::vector<PlasticBarData>::iterator begin() { return m_Data.begin(); };
    inline unsigned int size() { return m_Data.size(); };
    PlasticBarData* operator[](const unsigned int& i) { return &m_Data[i]; };

    inline void Set(const double& Energy, const double& Light, const double& PositionY, const double& Time,
                    const std::vector<unsigned int>& Nesting, const double& PositionX) {
      m_Data.push_back(PlasticBarData(Energy, Light, PositionY, Time, Nesting, PositionX));
    };

  };

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  class PS_PlasticBar : public G4VPrimitiveScorer {

   public: // with description
    PS_PlasticBar(G4String name, std::vector<G4int> NestingLevel, G4int depth = 0);
    ~PS_PlasticBar();

   protected: // with description
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    double EnergyToLight(double Energy, int Z);
    void Compute_Light();

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
    std::vector<G4int> m_NestingLevel;

   private:
    PlasticBarDataVector m_Data;

   private:
    std::vector<std::array<double, 7>> t_Energy_by_ChargeNumber;
    std::vector<std::array<double, 7>> t_Light_by_ChargeNumber;
    inline void AddEntry() {
      std::array<double, 7> nil_array = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
      t_Energy_by_ChargeNumber.push_back(nil_array);
      t_Light_by_ChargeNumber.push_back(nil_array);
    };
    std::array<double, 7> t_TotalEnergy_by_Z = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    std::array<double, 7> t_TotalLight_by_Z = {{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}};

   private:
    double t_Energy;
    double t_Light;
    double t_PositionX;
    double t_PositionY;
    double t_Time;
    std::vector<unsigned int> t_Level;

   private:
    bool t_DoConversion = 1;

    const double t_Convertor[7][4] = {
        {1.0, 0.0, 0.0, 0.0},      // electron
        {0.90, 7.55, 0.099, 0.74}, // proton
        {0.78, 39.3, 0.022, 0.91}, // alpha
        {0.61, 57.1, 0.012, 0.95}, // 7Li
        {0.44, 45.5, 0.010, 0.92}, // 9Be
        {0.35, 34.5, 0.011, 0.91}, // 11B
        {0.30, 25.6, 0.013, 0.91}, // 12C
    };

   public:
    inline unsigned int GetMult() { return m_Data.size(); };
    inline double GetEnergy(const unsigned int& i) { return m_Data[i]->GetEnergy(); };
    inline double GetLight(const unsigned int& i) { return m_Data[i]->GetLight(); };
    inline double GetPositionX(const unsigned int& i) { return m_Data[i]->GetPositionX(); };
    inline double GetPositionY(const unsigned int& i) { return m_Data[i]->GetPositionY(); };
    inline double GetTime(const unsigned int& i) { return m_Data[i]->GetTime(); };
    inline std::vector<unsigned int> GetLevel(const unsigned int& i) { return m_Data[i]->GetLevel(); };

  };
} // namespace PlasticBar

#endif
