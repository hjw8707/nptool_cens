#ifndef WireScorer_h
#define WireScorer_h 1
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Audrey ANNE contact address: anne@lpccaen.in2p3.fr       *
 *                                                                           *
 * Creation Date  : July 2024                                                *
 * Last update    : January 2025                                             *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *                                                                           *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/
#include "G4VPrimitiveScorer.hh"
#include "NPSHitsMap.hh"
// #include "NPSecondaries.hh"

#include <array>
#include <map>
using namespace CLHEP;

namespace WireScorers {

 
  /////////////////////////////////////////////////////////////////////
  //////////////////////// class WireData ////////////////////////
  /////////////////////////////////////////////////////////////////////
  // Hold One hit info
  class WireData {
   public:  
    WireData(const double& Time, const std::vector<unsigned int>& Nesting, const unsigned int& LayerNumber, const unsigned int& WireNumber, const unsigned int& Edge, const double& DriftLength, const int& DetectorNumber)
    {
      m_Index = CalculateIndex(Nesting);
      m_Level = Nesting;
      m_Time = Time;
      m_LayerNumber = LayerNumber;
      m_WireNumber = WireNumber;
      m_Edge = Edge;
      m_DriftLength = DriftLength;
      m_DetectorNumber = DetectorNumber;
    };
    ~WireData(){};

   private:
    unsigned int m_Index;
    std::vector<unsigned int> m_Level;
    double m_PositionX;
    double m_Time;
    unsigned int m_LayerNumber;
    unsigned int m_WireNumber;
    unsigned int m_Edge;
    double m_DriftLength;
    int m_DetectorNumber;
    
   public:
    static unsigned int CalculateIndex(const std::vector<unsigned int>& Nesting);

   public:
    inline unsigned int GetIndex() const { return m_Index; };
    inline std::vector<unsigned int> GetLevel() const { return m_Level; };
    inline double GetTime() const { return m_Time; };
    inline unsigned int GetLayerNumber() const { return m_LayerNumber; };
    inline unsigned int GetWireNumber() const { return m_WireNumber; };
    inline unsigned int GetEdge() const { return m_Edge; };
    inline double GetDriftLength() const { return m_DriftLength; };
    inline unsigned int GetDetectorNumber() const {return m_DetectorNumber; };
    
  };
  
  /////////////////////////////////////////////////////////////////////
  //////////////////////// class WireDataVector////////////////////////
  /////////////////////////////////////////////////////////////////////
  
  // Manage a vector of Wire hit
  class WireDataVector {
   public:
    WireDataVector(){};
    ~WireDataVector(){};

   private:
    std::vector<WireData> m_Data;

   public:
    std::vector<WireData>::iterator find(const unsigned int& index);
    inline void clear() { m_Data.clear(); };
    inline std::vector<WireData>::iterator end() { return m_Data.end(); };
    inline std::vector<WireData>::iterator begin() { return m_Data.begin(); };
    inline unsigned int size() { return m_Data.size(); };
    WireData* operator[](const unsigned int& i) { return &m_Data[i]; };

    inline void Set(const double& Time,const std::vector<unsigned int>& Nesting,  const unsigned int& LayerNumber, const unsigned int& WireNumber, const unsigned int& Edge, const double& DriftLength, const unsigned int& DetectorNumber) {
      m_Data.push_back(WireData(Time, Nesting, LayerNumber, WireNumber, Edge, DriftLength, DetectorNumber));
    };

  };
  
  //////////////////////////////////////////////////////
  //////////////// class PS_Wire ///////////////////////
  //////////////////////////////////////////////////////
  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  class PS_Wire : public G4VPrimitiveScorer {

   public: // with description
    PS_Wire(G4String name, std::vector<G4int> NestingLevel, G4int TotalNumberLayer, G4int NumberWireByLayer, G4double DriftSpeed, std::map<unsigned int , double> map_WireAngle,  G4int depth = 0 );
    ~PS_Wire();

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
    std::vector<G4int> m_NestingLevel;
    G4int m_TotalNumberLayer;
    G4int m_NumberWireByLayer;
    G4double m_DriftSpeed;
    std::map<unsigned int , double> m_WireAngle;

   private:
    WireDataVector m_Data;

   private:
    double t_Time;
    std::vector<unsigned int> t_Level;
    unsigned int t_LayerNumber;
    unsigned int t_WireNumber;
    unsigned int t_Edge;
    double t_DrifLength;
    unsigned int t_DetectorNumber;

    std::vector<double> t_PosInWireX;
    std::vector<double> t_PosInWireY;
    std::vector<double> t_PosInWireZ;
    
   public:
    inline unsigned int GetMult() { return m_Data.size(); };
    inline double GetTime(const unsigned int& i) { return m_Data[i]->GetTime(); };
    inline std::vector<unsigned int> GetLevel(const unsigned int& i) { return m_Data[i]->GetLevel(); };

    inline unsigned int GetLayerNumber (const unsigned int& i) {return m_Data[i]->GetLayerNumber(); };
    inline unsigned int GetWireNumber(const unsigned int& i) { return m_Data[i]->GetWireNumber(); };
    inline unsigned int GetEdge(const unsigned int& i) { return m_Data[i]->GetEdge(); };

    inline double GetDriftLength(const unsigned int& i ) {return m_Data[i]->GetDriftLength();};

    inline double GetDetectorNumber(const unsigned int& i) {return m_Data[i]->GetDetectorNumber();};
  };

  
} // namespace Wire

#endif
