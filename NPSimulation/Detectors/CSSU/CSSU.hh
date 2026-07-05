#ifndef CSSU_h
#define CSSU_h 1

#include <string>
#include <map>
#include <vector>

#include "G4LogicalVolume.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4ThreeVector.hh"
#include "NPInputParser.h"
#include "NPSVDetector.hh"
#include "TCSSUData.h"

class G4Material;

class CSSU : public NPS::VDetector {
 public:
  CSSU();
  virtual ~CSSU();

 public:
  void AddScintillatorBar(G4ThreeVector pos, double length, double width, double thickness, double pmtLength,
                          double pmtDiameter, double scintYield, double attenuationLength, double scintillatorRIndex);
  void AddGasBox(G4ThreeVector pos, double sizeX, double sizeY, double sizeZ, double pressureTorr, double stepLimit,
                 int segmentsZ, double segmentGapZ, double energyThreshold, double resoEnergy);

 public:
  void ReadConfiguration(NPL::InputParser);
  void ConstructDetector(G4LogicalVolume* world);
  void InitializeRootOutput();
  void ReadSensitive(const G4Event* event);
  void InitializeScorers();

 private:
  struct ScintillatorBar {
    G4ThreeVector pos;
    double length;
    double width;
    double thickness;
    double pmtLength;
    double pmtDiameter;
    double scintYield;
    double attenuationLength;
    double scintillatorRIndex;
  };

  struct GasBox {
    G4ThreeVector pos;
    double sizeX;
    double sizeY;
    double sizeZ;
    double pressureTorr;
    double stepLimit;
    int segmentsZ;
    double segmentGapZ;
    double energyThreshold;
    double resoEnergy;
  };

  struct GasResponse {
    double energyThreshold;
    double resoEnergy;
  };

 private:
  G4Material* BuildOpticalAir();
  G4Material* BuildOpticalVacuum();   // n=1 vacuum envelope (photons propagate, alpha unaffected)
  G4Material* BuildScintillatorMaterial(const ScintillatorBar& bar, int detectorNumber);
  G4LogicalVolume* BuildMotherVolume(const ScintillatorBar& bar, int detectorNumber);
  G4LogicalVolume* BuildScintillatorBar(const ScintillatorBar& bar, int detectorNumber);
  G4LogicalVolume* BuildPMT(const ScintillatorBar& bar, int detectorNumber);
  G4Material* BuildCF4Gas(const GasBox& gas, int detectorNumber);
  G4LogicalVolume* BuildGasMother(const GasBox& gas, int detectorNumber, G4Material* material);
  G4LogicalVolume* BuildGasSlice(const GasBox& gas, int detectorNumber, double sliceZ, G4Material* material);

 private:
  std::vector<ScintillatorBar> m_Bars;
  std::vector<GasBox> m_GasBoxes;
  std::map<int, int> m_GasSliceToOutputDetector;
  std::map<int, GasResponse> m_GasResponseByOutputDetector;
  TCSSUData* m_Event;

  G4MultiFunctionalDetector* m_CSSUScorer;
  G4MultiFunctionalDetector* m_PMTScorer;

  G4VisAttributes* m_VisMother;
  G4VisAttributes* m_VisScintillator;
  G4VisAttributes* m_VisPMT;
  G4VisAttributes* m_VisGasMother;
  G4VisAttributes* m_VisGas;

 public:
  static NPS::VDetector* Construct();
};

#endif
