#ifndef __SuperX3Physics__
#define __SuperX3Physics__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : january 2011                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *     This class holds the physics class for the SuperX3 detector from Micron    *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// C++ headers
#include <map>
#include <vector>
using namespace std;

// ROOT headers
#include "TCanvas.h"
#include "TH1.h"
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"

// NPTool headers
#include "NPCalibrationManager.h"
#include "NPInputParser.h"
#include "NPVDetector.h"
#include "TSuperX3Data.h"
#include "TSuperX3Spectra.h"

// forward declaration
class TSuperX3Spectra;

class TSuperX3Physics : public TObject, public NPL::VDetector {
 public: //   Constructor and Destructor
  TSuperX3Physics();
  ~TSuperX3Physics(){};

 public:
  void Clear();
  void Clear(const Option_t*){};

 private: // data obtained after BuildPhysicalEvent() and stored in ROOT output file
  vector<Int_t> fEventType;
  vector<Int_t> fDetectorNumber;
  vector<Double_t> fFrontEnergy;
  vector<Double_t> fBackEnergy;
  vector<Double_t> fHalfEnergy;
  vector<Double_t> fFrontTime;
  vector<Double_t> fBackTime;
  vector<Int_t> fFrontStrip;
  vector<Int_t> fBackStrip;

 public:
  // setters
  void SetEventType(Int_t evtType) { fEventType.push_back(evtType); }
  void SetDetectorNumber(Int_t moduleNbr) { fDetectorNumber.push_back(moduleNbr); }
  void SetFrontEnergy(Double_t ener) { fFrontEnergy.push_back(ener); }
  void SetBackEnergy(Double_t ener) { fBackEnergy.push_back(ener); }
  void SethalfEnergy(Double_t ener) { fHalfEnergy.push_back(ener); }
  void SetFrontTime(Double_t time) { fFrontTime.push_back(time); }
  void SetBackTime(Double_t time) { fBackTime.push_back(time); }
  void SetFrontStrip(Int_t x) { fFrontStrip.push_back(x); }
  void SetBackStrip(Int_t y) { fBackStrip.push_back(y); }

  // getters
  Int_t GetEventMultiplicity() { return fFrontEnergy.size(); }
  Int_t GetEventType(Int_t i) { return fEventType[i]; }
  Int_t GetDetectorNumber(Int_t i) { return fDetectorNumber[i]; }
  Double_t GetFrontEnergy(Int_t i) { return fFrontEnergy[i]; }
  Double_t GetBackEnergy(Int_t i) { return fBackEnergy[i]; }
  Double_t GetHalfEnergy(Int_t i) { return fHalfEnergy[i]; }
  Double_t GetFrontTime(Int_t i) { return fFrontTime[i]; }
  Double_t GetBackTime(Int_t i) { return fBackTime[i]; }
  Int_t GetFrontStrip(Int_t i) { return fFrontStrip[i]; }
  Int_t GetBackStrip(Int_t i) { return fBackStrip[i]; }

 public:
  Int_t m_nCounter;     //!
  Bool_t m_Counter[10]; //!

 public: // inherited from VDetector
  // Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
  void ReadConfiguration(NPL::InputParser);

  // Add parameters to the CalibrationManger
  void AddParameterToCalibrationManager();

  // Activate associated branches and link them to the private member object m_EventData
  void InitializeRootInputRaw();

  // Activate associated branches and link them to the private member m_EventPhysics
  void InitializeRootInputPhysics();

  // Create associated branches and associated private member m_EventPhysics
  void InitializeRootOutput();

  // This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract
  // physical parameter.
  void BuildPhysicalEvent();

  // Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but
  // less efficient ...). This method aimed to be used for analysis performed during experiment, when speed is
  // requiered. NB: This method can eventually be the same as BuildPhysicalEvent.
  void BuildSimplePhysicalEvent();

  // Same as above but for online analysis
  void BuildOnlinePhysicalEvent() { BuildPhysicalEvent(); };

  // Clear raw and physics data
  void ClearEventPhysics() { Clear(); }
  void ClearEventData() { m_EventData->Clear(); }

  // Methods related to the TSuperX3Spectra classes
  // Instantiate the TSuperX3Spectra class and the histograms
  void InitSpectra();
  // Fill the spectra defined in TSuperX3Spectra
  void FillSpectra();
  // Used for Online mainly, perform check on the histo and for example change their color if issues are found
  void CheckSpectra();
  // Used for Online only, clear all the spectra hold by the Spectra class
  void ClearSpectra();
  // Write Spectra to file
  void WriteSpectra();

 public: //   Specific to SuperX3
  // Remove bad channel, calibrate the data and apply thresholds
  void PreTreat();

  // Clear the pre treated object
  void ClearPreTreatedData() { m_PreTreatedData->Clear(); }

  // Return false if the channel is disabled by user
  // Frist argument is either "Front" or "Back"
  bool IsValidChannel(string Type, int detector, int channel);

  // Initialize the standard parameters for analysis, i.e.: all channel enable,
  // maximum multiplicity for strip = number of telescope
  void InitializeStandardParameters();

  //   Read the user configuration file; if no file found, load standard one
  void ReadAnalysisConfig();

  // Add detector using cartesian coordinates
  void AddDetector(TVector3 C_X1_Y1, TVector3 C_X16_Y1, TVector3 C_X1_Y16, TVector3 C_X16_Y16);

  // Add detector using spherical coordinates
  void AddDetector(double theta, double phi, double distance, double beta_u, double beta_v, double beta_w);

  // Give an external TSuperX3Data object to TSuperX3Physics. Needed for online analysis for example.
  void SetRawDataPointer(TSuperX3Data* rawDataPointer) { m_EventData = rawDataPointer; }

  // Use for reading Calibration Run, very simple methods; only apply calibration, no condition
  void ReadCalibrationRun(){};

 public: // Methods used for event treatement
  Int_t EventType();
  vector<TVector2> Match_Front_Back();

 private:                          // Data not written in the tree
  TSuperX3Data* m_EventData;       //!
  TSuperX3Data* m_PreTreatedData;  //!
  TSuperX3Physics* m_EventPhysics; //!

 public:
  TSuperX3Data* GetRawData() const { return m_EventData; }
  TSuperX3Data* GetPreTreatedData() const { return m_PreTreatedData; }

 private:                                      // Map of activated Channel
  map<int, vector<bool>> m_FrontChannelStatus; //!
  map<int, vector<bool>> m_BackChannelStatus;  //!

 private: // Parameters used in the analysis
  // If multiplicity is greater than m_MaximumStripMultiplicityAllowed
  // after PreTreat(), event is not treated
  int m_MaximumStripMultiplicityAllowed; //!

  // Tolerance for front / back energy match
  double m_StripEnergyMatchingSigma;         //!
  double m_StripEnergyMatchingNumberOfSigma; //!

  // Energy thresholds
  // Raw Threshold
  int m_FrontE_Raw_Threshold; //!
  int m_BackE_Raw_Threshold;  //!
  // Calibrated Threshold
  double m_FrontE_Calib_Threshold; //!
  double m_BackE_Calib_Threshold;  //!

 private:                                          // Spatial Position of Strip Calculated on bases of detector position
  int m_NumberOfDetectors;                         //!
  vector<vector<vector<double>>> m_StripPositionX; //!
  vector<vector<vector<double>>> m_StripPositionY; //!
  vector<vector<vector<double>>> m_StripPositionZ; //!

 public:
  double GetNumberOfDetectors() { return m_NumberOfDetectors; };
  double GetStripPositionX(int N, int Front, int Back) { return m_StripPositionX[N - 1][Front - 1][Back - 1]; };
  double GetStripPositionY(int N, int Front, int Back) { return m_StripPositionY[N - 1][Front - 1][Back - 1]; };
  double GetStripPositionZ(int N, int Front, int Back) { return m_StripPositionZ[N - 1][Front - 1][Back - 1]; };
  TVector3 GetPositionOfInteraction(int i);
  TVector3 GetDetectorNormal(int i);
  void DumpStrippingScheme(Int_t detecNumber);

 private:               // Geometry and strip number
  double m_SiliconFace; //!     // mm
  int m_NumberOfStrips; //!
  double m_StripPitch;  //!

 private:                     // Spectra Class
  TSuperX3Spectra* m_Spectra; // !

 public: // Spectra Getter
  map<string, TH1*> GetSpectra();

 public: // Static constructor to be passed to the Detector Factory
  static NPL::VDetector* Construct();

  ClassDef(TSuperX3Physics, 1) // TSuperX3Physics
};

namespace SuperX3_LOCAL {
  Double_t fSuperX3_Front_E(TSuperX3Data* EventData, Int_t i);
  Double_t fSuperX3_Front_T(TSuperX3Data* EventData, Int_t i);
  Double_t fSuperX3_Back_E(TSuperX3Data* EventData, Int_t i);
  Double_t fSuperX3_Back_T(TSuperX3Data* EventData, Int_t i);
} // namespace SuperX3_LOCAL

#endif
