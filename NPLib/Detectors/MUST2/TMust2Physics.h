#ifndef TMUST2PHYSICS_H
#define TMUST2PHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : febuary 2009                                             *
 * Last update    : July 2021
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold must2 treated data                                       *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <map>
#include <stdlib.h>
#include <vector>
// NPL
#include "NPCalibrationManager.h"
#include "NPDetectorFactory.h"
#include "NPEnergyLoss.h"
#include "NPInputParser.h"
#include "NPVDetector.h"
#include "NPVTreeReader.h"
#include "TMust2Data.h"
#include "TMust2PhysicsReader.h"
#include "TMust2Spectra.h"
// ROOT
#include "TGraphErrors.h"
#include "TH1.h"
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"

using namespace std;

// Forward Declaration
class TMust2Spectra;

class TMust2Physics : public TObject, public NPL::VDetector, public TMust2PhysicsReader {
 public:
  TMust2Physics();
  ~TMust2Physics();

 public:
  void Clear();
  void Clear(const Option_t*){};

 public:
  
  // Returns an array of good X,Y pairs (matching in energy)
  std::vector<std::pair<unsigned int, unsigned int>> Match_X_Y(); //!
  
  // Returns a pair <CsI,pixel> giving the corresponding pixel matching
  std::pair<unsigned int, unsigned int> GetPixel(std::pair<int,std::pair <unsigned int, unsigned int>> match);//!
  
  bool             Match_Si_CsI(int X, int Y, int CristalNbr, int DetectorNbr); //!
  std::vector<std::pair<int, std::pair<unsigned int, unsigned int>>> 
  Match_Si_CsI(std::vector<std::pair<unsigned int, unsigned int>> Array_of_good_XY_couples);//!
  bool             Match_Si_SiLi(int X, int Y, int PadNbr);
  bool             ResolvePseudoEvent();

 public:
  //   Provide Physical Multiplicity
  // Int_t EventMultiplicity;
  int EventMultiplicity;

  //   Provide a Classification of Event
  vector<int> EventType;

  // Telescope
  vector<int> TelescopeNumber;

  //   Si
  vector<double> Si_E;
  vector<double> Si_E_Raw;
  vector<double> Si_T;
  vector<int> Si_X;
  vector<int> Si_StripNumberX;
  vector<int> Si_Y;
  vector<int> Si_StripNumberY;

  // Use for checking purpose
  vector<double> Si_EX;
  vector<double> Si_TX;
  vector<double> Si_EY;
  vector<double> Si_TY;
  vector<int> TelescopeNumber_X;
  vector<int> TelescopeNumber_Y;

  //   Si(Li)
  vector<double> SiLi_E;
  vector<double> SiLi_T;
  vector<int> SiLi_N;

  //   CsI
  vector<double> CsI_E;
  vector<double> CsI_E_Raw;
  vector<double> CsI_T;
  vector<int> CsI_N;

  // Physical Value
  vector<double> TotalEnergy;

  std::vector<int> Pixel;

 public: //   Innherited from VDetector Class
  //   Read stream at ConfigFile to pick-up parameters of detector
  //   (Position,...) using Token
  void ReadConfiguration(NPL::InputParser parser);

  void ReadDoCalibration(NPL::InputParser parser);

  //   Add Parameter to the CalibrationManger
  void AddParameterToCalibrationManager();

  //   Activated associated Branches and link it to the private member
  //   DetectorData address
  //   In this method mother Branches (Detector) AND daughter leaf
  //   (fDetector_parameter) have to be activated
  void InitializeRootInputRaw();

  //   Activated associated Branches and link it to the private member
  //   DetectorPhysics address
  //   In this method mother Branches (Detector) AND daughter leaf (parameter)
  //   have to be activated
  void InitializeRootInputPhysics();

  //   Create associated branches and associated private member DetectorPhysics
  //   address
  void InitializeRootOutput();


  //   This method is called at each event read from the Input Tree. Aime is to
  //   build treat Raw dat in order to extract physical parameter.
  void BuildPhysicalEvent();

  //   Same as above, but only the simplest event and/or simple method are used
  //   (low multiplicity, faster algorythm but less efficient ...).
  //   This method aimed to be used for analysis performed during experiment,
  //   when speed is requiered.
  //   NB: This method can eventually be the same as BuildPhysicalEvent.
  void BuildSimplePhysicalEvent();

  // Same as above but for online analysis
  void BuildOnlinePhysicalEvent() { BuildPhysicalEvent(); };

  //   Those two method all to clear the Event Physics or Data
  void ClearEventPhysics() { Clear(); }
  void ClearEventData() { m_EventData->Clear(); }
  // Method related to the TSpectra classes, aimed at providing a framework for
  // online applications
  // Instantiate the Spectra class and the histogramm throught it

  void InitSpectra();
  // Fill the spectra hold by the spectra class
  void FillSpectra();
  // Used for Online mainly, perform check on the histo and for example change
  // their color if issues are found
  void CheckSpectra();
  // Used for Online only, clear all the spectra hold by the Spectra class
  void ClearSpectra();

  void SetTreeReader(TTreeReader* TreeReader);

 public: //   Specific to MUST2 Array
  //   Clear The PreTeated object
  void ClearPreTreatedData() { m_PreTreatedData->Clear(); }

  //   Remove bad channel, calibrate the data and apply threshold
  void PreTreat();

  //   Return false if the channel is disabled by user
  //   Frist argument is either 0 for X,1 Y,2 SiLi, 3 CsI
  bool IsValidChannel(const int& DetectorType, const int& telescope, const int& channel);

  //   Initialize the standard parameter for analysis
  //   ie: all channel enable, maximum multiplicity for strip = number of
  //   telescope
  void InitializeStandardParameter();

  //   Read the user configuration file; if no file found, load standard one
  void ReadAnalysisConfig();

  //   Add a Telescope using Corner Coordinate information
  void AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1, TVector3 C_X1_Y128, TVector3 C_X128_Y128);

  //   Add a Telescope using R Theta Phi of Si center information
  void AddTelescope(double theta, double phi, double distance, double beta_u, double beta_v, double beta_w);

  // Use for reading Calibration Run, very simple methods; only apply
  // calibration, no condition
  void ReadCalibrationRun();

  // Give and external TMustData object to TMust2Physics. Needed for online
  // analysis for example.
  void SetRawDataPointer(void* rawDataPointer) { m_EventData = (TMust2Data*)rawDataPointer; }
  // Retrieve raw and pre-treated data
  TMust2Data* GetRawData() const { // std::cout << "test" << std::endl;
    return m_EventData;
  }
  // TMust2Data* GetRawDataRead() const {std::cout << "test getrawdata" << std::endl; return r_EventData; }
  TMust2Data* GetPreTreatedData() const { return m_PreTreatedData; }
  // TMust2Physics* GetPhysicsData() const {return m_EventPhysics;}

  // Use to access the strip position
  double GetStripPositionX(const int N, const int X, const int Y) const {
    return m_StripPositionX[N - 1][X - 1][Y - 1];
  };
  double GetStripPositionY(const int N, const int X, const int Y) const {
    return m_StripPositionY[N - 1][X - 1][Y - 1];
  };
  double GetStripPositionZ(const int N, const int X, const int Y) const {
    return m_StripPositionZ[N - 1][X - 1][Y - 1];
  };

  double GetNumberOfTelescope() const { return m_NumberOfTelescope; };

  // To be called after a build Physical Event
  int GetEventMultiplicity() const { return EventMultiplicity; };

  double GetEnergyDeposit(const int i) const { return TotalEnergy[i]; };

  TVector3 GetPositionOfInteraction(const int i) const;
  TVector3 GetTelescopeNormal(const int i) const;
  bool GetCalPixel(){return Cal_Pixel;};
  unsigned int GetPixelSize(){return PixelSize;};

 private: //   Parameter used in the analysis
  // By default take EX and TY.
  bool m_Take_E_Y; //!
  bool m_Take_T_Y; //!

  // Size Container to be used in the different loop of the analysis (value must
  // be given locally)
  unsigned int m_StripXEMult; //!
  unsigned int m_StripYEMult; //!
  unsigned int m_StripXTMult; //!
  unsigned int m_StripYTMult; //!
  unsigned int m_SiLiEMult;   //!
  unsigned int m_SiLiTMult;   //!
  unsigned int m_CsIEMult;    //!
  unsigned int m_CsITMult;    //!

  //   Event over this value after pre-treatment are not treated / avoid long
  //   treatment time on spurious event
  unsigned int m_MaximumStripMultiplicityAllowed; //!
  //   Give the allowance in percent of the difference in energy between X and Y
  double m_StripEnergyMatchingSigma;         //!
  double m_StripEnergyMatchingNumberOfSigma; //!

  // Raw Threshold
  int m_Si_X_E_RAW_Threshold; //!
  int m_Si_Y_E_RAW_Threshold; //!
  int m_SiLi_E_RAW_Threshold; //!
  int m_CsI_E_RAW_Threshold;  //!

  // Calibrated Threshold
  double m_Si_X_E_Threshold; //!
  double m_Si_Y_E_Threshold; //!
  double m_SiLi_E_Threshold; //!
  double m_CsI_E_Threshold;  //!

  // Geometric Matching
  // size in strip of a pad
  int m_SiLi_Size; //!
  // center position of the pad on X
  vector<int> m_SiLi_MatchingX; //!
  // center position of the pad on Y
  vector<int> m_SiLi_MatchingY; //!
  // size in strip of a cristal
  unsigned int m_CsI_Size; //!
  // center position of the cristal on X
  vector<int> m_CsI_MatchingX;                                   //!
  std::map<int, std::pair<int, int>> m_ZeroDegree_CsI_MatchingX; //!
  // center position of the cristal on X
  vector<int> m_CsI_MatchingY;                                   //!
  std::map<int, std::pair<int, int>> m_ZeroDegree_CsI_MatchingY; //!

  // If set to true, all event that do not come in front of a cristal will be
  // ignore all time (crossing or not),
  // Warning, this option reduce statistic, however it help eliminating
  // unrealevent event that cross the DSSD
  // And go between pad or cristal.
  bool m_Ignore_not_matching_SiLi; //!
  bool m_Ignore_not_matching_CsI;  //!
  
  
  /////////////////////////// CALIBRATION RELATED METHODS ////////////////////////////////
  
  void DoCalibrationCSIPreTreat();//!

  void InitializeRootHistogramsCalib(); //!

  void InitializeRootHistogramsEnergyF(Int_t DetectorNumber); //!

  void InitializeRootHistogramsTimeF(Int_t DetectorNumber){}; //!

  void InitializeRootHistogramsCSIF(Int_t DetectorNumber); //!

  void FillHistogramsCalib(); //!

  void FillHistogramsCalibEnergyF(); //!

  void FillHistogramsCalibTimeF(){}; //!

  void FillHistogramsCalibCSIF(); //!

  void DoCalibration(); //!

  void DoCalibrationEnergyF(Int_t DetectorNumber); //!

  void DoCalibrationTimeF(Int_t DetectorNumber); //!

  void DoCalibrationCsIF(Int_t DetectorNumber); //!

  void MakeEnergyCalibFolders(); //!

  void MakeCSICalibFolders(); //!

  void CreateCalibrationEnergyFiles(unsigned int DetectorNumber, TString side, ofstream* calib_file,
                                    ofstream* dispersion_file); //!

  void CreateCalibrationCSIFiles(unsigned int DetectorNumber, ofstream* calib_file, TString ParticleType); //!

  void CloseCalibrationCSIFiles(ofstream* calib_file); //!

  void CloseCalibrationEnergyFiles(ofstream* calib_file, ofstream* dispersion_file); //!

  bool FindAlphas(TH1* CalibHist, TString side, unsigned int StripNb, unsigned int DetectorNumber); //!

  void FitLinearEnergy(TGraphErrors* FitHist, TString side, unsigned int StripNb, unsigned int DetectorNumber,
                       double* a, double* b, std::vector<double> Source_E); //!

  std::vector<double> SlowSource(double AlThickness);//!

  double FindMeanExtrapolation(TGraphErrors* DispersionHist);//!

  void WriteHistogramsCalib(); //!

  void WriteHistogramsEnergyF(); //!

  void WriteHistogramsCSIF(); //!

  void WriteHistogramsTimeF(){}; //!

  static Double_t source_Pu(Double_t* x, Double_t* par); //!
  static Double_t source_Am(Double_t* x, Double_t* par); //!
  static Double_t source_Cm(Double_t* x, Double_t* par); //!

  void DefineCalibrationSource(); //!

 private:                        //   Root Input and Output tree classes
  TMust2Data* m_EventData;       //!
  TMust2Data* m_PreTreatedData;  //!
  TMust2Physics* m_EventPhysics; //!
  // TMust2Data* r_EventData;

 private:                                     //   Map of activated channel
  map<int, vector<bool>> m_XChannelStatus;    //!
  map<int, vector<bool>> m_YChannelStatus;    //!
  map<int, vector<bool>> m_SiLiChannelStatus; //!
  map<int, vector<bool>> m_CsIChannelStatus;  //!

 private:
  int m_NumberOfTelescope; //!

  vector<vector<vector<double>>> m_StripPositionX; //!
  vector<vector<vector<double>>> m_StripPositionY; //!
  vector<vector<vector<double>>> m_StripPositionZ; //!

 public:
  // Prevent to treat event with ambiguous matching beetween X and Y
  bool m_multimatch;             //!
  vector<int> m_match_type;      //!
  map<int, int> m_NMatchDet;     //!
  map<int, int> m_StripXMultDet; //!
  map<int, int> m_StripYMultDet; //!
  map<int, int> m_NMatchX;       //!
  map<int, int> m_NMatchY;       //!

 private:
  map<int, bool> m_CsIPresent;  //!
  map<int, bool> m_SiLiPresent; //!
  map<int, bool> m_CsIOffset;   //!
  
  
  /////////////////////////// CALIBRATION RELATED VARIABLES ////////////////////////////////
 private:
  map<int, bool> DoCalibrationEnergy;                               //!
  map<int, bool> DoCalibrationTime;                                 //!
  map<int, bool> DoCalibrationCsI;                                  //!
  bool IsCalibCSI = false;                                          //!
  bool IsCalibEnergy = false;                                       //!
  bool Cal_Pixel; //!
  unsigned int PixelSize; //!
  std::map<TString, std::map<unsigned int, unsigned int>> BadStrip; //!
  std::map<unsigned int,std::vector<double>> AlphaSigma;                                    //!
  std::map<unsigned int,std::vector<double>> AlphaMean;                                    //!
  std::vector<TString> Source_isotope;                              //!
  std::vector<double> Source_E;                                     //!
  std::vector<double> Source_Sig;                                   //!
  std::vector<double> Source_branching_ratio;                       //!
  // ofstream peaks_file, calib_file, dispersion_file , calib_online_file, latex_file;//!

  ///////// Calib parameters for Si detectors
  map<int, double> EnergyXThreshold;  //!
  map<int, double> EnergyYThreshold;  //!
  map<int, std::string> AlphaFitType; //!

  ///////// Calib parameters for CsI detectors
  map<int, double> CSIEnergyXThreshold;             //!
  map<int, double> CSIEnergyYThreshold;             //!
  map<int, double> CSIEThreshold;                   //!
  map<int, double> SiThickness;                   //!
  map<int, double> AlThickness;                   //!
  TTreeReaderValue<std::vector<unsigned int>>* GATCONF_; //!
  bool DoCSIFit;                                    //!
  std::map<TString, NPL::EnergyLoss*> ParticleSi;   //!
  std::map<TString, NPL::EnergyLoss*> ParticleAl;   //!
  std::vector<string> ParticleType{"proton", "deuteron", "triton", "alpha"}; //!
  std::vector<string> ParticleTypePixel{"proton", "alpha"}; //!
  
  NPL::EnergyLoss* Alpha_Al = nullptr;//!

 private:                   // Spectra Class
  TMust2Spectra* m_Spectra; //!

 public:
  void WriteSpectra(); //!

 public: // Spectra Getter
  map<string, TH1*> GetSpectra();

 public: // Static constructor to be passed to the Detector Factory
  static NPL::VDetector* Construct();
  static NPL::VTreeReader* ConstructReader();
  ClassDef(TMust2Physics, 1) // Must2Physics structure
};

namespace MUST2_LOCAL {
  //   DSSD
  //   X
  double fSi_X_E(const TMust2Data* Data, const int& i);
  double fSi_X_T(const TMust2Data* Data, const int& i);

  //   Y
  double fSi_Y_E(const TMust2Data* Data, const int& i);
  double fSi_Y_T(const TMust2Data* Data, const int& i);

  //   SiLi
  double fSiLi_E(const TMust2Data* Data, const int& i);
  double fSiLi_T(const TMust2Data* Data, const int& i);

  //   CsI
  double fCsI_E(const TMust2Data* Data, const int& i);
  double fCsI_T(const TMust2Data* Data, const int& i);
} // namespace MUST2_LOCAL

#endif
