#ifndef TMUSETTPHYSICS_H
#define TMUSETTPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : June 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold MUSETT Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <map>
#include <vector>
#include <string>
// NPL
#include "NPCalibrationManager.h"
#include "NPInputParser.h"
#include "NPVDetector.h"
#include "TMUSETTData.h"
 #include "TMUSETTSpectra.h"
// ROOT
#include "TH1.h"
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom3.h"

using namespace std;

// Forward Declaration
 class TMUSETTSpectra;

class TMUSETTPhysics : public TObject, public NPL::VDetector {
  public:
  TMUSETTPhysics();
  ~TMUSETTPhysics();

  public:
  void Clear();
  void Clear(const Option_t*){};

  public:
  vector<TVector2> Match_X_Y();
  int  CheckEvent();

  TRandom3* m_random; //!

  public:
  //   Provide Physical Multiplicity
  Int_t EventMultiplicity;

  Int_t EventTrigger;


  // Detector
  vector<int> DetectorNumber;
   //   DSSD
  vector<double> DSSD_E;
  vector<double> DSSD_E_X;
  vector<double> DSSD_E_Y;
  vector<double> DSSD_T;
  vector<double> DSSD_T_X;
  vector<double> DSSD_T_Y;
  vector<int>    DSSD_X;
  vector<int>    DSSD_Y;

  //Garder E_x, E_y, T_x, T_y


  //Rajouter strip number
  vector<double> PosX;
  vector<double> PosY;
  vector<double> PosZ;
  vector<double> Theta;//!





  private:


  public: //   Innherited from VDetector Class
  //   Read stream at ConfigFile to pick-up parameters of detector
  //   (Position,...) using Token
  void ReadConfiguration(NPL::InputParser parser);

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
  // void InitSpectra();
  // Fill the spectra hold by the spectra class
  // void FillSpectra();
  // Used for Online mainly, perform check on the histo and for example change
  // their color if issues are found
  // void CheckSpectra();
  // Used for Online only, clear all the spectra hold by the Spectra class
  // void ClearSpectra();

  public: //   Specific to MUSETT Array
  //   Clear The PreTeated object
  void ClearPreTreatedData() { m_PreTreatedData->Clear(); }

  //   Remove bad channel, calibrate the data and apply threshold
  void PreTreat();

  //   Return false if the channel is disabled by user
  //   Frist argument is either 0 for X,1 Y,2 SecondLayer 3
  bool IsValidChannel(const int& Type, const int& telescope,
                      const int& channel);

  //   Initialize the standard parameter for analysis
  //   ie: all channel enable, maximum multiplicity for strip = number of
  //   telescope
  void InitializeStandardParameter();

  //   Read the user configuration file; if no file found, load standard one
  void ReadAnalysisConfig();

  //   Add a Detector using Corner Coordinate information
  void AddDetector(TVector3 C_X1_Y1, TVector3 C_X128_Y1, TVector3 C_X1_Y128,
                    TVector3 C_X128_Y128);

  //   Add a Detector using R Theta Phi of Si center information
  void AddDetector(double theta, double phi, double distance, double beta_u,
                    double beta_v, double beta_w);

 //   Special Method for Annular S1
  void AddDetector(TVector3 C_Center);

  // Give and external TMustData object to TMUSETTPhysics. Needed for online
  // analysis for example.
  void SetRawDataPointer(void* rawDataPointer) {
    m_EventData = (TMUSETTData*)rawDataPointer;
  }
  // Retrieve raw and pre-treated data
  TMUSETTData* GetRawData() const { return m_EventData; }
  TMUSETTData* GetPreTreatedData() const { return m_PreTreatedData; }


  // Use to access the strip position
  double GetStripPositionX(const int N, const int X, const int Y) {
    // if (N==9)
    // cout << N << " " << X << " " << Y << " " << m_DetectorNumberIndex[N] << " " << m_StripPositionX[ m_DetectorNumberIndex[N] - 1][X - 1][Y - 1] << endl;
    return m_StripPositionX[N][X][Y];
  };
  double GetStripPositionY(const int N, const int X, const int Y) {
    return m_StripPositionY[N][X][Y];
  };
  double GetStripPositionZ(const int N, const int X, const int Y) {
    return m_StripPositionZ[N][X][Y];
  };

  double GetNumberOfDetector() const { return m_NumberOfDetector; };

  // To be called after a build Physical Even
  int GetEventMultiplicity() const { return EventMultiplicity; };

  int GetEventTrigger() const {return EventTrigger;};


  //double GetEnergyDeposit(const int i) const { return TotalEnergy[i]; };

  TVector3 GetPositionOfInteraction(const int i,bool random=false) ;
  TVector3 GetDetectorNormal(const int i) ;
  double GetRelativeAngle();
  double GetS(const unsigned int A, const unsigned int Q2A);

  private: //   Parameter used in the analysis
  // Shape of the detector Trapezoid or Square

  // By default take EX and TY.
  bool m_Take_E_Y; //!
  bool m_Take_T_Y; //!

  UShort_t StripLimit = 0; //!

  //   Event over this value after pre-treatment are not treated / avoid long
  //   treatment time on spurious event
  unsigned int m_MaximumStripMultiplicityAllowed; //!
  //   Give the allowance in percent of the difference in energy between X and Y
  double m_StripEnergyMatching; //!

  // Raw Threshold
  int m_DSSD_X_E_RAW_Threshold; //!
  int m_DSSD_Y_E_RAW_Threshold; //!
  int m_SecondLayer_E_RAW_Threshold; //!

  // Calibrated Threshold
  double m_DSSD_X_E_Threshold; //!
  double m_DSSD_Y_E_Threshold; //!
  double m_SecondLayer_E_Threshold; //!


  private: //   Root Input and Output tree classes
  TMUSETTData*    m_EventData; //!
  TMUSETTData*    m_PreTreatedData; //!
  TMUSETTPhysics* m_EventPhysics; //!

  private: //   Map of activated channel
  map<int, vector<bool>> m_XChannelStatus; //!
  map<int, vector<bool>> m_YChannelStatus; //!

  private: // Spatial Position of Strip Calculated on bases of detector position
  int m_NumberOfDetector; //!

  vector<vector<vector<double>>> m_StripPositionX; //!
  vector<vector<vector<double>>> m_StripPositionY; //!
  vector<vector<vector<double>>> m_StripPositionZ; //!

  private:
  map<int, int>    m_HitDSSDX; //!
  map<int, int>    m_HitDSSDY; //!

  private: // Spectra Class
  // TMUSETTSpectra* m_Spectra; //!

  public:
  // void WriteSpectra(); //!

  public: // Spectra Getter
  // map<string, TH1*> GetSpectra();

  public: // Static constructor to be passed to the Detector Factory
  static NPL::VDetector* Construct();
  ClassDef(TMUSETTPhysics, 1) // MUSETTPhysics structure
};

namespace MUSETT_LOCAL {
  //   DSSD
  //   X
  double fDSSD_X_E(const TMUSETTData* Data, const int& i);
  double fDSSD_X_T(const TMUSETTData* Data, const int& i);

  //   Y
  double fDSSD_Y_E(const TMUSETTData* Data, const int& i);
  double fDSSD_Y_T(const TMUSETTData* Data, const int& i);

  //  Second Layer
  double fSecondLayer_E(const TMUSETTData* Data, const int& i);
  double fSecondLayer_T(const TMUSETTData* Data, const int& i);
}

#endif
