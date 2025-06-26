#ifndef TISSPHYSICS_H
#define TISSPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2019    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 * Author: M. Labiche                     address: marc.labiche@stfc.ac.uk   *
 *                                                                           *
 * Creation Date  : july 2019                                                *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Iss treated data                                         *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <map>
#include <vector>
// NPL
#include "NPCalibrationManager.h"
#include "NPInputParser.h"
#include "NPVDetector.h"
#include "TIssData.h"
#include "TIssSpectra.h"

// ROOT
#include "TH1.h"
#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom3.h"

using namespace std;

// Forward Declaration
class TIssSpectra;

class TIssPhysics : public TObject, public NPL::VDetector {
public:
  TIssPhysics();
  ~TIssPhysics();

public:
  void Clear();
  void Clear(const Option_t*){};

public:
    vector < TVector2 > Match_Front_Back() ;
    int  CheckEvent();

/* from MUST2
  vector<TVector2> Match_X_Y();
  int  CheckEvent(int N);
  bool Match_Si_CsI(int X, int Y, int CristalNbr);
  bool Match_Si_SiLi(int X, int Y, int PadNbr);
  bool ResolvePseudoEvent();
*/
public:
  //   Provide Physical Multiplicity
  // Int_t EventMultiplicity;
  int EventMultiplicity;

  //   Provide a Classification of Event
  vector<int> EventType;

    // Detector
    vector<int> DetectorNumber ;

    //   DSSD
    vector<double> Strip_E ;
    vector<double> Strip_T ;
    vector<double> StripFront_E ;
    vector<double> StripFront_T ;
    vector<double> StripBack_E ;
    vector<double> StripBack_T ;
    vector<int>    Strip_Front ;
    vector<int>    Strip_Back ;

    // Used to apply Pixel Cal
    vector<double> StripFront_OriginalE; //!
    vector<double> StripBack_OriginalE; //!
    vector<double> DeadLayer; //!    
    // Used for Calibration
    vector<double> Strip_Front_RawE;
    vector<double> Strip_Back_RawE;


/* from MUST2
  // Telescope
  vector<int> TelescopeNumber;

  //   Si
  vector<double> Si_E;
  vector<double> Si_T;
  vector<int>    Si_X;
  vector<int>    Si_Y;

  // Use for checking purpose
  vector<double> Si_EX;
  vector<double> Si_TX;
  vector<double> Si_EY;
  vector<double> Si_TY;
  vector<int>    TelescopeNumber_X;
  vector<int>    TelescopeNumber_Y;


  // Physical Value
  vector<double> TotalEnergy;
*/


public: //   Innherited from VDetector Class

    //   Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
    void ReadConfiguration(NPL::InputParser) ;


    //   Add Parameter to the CalibrationManger
    void AddParameterToCalibrationManager() ;      

    //   Activated associated Branches and link it to the private member DetectorData address
    //   In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw() ;

    //   Activated associated Branches and link it to the private member DetectorPhysics address
    //   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ;

    //   Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput() ;

    //   This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
    void BuildPhysicalEvent() ;

    //   Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
    //   This method aimed to be used for analysis performed during experiment, when speed is requiered.
    //   NB: This method can eventually be the same as BuildPhysicalEvent.
    void BuildSimplePhysicalEvent() ;

    // Same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();}

    //   Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra();
    // Fill the spectra hold by the spectra class
    void FillSpectra();
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra();
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra();
    // Write Spectra to file
    void WriteSpectra();



/* from MUST2

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
  void InitSpectra();
  // Fill the spectra hold by the spectra class
  void FillSpectra();
  // Used for Online mainly, perform check on the histo and for example change
  // their color if issues are found
  void CheckSpectra();
  // Used for Online only, clear all the spectra hold by the Spectra class
  void ClearSpectra();
*/


public: //   Specific to Iss Array

   //   Clear The PreTeated object
    void ClearPreTreatedData() {m_PreTreatedData->Clear();}

    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    //   Return false if the channel is disabled by user
    //   Frist argument is either "X","Y","SiLi","CsI"
    bool IsValidChannel(const string& DetectorType, const int& telescope , const int& channel);

    //   Initialize the standard parameter for analysis
    //   ie: all channel enable, maximum multiplicity for strip = number of telescope
    void InitializeStandardParameter();

    //   Read the user configuration file; if no file found, load standard one
    void ReadAnalysisConfig();

 
    //   Add a Detector
    void AddBoxDetector( double Z);
//    void AddQQQDetector( double R,double Phi,double Z);
    //   Add a Telescope using Corner Coordinate information
    void AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1, TVector3 C_X1_Y128, TVector3 C_X128_Y128);

    //   Add a Telescope using R Theta Phi of Si center information
    void AddTelescope(double theta, double phi, double distance, double beta_u, double beta_v, double beta_w);

    double GetNominalMField() { return m_NominalField; }


    // Give and external TMustData object to TSharcPhysics. Needed for online analysis for example.
    void SetRawDataPointer(TIssData* rawDataPointer) {m_EventData = rawDataPointer;}
    // Retrieve raw and pre-treated data
    TIssData* GetRawData()        const {return m_EventData;}
    TIssData* GetPreTreatedData() const {return m_PreTreatedData;}

    // Use to access the strip position
    inline double GetStripPositionX( const int& N , const int& Front , const int& Back )   const{ return m_StripPositionX[N-1][Front-1][Back-1] ; }  ;
    inline double GetStripPositionY( const int& N , const int& Front , const int& Back )   const{ return m_StripPositionY[N-1][Front-1][Back-1] ; }  ;
    inline double GetStripPositionZ( const int& N , const int& Front , const int& Back )   const{ return m_StripPositionZ[N-1][Front-1][Back-1] ; }  ;

    inline double GetNumberOfDetector() const { return m_NumberOfDetector; };

    // To be called after a build Physical Event 
    inline int GetEventMultiplicity() const { return EventMultiplicity; };

    TVector3 GetPositionOfInteraction(const int& i, bool random=false) const;   
    TVector3 GetDetectorNormal(const int& i) const;
    double   GetDeadLayer(const int& i) const;



/*from MUST2

  //   Clear The PreTeated object
  void ClearPreTreatedData() { m_PreTreatedData->Clear(); }

  //   Remove bad channel, calibrate the data and apply threshold
  void PreTreat();

  //   Return false if the channel is disabled by user
  //   Frist argument is either 0 for X,1 Y,2 SiLi, 3 CsI
  bool IsValidChannel(const int& DetectorType, const int& telescope,
                      const int& channel);

  //   Initialize the standard parameter for analysis
  //   ie: all channel enable, maximum multiplicity for strip = number of
  //   telescope
  void InitializeStandardParameter();

  //   Read the user configuration file; if no file found, load standard one
  void ReadAnalysisConfig();

  //   Add a Telescope using Corner Coordinate information
  void AddTelescope(TVector3 C_X1_Y1, TVector3 C_X128_Y1, TVector3 C_X1_Y128,
                    TVector3 C_X128_Y128);

  //   Add a Telescope using R Theta Phi of Si center information
  void AddTelescope(double theta, double phi, double distance, double beta_u,
                    double beta_v, double beta_w);

  // Use for reading Calibration Run, very simple methods; only apply
  // calibration, no condition
  void ReadCalibrationRun();

  // Give and external TMustData object to TIssPhysics. Needed for online
  // analysis for example.
  void SetRawDataPointer(void* rawDataPointer) {
    m_EventData = (TIssData*)rawDataPointer;
  }
  // Retrieve raw and pre-treated data
  TIssData* GetRawData() const { return m_EventData; }
  TIssData* GetPreTreatedData() const { return m_PreTreatedData; }

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


*/


private: //   Parameter used in the analysis

   // By default take EX and TY.
    bool m_Take_E_Front;//!
    bool m_Take_T_Back;//!

    //   Event over this value after pre-treatment are not treated / avoid long treatment time on spurious event   
    unsigned int m_MaximumStripMultiplicityAllowed  ;//!
    //   Give the allowance in percent of the difference in energy between X and Y
    double m_StripEnergyMatchingSigma  ; //!
    double m_StripEnergyMatchingNumberOfSigma  ; //!

    //  Threshold
    double m_StripFront_E_RAW_Threshold ;//!
    double m_StripFront_E_Threshold ;//!
    double m_StripBack_E_RAW_Threshold ;//!
    double m_StripBack_E_Threshold ;//!
 



/*from MUST2

  // By default take EX and TY.
  bool m_Take_E_Y; //!
  bool m_Take_T_Y; //!

  // Size Container to be used in the different loop of the analysis (value must
  // be given locally)
  unsigned int m_StripXEMult; //!
  unsigned int m_StripYEMult; //!
  unsigned int m_StripXTMult; //!
  unsigned int m_StripYTMult; //!
  unsigned int m_SiLiEMult; //!
  unsigned int m_SiLiTMult; //!
  unsigned int m_CsIEMult; //!
  unsigned int m_CsITMult; //!

  //   Event over this value after pre-treatment are not treated / avoid long
  //   treatment time on spurious event
  unsigned int m_MaximumStripMultiplicityAllowed; //!
  //   Give the allowance in percent of the difference in energy between X and Y
  double m_StripEnergyMatchingSigma; //!
  double m_StripEnergyMatchingNumberOfSigma; //!

  // Raw Threshold
  int m_Si_X_E_RAW_Threshold; //!
  int m_Si_Y_E_RAW_Threshold; //!
  int m_SiLi_E_RAW_Threshold; //!
  int m_CsI_E_RAW_Threshold; //!

  // Calibrated Threshold
  double m_Si_X_E_Threshold; //!
  double m_Si_Y_E_Threshold; //!
  double m_SiLi_E_Threshold; //!
  double m_CsI_E_Threshold; //!

  // Geometric Matching
  // size in strip of a pad
  int m_SiLi_Size; //!
  // center position of the pad on X
  vector<int> m_SiLi_MatchingX; //!
  // center position of the pad on Y
  vector<int> m_SiLi_MatchingY; //!
  // size in strip of a cristal
  int m_CsI_Size; //!
  // center position of the cristal on X
  vector<int> m_CsI_MatchingX; //!
  // center position of the cristal on X
  vector<int> m_CsI_MatchingY; //!

  // If set to true, all event that do not come in front of a cristal will be
  // ignore all time (crossing or not),
  // Warning, this option reduce statistic, however it help eliminating
  // unrealevent event that cross the DSSD
  // And go between pad or cristal.
  bool m_Ignore_not_matching_SiLi; //!
  bool m_Ignore_not_matching_CsI; //!

*/


private: //   Root Input and Output tree classes
  TIssData*    m_EventData; //!
  TIssData*    m_PreTreatedData; //!
  TIssPhysics* m_EventPhysics; //!

private: //   Map of activated channel

   map< int, vector<bool> > m_FrontChannelStatus;//!
   map< int, vector<bool> > m_BackChannelStatus;//! 


/* from MUST2
  map<int, vector<bool>> m_XChannelStatus; //!
  map<int, vector<bool>> m_YChannelStatus; //!
  map<int, vector<bool>> m_SiLiChannelStatus; //!
  map<int, vector<bool>> m_CsIChannelStatus; //!
*/


  private:   //   Spatial Position of Strip Calculated on bases of detector position

    int m_NumberOfDetector;//!
    vector< vector < vector < double > > >   m_StripPositionX;//!
    vector< vector < vector < double > > >   m_StripPositionY;//!
    vector< vector < vector < double > > >   m_StripPositionZ;//!
    vector< TVector3 > m_DetectorNormal;//!
    vector< TVector3 > m_U;//!
    vector< TVector3 > m_V;//!
    TRandom3* m_Rand;//!
    double m_BoxPitchBack ;//!
    double m_BoxPitchFront;//!

    double m_NominalField;



/* from MUST2
private:
  int m_NumberOfTelescope; //!

  vector<vector<vector<double>>> m_StripPositionX; //!
  vector<vector<vector<double>>> m_StripPositionY; //!
  vector<vector<vector<double>>> m_StripPositionZ; //!


public:
  // Prevent to treat event with ambiguous matching beetween X and Y
  bool          m_multimatch; //!
  vector<int>   m_match_type; //!
  map<int, int> m_NMatchDet; //!
  map<int, int> m_StripXMultDet; //!
  map<int, int> m_StripYMultDet; //!

private:
  map<int, bool> m_CsIPresent; //!
  map<int, bool> m_SiLiPresent; //!
*/


private: // Spectra Class
  TIssSpectra* m_Spectra; //!


public: // Spectra Getter
  map<string, TH1*> GetSpectra();

public: // Static constructor to be passed to the Detector Factory
  static NPL::VDetector* Construct();
  ClassDef(TIssPhysics, 1) // IssPhysics structure
};

namespace ISS_LOCAL {
//   DSSD

  //   Front
  double fStrip_Front_E(const TIssData* Data, const int& i);
  double fStrip_Front_T(const TIssData* Data, const int& i);

  //   Back   
  double fStrip_Back_E(const TIssData* Data, const int& i);
  double fStrip_Back_T(const TIssData* Data, const int& i);


/* from MUST2
//   X
double fSi_X_E(const TIssData* Data, const int& i);
double fSi_X_T(const TIssData* Data, const int& i);

//   Y
double fSi_Y_E(const TIssData* Data, const int& i);
double fSi_Y_T(const TIssData* Data, const int& i);

//   SiLi
double fSiLi_E(const TIssData* Data, const int& i);
double fSiLi_T(const TIssData* Data, const int& i);

//   CsI
double fCsI_E(const TIssData* Data, const int& i);
double fCsI_T(const TIssData* Data, const int& i);
*/
} // namespace ISS_LOCAL

#endif
