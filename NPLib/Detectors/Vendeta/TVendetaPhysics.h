#ifndef TVendetaPHYSICS_H
#define TVendetaPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace@cea.fr                        *
 *                                                                           *
 * Creation Date  : February 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Vendeta Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// C++ headers 
#include <vector>
#include <map>
#include <string>
using namespace std;

// ROOT headers
#include "TObject.h"
#include "TH1.h"
#include "TVector3.h"
// NPTool headers
#include "TVendetaData.h"
#include "TVendetaSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TVendetaSpectra;



class TVendetaPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TVendetaPhysics();
    ~TVendetaPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      LG_DetectorNumber;
    vector<double>   LG_Q1;
    vector<double>   LG_Q2;
    vector<double>   LG_Time;
    vector<double>   LG_Qmax;
    vector<int>      LG_isSat;
    vector<int>      LG_isGamma;

    vector<int>      HG_DetectorNumber;
    vector<double>   HG_Q1;
    vector<double>   HG_Q2;
    vector<double>   HG_Time;
    vector<double>   HG_Qmax;
    vector<int>      HG_isSat;
    vector<int>      HG_isGamma;


  /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 POS); 
  void AddDetector(double R, double Theta, double Phi); 
 
  double GetDistanceFromTarget(int DetNbr) {return m_DetectorPosition[DetNbr-1].Mag();}
  TVector3 GetVectorDetectorPosition(int DetNbr) {return m_DetectorPosition[DetNbr-1];}

  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    // add parameters to the CalibrationManger
    void AddParameterToCalibrationManager();

    // method called event by event, aiming at extracting the 
    // physical information from detector
    void BuildPhysicalEvent();

    // same as BuildPhysicalEvent() method but with a simpler
    // treatment
    void BuildSimplePhysicalEvent();

    // same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // activate raw data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputRaw();

    // activate physics data object and branches from input TChain
    // in this method mother branches (Detector) AND daughter leaves 
    // (fDetector_parameter) have to be activated
    void InitializeRootInputPhysics();

    // create branches of output ROOT file
    void InitializeRootOutput();

    // clear the raw and physical data objects event by event
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   

    // methods related to the TVendetaSpectra class
    // instantiate the TVendetaSpectra class and 
    // declare list of histograms
    void InitSpectra();

    // fill the spectra
    void FillSpectra();

    // used for Online mainly, sanity check for histograms and 
    // change their color if issues are found, for example
    void CheckSpectra();

    // used for Online only, clear all the spectra
    void ClearSpectra();

    // write spectra to ROOT output file
    void WriteSpectra();


  //////////////////////////////////////////////////////////////
  // specific methods to Vendeta array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TVendetaData object to TVendetaPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TVendetaData* rawDataPointer) {m_EventData = rawDataPointer;}
   
    void SetAnodeNumber(int AnodeNbr) {m_AnodeNumber = AnodeNbr;}; //!

  // objects are not written in the TTree
  private:
    TVendetaData*         m_EventData;        //!
    TVendetaData*         m_PreTreatedData;   //!
    TVendetaPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TVendetaData* GetRawData()        const {return m_EventData;}
    TVendetaData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!

  private:
    // Anode number from fission chamber for time offset calibration
    int m_AnodeNumber; //!
  // number of detectors
  private:
    int m_NumberOfDetectors;  //!
    vector<TVector3> m_DetectorPosition; //!

  // spectra class
  private:
    TVendetaSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TVendetaPhysics,1)  // VendetaPhysics structure
};
#endif
