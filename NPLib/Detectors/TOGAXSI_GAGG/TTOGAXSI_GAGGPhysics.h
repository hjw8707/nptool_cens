#ifndef TTOGAXSI_GAGGPHYSICS_H
#define TTOGAXSI_GAGGPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_GAGG Treated data                                *
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
#include "TTOGAXSI_GAGGData.h"
#include "TTOGAXSI_GAGGSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TTOGAXSI_GAGGSpectra;



class TTOGAXSI_GAGGPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TTOGAXSI_GAGGPhysics();
    ~TTOGAXSI_GAGGPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    Int_t EventMultiplicity;
    vector<int>      RecoilDetectorNumber;
    vector<double>   RecoilE;
    vector<double>   RecoilT;

    vector<int>      ClusterDetectorNumber;
    vector<double>   ClusterE;
    vector<double>   ClusterT;

  /// A usefull method to bundle all operation to add a detector
  //void AddDetector(TVector3 POS, string shape); 
  //void AddDetector(double R, double Theta, double Phi, string shape); 
  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    // A useful method to bundle all operation to add a detector
    void AddRecoilDetector(TVector3 Pos, double Phi, TVector3 Ref);
    void AddClusterDetector(TVector3 Pos, double Phi, TVector3 Ref);

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

    // methods related to the TTOGAXSI_GAGGSpectra class
    // instantiate the TTOGAXSI_GAGGSpectra class and 
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
  // specific methods to TOGAXSI_GAGG array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TTOGAXSI_GAGGData object to TTOGAXSI_GAGGPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TTOGAXSI_GAGGData* rawDataPointer) {m_EventData = rawDataPointer;}
  
    double GetNumberOfRecoilDetector() const {return m_NumberOfRecoilDetectors;}
    double GetNumberOfClusterDetector() const {return m_NumberOfClusterDetectors;}
    int GetEventMultiplicity() const {return EventMultiplicity;}

  // objects are not written in the TTree
  private:
    TTOGAXSI_GAGGData*         m_EventData;        //!
    TTOGAXSI_GAGGData*         m_PreTreatedData;   //!
    TTOGAXSI_GAGGPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TTOGAXSI_GAGGData* GetRawData()        const {return m_EventData;}
    TTOGAXSI_GAGGData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    int m_NumberOfRecoilDetectors; //!
    int m_NumberOfClusterDetectors; //!

    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!

  // spectra class
  private:
    TTOGAXSI_GAGGSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TTOGAXSI_GAGGPhysics,1)  // TOGAXSI_GAGGPhysics structure

  private: //geometry

  double Crystal_Length;
  double Crystal_Width;
  double Crystal_Height;

};
#endif
