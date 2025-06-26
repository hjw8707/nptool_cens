#ifndef TZDDPHYSICS_H
#define TZDDPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2022   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : October 2022                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold ZDD Treated data                                *
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
#include "TZDDData.h"
#include "TZDDCalData.h"
#include "TZDDSpectra.h"
#include "TZDDPhysicsReader.h"
#include "NPCalibrationManager.h"
#include "NPTimeStamp.h"
#include "NPVDetector.h"
#include "NPVTreeReader.h"
#include "NPInputParser.h"
// forward declaration
class TZDDSpectra;



class TZDDPhysics : public TObject, public NPL::VDetector, public TZDDPhysicsReader {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TZDDPhysics();
    ~TZDDPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    vector<int>      IC_Nbr;
    vector<double>   IC_E;
    vector<unsigned long long>   IC_TS;
    double ICSum;

    vector<int>      PL_Nbr;
    vector<double>   PL_E;
    vector<unsigned long long>   PL_TS;

    vector<int>      DC_DetectorNumber;
    vector<double>   DC_DriftTime;
  
  /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 POS, string shape); 
  void AddDetector(double R, double Theta, double Phi, string shape); 
  
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

    // methods related to the TZDDSpectra class
    // instantiate the TZDDSpectra class and 
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
  
    void SetTreeReader(TTreeReader* TreeReader);

    void ReadConfigurationTS(); 
  //////////////////////////////////////////////////////////////
  // specific methods to ZDD array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();


    void Treat_DC();
    // Matching IC
    void Match_IC();
    void Match_IC1();
    
    void Match_PL();
    // bool CheckGoodEvent();

    // PreTreating Energy for IC and Plastic 
    void PreTreatEnergy(std::string Detector, CalibrationManager* Cal);
    
    // Same for time
    void PreTreatTime(std::string Detector, CalibrationManager* Cal);

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    void ClearTSRef(){RefTS_Name=""; RefTS=0;}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TZDDData object to TZDDPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TZDDData* rawDataPointer) {m_EventData = rawDataPointer;}

    // Setting TS Reference
    void SetRefTS(std::string TSRef_Name, unsigned long long TSRef);

    void Add_ZDD(double R, double Theta) {
        m_R     = R;
        m_Theta = Theta;
        m_NumberOfDetectors++;
    }
    bool fSortIC(int i, int j);

  private:
    void ClaimReaderData();

  // objects are not written in the TTree
  private:
    TZDDData*         m_EventData;        //!
    TZDDCalData*      m_PreTreatedData;   //!
    TZDDPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TZDDData* GetRawData()        const {return m_EventData;}
    TZDDCalData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    // thresholds
    double m_IC_E_RAW_Threshold; //!
    double m_PL_E_RAW_Threshold; //!
    double m_DC_E_RAW_Threshold; //!
    double m_EXO_E_RAW_Threshold; //!
    double m_E_Threshold;     //!
    double m_R;//!
    double m_Theta;//!
    unsigned int m_IC_Mult; //!
    unsigned int m_PL_Mult; //!
    unsigned int m_DC_Mult; //!
    unsigned int m_EXO_Mult; //!
    std::map<unsigned int, unsigned int> Map_IC;//!
    std::map<unsigned int,std::pair<unsigned int, unsigned long long>> SortIC;//!
    std::map<unsigned int,std::pair<unsigned int, unsigned long long>> SortPL;//!

  private:
    int m_NumberOfDetectors;  //!

  private:
    TimeStamp* TSEvent;//!
    std::string RefTS_Name = "";//!
    unsigned long long RefTS = 0;//!

  public:

  // spectra class
  private:
    TZDDSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();
    static NPL::VTreeReader* ConstructReader();

    ClassDef(TZDDPhysics,1)  // ZDDPhysics structure
};
#endif
