#ifndef TTACPHYSICS_H
#define TTACPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Hugo Jacob  contact address: hjacob@ijclab.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TAC Treated data                                *
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
#include "TTACData.h"
#include "TTACSpectra.h"
#include "TTACPhysicsReader.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPTimeStamp.h"
// forward declaration
class TTACSpectra;



class TTACPhysics : public TObject, public NPL::VDetector, public TTACPhysicsReader {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TTACPhysics();
    ~TTACPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    std::vector<unsigned int> TAC_Time;
    std::vector<std::string> TAC_Name;
    std::vector<unsigned long long> TAC_TS;


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

    // methods related to the TTACSpectra class
    // instantiate the TTACSpectra class and 
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


  //////////////////////////////////////////////////////////////
  // specific methods to TAC array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    void Match_TAC();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TTACData object to TTACPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TTACData* rawDataPointer) {m_EventData = rawDataPointer;}
    
    void ClearTSRef(){RefTS_Name=""; RefTS=0;}
    
    // Setting TS Reference
    void SetRefTS(std::string TSRef_Name, unsigned long long TSRef);

    void ReadConfigurationTS(); 
  
  private:
    void ClaimReaderData();
    
  // objects are not written in the TTree
  private:
    TTACData*         m_EventData;        //!
    TTACData*         m_PreTreatedData;   //!
    TTACPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TTACData* GetRawData()        const {return m_EventData;}
    TTACData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    // thresholds
    double m_TAC_Time_RAW_Threshold = 0; //!
    unsigned int m_TAC_Mult; //!
    std::map<std::string,std::pair<unsigned int, unsigned long long>> SortTAC;//!

  // number of detectors
  private:
    int m_NumberOfDetectors;  //!

  // spectra class
  private:
    TTACSpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();
    static NPL::VTreeReader* ConstructReader();
  
  private:
    TimeStamp* TSEvent;//!
    std::string RefTS_Name = "";//!
    unsigned long long RefTS = 0;//!

    ClassDef(TTACPhysics,1)  // TACPhysics structure
};
#endif
