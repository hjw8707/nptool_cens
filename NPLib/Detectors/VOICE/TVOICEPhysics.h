#ifndef TVOICEPHYSICS_H
#define TVOICEPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2023   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sunghan Bae  contact address: shbae2703@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : July 2023                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold VOICE Treated data                                *
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
#include "TVOICEData.h"
//#include "TVOICEGasData.h"
//#include "TVOICESpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
//class TVOICESpectra;



class TVOICEPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TVOICEPhysics();
    ~TVOICEPhysics();


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
  /// A usefull method to bundle all operation to add a detector
  void AddDetector(string Type, double zpos, double Tilt, double Rot, double radi, double pitch); 
  //G4Material* SetGasType(string gas, double pressure);
  
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

    // methods related to the TVOICESpectra class
    // instantiate the TVOICESpectra class and 
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
  // specific methods to VOICE array
  public:
    // remove bad channels, calibrate the data and apply thresholds
//    void PreTreat();

    // clear the pre-treated object
//    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
//    void ReadAnalysisConfig();

    // give and external TVOICEData object to TVOICEPhysics. 
    // needed for online analysis for example
//    void SetRawDataPointer(TVOICEData* rawDataPointer) {m_EventData = rawDataPointer;}
    
  // objects are not written in the TTree
  private:
    TVOICEData*         m_EventData;        //!
//    TVOICEGasData*      m_EventGasData;        //!
    TVOICEData*         m_PreTreatedData;   //!
    TVOICEPhysics*      m_EventPhysics;     //!
					    

    //Configuration reading array
    vector<string> m_Type;
    vector<double> m_zpos;
    vector<double> m_Tilt;
    vector<double> m_Rot;
    vector<double> m_radi;
    vector<double> m_pitch;

    string m_gas;
    double m_pressure;

    double m_AC_gap;
    double m_AFG_gap;

    //tree branches
    Int_t nhit;
    Int_t type[50];
    Double_t E[50];
    Double_t X[50];
    Double_t Y[50];
    Double_t Z[50];
    Double_t Time[50];
    Double_t Theta[50];
    Double_t Phi[50];
    Int_t Atom_Z[50];
    Int_t A[50];
    Int_t PDG[50];
    
    Double_t E_f[50];
    Double_t X_f[50];
    Double_t Y_f[50];
    Double_t Z_f[50];
    Double_t Time_f[50];
    Double_t Theta_f[50];
    Double_t Phi_f[50];
					    

  // getters for raw and pre-treated data object
//  public:
//    TVOICEData* GetRawData()        const {return m_EventData;}
  //  TVOICEData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  //private:
    // thresholds
    //double m_E_RAW_Threshold; //!
    //double m_E_Threshold;     //!

  // number of detectors
//  private:
  //  int m_NumberOfDetectors;  //!

  // spectra class
//  private:
  //  TVOICESpectra* m_Spectra; // !

  // spectra getter
 // public:
   // map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TVOICEPhysics,1)  // VOICEPhysics structure
};
#endif
