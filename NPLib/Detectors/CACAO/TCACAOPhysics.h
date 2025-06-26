#ifndef TCACAOPHYSICS_H
#define TCACAOPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: jwhwang@ibs.re.kr                        *
 *                                                                           *
 * Creation Date  : 4월 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CACAO Treated data                                *
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
#include "TRotation.h"

// NPTool headers
#include "TCACAOData.h"
#include "TCACAOSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TCACAOSpectra;



class TCACAOPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
public:
  TCACAOPhysics();
  ~TCACAOPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
public: 
  void Clear();   
  void Clear(const Option_t*) {};

  ////////////////////////////////////////////////////////////
  // 
  void AddDetector(TVector3 Pos, TRotation Rot,
		   TVector3 Dim, Double_t ShieldThickness);

  TVector3 GetDetPosition(int i) { return m_Pos[i]; }
  inline double GetNumberOfDetectors() { return m_NumberOfDetectors; }

  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
public:
  // read stream from ConfigFile to pick-up detector parameters
  void ReadConfiguration(NPL::InputParser);

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

  // methods related to the TCACAOSpectra class
  // instantiate the TCACAOSpectra class and 
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
  // specific methods to CACAO array
public:
  // remove bad channels, calibrate the data and apply thresholds
  void PreTreat();

  // clear the pre-treated object
  void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

  // read the user configuration file. If no file is found, load standard one
  void ReadAnalysisConfig();

  // give and external TCACAOData object to TCACAOPhysics. 
  // needed for online analysis for example
  void SetRawDataPointer(TCACAOData* rawDataPointer) {m_EventData = rawDataPointer;}
    
  // objects are not written in the TTree
private:
  TCACAOData*         m_EventData;        //!
  TCACAOData*         m_PreTreatedData;   //!
  TCACAOPhysics*      m_EventPhysics;     //!

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
public:
  Int_t nhit;
  Int_t detN[20];
  Double_t E[20], T[20];

  
  // getters for raw and pre-treated data object
public:
  TCACAOData* GetRawData()        const {return m_EventData;}
  TCACAOData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
private:
  // thresholds
  double m_E_RAW_Threshold; //!
  double m_E_Threshold;     //!

private:
  vector<TVector3> m_Pos;
  vector<TRotation> m_Rot;
  vector<TVector3> m_Dim;
  vector<double> m_ShieldThickness;
  
  // number of detectors
private:
  int m_NumberOfDetectors;  //!

  // spectra class
private:
  TCACAOSpectra* m_Spectra; // !

  // spectra getter
public:
  map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
public:
  static NPL::VDetector* Construct();

  ClassDef(TCACAOPhysics,1)  // CACAOPhysics structure
};
#endif
