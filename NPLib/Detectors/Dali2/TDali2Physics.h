#ifndef TDali2PHYSICS_H
#define TDali2PHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Elidiano Tronchin  contact address: elidiano.tronchin@studenti.unipd.it                        *
 *                                                                           *
 * Creation Date  : septembre 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Dali2 Treated data                                *
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
#include "TDali2Data.h"
#include "TDali2Spectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TDali2Spectra;


class TDali2Physics: public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
 public:
  TDali2Physics();
  ~TDali2Physics() {};


 public:
  TVector3 GetDALIPosition(int DALINbr) {
    static TVector3 temp;
    temp.SetXYZ(m_R[DALINbr]*cos(m_Alpha[DALINbr]),
		m_R[DALINbr]*sin(m_Alpha[DALINbr]),
		m_Zeta[DALINbr]);
    return temp; }
  
  inline int GetNumberOfDetectors() { return m_NumberOfDetectors; }

  double GetDopplerCorrectedEnergy(double energy,
  				   TVector3 position,
  				   TVector3 beta);
  //!
  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
 public: 
  void Clear();   
  void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
 public:
  vector<int>      DetectorNumber;
  vector<double>   ADC;
  vector<double>   Energy;
  vector<double>   TDC;
  vector<double>   Time;
  /* vector<int>   ParticleID; */

  /// A usefull method to bundle all operation to add a detector
  //void AddDetector(TVector3 POS, string shape);
  void AddDetector(TVector3 POS); 
  void AddDetector(double R, double Theta, double Phi);//, string shape);
  void AddDetector(double R, double Alpha, double Zeta, int Ring); 
  
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

  // methods related to the TDaliSpectra class
  // instantiate the TDaliSpectra class and 
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
  // specific methods to Dali array
 public:
  // remove bad channels, calibrate the data and apply thresholds
  void PreTreat();

  // clear the pre-treated object
  void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

  // read the user configuration file. If no file is found, load standard one
  void ReadAnalysisConfig();

  // give and external TDaliData object to TDaliPhysics. 
  // needed for online analysis for example
  void SetRawDataPointer(TDali2Data* rawDataPointer) {m_EventData = rawDataPointer;}
    
  // objects are not written in the TTree
 private:
  TDali2Data*         m_EventData;        //!
  TDali2Data*         m_PreTreatedData;   //!
  TDali2Physics*      m_EventPhysics;     //!
  
  // getters for raw and pre-treated data object
 public:
  TDali2Data* GetRawData()        const {return m_EventData;}
  TDali2Data* GetPreTreatedData() const {return m_PreTreatedData;}

  // Use to access energies et det number

  int GetDetNumber(const int &i) const {return DetectorNumber[i];};
  int GetMultEnergy() const {return Energy.size();};
  double GetEnergy(const int &i) const {return Energy[i];};
  /* int GetParticleID(const int &i) const {return ParticleID[i];}; */

  // for tree output
  Int_t nhit;
  Int_t detN[200];
  Double_t energy[200];
  Double_t energyDC[200];
    

  // parameters used in the analysis
 private:
  // thresholds
  double m_E_RAW_Threshold; //!
  double m_E_Threshold;     //!

  // number of detectors
  private:
  int m_NumberOfDetectors;  //!

  // detector geometry
 public:
  vector<double>  m_Zeta;
  vector<double>  m_R; 
  vector<double>  m_Alpha; 
  vector<int>  m_Ring; 
    
  // spectra class
 private:
  TDali2Spectra* m_Spectra; // !

  // spectra getter
 public:
  map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
 public:
  static NPL::VDetector* Construct();

  ClassDef(TDali2Physics,1)  // DaliPhysics structure
    };
#endif
