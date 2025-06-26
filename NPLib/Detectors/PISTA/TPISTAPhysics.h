#ifndef TPISTAPHYSICS_H
#define TPISTAPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Pierre Morfouace  contact address: pierre.morfouace2@cea.fr                        *
 *                                                                           *
 * Creation Date  : May 2020                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold PISTA Treated data                                *
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
#include "TPISTAData.h"
#include "TPISTASpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TPISTASpectra;



class TPISTAPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TPISTAPhysics();
    ~TPISTAPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    vector<TVector2> Match_X_Y();
    int CheckEvent();

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    Int_t EventMultiplicity;
    Int_t m_E_NoMatching;
    Int_t m_DE_NoMatching;
    Int_t m_DE_E_NoMatching;
    vector<int>     DetectorNumber;
    vector<double>  E;
    vector<double>  DE;
    vector<double>  back_DE;
    vector<double>  back_E;
    vector<int>     DE_StripNbr;
    vector<int>     E_StripNbr;
    vector<double>  DE_Time;
    vector<double>  back_DE_Time;
    vector<double>  E_Time;
    vector<double>  back_E_Time;

    vector<double> PosX;
    vector<double> PosY;
    vector<double> PosZ;
    vector<int> mult_DE;
    vector<int> mult_E;

  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    /// A usefull method to bundle all operation to add a detector
    void AddDetector(TVector3 A, TVector3 B, TVector3 C, TVector3 D, int version); 
    void AddDetector(double R, double Theta, double Phi, int version); 
 
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

    // methods related to the TPISTASpectra class
    // instantiate the TPISTASpectra class and 
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
  // specific methods to PISTA array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();
  
    void InitializeStandardParameter();

    bool IsValidChannel(const int& DetectorType, const int& telescope, const int& channel);
    // give and external TPISTAData object to TPISTAPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TPISTAData* rawDataPointer) {m_EventData = rawDataPointer;}
 
    double GetNumberOfTelescope() const {return m_NumberOfDetectors;}
    int GetEventMultiplicity() const {return EventMultiplicity;}

    double GetStripPositionX(const int N, const int X, const int Y){
      return m_StripPositionX[N-1][X-1][Y-1];
    };
    double GetStripPositionY(const int N, const int X, const int Y){
      return m_StripPositionY[N-1][X-1][Y-1];
    };
    double GetStripPositionZ(const int N, const int X, const int Y){
      return m_StripPositionZ[N-1][X-1][Y-1];
    };


    TVector3 GetPositionOfInteraction(int DetectorNumber, int StripE, int StripDE);
    TVector3 GetDetectorNormal(const int i);
    
  // objects are not written in the TTree
  private:
    TPISTAData*         m_EventData;        //!
    TPISTAData*         m_PreTreatedData;   //!
    TPISTAPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TPISTAData* GetRawData()        const {return m_EventData;}
    TPISTAData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    // Map to activate/deactivate strips
    map<int, vector<bool>> m_DEChannelStatus; //!
    map<int, vector<bool>> m_EChannelStatus; //!

    int m_NumberOfDetectors; //!
    int m_NumberOfStripsX; //!
    int m_NumberOfStripsY; //!
    vector<TVector3> m_A; //!
    vector<TVector3> m_B; //!
    vector<TVector3> m_C; //!
    vector<TVector3> m_D; //!
    vector<vector<vector<double>>> m_StripPositionX; //!
    vector<vector<vector<double>>> m_StripPositionY; //!
    vector<vector<vector<double>>> m_StripPositionZ; //!

    int m_DistanceBetweenDEandE;

    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!
    double m_Back_E_Time_min;     //!
    double m_Back_E_Time_max;     //!
    int m_CALIB_BACK_E_PER_STRIP; //!

 private:
    unsigned int m_MaximumStripMultiplicityAllowed;//
    double m_StripEnergyMatching;//


  // spectra class
  private:
    TPISTASpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TPISTAPhysics,1)  // PISTAPhysics structure
};
#endif
