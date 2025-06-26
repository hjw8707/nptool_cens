#ifndef TICPHYSICS_H
#define TICPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2020   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: P. Morfouace  contact address: pierre.morfouace@cea.fr   *
 *                                                                           *
 * Creation Date  : October 2023                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold IC Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

// C++ headers 
#include <vector>
#include <map>
#include <set>
#include <string>
using namespace std;

// ROOT headers
#include "TObject.h"
#include "TH1.h"
#include "TVector3.h"
#include "TSpline.h"
#include "TFile.h"
// NPTool headers
#include "TICData.h"
#include "TTimeData.h"
#include "TProfileEvaluator.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"



class TICPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TICPhysics();
    ~TICPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};


  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    double DE;
    double Eres;
    double Etot;
    double EtotInit;
    double Chio_Z;
    double Chio_ZRaw;
    
    vector<long> fIC_TS;

    double fIC[11];
    double fIC_raw[11]; 
    double fIC_PID[11]; 
    double fIC_Init[11];//! 
    
    // For AoQ the order in ToF is 13 23 14 24
    double EtotAoQ[4] ;
    vector<array<double,11>> fIC_AoQ =vector<array<double,11>>(4);//!
  private:

    /// A usefull method to bundle all operation to add a detector
  void AddDetector(TVector3 Pos); 
  void AddDetector(double R, double Theta, double Phi); 
  
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


  //////////////////////////////////////////////////////////////
  // specific methods to IC array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TICData object to TICPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TICData* rawDataPointer) {m_EventData = rawDataPointer;}

    // Setter and getter for section of the FPMW.
    void SetFPMWSection(int section) {m_FPMW_Section = section;}
    int GetFPMWSection() {return m_FPMW_Section;}
    
    //Setter and Getter for TTimeData used for drift time calculation and XY
    void SetTTimeData(TTimeData *newTime) {m_TimeData = newTime;}
    TTimeData* GetTTimeData() {return m_TimeData;}
    
    void SetX(double PosX) {m_X = PosX;}
    double GetX() {return m_X;}
    void SetY(double PosY) {m_Y = PosY;}
    double GetY() {return m_Y;}
    void SetThetaf(double Thetaf) {m_Thetaf = Thetaf;}
    double GetThetaf() {return m_Thetaf;}


    bool LoadSpline(vector<TSpline3*> &iSpline, int &NSpline, string Path);
    template <typename T> bool LoadVector(vector<T> &vec, const string &Path);
    double ApplyZSpline();
    void ApplyXYCorrections();

  // objects are not written in the TTree
  private:
    TICData*         m_EventData;        //!
    TICData*         m_PreTreatedData;   //!
    TICPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TICData* GetRawData()        const {return m_EventData;}
    TICData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    int m_NumberOfDetectors;  //!
    int m_FPMW_Section; //!
    int m_number_zspline; //!
                          
    string m_Z_SPLINE_PATH; //!
    string m_Z_SPLINE_EVAL_PATH; //!
                                 
    bool m_Z_SPLINE_CORRECTION = false; //!
    bool m_Z_SPLINE_EVAL = false; //!
                                  
    vector<TSpline3*>  m_Zspline = vector<TSpline3*>(50); //!
    vector<Double_t>  m_Z_spline_eval ; //!
                                        
    double m_Eres_Threshold; //!

    //Time temporary variable
    TTimeData*      m_TimeData;  //!
    //Correction in IC
    string m_Y_SPLINE_PATH; //!
    string m_XY0_PROFILE_PATH; //!
    string m_DE_SPLINE_PATH; //!
    string m_Z_THETA_SPLINE_PATH; //!

    bool m_Y_SPLINE_CORRECTION = false; //!
    bool m_DE_SPLINE_CORRECTION = false; //!
    bool m_XY0_SPLINE_CORRECTION = false; //!
    bool m_Z_THETA_CORRECTION = false; //!

    int m_number_Y_spline; //!
    int m_number_XY0_spline; //!
    int m_number_DE_spline; //!
    int m_number_Z_THETA_spline; //!

    double m_X; //!
    double m_Y; //!
    double m_Thetaf; //!
    // Be careful this array contains only the correction from IC1 to 4
    vector<TSpline3*> m_Yspline = vector<TSpline3*>(11); //! 
    vector<TSpline3*> m_DEspline = vector<TSpline3*>(11); //! 
    vector<TSpline3*> m_Z_Theta_spline = vector<TSpline3*>(2); //! 
    // The following contains correction for IC0
    ProfileEvaluator m_IC0_Profile;//!

    // ToF Name for loading
    vector<string> ToFName {"13","23","14","24"}; //!
  private:
    //Year of acquisition
    double m_Data_Year; //!
  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TICPhysics,1)  // ICPhysics structure
};
#endif
