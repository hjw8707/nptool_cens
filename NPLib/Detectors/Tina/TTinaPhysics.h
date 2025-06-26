#ifndef TTinaPHYSICS_H
#define TTinaPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2018   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Serge Franchoo  contact address: franchoo@ipno.in2p3.fr                        *
 *                                                                           *
 * Creation Date  : February 2018                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Tina Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *   
 *                                                                           *
 *****************************************************************************/

#include <vector>
#include <map>
#include <string>
#include "TObject.h"
#include "TH1.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TTinaData.h"
#include "TTinaSpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"

using namespace std;

// forward declaration
class TTinaSpectra;

class TTinaPhysics : public TObject, public NPL::VDetector {
  // constructor and destructor
  public:
    TTinaPhysics();
    ~TTinaPhysics();

  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    vector<TVector2> Match_X_Y();
    void             CheckEvent(int N);

  // data obtained after BuildPhysicalEvent() and stored in output Root file
  public:
    int EventMultiplicity;
    vector<int>    EventType;
    vector<int>    SquareTelescopeNumber;
    vector<int>    TrapezTelescopeNumber;
    vector<double> TTT_E;
    vector<double> TTT_T;
    vector<int>    TTT_X;
    vector<int>    TTT_Y;
    //vector<double> TTT_EX;
    //vector<double> TTT_TX;
    //vector<double> TTT_EY;
    //vector<double> TTT_TY;
    //vector<int>    TelescopeNumber_X;
    //vector<int>    TelescopeNumber_Y;
    vector<double> Pad_E;
    vector<double> Pad_T;
    vector<int>    Pad_N;

    vector<double> YY1_E;
    vector<double> YY1_T;
    vector<int>    YY1_R;
    vector<int>    YY1_S;
    vector<double> CsI_E;
    vector<double> CsI_T;
    vector<int>    CsI_N;
    
    // physical value
    vector<double> Energy;
    vector<double> Time;
    
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
    // create branches of output Root file
    void InitializeRootOutput();
    // clear the raw and physical data objects event by event
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}   
    // methods related to the TTinaSpectra class
    // instantiate the TTinaSpectra class and declare list of histograms
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

  // specific methods to Tina
  public:
    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();
    // return false if the channel is disabled by user
    // first argument is either 0 for X,1 Y,2 SiLi, 3 CsI
    bool IsValidChannel(const int& DetectorType, const int& telescope,const int& channel);
    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();
    // method to bundle all operation to add a detector
    void AddDetector(TVector3 POS, string shape, double Translation);
    void AddDetector(double R, double Theta, double Phi, double Translation, string shape); 
    // give and external TTinaData object to TTinaPhysics
    // needed for online analysis for example
    void SetRawDataPointer(TTinaData* rawDataPointer) {m_EventData = rawDataPointer;}
    // getters for raw and pre-treated data object
    TTinaData* GetRawData()        const {return m_EventData;}
    TTinaData* GetPreTreatedData() const {return m_PreTreatedData;}
    // use to access the strip position
    double GetStripPositionX(const int N, const int X, const int Y) const {
      return m_StripPositionX[N - 1][X - 1][Y - 1];
    };
    double GetStripPositionY(const int N, const int X, const int Y) const {
      return m_StripPositionY[N - 1][X - 1][Y - 1];
    };
    double GetStripPositionZ(const int N, const int X, const int Y) const {
      return m_StripPositionZ[N - 1][X - 1][Y - 1];
    };
    double GetRingPositionX(const int N, const int R) const {
        return m_RingPositionX[N - m_NumberOfTTTPad - 1][R - 1];
    };
    double GetRingPositionY(const int N, const int R) const {
        return m_RingPositionY[N - m_NumberOfTTTPad - 1][R - 1];
    };
    double GetRingPositionZ(const int N, const int R) const {
        return m_RingPositionZ[N - m_NumberOfTTTPad - 1][R - 1];
    };
    int GetNumberOfTelescope() const { return m_NumberOfTelescopes; };
    int GetNumberOfTTTPad() const { return m_NumberOfTTTPad; };
    // to be called after a BuildPhysicalEvent
    int GetEventMultiplicity() const { return EventMultiplicity; };
    double GetEnergyDeposit(const int i) const { return Energy[i]; };
    TVector3 GetPositionOfInteraction(const int i) const;
    TVector3 GetTelescopeNormal(const int i) const;
    TVector3 GetRingPositionOfInteraction(const int i) const;
    TVector3 GetRingNormal(const int i) const;

  // parameters used in the analysis
  private: 
    // container size to be used in the analysis loops (value must be given locally)
    unsigned int m_StripXEMult; //!
    unsigned int m_StripYEMult; //!
    unsigned int m_StripXTMult; //!
    unsigned int m_StripYTMult; //!
    unsigned int m_PadEMult; //!
    unsigned int m_PadTMult; //!
    
    unsigned int m_YY1RingEMult; //!
    unsigned int m_YY1SectorEMult; //!
    unsigned int m_YY1RingTMult; //!
    unsigned int m_YY1SectorTMult; //!
    unsigned int m_CsIEMult; //!
    unsigned int m_CsITMult; //!
    // events over this value after pre-treatment are not treated to avoid long
    // treatment times of spurious events
    unsigned int m_MaximumStripMultiplicityAllowed; //!
    // allowance in percent of the difference in energy between X and Y
    double m_StripEnergyMatchingSigma; //!
    double m_StripEnergyMatchingNumberOfSigma; //!
    int m_TTT_XE_RAW_Threshold; //!
    int m_TTT_YE_RAW_Threshold; //!
    int m_Pad_E_RAW_Threshold; //!
    double m_TTT_XE_CAL_Threshold; //!
    double m_TTT_YE_CAL_Threshold; //!
    double m_Pad_E_CAL_Threshold; //!
    
    int m_YY1_RE_RAW_Threshold; //!
    int m_YY1_SE_RAW_Threshold; //!
    int m_CsI_E_RAW_Threshold; //!
    double m_YY1_RE_CAL_Threshold; //!
    double m_YY1_SE_CAL_Threshold; //!
    double m_CsI_E_CAL_Threshold; //!
    // geometric matching ignored for the time being, should come here
    // if true, all events that do not come in front of a crystal will be ignored (crossing or not)
    // this option reduces statistics, however it helps eliminating
    // unwanted events that cross the DSSD and pass between the crystals
    // bool m_Ignore_not_matching_CsI; //!
    // by default take EX and TY
    bool m_Take_YE; //!
    bool m_Take_YT; //!
    bool m_Take_SE; //!
    bool m_Take_ST; //!
    int m_Pad_Size; //! 
    int m_NumberOfStrips; //!
    int m_NumberOfTTTPad;
    double m_Face; //!
    int m_NumberOfRings; //!
    int m_NumberOfYY1CsI;
    double m_Face_Ring; //!

  private:
   // Root input and output tree classes
    TTinaData*         m_EventData;        //!
    TTinaData*         m_PreTreatedData;   //!

    TTinaPhysics*      m_EventPhysics;     //!

  private:
    int m_NumberOfTelescopes;  //!
    vector<vector<vector<double>>> m_StripPositionX; //!
    vector<vector<vector<double>>> m_StripPositionY; //!
    vector<vector<vector<double>>> m_StripPositionZ; //!
    vector<vector<double>> m_RingPositionX; //!
    vector<vector<double>> m_RingPositionY; //!
    vector<vector<double>> m_RingPositionZ; //!

  public:
    // prevent to treat event with ambiguous matching between X and Y
    // bool          m_multimatch; //!
    int           m_OrderMatch; //!
    vector<int>   m_match_type; //!
    map<int, int> m_NMatchDet; //!
    map<int, int> m_StripXMultDet; //!
    map<int, int> m_StripYMultDet; //!

  private:
    TTinaSpectra* m_Spectra; // !

  public:
    map<string,TH1*> GetSpectra(); 

  public: 
    // static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    // TinaPhysics structure
    ClassDef(TTinaPhysics,1)  
};

namespace TINA_LOCAL {

  // calibration functions
  double fTTT_XE(const TTinaData* Data, const int& i);
  double fTTT_XT(const TTinaData* Data, const int& i);
  double fTTT_YE(const TTinaData* Data, const int& i);
  double fTTT_YT(const TTinaData* Data, const int& i);
  double fPad_E(const TTinaData* Data, const int& i);
  double fPad_T(const TTinaData* Data, const int& i);
   
  double fYY1_RE(const TTinaData* Data, const int& i);
  double fYY1_RT(const TTinaData* Data, const int& i);
  double fYY1_SE(const TTinaData* Data, const int& i);
  double fYY1_ST(const TTinaData* Data, const int& i);
  double fCsI_E(const TTinaData* Data, const int& i);
  double fCsI_T(const TTinaData* Data, const int& i);
} 

#endif
