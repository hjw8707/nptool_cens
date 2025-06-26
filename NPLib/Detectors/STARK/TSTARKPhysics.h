#ifndef TSTARKPHYSICS_H
#define TSTARKPHYSICS_H

/*****************************************************************************
 * Original Author: Jongwon Hwang    contact address: jwhwang@ibs.re.kr      * 
 * Creation Date  : May 2023                                                 *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold STARK Treated data                                *
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
#include "TSTARKData.h"
#include "TClonesArray.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"

using namespace std;


class TSTARKPhysics: public TObject, public NPL::VDetector {
  // constructor and destructor
  public:
    TSTARKPhysics();
    ~TSTARKPhysics();

  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    void AddDetector(string Type, TVector3 Pos, int Flip, int Rev, double Beta);
    void AddStripPosition(string Type, TVector3 Pos, int Flip, int Rev, double Beta);

    // data obtained after BuildPhysicalEvent() and stored in output Root file
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
    void BuildSimplePhysicalEvent() {BuildPhysicalEvent();};
    // same as above but for online analysis
    void BuildOnlinePhysicalEvent() {BuildPhysicalEvent();};
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

    // methods related to the TSTARKSpectra class
    // instantiate the TSTARKSpectra class and declare list of histograms
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

  // parameters used in the analysis
  private:
   // Root input and output tree classes
    TSTARKData*         m_EventData;        //!
    TSTARKData*         m_PreTreatedData;   //!
    TSTARKPhysics*      m_EventPhysics;     //!

    ////////////////////////////////////////////////////////////
    // detector geometry parameters
    //
    // for X6
    double X6_activeX, X6_activeY; // mm
    int X6_NFrontStrips, X6_NBackStrips;
    //
    // for BB10
    double BB10_activeX, BB10_activeY; // mm
    int BB10_NFrontStrips, BB10_NBackStrips;
    //
    // for QQQ5
    double QQQ5_activeInR, QQQ5_activeOutR; // mm
    int QQQ5_NRStrips, QQQ5_NAStrips;

    ////////////////////////////////////////////////////////////
    // configuration
    vector<string> m_Type; // X6, BB10, ...
    vector<TVector3> m_DetPos; // detector center position
    vector<int> m_Flip; // flip of the surface (for all)
    vector<int> m_Rev; // flip of the long-side direction (only for X6 & BB10)
    vector<double> m_Beta; // rotation along the beam axis (only for QQQ5)
    vector<TVector3> m_DetOri; // detector Y axis direction (only for X6)
    ////////////////////////////////////////////////////////////
    // strip center position
    // m_StripPos[detN][frontStripN][backStripN]
    // -> Pad case: only [0] exist (but type should be std::vector)
    vector<vector<vector<TVector3>>> m_StripPos; 
    ////////////////////////////////////////////////////////////

 public:
     ////////////////////////////////////////////////////////////
    // tree branches
    Int_t nhit;
    Int_t type[20]; // type = 0 (X6), 1 (BB10)
    Int_t detN[20], fStrN[20], bStrN[20];
    Double_t uppE[20], dwnE[20], sumE[20];
    std::vector<TVector3> sPosArr, hPosArr;
    ////////////////////////////////////////////////////////////
    
  public: 
    // static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    // STARKPhysics structure
    ClassDef(TSTARKPhysics,1)  
};

#endif
