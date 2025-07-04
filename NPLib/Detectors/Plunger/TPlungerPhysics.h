#ifndef TPlungerPHYSICS_H
#define TPlungerPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2025   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Jongwon Hwang  contact address: hjw8707@gmail.com                        *
 *                                                                           *
 * Creation Date  : 7월 2025                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold Plunger Treated data                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// C++ headers
#include <map>
#include <string>
#include <vector>
using namespace std;

// ROOT headers
#include "TH1.h"
#include "TObject.h"
#include "TVector3.h"
// NPTool headers
#include "NPCalibrationManager.h"
#include "NPInputParser.h"
#include "NPVDetector.h"
#include "TPlungerData.h"
// forward declaration
class TPlungerSpectra;

class TPlungerPhysics : public TObject, public NPL::VDetector {
    //////////////////////////////////////////////////////////////
    // constructor and destructor
   public:
    TPlungerPhysics();
    ~TPlungerPhysics() {};

    //////////////////////////////////////////////////////////////
    // Inherited from TObject and overriden to avoid warnings
   public:
    void Clear();
    void Clear(const Option_t*) {};

    //////////////////////////////////////////////////////////////
    // data obtained after BuildPhysicalEvent() and stored in
    // output ROOT file
   public:
    vector<int> DetectorNumber;
    vector<double> Energy;
    vector<double> Time;

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
    void BuildOnlinePhysicalEvent() { BuildPhysicalEvent(); };

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
    void ClearEventPhysics() { Clear(); }
    void ClearEventData() { m_EventData->Clear(); }

    //////////////////////////////////////////////////////////////
    // specific methods to Plunger array
   public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData() { m_PreTreatedData->Clear(); }

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TPlungerData object to TPlungerPhysics.
    // needed for online analysis for example
    void SetRawDataPointer(TPlungerData* rawDataPointer) { m_EventData = rawDataPointer; }

    // objects are not written in the TTree
   private:
    TPlungerData* m_EventData;        //!
    TPlungerData* m_PreTreatedData;   //!
    TPlungerPhysics* m_EventPhysics;  //!

    // getters for raw and pre-treated data object
   public:
    TPlungerData* GetRawData() const { return m_EventData; }
    TPlungerData* GetPreTreatedData() const { return m_PreTreatedData; }

    // parameters used in the analysis
   private:
    // thresholds
    double m_E_RAW_Threshold;  //!
    double m_E_Threshold;      //!

    // number of detectors
   private:
    int m_NumberOfDetectors;  //!

    // parameters for the target
   private:
    bool m_TargetFound;        //!
    double m_TargetR;          //!
    double m_TargetThickness;  //!
    double m_TargetPosZ;       //!
    string m_TargetMaterial;   //!

    // parameters for the stopper
   private:
    bool m_StopperFound;        //!
    double m_StopperR;          //!
    double m_StopperThickness;  //!
    double m_StopperPosZ;       //!
    string m_StopperMaterial;   //!

    // parameters for the chamber
   private:
    bool m_ChamberFound;        //!
    double m_ChamberR;          //!
    double m_ChamberThickness;  //!
    string m_ChamberMaterial;   //!
    double m_ChamberPipeR;      //!
    double m_ChamberPipeZ0;     //!
    double m_ChamberPipeZ1;     //!

    // Static constructor to be passed to the Detector Factory
   public:
    static NPL::VDetector* Construct();

    ClassDef(TPlungerPhysics, 1)  // PlungerPhysics structure
};
#endif
