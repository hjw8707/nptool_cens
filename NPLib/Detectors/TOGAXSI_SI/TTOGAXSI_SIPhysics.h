#ifndef TTOGAXSI_SIPHYSICS_H
#define TTOGAXSI_SIPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2024   this file is part of the NPTool Project       *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Thomas Pohl  contact address: thomas.pohl@riken.jp                        *
 *                                                                           *
 * Creation Date  : February 2024                                           *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold TOGAXSI_SI Treated data                                *
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
#include "TTOGAXSI_SIData.h"
#include "TTOGAXSI_SISpectra.h"
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPInputParser.h"
// forward declaration
class TTOGAXSI_SISpectra;



class TTOGAXSI_SIPhysics : public TObject, public NPL::VDetector {
  //////////////////////////////////////////////////////////////
  // constructor and destructor
  public:
    TTOGAXSI_SIPhysics();
    ~TTOGAXSI_SIPhysics() {};


  //////////////////////////////////////////////////////////////
  // Inherited from TObject and overriden to avoid warnings
  public: 
    void Clear();   
    void Clear(const Option_t*) {};

  public:
    vector<TVector2> MatchInner();
    vector<TVector2> MatchOuter();
    vector<TVector2> MatchCluster1();
    vector<TVector2> MatchCluster2();
    int CheckEvent();

  //////////////////////////////////////////////////////////////
  // data obtained after BuildPhysicalEvent() and stored in
  // output ROOT file
  public:
    Int_t EventMultiplicity;
    Int_t EventMultiplicity_cluster;
    vector<int>      DetectorNumber;
    vector<int>      DetectorNumberInnerX;
    vector<int>      DetectorNumberInnerZ;
    vector<int>      DetectorNumberOuterX;
    vector<int>      DetectorNumberOuterZ;
    vector<int>      DetectorNumberClusterX1;
    vector<int>      DetectorNumberClusterY1;
    vector<int>      DetectorNumberClusterX2;
    vector<int>      DetectorNumberClusterY2;

    vector<int>	     InnerXStrip;
    vector<double>   InnerXE;
    vector<int>	     InnerZStrip;
    vector<double>   InnerZE;
    vector<int>	     OuterXStrip;
    vector<double>   OuterXE;
    vector<int>	     OuterZStrip;
    vector<double>   OuterZE;
    vector<int>	     ClusterInnerStrip;
    vector<double>   ClusterInnerE;
    vector<int>	     ClusterX1Strip;
    vector<double>   ClusterX1E;
    vector<int>	     ClusterY1Strip;
    vector<double>   ClusterY1E;
    vector<int>	     ClusterX2Strip;
    vector<double>   ClusterX2E;
    vector<int>	     ClusterY2Strip;
    vector<double>   ClusterY2E;
    
    vector<double>   InnerXPosX;
    vector<double>   InnerXPosY;
    vector<double>   InnerXPosZ;
    vector<double>   InnerZPosX;
    vector<double>   InnerZPosY;
    vector<double>   InnerZPosZ;
    vector<double>   OuterXPosX;
    vector<double>   OuterXPosY;
    vector<double>   OuterXPosZ;
    vector<double>   OuterZPosX;
    vector<double>   OuterZPosY;
    vector<double>   OuterZPosZ;

    vector<double>   ClusterInnerPosX;
    vector<double>   ClusterInnerPosY;
    vector<double>   ClusterInnerPosZ;
    vector<double>   ClusterX1PosX;
    vector<double>   ClusterX1PosY;
    vector<double>   ClusterX1PosZ;
    vector<double>   ClusterY1PosX;
    vector<double>   ClusterY1PosY;
    vector<double>   ClusterY1PosZ;
    vector<double>   ClusterX2PosX;
    vector<double>   ClusterX2PosY;
    vector<double>   ClusterX2PosZ;
    vector<double>   ClusterY2PosX;
    vector<double>   ClusterY2PosY;
    vector<double>   ClusterY2PosZ;

  
  //////////////////////////////////////////////////////////////
  // methods inherited from the VDetector ABC class
  public:
    // read stream from ConfigFile to pick-up detector parameters
    void ReadConfiguration(NPL::InputParser);

    //A usefull method to bundle all operation to add a detector
    void AddInnerXDetector(double R, double Z, double Phi, TVector3 Ref);
    void AddInnerZDetector(double R, double Z, double Phi, TVector3 Ref);
    void AddOuterXDetector(double R, double Z, double Phi, TVector3 Ref);
    void AddOuterZDetector(double R, double Z, double Phi, TVector3 Ref);
    void AddClusterInnerDetector(double R, double Z, double Phi, TVector3 Ref);
    void AddClusterX1Detector(double R, double Z, double Phi, TVector3 Ref);
    void AddClusterY1Detector(double R, double Z, double Phi, TVector3 Ref);
    void AddClusterX2Detector(double R, double Z, double Phi, TVector3 Ref);
    void AddClusterY2Detector(double R, double Z, double Phi, TVector3 Ref);

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

    // methods related to the TTOGAXSI_SISpectra class
    // instantiate the TTOGAXSI_SISpectra class and 
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
  // specific methods to TOGAXSI_SI array
  public:
    // remove bad channels, calibrate the data and apply thresholds
    void PreTreat();

    // clear the pre-treated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}

    // read the user configuration file. If no file is found, load standard one
    void ReadAnalysisConfig();

    // give and external TTOGAXSI_SIData object to TTOGAXSI_SIPhysics. 
    // needed for online analysis for example
    void SetRawDataPointer(TTOGAXSI_SIData* rawDataPointer) {m_EventData = rawDataPointer;}
   
    double GetNumberOfInnerXDetector() const {return m_NumberOfInnerXDetectors;}
    double GetNumberOfInnerZDetector() const {return m_NumberOfInnerZDetectors;}
    double GetNumberOfOuterXDetector() const {return m_NumberOfOuterXDetectors;}
    double GetNumberOfOuterZDetector() const {return m_NumberOfOuterZDetectors;}
    double GetNumberOfClusterInnerDetector() const {return m_NumberOfClusterInnerDetectors;}
    double GetNumberOfClusterX1Detector() const {return m_NumberOfClusterX1Detectors;}
    double GetNumberOfClusterY1Detector() const {return m_NumberOfClusterY1Detectors;}
    double GetNumberOfClusterX2Detector() const {return m_NumberOfClusterX2Detectors;}
    double GetNumberOfClusterY2Detector() const {return m_NumberOfClusterY2Detectors;}


    int GetEventMultiplicity() const {return EventMultiplicity;}
    int GetEventMultiplicity_cluster() const {return EventMultiplicity_cluster;}
    int GetDetNbrSize() const {return DetectorNumber.size();}
    TVector3 GetInnerXPos(const int i) const {return TVector3(InnerXPosX[i],InnerXPosY[i],InnerXPosZ[i]);}
    TVector3 GetInnerZPos(const int i) const {return TVector3(InnerZPosX[i],InnerZPosY[i],InnerZPosZ[i]);}
    TVector3 GetOuterXPos(const int i) const {return TVector3(OuterXPosX[i],OuterXPosY[i],OuterXPosZ[i]);}
    TVector3 GetOuterZPos(const int i) const {return TVector3(OuterZPosX[i],OuterZPosY[i],OuterZPosZ[i]);}
    TVector3 GetClusterInnerPos(const int i) const {return TVector3(ClusterInnerPosX[i],ClusterInnerPosY[i],ClusterInnerPosZ[i]);}
    TVector3 GetClusterX1Pos(const int i) const {return TVector3(ClusterX1PosX[i],ClusterX1PosY[i],ClusterX1PosZ[i]);}
    TVector3 GetClusterY1Pos(const int i) const {return TVector3(ClusterY1PosX[i],ClusterY1PosY[i],ClusterY1PosZ[i]);}
    TVector3 GetClusterX2Pos(const int i) const {return TVector3(ClusterX2PosX[i],ClusterX2PosY[i],ClusterX2PosZ[i]);}
    TVector3 GetClusterY2Pos(const int i) const {return TVector3(ClusterY2PosX[i],ClusterY2PosY[i],ClusterY2PosZ[i]);}

    double GetInnerXStripPositionX(const int N, const int X, const int Y) {
      return m_InnerXStripPositionX[N-1][X-1][Y-1];
    };
    double GetInnerXStripPositionY(const int N, const int X, const int Y) {
      return m_InnerXStripPositionY[N-1][X-1][Y-1];
    };
    double GetInnerXStripPositionZ(const int N, const int X, const int Y) {
      return m_InnerXStripPositionZ[N-1][X-1][Y-1];
    };

    double GetInnerZStripPositionX(const int N, const int X, const int Y) {
      return m_InnerZStripPositionX[N-1][X-1][Y-1];
    };
    double GetInnerZStripPositionY(const int N, const int X, const int Y) {
      return m_InnerZStripPositionY[N-1][X-1][Y-1];
    };
    double GetInnerZStripPositionZ(const int N, const int X, const int Y) {
      return m_InnerZStripPositionZ[N-1][X-1][Y-1];
    };


    double GetOuterXStripPositionX(const int N, const int X, const int Y) {
      return m_OuterXStripPositionX[N-1][X-1][Y-1];
    };
    double GetOuterXStripPositionY(const int N, const int X, const int Y) {
      return m_OuterXStripPositionY[N-1][X-1][Y-1];
    };
    double GetOuterXStripPositionZ(const int N, const int X, const int Y) {
      return m_OuterXStripPositionZ[N-1][X-1][Y-1];
    };

    double GetOuterZStripPositionX(const int N, const int X, const int Y) {
      return m_OuterZStripPositionX[N-1][X-1][Y-1];
    };
    double GetOuterZStripPositionY(const int N, const int X, const int Y) {
      return m_OuterZStripPositionY[N-1][X-1][Y-1];
    };
    double GetOuterZStripPositionZ(const int N, const int X, const int Y) {
      return m_OuterZStripPositionZ[N-1][X-1][Y-1];
    };


    double GetClusterInnerStripPositionX(const int N, const int X, const int Y) {
      return m_ClusterInnerStripPositionX[N-1][X-1][Y-1];
    };
    double GetClusterInnerStripPositionY(const int N, const int X, const int Y) {
      return m_ClusterInnerStripPositionY[N-1][X-1][Y-1];
    };
    double GetClusterInnerStripPositionZ(const int N, const int X, const int Y) {
      return m_ClusterInnerStripPositionZ[N-1][X-1][Y-1];
    };

    double GetClusterX1StripPositionX(const int N, const int X, const int Y) {
      return m_ClusterX1StripPositionX[N-1][X-1][Y-1];
    };
    double GetClusterX1StripPositionY(const int N, const int X, const int Y) {
      return m_ClusterX1StripPositionY[N-1][X-1][Y-1];
    };
    double GetClusterX1StripPositionZ(const int N, const int X, const int Y) {
      return m_ClusterX1StripPositionZ[N-1][X-1][Y-1];
    };

    double GetClusterY1StripPositionX(const int N, const int X, const int Y) {
      return m_ClusterY1StripPositionX[N-1][X-1][Y-1];
    };
    double GetClusterY1StripPositionY(const int N, const int X, const int Y) {
      return m_ClusterY1StripPositionY[N-1][X-1][Y-1];
    };
    double GetClusterY1StripPositionZ(const int N, const int X, const int Y) {
      return m_ClusterY1StripPositionZ[N-1][X-1][Y-1];
    };

    double GetClusterX2StripPositionX(const int N, const int X, const int Y) {
      return m_ClusterX2StripPositionX[N-1][X-1][Y-1];
    };
    double GetClusterX2StripPositionY(const int N, const int X, const int Y) {
      return m_ClusterX2StripPositionY[N-1][X-1][Y-1];
    };
    double GetClusterX2StripPositionZ(const int N, const int X, const int Y) {
      return m_ClusterX2StripPositionZ[N-1][X-1][Y-1];
    };

    double GetClusterY2StripPositionX(const int N, const int X, const int Y) {
      return m_ClusterY2StripPositionX[N-1][X-1][Y-1];
    };
    double GetClusterY2StripPositionY(const int N, const int X, const int Y) {
      return m_ClusterY2StripPositionY[N-1][X-1][Y-1];
    };
    double GetClusterY2StripPositionZ(const int N, const int X, const int Y) {
      return m_ClusterY2StripPositionZ[N-1][X-1][Y-1];
    };


    TVector3 GetInnerXPositionOfInteraction(const int i);
    TVector3 GetInnerZPositionOfInteraction(const int i);
    TVector3 GetOuterXPositionOfInteraction(const int i);
    TVector3 GetOuterZPositionOfInteraction(const int i);
    TVector3 GetClusterInnerPositionOfInteraction(const int i);
    TVector3 GetClusterX1PositionOfInteraction(const int i);
    TVector3 GetClusterY1PositionOfInteraction(const int i);
    TVector3 GetClusterX2PositionOfInteraction(const int i);
    TVector3 GetClusterY2PositionOfInteraction(const int i);

    TVector3 GetDetectorNormal(const int i);

  // objects are not written in the TTree
  private:
    TTOGAXSI_SIData*         m_EventData;        //!
    TTOGAXSI_SIData*         m_PreTreatedData;   //!
    TTOGAXSI_SIPhysics*      m_EventPhysics;     //!

  // getters for raw and pre-treated data object
  public:
    TTOGAXSI_SIData* GetRawData()        const {return m_EventData;}
    TTOGAXSI_SIData* GetPreTreatedData() const {return m_PreTreatedData;}

  // parameters used in the analysis
  private:
    int m_NumberOfInnerXDetectors; //!
    int m_NumberOfInnerZDetectors; //!
    int m_NumberOfOuterXDetectors; //!
    int m_NumberOfOuterZDetectors; //!
    int m_NumberOfClusterInnerDetectors; //!
    int m_NumberOfClusterX1Detectors; //!
    int m_NumberOfClusterY1Detectors; //!
    int m_NumberOfClusterX2Detectors; //!
    int m_NumberOfClusterY2Detectors; //!

    vector<vector<vector<double>>> m_InnerXStripPositionX; //!
    vector<vector<vector<double>>> m_InnerXStripPositionY; //!
    vector<vector<vector<double>>> m_InnerXStripPositionZ; //!

    vector<vector<vector<double>>> m_InnerZStripPositionX; //!
    vector<vector<vector<double>>> m_InnerZStripPositionY; //!
    vector<vector<vector<double>>> m_InnerZStripPositionZ; //!

    vector<vector<vector<double>>> m_OuterXStripPositionX; //!
    vector<vector<vector<double>>> m_OuterXStripPositionY; //!
    vector<vector<vector<double>>> m_OuterXStripPositionZ; //!

    vector<vector<vector<double>>> m_OuterZStripPositionX; //!
    vector<vector<vector<double>>> m_OuterZStripPositionY; //!
    vector<vector<vector<double>>> m_OuterZStripPositionZ; //!

    vector<vector<vector<double>>> m_ClusterInnerStripPositionX; //!
    vector<vector<vector<double>>> m_ClusterInnerStripPositionY; //!
    vector<vector<vector<double>>> m_ClusterInnerStripPositionZ; //!

    vector<vector<vector<double>>> m_ClusterX1StripPositionX; //!
    vector<vector<vector<double>>> m_ClusterX1StripPositionY; //!
    vector<vector<vector<double>>> m_ClusterX1StripPositionZ; //!

    vector<vector<vector<double>>> m_ClusterY1StripPositionX; //!
    vector<vector<vector<double>>> m_ClusterY1StripPositionY; //!
    vector<vector<vector<double>>> m_ClusterY1StripPositionZ; //!

    vector<vector<vector<double>>> m_ClusterX2StripPositionX; //!
    vector<vector<vector<double>>> m_ClusterX2StripPositionY; //!
    vector<vector<vector<double>>> m_ClusterX2StripPositionZ; //!

    vector<vector<vector<double>>> m_ClusterY2StripPositionX; //!
    vector<vector<vector<double>>> m_ClusterY2StripPositionY; //!
    vector<vector<vector<double>>> m_ClusterY2StripPositionZ; //!

    // thresholds
    double m_E_RAW_Threshold; //!
    double m_E_Threshold;     //!
  
  private:
    unsigned int m_MaximumStripMultiplicityAllowed; //!
    double m_StripEnergyMatching; //!
  
  // spectra class
  private:
    TTOGAXSI_SISpectra* m_Spectra; // !

  // spectra getter
  public:
    map<string, TH1*>   GetSpectra(); 

  // Static constructor to be passed to the Detector Factory
  public:
    static NPL::VDetector* Construct();

    ClassDef(TTOGAXSI_SIPhysics,1)  // TOGAXSI_SIPhysics structure
  
  private: //geometry
    //InnerX Detector
    //Wafer parameter
    double InnerX_Wafer_Length;
    double InnerX_Wafer_Width;
    double InnerX_Wafer_Thickness;
    double InnerX_Wafer_LongitudinalStrips;
    double InnerX_Wafer_TransverseStrips;

    // PCB parameter
    double InnerX_PCB_Thickness;

    double InnerX_PCB_Width;
    double InnerX_PCB_Length;

    //InnerZ Detector
    //Wafer parameter
    double InnerZ_Wafer_Length;
    double InnerZ_Wafer_Width;
    double InnerZ_Wafer_Thickness;
    double InnerZ_Wafer_LongitudinalStrips;
    double InnerZ_Wafer_TransverseStrips;

    // PCB parameter
    double InnerZ_PCB_Thickness;

    double InnerZ_PCB_Width;
    double InnerZ_PCB_Length;


    //OuterX Detector
    //Wafer parameter
    double OuterX_Wafer_Length;
    double OuterX_Wafer_Width;
    double OuterX_Wafer_Thickness;
    double OuterX_Wafer_LongitudinalStrips;
    double OuterX_Wafer_TransverseStrips;

    // PCB parameter
    double OuterX_PCB_Thickness;

    double OuterX_PCB_Width;
    double OuterX_PCB_Length;
    double OuterX_PCB_gap;

    //OuterZ Detector
    //Wafer parameter
    double OuterZ_Wafer_Length;
    double OuterZ_Wafer_Width;
    double OuterZ_Wafer_Thickness;
    double OuterZ_Wafer_LongitudinalStrips;
    double OuterZ_Wafer_TransverseStrips;

    // PCB parameter
    double OuterZ_PCB_Thickness;

    double OuterZ_PCB_Width;
    double OuterZ_PCB_Length;
    double OuterZ_PCB_gap;


    //ClusterInner Detector
    //Wafer parameter
    double ClusterInner_Wafer_Height;
    double ClusterInner_Wafer_Base;
    double ClusterInner_Wafer_Top;
    double ClusterInner_Wafer_Thickness;

    double ClusterInner_ActiveWafer_Height;
    double ClusterInner_ActiveWafer_Base;
    double ClusterInner_ActiveWafer_Top;
    double ClusterInner_ActiveWafer_Thickness;

    double ClusterInner_Wafer_LongitudinalStrips;
    double ClusterInner_Wafer_TransverseStrips;


    //ClusterX Detector
    //Wafer parameter
    double ClusterX_Wafer_Length;
    double ClusterX_Wafer_Width;
    double ClusterX_Wafer_Thickness;
    double ClusterX_Wafer_LongitudinalStrips;
    double ClusterX_Wafer_TransverseStrips;

    // PCB parameter
    double ClusterX_PCB_Thickness;
    double ClusterX_PCB_Width;
    double ClusterX_PCB_Length;
    double ClusterX_PCB_gap;

    //ClusterY Detector
    //Wafer parameter
    double ClusterY_Wafer_Length;
    double ClusterY_Wafer_Width;
    double ClusterY_Wafer_Thickness;
    double ClusterY_Wafer_LongitudinalStrips;
    double ClusterY_Wafer_TransverseStrips;

    // PCB parameter
    double ClusterY_PCB_Thickness;
    double ClusterY_PCB_Width;
    double ClusterY_PCB_Length;
    double ClusterY_PCB_gap;

};
#endif
