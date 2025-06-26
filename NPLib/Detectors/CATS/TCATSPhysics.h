#ifndef __CATSCAL__
#define __CATSCAL__
/*****************************************************************************
 * Copyright (C) 2009-2016   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Sandra Giron  contact address: giron@ipno.in2p3.fr       *
 *                                                                           *
 * Creation Date  : febuary 2010                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold CATS treated data                                        *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/

//   STL
#include <vector>
#include <set>
//   ROOT
#include "TObject.h"
#include "TVector3.h"
#include <TRandom1.h>
#include <TRandom2.h>
#include <TRandom3.h>
//   NPLib
#include "TCATSPhysicsReader.h"
#include "TCATSSpectra.h"
#include "NPVTreeReader.h"
#include "NPVDetector.h"
#include "NPCalibrationManager.h"
#include "NPDetectorFactory.h"
#include "NPInputParser.h"

// forward declaration
class TCATSSpectra;
using namespace std ;

class TCATSPhysics : public TObject, public NPL::VDetector, public TCATSPhysicsReader
{

  public:   //   Constructor and Destructor
    TCATSPhysics();
    ~TCATSPhysics();

  private:   //   Root Input and Output tree classes

    TCATSData*       m_EventData;//!
    TCATSData*       m_PreTreatedData;//!
    TCATSPhysics*    m_EventPhysics;//!

  public :
    vector<UShort_t>    DetNumber;
    vector<double>      PositionX;
    vector<double>      PositionY;
    vector<double>      PositionZ;
    vector<double>      StripNumberX;
    vector<double>      StripNumberY;
    vector<double>      ChargeX;
    vector<double>      ChargeY;
    double           PositionOnTargetX;
    double           PositionOnTargetY;
    double           PositionOnMask1X;
    double           PositionOnMask1Y;
    double           PositionOnMask2X;
    double           PositionOnMask2Y;
    double           m_Zproj;

    TVector3         BeamDirection;//!

    

  private :
    std::map<UShort_t, std::vector<std::pair<UShort_t,UShort_t>>> MapX;//!
    std::map<UShort_t, std::vector<std::pair<UShort_t,UShort_t>>> MapY;//!
    std::map<UShort_t,std::pair<UShort_t,UShort_t>> MaxQX;//!
    std::map<UShort_t,std::pair<UShort_t,UShort_t>> MaxQY;//!
    std::map<UShort_t,Double_t> QSumX;//!
    std::map<UShort_t,Double_t> QSumY;//!
    std::map<UShort_t, std::pair<Double_t,std::pair<Double_t,Double_t>>> Positions;//!
    std::map<UShort_t,vector< vector<double> > >   StripPositionX;//!
    std::map<UShort_t,vector< vector<double> > >   StripPositionY;//!
    std::map<UShort_t,double>             StripPositionZ;//!
    int m_NumberOfCATS; //!
    double m_TargetAngle; //!
    double m_TargetThickness; //!

    // Those two vector contain a pointer to the reconstruction function that should be used for each detector 
    // Method take as argument a vector<double> representing the Charge distribution and uint giving the Strip with Max Q
    // It returns a double, corresponding to the position in Strip unit
    // Methods are declared in the CATS namespace outside the class
    // CINT can't handle such complex vector, so the ifndef part is there to hide it away 
    #ifndef __CINT__
    std::map<UShort_t,double(*)(std::pair<UShort_t,UShort_t>&,std::vector<std::pair<UShort_t,UShort_t>>&,Double_t)> ReconstructionFunctionY;//!
    std::map<UShort_t,double(*)(std::pair<UShort_t,UShort_t>&,std::vector<std::pair<UShort_t,UShort_t>>&,Double_t)> ReconstructionFunctionX;//!
    #endif /* __CINT __ */
  
  public:
    // Set the reconstruction Method used for the CATS plane
    void SetReconstructionMethod(unsigned int CATSNumber, string XorY, string MethodName); //!

  private : 
    //   Map of activated channel
    map< int, vector<bool> > m_XChannelStatus;//!
    map< int, vector<bool> > m_YChannelStatus;//! 
    //   Map of inverted channel
    map< int, vector<int> > m_CATSXInversion;//!
    map< int, vector<int> > m_CATSYInversion;//! 

  public:   // Output data of interest
    //   for a CATS
    void SetTargetAngle(double m_TargetAngle) {m_TargetAngle = m_TargetAngle;} //!
    void SetTargetThickness(double m_TargetThickness) {m_TargetThickness = m_TargetThickness;} //!


    //   Remove bad channel, calibrate the data and apply threshold
    void PreTreat();

    //   Activated associated Branches and link it to the private member DetectorData address
    //   In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw() ;

    //   Activated associated Branches and link it to the private member DetectorPhysics address
    //   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ;

    //   Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput() ;

    //   Clear The PreTeated object
    void ClearPreTreatedData()   {m_PreTreatedData->Clear();}
    TCATSData* GetPreTreatedData()   {return m_PreTreatedData;}

    void BuildPhysicalEvent();

    void BuildSimplePhysicalEvent();

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra();
    // Fill the spectra hold by the spectra class
    void FillSpectra();
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra();
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra();

    void TestCalibration();

    void SetTreeReader(TTreeReader* TreeReader);

    //   Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {m_EventData->Clear();}    

    void  Clear();
    void  Clear(const Option_t*) {};  

    // Give and external TCATSData object to TCATSPhysics, needed for online analysis
    void SetRawDataPointer(void* rawDataPointer) {m_EventData = (TCATSData*)rawDataPointer;}

    //   Return false if the channel is disabled by user
    bool IsValidChannel(const string& DetectorType, const int& Detector , const int& channel);
    void InitializeStandardParameter();
    void AddParameterToCalibrationManager();
    void ReadAnalysisConfig();
    void ReadConfiguration(NPL::InputParser);
    void AddCATS(TVector3 C_X1_Y1, TVector3 C_X28_Y1, TVector3 C_X1_Y28, TVector3 C_X28_Y28, UShort_t N);
    void AddMask(Double_t Z, UShort_t N);

  public:
    TVector3 GetBeamDirection();
    TVector3 GetPositionOnTarget();

    double GetPositionOnTargetX()  {return PositionOnTargetX;}  
    double GetPositionOnTargetY()  {return PositionOnTargetY;}

  public:
    /////////////////////           GETTERS           ////////////////////////
    // X
    inline UShort_t	GetCATSMult()		              const {return DetNumber.size();}
    inline UShort_t	GetCATSDet(const Int_t& i)		const {return DetNumber[i];}
    
    inline UShort_t	GetCATSStripX(const Int_t& i)   const {return StripNumberX[i];}
    inline double	GetCATSChargeX(const Int_t& i)	const {return ChargeX[i];}
    inline double	GetCATSPosX(const Int_t& i)	const {return PositionX[i];} 
    
    inline UShort_t	GetCATSStripY(const Int_t& i)   const {return StripNumberY[i];}
    inline double	GetCATSChargeY(const Int_t& i)	const {return ChargeY[i];}
    inline double	GetCATSPosY(const Int_t& i)	const {return PositionY[i];} 
    


  private:
    unsigned int sizeY;//! 
    unsigned int sizeX;//! 
    unsigned int sizeDet;//!
    Double_t Mask1_Z = 0;//!
    Double_t Mask2_Z = 0;//!

    int counter = 0;//!
    UShort_t StrX;//! 
    UShort_t NX;//! 
    UShort_t CATS_X_Q;//! 
    UShort_t StrY;//! 
    UShort_t NY;//! 
    UShort_t CATS_Y_Q;//! 
  // Look how many CATS were fired
  // use a set to identify which detector has been hit
  std::set<int> DetectorHitX;//! 
  std::set<int> DetectorHit;//! 
  std::map<UShort_t, std::vector<double>>	 QsumSample;//!

  //Debugging
  //unsigned long counter;//!
  unsigned long long time_elapsed1;//!
  unsigned long long time_elapsed2;//!
  unsigned long long time_elapsed3;//!
  
  private: // Spectra Class   
    TCATSSpectra*      m_Spectra;//! 

  public: 
    void WriteSpectra() ;//!

  public: // Spectra Getter
    map< string , TH1*> GetSpectra();

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    static NPL::VTreeReader* ConstructReader();
     ClassDef(TCATSPhysics,1)  // CATSPhysics structure
};


namespace CATS_LOCAL{
  // The reconstruction methods need to be outside the class for practical purpose
  // of dealing with the function pointer
  double AnalyticHyperbolicSecant(std::pair<UShort_t,UShort_t>&,std::vector<std::pair<UShort_t,UShort_t>>&,Double_t);
  double FittedHyperbolicSecant(std::pair<UShort_t,UShort_t>&,std::vector<std::pair<UShort_t,UShort_t>>&,Double_t);
  //double AnalyticGaussian(std::pair<UShort_t,UShort_t>&,std::vector<std::pair<UShort_t,UShort_t>>&,Double_t,  Double_t, UShort_t);
  //double Centroide(std::pair<UShort_t,UShort_t>&,std::vector<std::pair<UShort_t,UShort_t>>&,Double_t,  Double_t, UShort_t);


  //   tranform an integer to a string
  string itoa(int value);

  double fCATS_X_Q(const TCATSData* Data, const int& i);
  double fCATS_Y_Q(const TCATSData* Data, const int& i);
  bool fCATS_Threshold_X(const TCATSData* Data, const int& i);
  bool fCATS_Threshold_Y(const TCATSData* Data, const int& i);
  double fCATS_Ped_X(const TCATSData* m_EventData, const int& i);
  double fCATS_Ped_Y(const TCATSData* m_EventData, const int& i);
}

#endif
