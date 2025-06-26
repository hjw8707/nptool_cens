#ifndef TEXOGAMPHYSICS_H
#define TEXOGAMPHYSICS_H
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: S. Giron   contact address: giron@ipno.in2p3.fr          *
 *                  B. Le Crom                  lecrom@ipno.in2p3.fr         *
 * Creation Date  : march 2014                                               *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold exogam treated data                                      *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *****************************************************************************/


// STL
#include <vector>
#include <map>
#include <algorithm>

// NPL
#include "NPCalibrationManager.h"
#include "NPVDetector.h"
#include "NPVTreeReader.h"
#include "TExogamGeo.h"
#include "TExogamData.h"
#include "TExogamCalData.h"
#include "TExogamSpectra.h"
#include "TExogamPhysicsReader.h"
#include "NPInputParser.h"
#include "RootHistogramsCalib.h"
#include "TSpectrum.h"
#include "NPTimeStamp.h"
#include "NPBeam.h"
// ROOT 
#include "TVector2.h" 
#include "TVector3.h" 
#include "TObject.h"
#include "Math/MinimizerOptions.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TF1.h"
//Cubix
#if CUBIX
#include "CXRecalEnergy.h"
#endif

using namespace std ;
// Forward Declaration
class TExogamSpectra;

class TExogamPhysics : public TObject, public NPL::VDetector, public TExogamPhysicsReader{
 public:
  TExogamPhysics()	;
  ~TExogamPhysics();


  TExogamPhysics& operator=(const TExogamPhysics& other) {
        if (this == &other)
            return *this;


        //FIXME For a very weird reason, when using a treereader to read TExogamPhysics
        // deleting TSevent causes a memory leak crash
        // I couldnt find why, when I'm commenting it it works ok, but it makes no sense
        delete m_PreTreatedData;
        delete m_EventData;
        // delete TSEvent;
        E = other.E;
        EHG = other.EHG;
        Outer1 = other.Outer1;
        Outer2 = other.Outer2;
        Outer3 = other.Outer3;
        Outer4 = other.Outer4;
        Flange = other.Flange;
        Crystal = other.Crystal;
        TDC = other.TDC;
        TS = other.TS;
        E_AB = other.E_AB;
        Size_AB = other.Size_AB;
        Flange_AB = other.Flange_AB;
        Crystal_AB = other.Crystal_AB;
        TDC_AB = other.TDC_AB;
        TS_AB = other.TS_AB;
        Outer_AB = other.Outer_AB;
        Theta = other.Theta;
        Phi = other.Phi;

        m_PreTreatedData = new TExogamCalData(*other.m_PreTreatedData);
        m_EventData = new TExogamData(*other.m_EventData);
        m_EventPhysics = this;

        return *this;
    }
 public: 
  void Clear()	              ;	
  void Clear(const Option_t*) {};

  // Only threshold and cal applied to exogam
  std::vector<double> E;
  std::vector<double> EHG;
  std::vector<double> Outer1;
  std::vector<double> Outer2;
  std::vector<double> Outer3;
  std::vector<double> Outer4;
  std::vector<unsigned int> Flange;
  std::vector<unsigned int> Crystal;
  std::vector<double> TDC;
  std::vector<unsigned long long> TS;
  
  // Previous treatment + Add_Back (size of vectors are not the same because of AB !)
  // Energy AddBack
  std::vector<double> E_AB;
  // Number of gammas added to form the Energy AddBack
  std::vector<unsigned int> Size_AB;
  // Flange (addback only done in a specific flange)
  std::vector<unsigned int> Flange_AB;
  // Crystal with highest E
  std::vector<unsigned int> Crystal_AB;
  // TDC AddBack
  std::vector<double> TDC_AB;
  // TS AddBack
  std::vector<unsigned long long> TS_AB;
  // Outer with highest E
  std::vector<int> Outer_AB;
  // Theta Outer with highest E
  std::vector<double> Theta;
  // Phi Outer with highest E
  std::vector<double> Phi;




 
  /* 
  TH1F*                 clover_mult                  ;  
  TH1F*                 cristal_mult                 ;  
  */

 public:		//	Innherited from VDetector Class
			
  //	Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
  void ReadConfiguration(NPL::InputParser) 				; //!
		

  //	Add Parameter to the CalibrationManger
  void AddParameterToCalibrationManager()	;		 //!
			
		
  //	Activated associated Branches and link it to the private member DetectorData address
  //	In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
  void InitializeRootInputRaw() 					; //!

    //   Activated associated Branches and link it to the private member DetectorPhysics address
    //   In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics() ; //!

  //	Create associated branches and associated private member DetectorPhysics address
  void InitializeRootOutput() 		 		; //!
		
  //	This method is called at each event read from the Input Tree. Aim is to build a tree of calibrated data.
  void PreTreat()			       ; //!

  //	This method is called at each event read from the Input Tree. 
  void BuildPhysicalEvent()					; //!
		
		
  //	Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
  //	This method aimed to be used for analysis performed during experiment, when speed is requiered.
  //	NB: This method can eventually be the same as BuildPhysicalEvent.
  void BuildSimplePhysicalEvent()	       ; //!

  bool TDCMatch(unsigned int event); //!

  //	Those two method all to clear the Event Physics or Data
  void ClearEventPhysics()		{Clear();}		 //!
  void ClearEventData()			{m_EventData->Clear();}	 //!
  void ClearPreTreatedData()	        {m_PreTreatedData->Clear();} //!

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra(); //!
    // Fill the spectra hold by the spectra class
    void FillSpectra(); //!
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra(); //!
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra(); //!
    
    void SetTreeReader(TTreeReader* TreeReader); //!

    double DopplerCorrection(double E, double Angle, double Beta);//!
    double DopplerCorrection(double E, double ThetaGamma, double PhiGamma, double ThetaPart, double PhiPart, double Beta);//!
  
  
    void InitializeRootHistogramsCalib(); //!
    
    void FillHistogramsCalib(); //!
    
    void DoCalibration();//!

    void InitializeRootHistogramsE_F(unsigned int Detector_Nbr); //!

    void InitializeRootHistogramsEHG_F(unsigned int Detector_Nbr); //!

    // void InitializeRootHistogramsT_F(unsigned int Detector_Nbr); //!
    
    void InitializeRootHistogramsOuter_F(unsigned int Detector_Nbr, unsigned int Outer_Nbr); //!
    
    void FillRootHistogramsCalib_F(); //!
  
    void DoCalibrationE_F(unsigned int  Detector_Nbr,std::string CalibType, ofstream* calib_file, ofstream* dispersion_file, unsigned int Threshold); //!
    
    void DoCalibrationEHG_F(unsigned int  Detector_Nbr, ofstream* calib_file, ofstream* dispersion_file){}; //!
    
    void DoCalibrationOuter_F(unsigned int  Detector_Nbr, unsigned int Outer_Nbr, ofstream* calib_file, ofstream* dispersion_file){}; //!
      
    void MakeFolder(std::string make_folder); //!
    
    void MakeECalibFolders(std::string make_folder); //!
    
    void MakeEHGCalibFolders(std::string make_folder); //!
    
    void MakeOuterCalibFolders(std::string make_folder); //!

    void DefineCalibrationSource(std::string source); //!

    void CreateCalibrationEFiles(ofstream* calib_file,ofstream* dispersion_file); //!
    
    void CreateCalibrationEHGFiles(ofstream* calib_file,ofstream* dispersion_file); //!
    
    void CreateCalibrationOuterFiles(ofstream* calib_file,ofstream* dispersion_file); //!

    void ReadDoCalibration(NPL::InputParser parser); //!

    void WriteHistogramsCalib(); //!
    
    void WriteHistogramsEfficiency(); //!
    
    void WriteHistogramsE(); //!

    
    void InitializeRootHistogramsEfficiency(); //!

    void Efficiency();
    
    void FitFunction(TH1* Histogram, int PeakId, double* Area);

    void InitFitParameters();

    void EvaluateEfficiency_F(std::string DetectorID, ofstream* efficiency_file, double Source_activity); //!
    
    void ReadEfficiency(NPL::InputParser parser); //!
    
    void CreateEfficiencyFile(ofstream* efficiency_file); //!

    unsigned long long FindFirstExoTS(); //!
    
    unsigned long long FindLastExoTS(); //!
    
    void FillHistogramsEfficiency(); //!

    void FindPeaks(TH1* Hist);//!

    int FindMatchingEnergy(double Energy);//!

    // Setting TS Reference
    void SetRefTS(std::string TSRef_Name, unsigned long long TSRef);//!

    void ReadConfigurationTS();//! 
   
    Vector3D ConvertVector(const TVector3& v){
    return Vector3D(v.X(), v.Y(), v.Z());
  };//!
  
  private:
    void ClaimReaderData(); //!

 private:	//	Root Input and Output tree classes

 				
  TExogamData*         m_EventData;        //!
  TExogamCalData*      m_PreTreatedData;   //!
  TExogamPhysics*      m_EventPhysics;     //!
  

 public:		//	Specific to EXOGAM Array
  //	Add a Clover
  // void AddClover(string AngleFile);
  void AddClover(int Board, int Flange, int Channel0, int Channel1); //!

 
  // Give and external TMustData object to TExogamPhysics. Needed for online analysis for example.
  void SetRawDataPointer(TExogamData* rawDataPointer) {m_EventData = rawDataPointer;} //!
  // Retrieve raw and pre-treated data
  TExogamData* GetRawData()        const {return m_EventData;}; //!
  TExogamCalData* GetPreTreatedData() const {return m_PreTreatedData;}; //!
  
  void ResetPreTreatVariable(); //!

  void ReadAnalysisConfig(); //!

  double ComputeMeanFreePath(double GammaEnergy); //!

  unsigned int GetFlangeNbr(unsigned int crystal_nbr); //!

  int GetMaxOuter(unsigned int EventId); //!
  
  double GetDoppler(double Energy, unsigned int Flange, unsigned int Crystal, unsigned int Outer); //!

  // Newton Raphson setter and getter

  void SetNRThreshold(const double &nrThreshold_){nrThreshold = nrThreshold_;}
  void SetNRMaxIter(const int &nrMaxIter_){nrMaxIter = nrMaxIter_;}
  void SetInitTheta(const double &initTheta_ ){initTheta = initTheta_;}

  private: // Variables for analysis

  unsigned int m_EXO_Mult;//!
  double m_EXO_OuterUp_RAW_Threshold;//!
  double m_EXO_E_RAW_Threshold;//!
  double m_EXO_E_Threshold;//!
  double m_EXO_EHG_RAW_Threshold;//!
  double m_EXO_TDC_RAW_Threshold;//!
  double m_ExoTDC_LowThreshold;//!
  double m_ExoTDC_HighThreshold;//!
  int m_NumberOfClovers;//!
  double EXO_E;//!
  double EXO_EHG;//!
  double EXO_TDC;//!
  double EXO_Outer1;//!
  double EXO_Outer2;//!
  double EXO_Outer3;//!
  double EXO_Outer4;//!
  unsigned int flange_nbr;//!
  unsigned int crystal_nbr;//!
  vector<TVector3> m_pos_segment[4][4];//!
  vector<int> m_flange;//!
  map<unsigned int, int> MapFlangeToCloverNumber;//!
  TExogamGeo* ExogamGeo;//!
  std::map<unsigned int, double> MapFlangeRadius;//!
  std::map<unsigned int, std::map<unsigned int,  std::map<unsigned int, Vector3D>>> MapFlangeCoordinates;//!
  
  std::map<unsigned int,std::pair<unsigned int,unsigned int>> MapCrystalFlangeCLover;//! Map key is raw crystal nbr, pair associated is flange nbr and crystal nbr in the flange
  
  std::vector<unsigned int> IgnoreTDC;//!

  double GeDensity = 0.005323; //! g/mm3
  std::map<double, double> Map_PhotonCS;//!
 
  map<int, bool> DoCalibrationE;                                    //!
  map<int, bool> DoCalibrationEHG;                                  //!
  map<int, bool> DoCalibrationT;                                    //!
  map<int, map<int,bool>> DoCalibrationOuter;                       //!

  unsigned int Threshold_E_Cal;//!
  unsigned int Threshold_EHG_Cal;//!
  unsigned int Threshold_Outers_Cal;//!
  
  std::vector<std::string> Source_isotope;                          //!
  std::vector<double> Source_E;                                     //!
  std::vector<double> Source_Sig;                                   //!
  std::vector<double> Source_branching_ratio;                       //!
  std::vector<double> Source_branching_ratio_err;                       //!
  std::string Source_name;//!
  double Source_activity;//!
  double Time;//!
  unsigned int FitPolOrder;//!
  
#if CUBIX
  CXRecalEnergy* CubixEnergyCal = new CXRecalEnergy();//!
#endif
  
  // Efficiency variables
  map<int, bool> DoEfficiency;//!
  std::vector<double> Energies;//!
  std::vector<int> SourceID;//!
  std::vector<std::pair<double, double>> MinMax;//!
  double sigma;//!
  double threshold;//!
  int peaks_nb;//!
  double minmax;//!

  //NR variables 

  double initTheta = 0.7;//!
  double nrThreshold = 0.0001;//!
  double nrMaxIter = 20;//!

  private:
  // Fit parameters
  std::string Minimizer;//!
  std::string Algorithm;//!
  double Tolerance;//!
  int PrintLevel;//!
  std::string FitOptions;//!

  Double_t DefFWHM;//!
  Double_t DefFWHM_min;//!
  Double_t DefFWHM_max;//!
  
  Double_t LeftTailVal;//!
  Double_t LeftTailValMin;//!
  Double_t LeftTailValMax;//!

  Double_t RightTailVal;//!
  Double_t RightTailValMin;//!
  Double_t RightTailValMax;//!

  Double_t StepVal;//!
  Double_t StepValMin;//!
  Double_t StepValMax;//!
    
  bool UseLT;//!
  bool UseRT;//!
  bool UseStep;//!

  TF1         *fFitFunction       = nullptr;//!
  TF1         *fBackFunction      = nullptr;//!
  TF1         *fResidue           = nullptr;//!
  TH1         *fHistogram       = nullptr;//!
  
  double DoubleTailedStepedGaussian(Double_t*xx,Double_t*pp);//!
  double StepedBackground(Double_t*xx,Double_t*pp);//!
  double PeakFunction(Double_t*xx,Double_t*pp);//!
  double Residue(Double_t*xx,Double_t*pp);//!
  
  std::string BackgroundType;//!
  
  TFitResultPtr FitResult= nullptr;//!

  private: // Spectra Class   
    TExogamSpectra*      m_Spectra;//! 

  public: // Spectra Getter
    map< string , TH1*> GetSpectra(); 		 //!

  public: // Static constructor to be passed to the Detector Factory
     static NPL::VDetector* Construct(); //!
     static NPL::VTreeReader* ConstructReader(); //!
  
  private:
    TimeStamp* TSEvent;//!
    std::string RefTS_Name = "";//!
    unsigned long long RefTS = 0;//!
  
  private:
    bool DataIsCal;//!
     
     ClassDef(TExogamPhysics,1)  // ExogamPhysics structure
    };

namespace EXOGAM_LOCAL
{
  double fEXO_E(const TExogamData* m_EventData, const unsigned int& i);
  double fEXO_EHG(const TExogamData* m_EventData, const unsigned int& i);
  double fEXO_T(const TExogamData* m_EventData, const unsigned int& i);
  double fEXO_Outer(const TExogamData* m_EventData, const unsigned int& i, const unsigned int OuterNumber);
   const double Threshold_ECC   = 50;
   const double Threshold_GOCCE = 0;
   const double RawThreshold_ECC   = 0;
   const double RawThreshold_GOCCE = 0;
   
   
   //	tranform an integer to a string
   string itoa(int value);
        

 
}


#endif
