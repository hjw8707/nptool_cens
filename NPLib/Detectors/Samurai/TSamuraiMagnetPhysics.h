#ifndef TSamuraiMagnetPhysics_H
#define TSamuraiMagnetPhysics_H
/*****************************************************************************
 * Copyright (C) 2009-2020    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: matta@lpccaen.in2p3.fr    *
 *                                                                           *
 * Creation Date  : October 2021                                             *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class is a front to load the lib when the magnet is called, nothing *
 *  is done                                                                  *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *  
 *                                                                           *
 *                                                                           *
 *****************************************************************************/
// STL
#include <vector>
#include <map>
#include <string>

// NPL
#include "NPVDetector.h"
#include "NPInputParser.h"
#include "NPXmlParser.h"

#if __cplusplus > 199711L && NPMULTITHREADING 
#include "NPDCReconstructionMT.h"
#else
#include "NPDCReconstruction.h"
#endif

// ROOT 
#include "TVector3.h" 

// Forward declaration
//class TSamuraiMagnetSpectra;

class TSamuraiMagnetPhysics : public TObject, public NPL::VDetector{
  public:
    TSamuraiMagnetPhysics();
    ~TSamuraiMagnetPhysics() {};

  public: 
    void Clear();   
    void Clear(const Option_t*) {};

    // Read stream at ConfigFile to pick-up parameters of detector (Position,...) using Token
    void ReadConfiguration(NPL::InputParser) ;


    // Add Parameter to the CalibrationManger
    void AddParameterToCalibrationManager() {;};      

    // Activated associated Branches and link it to the private member DetectorData address
    // In this method mother Branches (Detector) AND daughter leaf (fDetector_parameter) have to be activated
    void InitializeRootInputRaw(){;} ;

    // Activated associated Branches and link it to the private member DetectorPhysics address
    // In this method mother Branches (Detector) AND daughter leaf (parameter) have to be activated
    void InitializeRootInputPhysics(){;} ;

    // Create associated branches and associated private member DetectorPhysics address
    void InitializeRootOutput(){;} ;

    // This method is called at each event read from the Input Tree. Aime is to build treat Raw dat in order to extract physical parameter. 
    void BuildPhysicalEvent(){;} ;

    // Same as above, but only the simplest event and/or simple method are used (low multiplicity, faster algorythm but less efficient ...).
    // This method aimed to be used for analysis performed during experiment, when speed is requiered.
    // NB: This method can eventually be the same as BuildPhysicalEvent.
    void BuildSimplePhysicalEvent() {;};

    // Same as above but for online analysis
    void BuildOnlinePhysicalEvent()  {BuildPhysicalEvent();};

    // Those two method all to clear the Event Physics or Data
    void ClearEventPhysics() {Clear();}      
    void ClearEventData()    {Clear();}   

    // Method related to the TSpectra classes, aimed at providing a framework for online applications
    // Instantiate the Spectra class and the histogramm throught it
    void InitSpectra(){;};
    // Fill the spectra hold by the spectra class
    void FillSpectra(){;};
    // Used for Online mainly, perform check on the histo and for example change their color if issues are found
    void CheckSpectra(){;};
    // Used for Online only, clear all the spectra hold by the Spectra class
    void ClearSpectra(){;};
    // Write Spectra to file
    void WriteSpectra(){;};

  private: // Spectra Class
    // TSamuraiMagnetSpectra* m_Spectra; // !

  public: // Spectra Getter
    std::map< std::string , TH1*> GetSpectra(); 

  public: // Static constructor to be passed to the Detector Factory
    static NPL::VDetector* Construct();
    ClassDef(TSamuraiMagnetPhysics,1)  // SamuraiMagnetPhysics structure
};

#endif
