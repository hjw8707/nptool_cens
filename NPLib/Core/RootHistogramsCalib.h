#ifndef ROOTHISTOCALIB_HH
#define ROOTHISTOCALIB_HH
/*****************************************************************************
 * Copyright (C) 2009-2021   this file is part of the NPTool Project         *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: N. de Sereville  contact address: deserevi@ipno.in2p3.fr *
 *                                                                           *
 * Creation Date  : 21/07/09                                                 *
 * Last update    : 10/05/21                                                 *
 *---------------------------------------------------------------------------*
 * Decription: This class is a singleton class which deals with the ROOT     *
 *             output file and tree both for NPSimulation and NPAnalysis.    *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *   + 03/02/11: Add support for TAsciiFile objects (N. de Sereville)        *
 *   + 10/05/21: Add support for split tree output (A. Matta)                *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <map>

// NPL headers
#include "TAsciiFile.h"

// ROOT headers
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TCutG.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include <string>

class RootHistogramsCalib{
public:
  // The analysis class is designed to be a singleton (i.e. only one instance
  // can exist). A member function called Instance is defined, which allows
  // the user to get a pointer to the existing instance or to create it if
  // it does not yet exist:
  // (see the constructor for an explanation of the arguments)
  static RootHistogramsCalib* getInstance(std::string fileNameBase = "Calibration");

  // The analysis class instance can be deleted by calling the Destroy
  // method (NOTE: The class destructor is protected, and can thus not be
  // called directly):
  static void Destroy();
  
protected:
  // Constructor (protected)
  RootHistogramsCalib(std::string fileNameBase);
 
  // Destructor (protected)
  virtual ~RootHistogramsCalib();
  
  // Prevent copying
  RootHistogramsCalib(const RootHistogramsCalib& only);
  const RootHistogramsCalib& operator=(const RootHistogramsCalib& only);
  
private:
  // The static instance of the RootHistogramsCalib class:
  static RootHistogramsCalib* instance;
  
private:
  void InitAsciiFiles();
  void CreateFile(std::string name);
public:
  TFile*      GetFile(std::string name="global") ; 
  std::map<TString,std::map<TString,TH1*>>*  GetTH1Map(){return TH1Map;}; 
  std::map<TString,std::map<TString,TH2*>>*  GetTH2Map(){return TH2Map;}; 
  std::map<TString,std::map<TString,TGraphErrors*>>*  GetTGraphMap(){return TGraphMap;}; 
  std::map<TString,std::map<TString,TCutG*>>*  GetTCutGMap(){return TCutGMap;}; 
  std::map<TString,std::map<TString,TFile*>>*  GetTFileMap(){return TFileMap;}; 
  std::map<TString,std::map<TString,TCanvas*>>* GetTCanvasMap(){return TCanvasMap;};
  void        Fill();
  
private:
  std::string                   pBaseName;
  std::string                   pMasterFile;
  TDirectory*                   pCurrentDirectory;
  TFile* pRootFile = nullptr;
  std::map<TString,std::map<TString,TH1*>>* TH1Map = new std::map<TString,std::map<TString,TH1*>>;
  std::map<TString,std::map<TString,TH2*>>* TH2Map = new std::map<TString,std::map<TString,TH2*>>;
  std::map<TString,std::map<TString,TGraphErrors*>>* TGraphMap = new std::map<TString,std::map<TString,TGraphErrors*>>;
  std::map<TString,std::map<TString,TCutG*>>* TCutGMap = new std::map<TString,std::map<TString,TCutG*>>;
  std::map<TString,std::map<TString,TFile*>>* TFileMap = new std::map<TString,std::map<TString,TFile*>>;
  std::map<TString,std::map<TString,TCanvas*>>* TCanvasMap = new std::map<TString,std::map<TString,TCanvas*>>;
};

#endif // ROOTHISTOCALIB_HH
